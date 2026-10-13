# Step 00 - Data audit.
# 1. Download GEO series matrices + raw CEL files.
# 2. Tabulate what clinical metadata each series actually provides (age, sex, platelets, treatment,
#    timing) - this decides what the benchmark can adjust for.
# 3. Record platform technology: none of these are sequence data, so none are variant-callable.
# 4. Checksum every raw array: GSE46922 and GSE56232 (same lab) share physical array files under
#    different patient numbers. Duplicates are collapsed before any pooled analysis.

source("R/utils.R")
suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })

geo_dir <- file.path(DATA, "geo")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
all_gse <- c(cfg$geo$gpl570, cfg$geo$cdna)

get_gse <- function(id) {
  f <- file.path(geo_dir, paste0(id, "_eset.rds"))
  if (file.exists(f)) return(readRDS(f))
  msg("Fetching series matrix ", id)
  es <- getGEO(id, GSEMatrix = TRUE, destdir = geo_dir, getGPL = TRUE)[[1]]
  saveRDS(es, f)
  es
}
esets <- setNames(lapply(all_gse, get_gse), all_gse)

platform_info <- function(gpl) {
  url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief", gpl)
  txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) character())
  field <- function(k) { v <- grep(paste0("^!Platform_", k, " = "), txt, value = TRUE)
    if (length(v)) sub(".* = ", "", v[1]) else NA_character_ }
  data.table(platform = gpl, title = field("title"), technology = field("technology"))
}

# ---------------------------------------------------------------- metadata completeness
kw <- list(
  age = "\\bage\\b|years|yrs", sex = "\\bsex\\b|gender|\\bmale\\b|female",
  platelet = "platelet|plt", treatment = "steroid|prednis|ivig|immunoglobulin|treat|therapy|rituximab|eltrombopag",
  time_since_dx = "days|weeks|months? (of|after|since)|within \\d+ months|after \\d+ months",
  followup_outcome = "recover|remission|resolv|did not recover|follow")
audit <- rbindlist(lapply(all_gse, function(id) {
  pd <- pData(esets[[id]])
  char_cols <- grep("characteristics|description|source_name", names(pd), value = TRUE)
  text <- tolower(unlist(pd[, char_cols, drop = FALSE]))
  gpl <- annotation(esets[[id]])
  pi <- platform_info(gpl)
  data.table(
    gse = id, platform = gpl, platform_title = pi$title, technology = pi$technology,
    n_samples = nrow(pd),
    characteristics_fields = paste(unique(sub(":.*", "", unlist(pd[, grep("characteristics", names(pd)), drop = FALSE]))), collapse = "; "),
    has_age = any(grepl(kw$age, text)), has_sex = any(grepl(kw$sex, text)),
    has_platelet_count = any(grepl(kw$platelet, text)), has_treatment = any(grepl(kw$treatment, text)),
    has_timing = any(grepl(kw$time_since_dx, text)), has_outcome = any(grepl(kw$followup_outcome, text)),
    sequence_reads = FALSE,
    gatk_variant_callable = FALSE,
    reason = "intensity/ratio microarray: measures mRNA abundance, never reads DNA sequence")
}))
write_tab(audit, "audit_metadata_completeness")

# ---------------------------------------------------------------- GPL570 sample sheet + raw files
parse_group <- function(x) {
  x <- tolower(x)
  fifelse(grepl("remission", x), "remission",
          fifelse(grepl("newly|new itp", x), "ND",
                  fifelse(grepl("chronic", x), "chronic", NA_character_)))
}
cel_dir <- file.path(geo_dir, "cel")
dir.create(cel_dir, showWarnings = FALSE)

sheet <- rbindlist(lapply(cfg$geo$gpl570, function(id) {
  pd <- pData(esets[[id]])
  data.table(gsm = pd$geo_accession, gse = id, title = pd$title,
             group = parse_group(pd$title),
             patient = sub(".*?(\\d+)\\s*$", "\\1", pd$title),
             url = sub("^ftp://", "https://", pd$supplementary_file))
}))
stopifnot(!anyNA(sheet$group))
# Subjects: GSE46922 ND and chronic are different children; GSE56232 is paired ND -> remission.
sheet[, subject := fifelse(gse == "GSE56232", paste0("GSE56232_P", patient),
                           paste0(gse, "_", group, "_", patient))]

sheet[, file := file.path(cel_dir, paste0(gsm, ".CEL"))]
for (i in seq_len(nrow(sheet))) {
  if (file.exists(sheet$file[i])) next
  gz <- download_once(sheet$url[i], paste0(sheet$file[i], ".gz"))
  R.utils::gunzip(gz, destname = sheet$file[i], remove = TRUE, overwrite = TRUE)
}
sheet[, md5 := unname(tools::md5sum(file))]

# ---------------------------------------------------------------- duplicate detection
sheet[, n_copies := .N, by = md5]
sheet[, canonical_gsm := gsm[1], by = md5]       # series order puts GSE46922 first
sheet[, is_duplicate := gsm != canonical_gsm]
dups <- sheet[n_copies > 1, .(gsm, gse, title, group, md5 = substr(md5, 1, 12), canonical_gsm)]
dup_pairs <- merge(sheet[is_duplicate == TRUE, .(gsm, gse, title, canonical_gsm)],
                   sheet[, .(canonical_gsm = gsm, canonical_gse = gse, canonical_title = title)],
                   by = "canonical_gsm")
write_tab(dup_pairs, "audit_duplicate_arrays")
msg(nrow(dup_pairs), " raw array files are byte-identical copies of arrays in another series")

# A duplicated ND array means the same child appears in both series: give it one subject ID so the
# GSE56232 remission sample of that child stays paired correctly in pooled analyses.
alias <- unique(merge(sheet[, .(gsm, canonical_gsm, subject)],
                      sheet[, .(canonical_gsm = gsm, canonical_subject = subject)], by = "canonical_gsm")
                [subject != canonical_subject, .(subject, canonical_subject)])
sheet[, subject_unified := subject]
if (nrow(alias)) sheet[alias, subject_unified := i.canonical_subject, on = "subject"]

save_rds(sheet, "gpl570_sample_sheet")
write_tab(sheet[, .(gsm, gse, title, group, patient, subject, subject_unified, md5, is_duplicate, canonical_gsm)],
          "gpl570_sample_sheet")

# ---------------------------------------------------------------- GSE23754 sample sheet
pd <- pData(esets[[cfg$geo$cdna]])
state <- tolower(pd[["disease state:ch1"]])
cdna_sheet <- data.table(
  gsm = pd$geo_accession, title = pd$title,
  group = fifelse(grepl("progress", state), "progression",
                  fifelse(grepl("acute|self", state), "acute",
                          fifelse(grepl("chronic", state), "chronic", NA_character_))),
  patient = pd[["patient:ch1"]], description = pd$description)
stopifnot(!anyNA(cdna_sheet$group))
# "eventually chronic" = progression-to-chronic (sampled early) + chronic (sampled late)
cdna_sheet[, eventual_chronic := as.integer(group != "acute")]
cdna_sheet[, early_sample := group %in% c("acute", "progression")]
save_rds(cdna_sheet, "gse23754_sample_sheet")
write_tab(cdna_sheet, "gse23754_sample_sheet")

summary_tab <- rbind(
  sheet[, .(n = .N, n_unique_arrays = uniqueN(md5)), by = .(gse, group)],
  cdna_sheet[, .(gse = cfg$geo$cdna, n = .N, n_unique_arrays = .N), by = group])
write_tab(summary_tab, "audit_group_sizes")
print(summary_tab)
msg("Audit complete")
