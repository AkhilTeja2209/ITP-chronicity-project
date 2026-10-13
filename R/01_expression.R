# Step 01 - Expression preprocessing, re-annotation and confounder estimation.
#
# GPL570 (GSE46922 + GSE56232, sorted T cells): RMA on the de-duplicated raw CEL files, array QC,
#   probe re-annotation against the current Ensembl release (probe-sequence alignments are redone by
#   Ensembl each release) compared with the legacy hgu133plus2.db mapping.
# GPL10845 (GSE23754, whole blood, two-colour): GEO-processed log2(Cy5/Cy3) ratios; the 2006-era
#   clone annotation is updated through Entrez ID / HGNC alias tables. The raw .gpr files cannot be
#   reprocessed because their spot identifiers are blank (no GAL file was deposited).
# Both: expression-derived confounder scores (platelet, glucocorticoid, interferon, cell types, sex).

source("R/utils.R")
suppressPackageStartupMessages({
  library(Biobase); library(affy); library(limma)
  library(AnnotationDbi); library(org.Hs.eg.db); library(hgu133plus2.db)
})

ann_dir <- file.path(DATA, "annotation")
dir.create(ann_dir, showWarnings = FALSE)

# ---------------------------------------------------------------- Ensembl annotation (cached)
biomart_query <- function(attrs, cache) {
  f <- file.path(ann_dir, cache)
  if (file.exists(f)) return(fread(f))
  for (m in cfg$annotation$ensembl_mirror) {
    res <- tryCatch({
      msg("Querying Ensembl BioMart (", m, ") for ", cache)
      mart <- if (m == "www") biomaRt::useEnsembl("genes", "hsapiens_gene_ensembl")
              else biomaRt::useEnsembl("genes", "hsapiens_gene_ensembl", mirror = m)
      as.data.table(biomaRt::getBM(attrs, mart = mart))
    }, error = function(e) { msg("  BioMart failed: ", conditionMessage(e)); NULL })
    if (!is.null(res) && nrow(res) > 0) { fwrite(res, f, sep = "\t"); return(res) }
  }
  NULL
}

genes_bm <- biomart_query(c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                            "chromosome_name", "start_position", "end_position"), "ensembl_genes.tsv")
if (is.null(genes_bm)) {
  msg("BioMart unavailable - falling back to EnsDb.Hsapiens.v86 gene coordinates")
  suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
  g <- ensembldb::genes(EnsDb.Hsapiens.v86)
  genes_bm <- data.table(ensembl_gene_id = g$gene_id, external_gene_name = g$gene_name,
                         gene_biotype = g$gene_biotype, chromosome_name = as.character(GenomicRanges::seqnames(g)),
                         start_position = GenomicRanges::start(g), end_position = GenomicRanges::end(g))
  annotation_source <- "EnsDb.Hsapiens.v86 (fallback)"
} else annotation_source <- "Ensembl BioMart (current release)"
gene_lookup <- unique(genes_bm[chromosome_name %in% c(1:22, "X", "Y", "MT")],
                      by = "ensembl_gene_id")
setnames(gene_lookup, c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name",
                        "start_position", "end_position"),
         c("gene_id", "symbol", "biotype", "chr", "start", "end"))
gene_lookup[symbol == "", symbol := gene_id]
save_rds(gene_lookup, "gene_lookup")

# Symbol-keyed copy of an Ensembl-keyed matrix (for marker-set scoring by symbol).
to_symbol <- function(mat) {
  sym <- gene_lookup$symbol[match(rownames(mat), gene_lookup$gene_id)]
  collapse_to_genes(mat, sym)
}

# ---------------------------------------------------------------- confounder scores
confounder_scores <- function(mat_sym) {
  sets <- setdiff(names(gene_sets), c("sex_female", "sex_male"))
  sc <- as.data.table(lapply(setNames(sets, sets), function(s) score_set(mat_sym, gene_sets[[s]])))
  sex_score <- score_set(mat_sym, gene_sets$sex_male, 1) - score_set(mat_sym, gene_sets$sex_female, 1)
  sc[, sex_score := sex_score]
  sc[, sex := NA_character_]
  ok <- is.finite(sex_score)
  if (sum(ok) >= 4) {
    km <- kmeans(sex_score[ok], centers = range(sex_score[ok]))
    sep <- abs(diff(km$centers[, 1])) / sqrt(mean(km$withinss / pmax(km$size - 1, 1)))
    lab <- ifelse(km$cluster == which.max(km$centers[, 1]), "M", "F")
    # Two well-separated clusters -> call sex; otherwise (e.g. one sex only) leave unresolved.
    if (is.finite(sep) && sep > 3) sc[ok, sex := lab]
    attr(sc, "sex_separation") <- sep
  }
  sc
}

# ================================================================= GPL570
sheet <- load_rds("gpl570_sample_sheet")
uniq <- sheet[is_duplicate == FALSE]
msg("RMA on ", nrow(uniq), " unique GPL570 arrays (", sum(sheet$is_duplicate), " duplicates removed)")
ab <- ReadAffy(filenames = uniq$file, sampleNames = uniq$gsm)
parse_scan <- function(x) {
  d <- as.Date(x, format = "%m/%d/%y")
  d[is.na(d)] <- as.Date(x[is.na(d)], format = "%Y-%m-%d")
  as.character(d)
}
uniq[, scan_date := parse_scan(protocolData(ab)$ScanDate)]
probe_mat <- exprs(rma(ab))
rm(ab); invisible(gc())

# QC: relative log expression (RLE) and near-duplicate RNA (distinct files, ~identical signal)
rle <- probe_mat - apply(probe_mat, 1, median)
qc <- data.table(gsm = colnames(probe_mat), rle_median = apply(rle, 2, median),
                 rle_iqr = apply(rle, 2, IQR))
lim <- quantile(qc$rle_iqr, 0.75) + 1.5 * IQR(qc$rle_iqr)
qc[, outlier := abs(rle_median) > 0.1 | rle_iqr > lim]
write_tab(merge(uniq[, .(gsm, gse, group, scan_date)], qc, by = "gsm"), "gpl570_array_qc")
r <- cor(probe_mat)
near <- which(r > cfg$qc$near_duplicate_r & upper.tri(r), arr.ind = TRUE)
near_dup <- data.table(gsm_a = rownames(r)[near[, 1]], gsm_b = colnames(r)[near[, 2]], r = r[near])
write_tab(near_dup, "gpl570_near_duplicates")
if (isTRUE(cfg$qc$drop_outliers) && any(qc$outlier)) {
  probe_mat <- probe_mat[, !qc$outlier]; uniq <- uniq[gsm %in% colnames(probe_mat)]
}

# ---- probe re-annotation: legacy Bioconductor package vs current Ensembl alignment
probes <- rownames(probe_mat)
legacy_sym <- mapIds(hgu133plus2.db, probes, "SYMBOL", "PROBEID", multiVals = "first")
legacy_ens <- mapIds(hgu133plus2.db, probes, "ENSEMBL", "PROBEID", multiVals = "first")
affy_bm <- biomart_query(c("affy_hg_u133_plus_2", "ensembl_gene_id"), "ensembl_affy_hg_u133_plus_2.tsv")
if (!is.null(affy_bm)) {
  affy_bm <- affy_bm[affy_hg_u133_plus_2 != "" & ensembl_gene_id %in% gene_lookup$gene_id]
  n_genes_per_probe <- affy_bm[, .(n = uniqueN(ensembl_gene_id), gene = ensembl_gene_id[1]), by = affy_hg_u133_plus_2]
  current_ens <- setNames(rep(NA_character_, length(probes)), probes)
  uni <- n_genes_per_probe[n == 1]
  current_ens[uni$affy_hg_u133_plus_2[uni$affy_hg_u133_plus_2 %in% probes]] <-
    uni$gene[uni$affy_hg_u133_plus_2 %in% probes]
  ambiguous <- n_genes_per_probe[n > 1, affy_hg_u133_plus_2]
  reannotated <- TRUE
} else {
  msg("Current probe mapping unavailable - using legacy hgu133plus2.db mapping")
  current_ens <- legacy_ens; ambiguous <- character(); reannotated <- FALSE
}
probe_ann <- data.table(probe = probes, legacy_symbol = legacy_sym, legacy_gene = legacy_ens,
                        current_gene = current_ens)
probe_ann[, current_symbol := gene_lookup$symbol[match(current_gene, gene_lookup$gene_id)]]
probe_ann[, current_biotype := gene_lookup$biotype[match(current_gene, gene_lookup$gene_id)]]
probe_ann[, status := fcase(
  is.na(legacy_gene) & is.na(current_gene) & probe %in% ambiguous, "ambiguous_now (multi-gene)",
  is.na(legacy_gene) & is.na(current_gene), "unmapped_both",
  is.na(legacy_gene) & !is.na(current_gene), "recovered (new mapping)",
  !is.na(legacy_gene) & is.na(current_gene) & probe %in% ambiguous, "lost: now multi-gene",
  !is.na(legacy_gene) & is.na(current_gene), "lost: no longer maps",
  legacy_gene == current_gene, "unchanged",
  default = "reassigned to different gene")]
write_tab(probe_ann, "gpl570_probe_reannotation")
reannot_summary <- probe_ann[, .N, by = status][order(-N)]
legacy_genes <- unique(na.omit(probe_ann$legacy_gene)); current_genes <- unique(na.omit(probe_ann$current_gene))
new_genes <- setdiff(current_genes, legacy_genes)
biotype_new <- gene_lookup[gene_id %in% new_genes, .N, by = biotype][order(-N)]
write_tab(reannot_summary, "gpl570_reannotation_summary")
write_tab(biotype_new, "gpl570_reannotation_new_gene_biotypes")
msg("Re-annotation: ", length(legacy_genes), " legacy genes -> ", length(current_genes),
    " current genes; ", length(new_genes), " genes newly measurable")

gpl570_expr <- collapse_to_genes(probe_mat, current_ens)
gpl570_sc <- confounder_scores(to_symbol(gpl570_expr))
gpl570_samples <- cbind(uniq[, .(gsm, gse, group, subject = subject_unified, scan_date)], gpl570_sc)

# Per-series views include duplicate GSMs (identical columns) so each series keeps its own design.
series_view <- function(id) {
  s <- sheet[gse == id]
  m <- gpl570_expr[, s$canonical_gsm, drop = FALSE]; colnames(m) <- s$gsm
  sc <- gpl570_samples[match(s$canonical_gsm, gsm), !c("gsm", "gse", "group", "subject")]
  list(expr = m, samples = cbind(s[, .(gsm, gse, group, patient, subject, subject_unified)], sc))
}
per_series <- setNames(lapply(cfg$geo$gpl570, series_view), cfg$geo$gpl570)

# ---- pooled matrix (unique arrays) with batch handling
mod <- model.matrix(~ group, gpl570_samples)
bat <- model.matrix(~ gse, gpl570_samples)[, -1, drop = FALSE]
confounded <- qr(cbind(mod, bat))$rank < ncol(mod) + ncol(bat)
if (!confounded && uniqueN(gpl570_samples$gse) > 1) {
  pooled <- sva::ComBat(gpl570_expr, batch = gpl570_samples$gse, mod = mod)
  batch_note <- "ComBat applied (batch = series, disease group protected)"
} else {
  pooled <- gpl570_expr
  batch_note <- paste("ComBat NOT applied: after de-duplication series is fully confounded with group",
                      "(each group occurs in only one series), so batch cannot be separated from biology.")
}
msg(batch_note)

pca_plot <- function(m, s, title) {
  v <- apply(m, 1, var); top <- m[order(-v)[1:2000], ]
  p <- prcomp(t(top)); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  d <- cbind(s[, .(gse, group, scan_date)], PC1 = p$x[, 1], PC2 = p$x[, 2])
  ggplot(d, aes(PC1, PC2, colour = group, shape = gse)) + geom_point(size = 3) +
    labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"))
}
save_fig(pca_plot(gpl570_expr, gpl570_samples, "GPL570 unique arrays, RMA (before batch handling)"), "gpl570_pca_raw")
save_fig(pca_plot(pooled, gpl570_samples, "GPL570 pooled matrix used for WGCNA"), "gpl570_pca_pooled")

save_rds(list(expr = gpl570_expr, pooled = pooled, samples = gpl570_samples, per_series = per_series,
              batch_note = batch_note, reannotated = reannotated, annotation_source = annotation_source,
              sex_separation = attr(gpl570_sc, "sex_separation")), "gpl570")

# ================================================================= GSE23754 (two-colour cDNA)
es <- readRDS(file.path(DATA, "geo", paste0(cfg$geo$cdna, "_eset.rds")))
cs <- load_rds("gse23754_sample_sheet")
ratio <- exprs(es)[, cs$gsm]
fd <- fData(es)
parts <- strsplit(as.character(fd$NAME), "\\|\\|")
old_symbol <- trimws(vapply(parts, function(p) if (length(p) >= 3) p[3] else "", ""))
entrez <- trimws(vapply(parts, function(p) if (length(p) >= 7) p[7] else "", ""))
# current symbol via Entrez ID (survives renames); fall back to HGNC alias resolution
sym_from_eg <- suppressMessages(mapIds(org.Hs.eg.db, ifelse(entrez == "", "0", entrez), "SYMBOL", "ENTREZID",
                                       multiVals = "first"))
sym_from_alias <- tryCatch(suppressWarnings(alias2SymbolTable(ifelse(old_symbol == "", "NA_none", old_symbol), "Hs")),
                           error = function(e) rep(NA_character_, length(old_symbol)))
cur_symbol <- ifelse(!is.na(sym_from_eg), sym_from_eg, sym_from_alias)
cur_gene <- gene_lookup$gene_id[match(cur_symbol, gene_lookup$symbol)]
cdna_ann <- data.table(spot = rownames(ratio), old_symbol, entrez, current_symbol = cur_symbol, gene_id = cur_gene)
cdna_ann[, status := fcase(old_symbol == "" & is.na(gene_id), "unannotated clone",
                           old_symbol == "" & !is.na(gene_id), "recovered",
                           is.na(gene_id), "lost / withdrawn",
                           old_symbol == current_symbol, "unchanged",
                           default = "renamed")]
write_tab(cdna_ann[, .N, by = status][order(-N)], "gse23754_reannotation_summary")

# Average spots per gene (ratios are gene-centred, so "highest mean" probe selection is meaningless).
keep <- !is.na(cur_gene)
gene_mean <- rowsum(ifelse(is.na(ratio[keep, ]), 0, ratio[keep, ]), cur_gene[keep])
gene_n <- rowsum(1 * !is.na(ratio[keep, ]), cur_gene[keep])
cdna_expr <- gene_mean / gene_n
cdna_expr[gene_n == 0] <- NA
cdna_expr <- cdna_expr[rowMeans(!is.na(cdna_expr)) >= 0.7, ]
msg("GSE23754: ", nrow(cdna_expr), " genes with >=70% observed values")

cdna_sc <- confounder_scores(to_symbol(cdna_expr))
cdna_samples <- cbind(cs, cdna_sc)
save_rds(list(expr = cdna_expr, samples = cdna_samples,
              sex_separation = attr(cdna_sc, "sex_separation")), "gse23754")

# ================================================================= confounder-vs-group check
conf_names <- c(setdiff(names(gene_sets), c("sex_female", "sex_male")), "sex_score")
assoc <- function(samples, cohort) rbindlist(lapply(conf_names, function(v) {
  x <- samples[[v]]; g <- samples$group
  if (all(is.na(x)) || uniqueN(g[!is.na(x)]) < 2) return(NULL)
  med <- as.list(tapply(x, g, median, na.rm = TRUE))
  names(med) <- paste0("median_", names(med))
  cbind(data.table(cohort = cohort, score = v, kruskal_p = kruskal.test(x ~ factor(g))$p.value),
        as.data.table(med))
}), fill = TRUE)
conf_assoc <- rbind(assoc(per_series$GSE46922$samples, "GSE46922"),
                    assoc(per_series$GSE56232$samples, "GSE56232"),
                    assoc(cdna_samples, "GSE23754"), fill = TRUE)
write_tab(conf_assoc, "confounder_group_association")

long <- rbind(
  melt(per_series$GSE46922$samples[, c("group", conf_names), with = FALSE], id.vars = "group")[, cohort := "GSE46922"],
  melt(per_series$GSE56232$samples[, c("group", conf_names), with = FALSE], id.vars = "group")[, cohort := "GSE56232"],
  melt(cdna_samples[, c("group", conf_names), with = FALSE], id.vars = "group")[, cohort := "GSE23754"])
p <- ggplot(long[is.finite(value)], aes(group, value, fill = group)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8) + facet_grid(variable ~ cohort, scales = "free") +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = "score (mean z)", title = "Expression-derived confounder scores by disease group")
save_fig(p, "confounder_scores_by_group", width = 9, height = 14)
msg("Expression preprocessing complete")
