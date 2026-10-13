# Step 12 - Colocalisation of ITP GWAS loci with blood / T-cell eQTLs (coloc.abf).
# TWAS association can arise when a risk variant and an eQTL variant are merely in LD.
# coloc estimates the posterior that both signals share ONE causal variant (H4) versus two
# distinct variants (H3). eQTL summary statistics are fetched per region from the eQTL Catalogue
# FTP with tabix, so the multi-GB files are never downloaded whole.

source("R/utils.R")
suppressPackageStartupMessages(library(coloc))

gw <- load_rds("gwas")
gene_lookup <- load_rds("gene_lookup")
loci <- gw$loci
eqtl_dir <- normalizePath(file.path(DATA, "eqtl"), mustWork = FALSE)
dir.create(eqtl_dir, recursive = TRUE, showWarnings = FALSE)
helper <- normalizePath("R/tabix_query.py")
python <- Sys.getenv("PYTHON", "python3")
chk <- suppressWarnings(system2(python, c("-c", shQuote("import sys; print(sys.path); import pysam; print(pysam.__version__)")),
                                stdout = TRUE, stderr = TRUE))
chk <- rev(chk)
if (!is.null(attr(chk, "status"))) stop("Python with pysam is required (pip install pysam); ", python, " said: ", paste(chk, collapse = " "))
msg("Using ", python, " with pysam ", chk[1])

# Neither file marks its header with '#', so `tabix -h` does not return it: read the column names
# from the first line once and pass them explicitly. CRs are stripped (the GWAS file is CRLF).
header_cache <- list()
file_header <- function(src) {
  if (is.null(header_cache[[src]])) {
    cmd <- if (grepl("^https?://", src)) sprintf("curl -s '%s' | gzip -dc 2>/dev/null | head -1", src)
           else sprintf("gzip -dc '%s' | head -1", src)
    header_cache[[src]] <<- strsplit(sub("\r$", "", system(cmd, intern = TRUE)[1]), "\t")[[1]]
  }
  header_cache[[src]]
}
tabix <- function(src, region, cwd = getwd()) {
  err <- tempfile()
  out <- withr::with_dir(cwd, suppressWarnings(system2(python, c(shQuote(helper), shQuote(src), region),
                                                       stdout = TRUE, stderr = err)))
  if (length(readLines(err, warn = FALSE))) msg("  tabix stderr: ", paste(head(readLines(err, warn = FALSE), 3), collapse = " | "))
  out <- sub("\r$", "", out)
  if (!length(out)) return(NULL)
  fread(text = out, sep = "\t", header = FALSE, col.names = file_header(src), na.strings = c("NA", ""))
}

s_cases <- cfg$gwas$n_cases / cfg$gwas$n_total
results <- list(); plots <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]
  region <- sprintf("%s:%d-%d", L$chr, as.integer(L$start), as.integer(L$end))
  gr <- tabix(gw$file, region)
  if (is.null(gr)) { msg(L$locus, ": no GWAS rows returned by tabix"); next }
  gr <- gr[is.finite(beta) & is.finite(standard_error) & standard_error > 0]
  gr <- gr[, .(pos = base_pair_location, ea = effect_allele, oa = other_allele, beta, se = standard_error, gp = p_value)]
  for (ds in cfg$eqtl$datasets) {
    url <- sprintf("%s/%s/%s/%s.all.tsv.gz", cfg$eqtl$base_url, ds$study, ds$dataset, ds$dataset)
    msg(L$locus, " ", region, " x ", ds$label)
    eq <- tryCatch(tabix(url, region, eqtl_dir), error = function(e) { msg("  eQTL query error: ", conditionMessage(e)); NULL })
    if (is.null(eq) || !nrow(eq)) { msg("  no eQTL rows returned"); next }
    eq <- eq[is.finite(beta) & is.finite(se) & se > 0 & type == "SNP"]
    genes <- eq[, .(minp = min(pvalue)), by = gene_id][minp < cfg$eqtl$min_eqtl_p, gene_id]
    msg("  GWAS variants ", nrow(gr), "; eQTL rows ", nrow(eq), "; genes with eQTL p<", cfg$eqtl$min_eqtl_p, ": ", length(genes))
    for (gid in genes) {
      e <- eq[gene_id == gid, .(pos = position, ref, alt, eb = beta, ese = se, maf, an)]
      e <- unique(e, by = c("pos", "ref", "alt"))
      # align: eQTL effect is per alt allele; GWAS per effect_allele
      j <- merge(e, gr, by = "pos", allow.cartesian = TRUE)
      j[, sgn := fifelse(ea == alt & oa == ref, 1, fifelse(ea == ref & oa == alt, -1, NA_real_))]
      j <- j[!is.na(sgn) & maf > 0 & maf < 1]
      j[, snp := paste(pos, ref, alt, sep = "_")]
      j <- unique(j, by = "snp")
      if (nrow(j) < 50) { msg("    ", gid, ": only ", nrow(j), " shared SNPs, skipped"); next }
      d1 <- list(beta = j$beta * j$sgn, varbeta = j$se^2, type = "cc", s = s_cases, N = cfg$gwas$n_total, snp = j$snp)
      d2 <- list(beta = j$eb, varbeta = j$ese^2, type = "quant", N = round(median(j$an, na.rm = TRUE) / 2),
                 MAF = j$maf, snp = j$snp)
      cr <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = cfg$coloc$p1, p2 = cfg$coloc$p2, p12 = cfg$coloc$p12)),
                     error = function(e) NULL)
      if (is.null(cr)) next
      top <- as.data.table(cr$results)[order(-SNP.PP.H4)][1]
      results[[length(results) + 1]] <- data.table(
        locus = L$locus, chr = L$chr, lead_rsid = L$lead_rsid, lead_p = L$p, mhc = L$mhc,
        eqtl_dataset = ds$label, gene_id = gid, nsnps = nrow(j),
        PP.H0 = cr$summary[["PP.H0.abf"]], PP.H1 = cr$summary[["PP.H1.abf"]], PP.H2 = cr$summary[["PP.H2.abf"]],
        PP.H3 = cr$summary[["PP.H3.abf"]], PP.H4 = cr$summary[["PP.H4.abf"]],
        top_snp = top$snp, top_snp_PP4 = top$SNP.PP.H4,
        # direction: does the risk-increasing allele raise expression?
        risk_allele_raises_expression = {
          k <- j[snp == top$snp]; sign(k$beta * k$sgn) == sign(k$eb) })
      plots[[paste(L$locus, ds$label, gid)]] <- j[, .(locus = L$locus, eqtl = ds$label, gene_id = gid,
                                                      pos, gwas = -log10(gp), eqtl_p = -log10(2 * pnorm(-abs(eb / ese))))]
    }
  }
}

if (!length(results)) {
  res <- data.table(locus = character(), gene_id = character(), eqtl_dataset = character(), PP.H4 = numeric())
  msg("No locus-gene pairs could be tested for colocalisation")
} else {
  res <- rbindlist(results)
  res[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]
  setorder(res, -PP.H4)
  msg(nrow(res), " locus-gene-tissue tests; ", sum(res$PP.H4 > cfg$coloc$pp4), " with PP.H4 > ", cfg$coloc$pp4)
  top <- head(res[PP.H4 > 0.2], 6)
  if (nrow(top)) {
    pd <- rbindlist(lapply(seq_len(nrow(top)), function(k)
      plots[[paste(top$locus[k], top$eqtl_dataset[k], top$gene_id[k])]][, panel :=
        sprintf("%s %s (%s) PP4=%.2f", top$locus[k], top$symbol[k], top$eqtl_dataset[k], top$PP.H4[k])]))
    p <- ggplot(melt(pd, measure.vars = c("gwas", "eqtl_p")), aes(pos / 1e6, value, colour = variable)) +
      geom_point(size = 0.5, alpha = 0.6) + facet_wrap(~ panel, scales = "free") +
      scale_colour_manual(values = c(gwas = "firebrick", eqtl_p = "steelblue"), labels = c("ITP GWAS", "eQTL")) +
      labs(x = "position (Mb, GRCh38)", y = expression(-log[10](p)), colour = NULL, title = "Top colocalisation signals")
    save_fig(p, "coloc_top_loci", width = 11, height = 6)
  }
}
write_tab(res, "coloc_results")
save_rds(res, "coloc")
