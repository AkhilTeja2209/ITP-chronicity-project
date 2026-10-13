# Step 20 - Integrate the variant layer with the expression layer.
# Every gene gets four independent filters:
#   F1  GSE23754  differs early (progression-to-chronic vs self-limited, both <6 months from diagnosis)
#   F2  GSE46922  differs in chronic vs newly diagnosed T cells
#   F3  GSE56232  reverses when the same children remit (opposite sign to F2)
#   F4  genetics  colocalised eQTL (PP.H4 > threshold) or TWAS FDR significant
# and the overlap between genetic and expression evidence is tested formally, so "a few genes overlap"
# is never mistaken for convergence when it is what chance alone would give.

source("R/utils.R")
suppressPackageStartupMessages(library(limma))

gene_lookup <- load_rds("gene_lookup")
gpl <- load_rds("gpl570")
de <- load_rds("de_results")[model == "adjusted"]
wg <- if (has_rds("wgcna")) load_rds("wgcna") else NULL
twas <- if (has_rds("twas")) load_rds("twas") else NULL
coloc <- if (has_rds("coloc")) load_rds("coloc") else NULL
loci <- if (has_rds("gwas")) load_rds("gwas")$loci else NULL
pn <- cfg$de$p_nominal

pick <- function(co, cn, tag) {
  d <- de[cohort == co & contrast == cn, .(gene_id, logFC, t, P.Value)]
  setnames(d, c("logFC", "t", "P.Value"), paste0(c("logFC_", "t_", "p_"), tag)); d
}
ev <- data.table(gene_id = rownames(gpl$expr))
ev <- merge(ev, pick("GSE23754", "progression_vs_acute", "F1"), all.x = TRUE)
ev <- merge(ev, pick("GSE46922", "chronic_vs_ND", "F2"), all.x = TRUE)
ev <- merge(ev, pick("GSE56232", "remission_vs_ND", "F3"), all.x = TRUE)
ev[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]

# F1 counts only when the early whole-blood change points the same way as the chronic T-cell change;
# a gene up in chronic T cells but down in early progressors is not converging evidence.
ev[, F1 := !is.na(p_F1) & p_F1 < pn & (is.na(logFC_F2) | sign(logFC_F1) == sign(logFC_F2))]
ev[, F2 := !is.na(p_F2) & p_F2 < pn]
ev[, F3 := !is.na(p_F3) & p_F3 < pn & sign(logFC_F3) == -sign(logFC_F2)]
ev[, F1_same_direction_as_F2 := sign(logFC_F1) == sign(logFC_F2)]

if (!is.null(wg)) {
  ev[wg$membership, module := i.module, on = "gene_id"]
  ev[wg$module_tests, module_concordant := i.concordant, on = "module"]
}
if (!is.null(twas)) ev[twas, `:=`(twas_z = i.z, twas_p = i.p, twas_fdr = i.fdr), on = "gene_id"]
if (!is.null(coloc) && nrow(coloc)) {
  cbest <- coloc[credible == TRUE][order(-PP.H4)][!duplicated(gene_id), .(gene_id, coloc_PP4 = PP.H4, coloc_dataset = eqtl_dataset,
                                                         coloc_locus = locus, risk_allele_raises_expression)]
  ev <- merge(ev, cbest, by = "gene_id", all.x = TRUE)
}
for (v in c("twas_fdr", "coloc_PP4")) if (!v %in% names(ev)) ev[, (v) := NA_real_]
if (!is.null(loci) && nrow(loci)) {
  gl <- gene_lookup[gene_id %in% ev$gene_id]
  in_locus <- unique(unlist(lapply(seq_len(nrow(loci)), function(i)
    gl[chr == loci$chr[i] & end >= loci$start[i] & start <= loci$end[i], gene_id])))
  ev[, in_gwas_locus := gene_id %in% in_locus]
}
ev[, F4 := (!is.na(coloc_PP4) & coloc_PP4 > cfg$coloc$pp4) | (!is.na(twas_fdr) & twas_fdr < cfg$integration$twas_fdr)]
ev[, evidence_score := F1 + F2 + F3 + F4]
# genetics-expression direction: TWAS z>0 = higher predicted expression raises ITP risk
if ("twas_z" %in% names(ev)) ev[, genetic_expression_direction_agree := sign(twas_z) == sign(logFC_F2)]
fisher_p <- function(...) { p <- c(...); p <- p[is.finite(p)]
  if (!length(p)) NA_real_ else pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE) }
ev[, combined_p := mapply(fisher_p, p_F1, p_F2, p_F3, if ("twas_p" %in% names(ev)) twas_p else NA_real_)]
setorder(ev, -evidence_score, combined_p)
write_tab(ev, "integration_gene_evidence")
cand <- ev[evidence_score >= 2]
write_tab(cand, "integration_candidates")
msg(nrow(cand), " genes pass >=2 filters; ", sum(ev$evidence_score >= 3), " pass >=3; ",
    sum(ev$evidence_score == 4), " pass all 4")

# ---------------------------------------------------------------- formal convergence tests
tests <- list()
universe <- if ("twas_p" %in% names(ev)) ev[!is.na(twas_p) & !is.na(p_F2)] else ev[0]
if (nrow(universe)) {
  gen_set <- universe[F4 == TRUE, gene_id]
  gen_label <- "F4 genes (coloc or TWAS FDR)"
  if (length(gen_set) < 10) {   # too few for a test: use top 1% by TWAS p, labelled as such
    gen_set <- universe[order(twas_p)][seq_len(max(10, ceiling(0.01 * .N))), gene_id]
    gen_label <- "top 1% TWAS genes (too few F4 genes)"
  }
  for (f in c("F1", "F2", "F3")) {
    tab <- table(universe$gene_id %in% gen_set, universe[[f]])
    if (all(dim(tab) == 2)) {
      ft <- fisher.test(tab, alternative = "greater")
      tests[[length(tests) + 1]] <- data.table(test = paste("Fisher overlap:", gen_label, "vs", f),
        n_universe = nrow(universe), n_set = length(gen_set), n_overlap = tab["TRUE", "TRUE"],
        expected = round(length(gen_set) * mean(universe[[f]]), 2), odds_ratio = unname(ft$estimate), p = ft$p.value)
    }
  }
  # Competitive pre-ranked test (accounts for inter-gene correlation) of the genetic set in each contrast
  for (f in c("F1", "F2", "F3")) {
    tv <- universe[[paste0("t_", f)]]; ok <- is.finite(tv)
    if (sum(ok) < 100) next
    cr <- cameraPR(setNames(tv[ok], universe$gene_id[ok]), list(genetic = gen_set))
    tests[[length(tests) + 1]] <- data.table(test = paste("cameraPR:", gen_label, "in", f),
      n_universe = sum(ok), n_set = cr$NGenes, n_overlap = NA, expected = NA,
      odds_ratio = NA, p = cr$PValue, direction = cr$Direction)
  }
  # Genome-wide: does genetically predicted expression track the chronic-vs-ND difference?
  for (f in c("F1", "F2", "F3")) {
    tv <- universe[[paste0("t_", f)]]
    ct <- suppressWarnings(cor.test(universe$twas_z, tv, method = "spearman"))
    tests[[length(tests) + 1]] <- data.table(test = paste("Spearman TWAS z vs t", f),
      n_universe = sum(is.finite(tv)), n_set = NA, n_overlap = NA, expected = NA,
      odds_ratio = unname(ct$estimate), p = ct$p.value)
  }
  # Hypothesis-driven: test F2 only inside the genetic candidate set (small multiple-testing burden)
  sub <- universe[gene_id %in% gen_set]
  sub[, fdr_within := p.adjust(p_F2, "BH")]
  tests[[length(tests) + 1]] <- data.table(test = paste("A priori: F2 BH-FDR<0.1 within", gen_label),
    n_universe = nrow(sub), n_set = nrow(sub), n_overlap = sum(sub$fdr_within < 0.1), expected = NA,
    odds_ratio = NA, p = suppressWarnings(min(sub$fdr_within)))
  write_tab(sub[order(p_F2)], "integration_apriori_genetic_set")
}
tests <- rbindlist(tests, fill = TRUE)
write_tab(tests, "integration_tests")
print(tests)

if (nrow(universe)) {
  lab <- universe[order(-(F4 * 2 + F2))][1:min(15, .N)]
  p <- ggplot(universe, aes(twas_z, t_F2)) + geom_point(size = 0.4, colour = "grey60") +
    geom_point(data = universe[F4 == TRUE], colour = "darkorange", size = 1.2) +
    geom_point(data = universe[F4 == TRUE & F2 == TRUE], colour = "firebrick", size = 2) +
    geom_text(data = lab, aes(label = symbol), size = 2.6, vjust = -0.7) +
    labs(x = "TWAS z (ITP risk, whole blood)", y = "t (chronic vs ND, GSE46922 T cells)",
         title = "Inherited risk vs chronicity expression", subtitle = "orange = F4 genetic evidence; red = F4 and F2")
  save_fig(p, "integration_twas_vs_chronic", width = 7, height = 6)
}
fc <- melt(ev[, .(F1, F2, F3, F4)][, lapply(.SD, sum)], measure.vars = 1:4)
p <- ggplot(ev[, .N, by = evidence_score], aes(factor(evidence_score), N)) + geom_col(fill = "steelblue") +
  geom_text(aes(label = N), vjust = -0.3) + scale_y_log10() +
  labs(x = "number of filters passed (F1-F4)", y = "genes (log scale)", title = "Converging evidence per gene")
save_fig(p, "integration_evidence_counts", width = 6, height = 4)
write_tab(fc, "integration_filter_counts")
save_rds(list(evidence = ev, tests = tests), "integration")
