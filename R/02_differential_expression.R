# Step 02 - Differential expression (limma), per cohort, each fitted three ways:
#   unadjusted   : ~ group
#   adjusted     : ~ group + cohort-specific confounder scores (config de_covariates)
#   sensitivity  : adjusted + glucocorticoid score (does a steroid signature explain the hits?)
# Contrasts:
#   GSE23754  progression_vs_acute  (both sampled <6 months from diagnosis: the only *predictive* contrast)
#             chronic_vs_acute      (chronic sampled later: confounded with disease duration)
#   GSE46922  chronic_vs_ND         (sorted T cells, cross-sectional)
#   GSE56232  remission_vs_ND       (sorted T cells, same children 12 months later: paired)

source("R/utils.R")
suppressPackageStartupMessages(library(limma))

gene_lookup <- load_rds("gene_lookup")
gpl <- load_rds("gpl570")
cdna <- load_rds("gse23754")

fit_contrasts <- function(expr, samples, ref, covars, contrasts, block = NULL) {
  samples <- copy(samples)
  samples[, group := relevel(factor(group), ref = ref)]
  covars <- covars[vapply(covars, function(v) {
    x <- samples[[v]]; !all(is.na(x)) && uniqueN(na.omit(x)) > 1 }, TRUE)]
  terms <- c(block, "group", covars)
  ok <- complete.cases(samples[, terms, with = FALSE])
  s <- samples[ok]; m <- expr[, ok, drop = FALSE]
  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), s)
  design <- design[, qr(design)$pivot[seq_len(qr(design)$rank)], drop = FALSE]
  fit <- eBayes(lmFit(m, design))
  rbindlist(lapply(names(contrasts), function(cn) {
    coef <- paste0("group", contrasts[[cn]])
    if (!coef %in% colnames(design)) return(NULL)
    tt <- as.data.table(topTable(fit, coef = coef, number = Inf, sort.by = "none"), keep.rownames = "gene_id")
    tt[, `:=`(contrast = cn, covariates = paste(c(block, covars), collapse = "+"),
              n = nrow(s), df_resid = fit$df.residual[1])]
  }))
}

specs <- list(
  GSE23754 = list(expr = cdna$expr, samples = cdna$samples, ref = "acute", block = NULL,
                  contrasts = list(progression_vs_acute = "progression", chronic_vs_acute = "chronic")),
  GSE46922 = list(expr = gpl$per_series$GSE46922$expr, samples = gpl$per_series$GSE46922$samples,
                  ref = "ND", block = NULL, contrasts = list(chronic_vs_ND = "chronic")),
  GSE56232 = list(expr = gpl$per_series$GSE56232$expr, samples = gpl$per_series$GSE56232$samples,
                  ref = "ND", block = "subject", contrasts = list(remission_vs_ND = "remission")))

de <- rbindlist(lapply(names(specs), function(id) {
  sp <- specs[[id]]
  adj <- unlist(cfg$de_covariates[[id]])
  models <- list(unadjusted = character(), adjusted = adj, sensitivity_glucocorticoid = unique(c(adj, "glucocorticoid")))
  rbindlist(lapply(names(models), function(mn) {
    msg(id, " - ", mn)
    r <- fit_contrasts(sp$expr, sp$samples, sp$ref, models[[mn]], sp$contrasts, sp$block)
    r[, `:=`(cohort = id, model = mn)]
  }))
}))
de[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]
de[, biotype := gene_lookup$biotype[match(gene_id, gene_lookup$gene_id)]]
setcolorder(de, c("cohort", "contrast", "model", "gene_id", "symbol", "biotype", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val"))
save_rds(de, "de_results")
write_tab(de[model == "adjusted"][order(cohort, contrast, P.Value)], "de_adjusted_all")

de_summary <- de[, .(n_genes = .N, n_p05 = sum(P.Value < cfg$de$p_nominal),
                     n_fdr = sum(adj.P.Val < cfg$de$fdr), n = n[1], df_resid = df_resid[1],
                     covariates = covariates[1]),
                 by = .(cohort, contrast, model)]
# Expected false positives at p<0.05 is 5% of tested genes: report the ratio so the reader can tell
# a real signal from the null.
de_summary[, ratio_vs_null := round(n_p05 / (0.05 * n_genes), 2)]
write_tab(de_summary, "de_summary")
print(de_summary)

# ---------------------------------------------------------------- platelet / steroid artefact check
artefact <- rbindlist(lapply(names(specs), function(id) {
  sp <- specs[[id]]
  plt <- sp$samples$platelet; gc <- sp$samples$glucocorticoid
  rbindlist(lapply(unique(de[cohort == id, contrast]), function(cn) {
    rbindlist(lapply(c("unadjusted", "adjusted"), function(mn) {
      top <- de[cohort == id & contrast == cn & model == mn][order(P.Value)][1:100, gene_id]
      top <- intersect(top, rownames(sp$expr))
      cor_with <- function(score) if (all(is.na(score))) NA_real_ else
        mean(abs(apply(sp$expr[top, , drop = FALSE], 1, cor, y = score, use = "pairwise")) > 0.6, na.rm = TRUE)
      data.table(cohort = id, contrast = cn, model = mn,
                 top100_frac_platelet_correlated = cor_with(plt),
                 top100_frac_steroid_correlated = cor_with(gc),
                 platelet_marker_genes_in_top100 = sum(de[cohort == id & contrast == cn & model == mn &
                   gene_id %in% top, symbol] %in% gene_sets$platelet))
    }))
  }))
}))
write_tab(artefact, "de_artefact_check")

# ---------------------------------------------------------------- volcano plots
main <- de[model == "adjusted"]
main[, label := ifelse(rank(P.Value) <= 10, symbol, NA), by = .(cohort, contrast)]
p <- ggplot(main, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = adj.P.Val < cfg$de$fdr), size = 0.5, alpha = 0.6) +
  geom_text(aes(label = label), size = 2.5, vjust = -0.5, na.rm = TRUE, check_overlap = TRUE) +
  geom_hline(yintercept = -log10(cfg$de$p_nominal), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), name = paste0("FDR<", cfg$de$fdr)) +
  facet_wrap(~ paste(cohort, contrast), scales = "free") +
  labs(title = "Differential expression (confounder-adjusted limma)")
save_fig(p, "de_volcano", width = 10, height = 7)

# p-value histograms: flat = no signal, spike at 0 = signal
p <- ggplot(main, aes(P.Value)) + geom_histogram(bins = 40, boundary = 0, fill = "steelblue") +
  facet_wrap(~ paste(cohort, contrast), scales = "free_y") + labs(title = "p-value distributions")
save_fig(p, "de_pvalue_histograms", width = 10, height = 6)
msg("Differential expression complete")
