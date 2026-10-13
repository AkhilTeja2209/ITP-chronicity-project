# Step 21 - Benchmark: does a transcriptomic signature beat the baseline?
# Evaluation cohort: GSE23754 (whole blood, cDNA platform). Signatures are *derived only from the
# GPL570 T-cell cohorts*, so evaluation here is a true leave-one-cohort-out test with no leakage.
#   T1 (primary, predictive): progression-to-chronic vs self-limited acute; both sampled <6 months
#   T2 (secondary): eventually chronic (progression + chronic) vs acute; chronic samples were drawn
#      later, so T2 is confounded by disease duration and is reported for context only.
# Baseline: GEO provides no age, sex or platelet count for GSE23754, so the baseline uses their
# expression-derived surrogates (platelet transcript score, inferred sex). This is stated in the report.

source("R/utils.R")
suppressPackageStartupMessages({ library(pROC); library(glmnet) })

cdna <- load_rds("gse23754")
de <- load_rds("de_results")[model == "adjusted"]
integ <- if (has_rds("integration")) load_rds("integration")$evidence else NULL
wg <- if (has_rds("wgcna")) load_rds("wgcna") else NULL
X <- cdna$expr; S <- cdna$samples
N <- cfg$integration$top_n_signature
nperm <- cfg$benchmark$n_permutations

sig_score <- function(up, down) {
  u <- score_set(X, up, 3); d <- score_set(X, down, 3)
  if (all(is.na(u)) && all(is.na(d))) return(rep(NA_real_, ncol(X)))
  ifelse(is.na(u), 0, u) - ifelse(is.na(d), 0, d)
}
auc_of <- function(y, s) { ok <- is.finite(s); if (length(unique(y[ok])) < 2) return(NA_real_)
  as.numeric(auc(roc(y[ok], s[ok], levels = c(0, 1), direction = "<", quiet = TRUE))) }

# ---- fixed (unfitted) scores
t2 <- de[cohort == "GSE46922" & contrast == "chronic_vs_ND" & gene_id %in% rownames(X)][order(P.Value)]
tc_up <- head(t2[logFC > 0, gene_id], N); tc_down <- head(t2[logFC < 0, gene_id], N)
scores <- list(
  platelet_score_only = -S$platelet,   # lower platelet transcripts = lower count = worse prognosis
  Tcell_signature_GSE46922 = sig_score(tc_up, tc_down))
sig_sizes <- list(Tcell_signature_GSE46922 = c(length(tc_up), length(tc_down)))

if (!is.null(integ)) {
  # F2+F3 or F2+F4 genes, direction from F2; F1 is excluded from selection (it is GSE23754 itself)
  ig <- integ[F2 == TRUE & (F3 == TRUE | F4 == TRUE) & gene_id %in% rownames(X)]
  if (nrow(ig) >= 3) {
    scores$integrated_signature <- sig_score(ig[logFC_F2 > 0, gene_id], ig[logFC_F2 < 0, gene_id])
    sig_sizes$integrated_signature <- c(sum(ig$logFC_F2 > 0), sum(ig$logFC_F2 < 0))
  }
}
if (!is.null(wg)) {
  mt <- wg$module_tests[order(GSE46922_p)][1]   # chosen on GPL570 only
  g <- wg$membership[module == mt$module, gene_id]
  scores[[paste0("module_", mt$module)]] <- sign(mt$GSE46922_chronic_minus_ND) * score_set(X, g, 10)
}

tasks <- list(
  T1_progression_vs_acute = list(idx = S$group %in% c("progression", "acute"), y = as.integer(S$group == "progression")),
  T2_eventual_chronic_vs_acute = list(idx = rep(TRUE, nrow(S)), y = S$eventual_chronic))

loo_glm <- function(df) {
  vapply(seq_len(nrow(df)), function(i) {
    f <- suppressWarnings(glm(y ~ ., data = df[-i, , drop = FALSE], family = binomial))
    suppressWarnings(predict(f, df[i, , drop = FALSE], type = "response"))
  }, 0)
}
loo_enet <- function(x, y) {      # nested: feature preselection + lambda choice inside each fold
  vapply(seq_along(y), function(i) {
    xt <- x[-i, , drop = FALSE]; yt <- y[-i]
    tt <- abs(apply(xt, 2, function(v) tryCatch(t.test(v[yt == 1], v[yt == 0])$statistic, error = function(e) 0)))
    keep <- order(-tt)[1:min(200, ncol(xt))]
    cv <- tryCatch(suppressWarnings(cv.glmnet(xt[, keep], yt, family = "binomial", alpha = 0.5,
                                              nfolds = min(5, length(yt)), grouped = FALSE)), error = function(e) NULL)
    if (is.null(cv)) return(mean(yt))
    as.numeric(predict(cv, x[i, keep, drop = FALSE], s = "lambda.min", type = "response"))
  }, 0)
}

bench <- list(); roc_pts <- list()
for (tn in names(tasks)) {
  tk <- tasks[[tn]]; y <- tk$y[tk$idx]; s <- S[tk$idx]
  preds <- lapply(scores, function(v) v[tk$idx])
  base_df <- data.frame(y = y, platelet = s$platelet)
  if (!all(is.na(s$sex)) && uniqueN(na.omit(s$sex)) > 1) base_df$sex <- factor(s$sex)
  base_df <- base_df[, colSums(is.na(base_df)) == 0, drop = FALSE]
  preds$baseline_logistic_LOOCV <- loo_glm(base_df)
  if (!is.null(preds$Tcell_signature_GSE46922)) {
    preds$baseline_plus_Tcell_signature_LOOCV <- loo_glm(cbind(base_df, sig = preds$Tcell_signature_GSE46922))
  }
  xm <- t(X[, tk$idx]); xm <- xm[, colSums(is.na(xm)) == 0, drop = FALSE]
  set.seed(cfg$seed)
  preds$within_cohort_elastic_net_nestedLOOCV <- loo_enet(xm, y)

  for (mn in names(preds)) {
    p <- preds[[mn]]; a <- auc_of(y, p)
    ok <- is.finite(p)
    ci <- tryCatch(as.numeric(ci.auc(roc(y[ok], p[ok], levels = c(0, 1), direction = "<", quiet = TRUE), method = "delong")),
                   error = function(e) c(NA, NA, NA))
    perm <- if (is.finite(a)) mean(replicate(nperm, auc_of(sample(y), p)) >= a) else NA_real_
    # random gene-set null for unfitted signatures: same sizes, random genes from the platform
    rnd <- NA_real_
    if (mn %in% names(sig_sizes) && is.finite(a)) {
      sz <- sig_sizes[[mn]]; pool <- rownames(X)
      null <- replicate(500, { g <- sample(pool, sum(sz)); auc_of(y, sig_score(g[seq_len(sz[1])], g[-seq_len(sz[1])])[tk$idx]) })
      rnd <- mean(null >= a, na.rm = TRUE)
    }
    bench[[length(bench) + 1]] <- data.table(task = tn, model = mn, n = length(y), n_cases = sum(y),
      AUC = a, CI_low = ci[1], CI_high = ci[3], permutation_p = perm, random_geneset_p = rnd)
    if (is.finite(a)) {
      r <- roc(y[ok], p[ok], levels = c(0, 1), direction = "<", quiet = TRUE)
      roc_pts[[length(roc_pts) + 1]] <- data.table(task = tn, model = mn, fpr = 1 - r$specificities, tpr = r$sensitivities)
    }
  }
}
bench <- rbindlist(bench)
bench[, beats_chance := is.finite(CI_low) & CI_low > 0.5]
write_tab(bench, "benchmark_auc")
print(bench)

rp <- rbindlist(roc_pts)
if (nrow(rp)) {
  p <- ggplot(rp[order(fpr, tpr)], aes(fpr, tpr, colour = model)) + geom_step(linewidth = 0.8) +
    geom_abline(linetype = 2, colour = "grey50") + facet_wrap(~ task) + coord_equal() +
    labs(x = "1 - specificity", y = "sensitivity", title = "GSE23754 evaluation (signatures trained on GPL570 only)") +
    theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(p, "benchmark_roc", width = 10, height = 6.5)
}
save_rds(bench, "benchmark")
