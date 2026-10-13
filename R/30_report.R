# Step 30 - Assemble results/REPORT.md from whatever upstream steps succeeded.

source("R/utils.R")

tbl <- function(name, ...) {
  f <- file.path(RES, "tables", paste0(name, ".tsv"))
  if (!file.exists(f)) return("_(not available: upstream step did not run)_")
  md_table(fread(f), ...)
}
fig <- function(name, caption = name) {
  if (file.exists(file.path(RES, "figures", paste0(name, ".png"))))
    sprintf("![%s](figures/%s.png)", caption, name) else ""
}
opt <- function(name) if (has_rds(name)) load_rds(name) else NULL
status <- if (file.exists(file.path(RES, "pipeline_status.tsv")))
  fread(file.path(RES, "pipeline_status.tsv"), colClasses = "character", na.strings = NULL) else NULL

gpl <- opt("gpl570"); gw <- opt("gwas"); twas <- opt("twas"); coloc <- opt("coloc")
integ <- opt("integration"); bench <- opt("benchmark"); wg <- opt("wgcna")
sha <- Sys.getenv("GITHUB_SHA", unset = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "local"))

L <- c(
  "# Pediatric ITP chronicity: integrated variant x expression analysis - results",
  "",
  sprintf("_Generated %s by the automated pipeline (commit `%s`, %s)._", format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"),
          substr(sha, 1, 7), R.version.string),
  "",
  "> Everything below is produced by code in this repository from public data. Numbers are recomputed",
  "> on every run; interpretation sentences marked **[auto]** are generated from the numbers, not written by hand.",
  "",
  "## 0. Pipeline status", "",
  if (!is.null(status)) md_table(status) else "_(status unavailable)_",
  "",
  "## 1. Data audit", "",
  "### Metadata completeness and variant-callability", "",
  tbl("audit_metadata_completeness"),
  "",
  "No dataset contains sequencing reads, so none is variant-callable with GATK; the variant layer therefore",
  "uses published GWAS summary statistics.", "",
  "### Group sizes", "", tbl("audit_group_sizes"), "",
  "### Byte-identical raw arrays shared between series", "",
  tbl("audit_duplicate_arrays"), "",
  "**[auto]** ", if (file.exists(file.path(RES, "tables", "audit_duplicate_arrays.tsv")))
    sprintf("%d raw CEL files deposited in one series are exact copies of files in another, under different patient numbers. They are counted once in every pooled analysis; papers that merged GSE46922 with GSE56232 without checking would have double-counted these children.",
            nrow(fread(file.path(RES, "tables", "audit_duplicate_arrays.tsv")))) else "",
  "",
  "## 2. Expression layer", "",
  sprintf("### Probe re-annotation (GEO GPL570 table vs exact probe-sequence matches to GENCODE v%d)", cfg$annotation$gencode_release), "",
  tbl("gpl570_reannotation_summary"), "",
  "Biotypes of genes that become measurable only under the current annotation:", "",
  tbl("gpl570_reannotation_new_gene_biotypes", max_rows = 10), "",
  "GSE23754 clone re-annotation:", "", tbl("gse23754_reannotation_summary"), "",
  if (!is.null(gpl)) paste0("**Batch handling:** ", gpl$batch_note) else "", "",
  fig("gpl570_pca_raw"), "",
  "### Confounders estimated from expression", "",
  "Kruskal-Wallis test of each score against disease group. A small p here means that score differs between",
  "groups and any gene correlated with it will look 'differentially expressed' for that reason alone.", "",
  tbl("confounder_group_association", max_rows = 60), "",
  fig("confounder_scores_by_group"), "",
  "### Differential expression", "",
  "`ratio_vs_null` = genes at p<0.05 divided by the number expected by chance (5% of genes tested).",
  "Values near 1 mean the contrast is indistinguishable from noise at this sample size.", "",
  tbl("de_summary", max_rows = 30), "",
  "Artefact check (fraction of the top-100 genes strongly correlated, |r|>0.6, with the platelet or steroid score):", "",
  tbl("de_artefact_check"), "",
  fig("de_pvalue_histograms"), "", fig("de_volcano"), "",
  "### Co-expression modules", "",
  if (!is.null(wg)) sprintf("%d modules (soft power %d). Cross-cohort behaviour of each module:", length(wg$modules), wg$power) else "",
  "", tbl("wgcna_module_cross_cohort"), "",
  fig("wgcna_module_trait"), "",
  "## 3. Variant layer", "",
  if (!is.null(gw)) sprintf("GWAS %s: lambda_GC = %.3f; %d variants at p<5e-8; loci defined as **%s**.",
                            cfg$gwas$id, gw$lambda_gc, gw$n_gws, gw$locus_mode) else "_GWAS step did not run._",
  "", fig("gwas_manhattan"), "", tbl("gwas_loci"), "",
  "### TWAS (whole blood)", "",
  if (!is.null(twas)) sprintf("%d genes tested; %d at FDR<%.2f; %d Bonferroni-significant.", nrow(twas),
                              sum(twas$fdr < cfg$integration$twas_fdr), cfg$integration$twas_fdr, sum(twas$bonferroni_sig)) else "",
  "", tbl("twas_whole_blood", max_rows = 20), "", fig("twas_manhattan"), "",
  "### Colocalisation", "",
  if (!is.null(coloc) && nrow(coloc)) {
    hits <- coloc[PP.H4 > cfg$coloc$pp4]
    paste(sprintf("%d locus-gene-tissue tests; %d with PP.H4 > %.1f (shared causal variant)%s.",
                  nrow(coloc), nrow(hits), cfg$coloc$pp4,
                  if (nrow(hits)) paste0(": ", paste(unique(sprintf("%s (%s)", hits$symbol, hits$eqtl_dataset)), collapse = ", ")) else ""),
          if (nrow(hits) && all(hits$mhc)) "**[auto]** All of them lie in the MHC, where extreme LD and HLA allelic diversity make single-variant colocalisation unreliable; they are reported but should not be read as causal-gene evidence without HLA-specific fine-mapping." else "")
  } else "",
  "", tbl("coloc_results", max_rows = 20), "", fig("coloc_top_loci"), "",
  "## 4. Integration", "",
  "Formal tests of whether genetic evidence and chronicity expression converge more than chance:", "",
  tbl("integration_tests"), "",
  if (!is.null(integ)) {
    t <- integ$tests
    sig <- t[is.finite(p) & p < 0.05]
    bonf <- 0.05 / nrow(t)
    paste("**[auto]**", if (nrow(sig)) sprintf("%d of %d convergence tests reach nominal p<0.05 (%s); %s survive Bonferroni correction for %d tests (p<%.4f). Treat the convergence as %s.",
                                               nrow(sig), nrow(t), paste(sig$test, collapse = "; "),
                                               if (any(t$p < bonf, na.rm = TRUE)) sum(t$p < bonf, na.rm = TRUE) else "none",
                                               nrow(t), bonf, if (any(t$p < bonf, na.rm = TRUE)) "supported" else "not established")
          else sprintf("None of the %d convergence tests reaches p<0.05: at current sample sizes, the overlap between inherited ITP risk genes and chronicity-associated expression is not distinguishable from chance. Candidate genes below are leads, not findings.", nrow(t)))
  } else "", "",
  "Top candidate genes (ranked by number of filters passed, then combined p):", "",
  tbl("integration_candidates", max_rows = 30), "",
  fig("integration_twas_vs_chronic"), "", fig("integration_evidence_counts"), "",
  "## 5. Benchmark (GSE23754, signatures trained on GPL570 only)", "",
  tbl("benchmark_auc", max_rows = 30), "",
  if (!is.null(bench)) {
    t1 <- bench[task == "T1_progression_vs_acute"]
    b <- t1[model == "baseline_logistic_LOOCV", AUC]; best <- t1[which.max(AUC)]
    paste("**[auto]** Primary task T1 has", t1$n[1], "samples (", t1$n_cases[1], "progressors).",
          sprintf("Baseline AUC = %.2f; best model = %s (AUC %.2f, 95%% CI %.2f-%.2f, permutation p = %.3f).",
                  b, best$model, best$AUC, best$CI_low, best$CI_high, best$permutation_p),
          if (any(t1$beats_chance)) "At least one model's CI excludes 0.5." else
            "No model's confidence interval excludes 0.5; with 4 progressors this benchmark cannot distinguish a real signature from chance, which is the expected and reportable result at this n.",
          {
            t2 <- bench[task == "T2_eventual_chronic_vs_acute"]
            w <- t2[beats_chance == TRUE]
            if (nrow(w)) sprintf("In secondary task T2 (confounded by disease duration), only %s beat%s chance (%s); transcriptomic signatures add nothing beyond it.",
                                 paste(w$model, collapse = ", "), if (nrow(w) == 1) "s" else "",
                                 paste(sprintf("AUC %.2f, CI %.2f-%.2f", w$AUC, w$CI_low, w$CI_high), collapse = "; "))
            else "In secondary task T2 no model's CI excludes 0.5."
          })
  } else "", "",
  fig("benchmark_roc"), "",
  "## 6. Limitations stated by the design", "",
  "- The GWAS phenotype is ITP *susceptibility* in mostly adult biobank participants, not pediatric chronicity.",
  "- GSE23754 used a pre-2009 6-month definition of chronic ITP; 'progression' samples are the only ones drawn before outcome.",
  "- GEO provides no per-sample age, sex, platelet count or treatment for these series; sex and platelet burden are inferred from expression.",
  "- GSE46922/GSE56232 are sorted T cells; GSE23754 is whole blood on a different platform: cross-cohort agreement is tested, not assumed.",
  "- Time-to-remission (survival) modelling would be stronger but no public dataset links transcriptomes to follow-up time.",
  "", "See `docs/PROJECT_PLAN.md` for the full rationale.")

writeLines(L, file.path(RES, "REPORT.md"))
msg("Report written to ", file.path(RES, "REPORT.md"))
