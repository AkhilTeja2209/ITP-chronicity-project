# Step 03 - Co-expression modules (WGCNA) and pathway enrichment.
# Modules are built on the pooled, de-duplicated GPL570 T-cell matrix, then scored in every cohort
# (including the cDNA whole-blood GSE23754) as mean z of member genes. Modules travel across
# platforms far better than single genes, which matters at these sample sizes.

source("R/utils.R")
suppressPackageStartupMessages({ library(WGCNA); library(fgsea) })
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))

gene_lookup <- load_rds("gene_lookup")
gpl <- load_rds("gpl570")
cdna <- load_rds("gse23754")
de <- load_rds("de_results")

# ---------------------------------------------------------------- network
m <- gpl$pooled
m <- m[order(-apply(m, 1, var))[seq_len(min(cfg$wgcna$n_genes, nrow(m)))], ]
datExpr <- t(m)
n <- nrow(datExpr)
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)), networkType = "signed", verbose = 0)
fit_ok <- sft$fitIndices[sft$fitIndices$SFT.R.sq > 0.8 & sft$fitIndices$slope < 0, "Power"]
# WGCNA FAQ defaults for signed networks when scale-free fit is not reached at small n
power <- if (length(fit_ok)) min(fit_ok) else if (n < 20) 18 else if (n < 30) 16 else 14
msg("WGCNA: n = ", n, " samples, soft power = ", power, if (!length(fit_ok)) " (FAQ default; scale-free fit not reached)")
write_tab(as.data.table(sft$fitIndices), "wgcna_soft_threshold")

net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = cfg$wgcna$min_module_size, mergeCutHeight = cfg$wgcna$merge_cut_height,
                        maxBlockSize = ncol(datExpr) + 1, numericLabels = FALSE, pamRespectsDendro = FALSE,
                        verbose = 0)
colors <- net$colors
modules <- setdiff(unique(colors), "grey")
if (!length(modules)) stop("WGCNA found no modules (all genes unassigned); lower min_module_size in config")
kme <- signedKME(datExpr, net$MEs)
membership <- data.table(gene_id = colnames(datExpr), module = colors)
membership[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]
membership[, kME := mapply(function(g, mod) if (mod == "grey") NA_real_ else kme[g, paste0("kME", mod)],
                           gene_id, module)]
write_tab(membership[order(module, -kME)], "wgcna_module_membership")
with_png("wgcna_dendrogram", plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]],
         "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
         main = "GPL570 T-cell co-expression modules"))

# ---------------------------------------------------------------- module-trait (pooled)
s <- gpl$samples
traits <- data.table(chronic_vs_ND = fifelse(s$group == "chronic", 1, fifelse(s$group == "ND", 0, NA_real_)),
                     remission_vs_ND = fifelse(s$group == "remission", 1, fifelse(s$group == "ND", 0, NA_real_)),
                     s[, intersect(c("platelet", "glucocorticoid", "interferon", "monocyte", "b_cell", "nk_cell",
                                     "sex_score"), names(s)), with = FALSE])
MEs <- orderMEs(net$MEs)
mt <- rbindlist(lapply(names(traits), function(tr) rbindlist(lapply(names(MEs), function(me) {
  ok <- is.finite(traits[[tr]])
  ct <- suppressWarnings(cor.test(MEs[[me]][ok], traits[[tr]][ok]))
  data.table(module = sub("^ME", "", me), trait = tr, r = unname(ct$estimate), p = ct$p.value)
}))))
write_tab(mt, "wgcna_module_trait")
p <- ggplot(mt[module != "grey"], aes(trait, module, fill = r)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f\n(%.1e)", r, p)), size = 2.3) +
  scale_fill_gradient2(low = "steelblue", high = "firebrick", limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + labs(x = NULL, y = NULL, title = "Module-trait correlation (pooled GPL570)")
save_fig(p, "wgcna_module_trait", width = 9, height = 0.35 * length(modules) + 3)

# ---------------------------------------------------------------- module scores in every cohort
mod_genes <- split(membership$gene_id, membership$module)[modules]
score_cohort <- function(expr) as.data.table(lapply(mod_genes, function(g) score_set(expr, g, min_genes = 10)))
w <- function(x, y) if (sum(is.finite(x)) >= 2 && sum(is.finite(y)) >= 2)
  wilcox.test(x, y, exact = FALSE)$p.value else NA_real_

s46 <- gpl$per_series$GSE46922; s56 <- gpl$per_series$GSE56232
sc46 <- score_cohort(s46$expr); sc56 <- score_cohort(s56$expr); sc23 <- score_cohort(cdna$expr)
g23 <- cdna$samples$group
pair56 <- dcast(cbind(s56$samples[, .(patient, group)], sc56) |> melt(id.vars = c("patient", "group")),
                patient + variable ~ group)
mod_tests <- rbindlist(lapply(modules, function(mo) {
  pr <- pair56[variable == mo & is.finite(ND) & is.finite(remission)]
  data.table(
    module = mo, n_genes = length(mod_genes[[mo]]),
    n_genes_on_cdna = sum(mod_genes[[mo]] %in% rownames(cdna$expr)),
    GSE46922_chronic_minus_ND = median(sc46[[mo]][s46$samples$group == "chronic"]) - median(sc46[[mo]][s46$samples$group == "ND"]),
    GSE46922_p = w(sc46[[mo]][s46$samples$group == "chronic"], sc46[[mo]][s46$samples$group == "ND"]),
    GSE56232_remission_minus_ND = median(pr$remission - pr$ND),
    GSE56232_paired_p = if (nrow(pr) >= 3) wilcox.test(pr$remission, pr$ND, paired = TRUE, exact = FALSE)$p.value else NA_real_,
    GSE23754_progression_minus_acute = median(sc23[[mo]][g23 == "progression"], na.rm = TRUE) - median(sc23[[mo]][g23 == "acute"], na.rm = TRUE),
    GSE23754_progression_p = w(sc23[[mo]][g23 == "progression"], sc23[[mo]][g23 == "acute"]),
    GSE23754_chronic_minus_acute = median(sc23[[mo]][g23 == "chronic"], na.rm = TRUE) - median(sc23[[mo]][g23 == "acute"], na.rm = TRUE),
    GSE23754_chronic_p = w(sc23[[mo]][g23 == "chronic"], sc23[[mo]][g23 == "acute"]))
}))
# A module is "concordant" when it moves up (or down) in chronic T cells, moves the opposite way on
# remission, and moves the same way in the early progression-to-chronic whole-blood samples.
mod_tests[, concordant := sign(GSE46922_chronic_minus_ND) == -sign(GSE56232_remission_minus_ND) &
                          sign(GSE46922_chronic_minus_ND) == sign(GSE23754_progression_minus_acute)]
write_tab(mod_tests[order(GSE46922_p)], "wgcna_module_cross_cohort")

# ---------------------------------------------------------------- pathway enrichment
gmt_dir <- file.path(DATA, "msigdb")
pathways <- unlist(lapply(cfg$enrichment$collections, function(f)
  gmtPathways(download_once(file.path(cfg$enrichment$base_url, f), file.path(gmt_dir, f)))), recursive = FALSE)

# Module over-representation (universe = network genes)
universe <- unique(na.omit(membership$symbol))
ora <- rbindlist(lapply(modules, function(mo) {
  r <- fora(pathways, unique(na.omit(membership[module == mo, symbol])), universe,
            minSize = cfg$enrichment$min_size, maxSize = cfg$enrichment$max_size)
  if (nrow(r)) cbind(module = mo, as.data.table(r)[, .(pathway, pval, padj, overlap, size)]) else NULL
}))
write_tab(ora[order(module, pval)], "wgcna_module_enrichment")

# GSEA on each adjusted DE contrast, ranked by moderated t
gsea <- rbindlist(lapply(split(de[model == "adjusted" & !is.na(symbol)], by = c("cohort", "contrast")), function(d) {
  d <- d[is.finite(t)][order(-abs(t))][!duplicated(symbol)]
  r <- suppressWarnings(fgseaMultilevel(pathways, setNames(d$t, d$symbol),
                                        minSize = cfg$enrichment$min_size, maxSize = cfg$enrichment$max_size))
  cbind(cohort = d$cohort[1], contrast = d$contrast[1],
        as.data.table(r)[, .(pathway, NES, pval, padj, size,
                             leading_edge = vapply(leadingEdge, function(x) paste(head(x, 15), collapse = ","), ""))])
}))
write_tab(gsea[order(cohort, contrast, pval)], "gsea_de_contrasts")

save_rds(list(membership = membership, modules = modules, power = power, module_tests = mod_tests,
              module_trait = mt), "wgcna")
msg("Modules and enrichment complete: ", length(modules), " modules")
