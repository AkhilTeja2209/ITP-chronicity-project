# Step 04 - Genes that only the GENCODE re-annotation can measure.
# The chronic-vs-newly-diagnosed T-cell difference itself was reported by Jernås et al. 2013
# (GSE46922). What the original GEO annotation could not report are genes no GPL570 probe set was
# assigned to. This step isolates those genes (excluding mere renames), asks whether any are
# differentially expressed, and checks each hit for robustness:
#   - probe support (how many of the probe set's probes match the gene uniquely)
#   - expression level (low-intensity probe sets on GPL570 are mostly noise)
#   - steroid sensitivity (does the hit survive the glucocorticoid-score model?)
#   - reversal (does a chronic-vs-ND change move the opposite way on remission?)

source("R/utils.R")

gene_lookup <- load_rds("gene_lookup")
de <- load_rds("de_results")
ann <- fread(file.path(RES, "tables", "gpl570_probe_reannotation.tsv"))
pmap <- fread(file.path(DATA, "gencode", "gpl570_probeset_to_gene.tsv"))

# ---------------------------------------------------------------- gene-level novelty
# A gene was measurable before if any probe set mapping to it today was already annotated to it
# (same symbol or a renamed symbol). Otherwise it is newly measurable.
g <- ann[!is.na(current_gene) & current_gene != ""]
g <- merge(g, pmap[, .(probe = probeset, n_probes, top)], by = "probe", all.x = TRUE)
gene_status <- g[, .(previously_measurable = any(status %in% c("unchanged", "same gene, renamed")),
                     n_probesets = .N, best_probe_support = max(top / n_probes, na.rm = TRUE),
                     statuses = paste(sort(unique(status)), collapse = "; ")), by = .(gene_id = current_gene)]
gene_status[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]
gene_status[, biotype := gene_lookup$biotype[match(gene_id, gene_lookup$gene_id)]]
gene_status[, newly_measurable := !previously_measurable]
write_tab(gene_status[order(-newly_measurable, biotype)], "novel_gene_status")
msg(sum(gene_status$newly_measurable), " genes newly measurable (",
    gene_status[newly_measurable == TRUE & biotype == "lncRNA", .N], " lncRNA) out of ", nrow(gene_status))

# ---------------------------------------------------------------- DE hits among new genes
pick <- function(co, ct, mdl) de[cohort == co & contrast == ct & model == mdl,
                                 .(gene_id, logFC, AveExpr, P.Value, adj.P.Val)]
contrasts <- list(chronic_vs_ND = c("GSE46922", "chronic_vs_ND"), remission_vs_ND = c("GSE56232", "remission_vs_ND"))
expr_floor <- quantile(de[cohort == "GSE46922" & model == "adjusted", AveExpr], 0.25, na.rm = TRUE)

hits <- rbindlist(lapply(names(contrasts), function(cn) {
  co <- contrasts[[cn]][1]; ct <- contrasts[[cn]][2]
  a <- pick(co, ct, "adjusted"); st <- pick(co, ct, "sensitivity_glucocorticoid")
  d <- merge(a, st[, .(gene_id, fdr_steroid_model = adj.P.Val, logFC_steroid_model = logFC)], by = "gene_id")
  d <- merge(d, gene_status[, .(gene_id, symbol, biotype, newly_measurable, best_probe_support, n_probesets)], by = "gene_id")
  d[, contrast := cn]
  d
}))
# reversal: chronic change opposite to remission change
rev <- dcast(hits[, .(gene_id, contrast, logFC)], gene_id ~ contrast, value.var = "logFC")
hits[rev, reverses_on_remission := sign(i.chronic_vs_ND) == -sign(i.remission_vs_ND), on = "gene_id"]
hits[, `:=`(sig = adj.P.Val < cfg$de$fdr,
            robust = adj.P.Val < cfg$de$fdr & fdr_steroid_model < cfg$de$fdr & AveExpr > expr_floor &
                     best_probe_support >= 0.5)]

summary <- hits[, .(genes_tested = .N, fdr_hits = sum(sig),
                    new_genes_tested = sum(newly_measurable), new_gene_hits = sum(sig & newly_measurable),
                    new_lncRNA_hits = sum(sig & newly_measurable & biotype == "lncRNA"),
                    robust_new_hits = sum(robust & newly_measurable),
                    robust_new_hits_reversing = sum(robust & newly_measurable & reverses_on_remission %in% TRUE),
                    fisher_p_new_enriched = fisher.test(table(newly_measurable, sig), alternative = "greater")$p.value,
                    odds_ratio = unname(fisher.test(table(newly_measurable, sig))$estimate)),
                by = contrast]
write_tab(summary, "novel_gene_summary")
print(summary)

new_hits <- hits[sig & newly_measurable][order(contrast, P.Value)]
write_tab(new_hits[, .(contrast, symbol, gene_id, biotype, logFC, AveExpr, P.Value, adj.P.Val, fdr_steroid_model,
                       best_probe_support, reverses_on_remission, robust)], "novel_gene_hits")

# biotype breakdown of new hits vs all new genes tested
bt <- hits[newly_measurable == TRUE, .(tested = .N, hits = sum(sig)), by = .(contrast, biotype)][order(contrast, -hits)]
write_tab(bt, "novel_gene_biotypes")

# ---------------------------------------------------------------- figure
pd <- hits[, .(contrast, logFC, mlp = -log10(P.Value),
               class = fifelse(!newly_measurable, "previously measurable",
                               fifelse(robust, "new: robust hit", fifelse(sig, "new: FDR hit only", "new: not significant"))))]
pd[, class := factor(class, c("previously measurable", "new: not significant", "new: FDR hit only", "new: robust hit"))]
lab <- new_hits[robust == TRUE][, head(.SD, 6), by = contrast]
p <- ggplot(pd[order(class)], aes(logFC, mlp, colour = class)) + geom_point(size = 0.6, alpha = 0.7) +
  geom_text(data = lab, aes(logFC, -log10(P.Value), label = symbol), inherit.aes = FALSE, size = 2.6, vjust = -0.6) +
  scale_colour_manual(values = c("grey80", "#B9BEC6", "#C98A2B", "#9E1B32"), name = NULL) +
  facet_wrap(~ contrast, scales = "free") +
  labs(x = "log2 fold change", y = expression(-log[10](p)),
       title = "Genes measurable only after GENCODE re-annotation", subtitle = "robust = FDR < 0.1 in adjusted and steroid models, expressed, >= 50% probe support")
save_fig(p, "novel_genes_volcano", width = 10, height = 5)
save_rds(list(status = gene_status, hits = hits, summary = summary), "novel_genes")
