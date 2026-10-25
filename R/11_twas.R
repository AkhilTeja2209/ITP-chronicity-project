# Step 11 - Transcriptome-wide association (S-PrediXcan formulation).
# PredictDB GTEx v8 MASHR whole-blood models give, per gene, SNP weights that predict expression.
# Combined with GWAS z-scores and the model's SNP covariance:
#     z_gene = sum_l  w_l * sigma_l * z_l / sigma_gene,   sigma_gene^2 = w' Sigma w
# (Barbeira et al. 2018 Nat Commun). Implemented directly to avoid the MetaXcan Python 2-era stack.

source("R/utils.R")
source("R/twas_lib.R")

gw <- load_rds("gwas")
gene_lookup <- load_rds("gene_lookup")
models <- load_twas_models()
weights <- models$weights; extra <- models$extra; twas_dir <- models$dir
msg("TWAS models: ", uniqueN(weights$gene), " genes, ", uniqueN(weights$varID), " SNPs")

# Extract only model SNP positions from the full GWAS file
keys <- file.path(twas_dir, "model_positions.txt")
fwrite(unique(weights[, .(chr, pos)]), keys, sep = "\t", col.names = FALSE)
hdr <- gw$header
ci <- match(c("chromosome", "base_pair_location"), hdr)
cmd <- sprintf("gzip -dc '%s' | awk -F'\\t' '{sub(/\\r$/, \"\")} NR==FNR{k[$1\":\"$2]=1; next} FNR==1 || (($%d\":\"$%d) in k)' '%s' -",
               gw$file, ci[1], ci[2], keys)
g <- fread(cmd = cmd, sep = "\t", header = TRUE, na.strings = c("NA", ""))
g <- g[is.finite(beta) & is.finite(standard_error) & standard_error > 0]
g[, `:=`(chr = as.character(chromosome), pos = base_pair_location, z = beta / standard_error)]

res <- twas_gene_z(models, g)
n_model <- weights[, .N, by = gene]
res[, gene_id := sub("\\..*", "", model_gene)]
res[n_model, n_snps_model := i.N, on = c(model_gene = "gene")]
res[, p := 2 * pnorm(-abs(z))]
res[, fdr := p.adjust(p, "BH")]
res[, bonferroni_sig := p < 0.05 / .N]
res[, symbol := gene_lookup$symbol[match(gene_id, gene_lookup$gene_id)]]
res[is.na(symbol), symbol := extra$genename[match(model_gene[is.na(symbol)], extra$gene)]]
res[, c("chr", "start") := gene_lookup[match(res$gene_id, gene_id), .(chr, start)]]
setorder(res, p)
write_tab(res, "twas_whole_blood")
msg("TWAS: ", nrow(res), " genes tested; ", sum(res$fdr < cfg$integration$twas_fdr), " at FDR<",
    cfg$integration$twas_fdr, "; ", sum(res$bonferroni_sig), " Bonferroni")

pl <- res[!is.na(chr) & chr %in% c(1:22, "X")]
pl[, chr_f := factor(chr, c(1:22, "X"))]
p <- ggplot(pl, aes(as.integer(chr_f) + start / 3e8, -log10(p), colour = as.integer(chr_f) %% 2 == 0)) +
  geom_point(size = 0.7) + scale_colour_manual(values = c("grey30", "darkorange"), guide = "none") +
  geom_hline(yintercept = -log10(0.05 / nrow(res)), linetype = 2, colour = "firebrick") +
  geom_text(data = head(pl, 12), aes(label = symbol), size = 2.5, vjust = -0.6, colour = "black") +
  scale_x_continuous(breaks = 1:23, labels = c(1:22, "X")) +
  labs(x = "chromosome", y = expression(-log[10](p)), title = "TWAS (S-PrediXcan, GTEx v8 whole blood) for ITP risk") +
  theme(panel.grid = element_blank())
save_fig(p, "twas_manhattan", width = 12, height = 4.5)
save_rds(res, "twas")
