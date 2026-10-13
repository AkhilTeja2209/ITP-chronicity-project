# Step 11 - Transcriptome-wide association (S-PrediXcan formulation).
# PredictDB GTEx v8 MASHR whole-blood models give, per gene, SNP weights that predict expression.
# Combined with GWAS z-scores and the model's SNP covariance:
#     z_gene = sum_l  w_l * sigma_l * z_l / sigma_gene,   sigma_gene^2 = w' Sigma w
# (Barbeira et al. 2018 Nat Commun). Implemented directly to avoid the MetaXcan Python 2-era stack.

source("R/utils.R")
suppressPackageStartupMessages({ library(DBI); library(RSQLite) })

gw <- load_rds("gwas")
gene_lookup <- load_rds("gene_lookup")
twas_dir <- file.path(DATA, "twas")
tarball <- download_once(cfg$twas$url, file.path(twas_dir, "mashr_eqtl.tar"))
members <- untar(tarball, list = TRUE)
want <- grep(paste0("mashr_", cfg$twas$tissue, "\\.(db|txt\\.gz)$"), members, value = TRUE)
if (length(want) != 2) stop("Tissue model files not found in tarball: ", cfg$twas$tissue)
if (!all(file.exists(file.path(twas_dir, want)))) untar(tarball, files = want, exdir = twas_dir)
db_file <- file.path(twas_dir, grep("\\.db$", want, value = TRUE))
cov_file <- file.path(twas_dir, grep("\\.txt\\.gz$", want, value = TRUE))

con <- dbConnect(SQLite(), db_file)
weights <- as.data.table(dbGetQuery(con, "SELECT gene, varID, ref_allele, eff_allele, weight FROM weights"))
extra <- as.data.table(dbGetQuery(con, "SELECT * FROM extra"))
dbDisconnect(con)
weights[, c("chr", "pos") := tstrsplit(sub("^chr", "", varID), "_", keep = 1:2)]
weights[, pos := as.integer(pos)]
msg("TWAS models: ", uniqueN(weights$gene), " genes, ", uniqueN(weights$varID), " SNPs")

# Extract only model SNP positions from the full GWAS file
keys <- file.path(twas_dir, "model_positions.txt")
fwrite(unique(weights[, .(chr, pos)]), keys, sep = "\t", col.names = FALSE)
hdr <- gw$header
ci <- match(c("chromosome", "base_pair_location"), hdr)
cmd <- sprintf("gzip -dc '%s' | awk -F'\\t' 'NR==FNR{k[$1\":\"$2]=1; next} FNR==1 || (($%d\":\"$%d) in k)' '%s' -",
               gw$file, ci[1], ci[2], keys)
g <- fread(cmd = cmd, sep = "\t", header = TRUE, na.strings = c("NA", ""))
g <- g[is.finite(beta) & is.finite(standard_error) & standard_error > 0]
g[, `:=`(chr = as.character(chromosome), pos = base_pair_location, z = beta / standard_error)]

# Allele alignment: model eff_allele must match GWAS effect_allele (flip sign if swapped)
m <- merge(weights, g[, .(chr, pos, effect_allele, other_allele, z)], by = c("chr", "pos"), allow.cartesian = TRUE)
m[, sign := fifelse(effect_allele == eff_allele & other_allele == ref_allele, 1,
                    fifelse(effect_allele == ref_allele & other_allele == eff_allele, -1, NA_real_))]
m <- m[!is.na(sign)][, z := z * sign]
m <- unique(m, by = c("gene", "varID"))

covs <- fread(cov_file)
setnames(covs, c("GENE", "RSID1", "RSID2", "VALUE"))
covs <- covs[GENE %in% unique(m$gene)]
cov_by_gene <- split(covs, by = "GENE", keep.by = FALSE)
m_by_gene <- split(m, by = "gene")
n_model <- weights[, .N, by = gene]

res <- rbindlist(lapply(names(m_by_gene), function(gn) {
  d <- m_by_gene[[gn]]; cv <- cov_by_gene[[gn]]
  if (is.null(cv)) return(NULL)
  ids <- d$varID
  S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  cv <- cv[RSID1 %in% ids & RSID2 %in% ids]
  S[cbind(cv$RSID1, cv$RSID2)] <- cv$VALUE
  S[cbind(cv$RSID2, cv$RSID1)] <- cv$VALUE
  sd_l <- sqrt(pmax(diag(S), 0)); var_g <- drop(t(d$weight) %*% S %*% d$weight)
  if (!is.finite(var_g) || var_g <= 0) return(NULL)
  zg <- sum(d$weight * sd_l * d$z) / sqrt(var_g)
  data.table(model_gene = gn, z = zg, n_snps_used = nrow(d))
}))
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
