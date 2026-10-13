# Step 10 - ITP GWAS summary statistics: QC, genomic inflation, independent loci.
# Input is already-called variants (summary statistics), so no read processing / GATK is involved.
# The file is ~650 MB with ~25M variants; it is streamed through awk instead of loaded whole.

source("R/utils.R")

gwas_dir <- file.path(DATA, "gwas")
gwas_file <- download_once(cfg$gwas$url, file.path(gwas_dir, basename(cfg$gwas$url)))
download_once(paste0(cfg$gwas$url, ".tbi"), paste0(gwas_file, ".tbi"))

con <- gzfile(gwas_file); hdr <- strsplit(readLines(con, n = 1), "\t")[[1]]; close(con)
col <- function(name) { i <- match(name, hdr); if (is.na(i)) stop("GWAS column missing: ", name); i }
pcol <- col("p_value")

stream <- function(awk_body) {
  cmd <- sprintf("gzip -dc '%s' | awk -F'\\t' -v OFS='\\t' 'BEGIN{srand(%d)} %s'", gwas_file, cfg$seed, awk_body)
  fread(cmd = cmd, sep = "\t", header = TRUE, na.strings = c("NA", ""))
}

msg("Streaming GWAS: significant + thinned background variants")
# keep p<0.01 for plotting and a random 1% of all variants (flag=1) for unbiased inflation estimate
thin <- stream(sprintf("NR==1{print $0,\"random\"; next} {r=(rand()<0.01); if ($%d<0.01 || r) print $0,r}", pcol))
thin <- thin[is.finite(p_value) & is.finite(beta) & is.finite(standard_error)]
thin[, maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]
thin <- thin[is.na(maf) | maf >= cfg$gwas$min_maf]

rnd <- thin[random == 1]
lambda_gc <- median(qchisq(rnd$p_value, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
msg("Genomic inflation lambda_GC = ", round(lambda_gc, 3), " (", nrow(rnd), " random variants)")

# ---------------------------------------------------------------- independent loci (distance clumping)
clump <- function(d, window) {
  d <- d[order(p_value)]; leads <- d[0]
  while (nrow(d)) {
    top <- d[1]; leads <- rbind(leads, top)
    d <- d[!(chromosome == top$chromosome & abs(base_pair_location - top$base_pair_location) < window)]
  }
  leads
}
w <- cfg$gwas$window_kb * 1000
sig <- thin[p_value < cfg$gwas$lead_p]
leads <- clump(sig, w)
locus_mode <- sprintf("independent leads with p < %g", cfg$gwas$lead_p)
if (nrow(leads) == 0) {
  leads <- head(clump(thin[p_value < cfg$gwas$fallback_p], w), cfg$gwas$fallback_top_loci)
  locus_mode <- sprintf("FALLBACK: no variant reached p < %g; top %d suggestive loci (p < %g) used, hypothesis-generating only",
                        cfg$gwas$lead_p, nrow(leads), cfg$gwas$fallback_p)
}
msg(nrow(leads), " loci: ", locus_mode)
gene_lookup <- load_rds("gene_lookup")
loci <- leads[, .(locus = paste0("L", seq_len(.N)), chr = as.character(chromosome), lead_pos = base_pair_location,
                  lead_rsid = rsid, effect_allele, other_allele, beta, se = standard_error, p = p_value, maf,
                  start = pmax(1, base_pair_location - w), end = base_pair_location + w)]
loci[, mhc := chr == "6" & end > 25e6 & start < 34e6]
loci[, genes_in_window := vapply(seq_len(.N), function(i) {
  g <- gene_lookup[chr == loci$chr[i] & end >= loci$start[i] & start <= loci$end[i] &
                     biotype %in% c("protein_coding", "lncRNA")]
  paste(head(g[order(abs((start + end) / 2 - loci$lead_pos[i])), symbol], 15), collapse = ",")
}, "")]
write_tab(loci, "gwas_loci")

# ---------------------------------------------------------------- plots
chr_levels <- c(as.character(1:22), "X")
thin[, chr_f := factor(as.character(chromosome), chr_levels)]
thin <- thin[!is.na(chr_f)]
offs <- thin[, .(len = max(base_pair_location)), by = chr_f][order(chr_f)][, off := cumsum(as.numeric(len)) - len]
thin[offs, cum := base_pair_location + i.off, on = "chr_f"]
axis_pos <- thin[, .(mid = (min(cum) + max(cum)) / 2), by = chr_f]
p <- ggplot(thin[p_value < 0.01 | random == 1], aes(cum, -log10(p_value), colour = as.integer(chr_f) %% 2 == 0)) +
  geom_point(size = 0.4) + scale_colour_manual(values = c("grey30", "steelblue"), guide = "none") +
  geom_hline(yintercept = -log10(5e-8), colour = "firebrick", linetype = 2) +
  geom_hline(yintercept = -log10(cfg$gwas$lead_p), colour = "orange", linetype = 3) +
  scale_x_continuous(breaks = axis_pos$mid, labels = axis_pos$chr_f) +
  labs(x = "chromosome", y = expression(-log[10](p)),
       title = sprintf("ITP GWAS %s (%d cases / %d total); lambda_GC = %.3f", cfg$gwas$id,
                       cfg$gwas$n_cases, cfg$gwas$n_total, lambda_gc)) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6))
save_fig(p, "gwas_manhattan", width = 12, height = 4.5)

qq <- rnd[order(p_value)][, .(obs = -log10(p_value), exp = -log10(ppoints(.N)))]
p <- ggplot(qq, aes(exp, obs)) + geom_point(size = 0.5) + geom_abline(colour = "firebrick") +
  labs(x = "expected -log10 p", y = "observed -log10 p", title = sprintf("QQ (random 1%% of variants), lambda = %.3f", lambda_gc))
save_fig(p, "gwas_qq", width = 5, height = 5)

save_rds(list(file = gwas_file, header = hdr, loci = loci, lambda_gc = lambda_gc, locus_mode = locus_mode,
              n_gws = sum(thin$p_value < 5e-8), n_lead_p = nrow(sig)), "gwas")
msg("GWAS step complete")
