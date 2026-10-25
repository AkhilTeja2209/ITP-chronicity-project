# Step 13 - Replication of the variant layer in a second ITP GWAS (designed for FinnGen R12).
# Runs only when the replication file's URL is provided, as the GitHub secret named in
# config gwas.replication.url_env (FinnGen emails the download link after its access form).
#   - each discovery locus: the lead variant's effect in the replication GWAS (allele-aligned),
#     and the best variant in the same window
#   - TWAS recomputed with the replication GWAS using the same whole-blood models
# GCST90018865 already contains an earlier FinnGen release, so FinnGen R12 overlaps the discovery
# sample: agreement here is not independent replication, and the report says so.

source("R/utils.R")
source("R/twas_lib.R")

main <- function() {
  rc <- cfg$gwas$replication
  url <- Sys.getenv(rc$url_env, "")
  if (url == "") {
    msg("No replication GWAS provided (secret ", rc$url_env, " is empty); step skipped")
    write_tab(data.table(status = "not run", reason = paste("secret", rc$url_env, "not set")), "replication_status")
    return(invisible())
  }
  gw <- load_rds("gwas")
  disc_twas <- load_rds("twas")
  loci <- gw$loci
  rep_dir <- file.path(DATA, "replication")
  f <- download_once(url, file.path(rep_dir, sub("\\?.*$", "", basename(url))))

  con <- gzfile(f); hdr <- strsplit(sub("\r$", "", readLines(con, n = 1)), "\t")[[1]]; close(con)
  cols <- unlist(rc$columns)
  idx <- match(cols, hdr)
  if (anyNA(idx)) stop("Replication file lacks columns: ", paste(cols[is.na(idx)], collapse = ", "),
                       ". Its header is: ", paste(hdr, collapse = ", "))
  names(idx) <- names(cols)

  # one streaming pass: keep TWAS model positions and every discovery locus window
  models <- load_twas_models()
  keys <- file.path(rep_dir, "keys.txt")
  fwrite(unique(models$weights[, .(k = paste0(chr, ":", pos))]), keys, col.names = FALSE)
  wins <- file.path(rep_dir, "windows.txt")
  fwrite(loci[, .(chr, start, end)], wins, sep = "\t", col.names = FALSE)
  awk <- sprintf(paste0(
    "{sub(/\\r$/, \"\")} FILENAME==ARGV[1]{k[$1]=1; next} FILENAME==ARGV[2]{n++; wc[n]=$1; ws[n]=$2; we[n]=$3; next} ",
    "FNR==1{next} {c=$%d; sub(/^chr/, \"\", c); if (c==\"23\") c=\"X\"; p=$%d+0; keep=((c\":\"p) in k); ",
    "for (i=1; i<=n && !keep; i++) if (c==wc[i] && p>=ws[i] && p<=we[i]) keep=1; ",
    "if (keep) print c, p, $%d, $%d, $%d, $%d, $%d}"),
    idx["chromosome"], idx["position"], idx["effect_allele"], idx["other_allele"], idx["beta"], idx["se"], idx["p"])
  cmd <- sprintf("gzip -dc '%s' | awk -F'\\t' -v OFS='\\t' '%s' '%s' '%s' -", f, awk, keys, wins)
  g <- fread(cmd = cmd, sep = "\t", header = FALSE, na.strings = c("NA", ""),
             col.names = c("chr", "pos", "effect_allele", "other_allele", "beta", "se", "p"))
  g <- g[is.finite(beta) & is.finite(se) & se > 0]
  g[, `:=`(chr = as.character(chr), z = beta / se)]
  msg("Replication GWAS: ", nrow(g), " variants kept from ", rc$label)

  # ---- loci
  rep_loci <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
    L <- loci[i]
    w <- g[chr == L$chr & pos >= L$start & pos <= L$end]
    lead <- w[pos == L$lead_pos]
    lead[, sgn := fifelse(effect_allele == L$effect_allele & other_allele == L$other_allele, 1,
                          fifelse(effect_allele == L$other_allele & other_allele == L$effect_allele, -1, NA_real_))]
    lead <- lead[!is.na(sgn)][1]
    best <- w[which.min(p)]
    data.table(locus = L$locus, chr = L$chr, lead_rsid = L$lead_rsid, discovery_beta = L$beta, discovery_p = L$p,
               replication_beta = if (nrow(lead) && !is.na(lead$sgn)) lead$beta * lead$sgn else NA_real_,
               replication_p = if (nrow(lead) && !is.na(lead$sgn)) lead$p else NA_real_,
               best_window_p = if (nrow(best)) best$p else NA_real_,
               best_window_distance_kb = if (nrow(best)) round((best$pos - L$lead_pos) / 1000) else NA_real_,
               mhc = L$mhc)
  }))
  rep_loci[, same_direction := sign(replication_beta) == sign(discovery_beta)]
  write_tab(rep_loci, "replication_loci")

  # ---- TWAS with the replication GWAS
  rt <- twas_gene_z(models, g)
  rt[, p := 2 * pnorm(-abs(z))][, fdr := p.adjust(p, "BH")]
  cmp <- merge(disc_twas[, .(model_gene, symbol, z_discovery = z, p_discovery = p)],
               rt[, .(model_gene, z_replication = z, p_replication = p, fdr_replication = fdr)], by = "model_gene")
  setorder(cmp, p_discovery)
  write_tab(cmp, "replication_twas")

  found <- rep_loci[!is.na(replication_beta)]
  top <- head(cmp, 20)
  ct <- suppressWarnings(cor.test(cmp$z_discovery, cmp$z_replication, method = "spearman"))
  summary <- data.table(
    replication_gwas = rc$label, overlaps_discovery = isTRUE(rc$overlaps_discovery),
    loci = nrow(rep_loci), lead_found = nrow(found), lead_same_direction = sum(found$same_direction),
    lead_sign_test_p = if (nrow(found)) binom.test(sum(found$same_direction), nrow(found), 0.5, alternative = "greater")$p.value else NA_real_,
    lead_p05_same_direction = sum(found$same_direction & found$replication_p < 0.05),
    twas_genes_compared = nrow(cmp), twas_spearman_rho = unname(ct$estimate), twas_spearman_p = ct$p.value,
    twas_top20_same_direction = sum(sign(top$z_discovery) == sign(top$z_replication)),
    twas_replication_fdr_genes = sum(cmp$fdr_replication < cfg$integration$twas_fdr))
  write_tab(summary, "replication_summary")
  write_tab(data.table(status = "run", reason = rc$label), "replication_status")
  print(summary)
}
main()
