# TWAS helpers shared by step 11 (discovery GWAS) and step 13 (replication GWAS).
# S-PrediXcan statistic (Barbeira et al. 2018 Nat Commun):
#     z_gene = sum_l  w_l * sigma_l * z_l / sigma_gene,   sigma_gene^2 = w' Sigma w

suppressPackageStartupMessages({ library(DBI); library(RSQLite) })

# Download (once) and open the PredictDB MASHR model for the configured tissue.
load_twas_models <- function() {
  twas_dir <- file.path(DATA, "twas")
  tarball <- download_once(cfg$twas$url, file.path(twas_dir, "mashr_eqtl.tar"))
  members <- untar(tarball, list = TRUE)
  want <- grep(paste0("mashr_", cfg$twas$tissue, "\\.(db|txt\\.gz)$"), members, value = TRUE)
  if (length(want) != 2) stop("Tissue model files not found in tarball: ", cfg$twas$tissue)
  if (!all(file.exists(file.path(twas_dir, want)))) untar(tarball, files = want, exdir = twas_dir)
  con <- dbConnect(SQLite(), file.path(twas_dir, grep("\\.db$", want, value = TRUE)))
  on.exit(dbDisconnect(con))
  weights <- as.data.table(dbGetQuery(con, "SELECT gene, varID, ref_allele, eff_allele, weight FROM weights"))
  extra <- as.data.table(dbGetQuery(con, "SELECT * FROM extra"))
  weights[, c("chr", "pos") := tstrsplit(sub("^chr", "", varID), "_", keep = 1:2)]
  weights[, pos := as.integer(pos)]
  list(weights = weights, extra = extra, dir = twas_dir,
       cov_file = file.path(twas_dir, grep("\\.txt\\.gz$", want, value = TRUE)))
}

# g: data.table with chr (character), pos, effect_allele, other_allele, z.
# Returns one row per gene: model_gene, z, n_snps_used.
twas_gene_z <- function(models, g) {
  w <- models$weights
  m <- merge(w, g[, .(chr, pos, effect_allele, other_allele, z)], by = c("chr", "pos"), allow.cartesian = TRUE)
  # model eff_allele must match the GWAS effect_allele; flip the sign when the alleles are swapped
  m[, sign := fifelse(effect_allele == eff_allele & other_allele == ref_allele, 1,
                      fifelse(effect_allele == ref_allele & other_allele == eff_allele, -1, NA_real_))]
  m <- m[!is.na(sign)][, z := z * sign]
  m <- unique(m, by = c("gene", "varID"))
  covs <- fread(models$cov_file)
  setnames(covs, c("GENE", "RSID1", "RSID2", "VALUE"))
  covs <- covs[GENE %in% unique(m$gene)]
  cov_by_gene <- split(covs, by = "GENE", keep.by = FALSE)
  m_by_gene <- split(m, by = "gene")
  rbindlist(lapply(names(m_by_gene), function(gn) {
    d <- m_by_gene[[gn]]; cv <- cov_by_gene[[gn]]
    if (is.null(cv)) return(NULL)
    ids <- d$varID
    S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    cv <- cv[RSID1 %in% ids & RSID2 %in% ids]
    S[cbind(cv$RSID1, cv$RSID2)] <- cv$VALUE
    S[cbind(cv$RSID2, cv$RSID1)] <- cv$VALUE
    sd_l <- sqrt(pmax(diag(S), 0)); var_g <- drop(t(d$weight) %*% S %*% d$weight)
    if (!is.finite(var_g) || var_g <= 0) return(NULL)
    data.table(model_gene = gn, z = sum(d$weight * sd_l * d$z) / sqrt(var_g), n_snps_used = nrow(d))
  }))
}
