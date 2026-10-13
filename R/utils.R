# Shared helpers sourced by every pipeline step.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(ggplot2)
})

cfg <- yaml::read_yaml("config/config.yaml")
gene_sets <- yaml::read_yaml("config/gene_sets.yaml")
set.seed(cfg$seed)
options(timeout = 7200, stringsAsFactors = FALSE)

DATA <- cfg$paths$data
RES <- cfg$paths$results
for (d in c(DATA, file.path(RES, c("tables", "figures", "rds"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

write_tab <- function(x, name) {
  f <- file.path(RES, "tables", paste0(name, ".tsv"))
  fwrite(as.data.table(x), f, sep = "\t", na = "NA")
  invisible(f)
}
read_tab <- function(name) fread(file.path(RES, "tables", paste0(name, ".tsv")))

save_rds <- function(x, name) saveRDS(x, file.path(RES, "rds", paste0(name, ".rds")))
load_rds <- function(name) {
  f <- file.path(RES, "rds", paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing upstream output '", name, "'. Run the step that creates it first.")
  readRDS(f)
}
has_rds <- function(name) file.exists(file.path(RES, "rds", paste0(name, ".rds")))

save_fig <- function(plot, name, width = 7, height = 5) {
  ggsave(file.path(RES, "figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 150, bg = "white")
}
# For base-graphics figures (WGCNA dendrograms etc.)
with_png <- function(name, expr, width = 8, height = 6) {
  png(file.path(RES, "figures", paste0(name, ".png")), width = width, height = height,
      units = "in", res = 150)
  on.exit(dev.off())
  force(expr)
}

download_once <- function(url, dest, retries = 3) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")
  for (i in seq_len(retries)) {
    msg("Downloading ", url, if (i > 1) paste0(" (attempt ", i, ")"))
    ok <- tryCatch({
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl")
      TRUE
    }, error = function(e) { msg("  failed: ", conditionMessage(e)); FALSE })
    if (ok && file.size(tmp) > 0) { file.rename(tmp, dest); return(dest) }
    Sys.sleep(10 * i)
  }
  stop("Could not download ", url)
}

# Score each sample as the mean z-score of the set's genes. `mat` is genes x samples with
# rownames = gene symbols. Returns NA when fewer than `min_genes` genes are present.
score_set <- function(mat, genes, min_genes = 2) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < min_genes) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

# Collapse a probe x sample matrix to gene x sample keeping, per gene, the probe with highest mean.
collapse_to_genes <- function(mat, gene_of_probe) {
  keep <- !is.na(gene_of_probe) & gene_of_probe != ""
  mat <- mat[keep, , drop = FALSE]
  g <- gene_of_probe[keep]
  avg <- rowMeans(mat, na.rm = TRUE)
  ord <- order(g, -avg)
  first <- ord[!duplicated(g[ord])]
  out <- mat[first, , drop = FALSE]
  rownames(out) <- g[first]
  out
}

# Markdown table without extra dependencies.
md_table <- function(df, digits = 3, max_rows = 25) {
  df <- as.data.frame(df)
  if (nrow(df) == 0) return("_(no rows)_")
  df <- head(df, max_rows)
  fmt <- function(x) {
    if (is.numeric(x)) {
      ifelse(is.na(x), "NA",
             ifelse(abs(x) < 1e-3 & x != 0, formatC(x, format = "e", digits = 2),
                    formatC(x, format = "fg", digits = digits)))
    } else as.character(x)
  }
  body <- vapply(seq_len(nrow(df)), function(i)
    paste0("| ", paste(vapply(df[i, , drop = FALSE], fmt, ""), collapse = " | "), " |"), "")
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body) |> paste(collapse = "\n")
}

theme_set(theme_bw(base_size = 11))
