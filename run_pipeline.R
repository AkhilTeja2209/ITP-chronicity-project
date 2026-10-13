#!/usr/bin/env Rscript
# Runs the full pipeline. Usage:
#   Rscript run_pipeline.R            # all steps
#   Rscript run_pipeline.R 10 11 12   # only the listed steps (upstream outputs must already exist)
# The expression arm (00-03) and the variant arm (10-12) are independent: a failure in one is
# recorded and the other still runs. The report is always written; the exit code is non-zero if
# any step failed, so CI shows red.

steps <- list(
  "00" = list(file = "R/00_audit.R",                 needs = character()),
  "01" = list(file = "R/01_expression.R",            needs = "00"),
  "02" = list(file = "R/02_differential_expression.R", needs = "01"),
  "03" = list(file = "R/03_modules_enrichment.R",    needs = "02"),
  "10" = list(file = "R/10_gwas.R",                  needs = "01"),   # needs gene coordinates from 01
  "11" = list(file = "R/11_twas.R",                  needs = "10"),
  "12" = list(file = "R/12_coloc.R",                 needs = "10"),
  "20" = list(file = "R/20_integration.R",           needs = "02"),
  "21" = list(file = "R/21_benchmark.R",             needs = "02"))

args <- commandArgs(trailingOnly = TRUE)
run <- if (length(args)) intersect(names(steps), args) else names(steps)
status <- data.frame(step = character(), script = character(), status = character(),
                     minutes = numeric(), message = character())
failed <- character()

for (id in run) {
  st <- steps[[id]]
  blocked <- intersect(st$needs, failed)
  t0 <- Sys.time()
  if (length(blocked)) {
    res <- list(status = "skipped", message = paste("upstream failed:", paste(blocked, collapse = ",")))
    failed <- c(failed, id)
  } else {
    message("\n========== step ", id, ": ", st$file, " ==========")
    tb <- character()
    res <- tryCatch({
      withCallingHandlers(sys.source(st$file, envir = new.env(parent = globalenv())),
        error = function(e) tb <<- vapply(sys.calls(), function(cl) paste(deparse(cl, nlines = 1), collapse = ""), ""))
      list(status = "ok", message = "")
    }, error = function(e) {
      message("ERROR in step ", id, ": ", conditionMessage(e))
      message("Traceback (innermost last):\n  ", paste(tail(head(tb, -2), 8), collapse = "\n  "))
      list(status = "FAILED", message = gsub("[\r\n\t|]+", " ", substr(conditionMessage(e), 1, 300)))
    })
    if (res$status != "ok") failed <- c(failed, id)
  }
  status <- rbind(status, data.frame(step = id, script = st$file, status = res$status,
                                     minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                                     message = res$message))
  dir.create("results", showWarnings = FALSE)
  write.table(status, "results/pipeline_status.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

message("\n========== report ==========")
tryCatch(sys.source("R/30_report.R", envir = new.env(parent = globalenv())),
         error = function(e) { message("Report failed: ", conditionMessage(e)); failed <<- c(failed, "30") })
print(status)
if (length(failed)) quit(status = 1)
