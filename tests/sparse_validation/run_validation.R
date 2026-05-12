#!/usr/bin/env Rscript
# =============================================================================
# run_validation.R -- Shell entry point.
#
#   Rscript tests/sparse_validation/run_validation.R
#
# Exits 0 if all checks pass, 1 otherwise.
# Default sparse-fitter is the dense fastclogit (identity check / harness test).
# =============================================================================

.this_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grepl("^--file=", args)]
  if (length(fa) > 0) return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  for (i in seq_len(sys.nframe())) {
    fr <- sys.frame(i)
    if (!is.null(fr$ofile)) return(normalizePath(dirname(fr$ofile)))
  }
  normalizePath(file.path("tests", "sparse_validation"))
})()

source(file.path(.this_dir, "run_all.R"))

res <- run_sparse_validation(verbose = FALSE)

if (isTRUE(res$all_pass)) {
  quit(save = "no", status = 0L)
} else {
  quit(save = "no", status = 1L)
}
