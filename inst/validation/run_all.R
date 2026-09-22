# =============================================================================
# run_all.R — run every validation study and print one summary.
#
#   cd inst/validation && Rscript run_all.R
#
# Each study prints its own pre-specified checks as it goes; the summary at the
# end aggregates them and writes results_summary.csv.
#
# Total runtime is roughly 8-10 minutes, dominated by 04_offset.R.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
Sys.setenv(VC_DIR = vc_dir)

studies <- c(
  "01_dgp_equivalence.R",
  "02_design_matrix.R",
  "03_recovery_coverage.R",
  "04_offset.R",
  "05_cluster_robust.R",
  "06_khb.R",
  "07_benchmarks.R",
  "08_assumptions.R"
)

only <- commandArgs(trailingOnly = TRUE)
if (length(only)) studies <- studies[grepl(paste(only, collapse = "|"), studies)]

source(file.path(vc_dir, "00_helpers.R"))
t0 <- Sys.time()
for (s in studies) {
  # Each study is sourced into this session so vc_check() accumulates into one
  # shared results table.
  source(file.path(vc_dir, s), local = new.env(parent = globalenv()))
}
cat(sprintf("\n  total runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
vc_summary(file.path(vc_dir, "results_summary.csv"))
