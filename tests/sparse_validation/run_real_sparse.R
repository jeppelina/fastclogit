# Run the full validation suite with the REAL sparse fitter
# (fastclogit dispatching on dgCMatrix input).

suppressPackageStartupMessages({
  library(fastclogit); library(Matrix); library(survival); library(data.table)
})

source("tests/sparse_validation/gen_data.R")
source("tests/sparse_validation/compare.R")
source("tests/sparse_validation/reference_fits.R")
source("tests/sparse_validation/run_all.R")

# Real sparse fitter — passes dgCMatrix directly; fastclogit dispatches
# to clogit_fit_sparse_cpp via is(X, "sparseMatrix") check.
real_sparse_fit <- function(Xsp, choice, strata, cluster = NULL) {
  stopifnot(inherits(Xsp, "sparseMatrix"))
  fastclogit::fastclogit(
    X = Xsp, choice = choice, strata = strata, cluster = cluster,
    max_iter = 100L, tol = 1e-8, verbose = FALSE
  )
}

cat("\n", strrep("#", 72), "\n", sep="")
cat("# REAL sparse path validation: clogit_fit_sparse_cpp\n")
cat(strrep("#", 72), "\n", sep="")
results <- run_sparse_validation(
  sparse_fit_fn = real_sparse_fit,
  tol_coef = 1e-10, tol_se = 1e-8, tol_ll = 1e-6, verbose = TRUE
)
