###############################################################################
#### load_fastclogit.R — Load fastclogit on MONA without package install
####
#### Compiles the C++ files via Rcpp::sourceCpp() and sources the R files.
#### Runs a quick sanity check on simulated data to verify everything works.
####
#### Upload this file + the files below to MONA, then:
####   source("load_fastclogit.R")
####
#### Required files in the same directory:
####   clogit_newton.cpp      (Newton-Raphson with three-tier convergence)
####   clogit_sandwich.cpp    (Clustered sandwich variance estimation)
####   fastclogit.R           (Core: matrix-level fitting function)
####   fastclogit_methods.R   (S3 methods: summary, print, vcov, confint, tidy)
####   fclogit.R              (Formula interface: fclogit())
####   khb_decompose.R        (Generic KHB mediation decomposition)
####
#### Optional files (in archive/ subfolder):
####   simulate_clogit.R      (DGP for testing)
####
#### Dependencies: Rcpp, RcppArmadillo (both available on MONA / SCB servers)
###############################################################################

cat("=== Loading fastclogit (sourceCpp mode) ===\n")

# --- Locate the kernels -------------------------------------------------
# Try several candidates and take the first that DEMONSTRABLY holds the .cpp
# files, rather than trusting one guess.
#
# The old version trusted sys.frame(1)$ofile alone. sys.frame(1) is the
# OUTERMOST frame, not the calling one, so with nested source() calls it
# names whichever file started the chain. Sourced from a job script in a
# subdirectory (job -> helper -> wrapper -> engine -> here) it resolved to
# that subdirectory and the load died with "Cannot find:
# .../clogit_newton.cpp", 0.3 s into a job budgeted at four hours
# (Paper 3, 2026-09-10). It happened to work from a top-level source() only
# because the frame was then load_fastclogit.R's own.
#
# Checking for the file removes the guesswork: a candidate is either right or
# it is skipped. Set .this_script_dir or CODE_PATH before sourcing to pin it.
.fc_candidates <- c(
  if (exists(".this_script_dir")) .this_script_dir else NULL,
  if (exists("CODE_PATH")) CODE_PATH else NULL,
  if (exists("CODE_PATH")) file.path(CODE_PATH, "lib") else NULL,
  getwd(),
  file.path(getwd(), "lib"),
  tryCatch({
    sp <- sys.frame(1)$ofile
    if (!is.null(sp)) dirname(normalizePath(sp)) else NULL
  }, error = function(e) NULL)
)
.fc_found <- Filter(function(d)
  !is.null(d) && nzchar(d) && file.exists(file.path(d, "clogit_newton.cpp")),
  .fc_candidates)
if (!length(.fc_found))
  stop("load_fastclogit.R cannot locate clogit_newton.cpp. Looked in:\n  ",
       paste(unique(unlist(.fc_candidates)), collapse = "\n  "),
       "\nSet .this_script_dir to the directory holding the .cpp files ",
       "before sourcing.")
FASTCLOGIT_DIR <- .fc_found[[1]]
rm(.fc_candidates, .fc_found)

cat("  Directory:", FASTCLOGIT_DIR, "\n")

# --- Check dependencies ---
cat("  Checking Rcpp... ")
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp is not installed. Ask SCB to install it, or try: install.packages('Rcpp')\n",
       "Rcpp and RcppArmadillo are required on MONA.")
}
cat(as.character(packageVersion("Rcpp")), "\n")

cat("  Checking RcppArmadillo... ")
if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
  stop("RcppArmadillo is not installed. Ask SCB to install it, or try: install.packages('RcppArmadillo')\n",
       "Rcpp and RcppArmadillo are required on MONA.")
}
cat(as.character(packageVersion("RcppArmadillo")), "\n")

# --- Compile C++ files ---
cpp_newton          <- file.path(FASTCLOGIT_DIR, "clogit_newton.cpp")
cpp_newton_sparse   <- file.path(FASTCLOGIT_DIR, "clogit_newton_sparse.cpp")
cpp_sandwich        <- file.path(FASTCLOGIT_DIR, "clogit_sandwich.cpp")
cpp_sandwich_sparse <- file.path(FASTCLOGIT_DIR, "clogit_sandwich_sparse.cpp")
# (v0.4.1: no csr_matrix.h here — its contents are inlined into the
# sparse .cpp files to avoid include-path issues on UNC paths with '$'.)

if (!file.exists(cpp_newton))   stop("Cannot find: ", cpp_newton)
if (!file.exists(cpp_sandwich)) stop("Cannot find: ", cpp_sandwich)

# MONA's source-mode install: all four .cpp files are self-contained.
# The sparse kernels inline the CsrMatrix struct rather than relying on
# a csr_matrix.h header file — this avoids include-path failures on UNC
# paths that R/make may rewrite (e.g. paths containing '$').
#
# env = globalenv() makes compiled functions land in the global env
# regardless of whether load_fastclogit.R is sourced or run via Rscript.

cat("  Compiling clogit_newton.cpp (dense)... ")
Rcpp::sourceCpp(cpp_newton, env = globalenv())
cat("OK\n")

cat("  Compiling clogit_sandwich.cpp (dense)... ")
Rcpp::sourceCpp(cpp_sandwich, env = globalenv())
cat("OK\n")

# --- Compile sparse C++ files if present ---------------------------------
# The sparse path is optional and needs BOTH files: fastclogit() dispatches
# on inherits(X, "sparseMatrix") and would find only half a kernel if one
# were missing. Compile both or neither.
if (file.exists(cpp_newton_sparse) && file.exists(cpp_sandwich_sparse)) {
  cat("  Compiling clogit_newton_sparse.cpp... ")
  Rcpp::sourceCpp(cpp_newton_sparse, env = globalenv())
  cat("OK\n")
  cat("  Compiling clogit_sandwich_sparse.cpp... ")
  Rcpp::sourceCpp(cpp_sandwich_sparse, env = globalenv())
  cat("OK\n")

  # Smoke check: confirm both entry points actually reached the global env.
  if (!exists("clogit_fit_sparse_cpp", mode = "function"))
    stop("Sparse kernel compiled but clogit_fit_sparse_cpp not exported. ",
         "Check that the file's [[Rcpp::export]] attribute is intact.")
  if (!exists("clogit_sandwich_sparse_cpp", mode = "function"))
    stop("Sparse kernel compiled but clogit_sandwich_sparse_cpp not exported.")
} else {
  cat("  Sparse kernel files not found — sparse-X dispatch will be disabled.\n")
  cat("    Looked for: ", basename(cpp_newton_sparse),
      " and ", basename(cpp_sandwich_sparse), "\n", sep = "")
}

# --- Source R files ---
.source_if_exists <- function(filename, label = NULL) {
  path <- file.path(FASTCLOGIT_DIR, filename)
  if (file.exists(path)) {
    source(path)
    cat("  Sourced:", if (!is.null(label)) label else filename, "\n")
    return(TRUE)
  }
  return(FALSE)
}

# Core (required)
r_wrapper <- file.path(FASTCLOGIT_DIR, "fastclogit.R")
if (!file.exists(r_wrapper)) stop("Cannot find: ", r_wrapper)
source(r_wrapper)
cat("  Sourced: fastclogit.R (core fitting function)\n")

# S3 methods (required)
.source_if_exists("fastclogit_methods.R", "fastclogit_methods.R (S3 methods)")

# Formula interface
.source_if_exists("fclogit.R", "fclogit.R (formula interface)")

# DGP simulator (check archive/ subfolder too)
if (!.source_if_exists("simulate_clogit.R", "simulate_clogit.R (data simulator)")) {
  .source_if_exists("archive/simulate_clogit.R", "archive/simulate_clogit.R (data simulator)")
}

# KHB decomposition
.source_if_exists("khb_decompose.R", "khb_decompose.R (KHB mediation)")

# --- Quick sanity check ---
# Needs simulate_clogit.R, which is deliberately not deployed on MONA. Without
# this guard the check called a function that does not exist and printed
# "Sanity check: FAIL", which reads as a broken kernel when in fact all four
# .cpp kernels had just compiled cleanly. Say SKIPPED, and mean it.
cat("  Sanity check: ")
if (!exists("simulate_clogit_data")) {
  cat("SKIPPED (simulate_clogit.R not deployed; kernels compiled OK)\n")
} else
tryCatch({
  test_sim <- simulate_clogit_data(n_egos = 50, n_alts = 10, seed = 1)
  test_fit <- fastclogit(test_sim$X, test_sim$choice, test_sim$strata)
  stopifnot(test_fit$converged)

  # Also test formula interface if available
  if (exists("fclogit")) {
    test_fit2 <- fclogit(choice ~ x1 + x2,
                         data = test_sim$data, strata = "strata_id")
    stopifnot(test_fit2$converged)
    cat("PASS (matrix + formula, ", length(test_fit$coefficients), " coefs, loglik = ",
        round(test_fit$loglik, 2), ")\n", sep = "")
    rm(test_fit2)
  } else {
    cat("PASS (", length(test_fit$coefficients), " coefs, loglik = ",
        round(test_fit$loglik, 2), ")\n", sep = "")
  }
  rm(test_sim, test_fit)
}, error = function(e) {
  cat("FAIL: ", conditionMessage(e), "\n")
  warning("fastclogit sanity check failed — functions are loaded but may not work correctly")
})

cat("=== fastclogit ready ===\n")
cat("  Available functions:\n")
cat("    fastclogit()         - Matrix interface (X, choice, strata)\n")
if (exists("fclogit"))       cat("    fclogit()            - Formula interface (choice ~ x1 + x2)\n")
if (exists("khb_decompose")) cat("    khb_decompose()      - KHB mediation decomposition\n")
if (exists("tidy_fastclogit")) cat("    tidy_fastclogit()    - Broom-style tidy output\n")
if (exists("simulate_clogit_data")) cat("    simulate_clogit_data() - Data simulator for testing\n")
cat("\n")
