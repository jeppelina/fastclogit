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
####   clogit_newton.cpp      (Newton-Raphson with ridge regularization)
####   clogit_sandwich.cpp    (Clustered sandwich variance estimation)
####   fastclogit.R           (Core: matrix-level fitting function)
####   fastclogit_methods.R   (S3 methods: summary, print, vcov, confint, tidy)
####   fclogit.R              (Formula interface: fclogit())
####   simulate_clogit.R      (DGP for testing — optional)
####
#### Optional files:
####   khb_decompose.R        (Generic KHB mediation decomposition)
####   fastclogit_pure.R      (Pure-R fallback if Rcpp is not available)
####
#### Dependencies: Rcpp, RcppArmadillo (both available on MONA / SCB servers)
####
#### If Rcpp is NOT available, source fastclogit_pure.R instead of this file.
#### It provides the same interface but runs in pure R (slower, no compiler).
###############################################################################

cat("=== Loading fastclogit (sourceCpp mode) ===\n")

# --- Locate files relative to this script ---
if (exists(".this_script_dir")) {
  FASTCLOGIT_DIR <- .this_script_dir
} else {
  FASTCLOGIT_DIR <- tryCatch({
    script_path <- sys.frame(1)$ofile
    if (!is.null(script_path)) dirname(normalizePath(script_path))
    else getwd()
  }, error = function(e) getwd())
}

cat("  Directory:", FASTCLOGIT_DIR, "\n")

# --- Check dependencies ---
cat("  Checking Rcpp... ")
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp is not installed. Ask SCB to install it, or try: install.packages('Rcpp')\n",
       "Alternatively, source('fastclogit_pure.R') for a pure-R version (no Rcpp needed).")
}
cat(as.character(packageVersion("Rcpp")), "\n")

cat("  Checking RcppArmadillo... ")
if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
  stop("RcppArmadillo is not installed. Ask SCB to install it, or try: install.packages('RcppArmadillo')\n",
       "Alternatively, source('fastclogit_pure.R') for a pure-R version (no Rcpp needed).")
}
cat(as.character(packageVersion("RcppArmadillo")), "\n")

# --- Compile C++ files ---
cpp_newton   <- file.path(FASTCLOGIT_DIR, "clogit_newton.cpp")
cpp_sandwich <- file.path(FASTCLOGIT_DIR, "clogit_sandwich.cpp")

if (!file.exists(cpp_newton))   stop("Cannot find: ", cpp_newton)
if (!file.exists(cpp_sandwich)) stop("Cannot find: ", cpp_sandwich)

cat("  Compiling clogit_newton.cpp... ")
Rcpp::sourceCpp(cpp_newton)
cat("OK\n")

cat("  Compiling clogit_sandwich.cpp... ")
Rcpp::sourceCpp(cpp_sandwich)
cat("OK\n")

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

# DGP simulator
.source_if_exists("simulate_clogit.R", "simulate_clogit.R (data simulator)")

# KHB decomposition
.source_if_exists("khb_decompose.R", "khb_decompose.R (KHB mediation)")

# --- Quick sanity check ---
cat("  Sanity check: ")
tryCatch({
  test_sim <- simulate_clogit_data(n_egos = 50, n_alts = 10, seed = 1)
  test_fit <- fastclogit(test_sim$X, test_sim$choice, test_sim$strata)
  stopifnot(test_fit$converged)

  # Also test formula interface if available
  if (exists("fclogit")) {
    test_fit2 <- fclogit(choice ~ lnDist + n_years_same_cfar,
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
