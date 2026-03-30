###############################################################################
#### test_pure_vs_rcpp.R — Verify pure-R and Rcpp versions give same results
####
#### Run this on your Mac (where both Rcpp and pure R work):
####   cd ~/Research/fastclogit/mona
####   Rscript test_pure_vs_rcpp.R
####
#### Sources the pure-R version (fastclogit_pure.R) in its own environment,
#### then compares coefficients, SEs, and robust SEs against the installed
#### Rcpp package across 5 configurations (basic, offset, cluster, large, small).
####
#### NOTE: pure-R version uses explicit stats:: prefixes (var, pnorm, etc.)
#### to work correctly when loaded into restricted environments.
###############################################################################
setwd("/Users/jeppelina/Documents/R/Research/fastclogit/mona")

cat("=== Comparing pure-R vs Rcpp fastclogit ===\n\n")

# --- Load Rcpp version (installed package) ---
library(fastclogit)
cat("Rcpp version loaded: fastclogit", as.character(packageVersion("fastclogit")), "\n")

# Save Rcpp functions under different names
fastclogit_rcpp <- fastclogit::fastclogit
simulate_rcpp   <- fastclogit::simulate_clogit_data

# --- Load pure-R version into its own environment ---
pure_env <- new.env(parent = baseenv())
sys.source("fastclogit_pure.R", envir = pure_env)
fastclogit_pure <- pure_env$fastclogit
simulate_pure   <- pure_env$simulate_clogit_data
cat("Pure-R version loaded\n\n")

# --- Test configurations ---
tests <- list(
  list(label = "Basic (no offset, no cluster)",
       n_egos = 500, n_alts = 20, use_offset = FALSE, cluster_ratio = 1.0, seed = 100),
  list(label = "With offset",
       n_egos = 500, n_alts = 20, use_offset = TRUE, cluster_ratio = 1.0, seed = 200),
  list(label = "With offset + cluster",
       n_egos = 500, n_alts = 20, use_offset = TRUE, cluster_ratio = 1.5, seed = 300),
  list(label = "Larger (2K egos, 50 alts)",
       n_egos = 2000, n_alts = 50, use_offset = TRUE, cluster_ratio = 1.3, seed = 400),
  list(label = "Small groups (2 alts)",
       n_egos = 1000, n_alts = 2, use_offset = FALSE, cluster_ratio = 1.0, seed = 500)
)

all_pass <- TRUE

for (tt in tests) {
  cat("Test:", tt$label, "\n")

  # Simulate data (use Rcpp version — same DGP code, same seed)
  sim <- simulate_rcpp(
    n_egos = tt$n_egos, n_alts = tt$n_alts,
    use_offset = tt$use_offset, cluster_ratio = tt$cluster_ratio,
    seed = tt$seed
  )

  off <- if (tt$use_offset) sim$offset else NULL
  cl  <- if (tt$cluster_ratio > 1) sim$cluster else NULL

  # Fit with Rcpp
  t_rcpp <- system.time({
    fit_rcpp <- fastclogit_rcpp(sim$X, sim$choice, sim$strata,
                                offset = off, cluster = cl,
                                max_iter = 50, tol = 1e-8)
  })

  # Fit with pure R
  t_pure <- system.time({
    fit_pure <- fastclogit_pure(sim$X, sim$choice, sim$strata,
                                offset = off, cluster = cl,
                                max_iter = 50, tol = 1e-8)
  })

  # Compare
  max_coef_diff <- max(abs(fit_rcpp$coefficients - fit_pure$coefficients))
  ll_diff <- abs(fit_rcpp$loglik - fit_pure$loglik)
  max_se_diff <- max(abs(fit_rcpp$se - fit_pure$se) / fit_rcpp$se)

  cat("  Coefficients max|diff|:", format(max_coef_diff, digits = 3, scientific = TRUE), "\n")
  cat("  Log-lik |diff|:       ", format(ll_diff, digits = 3, scientific = TRUE), "\n")
  cat("  Model SE max|rel diff|:", format(max_se_diff, digits = 3, scientific = TRUE), "\n")

  if (!is.null(cl)) {
    max_rse_diff <- max(abs(fit_rcpp$se_robust - fit_pure$se_robust) / fit_rcpp$se_robust)
    cat("  Robust SE max|rel diff|:", format(max_rse_diff, digits = 3, scientific = TRUE), "\n")
  }

  cat("  Time: Rcpp =", round(t_rcpp[3], 3), "s | pure R =", round(t_pure[3], 3), "s",
      "| ratio =", round(t_pure[3] / max(t_rcpp[3], 0.001), 1), "x\n")

  pass <- max_coef_diff < 1e-8 && ll_diff < 1e-8
  cat("  Result:", if (pass) "PASS" else "*** FAIL ***", "\n\n")
  if (!pass) all_pass <- FALSE
}

cat("============================================================\n")
cat("  Overall:", if (all_pass) "ALL PASS" else "SOME FAILED", "\n")
cat("============================================================\n")
