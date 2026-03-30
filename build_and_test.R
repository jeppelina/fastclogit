#!/usr/bin/env Rscript
# ============================================================================
# build_and_test.R — Build fastclogit and run validation tests
#
# Run from the fastclogit/ directory:
#   Rscript build_and_test.R
#
# Or in RStudio: open the fastclogit.Rproj, then source this file.
#
# All output is logged to build_and_test_<timestamp>.log
# ============================================================================

# --- Logging setup ---
timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(getwd(), paste0("build_and_test_", timestamp_str, ".log"))

# Open log connection — tee to both console and file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)  # split = TRUE → console + file
sink(log_con, type = "message")                # capture warnings/errors too

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  # Print log location to console after sinks are closed
  message("\nLog saved to: ", log_file)
}, add = TRUE)

cat("=== fastclogit: Build & Test ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform: ", R.version$platform, "\n")
cat("Log file: ", log_file, "\n\n")

# --- 0. Check prerequisites ---
cat("Step 0: Checking prerequisites...\n")
for (pkg in c("Rcpp", "RcppArmadillo", "survival", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  cat("  ", pkg, ":", as.character(packageVersion(pkg)), "\n")
}
cat("  All prerequisites OK\n\n")

# --- 1. Quick compile test (sourceCpp) ---
cat("Step 1: Testing C++ compilation with Rcpp::sourceCpp()...\n")
tryCatch({
  Rcpp::sourceCpp(code = '
    #include <RcppArmadillo.h>
    // [[Rcpp::depends(RcppArmadillo)]]
    // [[Rcpp::export]]
    int test_armadillo() {
      arma::mat A(3, 3, arma::fill::eye);
      return A.n_rows;
    }
  ')
  stopifnot(test_armadillo() == 3)
  cat("  C++ compilation works!\n\n")
}, error = function(e) {
  cat("  ERROR: C++ compilation failed:\n  ", conditionMessage(e), "\n")
  cat("  Make sure you have a C++ compiler (g++ or clang++) installed.\n")
  cat("  On Mac: xcode-select --install\n")
  cat("  On Ubuntu: sudo apt install r-base-dev\n")
  cat("  On Windows: install Rtools from https://cran.r-project.org/bin/windows/Rtools/\n")
  stop("Cannot proceed without C++ compilation")
})

# --- 2. Build and install the package ---
cat("Step 2: Building and installing fastclogit...\n")
pkg_dir <- getwd()
if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  if (file.exists(file.path(dirname(pkg_dir), "DESCRIPTION"))) {
    pkg_dir <- dirname(pkg_dir)
  } else {
    stop("Cannot find DESCRIPTION. Run this script from the fastclogit/ directory.")
  }
}
cat("  Package dir:", pkg_dir, "\n")

install.packages(pkg_dir, repos = NULL, type = "source")
library(fastclogit)
cat("  Package installed and loaded! Version:", as.character(packageVersion("fastclogit")), "\n\n")

# --- 3. Quick smoke test ---
cat("Step 3: Smoke test...\n")
sim <- simulate_clogit_data(n_egos = 200, n_alts = 20, seed = 1)
fit <- fastclogit(sim$X, sim$choice, sim$strata,
                   offset = sim$offset, cluster = sim$cluster,
                   verbose = TRUE)
cat("\n")
print(fit)
cat("\n  Smoke test passed!\n\n")

# --- 4. Validation against survival::clogit ---
cat("Step 4: Validating against survival::clogit...\n\n")
library(survival)

run_comparison <- function(label, n_egos, n_alts, use_offset, cluster_ratio, seed) {
  cat("  Test: ", label, "\n")
  cat("    Config: n_egos=", n_egos, " n_alts=", n_alts,
      " offset=", use_offset, " cluster_ratio=", cluster_ratio,
      " seed=", seed, "\n")

  sim <- simulate_clogit_data(
    n_egos = n_egos, n_alts = n_alts,
    use_offset = use_offset, cluster_ratio = cluster_ratio, seed = seed
  )

  # fastclogit
  t_fast <- system.time({
    fit_fast <- fastclogit(sim$X, sim$choice, sim$strata,
                            offset = if (use_offset) sim$offset else NULL,
                            cluster = if (cluster_ratio > 1) sim$cluster else NULL,
                            max_iter = 50, tol = 1e-10)
  })

  # survival::clogit
  pred_names <- colnames(sim$X)
  fml_str <- paste("choice ~",
                    paste(pred_names, collapse = " + "))
  if (use_offset) fml_str <- paste(fml_str, "+ offset(correction)")
  fml_str <- paste(fml_str, "+ strata(strata_id)")
  if (cluster_ratio > 1) fml_str <- paste(fml_str, "+ cluster(cluster_id)")
  cat("    Formula:", fml_str, "\n")

  t_clogit <- system.time({
    fit_clogit <- clogit(as.formula(fml_str), data = sim$data, method = "efron")
  })

  # Compare coefficients
  beta_fast   <- fit_fast$coefficients
  beta_clogit <- coef(fit_clogit)
  common <- intersect(names(beta_fast), names(beta_clogit))

  max_coef_diff <- max(abs(beta_fast[common] - beta_clogit[common]))
  ll_diff <- abs(fit_fast$loglik - as.numeric(fit_clogit$loglik[2]))

  cat("    Coefficients: max|diff| =", format(max_coef_diff, digits = 3, scientific = TRUE), "\n")
  cat("    Log-lik:      |diff|    =", format(ll_diff, digits = 3, scientific = TRUE), "\n")
  cat("    Converged:   fastclogit =", fit_fast$converged,
      " (", fit_fast$iterations, " iters)\n")
  cat("    Time:         fastclogit =", round(t_fast[3], 3), "s",
      " | clogit =", round(t_clogit[3], 3), "s",
      " | speedup =", round(t_clogit[3] / max(t_fast[3], 0.001), 1), "x\n")

  # Compare SEs
  # When cluster() is used, coxph stores:
  #   $naive.var = model-based (information) variance
  #   $var       = robust (sandwich) variance
  # vcov() returns $var if it exists (i.e., the robust one)
  # Without cluster(), $var = model-based and $naive.var doesn't exist
  if (cluster_ratio > 1) {
    # Model-based SEs: use naive.var
    se_fast   <- fit_fast$se[common]
    se_clogit <- setNames(sqrt(diag(fit_clogit$naive.var)), names(coef(fit_clogit)))[common]
    max_se_diff <- max(abs(se_fast - se_clogit) / se_clogit)
    cat("    SEs (model):  max|rel diff| =", format(max_se_diff, digits = 3, scientific = TRUE), "\n")

    # Robust SEs: use $var
    se_fast_r   <- fit_fast$se_robust[common]
    se_clogit_r <- setNames(sqrt(diag(fit_clogit$var)), names(coef(fit_clogit)))[common]
    max_se_r_diff <- max(abs(se_fast_r - se_clogit_r) / se_clogit_r)
    cat("    SEs (robust): max|rel diff| =", format(max_se_r_diff, digits = 3, scientific = TRUE), "\n")

    # Log per-coefficient detail for robust SEs
    cat("    Per-coefficient robust SE comparison:\n")
    for (nm in common) {
      cat("      ", nm, ": fast=", round(se_fast_r[nm], 6),
          " clogit=", round(se_clogit_r[nm], 6),
          " rel_diff=", format(abs(se_fast_r[nm] - se_clogit_r[nm]) / se_clogit_r[nm],
                               digits = 3, scientific = TRUE), "\n")
    }
  } else {
    # No clustering: just compare model-based SEs
    se_fast   <- fit_fast$se[common]
    se_clogit <- setNames(sqrt(diag(vcov(fit_clogit))), names(coef(fit_clogit)))[common]
    max_se_diff <- max(abs(se_fast - se_clogit) / se_clogit)
    cat("    SEs (model):  max|rel diff| =", format(max_se_diff, digits = 3, scientific = TRUE), "\n")
  }

  # Per-coefficient detail
  cat("    Per-coefficient comparison:\n")
  for (nm in common) {
    cat("      ", nm, ": fast=", round(beta_fast[nm], 6),
        " clogit=", round(beta_clogit[nm], 6),
        " diff=", format(beta_fast[nm] - beta_clogit[nm], digits = 3, scientific = TRUE), "\n")
  }

  pass <- max_coef_diff < 1e-4 && ll_diff < 1e-3
  cat("    Result:      ", if (pass) "PASS" else "*** FAIL ***", "\n\n")
  pass
}

results <- c(
  run_comparison("Basic (no offset, no cluster)",
                 n_egos = 500, n_alts = 20, use_offset = FALSE,
                 cluster_ratio = 1.0, seed = 100),
  run_comparison("With offset",
                 n_egos = 500, n_alts = 20, use_offset = TRUE,
                 cluster_ratio = 1.0, seed = 200),
  run_comparison("With offset + cluster",
                 n_egos = 500, n_alts = 20, use_offset = TRUE,
                 cluster_ratio = 1.5, seed = 300),
  run_comparison("Larger (5K egos, 100 alts)",
                 n_egos = 5000, n_alts = 100, use_offset = TRUE,
                 cluster_ratio = 1.3, seed = 400),
  run_comparison("Small groups (2 alts)",
                 n_egos = 1000, n_alts = 2, use_offset = FALSE,
                 cluster_ratio = 1.0, seed = 500)
)

cat("=== Validation Results: ", sum(results), "/", length(results), " tests passed ===\n\n")

# --- 5. Coefficient recovery test ---
cat("Step 5: Coefficient recovery from known DGP...\n")
sim <- simulate_clogit_data(n_egos = 10000, n_alts = 100,
                             use_offset = FALSE, seed = 42)

t_recovery <- system.time({
  fit <- fastclogit(sim$X, sim$choice, sim$strata, max_iter = 50, tol = 1e-10)
})
cat("  Fit time:", round(t_recovery[3], 2), "s\n")

cat("  True vs Estimated coefficients:\n")
comparison <- data.frame(
  true  = sim$beta_true,
  est   = fit$coefficients[names(sim$beta_true)],
  se    = fit$se[names(sim$beta_true)],
  z_off = (fit$coefficients[names(sim$beta_true)] - sim$beta_true) / fit$se[names(sim$beta_true)]
)
comparison$within_3se <- abs(comparison$z_off) < 3
print(round(comparison, 4))
n_matched  <- sum(!is.na(comparison$within_3se))
n_within   <- sum(comparison$within_3se, na.rm = TRUE)
n_missing  <- sum(is.na(comparison$est))
recovery_pass <- n_missing == 0 && all(comparison$within_3se, na.rm = TRUE)
cat("\n  Matched:", n_matched, "/", nrow(comparison), "coefficients")
if (n_missing > 0) cat(" (", n_missing, " name mismatches — check beta_true names)")
cat("\n  All within 3 SEs?", recovery_pass, "\n\n")

# --- 6. Memory benchmark ---
cat("Step 6: Memory benchmark...\n")

# Helper to report memory state
log_memory <- function(label) {
  gc_info <- gc(verbose = FALSE)
  cat("  [Memory] ", label, ": ",
      round(gc_info[2, 2], 1), " MB (Vcells used)\n")
}

log_memory("Before simulation")

sim_med <- simulate_clogit_data(n_egos = 50000, n_alts = 100,
                                 use_offset = TRUE, cluster_ratio = 1.3,
                                 seed = 999)
cat("  Data size: X =", format(object.size(sim_med$X), units = "MB"),
    ", total sim =", format(object.size(sim_med), units = "MB"), "\n")
cat("  Rows:", format(nrow(sim_med$X), big.mark = ","),
    " | Cols:", ncol(sim_med$X), "\n")

log_memory("After simulation, before fit")

gc(reset = TRUE)
mem_before <- gc(verbose = FALSE)[2, 2]

t_fit <- system.time({
  fit_med <- fastclogit(sim_med$X, sim_med$choice, sim_med$strata,
                         offset = sim_med$offset, cluster = sim_med$cluster,
                         max_iter = 25, verbose = TRUE)
})

mem_after <- gc(verbose = FALSE)[2, 2]
log_memory("After fit")

cat("\n  Memory delta during fit: ", round(mem_after - mem_before, 1), " MB\n")
cat("  Time: ", round(t_fit[3], 2), " seconds\n")
cat("  Converged:", fit_med$converged, " (", fit_med$iterations, " iters)\n")

# Extrapolate to MONA scale
rows_med  <- nrow(sim_med$X)
rows_mona <- 200e6
scale_factor <- rows_mona / rows_med
cat("\n  --- MONA Extrapolation ---\n")
cat("  Simulation rows:", format(rows_med, big.mark = ","), "\n")
cat("  MONA rows:      ", format(rows_mona, big.mark = ","),
    " (", round(scale_factor, 1), "x)\n")
cat("  Estimated MONA X size: ",
    round(rows_mona * ncol(sim_med$X) * 8 / 1e9, 1), " GB\n")
cat("  Estimated MONA fit memory delta: ",
    round((mem_after - mem_before) * scale_factor / 1024, 1), " GB\n")
cat("  Estimated MONA fit time: ",
    round(t_fit[3] * scale_factor, 0), " seconds\n\n")

# --- 7. Session info ---
cat("Step 7: Session info (for reproducibility)...\n")
print(sessionInfo())
cat("\n")

# --- Summary ---
cat("\n")
cat("============================================================\n")
cat("  SUMMARY\n")
cat("============================================================\n")
cat("  Validation:          ", sum(results), "/", length(results), " passed\n")
cat("  Coefficient recovery:", if (recovery_pass) "PASS" else "FAIL", "\n")
cat("  Memory benchmark:    ", round(mem_after - mem_before, 1), " MB for ",
    format(rows_med, big.mark = ","), " rows\n")
cat("  Convergence:         ", fit_med$converged, "\n")
cat("  Log file:            ", log_file, "\n")
cat("============================================================\n\n")

cat("The package is ready. To install on MONA:\n")
cat("  1. R CMD build fastclogit/\n")
cat("  2. Upload fastclogit_0.1.0.tar.gz to MONA\n")
cat("  3. install.packages('fastclogit_0.1.0.tar.gz', repos = NULL, type = 'source')\n")
cat("\n=== All done! ===\n")
