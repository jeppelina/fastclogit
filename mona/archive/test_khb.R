###############################################################################
#### test_khb.R — Validation suite for khb_fastclogit.R
####
#### Tests the KHB decomposition end-to-end on simulated data:
####   1. Basic plumbing: khb_clogit() runs, returns correct structure
####   2. Known mediation: synthetic Z that mediates X→Y, verify indirect > 0
####   3. No mediation: Z is noise, verify indirect ≈ 0
####   4. Factor key variable: factor X with multiple levels decomposes
####   5. Coefficient names survive the full pipeline (names bug regression)
####   6. Offset auto-detect: verifies the global-check logic
####
#### Requirements: data.table, fastclogit (installed package)
####
#### Run from Mac:
####   cd ~/Research/fastclogit/mona
####   Rscript test_khb.R
###############################################################################

cat("=== KHB validation test ===\n\n")

library(data.table)
library(fastclogit)

# Source the KHB script (lives next to the pipeline scripts)
khb_path <- file.path(dirname(getwd()), "Scripts feb 11", "khb_fastclogit.R")
if (!file.exists(khb_path)) {
  # Try relative to this script's location
  khb_path <- file.path("..", "..", "Scripts feb 11", "khb_fastclogit.R")
}
if (!file.exists(khb_path)) {
  # Last resort: absolute path on Mac
  khb_path <- "~/Research/Scripts feb 11/khb_fastclogit.R"
}
if (!file.exists(khb_path)) {
  stop("Cannot find khb_fastclogit.R. Tried multiple paths.\n",
       "  Run this from ~/Research/fastclogit/mona/")
}
cat("  Sourcing:", khb_path, "\n\n")
source(khb_path)


# ============================================================================
# Helper: generate KHB test data
#
# Creates a data.table with choice-set structure mimicking the pipeline:
#   - CoupleId (strata), LopNrEgo (cluster), actualpartner (outcome)
#   - X: key variable(s), Z: mediator(s), C: controls
#   - Optional: correction column for McFadden/Manski offset
#
# DGP:
#   Z = gamma * X + noise           (mediation path: X -> Z)
#   P(chosen) ∝ exp(beta_x * X + beta_z * Z + beta_c * C + offset)
#
# If beta_z != 0 and gamma != 0, there is genuine mediation.
# ============================================================================

generate_khb_data <- function(n_egos = 500, n_alts = 30,
                               beta_x = 0.5, beta_z = 0.8,
                               beta_c = 0.3, gamma_xz = 0.6,
                               use_offset = TRUE, seed = 42) {
  set.seed(seed)
  n_total <- n_egos * n_alts

  # Strata and cluster structure
  CoupleId <- rep(seq_len(n_egos), each = n_alts)
  # Some clustering: ~70% unique clusters
  n_clusters <- max(1L, round(n_egos * 0.7))
  ego_cluster <- sample(seq_len(n_clusters), n_egos, replace = TRUE)
  LopNrEgo <- rep(ego_cluster, each = n_alts)

  # Predictors
  X_var <- rnorm(n_total, 0, 1)          # key variable (continuous)
  C_var <- rnorm(n_total, 0, 1)          # control
  Z_var <- gamma_xz * X_var + rnorm(n_total, 0, 1)  # mediator (correlated with X)
  Z_noise <- rnorm(n_total, 0, 1)        # pure noise (for no-mediation test)

  # Offset (McFadden/Manski correction)
  correction <- rep(0, n_total)
  if (use_offset) {
    strata_props <- c(0.1, 0.2, 0.3, 0.4)
    pop_sizes <- c(50, 200, 500, 5000)
    for (ego in seq_len(n_egos)) {
      rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
      stratum <- sample(seq_along(strata_props), n_alts, replace = TRUE,
                         prob = strata_props)
      for (s in seq_along(strata_props)) {
        s_rows <- rows[stratum == s]
        if (length(s_rows) > 0) {
          correction[s_rows] <- -log(length(s_rows) / pop_sizes[s])
        }
      }
    }
  }

  # Generate choices from the true model: Y ~ X + Z + C + offset
  eta <- beta_x * X_var + beta_z * Z_var + beta_c * C_var + correction
  actualpartner <- integer(n_total)
  for (ego in seq_len(n_egos)) {
    rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
    eta_j <- eta[rows]
    prob <- exp(eta_j - max(eta_j))
    prob <- prob / sum(prob)
    chosen <- sample.int(n_alts, 1, prob = prob)
    actualpartner[rows[chosen]] <- 1L
  }

  dt <- data.table(
    actualpartner = actualpartner,
    CoupleId = CoupleId,
    LopNrEgo = LopNrEgo,
    X_key = X_var,
    Z_mediator = Z_var,
    Z_noise = Z_noise,
    C_control = C_var,
    correction = correction
  )

  list(
    dt = dt,
    beta_true = c(X_key = beta_x, Z_mediator = beta_z, C_control = beta_c),
    gamma_xz = gamma_xz
  )
}


# ============================================================================
# TEST 1: Basic plumbing — khb_clogit runs and returns expected structure
# ============================================================================

cat("--- Test 1: Basic plumbing ---\n")

sim <- generate_khb_data(n_egos = 500, n_alts = 30, seed = 1)

# Set global for offset auto-detect
stratified_sampling <- TRUE

result <- tryCatch({
  khb_clogit(
    data = sim$dt,
    key_vars = "X_key",
    z_vars = "Z_mediator",
    concomitant = "C_control",
    group_var = "CoupleId",
    cluster_var = "LopNrEgo",
    outcome_var = "actualpartner",
    use_offset = NULL,   # test auto-detect
    verbose = TRUE
  )
}, error = function(e) {
  cat("\n  !!! ERROR:", e$message, "\n")
  cat("  !!! Call: "); print(e$call)
  return(list(success = FALSE, error = e$message))
})

# Check return structure
expected_fields <- c("success", "decomposition", "z_effects",
                      "coef_reduced", "coef_full",
                      "n_obs", "n_events", "key_vars", "z_vars", "concomitant")
has_fields <- all(expected_fields %in% names(result))
cat("  success:", result$success, "\n")
cat("  All expected fields present:", has_fields, "\n")

# Diagnostic: print coef matrix dimensions and names
if ("coef_reduced" %in% names(result)) {
  cat("  coef_reduced dim:", paste(dim(result$coef_reduced), collapse="x"), "\n")
  cat("  coef_reduced rownames:", paste(rownames(result$coef_reduced), collapse=", "), "\n")
  cat("  coef_reduced colnames:", paste(colnames(result$coef_reduced), collapse=", "), "\n")
}
if ("coef_full" %in% names(result)) {
  cat("  coef_full dim:", paste(dim(result$coef_full), collapse="x"), "\n")
  cat("  coef_full rownames:", paste(rownames(result$coef_full), collapse=", "), "\n")
}

if (!is.null(result$decomposition)) {
  cat("  Decomposition rows:", nrow(result$decomposition), "\n")
} else {
  cat("  Decomposition: NULL\n")
}
if (!is.null(result$z_effects)) {
  cat("  Z-effects rows:", nrow(result$z_effects), "\n")
} else {
  cat("  Z-effects: NULL\n")
}

pass1 <- isTRUE(result$success) && has_fields &&
         !is.null(result$decomposition) && nrow(result$decomposition) > 0
cat("  Test 1:", if (pass1) "PASS" else "FAIL", "\n\n")


# ============================================================================
# TEST 2: Known mediation — indirect effect should be meaningfully > 0
#
# DGP: Z = 0.6*X + noise, Y ~ 0.5*X + 0.8*Z + 0.3*C
# Total effect of X ≈ 0.5 + 0.8*0.6 = 0.98  (via Z path)
# Direct effect of X ≈ 0.5
# Indirect effect ≈ 0.48
#
# With 500 egos we won't nail exact values, but the sign/magnitude
# should be clearly positive and nontrivial.
# ============================================================================

cat("--- Test 2: Known mediation (indirect > 0) ---\n")

if (!pass1) {
  cat("  SKIPPED (Test 1 failed)\n\n")
  pass2 <- FALSE
} else {
  decomp <- result$decomposition

  cat("  Total effect (reduced model):  ", round(decomp$total_effect[1], 4), "\n")
  cat("  Direct effect (full model):    ", round(decomp$direct_effect[1], 4), "\n")
  cat("  Indirect effect (difference):  ", round(decomp$indirect_effect[1], 4), "\n")
  cat("  Confounding %:                 ", round(decomp$conf_pct[1], 1), "%\n")

  # Indirect should be positive (X→Z path has same-sign gamma and beta_z)
  indirect_positive <- decomp$indirect_effect[1] > 0
  # Indirect should be at least 10% of total (our true value is ~49%)
  indirect_nontrivial <- abs(decomp$conf_pct[1]) > 10
  # Total should be larger than direct (reduced model absorbs Z-mediation)
  total_gt_direct <- abs(decomp$total_effect[1]) > abs(decomp$direct_effect[1])

  pass2 <- indirect_positive && indirect_nontrivial && total_gt_direct
  cat("  Indirect positive:", indirect_positive, "\n")
  cat("  Indirect nontrivial (>10%):", indirect_nontrivial, "\n")
  cat("  |Total| > |Direct|:", total_gt_direct, "\n")
  cat("  Test 2:", if (pass2) "PASS" else "FAIL", "\n\n")
}


# ============================================================================
# TEST 3: No mediation — Z is pure noise, indirect ≈ 0
#
# Key: same X, but Z_noise is independent of X. The reduced model with
# Z_noise_resid should give ≈ same coefficient as full model with Z_noise.
# ============================================================================

cat("--- Test 3: No mediation (Z = noise) ---\n")

pass3 <- tryCatch({
  result_null <- khb_clogit(
    data = sim$dt,
    key_vars = "X_key",
    z_vars = "Z_noise",
    concomitant = "C_control",
    group_var = "CoupleId",
    cluster_var = "LopNrEgo",
    outcome_var = "actualpartner",
    use_offset = TRUE,
    verbose = FALSE
  )

  decomp_null <- result_null$decomposition

  cat("  Total effect:   ", round(decomp_null$total_effect[1], 4), "\n")
  cat("  Direct effect:  ", round(decomp_null$direct_effect[1], 4), "\n")
  cat("  Indirect effect:", round(decomp_null$indirect_effect[1], 4), "\n")
  cat("  Confounding %:  ", round(decomp_null$conf_pct[1], 1), "%\n")

  indirect_small <- abs(decomp_null$indirect_effect[1]) < 0.15
  conf_small <- abs(decomp_null$conf_pct[1]) < 15

  ok <- result_null$success && indirect_small && conf_small
  cat("  |Indirect| < 0.15:", indirect_small, "\n")
  cat("  |Conf %| < 15%:", conf_small, "\n")
  cat("  Test 3:", if (ok) "PASS" else "FAIL", "\n\n")
  ok
}, error = function(e) {
  cat("  !!! ERROR:", e$message, "\n")
  cat("  Test 3: FAIL (error)\n\n")
  FALSE
})


# ============================================================================
# TEST 4: Factor key variable — verify multi-level decomposition
#
# Creates a 3-level factor X (ref, level2, level3) and decomposes each.
# ============================================================================

cat("--- Test 4: Factor key variable ---\n")

set.seed(77)
n_egos_f <- 600
n_alts_f <- 25
n_total_f <- n_egos_f * n_alts_f

dt_fac <- data.table(
  CoupleId = rep(seq_len(n_egos_f), each = n_alts_f),
  LopNrEgo = rep(sample(seq_len(round(n_egos_f * 0.7)), n_egos_f, replace = TRUE),
                  each = n_alts_f),
  X_factor = factor(sample(c("ref", "level2", "level3"), n_total_f, replace = TRUE),
                     levels = c("ref", "level2", "level3")),
  Z_med = rnorm(n_total_f),
  C_ctrl = rnorm(n_total_f)
)

# Make Z correlated with factor levels
dt_fac[X_factor == "level2", Z_med := Z_med + 0.5]
dt_fac[X_factor == "level3", Z_med := Z_med + 1.0]

# Generate choices
beta_f <- c(level2 = 0.4, level3 = 0.8, Z = 0.6, C = 0.2)
dt_fac[, eta := (as.integer(X_factor == "level2") * beta_f["level2"] +
                  as.integer(X_factor == "level3") * beta_f["level3"] +
                  Z_med * beta_f["Z"] + C_ctrl * beta_f["C"])]
dt_fac[, actualpartner := 0L]
for (ego in seq_len(n_egos_f)) {
  rows <- ((ego - 1) * n_alts_f + 1):(ego * n_alts_f)
  eta_j <- dt_fac$eta[rows]
  prob <- exp(eta_j - max(eta_j))
  prob <- prob / sum(prob)
  chosen <- sample.int(n_alts_f, 1, prob = prob)
  dt_fac[rows[chosen], actualpartner := 1L]
}
dt_fac[, eta := NULL]

# No offset for this test
stratified_sampling <- FALSE

pass4 <- tryCatch({
  result_fac <- khb_clogit(
    data = dt_fac,
    key_vars = "X_factor",
    z_vars = "Z_med",
    concomitant = "C_ctrl",
    group_var = "CoupleId",
    cluster_var = "LopNrEgo",
    outcome_var = "actualpartner",
    use_offset = FALSE,
    verbose = FALSE
  )

  cat("  Success:", result_fac$success, "\n")
  cat("  Decomposition rows:", nrow(result_fac$decomposition), "\n")

  if (result_fac$success) {
    for (i in 1:nrow(result_fac$decomposition)) {
      r <- result_fac$decomposition[i]
      cat(sprintf("  %s: total=%.4f, direct=%.4f, indirect=%.4f (%.1f%%)\n",
                  r$coefficient, r$total_effect, r$direct_effect,
                  r$indirect_effect, r$conf_pct))
    }
  }

  has_two_levels <- nrow(result_fac$decomposition) == 2
  both_indirect_positive <- all(result_fac$decomposition$indirect_effect > 0)

  ok <- result_fac$success && has_two_levels && both_indirect_positive
  cat("  Two decomposition levels:", has_two_levels, "\n")
  cat("  Both indirects positive:", both_indirect_positive, "\n")
  cat("  Test 4:", if (ok) "PASS" else "FAIL", "\n\n")
  ok
}, error = function(e) {
  cat("  !!! ERROR:", e$message, "\n")
  cat("  Test 4: FAIL (error)\n\n")
  FALSE
})


# ============================================================================
# TEST 5: Coefficient names survive (regression test for the names bug)
#
# The bug: after unscaling, names(fit$coefficients) was NULL, causing
# tidy_fastclogit() to crash. Verify names are intact throughout.
# ============================================================================

cat("--- Test 5: Coefficient names survive pipeline ---\n")

if (!pass1) {
  cat("  SKIPPED (Test 1 failed)\n\n")
  pass5 <- FALSE
} else {
  # Re-run test 1 result (already computed)
  coef_r <- result$coef_reduced
  coef_f <- result$coef_full

  # Check that rownames exist on both coefficient matrices
  has_rownames_r <- !is.null(rownames(coef_r)) && length(rownames(coef_r)) > 0
  has_rownames_f <- !is.null(rownames(coef_f)) && length(rownames(coef_f)) > 0

  # Check that decomposition has non-empty coefficient names
  has_coef_names <- all(nchar(result$decomposition$coefficient) > 0)

  # Check that z_effects has names
  has_z_names <- all(nchar(result$z_effects$coefficient) > 0)

  cat("  Reduced model has rownames:", has_rownames_r, "\n")
  cat("  Full model has rownames:", has_rownames_f, "\n")
  cat("  Decomposition has coefficient names:", has_coef_names, "\n")
  cat("  Z-effects has coefficient names:", has_z_names, "\n")

  if (has_rownames_r) cat("  Reduced rownames: ", paste(rownames(coef_r), collapse = ", "), "\n")
  if (has_rownames_f) cat("  Full rownames:    ", paste(rownames(coef_f), collapse = ", "), "\n")

  pass5 <- has_rownames_r && has_rownames_f && has_coef_names && has_z_names
  cat("  Test 5:", if (pass5) "PASS" else "FAIL", "\n\n")
}


# ============================================================================
# TEST 6: Offset auto-detect logic
#
# Tests the auto-detect code path: use_offset = NULL should detect the
# global stratified_sampling variable and the correction column.
# ============================================================================

cat("--- Test 6: Offset auto-detect ---\n")

pass6 <- tryCatch({
  # Case A: stratified_sampling = TRUE + correction column exists → should use offset
  stratified_sampling <<- TRUE
  result_auto <- khb_clogit(
    data = sim$dt,
    key_vars = "X_key",
    z_vars = "Z_mediator",
    concomitant = "C_control",
    use_offset = NULL,
    verbose = FALSE
  )
  auto_on <- result_auto$success
  cat("  Case A (global=TRUE, column=present): success =", auto_on, "\n")

  # Case B: stratified_sampling = FALSE → should NOT use offset
  stratified_sampling <<- FALSE
  result_no_auto <- khb_clogit(
    data = sim$dt,
    key_vars = "X_key",
    z_vars = "Z_mediator",
    concomitant = "C_control",
    use_offset = NULL,
    verbose = FALSE
  )
  auto_off <- result_no_auto$success
  cat("  Case B (global=FALSE): success =", auto_off, "\n")

  # Case C: no correction column → should NOT use offset regardless
  dt_no_corr <- copy(sim$dt)
  dt_no_corr[, correction := NULL]
  stratified_sampling <<- TRUE
  result_no_col <- khb_clogit(
    data = dt_no_corr,
    key_vars = "X_key",
    z_vars = "Z_mediator",
    concomitant = "C_control",
    use_offset = NULL,
    verbose = FALSE
  )
  auto_no_col <- result_no_col$success
  cat("  Case C (global=TRUE, no column): success =", auto_no_col, "\n")

  # Results should be slightly different with/without offset (coefficients differ)
  coef_with <- result_auto$decomposition$total_effect[1]
  coef_without <- result_no_auto$decomposition$total_effect[1]
  offsets_differ <- abs(coef_with - coef_without) > 0.001
  cat("  Offset changes results:", offsets_differ,
      " (diff =", round(abs(coef_with - coef_without), 4), ")\n")

  ok <- auto_on && auto_off && auto_no_col
  cat("  Test 6:", if (ok) "PASS" else "FAIL", "\n\n")
  ok
}, error = function(e) {
  cat("  !!! ERROR:", e$message, "\n")
  cat("  Test 6: FAIL (error)\n\n")
  FALSE
})

# Clean up global
if (exists("stratified_sampling")) rm(stratified_sampling)


# ============================================================================
# SUMMARY
# ============================================================================

cat("=== Summary ===\n")
cat("  Test 1 (basic plumbing):       ", if (pass1) "PASS" else "FAIL", "\n")
cat("  Test 2 (known mediation):      ", if (pass2) "PASS" else "FAIL", "\n")
cat("  Test 3 (no mediation):         ", if (pass3) "PASS" else "FAIL", "\n")
cat("  Test 4 (factor key variable):  ", if (pass4) "PASS" else "FAIL", "\n")
cat("  Test 5 (coefficient names):    ", if (pass5) "PASS" else "FAIL", "\n")
cat("  Test 6 (offset auto-detect):   ", if (pass6) "PASS" else "FAIL", "\n")

all_pass <- all(pass1, pass2, pass3, pass4, pass5, pass6)
cat("\n  Overall:", if (all_pass) "ALL PASS" else "SOME FAILED", "\n")
cat("=== Done ===\n")
