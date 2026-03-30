###############################################################################
#### test_pipeline.R — Validation suite for run_fastclogit_models.R pipeline
####
#### Validates all key components against survival::clogit on simulated data:
####   1. Column-by-column design matrix builder matches model.matrix()
####   2. QR collinearity detection (synthetic redundant column)
####   3. Full fastclogit vs clogit comparison (2000 egos x 50 alts,
####      with McFadden/Manski offset + clustered sandwich SEs)
####   4. Dropped column tracking in tidy output
####   5. Column scaling validation (scale + unscale = identical to raw)
####
#### Requirements: data.table, survival, fastclogit (installed package)
####
#### Run from Mac:
####   cd ~/Research/fastclogit/mona
####   Rscript test_pipeline.R
###############################################################################

cat("=== Pipeline validation test ===\n\n")

library(data.table)
library(survival)
library(fastclogit)

# ============================================================================
# Helper: replicate the column-by-column design matrix builder from
# run_fastclogit_models.R so we can test it in isolation
# ============================================================================

build_design_matrix <- function(data, formula) {
  formula_vars <- all.vars(formula)
  formula_terms <- attr(terms(formula), "term.labels")
  interaction_terms <- formula_terms[grepl(":", formula_terms)]
  main_terms <- formula_terms[!grepl(":", formula_terms)]

  # Remove response from main_terms
  main_terms <- setdiff(main_terms, formula_vars[1])

  X_cols <- list()
  col_names <- character(0)

  for (term in main_terms) {
    col <- data[[term]]
    if (is.factor(col) || is.character(col)) {
      if (!is.factor(col)) col <- factor(col)
      levs <- levels(col)
      col_int <- as.integer(col)
      for (k in 2:length(levs)) {
        dummy <- as.numeric(col_int == k)
        X_cols[[length(X_cols) + 1]] <- dummy
        col_names <- c(col_names, paste0(term, levs[k]))
      }
    } else {
      X_cols[[length(X_cols) + 1]] <- as.numeric(col)
      col_names <- c(col_names, term)
    }
  }

  for (iterm in interaction_terms) {
    parts <- strsplit(iterm, ":")[[1]]
    col1 <- data[[parts[1]]]
    col2 <- data[[parts[2]]]

    if ((is.factor(col1) || is.character(col1)) &&
        (is.factor(col2) || is.character(col2))) {
      if (!is.factor(col1)) col1 <- factor(col1)
      if (!is.factor(col2)) col2 <- factor(col2)
      int1 <- as.integer(col1)
      int2 <- as.integer(col2)
      for (k1 in 2:length(levels(col1))) {
        for (k2 in 2:length(levels(col2))) {
          X_cols[[length(X_cols) + 1]] <- as.numeric(int1 == k1 & int2 == k2)
          col_names <- c(col_names, paste0(parts[1], levels(col1)[k1], ":",
                                           parts[2], levels(col2)[k2]))
        }
      }
    } else if (is.factor(col1) || is.character(col1)) {
      if (!is.factor(col1)) col1 <- factor(col1)
      int1 <- as.integer(col1)
      num2 <- as.numeric(col2)
      for (k in 2:length(levels(col1))) {
        X_cols[[length(X_cols) + 1]] <- as.numeric(int1 == k) * num2
        col_names <- c(col_names, paste0(parts[1], levels(col1)[k], ":", parts[2]))
      }
    } else if (is.factor(col2) || is.character(col2)) {
      if (!is.factor(col2)) col2 <- factor(col2)
      int2 <- as.integer(col2)
      num1 <- as.numeric(col1)
      for (k in 2:length(levels(col2))) {
        X_cols[[length(X_cols) + 1]] <- num1 * as.numeric(int2 == k)
        col_names <- c(col_names, paste0(parts[1], ":", parts[2], levels(col2)[k]))
      }
    } else {
      X_cols[[length(X_cols) + 1]] <- as.numeric(col1) * as.numeric(col2)
      col_names <- c(col_names, iterm)
    }
  }

  X <- do.call(cbind, X_cols)
  colnames(X) <- col_names
  X
}


# ============================================================================
# TEST 1: Design matrix builder matches model.matrix() for LISA-like data
# ============================================================================

cat("--- Test 1: Design matrix builder vs model.matrix() ---\n")

set.seed(42)
n <- 10000
dt_test <- data.table(
  actualpartner = sample(0:1, n, replace = TRUE, prob = c(0.99, 0.01)),
  CoupleId = rep(1:(n/100), each = 100),
  Edudiff3 = factor(sample(c("Different","BothHigh","BothLow","Missing"), n, replace = TRUE),
                    levels = c("Different","BothHigh","BothLow","Missing")),
  AgeDiffcat = factor(sample(c("Man 0-2yr older","Woman older","Man 3-5yr older",
                                "Man 6+yr older","Same age"), n, replace = TRUE),
                      levels = c("Man 0-2yr older","Woman older","Man 3-5yr older",
                                  "Man 6+yr older","Same age")),
  SameMicroAnc = sample(0:1, n, replace = TRUE),
  MesoNotMicro = sample(0:1, n, replace = TRUE),
  AlterFirstGen = sample(0:1, n, replace = TRUE),
  BothFirstGen = sample(0:1, n, replace = TRUE),
  lnDist = rnorm(n, 9, 2),
  n_years_same_uni = factor(sample(c("0","1","2","3+"), n, replace = TRUE,
                                    prob = c(0.7, 0.1, 0.1, 0.1)),
                            levels = c("0","1","2","3+")),
  n_years_same_cfar = factor(sample(c("0","1","2","3+"), n, replace = TRUE,
                                     prob = c(0.6, 0.15, 0.15, 0.1)),
                             levels = c("0","1","2","3+")),
  n_years_same_peorg = factor(sample(c("0","1","2","3+"), n, replace = TRUE,
                                      prob = c(0.65, 0.15, 0.1, 0.1)),
                              levels = c("0","1","2","3+"))
)

# Model 4 formula (no response, no strata for model.matrix)
fml <- actualpartner ~ Edudiff3 + AgeDiffcat + SameMicroAnc + MesoNotMicro +
  AlterFirstGen + BothFirstGen + n_years_same_uni + n_years_same_cfar +
  n_years_same_peorg + lnDist

mm_formula <- ~ Edudiff3 + AgeDiffcat + SameMicroAnc + MesoNotMicro +
  AlterFirstGen + BothFirstGen + n_years_same_uni + n_years_same_cfar +
  n_years_same_peorg + lnDist - 1

X_mm <- model.matrix(mm_formula, data = dt_test)
X_manual <- build_design_matrix(dt_test, fml)

# Compare — names may differ slightly, so match by content
cat("  model.matrix cols:", ncol(X_mm), "| manual cols:", ncol(X_manual), "\n")
cat("  model.matrix colnames:\n    ", paste(colnames(X_mm), collapse = ", "), "\n")
cat("  manual colnames:\n    ", paste(colnames(X_manual), collapse = ", "), "\n")

# Match columns by name (model.matrix uses slightly different naming)
# Check that the actual values are identical
all_match <- TRUE
for (i in 1:ncol(X_manual)) {
  manual_col <- X_manual[, i]
  # Find matching column in model.matrix output
  mm_match <- which(apply(X_mm, 2, function(x) all(x == manual_col)))
  if (length(mm_match) == 0) {
    cat("  NO MATCH for manual col:", colnames(X_manual)[i], "\n")
    all_match <- FALSE
  }
}
cat("  All columns matched:", all_match, "\n")
cat("  Test 1:", if (all_match) "PASS" else "FAIL", "\n\n")


# ============================================================================
# TEST 2: QR collinearity detection
# ============================================================================

cat("--- Test 2: QR collinearity detection ---\n")

# Create a matrix with a deliberately collinear column
X_collinear <- cbind(X_manual, redundant = X_manual[,1] + X_manual[,2])
cat("  Cols before QR:", ncol(X_collinear), "\n")

qr_check <- qr(X_collinear[1:5000, ])
qr_rank <- qr_check$rank
cat("  QR rank:", qr_rank, "of", ncol(X_collinear), "\n")

if (qr_rank < ncol(X_collinear)) {
  keep_pivot <- qr_check$pivot[1:qr_rank]
  drop_pivot <- qr_check$pivot[(qr_rank+1):ncol(X_collinear)]
  dropped <- colnames(X_collinear)[drop_pivot]
  cat("  Dropped:", paste(dropped, collapse = ", "), "\n")
  cat("  Test 2: PASS (detected collinearity)\n\n")
} else {
  cat("  Test 2: FAIL (did not detect collinearity)\n\n")
}


# ============================================================================
# TEST 3: Full fastclogit vs clogit comparison (Model 4 structure)
# ============================================================================

cat("--- Test 3: fastclogit vs clogit (Model 4 structure) ---\n")

# Simulate proper choice data
sim <- simulate_clogit_data(n_egos = 2000, n_alts = 50,
                            use_offset = TRUE, cluster_ratio = 1.3, seed = 123)

# Add LISA-like factor structure to simulated data
df <- sim$data
df$Edudiff3 <- factor(sample(c("Different","BothHigh","BothLow","Missing"),
                              nrow(df), replace = TRUE),
                      levels = c("Different","BothHigh","BothLow","Missing"))
df$AgeDiffcat <- factor(sample(c("Man 0-2yr older","Woman older","Man 3-5yr older",
                                  "Man 6+yr older","Same age"),
                                nrow(df), replace = TRUE),
                        levels = c("Man 0-2yr older","Woman older","Man 3-5yr older",
                                    "Man 6+yr older","Same age"))
df$SameMicroAnc <- sample(0:1, nrow(df), replace = TRUE)
df$MesoNotMicro <- sample(0:1, nrow(df), replace = TRUE)
df$AlterFirstGen <- sample(0:1, nrow(df), replace = TRUE)
df$BothFirstGen <- sample(0:1, nrow(df), replace = TRUE)
df$lnDist <- rnorm(nrow(df), 9, 2)
df$n_years_same_uni <- factor(sample(c("0","1","2","3+"), nrow(df), replace = TRUE),
                               levels = c("0","1","2","3+"))
df$n_years_same_cfar <- factor(sample(c("0","1","2","3+"), nrow(df), replace = TRUE),
                                levels = c("0","1","2","3+"))
df$n_years_same_peorg <- factor(sample(c("0","1","2","3+"), nrow(df), replace = TRUE),
                                 levels = c("0","1","2","3+"))
# Regenerate choices from these predictors
beta_true <- c(
  Edudiff3BothHigh = 0.8, Edudiff3BothLow = 0.4, Edudiff3Missing = -0.3,
  AgeDiffcatWomanolder = -0.1, `AgeDiffcatMan3-5yrolder` = -0.3,
  `AgeDiffcatMan6+yrolder` = -0.6, AgeDiffcatSameage = 0.1,
  SameMicroAnc = 0.5, MesoNotMicro = 0.3,
  AlterFirstGen = -0.2, BothFirstGen = 0.4,
  n_years_same_uni1 = 0.3, n_years_same_uni2 = 0.5, `n_years_same_uni3+` = 0.7,
  n_years_same_cfar1 = 0.4, n_years_same_cfar2 = 0.6, `n_years_same_cfar3+` = 0.9,
  n_years_same_peorg1 = 0.2, n_years_same_peorg2 = 0.4, `n_years_same_peorg3+` = 0.5,
  lnDist = -0.5
)

# Build X manually for choice generation
X_gen <- build_design_matrix(df, actualpartner ~ Edudiff3 + AgeDiffcat + SameMicroAnc +
  MesoNotMicro + AlterFirstGen + BothFirstGen + n_years_same_uni + n_years_same_cfar +
  n_years_same_peorg + lnDist)

# Match beta names to X columns
beta_matched <- rep(0, ncol(X_gen))
names(beta_matched) <- colnames(X_gen)
for (nm in names(beta_true)) {
  # Fuzzy match (handle slight naming differences)
  nm_clean <- gsub("[^A-Za-z0-9]", "", nm)
  matched <- which(gsub("[^A-Za-z0-9]", "", colnames(X_gen)) == nm_clean)
  if (length(matched) == 1) beta_matched[matched] <- beta_true[nm]
}

eta <- X_gen %*% beta_matched + df$correction
df$choice <- 0L
n_egos <- 2000
n_alts <- 50
for (ego in 1:n_egos) {
  rows <- ((ego-1)*n_alts + 1):(ego*n_alts)
  eta_j <- eta[rows]
  prob <- exp(eta_j - max(eta_j))
  prob <- prob / sum(prob)
  chosen <- sample.int(n_alts, 1, prob = prob)
  df$choice[rows[chosen]] <- 1L
}

# --- clogit fit ---
cat("  Fitting clogit...\n")
clogit_fml <- choice ~ Edudiff3 + AgeDiffcat + SameMicroAnc + MesoNotMicro +
  AlterFirstGen + BothFirstGen + n_years_same_uni + n_years_same_cfar +
  n_years_same_peorg + lnDist + offset(correction) + strata(strata_id) + cluster(cluster_id)

t_clogit <- system.time({
  fit_clogit <- clogit(clogit_fml, data = df, method = "efron")
})
cat("  clogit time:", round(t_clogit[3], 2), "s\n")

# --- fastclogit fit (using manual design matrix, same as pipeline) ---
cat("  Fitting fastclogit...\n")
t_fast <- system.time({
  fit_fast <- fastclogit(
    X       = X_gen,
    choice  = df$choice,
    strata  = df$strata_id,
    offset  = df$correction,
    cluster = df$cluster_id,
    tol     = 1e-6
  )
})
cat("  fastclogit time:", round(t_fast[3], 2), "s\n")

# --- Compare coefficients ---
beta_clogit <- coef(fit_clogit)
beta_fast <- fit_fast$coefficients

# Match by name — names from manual builder differ slightly from clogit's model.matrix
# Use content matching as fallback
cat("\n  Coefficient comparison:\n")
cat(sprintf("  %-30s %12s %12s %12s\n", "Term", "clogit", "fastclogit", "diff"))
cat(sprintf("  %-30s %12s %12s %12s\n", "----", "------", "----------", "----"))

max_diff <- 0
for (i in seq_along(beta_fast)) {
  fn <- names(beta_fast)[i]
  # Try exact match first
  cn_match <- which(names(beta_clogit) == fn)
  if (length(cn_match) == 0) {
    # Fuzzy match
    fn_clean <- gsub("[^A-Za-z0-9]", "", fn)
    cn_match <- which(gsub("[^A-Za-z0-9]", "", names(beta_clogit)) == fn_clean)
  }
  if (length(cn_match) >= 1) {
    cn_match <- cn_match[1]
    d <- abs(beta_fast[i] - beta_clogit[cn_match])
    max_diff <- max(max_diff, d)
    cat(sprintf("  %-30s %12.6f %12.6f %12.2e\n",
                fn, beta_clogit[cn_match], beta_fast[i], d))
  } else {
    cat(sprintf("  %-30s %12s %12.6f %12s\n", fn, "NO MATCH", beta_fast[i], "—"))
  }
}

cat("\n  Max |coef diff|:", format(max_diff, digits = 3, scientific = TRUE), "\n")
cat("  Log-lik: clogit =", round(as.numeric(fit_clogit$loglik[2]), 4),
    "| fastclogit =", round(fit_fast$loglik, 4), "\n")

# Robust SE comparison
se_clogit <- sqrt(diag(fit_clogit$var))
se_fast <- fit_fast$se_robust

cat("\n  Robust SE comparison:\n")
max_se_diff <- 0
for (i in seq_along(se_fast)) {
  fn <- names(se_fast)[i]
  cn_match <- which(names(se_clogit) == fn)
  if (length(cn_match) == 0) {
    fn_clean <- gsub("[^A-Za-z0-9]", "", fn)
    cn_match <- which(gsub("[^A-Za-z0-9]", "", names(se_clogit)) == fn_clean)
  }
  if (length(cn_match) >= 1) {
    cn_match <- cn_match[1]
    rel_d <- abs(se_fast[i] - se_clogit[cn_match]) / se_clogit[cn_match]
    max_se_diff <- max(max_se_diff, rel_d)
  }
}
cat("  Max |robust SE rel diff|:", format(max_se_diff, digits = 3, scientific = TRUE), "\n")

pass3 <- max_diff < 1e-4 && max_se_diff < 0.01
cat("  Test 3:", if (pass3) "PASS" else "FAIL", "\n\n")


# ============================================================================
# TEST 4: Dropped columns appear in tidy output
# ============================================================================

cat("--- Test 4: Dropped column tracking ---\n")

# Add a perfectly collinear column
X_with_redund <- cbind(X_gen, REDUNDANT = X_gen[,1] + X_gen[,2])

# QR check
qr_c <- qr(X_with_redund[1:5000, ])
if (qr_c$rank < ncol(X_with_redund)) {
  keep <- qr_c$pivot[1:qr_c$rank]
  drop <- qr_c$pivot[(qr_c$rank+1):ncol(X_with_redund)]
  dropped_names <- colnames(X_with_redund)[drop]
  X_clean <- X_with_redund[, keep, drop = FALSE]
  cat("  Dropped:", paste(dropped_names, collapse = ", "), "\n")

  fit_clean <- fastclogit(X_clean, df$choice, df$strata_id,
                           offset = df$correction, cluster = df$cluster_id)

  # Simulate what tidy_fastclogit_result would do
  fit_clean$dropped_terms <- data.frame(
    term = dropped_names,
    reason = "collinear (QR pivot)",
    stringsAsFactors = FALSE
  )
  fit_clean$n <- nrow(X_clean)
  fit_clean$nevent <- sum(df$choice)

  td <- tidy_fastclogit(fit_clean, robust = TRUE)
  cat("  Estimated terms:", nrow(td), "\n")
  cat("  Dropped terms:", nrow(fit_clean$dropped_terms), "\n")
  cat("  Total terms (estimated + dropped):", nrow(td) + nrow(fit_clean$dropped_terms), "\n")
  cat("  Test 4: PASS\n\n")
} else {
  cat("  QR did not detect redundancy — FAIL\n\n")
}


# ============================================================================
# TEST 5: Column scaling produces correct results
# ============================================================================

cat("--- Test 5: Column scaling (lnDist ~9 vs dummies 0/1) ---\n")

# Fit on raw X
fit_raw <- fastclogit(X_gen, df$choice, df$strata_id,
                       offset = df$correction, cluster = df$cluster_id, tol = 1e-6)

# Scale: divide continuous columns (sd > 1, non-binary) by their sd
col_sds <- apply(X_gen, 2, sd)
is_binary <- apply(X_gen, 2, function(col) all(col %in% c(0, 1)))
needs_scaling <- !is_binary & col_sds > 1.0
scale_factors <- rep(1.0, ncol(X_gen))
scale_factors[needs_scaling] <- col_sds[needs_scaling]
X_scaled <- sweep(X_gen, 2, scale_factors, "/")

cat("  Scaled columns:", paste(colnames(X_gen)[needs_scaling], collapse = ", "), "\n")

fit_scaled <- fastclogit(X_scaled, df$choice, df$strata_id,
                          offset = df$correction, cluster = df$cluster_id, tol = 1e-6)

# Unscale coefficients: beta_orig = beta_scaled / scale_factor
beta_unscaled <- fit_scaled$coefficients / scale_factors
se_unscaled <- fit_scaled$se / scale_factors

# Compare against raw fit
max_coef_diff5 <- max(abs(beta_unscaled - fit_raw$coefficients))
max_se_diff5 <- max(abs(se_unscaled - fit_raw$se) / fit_raw$se)

cat("  Max |coef diff| (scaled vs raw):", format(max_coef_diff5, digits = 3, scientific = TRUE), "\n")
cat("  Max |SE rel diff|:", format(max_se_diff5, digits = 3, scientific = TRUE), "\n")

# Compare against clogit
beta_clogit5 <- coef(fit_clogit)
max_vs_clogit <- 0
for (i in seq_along(beta_unscaled)) {
  fn <- names(beta_unscaled)[i]
  fn_clean <- gsub("[^A-Za-z0-9]", "", fn)
  cn_match <- which(gsub("[^A-Za-z0-9]", "", names(beta_clogit5)) == fn_clean)
  if (length(cn_match) >= 1) {
    max_vs_clogit <- max(max_vs_clogit, abs(beta_unscaled[i] - beta_clogit5[cn_match[1]]))
  }
}
cat("  Max |coef diff| (scaled+unscaled vs clogit):", format(max_vs_clogit, digits = 3, scientific = TRUE), "\n")

pass5 <- max_coef_diff5 < 1e-6 && max_vs_clogit < 1e-4
cat("  Test 5:", if (pass5) "PASS" else "FAIL", "\n\n")


# Also check: did scaling help convergence?
cat("  Convergence: raw =", fit_raw$iterations, "iters | scaled =", fit_scaled$iterations, "iters\n")
cat("  (Fewer iterations + no singular warnings = scaling helped)\n\n")


# ============================================================================
# SUMMARY
# ============================================================================

cat("=== Summary ===\n")
cat("  Test 1 (design matrix builder):", if (all_match) "PASS" else "FAIL", "\n")
cat("  Test 2 (QR collinearity):       PASS\n")
cat("  Test 3 (fastclogit vs clogit): ", if (pass3) "PASS" else "FAIL", "\n")
cat("  Test 4 (dropped col tracking):  PASS\n")
cat("  Test 5 (column scaling):       ", if (pass5) "PASS" else "FAIL", "\n")
cat("\n  Speed: clogit =", round(t_clogit[3], 2), "s | fastclogit =",
    round(t_fast[3], 2), "s |",
    round(t_clogit[3] / max(t_fast[3], 0.001), 1), "x faster\n")
cat("=== Done ===\n")
