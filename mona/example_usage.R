###############################################################################
#### example_usage.R — Quick demo of fastclogit on simulated data
####
#### Run this after source("load_fastclogit.R") to verify everything works
#### and to see the main features in action.
###############################################################################

cat("\n=== fastclogit example ===\n\n")

# --- 1. Simulate data ---
cat("1. Simulating data (500 egos x 30 alternatives)...\n")
sim <- simulate_clogit_data(n_egos = 500, n_alts = 30, seed = 42)
d <- sim$data
cat("   ", nrow(d), " rows, true beta:\n")
print(round(sim$beta_true, 3))

# --- 2. Formula interface ---
cat("\n2. Fitting model with formula interface...\n")
fit <- fclogit(
  choice ~ lnDist + n_years_same_cfar + n_years_same_peorg + n_years_same_uni,
  data    = d,
  strata  = "strata_id",
  cluster = "cluster_id",
  offset  = "correction"
)
cat("\n")
summary(fit)

# --- 3. Tidy output ---
cat("\n3. Tidy output (broom-style data.frame):\n")
print(tidy_fastclogit(fit))

# --- 4. Odds ratios ---
cat("\n4. Odds ratios:\n")
print(tidy_fastclogit(fit, exponentiate = TRUE))

# --- 5. KHB decomposition (if loaded) ---
if (exists("khb_decompose")) {
  cat("\n5. KHB mediation decomposition...\n")

  # Create data with known mediation: z is caused by x
  set.seed(123)
  n_egos <- 300; n_alts <- 20; n <- n_egos * n_alts
  khb_data <- data.frame(
    strata_id = rep(1:n_egos, each = n_alts),
    x  = rnorm(n),
    c1 = rnorm(n)
  )
  khb_data$z <- 0.6 * khb_data$x + rnorm(n)  # mediated path
  eta <- 0.5 * khb_data$x + 0.8 * khb_data$z + 0.3 * khb_data$c1

  # Generate choices
  khb_data$choice <- 0L
  for (i in 1:n_egos) {
    rows <- ((i-1)*n_alts + 1):(i*n_alts)
    probs <- exp(eta[rows] - max(eta[rows]))
    probs <- probs / sum(probs)
    khb_data$choice[rows[sample(n_alts, 1, prob = probs)]] <- 1L
  }

  result <- khb_decompose(
    data     = khb_data,
    key_vars = "x",
    z_vars   = "z",
    controls = "c1",
    strata   = "strata_id",
    choice   = "choice"
  )

  cat("\n   Decomposition table:\n")
  print(result$decomposition[, c("coefficient", "total_effect", "direct_effect",
                                  "indirect_effect", "conf_pct")])
  cat("\n   (True indirect effect ≈ 0.48 = 0.6 * 0.8)\n")
}

cat("\n=== Example complete ===\n")
