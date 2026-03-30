# test-basic.R — Validate fastclogit against survival::clogit
#
# These tests use small simulated data where clogit() is tractable,
# and verify that fastclogit produces matching results.

library(testthat)

test_that("fastclogit matches clogit on simple data (no offset, no cluster)", {
  skip_if_not_installed("survival")
  library(survival)

  sim <- simulate_clogit_data(
    n_egos = 500, n_alts = 20,
    use_offset = FALSE, cluster_ratio = 1.0, seed = 123
  )

  # fastclogit
  fit_fast <- fastclogit(sim$X, sim$choice, sim$strata,
                          max_iter = 50, tol = 1e-10)

  # survival::clogit
  # Build formula from column names
  pred_names <- colnames(sim$X)
  fml <- as.formula(paste("choice ~",
                          paste(pred_names, collapse = " + "),
                          "+ strata(strata_id)"))
  fit_clogit <- clogit(fml, data = sim$data, method = "efron")

  # Compare coefficients
  beta_fast   <- fit_fast$coefficients
  beta_clogit <- coef(fit_clogit)

  # Match names (clogit may reorder)
  common_names <- intersect(names(beta_fast), names(beta_clogit))
  expect_true(length(common_names) > 0)

  for (nm in common_names) {
    expect_equal(beta_fast[nm], beta_clogit[nm], tolerance = 1e-5,
                 label = paste("coefficient:", nm))
  }

  # Compare log-likelihood
  expect_equal(fit_fast$loglik, as.numeric(fit_clogit$loglik[2]),
               tolerance = 1e-4, label = "log-likelihood")

  # Compare model-based SEs
  se_fast   <- fit_fast$se[common_names]
  se_clogit <- sqrt(diag(vcov(fit_clogit)))[common_names]
  for (nm in common_names) {
    expect_equal(se_fast[nm], se_clogit[nm], tolerance = 1e-4,
                 label = paste("SE:", nm))
  }

  expect_true(fit_fast$converged)
})


test_that("fastclogit matches clogit WITH offset", {
  skip_if_not_installed("survival")
  library(survival)

  sim <- simulate_clogit_data(
    n_egos = 500, n_alts = 20,
    use_offset = TRUE, cluster_ratio = 1.0, seed = 456
  )

  # fastclogit
  fit_fast <- fastclogit(sim$X, sim$choice, sim$strata,
                          offset = sim$offset,
                          max_iter = 50, tol = 1e-10)

  # survival::clogit with offset
  pred_names <- colnames(sim$X)
  fml <- as.formula(paste("choice ~",
                          paste(pred_names, collapse = " + "),
                          "+ offset(correction) + strata(strata_id)"))
  fit_clogit <- clogit(fml, data = sim$data, method = "efron")

  # Compare coefficients
  beta_fast   <- fit_fast$coefficients
  beta_clogit <- coef(fit_clogit)
  common_names <- intersect(names(beta_fast), names(beta_clogit))

  for (nm in common_names) {
    expect_equal(beta_fast[nm], beta_clogit[nm], tolerance = 1e-5,
                 label = paste("offset coefficient:", nm))
  }

  # Compare log-likelihood
  expect_equal(fit_fast$loglik, as.numeric(fit_clogit$loglik[2]),
               tolerance = 1e-4, label = "offset log-likelihood")
})


test_that("fastclogit matches clogit with clustered SEs", {
  skip_if_not_installed("survival")
  library(survival)

  sim <- simulate_clogit_data(
    n_egos = 500, n_alts = 20,
    use_offset = TRUE, cluster_ratio = 1.5, seed = 789
  )

  # fastclogit with clustering
  fit_fast <- fastclogit(sim$X, sim$choice, sim$strata,
                          offset = sim$offset, cluster = sim$cluster,
                          max_iter = 50, tol = 1e-10)

  # survival::clogit with cluster()
  pred_names <- colnames(sim$X)
  fml <- as.formula(paste("choice ~",
                          paste(pred_names, collapse = " + "),
                          "+ offset(correction) + strata(strata_id) + cluster(cluster_id)"))
  fit_clogit <- clogit(fml, data = sim$data, method = "efron")

  # Compare coefficients (should be identical — clustering only affects SEs)
  beta_fast   <- fit_fast$coefficients
  beta_clogit <- coef(fit_clogit)
  common_names <- intersect(names(beta_fast), names(beta_clogit))

  for (nm in common_names) {
    expect_equal(beta_fast[nm], beta_clogit[nm], tolerance = 1e-5,
                 label = paste("cluster coefficient:", nm))
  }

  # Compare robust SEs
  # clogit stores robust variance when cluster() is used
  # Use names from the clogit variance matrix directly (may differ from fastclogit names)
  se_fast_robust   <- fit_fast$se_robust[common_names]
  clogit_var_names <- names(coef(fit_clogit))
  se_clogit_robust <- setNames(sqrt(diag(fit_clogit$var)), clogit_var_names)

  for (nm in common_names) {
    expect_equal(se_fast_robust[nm], se_clogit_robust[nm], tolerance = 1e-3,
                 label = paste("robust SE:", nm))
  }
})


test_that("coefficient recovery from known DGP", {
  sim <- simulate_clogit_data(
    n_egos = 10000, n_alts = 100,
    use_offset = FALSE, cluster_ratio = 1.0, seed = 42
  )

  fit <- fastclogit(sim$X, sim$choice, sim$strata,
                     max_iter = 50, tol = 1e-10)

  # With 10K egos and 100 alts, estimates should be close to truth
  beta_true <- sim$beta_true
  beta_hat  <- fit$coefficients

  for (nm in names(beta_true)) {
    se_nm <- fit$se[nm]
    # Should be within ~3 SEs of truth
    expect_true(abs(beta_hat[nm] - beta_true[nm]) < 3 * se_nm,
                label = paste("recovery:", nm,
                              "| true:", round(beta_true[nm], 3),
                              "| est:", round(beta_hat[nm], 3),
                              "| SE:", round(se_nm, 3)))
  }
})


test_that("edge case: groups of size 2", {
  sim <- simulate_clogit_data(
    n_egos = 1000, n_alts = 2,
    use_offset = FALSE, cluster_ratio = 1.0, seed = 101
  )

  fit <- fastclogit(sim$X, sim$choice, sim$strata, max_iter = 50)
  expect_true(fit$converged)
  expect_equal(fit$n_groups, 1000)
})


test_that("S3 methods work", {
  sim <- simulate_clogit_data(n_egos = 200, n_alts = 10, seed = 55)
  fit <- fastclogit(sim$X, sim$choice, sim$strata,
                     offset = sim$offset, cluster = sim$cluster)

  # coef
  expect_equal(length(coef(fit)), ncol(sim$X))

  # vcov
  v <- vcov(fit)
  expect_equal(nrow(v), ncol(sim$X))
  expect_equal(ncol(v), ncol(sim$X))

  # confint
  ci <- confint(fit)
  expect_equal(nrow(ci), ncol(sim$X))
  expect_equal(ncol(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))

  # logLik
  ll <- logLik(fit)
  expect_true(is.numeric(ll))
  expect_true(ll < 0)

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fastclogit")
  expect_equal(nrow(s$coef_table), ncol(sim$X))

  # tidy
  td <- tidy_fastclogit(fit)
  expect_true(is.data.frame(td))
  expect_equal(nrow(td), ncol(sim$X))
  expect_true(all(c("term", "estimate", "std.error", "p.value") %in% names(td)))

  # tidy with exponentiation (odds ratios)
  td_or <- tidy_fastclogit(fit, exponentiate = TRUE)
  expect_true(all(td_or$estimate > 0))
})
