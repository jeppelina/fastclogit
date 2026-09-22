# test-dense-sparse-consistency.R: Tier 2 of the test plan.
#
# Both paths computed FRESH in the same run and compared to each other. No
# cached reference, so nothing here can go stale.
#
# This is deliberately the live check for the open dense-path question:
# Paper 3 observed dense converging to a log-likelihood 271 units below sparse
# and survival at n=100 on real data, never reproduced locally
# (Paper 3/notes/technical/dense_kernel_local_simulation.md). If dense is
# genuinely broken, this is the test that should catch it, which is why the
# reference fits are generated from SPARSE and dense is checked against them
# live rather than being stored as the reference itself.

library(testthat)
library(fastclogit)

skip_if_no_matrix <- function() skip_if_not_installed("Matrix")

make_factor_heavy <- function(n_sets = 300, n_alts = 20, seed = 99) {
  set.seed(seed)
  n <- n_sets * n_alts
  d <- data.frame(
    strata = rep(seq_len(n_sets), each = n_alts),
    f1     = factor(sample.int(8L, n, TRUE)),
    f2     = factor(sample.int(6L, n, TRUE)),
    x      = rnorm(n)
  )
  # Drop the intercept column, NOT the reference levels. Keeping every level
  # of f1 would make its dummies sum to 1 in every row, which is a direction
  # constant within every stratum and therefore unidentified. The ridge then
  # makes vcov enormous along it and amplifies machine-epsilon differences
  # between the kernels into visible ones, so such a design cannot be used to
  # test cross-kernel agreement at all. See FASTCLOGIT_MERGE_MAP.md 5a.
  X <- model.matrix(~ f1 + f2 + x, data = d)[, -1, drop = FALSE]
  beta <- c(rnorm(ncol(X) - 1, 0, 0.4), -0.6)
  eta <- as.vector(X %*% beta)
  d$choice <- 0L
  for (s in split(seq_len(n), d$strata)) {
    d$choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  }
  list(X = X, choice = d$choice, strata = d$strata)
}

test_that("dense and sparse agree to machine epsilon", {
  skip_if_no_matrix()
  p <- make_factor_heavy()
  Xsp <- methods::as(p$X, "CsparseMatrix")

  fd <- fastclogit(p$X, p$choice, p$strata, max_iter = 100L)
  fs <- fastclogit(Xsp,  p$choice, p$strata, max_iter = 100L)

  expect_equal(unname(coef(fd)), unname(coef(fs)), tolerance = 1e-10)
  expect_equal(fd$loglik, fs$loglik, tolerance = 1e-10)
  expect_equal(unname(fd$se), unname(fs$se), tolerance = 1e-8)

  # The convergence ROUTE must agree too. Coefficients agreeing while the
  # routes differ means the two kernels are reaching the same place by
  # different ladders, which is the shape of the open dense/sparse question.
  expect_equal(fd$convergence_criterion, fs$convergence_criterion)
})

test_that("dense and sparse agree with a McFadden-Manski offset", {
  skip_if_no_matrix()
  p <- make_factor_heavy(seed = 123)
  n <- length(p$choice)
  offset <- c(2.9, 5.5, 8.1, 10.5)[sample.int(4L, n, TRUE)]

  fd <- fastclogit(p$X, p$choice, p$strata, offset = offset, max_iter = 100L)
  fs <- fastclogit(methods::as(p$X, "CsparseMatrix"), p$choice, p$strata,
                   offset = offset, max_iter = 100L)

  expect_equal(unname(coef(fd)), unname(coef(fs)), tolerance = 1e-9)
  expect_equal(fd$loglik, fs$loglik, tolerance = 1e-9)
})

test_that("dense and sparse agree on clustered sandwich SEs", {
  skip_if_no_matrix()
  p <- make_factor_heavy(seed = 321)
  cluster <- p$strata  # one cluster per choice set

  fd <- fastclogit(p$X, p$choice, p$strata, cluster = cluster, max_iter = 100L)
  fs <- fastclogit(methods::as(p$X, "CsparseMatrix"), p$choice, p$strata,
                   cluster = cluster, max_iter = 100L)

  expect_equal(unname(fd$se_robust), unname(fs$se_robust), tolerance = 1e-8)
})

test_that("the formula path and the matrix path give the same fit", {
  set.seed(8)
  sim <- simulate_clogit_data(n_egos = 400, n_alts = 15, seed = 8)
  fm <- fclogit(choice ~ x1 + x2 + x3,
                data = sim$data, strata = "strata_id")
  cols <- c("x1", "x2", "x3")
  fx <- fastclogit(as.matrix(sim$data[, cols]), sim$data$choice,
                   sim$data$strata_id)
  expect_equal(unname(coef(fm)[cols]), unname(coef(fx)[cols]), tolerance = 1e-10)
})
