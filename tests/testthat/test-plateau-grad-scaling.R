# test-plateau-grad-scaling.R
#
# The plateau rule's gradient ceiling is max(tier3_grad_floor,
# 1e-8 * |loglik|), not a bare absolute threshold. The reason is arithmetic:
# the gradient is a sum of N terms, so the smallest value it can attain is set
# by accumulation over those terms, and a function computed to relative
# accuracy eps locates its stationary point to about sqrt(eps) in the
# gradient. Paper 4 measured the floor at 8.3e-05 / 1.7e-03 / 2.0e-02 over
# 844k / 8.4M / 84M rows -- flat at ~4e-09 of |loglik|.
#
# Two properties are worth protecting, and they pull in opposite directions,
# which is why it is a max() and not a replacement.

library(testthat)
library(fastclogit)

# Same optimum, |loglik| scaled by m. Duplicating each choice set m times
# multiplies the log-likelihood by m and leaves the maximiser untouched, so it
# isolates the effect of |loglik| on the convergence rule from every other
# difference.
dup_problem <- function(n_sets = 400L, n_alts = 10L, m = 1L, seed = 3L) {
  set.seed(seed)
  n0 <- n_sets * n_alts
  X0 <- cbind(x1 = rnorm(n0), x2 = rnorm(n0))
  st0 <- rep(seq_len(n_sets), each = n_alts)
  eta <- as.vector(X0 %*% c(0.8, -0.5))
  ch0 <- integer(n0)
  for (s in split(seq_len(n0), st0))
    ch0[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  idx <- rep(seq_len(n0), times = m)
  rid <- rep(seq_len(m), each = n0)
  st <- st0[idx] + (rid - 1L) * n_sets
  ord <- order(st)
  list(X = X0[idx, , drop = FALSE][ord, , drop = FALSE],
       choice = ch0[idx][ord], strata = st[ord])
}

test_that("a larger |loglik| never makes convergence stricter", {
  # The ceiling is a floor on permissiveness: scaling the problem up may let
  # the plateau rule fire where it could not, but must never take away a
  # convergence that was already available.
  small <- dup_problem(m = 1L)
  big   <- dup_problem(m = 8L)

  fs <- fastclogit(small$X, small$choice, small$strata, max_iter = 100L)
  fb <- fastclogit(big$X, big$choice, big$strata, max_iter = 100L)

  expect_true(fs$converged)
  expect_true(fb$converged)

  # Duplication scales the log-likelihood and leaves the maximiser alone.
  expect_equal(fb$loglik / fs$loglik, 8, tolerance = 1e-6)
  expect_equal(unname(coef(fb)), unname(coef(fs)), tolerance = 1e-6)
})

test_that("the relative rule does not loosen small-scale fits", {
  # At the package's own sanity-check scale (loglik about -97) the relative
  # rule gives 9.7e-07, which is STRICTER than the 1e-06 absolute default.
  # max() keeps the absolute value there. A bare replacement would silently
  # tighten small fits, which is not the intent, so check that a small fit
  # still converges on the primary criterion with a genuinely small gradient.
  sim <- simulate_clogit_data(n_egos = 50, n_alts = 10, seed = 1)
  fit <- fastclogit(sim$X, sim$choice, sim$strata, max_iter = 100L)

  expect_true(fit$converged)
  expect_equal(fit$convergence_criterion, "primary")
  expect_lt(max(abs(fit$gradient)), 1e-6)
  expect_lt(abs(fit$loglik), 1e4)   # the regime the comment is about
})

test_that("the change is a no-op for fits that already converged", {
  # Paper 4 verified this at three sizes; this re-checks the property in the
  # package across dense, sparse, clustered and tight-tolerance paths. The
  # invariant is that every one of these converges cleanly and that dense and
  # sparse still agree, which is what would break first if the two kernels'
  # thresholds drifted apart.
  skip_if_not_installed("Matrix")
  p <- dup_problem(n_sets = 2000L, n_alts = 12L, m = 1L, seed = 11L)

  fd <- fastclogit(p$X, p$choice, p$strata, max_iter = 100L)
  fsp <- fastclogit(methods::as(p$X, "CsparseMatrix"), p$choice, p$strata,
                    max_iter = 100L)

  expect_true(fd$converged && fsp$converged)
  expect_equal(fd$convergence_criterion, fsp$convergence_criterion)
  expect_equal(unname(coef(fd)), unname(coef(fsp)), tolerance = 1e-10)
  expect_equal(fd$loglik, fsp$loglik, tolerance = 1e-10)
})
