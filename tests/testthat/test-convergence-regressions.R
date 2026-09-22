# test-convergence-regressions.R: Tier 3 of the test plan.
#
# One constructed problem per defect that has actually bitten us. These assert
# on BEHAVIOUR (convergence route, iteration trace, absence of an error), not
# on stored numbers, so they cannot go stale the way a cached reference fit
# does. A test that only checked coefficients would have passed throughout the
# entire line-search-defect era; that is the mistake these are written against.
#
# Background: Research/KNOWN_ISSUES_fastclogit.md and
# Research/FASTCLOGIT_MERGE_MAP.md.

library(testthat)
library(fastclogit)

# ---------------------------------------------------------------------------
# A problem shaped like the one that broke Paper 4: a McFadden-Manski offset
# with wide within-stratum spread, plus rare dummies carrying large true
# coefficients. This is what makes the opening Newton step enormous.
# ---------------------------------------------------------------------------
# Tuned so that the first Newton step is genuinely too long and the line search
# has to halve it. A parameter sweep showed the true coefficient magnitude is
# what controls this, not the offset spread: at beta_rare ~3.2 no seed halved
# at iteration 1, at ~4.5 every seed did. Keep these values unless you re-run
# that sweep, or the defect-1 assertion below quietly stops testing anything.
make_offset_problem <- function(n_sets = 400, n_alts = 25, seed = 42) {
  set.seed(seed)
  n <- n_sets * n_alts
  strata <- rep(seq_len(n_sets), each = n_alts)

  # Five sampling strata with very different inclusion probabilities, so the
  # correction spans ~8 log points inside a single choice set.
  stratum  <- sample.int(5L, n, replace = TRUE)
  offset   <- c(2.9, 4.8, 6.6, 8.7, 10.5)[stratum]

  # Rare dummies (the ones the offset design exists to make estimable).
  uni  <- rbinom(n, 1, 0.012)
  work <- rbinom(n, 1, 0.0084)
  dist <- rnorm(n)

  X <- cbind(uni = uni, work = work, dist = dist)
  beta_true <- c(uni = 3.375, work = 4.5, dist = -0.5)

  eta <- as.vector(X %*% beta_true) + offset
  choice <- integer(n)
  for (s in split(seq_len(n), strata)) {
    choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  }
  list(X = X, choice = choice, strata = strata, offset = offset,
       beta_true = beta_true)
}

test_that("the line search runs on the FIRST iteration (defect 1)", {
  # THE DIRECT SIGNATURE. Under the old code the opening Newton step was taken
  # at full length unconditionally, so the halving count for the step INTO
  # iteration 2 was structurally zero. iter_log row 2 describes that step.
  #
  # Whether a given draw needs halving is luck (seeds 11 and 42 do not, most
  # others do), so run several and require every one of them to halve. That
  # keeps the assertion strict without resting on a single lucky seed.
  seeds <- c(1L, 2L, 3L)
  for (sd in seeds) {
    p <- make_offset_problem(seed = sd)
    fit <- fastclogit(p$X, p$choice, p$strata, offset = p$offset,
                      max_iter = 100L)

    expect_true(fit$converged, info = paste("seed", sd))
    expect_true(fit$convergence_criterion %in%
                  c("primary", "secondary", "plateau", "flat_optimum"),
                info = paste("seed", sd))

    expect_gt(nrow(fit$iter_log), 1L)
    expect_gt(fit$iter_log$halving_count[2], 0L)

    # And it must land on the truth rather than in the saturated region. The
    # dummies are rare by design, so judge against the fit's own SEs rather
    # than a fixed tolerance: the old kernel ended up hundreds of log-odds
    # away, so a 4-SE band separates the two outcomes with room to spare.
    for (nm in names(p$beta_true)) {
      expect_lt(abs(coef(fit)[[nm]] - p$beta_true[[nm]]), 4 * fit$se[[nm]])
    }
  }
})

test_that("the optimiser never walks downhill (defect 2)", {
  p <- make_offset_problem(seed = 7)
  fit <- fastclogit(p$X, p$choice, p$strata, offset = p$offset, max_iter = 100L)

  ll <- fit$iter_log$loglik
  # Every accepted step must improve, up to the acceptance tolerance. The old
  # code assigned beta_new even when 20 halvings found nothing, losing ~0.14
  # of log-likelihood per iteration for 90 iterations.
  expect_true(all(diff(ll) >= -1e-6),
              info = paste("loglik decreased at iter",
                           paste(which(diff(ll) < -1e-6) + 1L, collapse = ", ")))

  # $loglik is the best-loglik value, not wherever the loop happened to stop.
  expect_equal(fit$loglik, max(ll), tolerance = 1e-10)
  expect_equal(fit$iter_log$loglik[fit$best_loglik_iter], max(ll),
               tolerance = 1e-10)
})

test_that("a flat optimum converges as code 6, not as a line-search failure", {
  # THE REGRESSION THAT COST A KHB RUN. On 2026-09-12 the first version of the
  # no-improvement guard broke unconditionally, so a fit sitting AT its maximum
  # (log-likelihood flat, gradient tiny) was reported FAILED and the men's
  # decomposition was skipped for the sex. Splitting on tier3_grad_floor is
  # what tells "finished" apart from "the Newton direction died".
  #
  # Forced by starting from a converged fit and demanding a tolerance no step
  # can improve on.
  p <- make_offset_problem(seed = 11)
  fit <- fastclogit(p$X, p$choice, p$strata, offset = p$offset,
                    max_iter = 100L, tol = 1e-14)

  expect_true(fit$converged)
  expect_false(fit$convergence_criterion == "line_search")
  # Whichever route it takes, a converged fit must have a small gradient.
  expect_lt(max(abs(fit$gradient)), 1e-2)
})

test_that("convergence_code decodes to every documented criterion name", {
  # Guards against a kernel gaining a code the R side reports as 'unknown',
  # which is what the package did for codes 5 and 6 until this merge.
  sim <- simulate_clogit_data(n_egos = 300, n_alts = 15, seed = 3)
  fit <- fastclogit(sim$X, sim$choice, sim$strata)
  expect_true(fit$convergence_criterion %in%
                c("primary", "secondary", "plateau", "iter_max",
                  "line_search", "flat_optimum"))
  expect_false(fit$convergence_criterion == "unknown")
  expect_true(is.character(fit$convergence_message))
  expect_gt(nchar(fit$convergence_message), 0L)
})

test_that("iter_log is well formed", {
  sim <- simulate_clogit_data(n_egos = 300, n_alts = 15, seed = 5)
  fit <- fastclogit(sim$X, sim$choice, sim$strata)

  expect_s3_class(fit$iter_log, "data.frame")
  expect_equal(names(fit$iter_log),
               c("iter", "loglik", "grad_max", "rel_ll_change",
                 "step_size", "halving_count", "tier_fired"))
  expect_equal(nrow(fit$iter_log), fit$iterations)
  expect_equal(fit$iter_log$iter, seq_len(fit$iterations))
  # First row has no preceding step, so these are NA by construction.
  expect_true(is.na(fit$iter_log$rel_ll_change[1]))
  expect_true(is.na(fit$iter_log$step_size[1]))
  expect_true(all(fit$iter_log$tier_fired %in%
                    c("none", "primary", "secondary", "plateau")))
})
