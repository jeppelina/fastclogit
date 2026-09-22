# test-reported-bugs.R: regressions for bugs reported from outside.
#
# Reported by Ben, 2026-09-22, against v0.3.0 from GitHub.

library(testthat)
library(fastclogit)

# A saturated assortative-mating design. The ego side (woman) is the stratum,
# so group_woman is constant WITHIN each stratum and is not identified in
# conditional logit. survival::clogit returns NA for those terms.
make_saturated_am <- function(n_sets = 200, n_alts = 20, A = 4, B = 4, seed = 1) {
  set.seed(seed)
  d <- data.frame(
    id_woman    = rep(seq_len(n_sets), each = n_alts),
    group_woman = factor(rep(sample(seq_len(A), n_sets, TRUE), each = n_alts)),
    group_man   = factor(sample(seq_len(B), n_sets * n_alts, TRUE))
  )
  eta <- 1.2 * (as.integer(d$group_woman) == as.integer(d$group_man))
  d$matched <- 0L
  for (s in split(seq_len(nrow(d)), d$id_woman)) {
    d$matched[s[sample.int(length(s), 1L, prob = exp(eta[s]))]] <- 1L
  }
  d
}

test_that("summary() prints when CI bounds fall outside [0, 1]", {
  # printCoefmat() takes the LAST column of the table as the p-value. The
  # confidence interval used to be appended after it, so symnum() was handed
  # the 97.5% bound and errored with "'x' must be between 0 and 1", or, when
  # the bounds happened to lie inside [0, 1], printed wrong stars in silence.
  d <- make_saturated_am()
  fit <- fclogit(matched ~ group_man + group_woman:group_man,
                 data = d, strata = "id_woman")

  s <- summary(fit)

  # The p-value must be last, and must be a probability.
  expect_equal(colnames(s$coef_table)[ncol(s$coef_table)], "Pr(>|z|)")
  pv <- s$coef_table[, ncol(s$coef_table)]
  expect_true(all(pv >= 0 & pv <= 1))

  # The interval is still available, just not in the printed table.
  expect_equal(ncol(s$conf_int), 2L)
  expect_true(any(s$conf_int > 1 | s$conf_int < 0),
              info = "test is only meaningful if some bound is outside [0,1]")

  expect_no_error(capture.output(print(s)))
})

test_that("a saturated a*b model with an ego-side factor does not error", {
  # The unidentified columns give an exactly singular Hessian. Before the ridge
  # fallback on the final vcov, inv_sympd() threw and the whole fit was lost
  # even though every identified coefficient was already correct.
  d <- make_saturated_am()
  fit <- NULL
  expect_warning(
    fit <- fclogit(matched ~ group_woman * group_man, data = d,
                   strata = "id_woman"),
    "not identified in conditional logit")
  expect_true(fit$converged)
})

test_that("identified coefficients match clogit exactly in the saturated model", {
  skip_if_not_installed("survival")
  d <- make_saturated_am()

  sv <- survival::clogit(matched ~ group_woman * group_man +
                           survival::strata(id_woman), data = d)
  fit <- suppressWarnings(
    fclogit(matched ~ group_woman * group_man, data = d, strata = "id_woman"))

  sv_coef <- coef(sv)
  identified <- names(sv_coef)[!is.na(sv_coef)]
  common <- intersect(identified, names(coef(fit)))
  expect_gt(length(common), 0L)

  expect_equal(unname(coef(fit)[common]), unname(sv_coef[common]),
               tolerance = 1e-6)

  # The terms clogit reports as NA must not be fitted here either. They are
  # dropped, named in the warning, and listed on the fit. Before v0.5.0 they
  # came back as ridge-determined values with meaningless standard errors and
  # no warning of any kind.
  aliased <- names(sv_coef)[is.na(sv_coef)]
  expect_gt(length(aliased), 0L)
  expect_false(any(aliased %in% names(coef(fit))))
  expect_setequal(fit$dropped_unidentified, aliased)
})
