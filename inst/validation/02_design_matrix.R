# =============================================================================
# 02_design_matrix.R — the interface studies.
#
# Two things nobody has looked at since they were written, and which produce a
# fit that converges beautifully to the WRONG answer when they are wrong.
#
#   A. fclogit builds its design matrix column-by-column by hand rather than
#      via model.matrix(), handling factors, reference levels and interactions
#      itself. Property-tested here over many random formulas by requiring the
#      formula path and the explicit-matrix path to produce the SAME FIT.
#      The log-likelihood is invariant to the basis, so it tests the span;
#      matching coefficients by name tests the parameterisation.
#
#   B. The collinearity dropper runs QR on a 50,000-row SUBSAMPLE and uses the
#      result to decide which columns to delete from a matrix that may have
#      tens of millions of rows. Rank on a subsample is not rank on the
#      population: a rare dummy absent from the subsample looks collinear.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

# -----------------------------------------------------------------------------
vc_header("A. fclogit formula path vs explicit model.matrix path")
# -----------------------------------------------------------------------------

make_random_case <- function(seed) {
  set.seed(seed)
  n_sets <- 150; n_alts <- 6; n <- n_sets * n_alts
  n_fac <- sample(1:3, 1); n_num <- sample(1:2, 1)
  d <- data.frame(sid = rep(seq_len(n_sets), each = n_alts))
  for (k in seq_len(n_fac))
    d[[paste0("f", k)]] <- factor(sample.int(sample(2:5, 1), n, TRUE))
  for (k in seq_len(n_num)) d[[paste0("v", k)]] <- round(rnorm(n), 3)

  tm <- c(paste0("f", seq_len(n_fac)), paste0("v", seq_len(n_num)))
  if (length(tm) >= 2 && runif(1) < 0.6)
    tm <- c(tm, paste(sample(tm, 2), collapse = ":"))
  fml <- stats::as.formula(paste("choice ~", paste(tm, collapse = " + ")))

  MM <- stats::model.matrix(stats::as.formula(paste("~", paste(tm, collapse = " + "))), data = d)
  MM <- MM[, colnames(MM) != "(Intercept)", drop = FALSE]
  keep <- apply(MM, 2, stats::var) > .Machine$double.eps
  MM <- MM[, keep, drop = FALSE]
  if (ncol(MM) < 1) return(NULL)

  b <- rnorm(ncol(MM), 0, 0.4)
  eta <- as.vector(MM %*% b)
  d$choice <- 0L
  for (s in split(seq_len(n), d$sid))
    d$choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  list(d = d, fml = fml, MM = MM)
}

n_trials <- 80L
bad <- character(0); checked <- 0L

for (i in seq_len(n_trials)) {
  cs <- make_random_case(i)
  if (is.null(cs)) next

  f_formula <- try(fclogit(cs$fml, data = cs$d, strata = "sid",
                           drop_collinear = FALSE, max_iter = 100L), silent = TRUE)
  f_matrix  <- try(fastclogit(cs$MM, cs$d$choice, cs$d$sid, max_iter = 100L),
                   silent = TRUE)
  if (inherits(f_formula, "try-error") || inherits(f_matrix, "try-error")) {
    bad <- c(bad, sprintf("%s -> error", deparse(cs$fml))); next
  }
  checked <- checked + 1L

  # Same span => same maximised log-likelihood, to machine precision.
  ll_ok <- abs(f_formula$loglik - f_matrix$loglik) <
    1e-8 * max(1, abs(f_matrix$loglik))

  # Same parameterisation => same coefficients under matching names.
  cf1 <- stats::coef(f_formula); cf2 <- stats::coef(f_matrix)
  common <- intersect(names(cf1), names(cf2))
  names_ok <- length(common) == length(cf2) && length(cf1) == length(cf2)
  coef_ok <- names_ok &&
    max(abs(cf1[common] - cf2[common])) < 1e-6

  if (!ll_ok || !coef_ok) {
    bad <- c(bad, sprintf("%s | dll=%.2e | ncol fc=%d mm=%d | matched=%d%s",
                          deparse(cs$fml),
                          abs(f_formula$loglik - f_matrix$loglik),
                          length(cf1), length(cf2), length(common),
                          if (names_ok) sprintf(" | dmax=%.2e",
                                                max(abs(cf1[common] - cf2[common]))) else ""))
  }
}

vc_note(sprintf("formulas exercised: %d of %d", checked, n_trials))
for (b in utils::head(bad, 12)) vc_note("MISMATCH: ", b)
if (length(bad) > 12) vc_note("... and ", length(bad) - 12, " more")

vc_check("design-matrix", "formula path == explicit matrix path",
         length(bad) == 0L,
         sprintf("%d mismatches over %d formulas", length(bad), checked))

# -----------------------------------------------------------------------------
vc_header("B. Subsample-QR collinearity detection on rare columns")
# -----------------------------------------------------------------------------
# fclogit: sample_n <- min(nrow(X), 50000L); qr(X[sample_idx, ]).
# With n > 50k the decision about which columns are collinear is made on a
# fraction of the data. A dummy whose few 1s all fall outside the subsample
# appears as a zero column, is called collinear, and is silently deleted --
# even though it is estimable in the full data.
#
# Pre-specified expectation: for a dummy with k ones in n rows and a subsample
# of m, the probability that every 1 is missed is about (1 - m/n)^k. At
# n = 200,000 and m = 50,000 that is 0.75^k: 32% at k = 4, 10% at k = 8.

n_sets <- 20000L; n_alts <- 10L; n <- n_sets * n_alts   # 200,000 rows
reps <- 25L
k_grid <- c(4L, 8L, 16L, 40L)

vc_note(sprintf("n = %s rows, QR subsample = 50,000 (%.0f%%), %d reps per cell",
                format(n, big.mark = ","), 100 * 50000 / n, reps))

drop_rate <- numeric(length(k_grid))
for (gi in seq_along(k_grid)) {
  k <- k_grid[gi]
  dropped <- 0L
  for (r in seq_len(reps)) {
    set.seed(1000L * gi + r)
    d <- data.frame(
      sid  = rep(seq_len(n_sets), each = n_alts),
      xnum = rnorm(n),
      rare = 0L
    )
    d$rare[sample.int(n, k)] <- 1L
    eta <- 0.7 * d$xnum + 2.0 * d$rare
    d$choice <- 0L
    for (s in split(seq_len(n), d$sid))
      d$choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L

    fit <- try(fclogit(choice ~ xnum + rare, data = d, strata = "sid",
                       drop_collinear = TRUE, max_iter = 100L), silent = TRUE)
    if (inherits(fit, "try-error")) { dropped <- dropped + 1L; next }
    if (!("rare" %in% names(stats::coef(fit)))) dropped <- dropped + 1L
  }
  drop_rate[gi] <- dropped / reps
  vc_note(sprintf("k = %2d ones : 'rare' wrongly dropped in %2d/%2d runs (%.0f%%)  [naive predicted %.0f%%]",
                  k, dropped, reps, 100 * drop_rate[gi], 100 * 0.75^k))
}

# Decision rule, fixed in advance: an estimable column must survive. Any
# non-zero drop rate for a column that IS identified in the full data is a
# defect, not a tuning matter.
vc_check("subsample-qr", "identified rare columns are never dropped",
         all(drop_rate == 0),
         sprintf("max drop rate %.0f%% (k = %d)",
                 100 * max(drop_rate), k_grid[which.max(drop_rate)]))

# Control: a genuinely collinear column MUST still be caught.
set.seed(99)
d <- data.frame(sid = rep(seq_len(2000L), each = 8L))
nn <- nrow(d)
d$a <- rnorm(nn); d$b <- rnorm(nn); d$dup <- d$a          # exact duplicate
eta <- 0.8 * d$a - 0.4 * d$b
d$choice <- 0L
for (s in split(seq_len(nn), d$sid))
  d$choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
fit <- fclogit(choice ~ a + b + dup, data = d, strata = "sid",
               drop_collinear = TRUE, max_iter = 100L)
vc_check("subsample-qr", "exactly collinear column is still detected",
         length(stats::coef(fit)) == 2L,
         sprintf("kept %d of 3 columns", length(stats::coef(fit))))
