# =============================================================================
# 00_helpers.R: shared machinery for the validation studies.
#
# Every study states its decision rule BEFORE it runs, as a call to
# vc_check(). The rules live in the scripts, in git, so the plan is its own
# pre-registration: changing a rule after seeing a result shows up as a diff.
#
# See notes/VALIDATION_PLAN.md for the reasoning behind each rule.
# =============================================================================

suppressPackageStartupMessages({
  library(fastclogit)
})

# One shared results table across every study, even when run_all.R sources each
# study into its own environment.
if (!exists("VC_RESULTS", envir = globalenv(), inherits = FALSE)) {
  .vc_env <- new.env(parent = emptyenv())
  .vc_env$rows <- list()
  assign("VC_RESULTS", .vc_env, envir = globalenv())
  rm(.vc_env)
}
VC_RESULTS <- get("VC_RESULTS", envir = globalenv())

vc_header <- function(txt) {
  cat("\n", strrep("=", 76), "\n", txt, "\n", strrep("=", 76), "\n", sep = "")
}

vc_note <- function(...) cat("  ", ..., "\n", sep = "")

# Record and print one pre-specified check.
vc_check <- function(study, what, pass, detail = "") {
  VC_RESULTS$rows[[length(VC_RESULTS$rows) + 1L]] <-
    data.frame(study = study, check = what, pass = isTRUE(pass),
               detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %-52s %s\n",
              if (isTRUE(pass)) "PASS" else "FAIL", what, detail))
  invisible(pass)
}

vc_summary <- function(path = NULL) {
  if (!length(VC_RESULTS$rows)) return(invisible(NULL))
  d <- do.call(rbind, VC_RESULTS$rows)
  vc_header("VALIDATION SUMMARY")
  for (s in unique(d$study)) {
    ds <- d[d$study == s, ]
    cat(sprintf("  %-34s %2d/%2d passed%s\n", s, sum(ds$pass), nrow(ds),
                if (all(ds$pass)) "" else "   <-- SEE ABOVE"))
  }
  cat(sprintf("\n  TOTAL: %d/%d checks passed\n", sum(d$pass), nrow(d)))
  if (!is.null(path)) {
    write.csv(d, path, row.names = FALSE)
    cat("  results ->", path, "\n")
  }
  invisible(d)
}

# --- Monte Carlo decision-rule helpers ---------------------------------------
# Reviewer C4: the SE-calibration band must be wider than its own MC error.
# sd() estimated from R draws has relative MC error ~ 1/sqrt(2(R-1)).
vc_se_ratio_band <- function(R, k = 3) 1 + c(-1, 1) * k / sqrt(2 * (R - 1))

# Reviewer C5: with p coefficients, per-coefficient 95% bands give a high
# family-wise false-alarm rate. Pre-specified joint rule: a study fails if TWO
# OR MORE coefficients breach, or if ANY breaches by more than 2x the band
# half-width. One marginal excursion is expected and is not a failure.
vc_joint_rule <- function(values, lo, hi) {
  half <- (hi - lo) / 2
  ctr  <- (hi + lo) / 2
  dev  <- abs(values - ctr)
  n_breach <- sum(dev > half)
  worst    <- max(dev / half)
  list(pass = (n_breach < 2L) && (worst <= 2),
       n_breach = n_breach, worst = worst,
       detail = sprintf("%d/%d outside band, worst = %.2fx half-width",
                        n_breach, length(values), worst))
}

# Coverage band for R replications at nominal level.
vc_cov_band <- function(R, level = 0.95, k = 2)
  level + c(-1, 1) * k * sqrt(level * (1 - level) / R)

# --- DGP: genuine within-cluster SCORE dependence, by duplication ------------
# Getting this design right took three attempts; the two that failed are worth
# recording because both fail SILENTLY, looking like a well-calibrated model.
#
#   (1) A random intercept shared within a choice set. Constant within the
#       stratum, so it cancels exactly in the conditional likelihood --
#       conditional logit conditions on the stratum. Induces NOTHING.
#
#   (2) Covariates correlated across an ego's several choice sets. Also
#       induces nothing, and this one is subtler: the per-stratum score is
#       x_chosen - E[x], which has conditional mean zero given X at the true
#       beta. So Cov(s_g1, s_g2) = E[Cov(s1,s2|X)] + Cov(E[s1|X], E[s2|X])
#       = 0 + Cov(0, 0) = 0. Correlating the DESIGN does not correlate the
#       scores, because the scores are a martingale difference sequence
#       whatever X does. Measured model-SE coverage under this design: 0.965
#       and 0.950, i.e. no clustering problem at all.
#
# What does work is exact duplication, which is also the cleanest possible
# test of a cluster-robust estimator. Each of `n_clusters` distinct choice
# sets is repeated `sets_per_cluster` times and the cluster is the original
# set. The log-likelihood is then m copies of the single-copy likelihood, so:
#   - beta_hat is IDENTICAL to the single-copy estimate, hence its true
#     sampling SD is the single-copy SD;
#   - the model-based SE shrinks by sqrt(m), because the fit believes it has
#     m times the information;
#   - a sandwich clustered on the original set must undo exactly that.
# The truth is known analytically, so the test is unambiguous.
vc_sim_clustered <- function(n_clusters = 800, sets_per_cluster = 4,
                             n_alts = 12, beta = c(x1 = 0.8, x2 = -0.5),
                             seed = 1) {
  set.seed(seed)
  p <- length(beta)
  n0 <- n_clusters * n_alts

  X0 <- matrix(rnorm(n0 * p), n0, p, dimnames = list(NULL, names(beta)))
  strata0 <- rep(seq_len(n_clusters), each = n_alts)
  eta <- as.vector(X0 %*% beta)
  choice0 <- integer(n0)
  for (s in split(seq_len(n0), strata0)) {
    choice0[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  }

  # Replicate each choice set m times; cluster = the original set.
  m <- sets_per_cluster
  idx     <- rep(seq_len(n0), times = m)
  rep_id  <- rep(seq_len(m), each = n0)
  strata  <- strata0[idx] + (rep_id - 1L) * n_clusters
  cluster <- strata0[idx]
  ord <- order(strata)

  list(X = X0[idx, , drop = FALSE][ord, , drop = FALSE],
       choice  = choice0[idx][ord],
       strata  = strata[ord],
       cluster = cluster[ord],
       beta_true = beta,
       m = m,
       single = list(X = X0, choice = choice0, strata = strata0))
}

# --- DGP: population choice set + stratified alternative sampling ------------
# Used by the offset study. Builds a POPULATION of alternatives per ego,
# draws the chosen one from the full population softmax, then samples
# alternatives at deliberately unequal rates across strata.
#
# The stratification MUST be correlated with the model covariates, or omitting
# the offset induces no bias and the study's falsification arm passes
# vacuously (reviewer C8).
vc_sim_sampled <- function(n_egos = 1000, n_pop = 2000, n_keep = 25,
                           beta = c(x = 1.0, z = -0.8),
                           keep_frac = c(0.60, 0.20, 0.08, 0.03, 0.01),
                           seed = 1) {
  set.seed(seed)
  S <- length(keep_frac)
  rows <- vector("list", n_egos)

  for (i in seq_len(n_egos)) {
    # Stratum membership, and covariates CORRELATED with stratum so that
    # dropping the offset must bite.
    stratum <- sample.int(S, n_pop, replace = TRUE)
    x <- rnorm(n_pop, mean = (stratum - mean(seq_len(S))) * 0.45)
    z <- rnorm(n_pop, mean = (stratum - mean(seq_len(S))) * -0.30)
    eta <- beta[["x"]] * x + beta[["z"]] * z
    chosen_idx <- sample.int(n_pop, 1L, prob = exp(eta - max(eta)))

    # Stratified sample of alternatives, chosen one always retained.
    # SAMPLING PROTOCOL. This matters more than it looks: the correction
    # formula is protocol-specific, and the standard formula is only valid for
    # one particular protocol.
    #
    # Here: draw exactly n_s alternatives from stratum s, and when the chosen
    # alternative lies in stratum s it OCCUPIES ONE OF THOSE n_s SLOTS (so we
    # draw n_s - 1 more from that stratum's pool). Every sampled choice set
    # then has the same stratum composition regardless of which alternative
    # was chosen.
    #
    # Under this protocol, with D the sampled set:
    #   q(D|i) = 1 / [ C(N_s - 1, n_s - 1) * prod_{t != s} C(N_t, n_t) ]
    # and since C(N-1, n-1) = C(N, n) * n/N,
    #   ln q(D|i) = const - ln(n_s / N_s)
    # so the correction is c_s = -log(n_s / N_s), the familiar formula.
    #
    # The first version of this helper instead forced the chosen in ON TOP of
    # n_s draws per stratum. That is a different protocol, for which the
    # familiar formula is simply wrong, and it biased the estimator by +0.10
    # on a true coefficient of 1.0 with 57% coverage. The bias came from the
    # simulation, not from fastclogit -- but it is a fair illustration of why
    # Paper 1 and Paper 4 both carry a patch script for exactly this
    # bookkeeping.
    N_s <- tabulate(stratum, nbins = S)
    n_s <- pmax(1L, round(n_keep * keep_frac / sum(keep_frac)))
    n_s <- pmin(n_s, N_s)
    s_chosen <- stratum[chosen_idx]

    keep <- chosen_idx
    for (s in seq_len(S)) {
      pool <- setdiff(which(stratum == s), chosen_idx)
      take <- if (s == s_chosen) n_s[s] - 1L else n_s[s]
      take <- min(take, length(pool))
      if (take > 0) keep <- c(keep, sample(pool, take))
    }

    corr_s <- -log(n_s / N_s)

    rows[[i]] <- data.frame(
      strata  = i,
      x       = x[keep],
      z       = z[keep],
      stratum = stratum[keep],
      choice  = as.integer(keep == chosen_idx),
      offset  = corr_s[stratum[keep]]
    )
  }
  d <- do.call(rbind, rows)
  list(data = d, beta_true = beta, n_strata_types = S)
}
