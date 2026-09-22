# =============================================================================
# 01_dgp_equivalence.R
#
# The repo contains two DGPs: R/simulate_clogit.R draws the chosen alternative
# by SOFTMAX SAMPLING, tests/sparse_validation/gen_data.R adds Gumbel(0,1) to
# eta and takes the ARGMAX. These are the same distribution, which is why
# nobody has noticed the difference -- but it is an unasserted assumption
# linking two test harnesses, and if either is subtly wrong (a scale error on
# the Gumbel, a dropped offset in one path) a whole family of tests is
# validating against a mis-specified truth.
#
# Assert it once, including with an offset, since the offset is where a
# discrepancy would actually matter.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("DGP equivalence: softmax sampling vs Gumbel argmax")

rgumbel <- function(n) -log(-log(stats::runif(n)))

check_equivalence <- function(eta, R, label, seed = 1) {
  K <- length(eta)
  set.seed(seed)
  # Softmax sampling
  pr <- exp(eta - max(eta)); pr <- pr / sum(pr)
  a <- tabulate(sample.int(K, R, replace = TRUE, prob = pr), nbins = K)
  # Gumbel argmax
  b <- tabulate(
    max.col(matrix(rep(eta, each = R), nrow = R) + matrix(rgumbel(R * K), nrow = R),
            ties.method = "first"),
    nbins = K)

  # Two-sample chi-square on the realised counts.
  tab <- rbind(a, b)
  keep <- colSums(tab) >= 5           # chi-square validity
  ts <- suppressWarnings(stats::chisq.test(tab[, keep, drop = FALSE]))

  vc_note(sprintf("%-22s K=%d  R=%s  chi2=%.2f  df=%d  p=%.3f",
                  label, K, format(R, big.mark = ","),
                  ts$statistic, ts$parameter, ts$p.value))
  # Also compare to the ANALYTIC probabilities, which is the stronger check:
  # both empirical distributions must match the softmax they claim to encode.
  dev_a <- max(abs(a / R - pr)); dev_b <- max(abs(b / R - pr))
  vc_note(sprintf("%-22s max |empirical - analytic softmax|: sampling %.4f, gumbel %.4f",
                  "", dev_a, dev_b))
  list(p = ts$p.value, dev_a = dev_a, dev_b = dev_b, se = sqrt(0.25 / R))
}

R <- 200000L

# Decision rule, fixed in advance: the two mechanisms must be statistically
# indistinguishable (chi-square p > 0.001, a deliberately lenient threshold so
# that only a real discrepancy trips it), AND each must track the analytic
# softmax to within 4 binomial standard errors.
r1 <- check_equivalence(c(0, 0.5, -0.3, 1.2, -1.0), R, "no offset")
r2 <- check_equivalence(c(0, 0.5, -0.3, 1.2, -1.0) + c(2.9, 4.8, 6.6, 8.7, 10.5),
                        R, "with wide offset", seed = 2)

for (nm in c("no offset", "with wide offset")) {
  r <- if (nm == "no offset") r1 else r2
  vc_check("dgp-equivalence", paste0("mechanisms agree (", nm, ")"),
           r$p > 0.001, sprintf("chi-square p = %.3f", r$p))
  vc_check("dgp-equivalence", paste0("both match analytic softmax (", nm, ")"),
           max(r$dev_a, r$dev_b) < 4 * r$se,
           sprintf("max dev %.4f vs 4 SE = %.4f", max(r$dev_a, r$dev_b), 4 * r$se))
}
