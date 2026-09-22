# =============================================================================
# 04_offset.R — does the McFadden-Manski correction actually work?
#
# THE HIGHEST-CONSEQUENCE STUDY HERE. The correction is the load-bearing
# assumption in Papers 1, 3 and 4, and nothing anywhere tests that a fit on a
# stratified sample of alternatives recovers the parameters of the POPULATION
# choice model. The existing tests pass an offset through and check we agree
# with survival on the same offset, which tests plumbing, not correctness.
#
# Four arms. Three of them are supposed to FAIL: a validation that cannot fail
# is not a validation.
#
#   (a) correct offset                  -> must recover the population beta
#   (b) no offset                       -> must be biased
#   (c) permuted (wrong) offset         -> must be biased
#   (d) correct offset for sampled alts, but the CHOSEN alternative assigned
#       to the wrong stratum           -> must be biased
#
# Arm (d) is the realistic bug, not a strawman: both Paper 1 and Paper 4 carry
# a patch script (patch_actual_partner_correction.R, job_p4_patch_correction.R)
# whose entire job is assigning the actual chosen partner to the right
# sampling stratum. It is smaller than (b) and (c) and therefore nastier.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("McFadden-Manski offset correctness")

R <- 300L
N_EGOS <- 600L; N_POP <- 1500L; N_KEEP <- 25L
BETA <- c(x = 1.0, z = -0.8)

vc_note(sprintf("R = %d, %d egos, population choice set %d, %d alternatives kept",
                R, N_EGOS, N_POP, N_KEEP))
vc_note("sampling fractions across 5 strata: 0.60 / 0.20 / 0.08 / 0.03 / 0.01")
vc_note("covariates are CORRELATED with stratum, so omitting the offset must bite")

arms <- c("a_correct", "b_none", "c_permuted", "d_chosen_misassigned")
est <- array(NA_real_, c(R, length(BETA), length(arms)),
             dimnames = list(NULL, names(BETA), arms))
ses <- est

t0 <- Sys.time()
for (r in seq_len(R)) {
  sim <- vc_sim_sampled(n_egos = N_EGOS, n_pop = N_POP, n_keep = N_KEEP,
                        beta = BETA, seed = 20000L + r)
  d <- sim$data
  X <- as.matrix(d[, c("x", "z")])

  # (a) correct
  f <- fastclogit(X, d$choice, d$strata, offset = d$offset, max_iter = 100L)
  est[r, , "a_correct"] <- stats::coef(f); ses[r, , "a_correct"] <- f$se

  # (b) no offset at all
  f <- fastclogit(X, d$choice, d$strata, max_iter = 100L)
  est[r, , "b_none"] <- stats::coef(f); ses[r, , "b_none"] <- f$se

  # (c) offsets permuted across strata types within each choice set
  off_c <- d$offset
  for (s in split(seq_len(nrow(d)), d$strata)) {
    u <- unique(d$stratum[s])
    if (length(u) > 1) {
      map <- stats::setNames(sample(u), u)
      off_c[s] <- d$offset[s][match(map[as.character(d$stratum[s])], d$stratum[s])]
      off_c[s][is.na(off_c[s])] <- d$offset[s][is.na(off_c[s])]
    }
  }
  f <- fastclogit(X, d$choice, d$strata, offset = off_c, max_iter = 100L)
  est[r, , "c_permuted"] <- stats::coef(f); ses[r, , "c_permuted"] <- f$se

  # (d) chosen alternative given the offset of the most heavily sampled
  #     stratum instead of its own -- the realistic misassignment
  off_d <- d$offset
  ch <- which(d$choice == 1L)
  for (s in split(seq_len(nrow(d)), d$strata)) {
    i <- intersect(s, ch)
    if (length(i)) {
      cand <- s[d$stratum[s] == 1L]
      if (length(cand)) off_d[i] <- d$offset[cand[1]]
    }
  }
  f <- fastclogit(X, d$choice, d$strata, offset = off_d, max_iter = 100L)
  est[r, , "d_chosen_misassigned"] <- stats::coef(f)
  ses[r, , "d_chosen_misassigned"] <- f$se
}
vc_note(sprintf("elapsed: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

z <- stats::qnorm(0.975)
summ <- do.call(rbind, lapply(arms, function(a) {
  e <- est[, , a, drop = TRUE]; s <- ses[, , a, drop = TRUE]
  data.frame(
    arm = a, term = names(BETA), truth = as.numeric(BETA),
    mean_est = colMeans(e),
    bias = colMeans(e) - BETA,
    mc_se = apply(e, 2, stats::sd) / sqrt(R),
    bias_in_mcse = (colMeans(e) - BETA) / (apply(e, 2, stats::sd) / sqrt(R)),
    coverage = colMeans(abs(e - matrix(BETA, R, length(BETA), byrow = TRUE)) <= z * s),
    row.names = NULL)
}))
print(format(summ, digits = 4))

# --- Decision rules, fixed in advance ---------------------------------------
ga <- summ[summ$arm == "a_correct", ]
bias_tol <- pmax(3 * ga$mc_se, 0.1 * apply(est[, , "a_correct"], 2, stats::sd))
vc_check("offset", "(a) correct offset recovers the population beta",
         all(abs(ga$bias) <= bias_tol),
         sprintf("max |bias| = %.4f vs tol %.4f", max(abs(ga$bias)), max(bias_tol)))

cb <- vc_cov_band(R)
jc <- vc_joint_rule(ga$coverage, cb[1], cb[2])
vc_check("offset", "(a) intervals cover at the nominal rate", jc$pass,
         sprintf("coverage %s, band [%.3f, %.3f]",
                 paste(sprintf("%.3f", ga$coverage), collapse = "/"), cb[1], cb[2]))

# Falsification arms. Each MUST be detectably biased; if one is not, the
# design has failed to create the problem the correction exists to solve and
# arm (a) passing means nothing.
for (a in c("b_none", "c_permuted", "d_chosen_misassigned")) {
  g <- summ[summ$arm == a, ]
  worst <- max(abs(g$bias_in_mcse))
  vc_check("offset", sprintf("(%s) is detectably biased [design check]",
                             substr(a, 1, 1)),
           worst > 5,
           sprintf("max |bias| = %.1f MC SE (need > 5)", worst))
}

saveRDS(summ, file.path(vc_dir, "results_04_offset.rds"))
