# =============================================================================
# 03_recovery_coverage.R: the statistical backbone.
#
# Nothing in the package's test suite checks that the estimator is unbiased,
# that the reported SEs match the sampling variability of the estimator, or
# that a nominal 95% interval covers 95% of the time. A systematically wrong
# variance formula passes every existing test, because survival::clogit would
# have to be wrong in the same way for the comparison to fail.
#
# Decision rules are calibrated to their own Monte Carlo error (reviewer C4)
# and applied jointly across coefficients (reviewer C5), because per-coefficient
# 95% bands on 11 coefficients produce a false-alarm rate near 43% on correct
# code.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("Recovery, SE calibration and coverage")

R <- 2000L
N_EGOS <- 2000L; N_ALTS <- 20L

sim0 <- simulate_clogit_data(n_egos = N_EGOS, n_alts = N_ALTS,
                             use_offset = FALSE, cluster_ratio = 1.0, seed = 1)
beta_true <- sim0$beta_true
p <- length(beta_true)

vc_note(sprintf("R = %d replications, %d strata x %d alts, p = %d",
                R, N_EGOS, N_ALTS, p))
vc_note("model-based SEs (no clustering in this DGP)")

est <- matrix(NA_real_, R, p, dimnames = list(NULL, names(beta_true)))
ses <- matrix(NA_real_, R, p, dimnames = list(NULL, names(beta_true)))
routes <- character(R)

t0 <- Sys.time()
for (r in seq_len(R)) {
  s <- simulate_clogit_data(n_egos = N_EGOS, n_alts = N_ALTS,
                            use_offset = FALSE, cluster_ratio = 1.0,
                            seed = 10000L + r)
  f <- fastclogit(s$X, s$choice, s$strata, max_iter = 100L)
  est[r, ] <- stats::coef(f)[names(beta_true)]
  ses[r, ] <- f$se[names(beta_true)]
  routes[r] <- f$convergence_criterion
}
vc_note(sprintf("elapsed: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
vc_note("convergence routes: ",
        paste(sprintf("%s=%d", names(table(routes)), table(routes)), collapse = "  "))

# --- Estimands ---------------------------------------------------------------
bias      <- colMeans(est) - beta_true
mc_se     <- apply(est, 2, stats::sd) / sqrt(R)
emp_sd    <- apply(est, 2, stats::sd)
se_ratio  <- colMeans(ses) / emp_sd
z         <- stats::qnorm(0.975)
covered   <- abs(est - matrix(beta_true, R, p, byrow = TRUE)) <= z * ses
coverage  <- colMeans(covered)

tab <- data.frame(
  term = names(beta_true), truth = as.numeric(beta_true),
  mean_est = colMeans(est), bias = bias, mc_se = mc_se,
  emp_sd = emp_sd, mean_se = colMeans(ses), se_ratio = se_ratio,
  coverage = coverage, row.names = NULL)
print(format(tab, digits = 4))

# --- Pre-specified decision rules -------------------------------------------
# Bias: |bias| < max(3 * MC SE, 0.1 * empirical SD). The second term keeps the
# rule valid as R grows: the conditional MLE carries finite-sample bias of
# order 1/n_strata, which becomes detectable -- and would reject a CORRECT
# estimator -- once the MC SE shrinks below it (reviewer C6).
bias_tol <- pmax(3 * mc_se, 0.1 * emp_sd)
bias_breach <- sum(abs(bias) > bias_tol)
vc_check("recovery", "coefficients are unbiased",
         bias_breach == 0L,
         sprintf("%d/%d exceed max(3 MC SE, 0.1 sd)", bias_breach, p))

# SE calibration: band set from the MC error of an SD estimated on R draws,
# which is ~1/sqrt(2(R-1)) in relative terms. A fixed +/-3% band would be
# NARROWER than its own noise at R = 500.
band <- vc_se_ratio_band(R, k = 3)
j <- vc_joint_rule(se_ratio, band[1], band[2])
vc_note(sprintf("SE-ratio band at R=%d: [%.4f, %.4f]", R, band[1], band[2]))
vc_check("recovery", "reported SEs match sampling variability",
         j$pass, j$detail)

# Coverage: joint rule over the same band logic.
cb <- vc_cov_band(R)
jc <- vc_joint_rule(coverage, cb[1], cb[2])
vc_note(sprintf("coverage band at R=%d: [%.4f, %.4f]", R, cb[1], cb[2]))
vc_check("recovery", "95% intervals cover at the nominal rate",
         jc$pass, jc$detail)

saveRDS(list(tab = tab, routes = routes, R = R),
        file.path(vc_dir, "results_03_recovery.rds"))
