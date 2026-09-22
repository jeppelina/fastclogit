# =============================================================================
# 05_cluster_robust.R: do the clustered sandwich SEs do their job?
#
# Nothing currently tests that the sandwich does anything useful. The kernels
# also apply a finite-sample correction C/(C-1) * (G-1)/G, described in the
# code as matching survival::coxph with cluster(); that claim is untested too.
#
# DESIGN NOTE. Two natural-looking designs induce NO dependence at all, and
# both fail silently -- see the long comment on vc_sim_clustered() in
# 00_helpers.R. The short version: a stratum-constant random effect cancels in
# the conditional likelihood, and correlating the covariates across an ego's
# choice sets does not correlate the scores, because the score has conditional
# mean zero given X whatever X does.
#
# This study therefore uses exact duplication, where the answer is known
# analytically. Each choice set is repeated m times and clustered on the
# original set. beta_hat is then identical to the single-copy estimate, so its
# true sampling SD is the single-copy SD; the model SE shrinks by sqrt(m);
# and the sandwich must undo exactly that factor.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("Cluster-robust standard errors")

R <- 400L
N_CL <- 400L; SETS <- 4L; N_ALTS <- 12L
BETA <- c(x1 = 0.8, x2 = -0.5)

vc_note(sprintf("R = %d, %d distinct choice sets each duplicated %d times, %d alts",
                R, N_CL, SETS, N_ALTS))
vc_note(sprintf("model SEs should be ~sqrt(%d) = %.2f times too small", SETS, sqrt(SETS)))

p <- length(BETA)
est <- matrix(NA_real_, R, p); se_m <- matrix(NA_real_, R, p)
se_r <- matrix(NA_real_, R, p)

t0 <- Sys.time()
for (r in seq_len(R)) {
  s <- vc_sim_clustered(n_clusters = N_CL, sets_per_cluster = SETS,
                        n_alts = N_ALTS, beta = BETA, seed = 30000L + r)
  f <- fastclogit(s$X, s$choice, s$strata, cluster = s$cluster, max_iter = 100L)
  est[r, ]  <- stats::coef(f)
  se_m[r, ] <- f$se
  se_r[r, ] <- f$se_robust
}
vc_note(sprintf("elapsed: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

z <- stats::qnorm(0.975)
truth <- matrix(BETA, R, p, byrow = TRUE)
emp_sd <- apply(est, 2, stats::sd)
cov_m <- colMeans(abs(est - truth) <= z * se_m)
cov_r <- colMeans(abs(est - truth) <= z * se_r)

tab <- data.frame(
  term = names(BETA), truth = as.numeric(BETA),
  bias = colMeans(est) - BETA,
  emp_sd = emp_sd,
  mean_se_model = colMeans(se_m), ratio_model = colMeans(se_m) / emp_sd,
  cov_model = cov_m,
  mean_se_robust = colMeans(se_r), ratio_robust = colMeans(se_r) / emp_sd,
  cov_robust = cov_r, row.names = NULL)
print(format(tab, digits = 4))

# --- Design check: the dependence must actually bite ------------------------
# If model SEs are already well calibrated there is no clustering problem to
# solve, and the robust arm passing would mean nothing.
vc_check("cluster-robust", "model SEs undercover [design check]",
         all(cov_m < 0.93),
         sprintf("model coverage %s (need < 0.93)",
                 paste(sprintf("%.3f", cov_m), collapse = "/")))

# --- The actual claim -------------------------------------------------------
# Band is wider than S1's: cluster-robust inference is asymptotic in the
# NUMBER OF CLUSTERS, not the number of observations.
vc_check("cluster-robust", "robust SEs restore nominal coverage",
         all(cov_r >= 0.93 & cov_r <= 0.97),
         sprintf("robust coverage %s", paste(sprintf("%.3f", cov_r), collapse = "/")))

vc_check("cluster-robust", "robust SEs match sampling variability",
         all(abs(colMeans(se_r) / emp_sd - 1) < 0.05),
         sprintf("ratio %s", paste(sprintf("%.3f", colMeans(se_r) / emp_sd), collapse = "/")))

# Analytic check: under exact m-fold duplication the model SE must be too
# small by exactly sqrt(m), and the sandwich must recover that factor.
infl <- colMeans(se_r) / colMeans(se_m)
vc_check("cluster-robust",
         sprintf("robust/model SE ratio equals sqrt(m) = %.3f", sqrt(SETS)),
         all(abs(infl / sqrt(SETS) - 1) < 0.05),
         sprintf("observed %s", paste(sprintf("%.3f", infl), collapse = "/")))

# --- The finite-sample correction claim -------------------------------------
# The kernels multiply by C/(C-1) * (G-1)/G and the comment says this matches
# survival::coxph with cluster(). Check it directly, on one dataset small
# enough for survival.
if (requireNamespace("survival", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survival))
  s <- vc_sim_clustered(n_clusters = 150L, sets_per_cluster = 3L, n_alts = 8L,
                        beta = BETA, seed = 77L)
  d <- data.frame(x1 = s$X[, 1], x2 = s$X[, 2], choice = s$choice,
                  sid = s$strata, cl = s$cluster)
  ff <- fastclogit(s$X, s$choice, s$strata, cluster = s$cluster, max_iter = 100L)
  sv <- clogit(choice ~ x1 + x2 + strata(sid), data = d, cluster = cl,
               method = "efron")
  se_sv <- sqrt(diag(stats::vcov(sv)))
  rel <- max(abs(ff$se_robust / se_sv - 1))
  vc_note(sprintf("fastclogit robust SE: %s", paste(sprintf("%.6f", ff$se_robust), collapse = " ")))
  vc_note(sprintf("survival  robust SE: %s", paste(sprintf("%.6f", se_sv), collapse = " ")))
  vc_check("cluster-robust", "robust SEs agree with survival::clogit",
           rel < 0.02, sprintf("max relative difference %.4f", rel))
}

saveRDS(tab, file.path(vc_dir, "results_05_cluster.rds"))
