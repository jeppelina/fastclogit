# =============================================================================
# 06_khb.R — does khb_decompose recover a known mediation structure?
#
# khb_decompose() is exported, is load-bearing for Papers 3 and 4, and is
# tested nowhere. Under a constructed DGP the decomposition has an analytic
# answer, so recovery is checkable.
#
# DGP:   z = a * x + e            (x mediated by z)
#        eta = b_x * x + b_z * z
#
# Substituting, eta = (b_x + b_z * a) * x + b_z * e, so on the scale of the
# FULL model:
#        direct   = b_x
#        indirect = b_z * a
#        total    = b_x + b_z * a
#
# The whole point of KHB is that these are all on the same scale, which naive
# cross-model comparison destroys through rescaling. If the package's
# implementation is right, the three quantities come back at their analytic
# values.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("KHB decomposition: recovery of a known mediation structure")

B_X <- 0.5; B_Z <- 0.8; A <- 0.6
DIRECT <- B_X; INDIRECT <- B_Z * A; TOTAL <- B_X + B_Z * A

R <- 200L
N_EGOS <- 800L; N_ALTS <- 15L

vc_note(sprintf("b_x = %.2f, b_z = %.2f, a = %.2f", B_X, B_Z, A))
vc_note(sprintf("analytic: direct = %.3f, indirect = %.3f, total = %.3f",
                DIRECT, INDIRECT, TOTAL))
vc_note(sprintf("R = %d, %d strata x %d alts", R, N_EGOS, N_ALTS))

make_med <- function(seed) {
  set.seed(seed)
  n <- N_EGOS * N_ALTS
  d <- data.frame(sid = rep(seq_len(N_EGOS), each = N_ALTS), x = rnorm(n))
  d$z <- A * d$x + rnorm(n)
  eta <- B_X * d$x + B_Z * d$z
  d$choice <- 0L
  for (s in split(seq_len(n), d$sid))
    d$choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  d
}

# One run first, to see what the function returns.
probe <- khb_decompose(make_med(1L), key_vars = "x", z_vars = "z",
                       strata = "sid", choice = "choice", verbose = FALSE)
cat("\n  khb_decompose() returns: ", paste(names(probe), collapse = ", "), "\n")
cat("  $decomposition columns:   ",
    paste(names(probe$decomposition), collapse = ", "), "\n\n")
print(probe$decomposition)

dec <- probe$decomposition
col_total    <- grep("^total$|total_", names(dec), value = TRUE)[1]
col_direct   <- grep("^direct$|direct_", names(dec), value = TRUE)[1]
col_indirect <- grep("^indirect$|indirect_", names(dec), value = TRUE)[1]

if (any(is.na(c(col_total, col_direct, col_indirect)))) {
  vc_check("khb", "decomposition columns found", FALSE,
           paste("got:", paste(names(dec), collapse = ", ")))
} else {
  res <- matrix(NA_real_, R, 3,
                dimnames = list(NULL, c("total", "direct", "indirect")))
  t0 <- Sys.time()
  for (r in seq_len(R)) {
    k <- try(khb_decompose(make_med(40000L + r), key_vars = "x", z_vars = "z",
                           strata = "sid", choice = "choice", verbose = FALSE),
             silent = TRUE)
    if (inherits(k, "try-error")) next
    row <- k$decomposition[1, ]
    res[r, ] <- c(row[[col_total]], row[[col_direct]], row[[col_indirect]])
  }
  vc_note(sprintf("elapsed: %.1f s (%d/%d runs succeeded)",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  sum(stats::complete.cases(res)), R))

  truth <- c(total = TOTAL, direct = DIRECT, indirect = INDIRECT)
  m  <- colMeans(res, na.rm = TRUE)
  sd_ <- apply(res, 2, stats::sd, na.rm = TRUE)
  mcse <- sd_ / sqrt(sum(stats::complete.cases(res)))

  tab <- data.frame(quantity = names(truth), analytic = as.numeric(truth),
                    mean_est = m, bias = m - truth, mc_se = mcse,
                    bias_in_mcse = (m - truth) / mcse, row.names = NULL)
  print(format(tab, digits = 4))

  # Pre-specified rule, same shape as the recovery study: bias within
  # max(3 MC SE, 0.1 sd). The 0.1-sd term keeps it valid as R grows.
  tol <- pmax(3 * mcse, 0.1 * sd_)
  vc_check("khb", "decomposition recovers the analytic values",
           all(abs(m - truth) <= tol),
           sprintf("max |bias| = %.4f vs tol %.4f", max(abs(m - truth)), max(tol)))

  # The identity must hold exactly in every run, not just on average.
  ident <- max(abs(res[, "total"] - res[, "direct"] - res[, "indirect"]),
               na.rm = TRUE)
  vc_check("khb", "total = direct + indirect holds exactly",
           ident < 1e-8, sprintf("max |total - direct - indirect| = %.2e", ident))

  saveRDS(tab, file.path(vc_dir, "results_06_khb.rds"))
}
