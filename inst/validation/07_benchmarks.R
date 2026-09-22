# =============================================================================
# 07_benchmarks.R — reproducible timing and memory, plus the tol-scaling
# observation.
#
# Every configuration runs in its OWN PROCESS under /usr/bin/time -l, so the
# reported peak RSS includes Armadillo's allocations. See 07_bench_worker.R
# for why gc() is the wrong instrument here.
#
# The README has long quoted "225 GB to 30 GB" and "95 minutes to 80 seconds"
# from Paper 3 production runs. Those are real, but a reader cannot reproduce
# them. This produces numbers a reader can reproduce; the production figures
# are cited separately as what they are.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("Benchmarks: wall time and peak RSS, measured per process")

worker <- file.path(vc_dir, "07_bench_worker.R")
have_time <- file.exists("/usr/bin/time")
if (!have_time) vc_note("/usr/bin/time not found; RSS will be reported as NA")

run_one <- function(engine, n_sets, n_alts, p_fac) {
  a <- c("--vanilla", worker, engine, n_sets, n_alts, p_fac)
  if (have_time) {
    out <- suppressWarnings(system2("/usr/bin/time",
                                    c("-l", "Rscript", a),
                                    stdout = TRUE, stderr = TRUE))
  } else {
    out <- suppressWarnings(system2("Rscript", a, stdout = TRUE, stderr = TRUE))
  }
  res <- grep("^RESULT\t", out, value = TRUE)
  if (!length(res)) return(NULL)
  f <- strsplit(res[1], "\t")[[1]]
  rss_line <- grep("maximum resident set size", out, value = TRUE)
  rss_gb <- if (length(rss_line))
    as.numeric(sub("^\\s*(\\d+).*", "\\1", rss_line[1])) / 2^30 else NA_real_
  data.frame(engine = f[2], n_rows = as.integer(f[3]), p = as.integer(f[4]),
             n_sets = as.integer(f[5]), n_alts = as.integer(f[6]),
             seconds = as.numeric(f[7]), loglik = as.numeric(f[8]),
             iters = f[9], peak_rss_gb = rss_gb, stringsAsFactors = FALSE)
}

grid <- list(
  list(n_sets =  10000L, n_alts = 10L, p_fac = 20L, survival = TRUE),
  list(n_sets =  50000L, n_alts = 10L, p_fac = 50L, survival = FALSE),
  list(n_sets =  50000L, n_alts = 20L, p_fac = 90L, survival = FALSE)
)

rows <- list()
for (g in grid) {
  for (eng in c("dense", "sparse", if (g$survival) "survival")) {
    vc_note(sprintf("running %-9s n_sets=%s n_alts=%d p_fac=%d ...",
                    eng, format(g$n_sets, big.mark = ","), g$n_alts, g$p_fac))
    r <- run_one(eng, g$n_sets, g$n_alts, g$p_fac)
    if (!is.null(r)) rows[[length(rows) + 1L]] <- r
  }
}
bench <- do.call(rbind, rows)
cat("\n")
print(format(bench, digits = 4))

# Cross-engine agreement is a correctness check riding along for free: the
# three engines must reach the same log-likelihood on the same data.
for (nr in unique(bench$n_rows)) {
  b <- bench[bench$n_rows == nr, ]
  if (nrow(b) > 1) {
    rel <- diff(range(b$loglik)) / abs(mean(b$loglik))
    vc_check("benchmarks",
             sprintf("engines agree on loglik at n = %s", format(nr, big.mark = ",")),
             rel < 1e-8, sprintf("relative spread %.2e across %s",
                                 rel, paste(b$engine, collapse = "/")))
  }
}

sp <- bench[bench$engine == "sparse", ]; dn <- bench[bench$engine == "dense", ]
mrg <- merge(sp, dn, by = "n_rows", suffixes = c("_sp", "_dn"))
if (nrow(mrg)) {
  cat("\n  sparse vs dense:\n")
  for (i in seq_len(nrow(mrg))) {
    cat(sprintf("    n = %-10s  time %5.2fx faster   peak RSS %5.2fx smaller\n",
                format(mrg$n_rows[i], big.mark = ","),
                mrg$seconds_dn[i] / mrg$seconds_sp[i],
                mrg$peak_rss_gb_dn[i] / mrg$peak_rss_gb_sp[i]))
  }
}
write.csv(bench, file.path(vc_dir, "results_07_benchmarks.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
vc_header("Tolerance semantics: how max|gradient| scales with n")
# -----------------------------------------------------------------------------
# `tol` applies to max|gradient|, and the gradient is a SUM over strata. What
# that costs is NOT that the optimum moves -- at the MLE the gradient is zero
# by definition -- but that the ATTAINABLE numerical floor rises: the summation
# error in accumulating millions of per-stratum terms puts a lower bound on how
# small max|grad| can be computed to be. The kernel comments already record
# this empirically ("with 75M+ rows, max|grad| of ~0.003 can be the numerical
# floor"), and it is the reason the convergence ladder needed tiers beyond the
# primary gradient test at all.
#
# Measured directly: fit at each scale with a tolerance far below anything
# achievable and record the max|grad| actually reached, plus the route taken.
# Reported, not acted on -- changing the convergence criterion under three live
# papers is not a pre-patch move.

scales <- c(2000L, 10000L, 50000L, 200000L)
floor_grad <- numeric(length(scales)); route <- character(length(scales))
iters <- integer(length(scales))
for (i in seq_along(scales)) {
  set.seed(7)
  ns <- scales[i]; na <- 10L; n <- ns * na
  X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  b <- c(0.8, -0.5)
  eta <- as.vector(X %*% b)
  ch <- integer(n); st <- rep(seq_len(ns), each = na)
  for (s in split(seq_len(n), st))
    ch[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  f <- fastclogit(X, ch, st, tol = 1e-14, max_iter = 200L)
  floor_grad[i] <- max(abs(f$gradient))
  route[i] <- f$convergence_criterion
  iters[i] <- f$iterations
}
tg <- data.frame(n_strata = scales, attainable_max_grad = floor_grad,
                 route = route, iters = iters)
print(format(tg, digits = 4))
# WHAT THIS ACTUALLY SHOWED, which is not what was expected.
#
# The hypothesis was that the attainable gradient floor RISES with n, making
# an absolute tol progressively harder to reach and explaining why the
# convergence ladder needed extra tiers. At these scales that is simply not
# true: the floor measured 1e-13 to 1e-9 over 2k to 200k strata, is not
# monotone in n, and sits far below the 1e-6 default. Fitting an exponent to
# four non-monotone points and extrapolating it to 37M strata would be noise
# presented as a trend, so no exponent is reported.
#
# The kernel comments record a floor of ~0.003 at 75M+ rows. That is three
# orders of magnitude beyond anything reproducible on this machine, so it can
# be neither confirmed nor refuted here, and the documentation should attribute
# it to production experience rather than to this study.
#
# What the run DID surface is a usability finding worth keeping: when tol is
# set below the attainable floor, the optimiser grinds to max_iter rather than
# recognising that it is finished. Two of the four fits ran all 200 iterations
# with final gradients of 3.8e-13 and 1.1e-11 and were reported as
# NOT CONVERGED. The plateau rule did not rescue them, because its
# side-condition requires evidence of struggling (halvings, or a tiny step)
# that a cleanly-converged fit never produces.
n_capped <- sum(route == "iter_max")
vc_note(sprintf("attainable floor spans %.1e to %.1e over %s to %s strata; not monotone",
                min(floor_grad), max(floor_grad),
                format(min(scales), big.mark = ","), format(max(scales), big.mark = ",")))
vc_note("hypothesis NOT supported at these scales: the floor is far below tol = 1e-6")
vc_note(sprintf("%d of %d fits ground to max_iter under an unreachable tol",
                n_capped, length(scales)))
vc_check("tol-semantics",
         "tol = 1e-6 is comfortably attainable up to 200k strata",
         all(floor_grad < 1e-6),
         sprintf("max attainable floor %.1e", max(floor_grad)))
vc_check("tol-semantics",
         "an unreachable tol is reported as failure, not recognised [finding]",
         n_capped == 0L,
         sprintf("%d/%d fits hit max_iter with max|grad| < 1e-10",
                 n_capped, length(scales)))

write.csv(tg, file.path(vc_dir, "results_07_tolscaling.csv"), row.names = FALSE)
