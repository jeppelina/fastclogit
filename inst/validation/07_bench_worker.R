# =============================================================================
# 07_bench_worker.R: ONE fit, in its own process, so peak RSS is attributable.
#
# Run via: /usr/bin/time -l Rscript 07_bench_worker.R <engine> <n_sets> <n_alts> <p_fac>
# where engine is one of: dense, sparse, survival
#
# Peak memory MUST be measured at the OS level. gc() reports R's heap;
# Armadillo allocates outside it, so the design matrix, the CSR view and every
# kernel workspace are invisible to gc(). Benchmarking this package's memory
# with gc() would measure everything except the thing being claimed.
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
engine <- args[1]
n_sets <- as.integer(args[2]); n_alts <- as.integer(args[3])
p_fac  <- as.integer(args[4])

suppressPackageStartupMessages({
  library(fastclogit); library(Matrix)
})

set.seed(42)
n <- n_sets * n_alts
d <- data.frame(f = factor(sample.int(p_fac, n, TRUE)), x = rnorm(n))
strata <- rep(seq_len(n_sets), each = n_alts)

X <- Matrix::sparse.model.matrix(~ f + x, data = d)[, -1, drop = FALSE]
beta <- c(rnorm(ncol(X) - 1, 0, 0.5), -0.4)
eta <- as.vector(X %*% beta)
choice <- integer(n)
for (s in split(seq_len(n), strata))
  choice[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L

t0 <- Sys.time()
if (engine == "dense") {
  Xd <- as.matrix(X)
  fit <- fastclogit(Xd, choice, strata, max_iter = 100L)
  ll <- fit$loglik; it <- fit$iterations
} else if (engine == "sparse") {
  fit <- fastclogit(X, choice, strata, max_iter = 100L)
  ll <- fit$loglik; it <- fit$iterations
} else if (engine == "survival") {
  suppressPackageStartupMessages(library(survival))
  dd <- as.data.frame(as.matrix(X)); dd$choice <- choice; dd$sid <- strata
  nm <- setdiff(names(dd), c("choice", "sid"))
  names(dd)[seq_along(nm)] <- paste0("V", seq_along(nm))
  fml <- stats::as.formula(paste("choice ~",
                                 paste(paste0("V", seq_along(nm)), collapse = " + "),
                                 "+ strata(sid)"))
  fit <- clogit(fml, data = dd, method = "efron")
  ll <- as.numeric(fit$loglik[2]); it <- NA_integer_
}
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("RESULT\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.6f\t%s\n",
            engine, n, ncol(X), n_sets, n_alts, el, ll, as.character(it)))
