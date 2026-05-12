# Stress-test sparse vs dense at Paper-3-like scale.
#
# Each fit runs in a fresh R subprocess so /usr/bin/time -l can measure
# its peak RSS independently. Driver collects results into a CSV.
#
# Usage: Rscript tests/sparse_validation/stress_bench.R [SCALE]
#   small   (default; finishes in <2 min, ≤4 GB peak per fit)
#   medium  (~5 min, ≤25 GB peak)
#   paper3  (the actual MONA target, may need 100+ GB on dense)

suppressPackageStartupMessages({
  library(fastclogit); library(Matrix); library(data.table)
})

# ---------------------------------------------------------------------------
# Data generator — defined at TOP so both driver and worker can use it.
# Same structure as gen_paper3_like (factor pair × decade × age × edu + dist).
# ---------------------------------------------------------------------------
gen_paper3_like_scalable <- function(n_strata, K = 30L, seed = 999L) {
  set.seed(seed)
  n <- n_strata * K
  pair_levels   <- paste0("P", sprintf("%02d", 1:25))
  decade_levels <- c("d1990", "d2000", "d2010")
  age_levels    <- c("a1", "a2", "a3", "a4")
  edu_levels    <- c("e1", "e2", "e3")

  pair_w <- c(0.70, rep(0.30 / 24, 24))  # 70% Sweden-Sweden (mimics endogamy)
  pair   <- sample(pair_levels, n, replace = TRUE, prob = pair_w)
  decade <- sample(decade_levels, n, replace = TRUE)
  age    <- sample(age_levels, n, replace = TRUE)
  edu    <- sample(edu_levels, n, replace = TRUE)
  dist   <- rnorm(n)

  df <- data.frame(pair = factor(pair, levels = pair_levels),
                   decade = factor(decade, levels = decade_levels),
                   age = factor(age, levels = age_levels),
                   edu = factor(edu, levels = edu_levels),
                   dist = dist)
  X <- stats::model.matrix(~ pair*decade + age*decade + edu*decade + dist*decade,
                            data = df)[, -1, drop = FALSE]
  p <- ncol(X)
  beta_true <- rnorm(p, sd = 0.3)
  util <- as.numeric(X %*% beta_true) - log(-log(stats::runif(n)))

  strata <- rep(seq_len(n_strata), each = K)
  choice <- integer(n)
  for (s in seq_len(n_strata)) {
    idx <- which(strata == s)
    choice[idx[which.max(util[idx])]] <- 1L
  }

  cluster <- rep(seq_len(n_strata %/% 5L + 1L), length.out = n_strata)[strata]
  Xsp <- as(X, "CsparseMatrix")
  list(X = X, Xsp = Xsp, choice = choice, strata = strata,
       cluster = cluster, beta_true = beta_true)
}

# ===========================================================================
# WORKER mode — runs a single fit and prints a CSV row to stdout
# ===========================================================================
if (Sys.getenv("FASTCLOGIT_BENCH_RUN") == "worker") {
  args     <- commandArgs(trailingOnly = TRUE)
  storage  <- args[1]
  n_strata <- as.integer(args[2])
  K        <- as.integer(args[3])
  seed     <- as.integer(args[4])

  d <- gen_paper3_like_scalable(n_strata = n_strata, K = K, seed = seed)
  Xfit   <- if (storage == "sparse") d$Xsp else d$X
  choice <- d$choice; strata <- d$strata; cluster <- d$cluster
  # Free the alternate-storage copy so peak RSS reflects only the chosen one
  if (storage == "sparse") d$X <- NULL else d$Xsp <- NULL
  d$beta_true <- NULL
  d$X <- NULL; d$Xsp <- if (storage == "sparse") Xfit else NULL
  rm(d); gc(verbose = FALSE); gc(verbose = FALSE)

  t0 <- Sys.time()
  fit <- fastclogit(X = Xfit, choice = choice, strata = strata,
                    cluster = cluster, max_iter = 100L, tol = 1e-8,
                    verbose = FALSE)
  t_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  density <- if (storage == "sparse")
    length(Xfit@x) / (nrow(Xfit) * ncol(Xfit)) else 1.0

  cat(sprintf(
    "RESULT,%s,%d,%d,%d,%.6f,%.3f,%.3e,%d,%s,%.6f\n",
    storage, n_strata, K, ncol(Xfit), density, t_secs,
    max(abs(fit$gradient)), fit$iterations,
    fit$converged, fit$loglik))
  q(save = "no")
}

# ===========================================================================
# DRIVER mode
# ===========================================================================
scale_arg <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "small"
scales <- switch(scale_arg,
  small  = list(n_strata = c(10000L, 50000L, 100000L), K = 30L),
  medium = list(n_strata = c(50000L, 200000L, 500000L), K = 30L),
  paper3 = list(n_strata = c(200000L, 500000L, 720000L), K = 50L),
  stop("Unknown scale: ", scale_arg))

cat(strrep("=", 70), "\n", sep="")
cat("Sparse-vs-dense stress benchmark (scale = ", scale_arg, ")\n", sep="")
cat(strrep("=", 70), "\n\n", sep="")

results <- list()
for (n_strata in scales$n_strata) {
  for (storage in c("sparse", "dense")) {
    cat(sprintf("[%s]  n_strata=%-6d  K=%2d  ... ", storage, n_strata, scales$K))
    cmd <- sprintf(
      "/usr/bin/time -l env FASTCLOGIT_BENCH_RUN=worker Rscript tests/sparse_validation/stress_bench.R %s %d %d 42 2>&1",
      storage, n_strata, scales$K)
    full_out <- system(cmd, intern = TRUE)
    result_line <- grep("^RESULT,", full_out, value = TRUE)
    rss_line    <- grep("maximum resident set size", full_out, value = TRUE)
    if (length(result_line) == 0L || length(rss_line) == 0L) {
      cat("FAILED\n")
      cat("  Tail of output:\n")
      cat(paste0("    ", tail(full_out, 8), collapse = "\n"), "\n")
      next
    }
    rss_bytes <- as.numeric(sub("^\\s*([0-9]+)\\s.*", "\\1", rss_line))
    rss_gb    <- rss_bytes / 1024^3
    pieces    <- strsplit(result_line, ",")[[1]]
    res <- list(storage   = pieces[2],
                n_strata  = as.integer(pieces[3]),
                K         = as.integer(pieces[4]),
                p         = as.integer(pieces[5]),
                density   = as.numeric(pieces[6]),
                time_sec  = as.numeric(pieces[7]),
                max_grad  = as.numeric(pieces[8]),
                iters     = as.integer(pieces[9]),
                converged = pieces[10] == "TRUE",
                loglik    = as.numeric(pieces[11]),
                peak_gb   = rss_gb)
    cat(sprintf("time=%6.1fs  peak=%5.2f GB  iters=%d  density=%.1f%%  ll=%.2f\n",
                res$time_sec, res$peak_gb, res$iters,
                100*res$density, res$loglik))
    results[[length(results) + 1L]] <- res
  }
  cat("\n")
}

if (length(results) == 0L) { cat("No successful runs.\n"); q(status = 1) }

dt <- rbindlist(results)
cat("--- Final results ---\n")
print(dt)

# Pivot for sparse vs dense comparison
if (length(unique(dt$storage)) == 2L) {
  cat("\n--- Comparison table ---\n")
  d_dense  <- dt[storage == "dense",  .(n_strata, p, time_dense=time_sec, peak_dense=peak_gb, ll_dense=loglik, iters_dense=iters)]
  d_sparse <- dt[storage == "sparse", .(n_strata, p, time_sparse=time_sec, peak_sparse=peak_gb, ll_sparse=loglik, iters_sparse=iters)]
  pivot <- merge(d_dense, d_sparse, by = c("n_strata", "p"), all = TRUE)
  pivot[, mem_savings_x := peak_dense / peak_sparse]
  pivot[, time_ratio    := time_sparse / time_dense]
  pivot[, ll_diff       := abs(ll_sparse - ll_dense)]
  print(pivot)
}

out_csv <- sprintf("tests/sparse_validation/bench_%s.csv", scale_arg)
fwrite(dt, out_csv)
cat("\nSaved CSV: ", out_csv, "\n", sep="")
