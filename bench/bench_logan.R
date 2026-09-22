#!/usr/bin/env Rscript
# ============================================================================
# bench_logan.R: Speed & memory benchmark: fastclogit vs survival::clogit
#
# Two benchmarks using the Logan occupational mobility data:
#   Part 1, SPEED: fit the same model N times, compare wall-clock time
#   Part 2, MEMORY: single fit on increasingly large data, compare peak RSS
#
# Usage:
#   Rscript bench/bench_logan.R              # default 1000 speed fits
#   Rscript bench/bench_logan.R 500          # 500 speed fits
#
# Requirements: fastclogit (installed or devtools::load_all), survival
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
})

# --- Set up logging: mirror all output to a timestamped file ----------------
# Find script directory robustly
script_dir <- tryCatch({
  # Works with Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
}, error = function(e) getwd())

log_file <- file.path(
  script_dir,
  paste0("bench_logan_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)   # stdout to file + console
sink(log_con, type = "message", append = TRUE)  # stderr to same file

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message("Log saved to: ", log_file)
}, add = TRUE)

cat("Log: ", log_file, "\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:    ", R.version.string, "\n\n")

# --- Try to load fastclogit -------------------------------------------------
tryCatch(library(fastclogit), error = function(e) {
  # Package not installed; try devtools::load_all from parent of script dir
  pkg_root <- dirname(script_dir)
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    stop("fastclogit not installed and devtools not available. ",
         "Install with: devtools::install() from the package root.")
  }
})

# ============================================================================
# 1. Prepare the Logan data in conditional logit (long) format
# ============================================================================
data(logan, package = "survival")

# logan is one-row-per-person: occupation, foession, education, race.
# Expand to choice format: each person gets one row per occupation level,
# case = 1 for their actual occupation, 0 for the alternatives.

logan_dt   <- as.data.frame(logan)
logan_dt$id <- seq_len(nrow(logan_dt))
occ_levels <- sort(unique(logan_dt$occupation))
n_alts     <- length(occ_levels)

# Build long format
rows <- vector("list", nrow(logan_dt))
for (i in seq_len(nrow(logan_dt))) {
  rows[[i]] <- data.frame(
    id        = logan_dt$id[i],
    alt_occ   = occ_levels,
    case      = as.integer(occ_levels == logan_dt$occupation[i]),
    education = logan_dt$education[i],
    race_black = as.integer(logan_dt$race[i] == "black"),
    stringsAsFactors = FALSE
  )
}
logan_long <- do.call(rbind, rows)

# Build identified design matrix for conditional logit.
# Person-level covariates only enter via interaction with alt dummies.
# Drop first alternative as reference.
ref_alt <- occ_levels[1]
non_ref <- occ_levels[-1]

mm_list <- list()
for (lev in non_ref) {
  d <- as.integer(logan_long$alt_occ == lev)
  tag <- gsub("[^A-Za-z0-9]", "", lev)
  mm_list[[paste0("d_", tag)]]    <- d
  mm_list[[paste0("edu_", tag)]]  <- d * logan_long$education
  mm_list[[paste0("race_", tag)]] <- d * logan_long$race_black
}
X      <- do.call(cbind, mm_list)
choice <- logan_long$case
strata <- logan_long$id

n_persons <- length(unique(strata))

cat("=== Logan Occupational Mobility Benchmark ===\n\n")
cat("Persons:      ", nrow(logan_dt), "\n", sep = "")
cat("Alternatives: ", n_alts, " (", paste(occ_levels, collapse = ", "), ")\n", sep = "")
cat("Reference:    ", ref_alt, "\n", sep = "")
cat("Long-format:  ", nrow(logan_long), " rows, ", n_persons, " choice sets\n", sep = "")
cat("Parameters:   ", ncol(X), "\n\n", sep = "")

# Also prepare a data.frame for survival::clogit
clogit_df <- data.frame(X, case = choice, id = strata)
clogit_fml <- as.formula(paste("case ~",
                               paste(colnames(X), collapse = " + "),
                               "+ strata(id)"))

# ============================================================================
# 2. Benchmark
# ============================================================================
args   <- commandArgs(trailingOnly = TRUE)
n_fits <- if (length(args) > 0) as.integer(args[1]) else 1000L

cat("Number of fits: ", n_fits, "\n\n", sep = "")

# --- Warm-up (compile, JIT, cache) -----------------------------------------
invisible(fastclogit(X, choice, strata, max_iter = 25, tol = 1e-8, verbose = FALSE))
invisible(clogit(clogit_fml, data = clogit_df, method = "efron"))

# --- fastclogit -------------------------------------------------------------
cat("Running fastclogit x", n_fits, "...\n")
# Helper: peak Vcells MB from gc(), R 4.1+ has 7 columns, older has 6
peak_mb <- function() {
  g <- gc(verbose = FALSE)
  # "max used (Mb)" is the last column, row 2 (Vcells)
  g[2, ncol(g)]
}

gc(reset = TRUE, verbose = FALSE)
mem_before_fc <- peak_mb()

t_fc <- system.time({
  for (i in seq_len(n_fits)) {
    fit_fc <- fastclogit(X, choice, strata, max_iter = 25, tol = 1e-8, verbose = FALSE)
  }
})

mem_after_fc <- peak_mb()

cat(sprintf("  Total:    %6.2f sec\n", t_fc[3]))
cat(sprintf("  Per fit:  %6.4f sec\n", t_fc[3] / n_fits))
cat(sprintf("  Memory:   %6.1f MB peak\n", max(mem_after_fc - mem_before_fc, 0)))
cat(sprintf("  Converged: %s (%d iters)\n\n", fit_fc$converged, fit_fc$iterations))

# --- survival::clogit -------------------------------------------------------
cat("Running survival::clogit x", n_fits, "...\n")
gc(reset = TRUE, verbose = FALSE)
mem_before_cl <- peak_mb()

t_cl <- system.time({
  for (i in seq_len(n_fits)) {
    fit_cl <- clogit(clogit_fml, data = clogit_df, method = "efron")
  }
})

mem_after_cl <- peak_mb()

cat(sprintf("  Total:    %6.2f sec\n", t_cl[3]))
cat(sprintf("  Per fit:  %6.4f sec\n", t_cl[3] / n_fits))
cat(sprintf("  Memory:   %6.1f MB peak\n\n", max(mem_after_cl - mem_before_cl, 0)))

# ============================================================================
# 3. Comparison
# ============================================================================
speedup <- t_cl[3] / max(t_fc[3], 0.001)

# Coefficient agreement
coefs_fc <- fit_fc$coefficients
coefs_cl <- coef(fit_cl)
common   <- intersect(names(coefs_fc), names(coefs_cl))
max_diff <- if (length(common) > 0) max(abs(coefs_fc[common] - coefs_cl[common])) else NA

cat("=== Results ===\n\n")
cat(sprintf("  Speedup:          %.1fx\n", speedup))
cat(sprintf("  Max |coef diff|:  %.2e\n", max_diff))
cat(sprintf("  fastclogit total: %.2f sec (%.4f sec/fit)\n", t_fc[3], t_fc[3] / n_fits))
cat(sprintf("  clogit total:     %.2f sec (%.4f sec/fit)\n", t_cl[3], t_cl[3] / n_fits))

# Print coefficients side by side
if (length(common) > 0) {
  cat("\n  Coefficients:\n")
  cat(sprintf("  %-20s %12s %12s %12s\n", "Parameter", "fastclogit", "clogit", "diff"))
  cat(sprintf("  %-20s %12s %12s %12s\n", "---------", "----------", "------", "----"))
  for (nm in common) {
    cat(sprintf("  %-20s %12.6f %12.6f %12.2e\n",
                nm, coefs_fc[nm], coefs_cl[nm], coefs_fc[nm] - coefs_cl[nm]))
  }
}

# ============================================================================
# 4. MEMORY BENCHMARK, single fit on large data
# ============================================================================
cat("\n=== Part 2: Memory Benchmark (large data) ===\n\n")

# Helper: replicate the logan long-format data N times with unique strata
replicate_data <- function(X_base, choice_base, strata_base, n_ids, n_reps) {
  n_base <- nrow(X_base)
  n_total <- n_base * n_reps

  X_big      <- matrix(0, nrow = n_total, ncol = ncol(X_base))
  choice_big <- integer(n_total)
  strata_big <- integer(n_total)
  colnames(X_big) <- colnames(X_base)

  for (r in seq_len(n_reps)) {
    idx <- ((r - 1L) * n_base + 1L):(r * n_base)
    X_big[idx, ]    <- X_base
    choice_big[idx] <- choice_base
    strata_big[idx] <- strata_base + (r - 1L) * n_ids
  }

  list(X = X_big, choice = choice_big, strata = strata_big,
       n_rows = n_total, n_groups = n_ids * n_reps)
}

# Scale levels: number of replications of the base 4190-row dataset
mem_levels <- c(10, 50, 100, 500)

cat(sprintf("Base: %d rows x %d cols = %.1f KB\n",
            nrow(X), ncol(X), object.size(X) / 1024))
cat(sprintf("Replication levels: %s\n\n", paste(mem_levels, collapse = ", ")))

mem_results <- list()

for (n_rep in mem_levels) {
  dat <- replicate_data(X, choice, strata, n_persons, n_rep)

  cat(sprintf("--- %dx: %s rows, %s choice sets ---\n",
              n_rep, format(dat$n_rows, big.mark = ","),
              format(dat$n_groups, big.mark = ",")))

  input_mb <- as.numeric(object.size(dat$X)) / 1024^2

  # --- fastclogit: measure memory overhead beyond the input data ---
  gc(verbose = FALSE); gc(reset = TRUE, verbose = FALSE)
  mem_before <- peak_mb()

  t_fc <- system.time({
    fit_fc <- fastclogit(dat$X, dat$choice, dat$strata,
                         max_iter = 200, tol = 1e-6, verbose = FALSE)
  })

  mem_after <- peak_mb()
  fc_overhead <- max(mem_after - mem_before, 0)

  cat(sprintf("  fastclogit:  %5.2f sec | overhead %6.1f MB | input %5.1f MB | conv: %s (%d it)\n",
              t_fc[3], fc_overhead, input_mb, fit_fc$converged, fit_fc$iterations))

  row <- data.frame(
    reps       = n_rep,
    n_rows     = dat$n_rows,
    n_groups   = dat$n_groups,
    input_mb   = round(input_mb, 1),
    fc_sec     = round(t_fc[3], 3),
    fc_over_mb = round(fc_overhead, 1),
    stringsAsFactors = FALSE
  )

  # --- survival::clogit ---
  bench_df <- data.frame(dat$X, case = dat$choice, id = dat$strata)

  gc(verbose = FALSE); gc(reset = TRUE, verbose = FALSE)
  mem_before <- peak_mb()

  t_cl <- system.time({
    fit_cl <- clogit(clogit_fml, data = bench_df, method = "efron")
  })

  mem_after <- peak_mb()
  cl_overhead <- max(mem_after - mem_before, 0)

  cat(sprintf("  clogit:      %5.2f sec | overhead %6.1f MB\n",
              t_cl[3], cl_overhead))

  row$cl_sec     <- round(t_cl[3], 3)
  row$cl_over_mb <- round(cl_overhead, 1)
  row$mem_ratio  <- round(cl_overhead / max(fc_overhead, 0.1), 1)
  row$speedup    <- round(t_cl[3] / max(t_fc[3], 0.001), 1)

  rm(bench_df, fit_cl)

  mem_results[[length(mem_results) + 1]] <- row

  rm(dat, fit_fc)
  gc(verbose = FALSE)
  cat("\n")
}

# ============================================================================
# 5. Summary
# ============================================================================
mem_dt <- do.call(rbind, mem_results)

cat("=== Memory Benchmark Summary ===\n\n")
cat("  'overhead' = peak memory allocated during fit, beyond the input data.\n")
cat("  'mem_ratio' = clogit overhead / fastclogit overhead.\n\n")
print(mem_dt, row.names = FALSE)

cat("\n=== Overall ===\n\n")
cat(sprintf("  Speed (Part 1, %d fits):   fastclogit %.1fx faster\n", n_fits, speedup))
has_mem <- !is.na(mem_dt$mem_ratio)
if (any(has_mem)) {
  cat(sprintf("  Memory (Part 2, largest):  clogit uses %.1fx more overhead\n",
              max(mem_dt$mem_ratio[has_mem])))
}

cat("\n=== Session Info ===\n\n")
print(sessionInfo())
cat("\nDone.\n")
