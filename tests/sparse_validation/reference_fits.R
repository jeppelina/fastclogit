# =============================================================================
# reference_fits.R -- Build reference fits with fastclogit (dense) and
# survival::clogit() for every simulation problem.
#
# Results are cached under tests/sparse_validation/reference/<name>.rds.
# Re-running is idempotent: existing RDS with matching seed are reused.
#
# Each cached RDS holds:
#   list(
#     name, seed, data,                   # data = the list from gen_data
#     fastclogit_fit = list(
#       coefficients, vcov, vcov_robust, se, se_robust,
#       loglik, gradient, iterations, converged,
#       time_sec, n_obs, n_groups, n_clusters
#     ),
#     survival_fit = list( ... ) or NULL   # NULL if survival skipped
#   )
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(survival)
  library(fastclogit)
})

# Locate this file's directory robustly (works for source() and Rscript)
.this_dir <- (function() {
  # 1. sourced via source(): use sys.frames
  for (i in seq_len(sys.nframe())) {
    fr <- sys.frame(i)
    if (!is.null(fr$ofile)) return(normalizePath(dirname(fr$ofile)))
  }
  # 2. Rscript: parse commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grepl("^--file=", args)]
  if (length(fa) > 0) return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  # 3. fallback to cwd-relative
  normalizePath(file.path("tests", "sparse_validation"))
})()

source(file.path(.this_dir, "gen_data.R"))

REF_DIR <- file.path(.this_dir, "reference")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

# ----- helpers ---------------------------------------------------------------

.section <- function(txt) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_h1(txt)
  } else {
    cat("\n", strrep("=", 70), "\n", txt, "\n", strrep("=", 70), "\n", sep = "")
  }
}
.note <- function(txt) {
  if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_info(txt)
  else cat("[i] ", txt, "\n", sep = "")
}

.fit_fastclogit <- function(data) {
  t0 <- Sys.time()
  fit <- fastclogit::fastclogit(
    X       = data$X,
    choice  = data$choice,
    strata  = data$strata,
    cluster = data$cluster,
    max_iter = 100L,
    tol     = 1e-8,
    verbose = FALSE
  )
  t1 <- Sys.time()
  list(
    coefficients = fit$coefficients,
    vcov         = fit$vcov,
    vcov_robust  = fit$vcov_robust,
    se           = fit$se,
    se_robust    = fit$se_robust,
    loglik       = fit$loglik,
    gradient     = fit$gradient,
    iterations   = fit$iterations,
    converged    = fit$converged,
    n_obs        = fit$n_obs,
    n_groups     = fit$n_groups,
    n_clusters   = fit$n_clusters,
    time_sec     = as.numeric(difftime(t1, t0, units = "secs"))
  )
}

.fit_survival <- function(data, time_budget_sec = 600) {
  # survival::clogit accepts a 0/1 response directly (no Surv() wrapper needed).
  # `robust = TRUE` + `cluster = .cluster` gives Lin-Wei sandwich SEs, BUT
  # survival::clogit forbids robust+exact; for robust we use method="approximate"
  # (the Prentice approximation). Without cluster we use method="exact".
  df <- as.data.frame(data$X)
  # Use safe response name that won't collide with X column names
  df$.choice <- data$choice
  df$.strata <- data$strata
  if (!is.null(data$cluster)) df$.cluster <- data$cluster

  rhs <- paste(sprintf("`%s`", colnames(data$X)), collapse = " + ")

  t0 <- Sys.time()
  if (!is.null(data$cluster)) {
    fml <- stats::as.formula(sprintf(
      ".choice ~ %s + strata(.strata)", rhs
    ))
    fit <- tryCatch(
      survival::clogit(fml, data = df, method = "approximate",
                       cluster = .cluster, robust = TRUE),
      error = function(e) {
        message("  survival::clogit ERROR: ", conditionMessage(e)); NULL
      }
    )
  } else {
    fml <- stats::as.formula(sprintf(
      ".choice ~ %s + strata(.strata)", rhs
    ))
    fit <- tryCatch(
      survival::clogit(fml, data = df, method = "exact"),
      error = function(e) {
        message("  survival::clogit ERROR: ", conditionMessage(e)); NULL
      }
    )
  }
  t1 <- Sys.time()
  if (is.null(fit)) return(NULL)

  cf  <- stats::coef(fit)
  # survival keeps backticks in names when columns are backtick-quoted in RHS;
  # strip them so names align with fastclogit's coefficient names.
  names(cf) <- gsub("`", "", names(cf), fixed = TRUE)
  nms <- names(cf)
  V_model <- if (is.null(data$cluster)) stats::vcov(fit) else (fit$naive.var %||% stats::vcov(fit))
  V_robust <- if (!is.null(data$cluster)) fit$var else NULL
  if (!is.null(V_model)) { rownames(V_model) <- colnames(V_model) <- nms }
  if (!is.null(V_robust)) { rownames(V_robust) <- colnames(V_robust) <- nms }
  list(
    coefficients = cf,
    vcov         = V_model,
    vcov_robust  = V_robust,
    se           = setNames(sqrt(diag(V_model)), nms),
    se_robust    = if (!is.null(V_robust)) setNames(sqrt(diag(V_robust)), nms) else NULL,
    loglik       = as.numeric(stats::logLik(fit)),
    iterations   = fit$iter[1L],
    converged    = isTRUE(fit$info < 1) || is.null(fit$info),
    time_sec     = as.numeric(difftime(t1, t0, units = "secs"))
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.refit_needed <- function(rds_path, seed) {
  if (!file.exists(rds_path)) return(TRUE)
  cached <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(cached)) return(TRUE)
  if (!identical(cached$seed, seed)) return(TRUE)
  FALSE
}

.build_one <- function(name, generator, seed, force = FALSE,
                       skip_survival_for = character(0)) {
  rds_path <- file.path(REF_DIR, paste0(name, ".rds"))

  if (!force && !.refit_needed(rds_path, seed)) {
    .note(sprintf("[%s] using cached reference (%s)", name, rds_path))
    return(invisible(readRDS(rds_path)))
  }

  .section(sprintf("Building reference for problem '%s' (seed=%d)", name, seed))
  .note("simulating data ...")
  data <- generator(seed)
  .note(sprintf("n_rows=%s  p=%d  density=%.3f%%",
                format(nrow(data$X), big.mark = ","), ncol(data$X),
                100 * data$meta$density))

  .note("fitting fastclogit (dense) ...")
  fc_fit <- .fit_fastclogit(data)
  .note(sprintf("  -> %.2fs, ll=%.4f, iter=%d, converged=%s",
                fc_fit$time_sec, fc_fit$loglik, fc_fit$iterations, fc_fit$converged))

  sv_fit <- NULL
  if (!(name %in% skip_survival_for)) {
    .note("fitting survival::clogit ...")
    sv_fit <- .fit_survival(data)
    if (!is.null(sv_fit)) {
      .note(sprintf("  -> %.2fs, ll=%.4f, iter=%s",
                    sv_fit$time_sec, sv_fit$loglik,
                    as.character(sv_fit$iterations)))
    } else {
      .note("  survival::clogit failed; recording NULL")
    }
  } else {
    .note("skipping survival::clogit for this problem (too large)")
  }

  out <- list(
    name = name, seed = seed,
    data = data,
    fastclogit_fit = fc_fit,
    survival_fit = sv_fit
  )
  saveRDS(out, rds_path)
  .note(sprintf("saved -> %s", rds_path))
  invisible(out)
}

# ----- main builder ----------------------------------------------------------

build_all_references <- function(force = FALSE) {
  # survival::clogit is O(n_alt^k) for the exact method; the paper3_like
  # problem (1M rows, 50 alters, ~100 cols) is infeasible. Skip it.
  skip_survival_for <- c("paper3_like")

  seeds <- c(small_dense = 1L, medium_factor = 2L,
             paper3_like = 3L, edge_rare_cells = 4L)

  results <- list()
  for (nm in names(ALL_GENERATORS)) {
    results[[nm]] <- .build_one(
      name = nm, generator = ALL_GENERATORS[[nm]],
      seed = seeds[[nm]], force = force,
      skip_survival_for = skip_survival_for
    )
  }

  # ----- summary table -----
  .section("Reference fits summary")
  rows <- list()
  for (nm in names(results)) {
    r <- results[[nm]]
    data <- r$data
    n_obs <- nrow(data$X); n_strata <- length(unique(data$strata)); p <- ncol(data$X)

    fc <- r$fastclogit_fit
    rows[[length(rows) + 1L]] <- data.frame(
      problem = nm, fitter = "fastclogit",
      n_obs = n_obs, n_strata = n_strata, p = p,
      time_sec = round(fc$time_sec, 3),
      logLik   = round(fc$loglik, 4),
      converged = fc$converged,
      stringsAsFactors = FALSE
    )
    sv <- r$survival_fit
    if (!is.null(sv)) {
      rows[[length(rows) + 1L]] <- data.frame(
        problem = nm, fitter = "survival",
        n_obs = n_obs, n_strata = n_strata, p = p,
        time_sec = round(sv$time_sec, 3),
        logLik   = round(sv$loglik, 4),
        converged = sv$converged,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        problem = nm, fitter = "survival",
        n_obs = n_obs, n_strata = n_strata, p = p,
        time_sec = NA_real_, logLik = NA_real_, converged = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  tbl <- do.call(rbind, rows)
  print(tbl, row.names = FALSE)
  invisible(list(results = results, summary = tbl))
}

# When sourced interactively from package root, allow `Rscript reference_fits.R`
if (sys.nframe() == 0L) {
  build_all_references(force = FALSE)
}
