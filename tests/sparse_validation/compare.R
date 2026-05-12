# =============================================================================
# compare.R -- Comparison utilities for two fitted clogit objects.
#
# Handles both fastclogit and survival::clogit outputs (and the lightweight
# "ref" list format saved by reference_fits.R).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ANSI helpers ---------------------------------------------------------------

.RED   <- "\033[31m"
.GREEN <- "\033[32m"
.YEL   <- "\033[33m"
.BLD   <- "\033[1m"
.OFF   <- "\033[0m"
.col   <- function(x, code) paste0(code, x, .OFF)

# Coefficient / vcov extraction ---------------------------------------------
#
# Accepts:
#   - a fastclogit object       -> $coefficients, $vcov [, $vcov_robust]
#   - a survival::clogit object -> coef(), vcov(), $var (robust)
#   - the slim list from reference_fits.R (named `coefficients`, `vcov`, etc.)

.extract_coef <- function(fit) {
  if (is.list(fit) && !is.null(fit$coefficients)) return(fit$coefficients)
  stats::coef(fit)
}
.extract_vcov <- function(fit) {
  if (is.list(fit) && !is.null(fit$vcov)) return(fit$vcov)
  stats::vcov(fit)
}
.extract_vcov_robust <- function(fit) {
  if (is.list(fit)) {
    if (!is.null(fit$vcov_robust)) return(fit$vcov_robust)
    if (!is.null(fit$var))         return(fit$var)
    return(NULL)
  }
  NULL
}
.extract_loglik <- function(fit) {
  if (is.list(fit) && !is.null(fit$loglik)) return(fit$loglik)
  as.numeric(stats::logLik(fit))
}

# Core comparison -----------------------------------------------------------

compare_fits <- function(fit_a, fit_b, name_a = "A", name_b = "B",
                         tol_coef = 1e-10, tol_se = 1e-8, tol_ll = 1e-6,
                         verbose = FALSE) {
  ca <- .extract_coef(fit_a); cb <- .extract_coef(fit_b)
  va <- .extract_vcov(fit_a); vb <- .extract_vcov(fit_b)
  vra <- .extract_vcov_robust(fit_a); vrb <- .extract_vcov_robust(fit_b)
  lla <- .extract_loglik(fit_a); llb <- .extract_loglik(fit_b)

  # Align by name (survival::clogit may report a different ordering)
  nms <- union(names(ca), names(cb))
  if (length(nms) == 0L) nms <- paste0("V", seq_along(ca))
  ca <- ca[nms]; cb <- cb[nms]

  # Prefer model-based SEs; if one side lacks them but both have robust SEs,
  # fall back to comparing robust SEs.
  use_robust <- (is.null(va) || is.null(vb)) && !is.null(vra) && !is.null(vrb)
  if (use_robust) {
    sea <- sqrt(diag(vra))[nms]; seb <- sqrt(diag(vrb))[nms]
    se_label <- "robust"
  } else {
    sea <- if (!is.null(va)) sqrt(diag(va))[nms] else rep(NA_real_, length(nms))
    seb <- if (!is.null(vb)) sqrt(diag(vb))[nms] else rep(NA_real_, length(nms))
    se_label <- "model"
  }

  diff_c <- abs(ca - cb)
  diff_s <- abs(sea - seb)
  pass_c <- !is.na(diff_c) & diff_c <= tol_coef
  # Treat fully-missing SE comparisons as "skip" (pass = NA), but if both
  # vectors are present then mismatches fail.
  both_se_present <- !is.na(sea) & !is.na(seb)
  pass_s <- ifelse(both_se_present, diff_s <= tol_se, NA)

  body <- data.table::data.table(
    term     = nms,
    est_a    = unname(ca),
    est_b    = unname(cb),
    diff     = unname(diff_c),
    pass_coef = pass_c,
    se_a     = unname(sea),
    se_b     = unname(seb),
    se_diff  = unname(diff_s),
    pass_se  = pass_s
  )

  header <- data.table::data.table(
    term     = "__logLik__",
    est_a    = lla,
    est_b    = llb,
    diff     = abs(lla - llb),
    pass_coef = abs(lla - llb) <= tol_ll,
    se_a     = NA_real_, se_b = NA_real_, se_diff = NA_real_, pass_se = NA
  )

  out <- rbind(header, body)
  attr(out, "name_a")    <- name_a
  attr(out, "name_b")    <- name_b
  attr(out, "tol_coef")  <- tol_coef
  attr(out, "tol_se")    <- tol_se
  attr(out, "tol_ll")    <- tol_ll
  attr(out, "se_label")  <- se_label

  if (verbose) format_comparison(out)
  out
}

compare_summary <- function(comparison_dt) {
  body <- comparison_dt[term != "__logLik__"]
  header <- comparison_dt[term == "__logLik__"]
  # NA pass_se = "skip" (one side missing SEs); only explicit FALSE counts.
  se_failures <- sum(body$pass_se == FALSE, na.rm = TRUE)
  list(
    all_pass       = all(body$pass_coef) && (se_failures == 0L) && isTRUE(header$pass_coef),
    max_coef_diff  = max(body$diff, na.rm = TRUE),
    max_se_diff    = suppressWarnings(max(body$se_diff, na.rm = TRUE)),
    ll_diff        = header$diff,
    n_failures     = sum(!body$pass_coef) + se_failures + (!isTRUE(header$pass_coef))
  )
}

format_comparison <- function(comparison_dt) {
  name_a <- attr(comparison_dt, "name_a") %||% "A"
  name_b <- attr(comparison_dt, "name_b") %||% "B"
  tol_c  <- attr(comparison_dt, "tol_coef") %||% 1e-10
  tol_s  <- attr(comparison_dt, "tol_se")   %||% 1e-8
  tol_ll <- attr(comparison_dt, "tol_ll")   %||% 1e-6

  cat(.col(sprintf("\n  Comparing %s vs %s\n", name_a, name_b), .BLD))
  cat(sprintf("  tol_coef=%.0e  tol_se=%.0e  tol_ll=%.0e\n", tol_c, tol_s, tol_ll))

  header <- comparison_dt[term == "__logLik__"]
  body   <- comparison_dt[term != "__logLik__"]

  ll_col <- if (isTRUE(header$pass_coef)) .col("PASS", .GREEN) else .col("FAIL", .RED)
  cat(sprintf("  logLik: %s  (%.10g vs %.10g, |diff|=%.3e)\n",
              ll_col, header$est_a, header$est_b, header$diff))

  fmt_row <- function(i) {
    row <- body[i]
    se_ok <- is.na(row$pass_se) || isTRUE(row$pass_se)
    tag <- if (row$pass_coef && se_ok) .col("OK  ", .GREEN) else .col("FAIL", .RED)
    pad <- format(row$term, width = 40)
    coef_str <- sprintf("est %s=%14.7g  %s=%14.7g  d=%9.2e",
                        name_a, row$est_a, name_b, row$est_b, row$diff)
    se_str   <- if (!is.na(row$se_a)) {
      sprintf("  se d=%9.2e", row$se_diff)
    } else ""
    cat(sprintf("  %s %s  %s%s\n", tag, pad, coef_str, se_str))
  }

  for (i in seq_len(nrow(body))) {
    se_fail <- !is.na(body$pass_se[i]) && !body$pass_se[i]
    if (!body$pass_coef[i] || se_fail) {
      fmt_row(i)
    }
  }

  s <- compare_summary(comparison_dt)
  status <- if (s$all_pass) .col("ALL PASS", .GREEN) else .col(sprintf("%d FAILURES", s$n_failures), .RED)
  cat(sprintf("  -> %s  (max |dCoef|=%.3e, max |dSE|=%.3e, |dLL|=%.3e)\n\n",
              status, s$max_coef_diff, s$max_se_diff, s$ll_diff))

  invisible(s)
}

# `%||%` polyfill
`%||%` <- function(a, b) if (is.null(a)) b else a
