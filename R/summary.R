#' @export
coef.fastclogit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.fastclogit <- function(object, robust = TRUE, ...) {
  if (robust && !is.null(object$vcov_robust)) {
    object$vcov_robust
  } else {
    object$vcov
  }
}

#' @export
logLik.fastclogit <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(object$coefficients)
  attr(ll, "nobs") <- object$n_obs
  class(ll) <- "logLik"
  ll
}

#' @export
confint.fastclogit <- function(object, parm = NULL, level = 0.95,
                                robust = TRUE, ...) {
  cf <- object$coefficients
  se <- if (robust && !is.null(object$se_robust)) object$se_robust else object$se
  z <- qnorm((1 + level) / 2)

  ci <- cbind(cf - z * se, cf + z * se)
  colnames(ci) <- paste0(round(100 * c((1 - level)/2, (1 + level)/2), 1), " %")
  rownames(ci) <- names(cf)

  if (!is.null(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}

#' @export
print.fastclogit <- function(x, ...) {
  cat("fastclogit fit\n")
  cat("  Observations:", format(x$n_obs, big.mark = ","), "\n")
  cat("  Choice sets: ", format(x$n_groups, big.mark = ","), "\n")
  if (!is.null(x$n_clusters)) {
    cat("  Clusters:    ", format(x$n_clusters, big.mark = ","), "\n")
  }
  cat("  Log-lik:     ", format(x$loglik, digits = 8), "\n")
  cat("  Converged:   ", x$converged, " (", x$iterations, " iterations)\n", sep = "")
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 6))
  invisible(x)
}

#' @export
summary.fastclogit <- function(object, robust = TRUE, ...) {
  cf <- object$coefficients
  se <- if (robust && !is.null(object$se_robust)) object$se_robust else object$se
  se_type <- if (robust && !is.null(object$se_robust)) "robust" else "model"

  z <- cf / se
  p <- 2 * pnorm(-abs(z))

  ci <- confint(object, robust = robust)

  # Note: cbind(Name = named_vec) can silently drop column names in R 4.5.x.
  coef_table <- cbind(cf, se, z, p, ci)
  colnames(coef_table)[1:4] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  out <- list(
    coef_table = coef_table,
    loglik     = object$loglik,
    n_obs      = object$n_obs,
    n_groups   = object$n_groups,
    n_clusters = object$n_clusters,
    converged  = object$converged,
    iterations = object$iterations,
    se_type    = se_type,
    call       = object$call
  )
  class(out) <- "summary.fastclogit"
  out
}

#' @export
print.summary.fastclogit <- function(x, digits = 4, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Conditional Logit (fastclogit)\n")
  cat("  Observations: ", format(x$n_obs, big.mark = ","), "\n")
  cat("  Choice sets:  ", format(x$n_groups, big.mark = ","), "\n")
  if (!is.null(x$n_clusters)) {
    cat("  Clusters:     ", format(x$n_clusters, big.mark = ","),
        " (", x$se_type, " SEs)\n", sep = "")
  }
  cat("  Log-likelihood: ", format(x$loglik, digits = 8), "\n")
  cat("  Converged:    ", x$converged, " (", x$iterations, " iterations)\n\n", sep = "")

  cat("Coefficients:\n")
  printCoefmat(x$coef_table, digits = digits, signif.stars = TRUE,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1:2)
  cat("---\n")
  cat("SE type:", x$se_type, "\n")
  invisible(x)
}


#' Extract results as a tidy data.frame
#'
#' Returns coefficients, SEs, z-values, p-values, and CIs in a data.frame.
#' Compatible with broom-style workflows.
#'
#' Includes a safety check: if \code{names(object$coefficients)} is
#' \code{NULL} (can happen when post-fit operations like unscaling strip
#' names), falls back to \code{object$terms}, then to V1, V2, ...
#'
#' @param object A \code{fastclogit} object.
#' @param robust Use robust (clustered sandwich) SEs if available? Default
#'   \code{TRUE}.
#' @param conf.level Confidence level for intervals. Default 0.95.
#' @param exponentiate Return odds ratios instead of log-odds? Default
#'   \code{FALSE}.
#' @return A data.frame with columns: \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic}, \code{p.value}, \code{conf.low},
#'   \code{conf.high}.
#' @export
tidy_fastclogit <- function(object, robust = TRUE, conf.level = 0.95,
                             exponentiate = FALSE) {
  cf <- object$coefficients
  se <- if (robust && !is.null(object$se_robust)) object$se_robust else object$se

  # Safety: ensure names are present (they can be lost during unscaling)
  if (is.null(names(cf))) {
    if (!is.null(object$terms)) {
      names(cf) <- object$terms
    } else {
      names(cf) <- paste0("V", seq_along(cf))
    }
  }
  if (is.null(names(se))) names(se) <- names(cf)

  z <- cf / se
  p <- 2 * pnorm(-abs(z))
  zq <- qnorm((1 + conf.level) / 2)

  out <- data.frame(
    term      = names(cf),
    estimate  = unname(cf),
    std.error = unname(se),
    statistic = unname(z),
    p.value   = unname(p),
    conf.low  = unname(cf - zq * se),
    conf.high = unname(cf + zq * se),
    stringsAsFactors = FALSE
  )

  if (exponentiate) {
    out$estimate  <- exp(out$estimate)
    out$conf.low  <- exp(out$conf.low)
    out$conf.high <- exp(out$conf.high)
  }

  out
}
