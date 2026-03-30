###############################################################################
#### fastclogit_pure.R — Pure-R conditional logit (NO Rcpp required)
####
#### Self-contained fallback for environments without Rtools/Rcpp.
#### Identical interface and results to the Rcpp version, but runs anywhere
#### R is installed — no compiler toolchain needed.
####
#### Includes all components in a single file:
####   - clogit_fit_r():        Newton-Raphson with LogSumExp, adaptive ridge
####                             regularization, and step-halving line search
####   - clogit_sandwich_r():   Clustered sandwich variance (H^-1 B H^-1)
####   - fastclogit():          Main wrapper (input validation, grouping, names)
####   - S3 methods:            summary, print, vcov, confint, tidy_fastclogit
####   - simulate_clogit_data(): DGP for testing (partner choice structure)
####
#### Memory-efficient: processes group-by-group, accumulates gradient/Hessian
#### in-place. Peak memory ≈ X matrix + O(p²) workspace.
####
#### NOTE: All stats functions use explicit stats:: namespace (var, pnorm,
#### qnorm, printCoefmat, rnorm, rpois) so the file works correctly when
#### sourced into restricted environments via sys.source().
####
#### Usage on MONA (if no Rtools):
####   source("fastclogit_pure.R")
####   fit <- fastclogit(X, choice, strata, offset = off, cluster = cl)
####   summary(fit)
####
#### Author: Jesper Lindmarker
#### License: MIT
###############################################################################

# =============================================================================
# 1. NEWTON-RAPHSON CORE (replaces clogit_newton.cpp)
# =============================================================================

clogit_fit_r <- function(X, choice, offset, group_start, group_size,
                         max_iter, tol, verbose) {
  n <- nrow(X)
  p <- ncol(X)
  G <- length(group_start)

  beta   <- rep(0.0, p)
  grad   <- rep(0.0, p)
  hess   <- matrix(0.0, p, p)
  loglik <- 0.0
  converged <- FALSE
  iter <- 0L

  for (it in seq_len(max_iter)) {
    iter <- it
    loglik_new <- 0.0
    grad[] <- 0.0
    hess[] <- 0.0

    for (j in seq_len(G)) {
      start <- group_start[j] + 1L   # convert 0-based to 1-based
      K     <- group_size[j]
      idx   <- start:(start + K - 1L)

      Xj <- X[idx, , drop = FALSE]
      oj <- offset[idx]

      # Linear predictor
      eta <- as.numeric(Xj %*% beta) + oj

      # LogSumExp trick
      max_eta <- max(eta)
      exp_eta <- exp(eta - max_eta)
      sum_exp <- sum(exp_eta)
      prob    <- exp_eta / sum_exp

      # Find chosen
      c_local <- which(choice[idx] == 1L)
      if (length(c_local) == 0L) next
      c_local <- c_local[1L]

      # Log-likelihood
      loglik_new <- loglik_new + eta[c_local] - (max_eta + log(sum_exp))

      # Weighted mean: xbar = X_j' %*% prob
      xbar <- as.numeric(crossprod(Xj, prob))   # p x 1

      # Gradient: x_chosen - E[x]
      grad <- grad + Xj[c_local, ] - xbar

      # Hessian: -(E[xx'] - E[x]E[x]')
      # Efficient: X_w = Xj * sqrt(prob), then E[xx'] = X_w' %*% X_w
      sqrt_prob <- sqrt(prob)
      Xj_w <- Xj * sqrt_prob                    # K x p, recycling sqrt_prob
      hess <- hess - (crossprod(Xj_w) - tcrossprod(xbar))
    }

    if (verbose) {
      cat(sprintf("Iter %d: loglik = %.8f, max|grad| = %.2e\n",
                  it, loglik_new, max(abs(grad))))
    }

    # --- Convergence checks (BEFORE Newton step) ---
    grad_max <- max(abs(grad))

    # Primary: gradient norm
    if (grad_max < tol) {
      loglik <- loglik_new
      converged <- TRUE
      break
    }

    # Secondary: relative LL change (catches flat regions near optimum)
    if (it > 1L) {
      rel_ll <- abs(loglik_new - loglik) / (abs(loglik) + 1e-10)
      if (rel_ll < tol * 0.01 && grad_max < tol * 100) {
        loglik <- loglik_new
        converged <- TRUE
        if (verbose) {
          cat(sprintf("  Converged on relative loglik (%.2e) + gradient (%.2e)\n",
                      rel_ll, grad_max))
        }
        break
      }
    }

    loglik <- loglik_new

    # --- Newton step: delta = -H^{-1} %*% grad ---
    # Adaptive ridge regularization for near-singular Hessians
    # (rare factor categories with little within-group variation)
    neg_hess <- -hess
    ridge <- max(1e-12, 1e-8 * max(abs(diag(neg_hess))))
    diag(neg_hess) <- diag(neg_hess) + ridge

    delta <- tryCatch(
      solve(neg_hess, grad),
      error = function(e) {
        # Fallback: larger ridge
        diag(neg_hess) <<- diag(neg_hess) + 1e-4 * max(abs(diag(neg_hess)))
        tryCatch(solve(neg_hess, grad), error = function(e2) {
          warning("Hessian singular at iteration ", it)
          NULL
        })
      }
    )
    if (is.null(delta)) break

    # --- Step-halving (skip on first iteration) ---
    step_size <- 1.0
    beta_new  <- beta + step_size * delta

    if (it > 1L) {
      max_halving <- 20L
      for (h in seq_len(max_halving)) {
        ll_cand <- 0.0
        for (j in seq_len(G)) {
          start <- group_start[j] + 1L
          K     <- group_size[j]
          idx   <- start:(start + K - 1L)
          Xj <- X[idx, , drop = FALSE]
          oj <- offset[idx]
          eta_c <- as.numeric(Xj %*% beta_new) + oj
          max_eta_c <- max(eta_c)
          lse <- max_eta_c + log(sum(exp(eta_c - max_eta_c)))
          c_local <- which(choice[idx] == 1L)
          if (length(c_local) == 0L) next
          ll_cand <- ll_cand + eta_c[c_local[1L]] - lse
        }
        if (ll_cand >= loglik - 1e-10) break
        step_size <- step_size * 0.5
        beta_new  <- beta + step_size * delta
      }
      if (h > 1L && verbose) {
        cat(sprintf("  Step-halving: %d halvings, step_size = %.4e\n",
                    h - 1L, step_size))
      }
    }

    beta <- beta_new
  }

  # --- Recompute Hessian at final beta for variance ---
  hess[] <- 0.0
  for (j in seq_len(G)) {
    start <- group_start[j] + 1L
    K     <- group_size[j]
    idx   <- start:(start + K - 1L)
    Xj <- X[idx, , drop = FALSE]
    oj <- offset[idx]
    eta <- as.numeric(Xj %*% beta) + oj
    max_eta <- max(eta)
    exp_eta <- exp(eta - max_eta)
    prob    <- exp_eta / sum(exp_eta)
    xbar    <- as.numeric(crossprod(Xj, prob))
    Xj_w   <- Xj * sqrt(prob)
    hess    <- hess - (crossprod(Xj_w) - tcrossprod(xbar))
  }

  vcov <- solve(-hess)

  list(
    coefficients = beta,
    vcov         = vcov,
    hessian      = hess,
    loglik       = loglik,
    iterations   = iter,
    converged    = converged,
    gradient     = grad
  )
}


# =============================================================================
# 2. SANDWICH SEs (replaces clogit_sandwich.cpp)
# =============================================================================

clogit_sandwich_r <- function(X, choice, offset, group_start, group_size,
                              cluster_id, beta, hess_inv) {
  p <- ncol(X)
  G <- length(group_start)
  n_clusters <- max(cluster_id) + 1L

  # Accumulate per-cluster score: U[, e] = sum of per-group scores for cluster e
  U <- matrix(0.0, nrow = p, ncol = n_clusters)

  for (j in seq_len(G)) {
    start <- group_start[j] + 1L
    K     <- group_size[j]
    idx   <- start:(start + K - 1L)

    Xj <- X[idx, , drop = FALSE]
    oj <- offset[idx]

    eta <- as.numeric(Xj %*% beta) + oj
    max_eta <- max(eta)
    exp_eta <- exp(eta - max_eta)
    prob    <- exp_eta / sum(exp_eta)
    xbar    <- as.numeric(crossprod(Xj, prob))

    c_local <- which(choice[idx] == 1L)
    if (length(c_local) == 0L) next

    score_j <- Xj[c_local[1L], ] - xbar

    cl <- cluster_id[j] + 1L   # 0-based → 1-based column
    U[, cl] <- U[, cl] + score_j
  }

  # Meat: B = U %*% t(U)
  B <- tcrossprod(U)

  # Small-sample correction: C/(C-1) * (G-1)/G  (matches survival::coxph)
  correction <- (n_clusters / (n_clusters - 1)) * ((G - 1) / G)

  # Sandwich: correction * H^{-1} B H^{-1}
  vcov_robust <- correction * (hess_inv %*% B %*% hess_inv)

  list(
    vcov_robust  = vcov_robust,
    bread        = hess_inv,
    meat         = B,
    n_clusters   = n_clusters,
    df_correction = correction
  )
}


# =============================================================================
# 3. MAIN WRAPPER (same interface as Rcpp version)
# =============================================================================

fastclogit <- function(X, choice, strata, offset = NULL, cluster = NULL,
                       max_iter = 25L, tol = 1e-6, verbose = FALSE) {

  cl <- match.call()

  # --- Input validation ---
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or data.frame coercible to numeric matrix")
  }

  n <- nrow(X)
  p <- ncol(X)

  if (length(choice) != n) stop("length(choice) must equal nrow(X)")
  if (length(strata) != n) stop("length(strata) must equal nrow(X)")

  choice <- as.integer(choice)
  if (!all(choice %in% c(0L, 1L))) stop("choice must be 0/1 or logical")

  if (is.null(offset)) {
    offset <- rep(0.0, n)
  } else {
    if (length(offset) != n) stop("length(offset) must equal nrow(X)")
    offset <- as.numeric(offset)
  }

  # --- Sort by strata ---
  ord <- order(strata)
  X       <- X[ord, , drop = FALSE]
  choice  <- choice[ord]
  strata_sorted <- strata[ord]
  offset  <- offset[ord]
  if (!is.null(cluster)) {
    cluster_orig <- cluster[ord]
  }

  # --- Compute group boundaries (0-based for consistency with C++ version) ---
  grp_rle     <- rle(as.character(strata_sorted))
  group_size  <- as.integer(grp_rle$lengths)
  group_start <- as.integer(c(0L, cumsum(group_size[-length(group_size)])))
  n_groups    <- length(group_size)

  # --- Validate: each group has exactly one chosen ---
  total_chosen <- sum(choice)
  if (total_chosen != n_groups) {
    chosen_per_group <- tapply(choice, rep(seq_along(group_size), group_size), sum)
    n_zero <- sum(chosen_per_group == 0)
    n_multi <- sum(chosen_per_group > 1)
    if (n_zero > 0) warning(n_zero, " group(s) have no chosen alternative")
    if (n_multi > 0) warning(n_multi, " group(s) have multiple chosen alternatives")
  }

  # --- Drop zero-variance columns ---
  col_vars <- apply(X, 2, stats::var)
  zero_var <- which(col_vars < .Machine$double.eps)
  if (length(zero_var) > 0) {
    warning("Dropping ", length(zero_var), " zero-variance column(s): ",
            paste(colnames(X)[zero_var], collapse = ", "))
    X <- X[, -zero_var, drop = FALSE]
    p <- ncol(X)
  }

  # --- Fit via pure R Newton-Raphson ---
  fit <- clogit_fit_r(X, choice, offset, group_start, group_size,
                      as.integer(max_iter), tol, verbose)

  # --- Attach column names ---
  cnames <- colnames(X)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(p))
  names(fit$coefficients) <- cnames
  rownames(fit$vcov) <- colnames(fit$vcov) <- cnames
  names(fit$gradient) <- cnames
  fit$se <- setNames(sqrt(diag(fit$vcov)), cnames)

  # --- Clustered sandwich SEs ---
  if (!is.null(cluster)) {
    cluster_fac <- as.integer(as.factor(cluster_orig)) - 1L
    group_cluster <- cluster_fac[group_start + 1L]

    sandwich <- clogit_sandwich_r(
      X, choice, offset, group_start, group_size,
      group_cluster, fit$coefficients, fit$vcov
    )

    fit$vcov_robust <- sandwich$vcov_robust
    rownames(fit$vcov_robust) <- colnames(fit$vcov_robust) <- cnames
    fit$se_robust <- setNames(sqrt(diag(fit$vcov_robust)), cnames)
    fit$n_clusters <- sandwich$n_clusters
  }

  # --- Metadata ---
  fit$n_obs    <- n
  fit$n_groups <- n_groups
  fit$call     <- cl
  fit$terms    <- cnames

  if (!fit$converged) {
    warning("fastclogit did not converge in ", fit$iterations, " iterations. ",
            "Consider increasing max_iter or checking for separation.")
  }

  class(fit) <- "fastclogit"
  fit
}


# =============================================================================
# 4. S3 METHODS: summary, print, vcov, confint, tidy
# =============================================================================

coef.fastclogit <- function(object, ...) object$coefficients

vcov.fastclogit <- function(object, robust = TRUE, ...) {
  if (robust && !is.null(object$vcov_robust)) {
    return(object$vcov_robust)
  }
  object$vcov
}

logLik.fastclogit <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(object$coefficients)
  attr(ll, "nobs") <- object$n_obs
  class(ll) <- "logLik"
  ll
}

print.fastclogit <- function(x, ...) {
  cat("fastclogit (pure R)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 6))
  cat("\nLog-likelihood:", round(x$loglik, 4), "\n")
  cat("Observations:", format(x$n_obs, big.mark = ","),
      " | Groups:", format(x$n_groups, big.mark = ","))
  if (!is.null(x$n_clusters)) {
    cat(" | Clusters:", format(x$n_clusters, big.mark = ","))
  }
  cat("\nConverged:", x$converged, "(", x$iterations, "iterations )\n")
  invisible(x)
}

summary.fastclogit <- function(object, robust = TRUE, ...) {
  coefs <- object$coefficients
  p <- length(coefs)

  if (robust && !is.null(object$se_robust)) {
    se <- object$se_robust
    se_type <- "robust (clustered)"
  } else {
    se <- object$se
    se_type <- "model-based"
  }

  z_val <- coefs / se
  p_val <- 2 * stats::pnorm(-abs(z_val))

  # Note: cbind(Name = named_vec) can silently drop column names in R 4.5.x.
  coef_table <- cbind(coefs, se, z_val, p_val)
  colnames(coef_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  structure(list(
    call       = object$call,
    coefficients = coef_table,
    loglik     = object$loglik,
    n_obs      = object$n_obs,
    n_groups   = object$n_groups,
    n_clusters = object$n_clusters,
    converged  = object$converged,
    iterations = object$iterations,
    se_type    = se_type
  ), class = "summary.fastclogit")
}

print.summary.fastclogit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("fastclogit (pure R)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Coefficients (SEs:", x$se_type, "):\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE,
               na.print = "NA", cs.ind = 1:2, tst.ind = 3, P.values = TRUE,
               has.Pvalue = TRUE)

  cat("\nLog-likelihood:", round(x$loglik, 4), "\n")
  cat("Observations:", format(x$n_obs, big.mark = ","),
      " | Groups:", format(x$n_groups, big.mark = ","))
  if (!is.null(x$n_clusters)) {
    cat(" | Clusters:", format(x$n_clusters, big.mark = ","))
  }
  cat("\nConverged:", x$converged, "(", x$iterations, "iterations )\n")
  invisible(x)
}

confint.fastclogit <- function(object, parm, level = 0.95, robust = TRUE, ...) {
  cf  <- coef(object)
  ses <- if (robust && !is.null(object$se_robust)) object$se_robust else object$se
  a   <- (1 - level) / 2
  z   <- stats::qnorm(1 - a)
  ci  <- cbind(cf - z * ses, cf + z * ses)
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

tidy_fastclogit <- function(object, conf.int = TRUE, conf.level = 0.95,
                            robust = TRUE, exponentiate = FALSE) {
  coefs <- object$coefficients
  if (robust && !is.null(object$se_robust)) {
    se <- object$se_robust
  } else {
    se <- object$se
  }

  # Safety: ensure names are present (they can be lost during unscaling)
  if (is.null(names(coefs))) {
    if (!is.null(object$terms)) {
      names(coefs) <- object$terms
    } else {
      names(coefs) <- paste0("V", seq_along(coefs))
    }
  }
  if (is.null(names(se))) names(se) <- names(coefs)

  z_val <- coefs / se
  p_val <- 2 * stats::pnorm(-abs(z_val))

  out <- data.frame(
    term      = names(coefs),
    estimate  = unname(coefs),
    std.error = unname(se),
    statistic = unname(z_val),
    p.value   = unname(p_val),
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    ci <- confint(object, level = conf.level, robust = robust)
    out$conf.low  <- ci[, 1]
    out$conf.high <- ci[, 2]
  }

  if (exponentiate) {
    out$estimate  <- exp(out$estimate)
    if (conf.int) {
      out$conf.low  <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }

  out
}


# =============================================================================
# 5. SIMULATE DATA (same DGP as Rcpp version)
# =============================================================================

simulate_clogit_data <- function(
  n_egos        = 5000L,
  n_alts        = 100L,
  beta_continuous = c(lnDist = -0.5,
                      n_years_same_cfar = 1.0,
                      n_years_same_peorg = 0.8,
                      n_years_same_uni = 0.6),
  beta_factor   = list(
    Edudiff3   = c(Edudiff3_BothHigh = 0.8,
                   Edudiff3_BothLow = 0.4,
                   Edudiff3_Missing = -0.3),
    AgeDiffcat = c(AgeDiffcat_2 = -0.2,
                   AgeDiffcat_3 = -0.5,
                   AgeDiffcat_4 = -0.8,
                   AgeDiffcat_5 = -1.2)
  ),
  use_offset    = TRUE,
  strata_props  = c(0.05, 0.10, 0.10, 0.60, 0.15),
  strata_popsizes = c(20, 100, 50, 5000, 50000),
  cluster_ratio = 1.3,
  seed          = 42L
) {
  set.seed(seed)

  n_total <- as.integer(n_egos) * as.integer(n_alts)
  n_continuous <- length(beta_continuous)

  beta_factor_vec <- unlist(beta_factor, use.names = FALSE)
  names(beta_factor_vec) <- unlist(lapply(beta_factor, names))
  beta_true <- c(beta_continuous, beta_factor_vec)
  p <- length(beta_true)

  strata_id <- rep(seq_len(n_egos), each = n_alts)

  X_cont <- matrix(0, nrow = n_total, ncol = n_continuous)
  colnames(X_cont) <- names(beta_continuous)

  for (i in seq_len(n_continuous)) {
    nm <- names(beta_continuous)[i]
    if (grepl("lnDist", nm)) {
      X_cont[, i] <- stats::rnorm(n_total, mean = 9, sd = 2)
    } else if (grepl("n_years", nm)) {
      X_cont[, i] <- stats::rpois(n_total, lambda = 0.3)
    } else {
      X_cont[, i] <- stats::rnorm(n_total, mean = 0, sd = 1)
    }
  }

  X_factor_list <- list()
  factor_raw <- list()

  for (fac_name in names(beta_factor)) {
    betas <- beta_factor[[fac_name]]
    n_levels <- length(betas) + 1
    fac_vals <- sample(seq_len(n_levels), n_total, replace = TRUE)
    factor_raw[[fac_name]] <- fac_vals

    dummy_mat <- matrix(0, nrow = n_total, ncol = length(betas))
    colnames(dummy_mat) <- names(betas)
    for (k in seq_along(betas)) {
      dummy_mat[, k] <- as.integer(fac_vals == (k + 1))
    }
    X_factor_list[[fac_name]] <- dummy_mat
  }

  X_factor <- do.call(cbind, X_factor_list)
  X <- cbind(X_cont, X_factor)

  offset_vec <- rep(0, n_total)
  sampling_stratum <- rep(NA_integer_, n_total)

  if (use_offset) {
    for (ego in seq_len(n_egos)) {
      rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
      stratum_assignment <- sample(
        seq_along(strata_props), size = n_alts,
        replace = TRUE, prob = strata_props
      )
      sampling_stratum[rows] <- stratum_assignment
      for (s in seq_along(strata_props)) {
        s_rows <- rows[stratum_assignment == s]
        if (length(s_rows) > 0) {
          n_s <- length(s_rows)
          N_s <- strata_popsizes[s]
          offset_vec[s_rows] <- -log(n_s / N_s)
        }
      }
    }
  }

  eta <- X %*% beta_true + offset_vec
  choice <- integer(n_total)

  for (ego in seq_len(n_egos)) {
    rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
    eta_j <- eta[rows]
    max_eta <- max(eta_j)
    prob <- exp(eta_j - max_eta)
    prob <- prob / sum(prob)
    chosen_idx <- sample.int(n_alts, size = 1, prob = prob)
    choice[rows[chosen_idx]] <- 1L
  }

  n_unique_clusters <- max(1L, round(n_egos / cluster_ratio))
  ego_cluster <- sample(seq_len(n_unique_clusters), n_egos, replace = TRUE)
  cluster_id <- rep(ego_cluster, each = n_alts)

  df <- data.frame(
    choice = choice, strata_id = strata_id, cluster_id = cluster_id,
    X, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (use_offset) df$correction <- offset_vec
  for (fac_name in names(factor_raw)) {
    df[[paste0(fac_name, "_raw")]] <- factor(factor_raw[[fac_name]])
  }

  list(
    X = X, choice = choice, strata = strata_id,
    offset = offset_vec, cluster = cluster_id,
    beta_true = beta_true, data = df,
    n_egos = n_egos, n_alts = n_alts
  )
}


cat("=== fastclogit (pure R) loaded ===\n")
cat("  No Rcpp/compiler needed. Same interface as Rcpp version.\n")
cat("  Usage: fit <- fastclogit(X, choice, strata, offset=..., cluster=...)\n\n")
