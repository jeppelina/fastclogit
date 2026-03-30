#' KHB Decomposition for Conditional Logit Models
#'
#' Implements the Kohler, Karlson & Holm (2011) method to decompose total
#' effects into direct and indirect (mediated) effects in conditional logit
#' models, while correctly accounting for rescaling bias.
#'
#' In nonlinear models like logit, adding mediator variables changes the scale
#' of all coefficients (because the residual variance changes). This means you
#' cannot simply compare coefficients across nested models to assess mediation.
#' The KHB method solves this by residualizing the mediators on the key
#' variables within strata, so that both models share the same residual
#' variance scale.
#'
#' @section Method:
#' \enumerate{
#'   \item Residualize each Z-variable on the key variables X and controls C
#'     within strata: \code{Z_resid = residuals(lm(Z ~ X + C | strata))}
#'   \item Fit a \strong{reduced} model: \code{Y ~ X + Z_resid + C} (total
#'     effect of X, on the same scale as the full model)
#'   \item Fit a \strong{full} model: \code{Y ~ X + Z + C} (direct effect of X)
#'   \item Indirect effect = Total - Direct
#'   \item Confounding ratio = Total / Direct
#' }
#'
#' @param data A data.frame or data.table.
#' @param key_vars Character vector of key variable names (X) whose effects to
#'   decompose. Supports factor variables and interaction terms (e.g.,
#'   \code{"edu:age"}).
#' @param z_vars Character vector of mediator variable names (Z). Must be
#'   simple column names (no interactions).
#' @param controls Character vector of control variable names (C), or
#'   \code{NULL}.
#' @param strata Character string naming the choice-set/strata column.
#' @param cluster Optional character string naming the cluster column for
#'   robust SEs. If \code{NULL}, model-based SEs are used.
#' @param choice Character string naming the binary choice/outcome column
#'   (default \code{"choice"}).
#' @param offset Optional character string naming an offset column (e.g.,
#'   McFadden/Manski correction). If \code{NULL}, no offset.
#' @param verbose Logical. Print progress? Default \code{TRUE}.
#'
#' @return A list with:
#'   \item{success}{Logical.}
#'   \item{decomposition}{A data.frame with columns: \code{variable},
#'     \code{coefficient}, \code{total_effect}, \code{total_se},
#'     \code{direct_effect}, \code{direct_se}, \code{indirect_effect},
#'     \code{indirect_se}, \code{conf_ratio}, \code{conf_pct}.}
#'   \item{z_effects}{A data.frame of mediator effects from the full model.}
#'   \item{coef_reduced}{Coefficient matrix from the reduced model.}
#'   \item{coef_full}{Coefficient matrix from the full model.}
#'   \item{n_obs}{Number of observations used.}
#'   \item{n_groups}{Number of choice sets.}
#'
#' @references
#' Kohler, U., Karlson, K. B. & Holm, A. (2011). Comparing coefficients of
#' nested nonlinear probability models. \emph{The Stata Journal}, 11(3),
#' 420--438.
#'
#' @examples
#' \dontrun{
#' # Generate data with a known mediation structure
#' set.seed(42)
#' n_egos <- 500; n_alts <- 30
#' n <- n_egos * n_alts
#' d <- data.frame(
#'   strata_id = rep(1:n_egos, each = n_alts),
#'   x = rnorm(n),
#'   c1 = rnorm(n)
#' )
#' d$z <- 0.6 * d$x + rnorm(n)  # z is mediated by x
#' eta <- 0.5 * d$x + 0.8 * d$z + 0.3 * d$c1
#'
#' # Generate choices via softmax within strata
#' d$choice <- 0L
#' for (i in 1:n_egos) {
#'   rows <- ((i-1)*n_alts + 1):(i*n_alts)
#'   probs <- exp(eta[rows] - max(eta[rows]))
#'   probs <- probs / sum(probs)
#'   d$choice[rows[sample(n_alts, 1, prob = probs)]] <- 1L
#' }
#'
#' result <- khb_decompose(d, key_vars = "x", z_vars = "z",
#'                          controls = "c1", strata = "strata_id",
#'                          choice = "choice")
#' result$decomposition
#' }
#'
#' @export
khb_decompose <- function(data,
                           key_vars,
                           z_vars,
                           controls = NULL,
                           strata,
                           cluster = NULL,
                           choice = "choice",
                           offset = NULL,
                           verbose = TRUE) {

  .msg <- function(...) if (verbose) message(...)

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  # MEMORY: We never copy the full dataset. Residuals are stored in a
  # separate list, and we build slim data frames (only needed columns)
  # for the fclogit() calls.

  has_interaction <- grepl(":", z_vars)
  if (any(has_interaction)) {
    stop("Interaction terms in z_vars are not supported. Found: ",
         paste(z_vars[has_interaction], collapse = ", "))
  }

  # Column existence checks
  .terms_to_cols <- function(terms) unique(unlist(strsplit(terms, ":")))
  all_term_vars <- unique(c(key_vars, z_vars, controls))
  all_raw_cols <- .terms_to_cols(all_term_vars)
  needed <- unique(c(choice, strata, all_raw_cols))
  if (!is.null(cluster)) needed <- c(needed, cluster)
  if (!is.null(offset))  needed <- c(needed, offset)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    stop("Column(s) not found in data: ", paste(missing, collapse = ", "))
  }

  # Validate choice sets: each stratum must have exactly 1 chosen
  chosen_per <- tapply(data[[choice]], data[[strata]], sum)
  bad <- sum(chosen_per != 1)
  if (bad > 0) stop(bad, " strata do not have exactly 1 chosen alternative")

  .msg("========== KHB Decomposition (fastclogit) ==========\n")
  .msg("Key variables (X):  ", paste(key_vars, collapse = ", "))
  .msg("Mediators (Z):      ", paste(z_vars, collapse = ", "))
  if (!is.null(controls)) .msg("Controls (C):       ", paste(controls, collapse = ", "))
  .msg("Strata:             ", strata)
  .msg("")

  # -------------------------------------------------------------------------
  # Step 1: Residualize Z on X + C within strata
  # -------------------------------------------------------------------------
  .msg("Step 1: Residualizing Z on X + C within strata")

  resid_predictors <- c(key_vars, controls)
  resid_columns <- .terms_to_cols(resid_predictors)

  # Convert characters to factors in a slim copy of needed columns only
  # (avoids triggering copy-on-modify on the full dataset)
  slim_cols <- unique(c(all_raw_cols, choice, strata))
  if (!is.null(cluster)) slim_cols <- c(slim_cols, cluster)
  if (!is.null(offset))  slim_cols <- c(slim_cols, offset)
  slim <- data[, slim_cols, drop = FALSE]
  if (inherits(slim, "data.table")) slim <- as.data.frame(slim)

  for (v in all_raw_cols) {
    if (is.character(slim[[v]])) slim[[v]] <- factor(slim[[v]])
  }

  # Build predictor design matrix for residualization
  mm_formula <- stats::as.formula(paste("~ 0 +", paste(resid_predictors, collapse = " + ")))
  X_dm <- stats::model.matrix(mm_formula, data = slim)

  # Within-strata demeaning
  strata_ids <- slim[[strata]]
  uid <- unique(strata_ids)
  strata_int <- match(strata_ids, uid)
  n_strata <- length(uid)
  grp_sums <- rowsum(X_dm, strata_int, reorder = FALSE, na.rm = TRUE)
  grp_n <- tabulate(strata_int, nbins = n_strata)
  grp_means <- grp_sums / grp_n
  for (j in seq_len(ncol(X_dm))) {
    X_dm[, j] <- X_dm[, j] - grp_means[strata_int, j]
  }
  rm(grp_sums, grp_means)

  # Completeness mask for residualization predictors
  n_rows <- nrow(slim)
  X_ok <- rep(TRUE, n_rows)
  for (v in resid_columns) X_ok <- X_ok & !is.na(slim[[v]])

  # Expand factor Z to dummies and residualize
  # MEMORY: residuals stored in a list, NOT added to data/slim
  z_expanded <- list()
  z_orig_map <- list()
  resid_store <- list()  # name -> numeric vector of residuals

  for (z in z_vars) {
    if (is.factor(slim[[z]])) {
      mm <- stats::model.matrix(~ 0 + slim[[z]])
      dummy_names <- make.names(paste0(z, "_", levels(slim[[z]])), unique = TRUE)
      colnames(mm) <- dummy_names
      # Store dummy columns in resid_store temporarily (will be overwritten with residuals)
      for (dn in dummy_names) resid_store[[dn]] <- mm[, dn]
      z_expanded[[z]] <- dummy_names
      for (dn in dummy_names) z_orig_map[[dn]] <- z
    } else {
      z_expanded[[z]] <- z
      z_orig_map[[z]] <- z
    }
  }

  z_all <- unlist(z_expanded, use.names = FALSE)
  z_resid_names <- paste0(z_all, "_resid")

  for (i in seq_along(z_all)) {
    z_col <- z_all[i]
    z_orig <- z_orig_map[[z_col]]
    z_resid <- z_resid_names[i]

    # Get raw values: from resid_store (factor dummies) or slim (numeric)
    z_vec <- if (!is.null(resid_store[[z_col]])) {
      as.numeric(resid_store[[z_col]])
    } else {
      as.numeric(slim[[z_col]])
    }
    z_dm <- z_vec - ave(z_vec, strata_ids, FUN = function(x) mean(x, na.rm = TRUE))
    ok <- X_ok & !is.na(slim[[z_orig]])

    if (sum(ok) == 0) {
      warning("No complete cases for residualization of '", z_col, "'")
      resid_store[[z_resid]] <- rep(NA_real_, n_rows)
      next
    }

    fit <- stats::lm.fit(x = X_dm[ok, , drop = FALSE], y = z_dm[ok])
    coefs <- fit$coefficients
    coefs[is.na(coefs)] <- 0

    resid_vec <- rep(NA_real_, n_rows)
    resid_vec[ok] <- as.numeric(z_dm[ok] - X_dm[ok, , drop = FALSE] %*% coefs)
    resid_store[[z_resid]] <- resid_vec

    if (verbose) {
      ss_res <- sum(resid_vec[ok]^2, na.rm = TRUE)
      ss_tot <- sum(z_dm[ok]^2, na.rm = TRUE)
      r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else 0
      .msg(sprintf("  %s -> %s (demeaned R^2 = %.4f)", z_col, z_resid, r2))
    }
  }

  rm(X_dm, X_ok)
  gc(verbose = FALSE)

  # Check for all-NA residual columns
  all_na <- vapply(z_resid_names, function(zr) all(is.na(resid_store[[zr]])), logical(1))
  if (any(all_na)) {
    stop("Residual column(s) entirely NA: ", paste(z_resid_names[all_na], collapse = ", "))
  }

  # -------------------------------------------------------------------------
  # Step 2: Build slim data frame for model fitting
  # -------------------------------------------------------------------------
  # MEMORY: Instead of adding residual columns to the (possibly huge) input
  # data, we build a slim data frame containing only the columns needed by
  # the two fclogit() calls.  This keeps peak memory to ~2x the columns
  # actually used, rather than copying the entire dataset.

  # Columns shared by both models
  meta_cols <- c(choice, strata)
  if (!is.null(cluster)) meta_cols <- c(meta_cols, cluster)
  if (!is.null(offset))  meta_cols <- c(meta_cols, offset)
  shared_pred_cols <- unique(c(.terms_to_cols(key_vars), .terms_to_cols(controls)))

  # Build reduced-model data: slim[meta + X + C] + residual columns from resid_store
  fit_data <- slim[, unique(c(meta_cols, shared_pred_cols, z_vars)), drop = FALSE]
  for (zr in z_resid_names) {
    fit_data[[zr]] <- resid_store[[zr]]
  }

  # Free resid_store and slim now — no longer needed
  rm(resid_store, slim)
  gc(verbose = FALSE)

  # -------------------------------------------------------------------------
  # Step 2b: Fit reduced model (X + Z_resid + C) -> total effect
  # -------------------------------------------------------------------------
  .msg("\nStep 2: Fitting reduced model (total effect)")

  reduced_terms <- c(key_vars, z_resid_names, controls)
  reduced_formula <- stats::as.formula(
    paste(choice, "~", paste(reduced_terms, collapse = " + "))
  )

  fit_reduced <- tryCatch(
    fclogit(reduced_formula, data = fit_data, strata = strata,
            cluster = cluster, offset = offset,
            max_iter = 100L, verbose = FALSE),
    error = function(e) {
      .msg("  ERROR: ", e$message)
      NULL
    }
  )

  if (is.null(fit_reduced)) {
    return(list(success = FALSE, error = "Reduced model failed to fit"))
  }

  .msg("  Converged: ", fit_reduced$converged,
       " (", fit_reduced$iterations, " iterations)")

  coef_reduced <- .build_khb_coef_matrix(fit_reduced)

  # -------------------------------------------------------------------------
  # Step 3: Fit full model (X + Z + C) -> direct effect
  # -------------------------------------------------------------------------
  .msg("\nStep 3: Fitting full model (direct effect)")

  full_terms <- c(key_vars, z_vars, controls)
  full_formula <- stats::as.formula(
    paste(choice, "~", paste(full_terms, collapse = " + "))
  )

  fit_full <- tryCatch(
    fclogit(full_formula, data = fit_data, strata = strata,
            cluster = cluster, offset = offset,
            max_iter = 100L, verbose = FALSE),
    error = function(e) {
      .msg("  ERROR: ", e$message)
      NULL
    }
  )

  # Free fit_data — both models are fitted
  rm(fit_data)
  gc(verbose = FALSE)

  if (is.null(fit_full)) {
    return(list(success = FALSE, error = "Full model failed to fit"))
  }

  .msg("  Converged: ", fit_full$converged,
       " (", fit_full$iterations, " iterations)")

  coef_full <- .build_khb_coef_matrix(fit_full)

  # -------------------------------------------------------------------------
  # Step 4: Decompose
  # -------------------------------------------------------------------------
  .msg("\nStep 4: Computing decomposition")

  results_list <- list()

  for (x in key_vars) {
    cn_reduced <- .match_khb_names(x, rownames(coef_reduced))
    cn_full    <- .match_khb_names(x, rownames(coef_full))
    cn <- intersect(cn_reduced, cn_full)

    if (length(cn) == 0) {
      .msg("  Warning: '", x, "' not found in model coefficients")
      next
    }

    for (nm in cn) {
      b_total  <- coef_reduced[nm, "coef"]
      se_total <- coef_reduced[nm, "se"]
      b_direct  <- coef_full[nm, "coef"]
      se_direct <- coef_full[nm, "se"]
      b_indirect <- b_total - b_direct
      se_indirect <- sqrt(se_total^2 + se_direct^2)  # conservative

      conf_ratio <- if (abs(b_direct) > 1e-10) b_total / b_direct else NA_real_
      conf_pct <- if (abs(b_total) > 1e-10) 100 * (b_total - b_direct) / b_total else NA_real_

      results_list[[nm]] <- data.frame(
        variable = x,
        coefficient = nm,
        total_effect = b_total,
        total_se = se_total,
        direct_effect = b_direct,
        direct_se = se_direct,
        indirect_effect = b_indirect,
        indirect_se = se_indirect,
        conf_ratio = conf_ratio,
        conf_pct = conf_pct,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) == 0) {
    return(list(success = FALSE, error = "No coefficients matched key variables"))
  }

  decomposition <- do.call(rbind, results_list)
  rownames(decomposition) <- NULL

  # Z effects from full model
  z_effects_list <- list()
  for (z in z_vars) {
    z_names <- .match_khb_names(z, rownames(coef_full))
    for (zn in z_names) {
      z_effects_list[[zn]] <- data.frame(
        variable = z, coefficient = zn,
        estimate = coef_full[zn, "coef"],
        se = coef_full[zn, "se"],
        z_stat = coef_full[zn, "z"],
        p_value = coef_full[zn, "p"],
        stringsAsFactors = FALSE
      )
    }
  }
  z_effects <- if (length(z_effects_list) > 0) {
    do.call(rbind, z_effects_list)
  } else {
    data.frame()
  }
  rownames(z_effects) <- NULL

  # -------------------------------------------------------------------------
  # Print summary
  # -------------------------------------------------------------------------
  if (verbose) {
    .msg("\n========== KHB Results ==========\n")
    for (i in seq_len(nrow(decomposition))) {
      r <- decomposition[i, ]
      .msg(sprintf("  %s:", r$coefficient))
      .msg(sprintf("    Total  (reduced): %8.4f (SE: %.4f)", r$total_effect, r$total_se))
      .msg(sprintf("    Direct (full):    %8.4f (SE: %.4f)", r$direct_effect, r$direct_se))
      .msg(sprintf("    Indirect (diff):  %8.4f (SE: %.4f*)", r$indirect_effect, r$indirect_se))
      .msg(sprintf("    Confounding ratio: %.4f", r$conf_ratio))
      .msg(sprintf("    Confounding %%:     %.2f%%", r$conf_pct))
    }
    .msg("\n  * Indirect SE assumes independence (conservative)")
  }

  list(
    success = TRUE,
    decomposition = decomposition,
    z_effects = z_effects,
    coef_reduced = coef_reduced,
    coef_full = coef_full,
    n_obs = fit_full$n_obs,
    n_groups = fit_full$n_groups
  )
}


# ---------------------------------------------------------------------------
# Internal helpers for KHB
# ---------------------------------------------------------------------------

#' @keywords internal
.build_khb_coef_matrix <- function(fit) {
  cf <- coef(fit)
  se <- if (!is.null(fit$se_robust)) fit$se_robust else fit$se
  z_val <- cf / se
  p_val <- 2 * stats::pnorm(-abs(z_val))
  mat <- cbind(cf, se, z_val, p_val)
  colnames(mat) <- c("coef", "se", "z", "p")
  rownames(mat) <- names(cf)
  mat
}

#' @keywords internal
.match_khb_names <- function(var_name, coef_names) {
  # Exact match
  if (var_name %in% coef_names) return(var_name)
  # Interaction match
  if (grepl(":", var_name)) {
    parts <- strsplit(var_name, ":")[[1]]
    matches <- vapply(coef_names, function(cn) {
      cn_parts <- strsplit(cn, ":")[[1]]
      if (length(cn_parts) != length(parts)) return(FALSE)
      all(mapply(function(p, cp) startsWith(cp, p), parts, cn_parts))
    }, logical(1))
    return(coef_names[matches])
  }
  # Prefix match (factor dummies)
  matches <- coef_names[startsWith(coef_names, var_name) & !grepl(":", coef_names)]
  matches
}
