#' Simulate Conditional Logit Data (Partner Choice DGP)
#'
#' Generates synthetic data mimicking the LISA partner choice pipeline:
#' stratified importance sampling with McFadden/Manski offsets,
#' factor and continuous predictors, and clustered egos.
#'
#' @param n_egos Integer. Number of choice sets (egos).
#' @param n_alts Integer. Alternatives per choice set (including chosen).
#' @param beta_continuous Numeric vector. True coefficients for continuous
#'   predictors. Length determines number of continuous predictors.
#'   Default mimics: lnDist(-0.5), n_years_same_cfar(1.0),
#'   n_years_same_peorg(0.8), n_years_same_uni(0.6).
#' @param beta_factor List of numeric vectors. True coefficients for factor
#'   dummy variables (excluding reference level). Names become factor names.
#'   Default mimics Edudiff3 (4 levels, ref=1) and AgeDiffcat (5 levels, ref=1).
#' @param use_offset Logical. Simulate McFadden/Manski stratified sampling?
#' @param strata_props Numeric vector of length 5. Sampling proportions per
#'   stratum (Cfar/PeOrg/Uni/County/Rest). Default: c(5, 10, 10, 60, 15)/100.
#' @param strata_popsizes Numeric vector of length 5. Typical population sizes
#'   per stratum for offset computation.
#' @param cluster_ratio Numeric > 0. Ratio of egos to unique cluster IDs.
#'   1.0 = each ego is unique; 1.5 = some egos appear ~1.5 times on average.
#' @param seed Integer. Random seed.
#'
#' @return A list with components:
#'   \item{X}{Numeric matrix (n_total x p) — design matrix}
#'   \item{choice}{Integer vector (n_total) — 1 for chosen, 0 otherwise}
#'   \item{strata}{Integer vector (n_total) — choice set ID}
#'   \item{offset}{Numeric vector (n_total) — McFadden/Manski correction (0 if !use_offset)}
#'   \item{cluster}{Integer vector (n_total) — cluster ID}
#'   \item{beta_true}{Named numeric vector — true coefficients}
#'   \item{data}{data.frame with all columns for easy clogit() comparison}
#'   \item{n_egos}{Number of choice sets}
#'   \item{n_alts}{Alternatives per set}
#'
#' @examples
#' sim <- simulate_clogit_data(n_egos = 500, n_alts = 50)
#' str(sim$X)
#' table(sim$choice)
#'
#' @export
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

  # --- True beta vector ---
  # Use recursive=FALSE to avoid prepending list names (e.g., "Edudiff3.Edudiff3_BothHigh")
  beta_factor_vec <- unlist(beta_factor, use.names = FALSE)
  names(beta_factor_vec) <- unlist(lapply(beta_factor, names))
  beta_true <- c(beta_continuous, beta_factor_vec)
  p <- length(beta_true)

  # --- Strata IDs ---
  strata_id <- rep(seq_len(n_egos), each = n_alts)

  # --- Continuous predictors ---
  # lnDist: log-normal-ish (shifted, so variation is realistic)
  # Institution overlaps: mostly 0, some positive (mimics count data)
  X_cont <- matrix(0, nrow = n_total, ncol = n_continuous)
  colnames(X_cont) <- names(beta_continuous)

  for (i in seq_len(n_continuous)) {
    nm <- names(beta_continuous)[i]
    if (grepl("lnDist", nm)) {
      # Log distance: centered around 8-10 (few km to hundreds of km)
      X_cont[, i] <- rnorm(n_total, mean = 9, sd = 2)
    } else if (grepl("n_years", nm)) {
      # Institutional overlap: mostly 0, some 1-3+
      X_cont[, i] <- rpois(n_total, lambda = 0.3)
    } else {
      X_cont[, i] <- rnorm(n_total, mean = 0, sd = 1)
    }
  }

  # --- Factor predictors (as dummies, excluding reference) ---
  X_factor_list <- list()
  factor_raw <- list()  # keep original factor for data.frame version

  for (fac_name in names(beta_factor)) {
    betas <- beta_factor[[fac_name]]
    n_levels <- length(betas) + 1  # +1 for reference
    # Draw random factor levels (1 = reference, 2..n_levels = dummies)
    fac_vals <- sample(seq_len(n_levels), n_total, replace = TRUE)
    factor_raw[[fac_name]] <- fac_vals

    # Create dummy matrix
    dummy_mat <- matrix(0, nrow = n_total, ncol = length(betas))
    colnames(dummy_mat) <- names(betas)
    for (k in seq_along(betas)) {
      dummy_mat[, k] <- as.integer(fac_vals == (k + 1))
    }
    X_factor_list[[fac_name]] <- dummy_mat
  }

  X_factor <- do.call(cbind, X_factor_list)

  # --- Full design matrix ---
  X <- cbind(X_cont, X_factor)

  # --- Offset (McFadden/Manski correction) ---
  offset_vec <- rep(0, n_total)
  sampling_stratum <- rep(NA_integer_, n_total)

  if (use_offset) {
    for (ego in seq_len(n_egos)) {
      rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
      # Assign each alternative to a sampling stratum
      stratum_assignment <- sample(
        seq_along(strata_props),
        size = n_alts,
        replace = TRUE,
        prob = strata_props
      )
      sampling_stratum[rows] <- stratum_assignment

      # Compute correction per stratum
      for (s in seq_along(strata_props)) {
        s_rows <- rows[stratum_assignment == s]
        if (length(s_rows) > 0) {
          n_s <- length(s_rows)
          N_s <- strata_popsizes[s]
          # correction = -log(n_s / N_s)
          offset_vec[s_rows] <- -log(n_s / N_s)
        }
      }
    }
  }

  # --- Generate choices from the model ---
  eta <- X %*% beta_true + offset_vec
  choice <- integer(n_total)

  for (ego in seq_len(n_egos)) {
    rows <- ((ego - 1) * n_alts + 1):(ego * n_alts)
    eta_j <- eta[rows]
    # Softmax with LogSumExp trick
    max_eta <- max(eta_j)
    prob <- exp(eta_j - max_eta)
    prob <- prob / sum(prob)
    # Draw chosen alternative
    chosen_idx <- sample.int(n_alts, size = 1, prob = prob)
    choice[rows[chosen_idx]] <- 1L
  }

  # --- Cluster IDs ---
  # cluster_ratio > 1 means some egos share a cluster (repeated observations)
  n_unique_clusters <- max(1L, round(n_egos / cluster_ratio))
  ego_cluster <- sample(seq_len(n_unique_clusters), n_egos, replace = TRUE)
  cluster_id <- rep(ego_cluster, each = n_alts)

  # --- Build data.frame for clogit() comparison ---
  df <- data.frame(
    choice = choice,
    strata_id = strata_id,
    cluster_id = cluster_id,
    X,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (use_offset) {
    df$correction <- offset_vec
  }

  # Add raw factor columns for clogit formula with factor()
  for (fac_name in names(factor_raw)) {
    df[[paste0(fac_name, "_raw")]] <- factor(factor_raw[[fac_name]])
  }

  list(
    X         = X,
    choice    = choice,
    strata    = strata_id,
    offset    = offset_vec,
    cluster   = cluster_id,
    beta_true = beta_true,
    data      = df,
    n_egos    = n_egos,
    n_alts    = n_alts
  )
}
