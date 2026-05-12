# =============================================================================
# gen_data.R -- Simulation generators for sparse-X validation harness.
#
# Each generator returns:
#   list(
#     X         = dense numeric matrix (n_rows x p),
#     Xsp       = Matrix::dgCMatrix, identical to X,
#     choice    = integer 0/1 vector (n_rows),
#     strata    = integer vector (n_rows), choice-set ID,
#     cluster   = integer vector (n_rows) or NULL,
#     beta_true = numeric vector (p),
#     meta      = list of bookkeeping info (density, formula etc.)
#   )
#
# Standard conditional-logit DGP:
#   eta = X beta_true + Gumbel(0,1) noise
#   per stratum, alternative with max eta gets choice == 1.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

# ----- helpers ---------------------------------------------------------------

.gumbel <- function(n) -log(-log(stats::runif(n)))

.pick_choice <- function(eta, strata) {
  # For each stratum, set choice=1 at the row with the max eta.
  n <- length(eta)
  choice <- integer(n)
  # split rows by strata, pick argmax within each
  ord <- order(strata)
  s_ord <- strata[ord]
  e_ord <- eta[ord]
  rle_s <- rle(s_ord)
  starts <- c(0L, cumsum(rle_s$lengths[-length(rle_s$lengths)]))
  for (g in seq_along(rle_s$lengths)) {
    rng <- (starts[g] + 1L):(starts[g] + rle_s$lengths[g])
    winner_local <- which.max(e_ord[rng])
    choice[ord[rng[winner_local]]] <- 1L
  }
  choice
}

.verify_sparse <- function(X, Xsp) {
  stopifnot(inherits(Xsp, "dgCMatrix"))
  stopifnot(nrow(X) == nrow(Xsp), ncol(X) == ncol(Xsp))
  same <- all(as.matrix(Xsp) == X)
  if (!isTRUE(same)) stop("Xsp does not round-trip to X (mismatch in dense/sparse).")
  invisible(TRUE)
}

# Build dgCMatrix efficiently from a dense matrix, preserving colnames.
.densify_to_sparse <- function(X) {
  nz <- which(X != 0, arr.ind = TRUE)
  Xsp <- Matrix::sparseMatrix(
    i = nz[, 1], j = nz[, 2], x = X[nz],
    dims = dim(X), dimnames = dimnames(X)
  )
  methods::as(Xsp, "CsparseMatrix")
}

# ----- 1. gen_small_dense ----------------------------------------------------
#
# 500 strata x 5 alters; 5 continuous covariates. Density ~100%.
# Sanity test that sparse path does not break on dense data.
# ----------------------------------------------------------------------------

gen_small_dense <- function(seed = 1L) {
  set.seed(seed)
  n_strata <- 500L
  n_alt    <- 5L
  p        <- 5L
  n        <- n_strata * n_alt

  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  beta_true <- stats::rnorm(p, sd = 0.5)

  strata <- rep(seq_len(n_strata), each = n_alt)
  eta    <- as.numeric(X %*% beta_true) + .gumbel(n)
  choice <- .pick_choice(eta, strata)

  Xsp <- .densify_to_sparse(X)
  .verify_sparse(X, Xsp)

  list(
    X = X, Xsp = Xsp,
    choice = choice, strata = strata, cluster = NULL,
    beta_true = beta_true,
    meta = list(
      name = "small_dense",
      formula = "~ x1+x2+x3+x4+x5",
      n_strata = n_strata, n_alt = n_alt, p = p, n_rows = n,
      density = mean(X != 0)
    )
  )
}

# ----- 2. gen_medium_factor --------------------------------------------------
#
# 5000 strata x 8 alters. Formula:  ~ x1*decade + factor_a*decade + factor_b
#   x1       continuous
#   decade   3 levels
#   factor_a 6 levels
#   factor_b 4 levels
# Density ~15%. 50 clusters.
# ----------------------------------------------------------------------------

gen_medium_factor <- function(seed = 2L) {
  set.seed(seed)
  n_strata <- 5000L
  n_alt    <- 8L
  n        <- n_strata * n_alt
  strata   <- rep(seq_len(n_strata), each = n_alt)

  # 50 clusters at the strata level -> per-row cluster
  n_clusters <- 50L
  cluster_per_stratum <- sample.int(n_clusters, n_strata, replace = TRUE)
  cluster <- cluster_per_stratum[strata]

  x1       <- stats::rnorm(n)
  decade   <- factor(sample(c("d1990", "d2000", "d2010"), n, replace = TRUE),
                     levels = c("d1990", "d2000", "d2010"))
  factor_a <- factor(sample(paste0("A", 1:6), n, replace = TRUE),
                     levels = paste0("A", 1:6))
  factor_b <- factor(sample(paste0("B", 1:4), n, replace = TRUE),
                     levels = paste0("B", 1:4))

  d <- data.frame(x1 = x1, decade = decade, factor_a = factor_a, factor_b = factor_b)
  # build design matrix the same way fclogit() expands -- via model.matrix
  # but drop intercept (clogit has none)
  mf <- stats::model.frame(~ x1 * decade + factor_a * decade + factor_b, data = d)
  X  <- stats::model.matrix(stats::terms(mf), mf)
  X  <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  p  <- ncol(X)

  beta_true <- stats::rnorm(p, sd = 0.3)
  eta       <- as.numeric(X %*% beta_true) + .gumbel(n)
  choice    <- .pick_choice(eta, strata)

  Xsp <- .densify_to_sparse(X)
  .verify_sparse(X, Xsp)

  list(
    X = X, Xsp = Xsp,
    choice = choice, strata = strata, cluster = cluster,
    beta_true = beta_true,
    meta = list(
      name = "medium_factor",
      formula = "~ x1*decade + factor_a*decade + factor_b",
      n_strata = n_strata, n_alt = n_alt, p = p, n_rows = n,
      density = mean(X != 0),
      n_clusters = n_clusters
    )
  )
}

# ----- 3. gen_paper3_like ----------------------------------------------------
#
# 20,000 strata x 50 alters (1M rows). Mirrors Paper 3 Step 4.
# Formula:  ~ pair*decade + age*decade + edu*decade + dist*decade
#   pair    25 levels (mimic pair_gen_meso); reference level is dominant
#           in 70% of choice events ("Sweden-Sweden" analog)
#   decade  3 levels
#   age     4 levels
#   edu     3 levels
#   dist    continuous
# Density ~5%. 4000 ego-clusters, ~5 obs/cluster.
# ----------------------------------------------------------------------------

gen_paper3_like <- function(seed = 3L) {
  set.seed(seed)
  n_strata <- 20000L
  n_alt    <- 50L
  n        <- n_strata * n_alt
  strata   <- rep(seq_len(n_strata), each = n_alt)

  # 4000 clusters with ~5 strata per cluster (= ~5 obs/cluster in the
  # ego sense -- one stratum per ego-time observation)
  n_clusters <- 4000L
  cluster_per_stratum <- sample.int(n_clusters, n_strata, replace = TRUE)
  cluster <- cluster_per_stratum[strata]

  # ---- pair: 25 levels, reference level dominant -----------------------
  # For 70% of strata, designate the "reference" level (pairP01) as the
  # chosen alternative *opportunity* by making it disproportionately likely
  # to be among the alternatives. We achieve this implicitly via beta_true:
  # set reference very attractive so it wins ~70% of strata. We still
  # sample alternatives uniformly.
  pair_levels <- paste0("pairP", sprintf("%02d", 1:25))
  pair <- factor(sample(pair_levels, n, replace = TRUE), levels = pair_levels)
  # Plant 1 alternative-per-stratum guaranteed to be the reference level
  # so the "Sweden-Sweden" can win most of the time.
  ref_row_per_stratum <- (strata - 1L) * n_alt + 1L
  pair[ref_row_per_stratum] <- pair_levels[1]

  decade <- factor(sample(c("d1990", "d2000", "d2010"), n, replace = TRUE),
                   levels = c("d1990", "d2000", "d2010"))
  age <- factor(sample(c("a1", "a2", "a3", "a4"), n, replace = TRUE),
                levels = c("a1", "a2", "a3", "a4"))
  edu <- factor(sample(c("e1", "e2", "e3"), n, replace = TRUE),
                levels = c("e1", "e2", "e3"))
  dist <- stats::rnorm(n)

  d <- data.frame(pair = pair, decade = decade, age = age, edu = edu, dist = dist)

  # Build piece by piece to avoid model.matrix memory blowups on 1M rows.
  # 1M x ~100 dense doubles = ~800 MB; tolerable but tight. Use model.matrix
  # column-by-column construction via sparse.model.matrix.
  X_sp_full <- Matrix::sparse.model.matrix(
    ~ pair * decade + age * decade + edu * decade + dist * decade,
    data = d
  )
  # drop intercept
  keep_cols <- colnames(X_sp_full) != "(Intercept)"
  X_sp_full <- X_sp_full[, keep_cols, drop = FALSE]
  Xsp <- methods::as(X_sp_full, "CsparseMatrix")
  X   <- as.matrix(Xsp)
  p   <- ncol(X)

  # beta_true: make pairP01 reference attractive -> it wins ~70% of strata.
  # Reference is excluded as baseline; positive coefficient on the *negation*
  # is implicit. Instead we draw all betas small except give "dist" a real
  # effect, and rely on Gumbel noise to give pair-reference its 70%.
  # Easier: keep betas modest and let geometry handle it.
  beta_true <- stats::rnorm(p, sd = 0.15)
  # Strong negative on all non-reference pair main effects -> reference wins
  pair_main_idx <- grep("^pair", colnames(X))
  pair_main_idx <- pair_main_idx[!grepl(":", colnames(X)[pair_main_idx])]
  beta_true[pair_main_idx] <- stats::runif(length(pair_main_idx), -2.0, -1.0)

  eta    <- as.numeric(X %*% beta_true) + .gumbel(n)
  choice <- .pick_choice(eta, strata)

  .verify_sparse(X, Xsp)

  list(
    X = X, Xsp = Xsp,
    choice = choice, strata = strata, cluster = cluster,
    beta_true = beta_true,
    meta = list(
      name = "paper3_like",
      formula = "~ pair*decade + age*decade + edu*decade + dist*decade",
      n_strata = n_strata, n_alt = n_alt, p = p, n_rows = n,
      density = mean(X != 0),
      n_clusters = n_clusters,
      ref_win_rate = mean(choice[ref_row_per_stratum])
    )
  )
}

# ----- 4. gen_edge_rare_cells ------------------------------------------------
#
# 2000 strata x 10 alters. Constructed so (rare_factor_level x decade=2010)
# has only ~30 chosen rows -- mimics the cell that triggers step-halving in
# Paper 3. Cluster present.
# ----------------------------------------------------------------------------

gen_edge_rare_cells <- function(seed = 4L) {
  set.seed(seed)
  n_strata <- 2000L
  n_alt    <- 10L
  n        <- n_strata * n_alt
  strata   <- rep(seq_len(n_strata), each = n_alt)

  n_clusters <- 200L
  cluster_per_stratum <- sample.int(n_clusters, n_strata, replace = TRUE)
  cluster <- cluster_per_stratum[strata]

  # rare_factor: 5 levels, level "rare" has prevalence ~3%
  rare_factor_probs <- c(rare = 0.03, common1 = 0.30, common2 = 0.30,
                         common3 = 0.27, common4 = 0.10)
  rare_factor <- factor(
    sample(names(rare_factor_probs), n, replace = TRUE, prob = rare_factor_probs),
    levels = names(rare_factor_probs)
  )
  decade <- factor(
    sample(c("d1990", "d2000", "d2010"), n, replace = TRUE),
    levels = c("d1990", "d2000", "d2010")
  )
  x1 <- stats::rnorm(n)

  d <- data.frame(rare_factor = rare_factor, decade = decade, x1 = x1)

  X_sp_full <- Matrix::sparse.model.matrix(
    ~ rare_factor * decade + x1,
    data = d
  )
  X_sp_full <- X_sp_full[, colnames(X_sp_full) != "(Intercept)", drop = FALSE]
  Xsp <- methods::as(X_sp_full, "CsparseMatrix")
  X   <- as.matrix(Xsp)
  p   <- ncol(X)

  beta_true <- stats::rnorm(p, sd = 0.4)
  eta       <- as.numeric(X %*% beta_true) + .gumbel(n)
  choice    <- .pick_choice(eta, strata)

  # Count the (rare, d2010) cell among chosen rows
  rare_2010_idx <- which(rare_factor == "rare" & decade == "d2010" & choice == 1L)
  n_rare_2010_chosen <- length(rare_2010_idx)

  .verify_sparse(X, Xsp)

  list(
    X = X, Xsp = Xsp,
    choice = choice, strata = strata, cluster = cluster,
    beta_true = beta_true,
    meta = list(
      name = "edge_rare_cells",
      formula = "~ rare_factor*decade + x1",
      n_strata = n_strata, n_alt = n_alt, p = p, n_rows = n,
      density = mean(X != 0),
      n_clusters = n_clusters,
      n_rare_2010_chosen = n_rare_2010_chosen
    )
  )
}

# ----- registry --------------------------------------------------------------

ALL_GENERATORS <- list(
  small_dense     = gen_small_dense,
  medium_factor   = gen_medium_factor,
  paper3_like     = gen_paper3_like,
  edge_rare_cells = gen_edge_rare_cells
)
