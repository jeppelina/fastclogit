# =============================================================================
# 08_assumptions.R: what the kernels assume, and what happens when it is false.
#
# The failure mode that matters is SILENCE: a violated assumption that produces
# a fit which converges cleanly to the wrong answer. Each probe below states
# the assumption, violates it, and records whether the package errors, warns,
# or proceeds quietly.
#
# A probe "passes" when the behaviour is safe (error or warning), and fails
# when the violation is accepted in silence. Some failures here are accepted
# and documented rather than fixed; the point is that the list is known.
# =============================================================================

vc_dir <- Sys.getenv("VC_DIR", unset = ".")
source(file.path(vc_dir, "00_helpers.R"))

vc_header("Assumption sweep")

base_data <- function(n_sets = 300L, n_alts = 8L, seed = 5L) {
  set.seed(seed)
  n <- n_sets * n_alts
  X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  st <- rep(seq_len(n_sets), each = n_alts)
  eta <- as.vector(X %*% c(0.8, -0.5))
  ch <- integer(n)
  for (s in split(seq_len(n), st))
    ch[s[sample.int(length(s), 1L, prob = exp(eta[s] - max(eta[s])))]] <- 1L
  list(X = X, choice = ch, strata = st, cluster = st)
}

# Classify what happened: "error", "warning", or "silent".
behaviour <- function(expr) {
  w <- NULL
  r <- withCallingHandlers(
    tryCatch(force(expr), error = function(e) structure("error", msg = conditionMessage(e))),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  if (identical(unclass(r)[1], "error")) return(list(kind = "error", msg = attr(r, "msg"), value = NULL))
  if (length(w)) return(list(kind = "warning", msg = w[1], value = r))
  list(kind = "silent", msg = "", value = r)
}

report <- function(label, b, safe_kinds = c("error", "warning"), note = "") {
  vc_note(sprintf("%-46s -> %-8s %s", label, b$kind,
                  substr(gsub("\\s+", " ", b$msg), 1, 60)))
  vc_check("assumptions", label, b$kind %in% safe_kinds,
           if (nzchar(note)) note else b$kind)
}

d <- base_data()

# --- 1. Exactly one chosen alternative per stratum ---------------------------
# The kernel breaks at the FIRST chosen == 1 it finds in a stratum, so a
# duplicated choice indicator fits silently on the first one and ignores the
# rest. A stratum with zero chosen is skipped via `continue`.
d2 <- d; d2$choice[which(d2$choice == 1L)[1] + 1L] <- 1L   # two chosen in set 1
report("two chosen alternatives in a stratum",
       behaviour(fastclogit(d2$X, d2$choice, d2$strata, max_iter = 50L)))

d3 <- d; d3$choice[which(d3$choice == 1L)[1]] <- 0L        # zero chosen in set 1
report("zero chosen alternatives in a stratum",
       behaviour(fastclogit(d3$X, d3$choice, d3$strata, max_iter = 50L)))

# --- 2. Cluster is constant within a stratum ---------------------------------
# The sandwich takes the cluster of each group's FIRST row. A stratum spanning
# two clusters is assigned wholly to one of them, so the robust SEs are
# computed for a clustering the user did not ask for. In the papers the
# stratum is an ego's choice set and the cluster is the ego, so this holds --
# but by convention, not by construction.
d4 <- d
d4$cluster <- rep(seq_len(length(d4$choice) / 4L), each = 4L)  # splits strata
report("cluster varies within a stratum",
       behaviour(fastclogit(d4$X, d4$choice, d4$strata, cluster = d4$cluster,
                            max_iter = 50L)))

# --- 3. Offset is finite -----------------------------------------------------
d5 <- d; off <- rep(0, length(d5$choice)); off[1] <- Inf
report("non-finite offset (Inf)",
       behaviour(fastclogit(d5$X, d5$choice, d5$strata, offset = off,
                            max_iter = 50L)))
off[1] <- NA_real_
report("NA in offset",
       behaviour(fastclogit(d5$X, d5$choice, d5$strata, offset = off,
                            max_iter = 50L)))

# --- 4. No NA in X -----------------------------------------------------------
d6 <- d; d6$X[1, 1] <- NA_real_
report("NA in the design matrix",
       behaviour(fastclogit(d6$X, d6$choice, d6$strata, max_iter = 50L)))

# --- 5. Singleton strata -----------------------------------------------------
# A stratum with one alternative contributes nothing to the gradient or
# Hessian, but is still counted in G for the finite-sample correction
# (G-1)/G. survival::clogit drops such strata entirely.
d7 <- base_data(n_sets = 300L)
sing <- which(d7$strata %in% 1:20)
keep <- c(sing[d7$choice[sing] == 1L], which(!(d7$strata %in% 1:20)))
d7 <- list(X = d7$X[keep, , drop = FALSE], choice = d7$choice[keep],
           strata = d7$strata[keep], cluster = d7$strata[keep])
b <- behaviour(fastclogit(d7$X, d7$choice, d7$strata, cluster = d7$cluster,
                          max_iter = 50L))
report("singleton strata (K = 1)", b, safe_kinds = c("error", "warning", "silent"),
       note = paste0(b$kind, "; n_groups reported = ",
                     if (!is.null(b$value)) b$value$n_groups else NA))

# --- 6. Scale equivariance ---------------------------------------------------
# A conditional logit is equivariant under rescaling of a covariate: scaling
# column j by c must divide beta_j by exactly c and leave every other
# coefficient and the log-likelihood unchanged. The Newton step uses a ridge
# of 1e-8 * max|diag(-H)|, which is relative and so should preserve this; if
# it does not, the ridge is doing something visible.
f1 <- fastclogit(d$X, d$choice, d$strata, max_iter = 100L, tol = 1e-10)
Xs <- d$X; Xs[, 1] <- Xs[, 1] * 1000
f2 <- fastclogit(Xs, d$choice, d$strata, max_iter = 100L, tol = 1e-10)
rel_b1 <- abs(stats::coef(f2)[1] * 1000 / stats::coef(f1)[1] - 1)
rel_b2 <- abs(stats::coef(f2)[2] / stats::coef(f1)[2] - 1)
rel_ll <- abs(f2$loglik / f1$loglik - 1)
vc_note(sprintf("rescale x1 by 1e3: beta1 rel err %.2e, beta2 rel err %.2e, loglik rel err %.2e",
                rel_b1, rel_b2, rel_ll))
vc_check("assumptions", "fit is equivariant under covariate rescaling",
         max(rel_b1, rel_b2, rel_ll) < 1e-6,
         sprintf("max relative error %.2e", max(rel_b1, rel_b2, rel_ll)))

# --- 7. Reproducibility ------------------------------------------------------
a1 <- fastclogit(d$X, d$choice, d$strata, max_iter = 100L)
a2 <- fastclogit(d$X, d$choice, d$strata, max_iter = 100L)
vc_check("assumptions", "identical input gives a bit-identical fit",
         identical(stats::coef(a1), stats::coef(a2)) &&
           identical(a1$loglik, a2$loglik) &&
           identical(a1$convergence_criterion, a2$convergence_criterion),
         "coefficients, loglik and route all identical")
