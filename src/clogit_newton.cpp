// clogit_newton.cpp — Newton-Raphson conditional logit with offset support
//
// Memory-efficient: works on pre-built design matrix, no copies.
// Each group's contribution to gradient/Hessian is accumulated in-place.
//
// Convergence ladder (added 2026-06-17 to bring fastclogit closer to
// survival::clogit on flat-ridge problems while staying strict enough to
// refuse clearly-unconverged fits):
//   PRIMARY:   max|grad| < tol                                   (code = 1)
//   SECONDARY: rel_ll_change < tol*0.01 AND grad_max < tol*10
//              AND prev_unhalved_step_norm < tol*1e3             (code = 2)
//   PLATEAU:   rel_ll_change < tier3_plateau_tol for K iters
//              AND grad_max < tier3_grad_floor
//              AND (prev_halving_count >= tier3_halving_floor OR
//                   prev_step_size < tier3_step_floor)           (code = 3)
//   iter_max:  none of the above before max_iter                 (code = 4)
//
// PLATEAU matches survival::clogit's rel-ll criterion pointwise (1e-9 default)
// but adds three guardrails (consecutive-iter persistence, grad floor,
// optimizer-struggling evidence) so it cannot fire on clearly-unconverged
// problems where the gradient is still large.
//
// SECONDARY synced with sparse kernel 2026-06-17: tightened from
// grad_max < tol*1e4 to grad_max < tol*10 and added the unhalved-Newton-step
// guard, so dense and sparse kernels converge under the same criteria.
//
// Best-loglik tracking (2026-06-17): tracks beta_at_max_loglik across iters
// and returns it as $coefficients, with $best_loglik_iter recording which
// iter produced it. Survival does the same; matters only when step-halving
// overshoots at the tail.
//
// Per-iter trace: returned as flat vectors (iter_log_*), assembled into a
// data.table on the R side.
//
// Other features:
//   - LogSumExp trick for numerical stability in softmax
//   - Adaptive ridge regularization (Levenberg-Marquardt) for near-singular
//     Hessians: ridge = 1e-8 * max|diag(-H)|, with 1e-4 fallback
//   - Step-halving line search (up to 20 halvings) to ensure monotone LL
//   - Recomputes Hessian at best-loglik beta for accurate variance estimation
//
// Called by: fastclogit() in fastclogit.R (via Rcpp::sourceCpp or package)
//
// Author: Jesper Lindmarker
// License: MIT

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List clogit_fit_cpp(
    const arma::mat& X,            // n x p design matrix
    const arma::ivec& chosen,      // n x 1 binary (1 = chosen alternative)
    const arma::vec& offset,       // n x 1 offset (McFadden/Manski correction)
    const arma::ivec& group_start, // G x 1 (0-based start index of each group)
    const arma::ivec& group_size,  // G x 1 (number of rows in each group)
    int max_iter,                  // maximum Newton-Raphson iterations
    double tol,                    // tier-1 / tier-2 tolerance on max|gradient|
    bool tier3_enable,             // enable loglik-plateau (tier-3) criterion
    double tier3_plateau_tol,      // rel-ll-change ceiling for plateau (default 1e-9)
    int tier3_plateau_iters,       // consecutive plateau iters required (default 3)
    double tier3_grad_floor,       // |grad| sanity ceiling for plateau (default 1e-2)
    int tier3_halving_floor,       // halvings in prev iter required (default 2)
    double tier3_step_floor,       // OR step_size < this in prev iter (default 1e-3)
    bool verbose                   // print iteration info
) {
    int p = X.n_cols;
    int G = group_start.n_elem;

    arma::vec beta(p, arma::fill::zeros);
    arma::vec grad(p, arma::fill::zeros);
    arma::mat hess(p, p, arma::fill::zeros);
    double loglik = 0.0;
    int iter = 0;
    bool converged = false;
    int convergence_code = 4;  // default = iter_max

    // Best-loglik beta tracking — returned as $coefficients and used for the
    // final Hessian recomputation.
    arma::vec best_beta = beta;
    double    best_loglik = -std::numeric_limits<double>::infinity();
    int       best_loglik_iter = 0;

    // Tier-3 state
    int    plateau_count        = 0;
    int    prev_halving_count   = 0;
    double prev_step_size       = 1.0;
    double prev_newton_step_norm = std::numeric_limits<double>::infinity();

    // Per-iter trace (flat vectors; assembled into a data.table on R side)
    std::vector<int>    log_iter;
    std::vector<double> log_loglik;
    std::vector<double> log_grad_max;
    std::vector<double> log_rel_ll_change;
    std::vector<double> log_step_size;
    std::vector<int>    log_halving_count;
    std::vector<int>    log_tier_fired;   // 0=none, 1=primary, 2=secondary, 3=plateau
    log_iter.reserve(max_iter);
    log_loglik.reserve(max_iter);
    log_grad_max.reserve(max_iter);
    log_rel_ll_change.reserve(max_iter);
    log_step_size.reserve(max_iter);
    log_halving_count.reserve(max_iter);
    log_tier_fired.reserve(max_iter);

    for (iter = 0; iter < max_iter; iter++) {
        double loglik_new = 0.0;
        grad.zeros();
        hess.zeros();

        for (int j = 0; j < G; j++) {
            int start = group_start(j);
            int K = group_size(j);

            const arma::mat Xj = X.rows(start, start + K - 1);
            const arma::vec oj = offset.subvec(start, start + K - 1);

            arma::vec eta = Xj * beta + oj;

            double max_eta = eta.max();
            arma::vec exp_eta = arma::exp(eta - max_eta);
            double sum_exp = arma::accu(exp_eta);
            arma::vec prob = exp_eta / sum_exp;

            int c_idx = -1;
            for (int k = 0; k < K; k++) {
                if (chosen(start + k) == 1) {
                    c_idx = k;
                    break;
                }
            }

            if (c_idx < 0) {
                Rcpp::warning("Group %d has no chosen alternative — skipping", j + 1);
                continue;
            }

            loglik_new += eta(c_idx) - (max_eta + std::log(sum_exp));

            arma::vec xbar = Xj.t() * prob;
            grad += Xj.row(c_idx).t() - xbar;
            arma::mat Xj_w = Xj.each_col() % arma::sqrt(prob);
            hess -= (Xj_w.t() * Xj_w - xbar * xbar.t());
        }

        double grad_max      = arma::abs(grad).max();
        double abs_ll_change = (iter > 0) ? std::abs(loglik_new - loglik) : 0.0;
        double rel_ll_change = (iter > 0)
            ? abs_ll_change / (std::abs(loglik) + 1e-10) : 0.0;

        if (verbose) {
            Rprintf("Iter %d: loglik = %.8f, max|grad| = %.2e\n",
                    iter + 1, loglik_new, grad_max);
        }

        // Best-loglik tracking. beta entered THIS iter, loglik_new is its
        // loglik. Always record on iter 0 to seed best_beta off the initial
        // beta=0 if subsequent iters somehow fail to improve.
        if (loglik_new > best_loglik) {
            best_loglik      = loglik_new;
            best_beta        = beta;
            best_loglik_iter = iter + 1;
        }

        int tier_fired = 0;

        // PRIMARY — strict gradient test
        if (grad_max < tol) {
            tier_fired = 1;
            convergence_code = 1;
            converged = true;
        }
        // SECONDARY — rel_ll AND grad AND step all tight (synced with sparse)
        else if (iter > 0 &&
                 rel_ll_change          < tol * 0.01 &&
                 grad_max               < tol * 10.0 &&
                 prev_newton_step_norm  < tol * 1e3) {
            tier_fired = 2;
            convergence_code = 2;
            converged = true;
            if (verbose) {
                Rprintf("  Converged via secondary: rel_ll=%.2e, grad=%.2e, prev_step=%.2e\n",
                        rel_ll_change, grad_max, prev_newton_step_norm);
            }
        }
        // PLATEAU — survival-style loglik plateau with side conditions
        else if (iter > 0 && tier3_enable &&
                 rel_ll_change < tier3_plateau_tol &&
                 grad_max      < tier3_grad_floor &&
                 (prev_halving_count >= tier3_halving_floor ||
                  prev_step_size      <  tier3_step_floor)) {
            plateau_count++;
            if (plateau_count >= tier3_plateau_iters) {
                tier_fired = 3;
                convergence_code = 3;
                converged = true;
                if (verbose) {
                    Rprintf("  Converged via plateau (tier-3): rel_ll=%.2e for %d iters, grad=%.2e, prev_halvings=%d\n",
                            rel_ll_change, plateau_count, grad_max,
                            prev_halving_count);
                }
            }
        } else {
            plateau_count = 0;
        }

        // Log this iter. step_size and halving_count describe the step that
        // produced the current beta (prev iter's Newton step). On iter 0 they
        // are 1.0 / 0 (no prior step taken).
        log_iter.push_back(iter + 1);
        log_loglik.push_back(loglik_new);
        log_grad_max.push_back(grad_max);
        log_rel_ll_change.push_back((iter > 0) ? rel_ll_change : NA_REAL);
        log_step_size.push_back((iter > 0) ? prev_step_size : NA_REAL);
        // -1, not NA_INTEGER. log_halving_count is returned as an arma::ivec,
        // and with ARMA_64BIT_WORD that is a 64-bit type which reaches R as a
        // double. NA_INTEGER is INT_MIN, which then sits outside R's integer
        // range, so as.integer() on the R side produced the right NA but also a
        // spurious "NAs introduced by coercion to integer range" warning on
        // every single fit. -1 round-trips cleanly and is mapped to NA in
        // fastclogit().
        log_halving_count.push_back((iter > 0) ? prev_halving_count : -1);
        log_tier_fired.push_back(tier_fired);

        if (converged) {
            loglik = loglik_new;
            iter++;  // returned iter count is number-of-iters-completed (1-based)
            break;
        }

        loglik = loglik_new;

        // Newton step: delta = -H^{-1} * grad
        // Adaptive ridge regularization for near-singular Hessians.
        arma::mat neg_hess = -hess;
        double ridge = 1e-8 * arma::abs(neg_hess.diag()).max();
        if (ridge < 1e-12) ridge = 1e-12;
        neg_hess.diag() += ridge;

        arma::vec delta;
        bool solve_ok = arma::solve(delta, neg_hess, grad,
                                     arma::solve_opts::likely_sympd);
        if (!solve_ok) {
            neg_hess.diag() += 1e-4 * arma::abs(neg_hess.diag()).max();
            solve_ok = arma::solve(delta, neg_hess, grad);
        }
        if (!solve_ok) {
            Rcpp::warning("Hessian is singular at iteration %d", iter + 1);
            break;
        }

        // Record UNHALVED Newton step magnitude — used by NEXT iter's
        // secondary convergence check.
        prev_newton_step_norm = arma::abs(delta).max();

        // Step-halving: ensure log-likelihood does not decrease
        double step_size = 1.0;
        arma::vec beta_new = beta + step_size * delta;
        int halving_count = 0;
        int max_halving = 20;

        // LINE SEARCH. Two changes, 2026-09-11, after models 2-4 of Paper 4 ran
        // all 100 iterations without converging.
        //
        // (1) It now runs on the FIRST iteration too. It used to be skipped
        //     there, so the opening Newton step was accepted unbounded however
        //     large. From beta = 0 with an offset spanning ~8 log points, that
        //     step had norm ~900 and drove the log-likelihood from -5.5e6 to
        //     -4.1e8, into a region where the softmax is numerically degenerate:
        //     the Hessian collapses toward zero, the gradient freezes, and the
        //     Newton direction stops being an ascent direction. Nothing after
        //     that can recover.
        //
        // (2) A step that improves nothing is no longer accepted. When the 20
        //     halvings were exhausted without finding an improvement, beta_new
        //     was assigned anyway, so the optimiser walked DOWNHILL by a tiny
        //     amount every iteration: 90 consecutive iterations each losing
        //     ~0.14 of log-likelihood, at roughly 70 seconds apiece. It now
        //     stops and reports the failure instead of grinding to the cap.
        bool ls_improved = false;
        {
            while (halving_count < max_halving) {
                double ll_candidate = 0.0;
                for (int j = 0; j < G; j++) {
                    int start = group_start(j);
                    int K = group_size(j);
                    const arma::mat Xj = X.rows(start, start + K - 1);
                    const arma::vec oj = offset.subvec(start, start + K - 1);
                    arma::vec eta_c = Xj * beta_new + oj;
                    double max_eta_c = eta_c.max();
                    double lse = max_eta_c + std::log(arma::accu(arma::exp(eta_c - max_eta_c)));

                    int c_idx = -1;
                    for (int k = 0; k < K; k++) {
                        if (chosen(start + k) == 1) { c_idx = k; break; }
                    }
                    if (c_idx >= 0) {
                        ll_candidate += eta_c(c_idx) - lse;
                    }
                }

                // Acceptance tolerance scaled to the magnitude of the
                // log-likelihood. A fixed 1e-10 is BELOW one unit in the last
                // place once |loglik| passes about 4.5e5: at Paper 4's
                // -4.26e6 an ulp is 9.5e-10, so near the optimum this test was
                // comparing rounding noise and could never register an
                // improvement. That is what produced 20 halvings on every
                // iteration in fits that were already finished, and, once the
                // 2026-09-11 guard was added, what made them report failure.
                const double ll_tol = 1e-10 + 8.0 * std::fabs(loglik) * 2.220446049250313e-16;
                if (ll_candidate >= loglik - ll_tol) {
                    ls_improved = true;
                    break;
                }

                step_size *= 0.5;
                beta_new = beta + step_size * delta;
                halving_count++;
            }

            if (halving_count > 0 && verbose) {
                Rprintf("  Step-halving: %d halvings, step_size = %.4e\n",
                        halving_count, step_size);
            }
        }

        // Save for next iter's tier-3 check and iter-log row.
        prev_halving_count = halving_count;
        prev_step_size     = step_size;

        if (!ls_improved) {
            // A line search that cannot improve is ambiguous, and the gradient
            // tells the two cases apart.
            //
            // AT a maximum, no step improves the log-likelihood: that is what
            // being at a maximum means. The first version of this guard
            // (2026-09-11) broke unconditionally and reported failure, which
            // mislabelled fits that were finished. On 2026-09-12 that cost a
            // men's KHB decomposition: M1 stopped here with the log-likelihood
            // flat to twelve significant figures and max|grad| 8.4e-04, was
            // renamed FAILED, and the decomposition was skipped for the sex.
            //
            // AWAY from a maximum the gradient is still large, the Newton
            // direction has stopped being an ascent direction, and continuing
            // only walks downhill. That is the Paper 4 M2 case: gradient frozen
            // at 3.95e+06.
            //
            // tier3_grad_floor is the same ceiling the plateau rule uses for
            // exactly this judgement, so the two agree by construction.
            if (grad_max < tier3_grad_floor) {
                converged        = true;
                convergence_code = 6;
                if (verbose) {
                    Rprintf("  Line search cannot improve and max|grad| = %.2e is below the "
                            "floor (%.2e).\n  At a flat optimum; treating as converged.\n",
                            grad_max, tier3_grad_floor);
                }
            } else {
                convergence_code = 5;
                if (verbose) {
                    Rprintf("  Line search exhausted %d halvings without improving the "
                            "log-likelihood, and max|grad| = %.2e is still large.\n  Stopping: "
                            "the Newton direction is not an ascent direction here.\n",
                            max_halving, grad_max);
                }
            }
            break;
        }

        beta = beta_new;
    }

    // If we exited via iter_max (loop counter reached max_iter), record code 4.
    if (!converged) {
        convergence_code = 4;
    }

    // Recompute Hessian (and grad, loglik) at best_beta for variance estimation.
    arma::vec beta_final = best_beta;
    grad.zeros();
    hess.zeros();
    double loglik_final = 0.0;
    for (int j = 0; j < G; j++) {
        int start = group_start(j);
        int K = group_size(j);
        const arma::mat Xj = X.rows(start, start + K - 1);
        const arma::vec oj = offset.subvec(start, start + K - 1);
        arma::vec eta = Xj * beta_final + oj;
        double max_eta = eta.max();
        arma::vec exp_eta = arma::exp(eta - max_eta);
        double sum_exp = arma::accu(exp_eta);
        arma::vec prob = exp_eta / sum_exp;

        int c_idx = -1;
        for (int k = 0; k < K; k++) {
            if (chosen(start + k) == 1) { c_idx = k; break; }
        }
        if (c_idx < 0) continue;

        loglik_final += eta(c_idx) - (max_eta + std::log(sum_exp));
        arma::vec xbar = Xj.t() * prob;
        grad += Xj.row(c_idx).t() - xbar;
        arma::mat Xj_w = Xj.each_col() % arma::sqrt(prob);
        hess -= (Xj_w.t() * Xj_w - xbar * xbar.t());
    }

    // Final variance estimate. Rank-deficient designs (aliased columns,
    // separation-on-rare-cells) produce a singular -H, so inv_sympd fails
    // even at the MLE. Mirror the in-iter solve's adaptive ridge here:
    // tiny ridge first, then a heavier ridge + pinv fallback. Cells where
    // the ridge actually load-bears are flagged via vcov_singular = true
    // so downstream can refuse to trust the corresponding SEs.
    arma::mat neg_H_final = -hess;
    double ridge_final = 1e-8 * arma::abs(neg_H_final.diag()).max();
    if (ridge_final < 1e-12) ridge_final = 1e-12;
    neg_H_final.diag() += ridge_final;
    arma::mat vcov;
    bool vcov_singular = false;
    if (!arma::inv_sympd(vcov, neg_H_final)) {
        vcov_singular = true;
        neg_H_final.diag() += 1e-4 * arma::abs(neg_H_final.diag()).max();
        if (!arma::inv_sympd(vcov, neg_H_final)) {
            vcov = arma::pinv(neg_H_final);
            Rcpp::warning("Final Hessian singular even after ridge; using pseudoinverse for vcov.");
        }
    }

    // Wrap trace vectors as Armadillo / NumericVector for return.
    arma::ivec log_iter_arma         = arma::conv_to<arma::ivec>::from(log_iter);
    arma::vec  log_loglik_arma       = arma::conv_to<arma::vec >::from(log_loglik);
    arma::vec  log_grad_max_arma     = arma::conv_to<arma::vec >::from(log_grad_max);
    arma::vec  log_rel_ll_change_arma = arma::conv_to<arma::vec>::from(log_rel_ll_change);
    arma::vec  log_step_size_arma    = arma::conv_to<arma::vec >::from(log_step_size);
    arma::ivec log_halving_count_arma = arma::conv_to<arma::ivec>::from(log_halving_count);
    arma::ivec log_tier_fired_arma   = arma::conv_to<arma::ivec>::from(log_tier_fired);

    return Rcpp::List::create(
        Rcpp::Named("coefficients")           = beta_final,
        Rcpp::Named("vcov")                   = vcov,
        Rcpp::Named("vcov_singular")          = vcov_singular,
        Rcpp::Named("hessian")                = hess,
        Rcpp::Named("loglik")                 = loglik_final,
        Rcpp::Named("iterations")             = iter,
        Rcpp::Named("converged")              = converged,
        Rcpp::Named("convergence_code")       = convergence_code,
        Rcpp::Named("gradient")               = grad,
        Rcpp::Named("best_loglik_iter")       = best_loglik_iter,
        Rcpp::Named("iter_log_iter")          = log_iter_arma,
        Rcpp::Named("iter_log_loglik")        = log_loglik_arma,
        Rcpp::Named("iter_log_grad_max")      = log_grad_max_arma,
        Rcpp::Named("iter_log_rel_ll_change") = log_rel_ll_change_arma,
        Rcpp::Named("iter_log_step_size")     = log_step_size_arma,
        Rcpp::Named("iter_log_halving_count") = log_halving_count_arma,
        Rcpp::Named("iter_log_tier_fired")    = log_tier_fired_arma
    );
}
