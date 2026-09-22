# fastclogit v0.5.0: validation, sweep and documentation plan

Drafted 2026-09-22, before pushing v0.5.0. Three phases: simulation runs, a
bug-and-assumptions sweep informed by what those runs show, then documentation
that cites the results rather than asserting them.

---

## 0. What the current tests do and do not establish

The suite as it stands (161 assertions) answers one question well: **does
fastclogit agree with `survival::clogit`?** It checks coefficients, log-likelihood
and SEs against survival on four small problems, plus dense/sparse agreement and
the defect regressions added in this merge.

That is a *software equivalence* check. Three things it cannot do:

1. **It cannot speak to regimes survival cannot reach.** Survival is skipped on
   the only problem above 40k rows. Every defect this package has had lived at
   scale, under an offset, with rare cells. Agreement at n=2,500 is weak evidence
   about n=37M.
2. **It never checks statistical properties.** Nothing tests that the estimator
   is unbiased, that the reported SEs match the sampling variability of the
   estimator, or that a nominal 95% interval covers 95% of the time. A
   systematically wrong SE would pass every current test, because survival would
   have to be wrong in the same way for the comparison to fail.
3. **It never tests the offset.** The McFadden-Manski correction is the single
   most load-bearing assumption in Papers 1, 3 and 4. There is no test anywhere
   that a fit on a stratified sample of alternatives with the correction
   recovers the parameters of the population choice model. The existing tests
   pass an offset through and check we agree with survival on the same offset,
   which tests plumbing, not correctness.

Point 3 is the most consequential gap in the package and phase 1 leads with it.

### A latent inconsistency worth resolving first

The repo contains **two different DGPs**:

- `R/simulate_clogit.R` draws the chosen alternative by **softmax sampling**.
- `tests/sparse_validation/gen_data.R` adds **Gumbel(0,1) noise to eta and takes
  the argmax**.

These are the same distribution (Gumbel-max is exactly categorical sampling from
the softmax), which is why nobody has noticed. But it is an unasserted
assumption linking the two harnesses, and if either is subtly wrong (a scale
error on the Gumbel, a missing offset in one path) then a whole family of tests
is validating against a mis-specified truth. Assert it, cheaply, once.

---

## Phase 1 — Simulation runs

Monte Carlo, not single fits. A single fit landing within 3 SEs of the truth is
consistent with both a correct estimator and a badly wrong one; only replication
separates them. All runs are scripted under `inst/validation/`, write results to
CSV, and are reproducible from a seed. Sizes are chosen so the whole phase runs
in well under an hour on this machine.

Every study states its **decision rule before it runs**. A simulation study whose
pass criterion is chosen after seeing the output is not evidence.

### S1. Recovery, SE calibration and coverage — the backbone

**Design.** R = 500 replications. n_egos = 2,000, n_alts = 20 (40k rows), the
default `simulate_clogit_data` DGP with continuous and factor predictors, no
offset, no clustering. Dense path.

**Estimands.** For each coefficient j: bias, the ratio of mean reported SE to the
Monte Carlo SD of the estimate, and the empirical coverage of the nominal 95%
interval.

**Decision rules, fixed in advance.**
- Bias: |mean(beta_hat_j) - beta_j| < 3 * MC standard error of the mean.
- SE calibration: mean(SE_j) / sd(beta_hat_j) within [0.97, 1.03].
- Coverage: within 0.95 +/- 2*sqrt(0.95*0.05/500) = [0.931, 0.969].

**What it catches that nothing currently does.** A wrong variance formula, a
missing finite-sample correction, a systematically biased estimator.

### S2. Does the McFadden-Manski offset actually work

**The highest-value study here.** Design mirrors what the papers do.

**Design.** Build a *population* choice set per ego: N_alt = 2,000 alternatives
drawn from five strata of known sizes. Compute the true population choice
probabilities under known beta and draw the chosen alternative from the full
population set. Then sample alternatives stratified, at deliberately unequal
rates (mirroring Paper 4: sampling fractions spanning ~8 log points), retaining
the chosen one. Fit three ways:

  (a) with the correct offset c_s = -log(n_s / N_s)
  (b) with **no** offset
  (c) with a deliberately **wrong** offset (strata permuted)

R = 300 replications, 1,000 egos, 25 sampled alternatives.

**Decision rules.**
- (a) must satisfy the S1 bias and coverage rules against the *population*
  beta.
- (b) must be detectably biased on at least one coefficient: |bias| > 5 * MC SE.
  If (b) is *not* biased the design has failed to create the problem the
  correction exists to solve, and the study is uninformative rather than
  reassuring — this is a design check, not a result.
- (c) must be biased. Same purpose: confirms the test is sensitive to the
  offset at all.

**Why (b) and (c) matter.** Without them, (a) passing could simply mean the
sampling design induced no bias. A validation that cannot fail is not a
validation.

### S3. Do the clustered sandwich SEs do their job

**Design.** Introduce genuine within-cluster dependence: each ego (cluster)
appears in several choice sets and carries a shared random intercept on a
covariate, so the per-stratum scores are positively correlated within cluster.
R = 400, 1,000 clusters, ~4 choice sets each.

**Decision rules.**
- Model-based SEs must **undercover**: coverage < 0.93. (If they do not, the
  dependence is too weak and the design is uninformative.)
- Robust SEs must recover nominal: coverage in [0.93, 0.97]. The band is wider
  than S1's because cluster-robust inference is only asymptotic in the number of
  clusters.
- Report the ratio mean(SE_robust)/sd(beta_hat) alongside.

**Secondary question this answers.** The kernels apply a finite-sample
correction `C/(C-1) * (G-1)/G`, described in the code as matching
`survival::coxph` with `cluster()`. That claim is currently untested. S3 checks
it empirically and against survival directly on a small case.

### S4. Dense versus sparse, and the open 271-loglik question

**The outstanding item from `FASTCLOGIT_MERGE_MAP.md` section 6.1.** Paper 3
observed the dense kernel converging to a log-likelihood 271 units *below*
sparse and survival at n=100 on real meso13 data, with the local reproduction
attempt failing at every scale up to 100k choice sets. The line-search fix is a
plausible explanation that was never tested, because the kernel-fix check fits
sparse only.

**Design.** A ladder of factor-heavy problems with the Paper-3 structure (a
dominant reference category holding ~70% of choice events, decade interactions,
rare cells) crossed with the Paper-4 offset structure, at n_alts in
{30, 50, 100} and scales up to the largest that fits locally. Both kernels on
identical input. Record coefficients, log-likelihood, convergence route, final
gradient, iteration count and the full `iter_log`.

**Decision rules.**
- Identified coefficients must agree to 1e-8 and log-likelihood to 1e-10
  relative.
- Convergence routes must match.
- Any divergence is recorded with its full trace and is a finding, not a
  failure to tune away.

**Honest statement of what this can and cannot settle.** The original
observation was on real data at n=100 with the *old* kernel, on Windows with a
different BLAS. If this reproduces the divergence, we have the cause. If it does
not, we have **not** cleared the dense path: we have only failed again to
reproduce it locally, exactly as in June. The note must say so, and the
conclusion for the papers stays "prefer sparse at production scale" either way.
Nothing in phase 1 can settle the question without a MONA run on real data.

### S5. Convergence-route audit

The papers rely on the claim that a non-primary convergence route is still
trustworthy. That claim has never been tested systematically.

**Design.** A stress grid designed to *force* each route: rare cells at varying
frequency, offsets of varying within-stratum spread, near-separation, and
deliberately ill-conditioned designs (near-collinear columns). For every fit,
record the route, then verify the fit is genuinely at a maximum by refitting from
a perturbed start with `tol = 1e-12` and comparing.

**Decision rule.** For every fit reporting `converged = TRUE` by any route, the
log-likelihood must be within 1e-6 (relative) of the high-precision refit, and
max|gradient| must be below `tier3_grad_floor`. Any route that fails this is a
route we should stop trusting, and that is a substantive finding about the
papers, not just about the package.

### S6. DGP equivalence

Assert softmax-sampling and Gumbel-argmax produce the same choice distribution,
at one configuration, with a chi-square goodness-of-fit test on realised choice
frequencies over R = 200,000 draws. Cheap, and it licenses comparing results
across the two harnesses.

### S7. Scale, memory and timing

Feeds phase 3's benchmark documentation. Dense and sparse across a grid of
n_rows, p and density; record wall time and peak memory (`gc()` based), plus
`survival::clogit` where it is feasible so the comparison is real rather than
quoted from memory. The README currently claims "225 GB to 30 GB" and "95
minutes to 80 seconds" from Paper 3 production runs; those are real but they are
not reproducible by a reader. Publish numbers a reader can reproduce, and cite
the production figures separately as what they are.

---

## Phase 2 — Bug and assumptions sweep

Run **after** phase 1, because the simulation output tells us where to look.
Two parts.

### 2a. Enumerate and test the implicit assumptions

The kernels assume a great deal that is never checked. Each item gets a test
that either confirms the guard exists or demonstrates the failure mode.

| # | Assumption | Failure mode if violated | Currently guarded? |
|---|---|---|---|
| 1 | Exactly one chosen alternative per stratum | kernel `continue`s on 0; 2 chosen is silently wrong (first wins) | partially |
| 2 | Strata contiguous after sorting | wrong group boundaries, silent | by construction, untested |
| 3 | Cluster is constant within stratum | kernel takes the cluster of each group's FIRST row; a stratum spanning clusters is silently miscounted | no |
| 4 | Offset is finite | NaN propagates to the whole fit | no |
| 5 | No NA in X | silent NaN | no |
| 6 | group_size >= 2 | singleton strata contribute nothing but are counted in G, distorting the finite-sample correction | no |
| 7 | nnz < 2^31 | overflow | yes, added in this merge |
| 8 | n * p does not overflow int | fixed in the verbose path only | partial |
| 9 | No stratum-constant column | not identified; ridge-determined value with meaningless SE | **no — known, documented** |

Three of these are worth singling out as candidate *real* bugs rather than
documentation gaps:

- **#3 (cluster within stratum).** If any stratum spans two clusters, the
  sandwich assigns the whole stratum to one of them. In the papers a stratum is
  an ego's choice set and the cluster is the ego, so this holds — but it holds
  by *convention*, not by construction, and nothing checks it. A user who
  clusters on something coarser than the stratum (say, county) gets silently
  wrong robust SEs.
- **#6 (singleton strata).** A stratum with one alternative contributes zero to
  the gradient and Hessian but is still counted in `G` for the correction
  `(G-1)/G`. `survival::clogit` drops such strata. If we count them, our robust
  SEs differ from survival's by a factor we have never measured.
- **#1 (two chosen).** The kernel breaks at the first `chosen == 1` it finds.
  Data with a duplicated choice indicator fits silently and wrongly.

### 2b. Scale invariance and tolerance semantics

Two methodological issues, not coding errors, which the simulation studies are
well placed to expose:

- **The ridge is scale-dependent.** `ridge = 1e-8 * max|diag(-H)|` is relative,
  which is right, but the *fallback* `1e-4 * max|diag|` is large enough to move
  a genuinely near-singular fit. Test: rescale a covariate by 1e3 and check the
  fitted coefficients transform exactly as they should. A conditional logit is
  equivariant under rescaling of X; if our fit is not, the ridge is doing
  something visible.
- **`tol` is on max|gradient|, which grows with n.** The gradient is a sum over
  strata, so at 37M rows a gradient of 1e-3 may be closer to the maximum than
  1e-6 is at 2,500 rows. This means `tol = 1e-6` is a *different* criterion at
  different scales, which is very likely part of why the convergence ladder
  needed three tiers in the first place. Test the relationship empirically
  across the S7 scale ladder, and document it. A gradient normalised by n (or by
  the number of strata) would be scale-free; whether to change the default is a
  decision for after the evidence, not now.

### 2c. Code review of the merge

Run the `code-review` skill over the v0.5.0 diff at high effort, as an
independent pass over what phases 1 and 2 do not reach.

---

## Phase 3 — Documentation

Written last, so it can cite phases 1 and 2 instead of asserting. Nothing goes
into the docs that a phase-1 script did not produce.

1. **README rewrite.** What it is, when to use it over `survival::clogit`, when
   *not* to, installation, a 10-line quick start, the reproducible benchmark
   table from S7, and an explicit limitations section (stratum-constant columns,
   the open dense question, tolerance semantics at scale).
2. **Complete options reference** in `?fastclogit`: every argument, what it
   does, when to change it, and what each convergence route means operationally.
   The six `tier3_*` arguments are currently documented individually but with no
   guidance on when anyone would touch them.
3. **Vignette 1, "Getting started"** — update the existing one for v0.5.0
   (convergence reporting, `$iter_log`, `$conf_int`).
4. **Vignette 2, "Convergence and diagnostics"** — NEW, and the one this package
   most needs. What each route means, how to read `$iter_log`, how to tell a
   finished fit from a stuck one, what to do about each. Built directly on the
   S5 stress grid, so every example is a real fit.
5. **Vignette 3, "Large-scale and restricted environments"** — the sparse path,
   when it pays, memory figures from S7, and the MONA source-mode workflow with
   the generated bundle.
6. **`inst/validation/README.md`** — how to re-run every study, what each
   decision rule was, and the results as of v0.5.0. This is what makes the
   claims checkable rather than promotional.

---

## Sequencing and stopping rules

1. S6 first (cheap, and licenses everything else).
2. S1, S2, S3 — the statistical backbone. **If S2(a) fails, stop and escalate.**
   That would mean the offset correction does not do what three papers assume it
   does, and no amount of documentation is the right next action.
3. S4, S5 — the kernel questions.
4. S7 — benchmarks.
5. Phase 2 sweep, informed by 1-4.
6. Phase 3 docs, citing 1-5.
7. Re-run the full testthat suite, then decide on the push.

**Explicit non-goal.** None of this settles the dense-path question on real
data; only a MONA run does. The plan is designed to say clearly what it has and
has not established.
