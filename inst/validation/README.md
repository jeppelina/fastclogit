# Validation studies

Simulation evidence for the claims this package makes. Every study states its
decision rule **before** it runs, as a call to `vc_check()`, and the rules live
in the scripts, in git. Changing a rule after seeing a result shows up as a
diff, so the plan is its own pre-registration.

Run everything:

```bash
cd inst/validation && Rscript run_all.R
```

or one study: `Rscript run_all.R 04_offset`. Total runtime is about 7.5 minutes.

## Why these exist

The `testthat` suite answers a *software* question: does fastclogit agree with
`survival::clogit`? That is necessary and not sufficient. It says nothing about
regimes survival cannot reach, it never checks statistical properties, and a
systematically wrong variance formula would pass every test in it, because
survival would have to be wrong in the same way for the comparison to fail.

These studies check the statistics.

## Results, v0.5.0 (2026-09-22): 35 of 36 checks pass

| Study | Checks | What it establishes |
|---|---|---|
| `01_dgp_equivalence` | 4/4 | The two DGPs in this repo (softmax sampling, Gumbel argmax) are the same distribution, with and without a wide offset |
| `02_design_matrix` | 3/3 | The formula path and the explicit-matrix path give identical fits over 80 random formulas; rare columns survive the collinearity screen |
| `03_recovery_coverage` | 3/3 | Unbiased; reported SEs match sampling variability; 95% intervals cover at 0.941–0.957 |
| `04_offset` | 5/5 | The McFadden-Manski correction recovers the population parameters, and three falsification arms are biased by 126–152 MC SE |
| `05_cluster_robust` | 5/5 | Model SEs undercover at 0.65–0.68 under clustering; robust SEs restore 0.935–0.955 and match survival to 0.2% |
| `06_khb` | 2/2 | The KHB decomposition recovers a known mediation structure; the total = direct + indirect identity holds exactly |
| `07_benchmarks` | 3/4 | Three engines agree on the log-likelihood exactly; sparse is 5–35x faster; one open finding, below |
| `08_assumptions` | 9/9 | Nine assumption violations all produce an error or a warning, none are silent |

### Headline numbers

**Coverage** (R = 2,000, 2,000 strata x 20 alternatives, 11 coefficients):
every coefficient covers within [0.941, 0.957] against a pre-specified band of
[0.940, 0.960]. Mean reported SE over empirical SD ranges 0.971 to 1.008.

**The offset** (R = 300, population choice set of 1,500, 25 alternatives
sampled at rates spanning 0.60 to 0.01):

| Arm | Bias on `x` (truth 1.0) | Coverage |
|---|---|---|
| (a) correct offset | +0.0016 (0.5 MC SE) | 0.970 |
| (b) no offset | +0.469 (135 MC SE) | 0.000 |
| (c) permuted offset | +0.678 (152 MC SE) | 0.000 |
| (d) chosen alternative misassigned | +0.765 (126 MC SE) | 0.000 |

Arms (b), (c) and (d) are supposed to fail. A validation that cannot fail is
not a validation, and without them arm (a) passing could just mean the sampling
design induced no bias in the first place.

**Cluster-robust SEs** under exact 4-fold duplication, where the answer is
known analytically: model SEs are too small by 1.998 and 1.996 against a true
factor of exactly 2.000, and the sandwich recovers it.

**Benchmarks** (this machine, per-process peak RSS via `/usr/bin/time -l`,
because `gc()` cannot see Armadillo's allocations):

| rows | p | dense | sparse | survival |
|---|---|---|---|---|
| 100,000 | 20 | 0.33 s / 0.30 GB | 0.05 s / 0.25 GB | 1.06 s / 0.40 GB |
| 500,000 | 50 | 4.32 s / 1.24 GB | 0.26 s / 0.41 GB | — |
| 1,000,000 | 90 | 21.34 s / 2.91 GB | 0.60 s / 0.64 GB | — |

All engines reach the same log-likelihood to a relative spread of exactly 0.

## What these studies found

Four defects, all fixed in v0.5.0:

1. **The collinearity screen deleted estimable columns.** `fclogit()` decided
   which columns to drop by running QR on a 50,000-row subsample. A dummy whose
   few 1s all fell outside the subsample looked like a zero column and was
   silently removed. Measured at 200,000 rows: a dummy with **4 ones was dropped
   in 40% of runs**, one with 8 ones in 4%. That is precisely the rare-cell
   regime these models exist for. Now 0% at every support level tested.

2. **`khb_decompose()` errored on its own documented default.** `controls`
   defaults to `NULL` and was passed straight to `strsplit()`. The roxygen
   example always supplies controls, so the default path had never been run.

3. **Three input violations failed deep inside the kernel.** An `NA` anywhere
   in `X` or the offset surfaced as `pinv(): svd failed`, preceded by Armadillo
   warnings about a non-symmetric matrix, which tells a user nothing about the
   one bad cell in their data. Now caught up front by name and count.

4. **Strata spanning several clusters were accepted silently.** The sandwich
   assigns each stratum to the cluster of its first row, so the robust SEs were
   for a clustering the caller never asked for. Now a warning.

## Open finding

`07_benchmarks` reports one deliberate failure. When `tol` is set below the
attainable numerical floor, the optimiser grinds to `max_iter` and reports
**not converged** at a gradient that is already tiny: at 10,000 strata with
`tol = 1e-14` it ran all 200 iterations and finished at `max|grad| = 3.8e-13`.
The plateau rule does not rescue it, because that rule's side-condition
requires evidence of struggling (halvings, or a tiny step) that a cleanly
converged fit never produces.

v0.5.0 makes the warning say so explicitly rather than changing the criterion.
Changing convergence behaviour under three live papers is not a patch-time
decision.

## A hypothesis these studies refuted

The plan predicted that the attainable gradient floor **rises** with n, making
an absolute `tol` progressively harder to reach and explaining why the
convergence ladder needed extra tiers. At the scales reachable here that is
simply not true: the floor measured 8e-14 to 3e-9 over 2,000 to 200,000 strata,
is not monotone in n, and sits far below the 1e-6 default.

The kernel comments record a floor of ~0.003 at 75M+ rows. That is three orders
of magnitude beyond anything reproducible on one machine, so it can be neither
confirmed nor refuted here, and the documentation attributes it to production
experience rather than to this study. No exponent is fitted to four
non-monotone points.

## Two simulation designs that were wrong, and why

Both produced results that looked fine, which is the reason they are recorded
rather than quietly corrected.

**The offset study was wrong twice.** The McFadden-Manski correction formula is
*protocol-specific*. The first version counted the chosen alternative, which is
retained with certainty rather than sampled, in the stratum's draw count; that
biased the estimator by +0.10 on a true coefficient of 1.0, with 57% coverage.
The second forced the chosen alternative in *on top of* `n_s` draws per
stratum, which is a different sampling protocol from the one `-log(n_s/N_s)` is
derived for. Only when the chosen alternative occupies one of its stratum's
`n_s` slots does the familiar formula apply. Both errors were in the
simulation, not the package — but Paper 1 and Paper 4 each carry a patch script
for exactly this bookkeeping, which is some evidence about how easy it is to
get wrong.

**Two natural-looking clustering designs induce no dependence at all.** A
random intercept shared within a choice set cancels exactly in the conditional
likelihood, because conditional logit conditions on the stratum. Correlating
covariates across an ego's several choice sets also does nothing, and this one
is subtler: the per-stratum score has conditional mean zero given `X`, so
`Cov(s_1, s_2) = E[Cov(s_1,s_2|X)] + Cov(0, 0) = 0`. Correlating the design
does not correlate the scores. Measured model-SE coverage under that design was
0.965 and 0.950 — no clustering problem to solve. The study now uses exact
duplication, where the truth is analytic.

## Deferred

- **The dense/sparse divergence at n=100.** Paper 3 observed the dense kernel
  converging to a log-likelihood 271 units below sparse and survival on real
  meso13 data. Local reproduction failed at every scale in June, and the
  line-search fix was never tested against it because the kernel-fix check
  fits sparse only. Nothing runnable on one machine can settle this; it needs a
  MONA run on real data. `tests/testthat/test-dense-sparse-consistency.R`
  guards agreement live in the meantime.
- **A full convergence-route audit.** Verifying that a plateau-converged fit is
  genuinely at the maximum requires an *independently coded* optimiser as
  arbiter. Checking it with the same kernel is circular, and checking
  `max|grad| < tier3_grad_floor` is doubly circular, since that is the
  threshold the plateau rule already required.
