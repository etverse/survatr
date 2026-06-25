# survatr (development version)

## Chunk 13: cluster-robust sandwich variance

**`contrast(cluster = "<column>")` — cluster-robust standard errors.** Pass a
person-period column that labels the independent sampling unit (site, household,
provider, repeated enrolment) and the sandwich variance becomes robust to
within-cluster correlation. The per-individual influence function is already the
variance currency, so clustering is a single aggregation step: sum each
individual's IF row into its cluster's row before the cross-product, keeping the
`n^2` divisor (`V = crossprod(IF_g) / n^2`, matching `sandwich::vcovCL(type =
"HC0", cadjust = FALSE)`). The point estimate is unchanged. The reduction is the
shared `causatr:::vcov_from_if()` primitive (causatr-reuse audit), so the
per-individual and clustered paths share one tested implementation.

`cluster = "<id-column>"` reproduces the per-individual sandwich exactly
(singleton clusters), the load-bearing regression invariant. Available across
gcomp, IPW, IPCW, longitudinal ICE-hazard, competing-risks (CIF + all-cause),
and the survival quantile; the bootstrap resamples whole clusters when `cluster`
is set. Validation aborts: `survatr_cluster_varies_within_id`, `survatr_cluster_na`,
`survatr_cluster_degenerate`, `survatr_bad_cluster`.

## Chunk 12: years of life lost (YLL) + all-cause RMST/RMTL on CR fits

**`contrast(type = "yll")` — per-cause years of life lost.** On a
competing-risks fit, `YLL^(j)(t*) = int_0^{t*} F^(j)(u) du` is the integral of
the cause-`j` cumulative incidence. It reuses the chunk-7 CIF influence-function
matrix, mapped through the same trapezoidal weight matrix RMST/RMTL use (the CIF
starts at `F(0) = 0`, so the quadrature is exactly `rmst_weights() %*% cif`),
and carries the `cause` dimension and the truncation-by-death caveat of the CIF
estimands. It satisfies the identity `sum_j YLL^(j)(t*) = RMTL(t*)` (verified to
machine precision), since `sum_j F^(j) = 1 - S`. Requires a competing-risks fit:
`survatr_yll_needs_cr` on a single-event fit.

**All-cause RMST / RMTL on a competing-risks fit.** `type = "rmst"` and
`"rmtl"` now also apply to a competing-risks fit (the integral of the all-cause
survival `S = prod(1 - sum_j h^(j))`), so the `sum_j YLL = RMTL` identity is
expressible on the same fit. The per-cause RMST and the `rmst_difference` /
`rmtl_difference` contrasts on a CR fit remain out of scope.

## Chunk 12: survival quantiles / median

**`contrast(type = "quantile", q = ...)` — survival-time quantiles / median.**
The `q`-quantile is `tau_q = inf{t : S^a(t) <= 1 - q}` (median is `q = 0.5`),
found by linear interpolation between the bracketing grid points of the
survival curve chunks 2/3 already compute. The sandwich SE is the
implicit-function delta method, `d tau_q = -dS(tau_q) / S'(tau_q)`, reading the
influence function at `tau_q` by interpolating the IF columns; bootstrap is the
documented fallback when the curve is near-flat at the crossing. `q` accepts a
vector (`q = c(0.25, 0.5, 0.75)` → one row per `(intervention, q)`), and the
result is **`q`-indexed** rather than time-indexed (`estimates`:
`intervention | q | tau_hat | se | ci_*`). A median *difference* contrast is
built automatically when there are two interventions; a single intervention is
accepted (a lone median is a valid request). Available across gcomp, IPW, IPCW,
longitudinal ICE-hazard, and competing risks (all-cause survival quantile). Aborts:
`survatr_quantile_unreached` (curve never crosses `1 - q` on the grid),
`survatr_bad_q` (`q` outside `(0, 1)`).

**Central estimand registry.** Type dispatch (point column, curve vs contrast,
contrast operator, result index, SE family, plot label, cause dimension, valid
fit kinds) now lives in one descriptor table (`R/estimand_registry.R`) that the
contrast, variance, bootstrap, and S3 paths look up, replacing the scattered
`switch(type, ...)` blocks. Adding an estimand is one registry row plus a new
SE family only when the math is genuinely new.

## Chunk 12: RMTL estimand

**`contrast(type = "rmtl")` / `"rmtl_difference"` — restricted mean time lost.**
RMTL up to t* is the complement of RMST within the same window,
`RMTL(t*) = t* - RMST(t*) = ∫₀^{t*} (1 - S^a(u)) du`. It is a smooth functional
of the survival curve chunks 2/3 already compute, so it reuses the existing
`n × |t-grid|` influence-function matrix with no refitting: the point estimate
is `t* - RMST`, and the sandwich SE is the **identical** RMST trapezoidal
quadratic form (the constant restriction time has zero gradient, so
`Var(RMTL) = Var(RMST)`). Available wherever RMST is — gcomp, IPW, IPCW, and
the longitudinal ICE path — via the shared SE filler, plus bootstrap, `tidy()`, and
`plot()`. The `rmtl_difference` contrast is the negative of `rmst_difference`.

## 2026-06-09 — Critical review: diagnose + matching-rejection fixes

**`$censoring` column renamed `n_at_risk` → `n_pp_rows`.** The old name was
misleading: the denominator in the censoring panel is the full person-period
row count, not the hazard-model at-risk set. Post-event rows are included so
the censoring event at the period it occurs is always visible. `$positivity`
still uses `n_at_risk` (risk-set rows), and the docstring now clarifies the
difference. The `print.survatr_diag` output is unchanged (wording was already
"person-periods").

**Missing snapshot for `method = "match"` rejection added.** The
`test-matching-rejection.R` test for `method = "match"` was missing
`expect_snapshot(error = TRUE)`; the snapshot file now covers all 5 rejection
routes.

## 2026-06-09 — Chunk 10: survival-aware `diagnose()`

**`diagnose.survatr_fit()` returns a `survatr_diag` with five per-period diagnostic
panels.** `$positivity` summarises the predicted hazard distribution at each
time point (model-predicted for gcomp/IPW, all-cause summed for competing risks,
empirical for ICE) and flags periods with h < 0.001 or h > 0.999.
`$balance` computes per-period standardised mean differences (or
treatment–confounder correlation for continuous treatment) across all supplied
confounders. `$weights` (IPW only) reports per-id effective sample size, maximum
weight, and top-5% weight share. `$censoring` (when a censoring column is
supplied) reports per-arm censoring incidence by period. `$competing` (competing-
risks fits only) computes the marginal per-cause CIF and checks the partition-of-
unity identity `Σ F^(j)(t) + S(t) = 1` (maximum observed deviation ≈ 0);
the truncation-by-death caveat is attached as a text attribute.

`print.survatr_diag` renders a compact per-panel banner. 40 new assertions in
`test-diagnose.R`. **Completes v1 (chunks 1–10).**

## 2026-06-09 — Chunk 8: full matching rejection surface

**All matching entry routes now hard-abort with `survatr_matching_rejected`.**
Previously only `estimator = "matching"` / `"match"` was caught. Chunk 8 closes
the remaining two routes: (1) `method = "matching"` or `"match"` in `...` (a
causatr-style API mis-call, intercepted before model dispatch); (2) a `matchit`
object (MatchIt output) passed directly as `data` (caught before any column
lookup). All three routes emit the same classed error and redirect to
`survival::coxph(..., weights = match_weights, cluster = subclass)`.

New test file `test-matching-rejection.R` with 9 assertions covering all routes
and snapshot-pinning the error message for each.

## 2026-06-09 — P3 oracle hardening: gfoRmula 🟡→🟢, survRM2 RMST sanity

**gfoRmula oracle hardened (🟡→🟢, `test-ice-survival-oracle.R`).** The
`gfoRmula::gformula_survival()` cross-check for longitudinal ICE-hazard survival was
previously a loose `expect_lt(max(abs(...)), 0.05)` with `nsimul = 30 000`,
flagged 🟡 due to Monte Carlo noise. `nsimul` is raised to 100 000 (reducing MC
noise to < 0.001 per estimate) and the assertion is tightened to an absolute
bound of 0.04 (matching the known ~0.02–0.04 systematic inter-estimator offset:
gfoRmula under-estimates risk on the feedback DGP while ICE tracks the analytic
forward-sim truth to ~0.001–0.013). A directional structural pin (`ICE ≥
gfoRmula − 0.01`) is added to document the bias direction. Row promoted 🟡→🟢
in `FEATURE_COVERAGE_MATRIX.md`.

**survRM2 RMST sanity check (🟢, `test-rmst-survRM2.R`).** On an unadjusted
constant-hazard DGP (n = 3 000, K = 20, h = 0.05), survatr's per-arm
pooled-logistic trapezoidal RMST at t = 20 agrees with the KM-based RMST from
`survRM2::rmst2()` within 0.05. Validates the RMST scale without conflating
estimator class (KM vs pooled-logistic); tight precision pins remain with the
closed-form oracle in `test-rmst.R` and the `delicatessen` M-estimator in
`test-gcomp-delicatessen.R`. `survRM2` added to `Suggests`.

## 2026-06-08 — Competing risks: cause-specific hazards + cumulative incidence

`surv_fit(competing = <cause-col>)` adds first-class competing-risks survival
via parallel **cause-specific hazards + cumulative incidence functions (CIF)**.
For `J` competing event types (a single multi-valued column: `0` = no event,
`1..J` = the cause), it fits `J` pooled-logistic cause-specific hazard models on
a **shared all-cause risk set** (which is exactly "treat the other causes as
censored at their event time"). `contrast()` gains a `cause` argument and the
estimands `type = "cif"` / `"cif_difference"` / `"cif_ratio"` (per cause), plus
all-cause `"survival"` / `"risk"` from the summed hazards. The
`sum_j F^(j)(t) + S(t) = 1` identity holds numerically.

The sandwich variance stacks the `J` cause-specific hazard scores
(block-diagonal bread) and propagates them through the CIF / survival cross-time
delta. The per-cause sensitivity carries the **all-cause** `(1 - H)` denominator
(`H = sum_j h^(j)`), which does not cancel `mu_eta = h(1 - h)` the way the
single-event case does; it reduces to the chunk-3 survival IF when `J = 1`. The
point estimates match a closed-form two-cause constant-hazard CIF and the
Aalen-Johansen estimator (`survival::survfit`); the sandwich SEs match an
independent `delicatessen` stacked-EE M-estimator to ~1e-4 and the empirical
bootstrap to within ~2%. The bootstrap, `print` / `tidy` / `plot` / `forrest`,
and the contrast surface are all cause-aware. Cause-specific CIF contrasts carry
a documented truncation-by-death caveat (printed for difference / ratio results,
emitted once per session, and covered in the vignette). Fine--Gray /
subdistribution hazards are out of scope (cause-specific only). IPW / ICE
competing risks, and per-cause RMST / years-of-life-lost, ship in later chunks.

## 2026-06-03 — Longitudinal ICE-hazard: longitudinal survival (`estimator = "ice"`)

`surv_fit(estimator = "ice")` adds longitudinal causal survival estimation for
a time-varying treatment via backward iterated conditional expectations on the
discrete-time hazard, with a survival-tail pseudo-outcome
`Ỹ_k = D_k + (1 − D_k) q_{k+1}`. New longitudinal-ICE arguments split baseline
confounders (`confounders`, never lagged) from time-varying ones
(`confounders_tv`, lag-expanded) and set the Markov lag order (`history`). The
curve, contrasts (risk / RR / RD / RMST), and the stacked-EE sandwich reuse the
existing `contrast()` surface. The survival influence-function chain reuses
causatr's single-model primitives but injects a `(1 − D_k)` failure
carry-forward factor that causatr's terminal-outcome chain omits; it is
validated to ~1e-5 against an independent `delicatessen` stacked-EE M-estimator,
and the point estimates against a forward-simulation g-formula truth, `lmtp`,
and `gfoRmula`. Entry (period-1) censoring is handled by standardising over the
at-risk-at-baseline population (consistent under MCAR), rather than returning an
all-`NA` curve.

## 2026-06-02 — `risk_ratio` sandwich `se` now reported on the natural scale

The `se` column for `type = "risk_ratio"` under `ci_method = "sandwich"` was
the SE of `log(RR)` (the natural delta-method quantity), while the bootstrap
path reported the SD of the RR replicates — so the two methods printed
different `se` numbers (off by roughly a factor of `RR`) for the same
estimand. The confidence intervals were correct in both paths (the sandwich
builds the interval on the log scale and exponentiates). The sandwich now
reports `se` on the natural RR scale (`RR * se(log RR)`) so the column always
means "SE of the reported estimand" and matches the bootstrap; the CI is
still log-based (strictly positive, not `se`-symmetric). Surfaced by the IPW
work but pre-existing in the gcomp sandwich.

## 2026-06-02 — IPW weighted hazard MSM (point-treatment, `estimator = "ipw"`)

`surv_fit(estimator = "ipw")` adds inverse-probability weighting as a
second, methodologically-distinct estimator of the counterfactual survival
curve (binary point treatment; `static` / `dynamic` interventions). It fits
a baseline treatment model on the original-row data, forms **stabilized**
density-ratio weights `w_i = f(A_i) / f(A_i | L_i)` at the observed
treatment, broadcasts the per-id weight onto every person-period row, and
fits a weighted marginal hazard MSM `logit h(t | A) = α(t) + β_A·A`
(`quasibinomial`). Confounding is handled by the weights, so the hazard
model carries `A` only. The `contrast()` curve / risk / RMST machinery is
reused unchanged. New arguments: `propensity_model_fn` (treatment-model
fitter) and `trim` (optional weight winsorization at a quantile).

The sandwich variance is a two-stage stacked estimating equation: the
chunk-3 cross-time delta on the weighted hazard MSM, plus a **subtracted**
treatment-model correction that propagates the propensity-score uncertainty
through a numeric cross-derivative. For stabilized weights this *narrows*
the SE relative to treating the weights as known (Robins 1999; Hernán et
al. 2000); the stacked sandwich matches the full two-stage bootstrap (which
re-estimates both models per replicate) to within 15%, and matches an
**independent** `delicatessen` stacked-M-estimation sandwich (Python) to ~1e-4
on shared data. The point estimate is validated against
`lmtp::lmtp_tmle(outcome_type = "survival")` and, degenerately, against
`causatr::causat(estimator = "ipw")`. All IPW
weight, density, truncation, bread, and correction primitives are reused
from `causatr`; survatr composes the stabilized weight and the cross-time
delta on top.

Continuous / categorical / count treatment types and `ipsi()` are deferred
to dedicated follow-up chunks (19 / 20) and abort with clear classed errors.
A constant treatment (no positivity contrast) aborts with
`survatr_ipw_no_treatment_variation`, and external weights with
`survatr_ipw_external_weights` (transport is chunk 21). Missing data is
rejected upfront by `check_no_na_in_predictors()` for the treatment-model
predictors too — survatr does not delegate NA handling to causatr's
row-dropping, so the influence-function row alignment is preserved.

## 2026-06-02 — GAM hazard models work with sandwich variance

`contrast(..., ci_method = "sandwich")` aborted with "non-conformable
arrays" for any hazard model fit with `mgcv::gam` — advertised as a
supported `model_fn`. Two independent causes:

- The counterfactual design matrix was built with the terms-based
  `model.matrix()`, which silently degrades a smooth `s(t)` term to a
  linear `t` (fewer columns than the penalized linear-predictor basis),
  leaving it non-conformable with the `model$Vp` bread that
  `causatr:::prepare_model_if()` correctly returns for a GAM. It is now
  built via `causatr:::iv_design_matrix()`, which uses the
  `predict(type = "lpmatrix")` basis for GAMs and a level-aware,
  aliased-column-dropping `model.matrix()` for GLMs — the latter also
  hardens single-level factor interventions and rank-deficient GLM
  designs.
- `predict.gam()` returns a 1-D array (carrying a `dim` attribute) where
  `predict.glm()` returns a plain vector, so the per-row hazard
  sensitivity scaling multiplied two arrays of different shape. The link
  and response predictions are now coerced to plain numeric.

The `Vp`-as-bread + lpmatrix strategy is mgcv's own default (`vcov.gam`)
and is justified for frequentist coverage by Marra & Wood (2012, Scand.
J. Stat. 39:53–74). Validated numerically: on a constant-hazard DGP the
GAM sandwich SE matches the analytically-anchored GLM sandwich SE
pointwise to under 2%, and tracks the bootstrap SE identically to the
GLM (both ~4% below at n = 1000, R = 1000 — a generic finite-sample
sandwich-vs-bootstrap gap, not GAM-specific). New tests in
`test-sandwich-gam.R`.

## 2026-06-02 — Smaller correctness and robustness fixes

- `forrest()` raised the same `survatr_bad_t_ref` class for two distinct
  failures: a malformed / off-grid `t_ref`, and a valid grid time with no
  pairwise contrast rows (the empty `contrasts` table a single-intervention
  result produces). The latter now aborts with
  `survatr_forrest_no_contrasts` so callers can tell the two apart.
- `compute_survival_curve()` now aborts (`survatr_hazard_misaligned`) when
  the per-row hazard vector does not align 1:1 with the person-period rows,
  rather than letting `:=` silently recycle a short vector into a corrupted
  cumulative product.
- The sandwich assembly asserts each influence-function matrix is
  `n_ids x |times|` (`survatr_if_failed`) before forming
  `crossprod(IF) / n_ids^2`, so a malformed IF surfaces as a clear error
  rather than a silently wrong covariance.
- `plot.survatr_result()` builds its per-group row mask in plain R instead of
  `tbl[get(group_col) == g]`, so a column literally named `g` can no longer
  shadow the grouping variable under data.table non-standard evaluation.
- Added an integration test that pins the contract (formal names + return
  shape) of the three unexported causatr functions survatr calls via `:::`
  (`apply_intervention`, `prepare_model_if`, `iv_design_matrix`). Because the
  causatr remote is an unpinned GitHub `main`, this turns an upstream
  signature drift into a loud CI failure instead of a cryptic runtime error.

## 2026-04-22 — Round-1 critical review: 9 fixes across chunks 1–4

A full adversarial review of the 214-test chunk-1/2/3/4 surface produced
four blocking correctness bugs, three required fixes, and two
suggestions. Each reproduced numerically against a focused script under
`/tmp/survatr_repro_*.R` before a fix landed; each has a classed
regression test colocated with the affected area.

**Blocking (numerical correctness)**

- `rmst_weights()` had off-by-one indexing + double-counted the `S(0) =
  1` prefix, inflating the sandwich SE for `rmst` and `rmst_difference`
  by up to 2× at `t = 1` and ~57% on irregular grids. Point estimate
  (`rmst_hat` via `trapezoidal_rmst`) was unaffected. Fix `dd873e5`.
- `na.action = na.omit` + NA in a confounder caused `prep$X_fit`
  (post-NA-drop) to misalign with `fit_idx` (pre-drop), crashing the
  sandwich IF chain with a subscript-OOB. Rejected upfront with
  `survatr_na_in_predictors` (NA in `censoring` still allowed). Fix
  `f7c7f95`.
- Ragged PP (ids dropped post-event) crashed the sandwich / bootstrap
  paths on a defensive assertion. Rejected upfront with
  `survatr_ragged_pp` and a padding recipe. The legacy
  single-row-id = wide heuristic retired in favor of the
  rectangularity check. Fix `d6a2765`.
- `print.survatr_result` used `show[seq_len(n)]` where `show` is a
  data.table whose `n` column shadowed the function argument via NSE,
  printing hundreds of phantom NA rows. Fixed by `utils::head()`.
  Fix `de4c5e0`.

**Required (latent but fragile)**

- Bootstrap `seed` was silently non-reproducible under `parallel !=
  "no"` — `mclapply` / `parLapply` ignore the serial RNG unless
  L'Ecuyer-CMRG is set first. Now sets and restores the RNG kind
  around the parallel call. Fix `9c2b0b6`.
- Non-binary outcome / censoring values silently walked the risk-set
  cumsums. Rejected upfront with `survatr_bad_indicator`. Fix
  `bcb9d95`.
- `validate_times()` rejected `Date` / `POSIXct` / `difftime` grids.
  Relaxed to accept the time-like classes. Fix `c8e37d6`.

**Suggestions**

- `.cf_hazard` / `.cf_surv` added to `SURVATR_RESERVED_COLS`.
- `build_contrasts()` merge slimmed to select only the target
  estimand column, stable under chunk-5+ IPW columns. Fix `4935adb`.

Full suite after all fixes: pass (see commit messages for per-chunk
coverage additions).

**CI infrastructure fix** (separate from the review): `.Rbuildignore`
pattern `^\.claude$` only matched the literal path and let R CMD build
recurse into dangling `.claude/skills/` symlinks on CI. Widened to
`^\.claude(/|$)`. Fix `915d5eb`.

## 2026-04-22 — Point-treatment g-computation bootstrap + S3 polish

Ship `ci_method = "bootstrap"` in `contrast.survatr_fit()` and the S3
method surface that turns a `survatr_result` into a plot, a tidy long
table, or a forest plot.

**Bootstrap.** Resamples **individuals** (all of each id's PP rows
together), refits the hazard model via `surv_fit()` per replicate,
recomputes the per-intervention curves and the requested contrasts,
and derives pointwise CIs from the replicate distribution. Per-id
resampling is load-bearing: row-level resampling would break the
within-id cumulative-product dependence and bias variance for longer
horizons. Two CI flavors:

- **Percentile** (default): `quantile(replicates, (1-conf)/2,
  (1+conf)/2)`. Transform-invariant, so the same code handles
  `risk_ratio` without a log-scale detour.
- **Wald**: `point +/- z * sd(replicates)`. Matches the sandwich
  shape.

New arguments on `contrast()`: `n_boot` (default `500L`), `boot_ci`
(`"percentile"` | `"wald"`), `parallel` (forwarded to `boot::boot()`),
`ncpus`, `seed` (when non-null, `set.seed(seed)` before the replicate
loop so the whole sequence is reproducible). When a replicate fails
to refit / contrast it is recorded as `NA`; if > 10% of replicates
fail, the call aborts with `survatr_boot_failed`.

New rejection surface: `survatr_bad_n_boot`, `survatr_bad_boot_ci`,
`survatr_bad_parallel`.

`ci_method` now accepts all three values (`"none"`, `"sandwich"`,
`"bootstrap"`); the chunk-3 `survatr_ci_not_available` placeholder
signal is retired.

**S3 polish.**

- `tidy.survatr_result()` -- long `data.frame` with columns
  `intervention | contrast | time | estimand | estimate | se |
  ci_lower | ci_upper`. `which = c("all", "estimates", "contrasts")`,
  `conf.int = TRUE`. The `tidy` generic is re-exported from
  `generics` (now an `Imports:` dependency).
- `plot.survatr_result()` -- base-R graphics with optional CI
  ribbons. Auto-dispatches curves for `survival` / `risk` / `rmst`
  and contrasts (with reference line at 0 or 1) for the three
  pairwise types.
- `forrest.survatr_result()` -- forest plot at a user-chosen
  `t_ref`. One row per pairwise contrast. Rejects curve-only types
  with `survatr_forrest_wrong_type` and out-of-grid `t_ref` with
  `survatr_bad_t_ref`.

Bug fix on the way: `build_contrasts()` previously returned a
zero-column data.table when the user passed a contrast-shaped `type`
with only a single intervention (everyone vs everyone -> empty). It
now returns the schema-complete empty-stub used by curve-only types,
so downstream code sees a stable shape.

Testing: 30+ new tests across `test-bootstrap-survival.R`,
`test-bootstrap-rejections.R`, `test-tidy-survatr_result.R`,
`test-plot-survatr_result.R`, `test-forrest-survatr_result.R`. The
bootstrap-vs-sandwich cross-check (skipped on CRAN) pins B = 500 at
n = 1500 to within 30% of the sandwich SE.

## 2026-04-22 — Point-treatment g-computation sandwich variance

Ship `ci_method = "sandwich"` in `contrast.survatr_fit()` — delta-method
cross-time influence function aggregation on the cumulative-product
survival curve. The chunk-2 placeholder `NA_real_` columns for `se` /
`ci_lower` / `ci_upper` are now filled with pointwise Wald CIs at the
user-supplied `conf_level` (default `0.95`).

Per-individual IF on `S^a(t)` decomposes into two pieces:

1. **Ch1** (target variability): `IF1_i(t) = S^a_i(t) - S^a(t)`.
2. **Ch2** (beta uncertainty): `IF2_i(t) = -n_ids * J_bar(t)' * B_inv * psi_i`,
   where `J_bar(t)` is the population gradient of `S^a(t)` wrt beta,
   `B_inv = (X'WX)^{-1}` from `causatr:::prepare_model_if()`, and `psi_i`
   is the per-individual sum of score contributions across the
   at-risk rows belonging to individual `i`.

The cross-time covariance `V = crossprod(IF_mat) / n_ids^2` yields
pointwise SEs on the diagonal; the RMST SE is the quadratic form
`w_j' V w_j` with trapezoidal weights from `rmst_weights()`. Risk
difference / RMST difference propagate via the contrast IF
`IF_ref - IF_a1`. Risk ratio is built on the log scale and the CI
endpoints are exponentiated so the reported bounds are strictly positive.

Wiring: the `survatr_ci_not_available` signal now triggers only on
`ci_method = "bootstrap"` (reserved for chunk 4). The `conf_level`
argument is validated at the boundary with the new
`survatr_bad_conf_level` class.

New rejection surface: `survatr_bad_conf_level`.

Testing: 26 new tests across `test-variance_if_survival.R`,
`test-sandwich-survival.R`, `test-sandwich-risk-difference.R`,
`test-sandwich-rmst.R`, `test-sandwich-risk-ratio.R`. Coverage pinned by
two 200-rep simulations (nominal 95% on `s_hat` and on `risk_difference`
around 0 under a no-effect DGP) achieving ≥ 88% coverage — skipped on
CRAN to keep the full-check runtime reasonable. A companion single-seed
check confirms the sandwich SE matches the empirical SD across 100
replicates at n = 3000 within 7% at t = 5.

Full suite: 173 passing / 0 failing / 1 skip.
`devtools::check()`: 0 errors / 0 warnings / 2 NOTEs (env timestamp +
unused `Imports: boot` reserved for chunk 4).

## 2026-04-22 — Point-treatment g-computation contrast path (no variance)

Ship `contrast.survatr_fit()` -- the curve-shaped entry point for
point-treatment g-computation. Given a `survatr_fit` and a named list of `causatr` interventions,
build the counterfactual person-period data under each intervention,
predict per-row hazards, cumulate within individual to `S^a_i(t)`, and
average across individuals to `S^a(t) = (1/n) sum_i S^a_i(t)`. Derive
risk / risk-difference / risk-ratio / RMST / RMST-difference contrasts
on a user-chosen time grid.

Returns a `survatr_result` S3: `$estimates` and `$contrasts` data.tables
indexed by `time`, plus `time_grid`, `type`, `reference`, `ci_method`,
and `call`. `se` / `ci_lower` / `ci_upper` columns are `NA_real_` at
this chunk; sandwich variance ships in chunk 3 and bootstrap in
chunk 4. `ci_method = "sandwich"` / `"bootstrap"` are rejected now with
class `survatr_ci_not_available`.

RMST is computed with a trapezoidal quadrature on the user's `times`
grid, prepending `S(0) = 1` to integrate from the origin; the
quadrature weights are exposed as `rmst_weights()` for chunk 3's
delta-method RMST SE.

`contrast()` is an S3 generic local to survatr (not a re-export of
`causatr::contrast`, which is a plain function there). When both
packages are attached, survatr's generic shadows the causatr function;
users who need the scalar-outcome path should call
`causatr::contrast()` explicitly.

Testing: 68 new tests across `test-contrast.R`, `test-contrast-survival.R`,
`test-contrast-rejections.R`, `test-predict_hazard.R`, `test-survival_curve.R`,
`test-rmst.R`, and a 🟡 `test-contrast-lmtp-oracle.R` that skips on CRAN
and on toy-DGP fit failures. Closed-form oracle: `S(t) = (1 - h)^t`
agrees with `contrast(type = "survival")` within 0.03 absolute at
n = 5000; trapezoidal RMST matches hand-expanded sum to 1e-12.

## 2026-04-22 — Point-treatment g-computation fit skeleton

First working entry point: `surv_fit()` fits the pooled-logistic
discrete-time hazard model `logit h(t | A, L) = alpha(t) + beta_A A + beta_L L`
on person-period data and returns a `survatr_fit` S3 object. The fit builds
the risk set by dropping rows at/after the first event per id and rows
at/after the first censor per id, and automatically switches the GLM family
to `quasibinomial()` when external weights are supplied.

Ports the pre-removal `causatr::causat_survival()` fit path into a dedicated
package with the naming / error-class surface tightened: reserved internal
columns are `.survatr_prev_event` / `.survatr_prev_cens`, classed errors use
the `survatr_*` prefix, and `model_fn` is a first-class argument (default
`stats::glm`, accepts any `mgcv::gam`-style fitter).

Rejection paths wired with classed errors: `survatr_matching_rejected`
(points to `survival::coxph(..., weights = , cluster = )`),
`survatr_bad_estimator`, `survatr_competing_misuse`, `survatr_bad_na_action`,
`survatr_reserved_col`, `survatr_bad_weights`, `survatr_col_not_found`,
`survatr_not_person_period`, `survatr_duplicate_pp_row`.

Testing: 79 tests across `test-checks.R`, `test-prepare_data.R`,
`test-risk_set.R`, `test-surv_fit.R`, `test-surv_fit-weighted.R`, and
`test-surv_fit-family-oracle.R`. The oracle test pins the logit intercept
against `qlogis(h)` on a closed-form constant-hazard DGP (n = 5000) for
both `binomial()` and `quasibinomial()` families.

Contrast, variance, bootstrap, IPW, ICE, competing risks, and
`diagnose()` ship in subsequent chunks.
