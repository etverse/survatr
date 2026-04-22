# survatr (development version)

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

## 2026-04-22 — Track A bootstrap + S3 polish

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

## 2026-04-22 — Track A sandwich variance

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

## 2026-04-22 — Track A contrast path (no variance)

Ship `contrast.survatr_fit()` -- the curve-shaped entry point for
Track A. Given a `survatr_fit` and a named list of `causatr` interventions,
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

## 2026-04-22 — Track A fit skeleton

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
