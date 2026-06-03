# Feature coverage matrix

Single source of truth for **what works** in survatr. Mirrors causatr's
`FEATURE_COVERAGE_MATRIX.md` convention: every PR that adds, removes, or
changes a feature MUST update this file and the corresponding tests.

## Legend

- 🟢 — supported, truth-based test pinned against an analytical or external
  reference (`lmtp::lmtp_tmle(outcome_type = "survival")`,
  `gfoRmula::gformula_survival()`, `survival::survfit` /
  `survival::coxph`, or a closed-form DGP).
- 🟡 — smoke-tested only (runs without error; point estimate / SE not
  pinned to a truth). Acceptable for combinations where no oracle exists,
  temporary during a multi-chunk feature rollout, or where the oracle is
  too expensive to run in CI.
- 🔴 — hard-rejected by a classed error with a regression test pinning the
  rejection.

<!-- Rows are added as features ship. No speculative ⚪ entries: the matrix
reflects **current** state, not planned scope. Planned scope lives in
`SURVIVAL_PACKAGE_HANDOFF.md` §10 (implementation chunks). -->

## Track A — Point survival, pooled-logistic hazard

### Fit path (`surv_fit()`)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| `estimator = "gcomp"`, `binomial()` family, unweighted | 🟢 | `test-surv_fit-family-oracle.R` | Closed-form constant-hazard DGP: β₀ → qlogis(h) at n = 5000, K = 10, h = 0.05. Tolerance 0.1 (~1.6 × MC-SE). |
| `estimator = "gcomp"`, `quasibinomial()` family, `weights ≡ 1` | 🟢 | `test-surv_fit-family-oracle.R`, `test-surv_fit-weighted.R` | Same closed-form DGP (intercept oracle); coefficient equivalence vs unweighted `binomial()` to 1e-10. |
| Risk-set construction (drop at/after first event; drop at/after first censor) | 🟢 | `test-risk_set.R` | Hand-verified 5-id, 3-period fixture (`fixture_small_pp()`): n_at_risk = 12 without censoring, 10 with. |
| `is_uncensored()` (NA \| 0 ⇒ at-risk) | 🟢 | `test-risk_set.R` | Hand-specified truth vector. |
| Reserved-column guard (`.survatr_prev_event`, `.survatr_prev_cens`) | 🔴 | `test-checks.R`, `test-surv_fit.R` | `survatr_reserved_col`. |
| External weights validation (NA / Inf / NaN / negative / mis-sized / non-numeric) | 🔴 | `test-checks.R`, `test-surv_fit-weighted.R` | `survatr_bad_weights`. Zero weights allowed. |
| `na.action = na.exclude` | 🔴 | `test-checks.R`, `test-surv_fit.R` | `survatr_bad_na_action`. Inherited rationale from causatr (residuals padding vs `model.matrix` drop misalignment). |
| `estimator = "matching"` / `"match"` | 🔴 | `test-surv_fit.R` | `survatr_matching_rejected`. Points to `survival::coxph(..., weights = match_weights, cluster = subclass)`. |
| `estimator = "ipw"` (weighted marginal hazard MSM) | 🟢 | `test-ipw-survival.R` | See the IPW section below. |
| `estimator = "ice"` (longitudinal ICE hazards) | 🟢 | `test-ice-survival.R` | See the Track B section below. `track = "B"`, `model = NULL` (per-step models fit lazily in `contrast()`). |
| `estimator = <unknown>` | 🔴 | `test-surv_fit.R` | `survatr_bad_estimator`. |
| `competing = <non-NULL>` | 🔴 | `test-surv_fit.R` | `survatr_competing_misuse`. Cause-specific + CIF path ships in chunk 7. |
| Missing column name in `data` | 🔴 | `test-prepare_data.R` | `survatr_col_not_found`. |
| Wide input (one row per id) | 🔴 | `test-prepare_data.R` | `survatr_not_person_period`. Points to `causatr::to_person_period()`. |
| Duplicated `(id, time)` rows | 🔴 | `test-prepare_data.R` | `survatr_duplicate_pp_row`. |
| Input mutation safety (`data.frame` / `data.table` both copied) | 🟢 | `test-prepare_data.R`, `test-surv_fit.R` | Post-fit name / row-count equality with pre-fit snapshot. |

### Contrast path (`contrast.survatr_fit()`)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| `type = "survival"` on constant-hazard DGP | 🟢 | `test-contrast-survival.R`, `test-survival_curve.R` | Closed-form `S(t) = (1 - h)^t` on n = 5000, absolute tolerance 0.03 (~4 MC-SE). |
| `type = "risk"` | 🟢 | `test-contrast-survival.R` | Derived from `s_hat`; returns curve-only (empty `contrasts` stub). |
| `type = "risk_difference"` | 🟢 | `test-contrast.R` | DGP with no treatment effect: RD ≈ 0 across time, tolerance 0.02 at n = 5000. |
| `type = "risk_ratio"` | 🟢 | `test-contrast.R` | DGP with no treatment effect: RR ≈ 1, tolerance 0.15. |
| `type = "rmst"` | 🟢 | `test-rmst.R`, `test-contrast.R` | Closed-form trapezoidal integral of `(1-h)^t` matched to 1e-12; curve-only shape verified. |
| `type = "rmst_difference"` | 🟢 | `test-contrast.R` | DGP with no effect: RMST-diff ≈ 0, tolerance 0.1. |
| Oracle cross-check vs `lmtp::lmtp_tmle(outcome_type = "survival")` | 🟢 | `test-contrast-lmtp-oracle.R` | gcomp `S^a(t)` (a1 and a0) on a confounded DGP (n = 2000) matches lmtp's TMLE survival within 0.05 at t ∈ {3, 5}. lmtp 1.5.3: one fit per horizon, `ife@x` estimate (the prior `folds`/`$theta` form silently skipped). |
| Per-individual cumulative product (Jensen-safe) | 🟢 | `test-survival_curve.R` | Cumulative product within id before averaging across ids; monotone non-increasing in t on random DGPs. |
| RMST trapezoidal quadrature (closed form) | 🟢 | `test-rmst.R` | Explicit sum of `(S(t_i) + S(t_{i+1}))/2 * (t_{i+1} - t_i)` reproduced to 1e-12. |
| `fit$pp_data` mutation safety | 🟢 | `test-contrast.R` | Names / nrow / treatment column identical pre- and post- `contrast()`. |
| Empty / unnamed / duplicate-named interventions | 🔴 | `test-contrast-rejections.R` | `survatr_bad_interventions`. |
| Non-`causatr_intervention` list elements | 🔴 | `test-contrast-rejections.R` | `survatr_bad_interventions`; error message points to `causatr::static()` et al. |
| `times` not numeric / empty / NA / outside `fit$time_grid` | 🔴 | `test-contrast-rejections.R` | `survatr_bad_times` (structural) / `survatr_time_extrapolation` (values outside grid). |
| Bad `reference` (not a name in `interventions`) | 🔴 | `test-contrast-rejections.R` | `survatr_bad_reference`. |
| `ci_method = "bootstrap"` | 🟢 | `test-bootstrap-survival.R` | Percentile and Wald CIs over B replicates; per-id resampling preserves within-id dependence. See the "Bootstrap variance" section below. |
| Unknown `ci_method` / `type` / `conf_level` | 🔴 | `test-contrast-rejections.R` | `survatr_bad_ci_method` / `match.arg` / `survatr_bad_conf_level`. |
| `se` / `ci_lower` / `ci_upper` columns when `ci_method = "none"` | 🟡 | `test-contrast.R` | All `NA_real_` by design (opt-in CI path). |

### Sandwich variance (`ci_method = "sandwich"`, delta-method cross-time IF)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| Per-individual IF matrix for `S^a(t)` (Ch1 + Ch2) | 🟢 | `test-variance_if_survival.R` | Column means ~ 0 across individuals on constant-hazard DGP; sandwich SE within 75% of hand-derived closed-form target (tolerance conservative because closed form ignores at-risk attrition). Also cross-checked: sandwich SE matches empirical SD of `s_hat` across 100 simulation replicates to ~7% at t = 5, n = 3000. |
| Sandwich CI for `s_hat` / `risk_hat` | 🟢 | `test-sandwich-survival.R` | 200-rep coverage simulation at n = 1000, K = 6, h = 0.08 achieves ≥ 88% nominal 95% coverage at t ∈ {1, 5}. Single-seed sanity test at n = 5000 covers truth. |
| Sandwich CI for `risk_difference` | 🟢 | `test-sandwich-risk-difference.R` | 200-rep coverage simulation on a no-effect DGP: RD CI covers 0 at ≥ 88% nominal 95%. |
| Sandwich CI for `risk_ratio` (log-scale) | 🟢 | `test-sandwich-risk-ratio.R` | Single-seed test: CI is strictly positive, `ci_upper > ci_lower`, and covers 1 on a no-effect DGP. |
| Sandwich CI for `rmst` (trapezoidal quadratic form) | 🟢 | `test-sandwich-rmst.R` | SE is non-negative, monotone non-decreasing in `t` (by construction of the cumulative trapezoidal integral of a positive IF), CI bounds finite. |
| Sandwich CI for `rmst_difference` | 🟢 | `test-sandwich-rmst.R` | CI at `t = 10` covers 0 on a no-effect DGP (n = 3000). |
| `conf_level` in (0, 1) | 🟢 | `test-contrast-rejections.R` | Rejects values outside the open interval with `survatr_bad_conf_level`. |
| `model_fn = mgcv::gam` (penalized `s(t)` baseline) × sandwich | 🟢 | `test-sandwich-gam.R` | Counterfactual design built on the gam `lpmatrix` basis via `causatr:::iv_design_matrix()` to match the `model$Vp` bread; `predict.gam` 1-D-array output coerced to plain numeric. GAM sandwich SE matches the analytically-anchored GLM sandwich SE within 2% on a constant-hazard DGP, and tracks the bootstrap SE identically to the GLM. `Vp`-as-bread justified for frequentist coverage by Marra & Wood (2012). A gam fit lacking `$Vp` still aborts in `causatr:::bread_inv()`. |

### Bootstrap variance (`ci_method = "bootstrap"`, resample individuals)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| Per-id resampling + per-replicate refit | 🟢 | `test-bootstrap-survival.R` | Smoke: CIs populated, point estimate in `[ci_lower, ci_upper]`. Cluster = id; each replicate draws n_ids ids with replacement, concatenates their PP blocks (renumbered), refits via `surv_fit()`, contrasts via `contrast(ci_method = "none")`. |
| Percentile CI | 🟢 | `test-bootstrap-survival.R` | Single-seed coverage of `(1-h)^t` at n = 2000, B = 300. Default `boot_ci = "percentile"` (transform-invariant; safer for ratios / RMST). |
| Wald CI | 🟢 | `test-bootstrap-survival.R` | Sample-SD × `z` bands around the observed point estimate. |
| Reproducibility with `seed` | 🟢 | `test-bootstrap-survival.R` | Two calls at the same `seed` return identical SEs and CI endpoints. |
| Bootstrap SE ≈ sandwich SE (cross-check) | 🟢 | `test-bootstrap-survival.R` | Skipped on CRAN. B = 500 at n = 1500: per-time SE agrees within 15%. |
| Empirical-SD oracle: sandwich + bootstrap ≈ sampling SD | 🟢 | `test-bootstrap-survival.R` | Skipped on CRAN. 100-replicate sim truth (n = 1000, h = 0.06, K = 6). Out-of-band 300-rep validation (2026-04-22) pinned both to within 1-2% of truth; the in-test 20% tolerance catches class-of-factor-n_ids scaling bugs like the one fixed in chunk 3 (`a3f79cb`). |
| `risk_ratio` via percentile CI | 🟢 | `test-bootstrap-survival.R` | Strictly positive; covers 1 on a no-effect DGP. |
| `risk_difference` with populated contrast CIs | 🟢 | `test-bootstrap-survival.R` | `contrasts$se` / `ci_*` non-NA, point ∈ CI. |
| Failure guard (>10% replicate failures) | 🔴 | — | `survatr_boot_failed`; point the user at sandwich or at a smaller / simpler DGP. |
| Parallel backend (`parallel` + `ncpus`) | 🟢 | `test-bootstrap-rejections.R` | Validated upfront; accepts `"no"`, `"multicore"`, `"snow"`; `ncpus` positive integer. |
| Bad `n_boot` / `boot_ci` / `parallel` | 🔴 | `test-bootstrap-rejections.R` | `survatr_bad_n_boot`, `survatr_bad_boot_ci`, `survatr_bad_parallel`. |

### S3 methods on `survatr_result`

| Surface | Status | Test file | Notes |
|---|---|---|---|
| `print()` | 🟢 | `test-surv_fit.R`, `test-contrast.R` | Snapshot-pinned. Shows type, reference, ci_method, time grid, head of contrasts (or estimates for curve-only). |
| `tidy()` | 🟢 | `test-tidy-survatr_result.R` | Long `data.frame` with `intervention`, `contrast`, `time`, `estimand`, `estimate`, `se`, `ci_lower`, `ci_upper`. `which` in `{"all", "estimates", "contrasts"}`; `conf.int = FALSE` drops CI columns. S3 method on the `generics::tidy` generic (re-exported). |
| `plot()` | 🟢 | `test-plot-survatr_result.R` | Base-R graphics: curves for `survival` / `risk` / `rmst`, contrasts with reference line at 0 / 1 for the three pairwise types. CI ribbons via `adjustcolor` when populated. Smoke-only (no `vdiffr`). |
| `forrest()` | 🟢 | `test-forrest-survatr_result.R` | Forest plot at a user-chosen `t_ref`. Aborts on curve-only types (`survatr_forrest_wrong_type`), on `t_ref` outside `time_grid` (`survatr_bad_t_ref`), and on a valid grid time with no pairwise contrast rows (`survatr_forrest_no_contrasts`). |
| causatr internal-API contract (`apply_intervention`, `prepare_model_if`, `iv_design_matrix`) | 🟢 | `test-causatr-integration.R` | Pins the formal names + return shape of the three `causatr:::` internals survatr depends on, so an upstream signature drift fails loudly in CI rather than as a cryptic runtime error. |
### IPW — weighted marginal hazard MSM (`estimator = "ipw"`, binary treatment)

Fit a baseline treatment model, form **stabilized** density-ratio weights
`w_i = f(A_i) / f(A_i | L_i)` at the observed treatment, broadcast onto the
person-period rows, and fit a weighted marginal MSM `logit h(t|A) = α(t) + β_A·A`
(quasibinomial). Confounding is handled by the weights, not the hazard model.

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| IPW survival curve `S^a(t)`, binary `static` interventions | 🟢 | `test-ipw-survival.R` | `lmtp::lmtp_tmle(outcome_type = "survival")` on a confounded DGP (n = 2000): per-arm S(t) at t ∈ {3, 5} within 0.05. Mid-test agreement ~0.01. |
| Degenerate single-period IPW ⇒ scalar point IPW | 🟢 | `test-ipw-survival.R` | `causatr::causat(estimator = "ipw")` on a one-period PP table: risk under each arm matches to 1e-2 (saturated binary MSM ⇔ Hájek point IPW). |
| Unconfounded (γ = 0) ⇒ unadjusted gcomp curve; weights ≈ 1 | 🟢 | `test-ipw-survival.R` | IPW curve matches the `confounders = ~1` gcomp curve to 0.02; mean weight ≈ 1, sd < 0.25. |
| Stabilized weight constant within id (broadcast invariant) | 🟢 | `test-ipw-survival.R` | One unique weight value per id; `fit$weights` length = `nrow(data)`. |
| `trim` winsorization (fixed cutoff for the sandwich) | 🟢 | `test-ipw-survival.R` | `trim = 0.9` caps the per-id weights at the resolved quantile; untrimmed records `NA` threshold. |
| Collinear confounder (aliased propensity coef dropped) | 🟢 | `test-ipw-survival.R` | causatr drops the aliased coef (warning); IPW fit + sandwich finite. |
| IPW stacked-EE sandwich (hazard block + treatment-model correction) | 🟢 | `test-ipw-survival-sandwich.R`, `test-ipw-delicatessen.R` | Stacked SE matches the full two-stage bootstrap SE (resamples ids, refits both models) within 15%; the treatment-model correction is strictly variance-reducing vs the naive weights-as-known SE (Robins 1999; Hernán et al. 2000). |
| IPW sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-ipw-delicatessen.R` | `S^a(t)`, `S^0(t)`, and `RD(t)` point + sandwich SE match a Python `delicatessen` stacked-EE sandwich to ~1e-4 on a shared person-period fixture. Reference generated by `data-raw/delicatessen_ipw_survival.py`; both read `fixtures/python/ipw_survival_data.csv`. |
| IPW sandwich RD CI coverage | 🟢 | `test-ipw-survival-sandwich.R` | 200-rep simulation, n = 800: RD CI covers the marginal-RD truth (∫ over L) at ≥ 88% nominal 95% at t ∈ {2, 4}. |
| IPW bootstrap (per-replicate dual refit) | 🟢 | `test-ipw-survival-bootstrap.R` | Weights re-estimated each replicate (treatment + hazard refit); CIs populated, point ∈ CI, reproducible under `seed`, SE ≈ sandwich. |
| IPW × `mgcv::gam` baseline hazard | 🟢 | `test-ipw-survival-sandwich.R` | lpmatrix-basis + `Vp` bread (as gcomp gam); IPW treatment correction composes via `prep$X_fit`. Stacked sandwich tracks the two-stage bootstrap within 15% on `s(t, k=4)`. |
| Continuous / categorical / count treatment IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_treatment_unsupported` (binary only this chunk; extended types → chunk 19). |
| Time-varying treatment under IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_time_varying_treatment` (point treatment only; longitudinal IPW out of scope). |
| External `weights` + `estimator = "ipw"` | 🔴 | `test-ipw-survival.R` | `survatr_ipw_external_weights` (survey/design-weight composition = transport; deferred → chunk 21). |
| `ipsi()` under IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_ipsi_deferred` (weight-path estimand → chunk 20). |
| `trim` validation | 🔴 | `test-checks.R` | `survatr_bad_trim` (NULL or scalar in (0, 1]). |
| causatr IPW internal-API contract (`fit_treatment_model`, `evaluate_density`, `truncate_weights`, `apply_model_correction`) | 🟢 | `test-causatr-integration.R` | Pins formals + return shapes of the IPW `causatr:::` internals survatr reuses. |

## Track B — Longitudinal survival (ICE hazards)

Time-varying treatment + time-to-event via backward iterated conditional
expectations on the discrete-time hazard, with a survival-tail pseudo-outcome
`Ỹ_k = D_k + (1−D_k)q_{k+1}`. survatr owns the backward loop + cross-step IF
cascade (with the `(1−D_k)` failure carry-forward), reusing causatr's
single-model primitives. Confounders split into `confounders` (baseline) +
`confounders_tv` (lagged); `history` Markov order.

| Feature combination | Status | Test file | Oracle / notes |
|---|---|---|---|
| Track B risk curve `R^d(t)`, static strategies | 🟢 | `test-ice-survival.R` | Forward-simulation Monte-Carlo g-formula truth on a treatment-confounder-feedback DGP (n = 4000, K = 4): `S^a(t)` within 0.02 at every t; protective treatment ⇒ `S^1 > S^0`. |
| Per-step link forcing (binomial at horizon, quasibinomial earlier) | 🟢 | `test-ice-survival.R` | Asserts `models[[K]]$family == "binomial"`, `models[[k<K]] == "quasibinomial"`. |
| Lag columns carry observed (not intervened) treatment | 🟢 | `test-ice-survival.R` | Under `shift(1)`: current `A` shifted, `lag1_A` unchanged. |
| Survival-aware stacked-EE sandwich (`ci_method = "sandwich"`) | 🟢 | `test-ice-survival.R`, `test-ice-survival-delicatessen.R` | SEs finite, CIs bracket the point; `risk_difference` / `risk_ratio` / `rmst_difference` contrasts well-formed. |
| ICE sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-ice-survival-delicatessen.R` | `S^1`, `S^0`, `RD` point + sandwich SE match a Python `delicatessen` stacked-EE survival-tail sandwich to ~1e-5 at t ∈ {1,2,3} on a shared fixture. Reference: `data-raw/delicatessen_ice_survival.py`; both read `fixtures/python/ice_survival_data.csv`. Pins the `(1−D_k)` failure carry-forward. |
| ICE sandwich ≈ empirical bootstrap | 🟢 | `test-ice-survival.R` | Per-time sandwich vs bootstrap SE within 20% (n = 1200, 400 reps). Guards against regressing to causatr's verbatim chain (which over-covers, growing in t). |
| Bootstrap variance (per-replicate ICE refit) | 🟢 | `test-ice-survival.R` | `bootstrap_survival()` refits Track B per replicate (lags rebuilt, `confounders_tv` / `history` threaded). |
| External point-estimate oracle: `lmtp::lmtp_tmle(survival)` | 🟢 | `test-ice-survival-oracle.R` | Static strategies: lmtp and ICE both within 0.03 of the forward-sim truth (`skip_if_not_installed`). |
| External oracle: `gfoRmula::gformula_survival()` | 🟡 | `test-ice-survival-oracle.R` | Cross-check when installed; `skip_if_not_installed` + defensive `tryCatch` (API-sensitive). |
| `estimator = "ice"` with constant-within-id treatment | 🟢 | `test-ice-survival.R` | Informs `survatr_ice_static_treatment` (Track A cheaper) but proceeds. |
| Time-varying treatment under Track A (gcomp/ipw) | 🟢 | `test-ice-survival.R` | Warns `survatr_tv_treatment_track_a`, points to `estimator = "ice"`. |
| External `weights` + `estimator = "ice"` | 🔴 | `test-ice-survival.R` | `survatr_ice_external_weights` (weighted / IPCW longitudinal → later chunk). |
| `ipsi()` / stochastic interventions under Track B | 🔴 | `test-ice-survival.R` | `survatr_ice_intervention_deferred` (weight-path / Monte-Carlo → later chunks). |
| causatr ICE primitive contract pins | 🟢 | `test-causatr-integration.R` | Signatures of `ice_fit_step`, `ice_predict_step`, `ice_build_formula`, `ice_apply_intervention_long`, `create_lag_vars`, `ice_if_setup`, `correct_model`, `new_causatr_fit`, etc. |

## Competing risks (cause-specific hazards + CIF)

Ships in chunk 7.
