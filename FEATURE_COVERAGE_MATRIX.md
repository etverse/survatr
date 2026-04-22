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
| `estimator ∈ {"ipw", "ice", <unknown>}` | 🔴 | `test-surv_fit.R` | `survatr_bad_estimator`. Tracks B / IPW ship in later chunks. |
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
| Oracle cross-check vs `lmtp::lmtp_tmle(outcome_type = "survival")` | 🟡 | `test-contrast-lmtp-oracle.R` | Point-estimate smoke test skipped on CRAN and skipped when the lmtp call can't fit the toy DGP. Kept to catch regressions on larger DGPs. |
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
| `model_fn` ≠ `stats::glm` (e.g. `mgcv::gam`) | 🔴 | — | Sandwich code uses `causatr:::prepare_model_if()` which abort-early on `mgcv::gam` without `$Vp`. Bootstrap (below) is the user path for non-GLM fitters. |

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
| `forrest()` | 🟢 | `test-forrest-survatr_result.R` | Forest plot at a user-chosen `t_ref`. Aborts on curve-only types (`survatr_forrest_wrong_type`) and on `t_ref` outside `time_grid` (`survatr_bad_t_ref`). |
### Bootstrap + S3 polish — ships in chunk 4.
### IPW weighted MSM — ships in chunk 5.

## Track B — Longitudinal survival (ICE hazards)

Ships in chunk 6.

## Competing risks (cause-specific hazards + CIF)

Ships in chunk 7.
