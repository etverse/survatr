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
| `ci_method = "sandwich"` / `"bootstrap"` | 🔴 | `test-contrast-rejections.R` | `survatr_ci_not_available`. Sandwich ships in chunk 3; bootstrap in chunk 4. |
| Unknown `ci_method` / `type` | 🔴 | `test-contrast-rejections.R` | `survatr_bad_ci_method` / `match.arg` error. |
| `se` / `ci_lower` / `ci_upper` columns | 🟡 | `test-contrast.R` | All `NA_real_` at chunk 2 by design; sandwich / bootstrap fill these in chunks 3 / 4. |

### Sandwich variance — ships in chunk 3.
### Bootstrap + S3 polish — ships in chunk 4.
### IPW weighted MSM — ships in chunk 5.

## Track B — Longitudinal survival (ICE hazards)

Ships in chunk 6.

## Competing risks (cause-specific hazards + CIF)

Ships in chunk 7.
