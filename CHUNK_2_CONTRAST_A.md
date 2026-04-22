# Chunk 2 -- Track A contrast path (no variance)

> **Status: ✅ done** — landed in the follow-on commit to `6e911d3`.
> 147 tests passing (1 skip on the `lmtp::lmtp_tmle` smoke oracle).
> `devtools::check()`: 0 errors / 0 warnings / 2 NOTEs (env + reserved imports).

> Per-chunk implementation guide for `survatr` chunk 2. Reads alongside
> [SURVIVAL_PACKAGE_HANDOFF.md](SURVIVAL_PACKAGE_HANDOFF.md) and
> [.claude/hard-rules.md](.claude/hard-rules.md).

## Goal

Ship the contrast path for Track A. Given a `survatr_fit` (from chunk 1) and
a list of interventions, produce per-intervention **survival curves**, then
derive risk, risk difference, risk ratio, and RMST contrasts indexed over
a user-chosen time grid.

Mechanics (non-negotiable, per `hard-rules.md`):

- Per-individual, per-period predicted hazard under intervention `a`:
  `h_hat_{i,k}^a = expit(X_{i,k}^a * beta_hat)`.
- Per-individual survival: `S_hat_i^a(t) = prod_{k <= t} (1 - h_hat_{i,k}^a)`.
- Population survival: `S_hat^a(t) = (1/n) sum_i S_hat_i^a(t)`.
- **Never average hazards before cumulating** -- Jensen's inequality on the
  nonlinear cumulative product would bias the estimate.

**No variance yet.** Sandwich + bootstrap ship in chunks 3 and 4. The
`survatr_result` returned by `contrast()` carries `se = NA_real_` columns so
downstream code (print / plot / tidy in chunk 4) reads the same shape.

## Deliverables

### New R files

| File | Contents |
|---|---|
| `R/contrast.R` | **Exported** `contrast()` S3 generic + `contrast.survatr_fit()` method. Accepts `interventions` list, `times` grid, `type` (`"survival"`, `"risk"`, `"risk_difference"`, `"risk_ratio"`, `"rmst"`, `"rmst_difference"`), `reference` for contrasts. |
| `R/apply_intervention.R` | Internal `apply_intervention_pp()` -- thin wrapper around `causatr:::apply_intervention` to operate on person-period data (broadcast baseline value to every period for point treatments, leave time-varying treatments untouched for future Track B use). Validates that the intervention yields a column of the right type. |
| `R/predict_hazard.R` | Internal `predict_hazard_pp()` -- for a given intervention, build the counterfactual PP data, call `stats::predict(fit$model, newdata, type = "response")` on **every** PP row (not just the at-risk subset), return a vector of length `nrow(pp_data)`. |
| `R/survival_curve.R` | Internal `compute_survival_curve()` -- cumulative product within id on the predicted hazards, then average across ids at each t in `times`. Returns a `data.table` with columns `intervention | time | s_hat | risk_hat | n`. |
| `R/rmst.R` | Internal `trapezoidal_rmst()` -- `RMST^a(t*) = int_0^{t*} S^a(u) du` via trapezoidal rule on the user's time grid. Also exposes `rmst_weights()` for reuse by the variance engine in chunk 3. |
| `R/contrast_types.R` | Internal dispatch for contrast `type`: assemble the `contrasts` data.table from the `estimates` data.table. Difference / ratio arithmetic per time, RMST contrast as the time integral of the difference. |

### Updated R files

| File | Change |
|---|---|
| `R/utils.R` | Add `new_survatr_result()` S3 constructor: holds `estimates`, `contrasts`, `time_grid`, `type`, `reference`, `call`. |
| `R/print.R` | Add a placeholder `print.survatr_result()` -- minimal banner + head of `contrasts`. Polished (plot / tidy / forrest) in chunk 4. |
| `R/survatr-package.R` | Add `@importFrom stats predict` if needed. |

### Tests (`tests/testthat/`)

| File | Asserts |
|---|---|
| `test-predict_hazard.R` | Predicted hazards on the NHEFS-like constant-hazard DGP under `static(0)` and `static(1)` agree with the closed-form `expit(qlogis(h) + 0)`. Length matches `nrow(pp_data)`. |
| `test-survival_curve.R` | On the constant-hazard DGP: `S_hat^a(t) -> (1 - h)^t` as n grows. Assert pointwise agreement at n = 5000, t in {1, 5, 10}, tolerance 0.01. Also assert monotone non-increasing in t. |
| `test-rmst.R` | Closed-form: with `S(t) = (1 - h)^t` on an integer grid `0:K`, `RMST(K) = sum_{t=0}^{K-1} ((1-h)^t + (1-h)^{t+1}) / 2`. Assert trapezoidal implementation matches to 1e-10. |
| `test-contrast.R` | **Happy path**: `contrast(fit, list(a1 = static(1), a0 = static(0)), times = 0:10, type = "risk_difference")` under a DGP with no treatment effect returns RD pointwise close to zero (tolerance 0.05 at n = 5000). Correct shape: `estimates` has 2 * 11 rows, `contrasts` has 11 rows with `contrast = "a1 vs a0"`. `se` columns all `NA_real_`. |
| `test-contrast-survival.R` | `type = "survival"` returns only `estimates`, empty `contrasts`. |
| `test-contrast-rejections.R` | `contrast(fit)` with empty `interventions` -> `survatr_bad_interventions`; `times` not numeric or outside observed grid -> `survatr_bad_times`; bad `type` -> `survatr_bad_contrast_type`; `reference` not in names(interventions) -> `survatr_bad_reference`. |
| `test-contrast-nhefs.R` | **Oracle smoke** against `lmtp::lmtp_tmle(outcome_type = "survival")` on a small simulated longitudinal PP dataset: point estimate of `S^a(t)` under `static(1)` and `static(0)` at t = 5 agrees within 0.05 (TMLE vs pooled-logistic gcomp differ by finite-sample bias and model-form). Marked 🟡. Skipped on CRAN. |

## API contract

```r
contrast(
  fit,                              # a survatr_fit
  interventions,                    # named list of causatr intervention objects
  times,                            # numeric vector; must be subset of fit$time_grid
  type = "risk_difference",         # "survival" | "risk" | "risk_difference" |
                                    # "risk_ratio" | "rmst" | "rmst_difference"
  reference = NULL,                 # name of reference intervention for contrasts
                                    # (default: first name in interventions)
  ci_method = "none",               # chunk 2 accepts only "none"; "sandwich" / "bootstrap"
                                    # ship in chunks 3 / 4 and are rejected here
  ...
)
# => survatr_result: list(estimates, contrasts, time_grid, type, reference, ci_method, call)
```

### Return shape

```r
result$estimates     # data.table: intervention | time | s_hat | risk_hat | se | ci_lower | ci_upper | n
result$contrasts     # data.table: contrast | time | estimate | se | ci_lower | ci_upper
result$time_grid     # numeric vector (the `times` input)
```

For `type = "survival"` and `type = "risk"`, `contrasts` is an empty
(0-row) data.table with the same columns. For `type = "rmst"`, `estimates`
has a single row per intervention (`time = t*`, `s_hat = RMST`); `contrasts`
empty. For `type = "rmst_difference"`, `estimates` has the per-intervention
RMSTs at every user time and `contrasts` has the RMST-difference at `t*`
(the max of `times`) -- or, if the user passes a scalar `times`, at that
single point.

## Behaviour rules

- **Counterfactual prediction on every PP row.** Do NOT restrict to
  at-risk rows when predicting; the cumulative product over individual `i`
  needs the hazard at every period `k` in `1..t`, whether or not `i` was
  still at risk in reality.
- **Intervention applied to the treatment column only.** Use
  `causatr:::apply_intervention()` (internal) on the copied PP data; never
  modify `fit$pp_data` in place.
- **Time grid must be a subset of the fit's observed grid.** Predictions
  at unseen times are extrapolation through `alpha(t)` -- allowed by the
  model but flagged with `rlang::warn()` (class `survatr_time_extrapolation`)
  so the user knows.
- **RMST uses trapezoidal integration on the user's `times` grid.** Users
  who want a finer grid pass more points; we don't resample.
- **Sandwich / bootstrap arguments are rejected** (`ci_method = "sandwich"`
  or `"bootstrap"`) with class `survatr_ci_not_available` pointing to the
  relevant future chunks. This is a placeholder shape -- the error message
  is informational, not "not yet implemented" phrasing.
- **IPW reweighting not handled in this chunk.** `fit$weights` is carried
  through on the result for downstream chunks but does not enter the
  population-average weighting at chunk 2 (unweighted mean across
  individuals). IPW-weighted curves ship in chunk 5.

## Non-goals (deferred)

- Sandwich variance -- chunk 3.
- Bootstrap + `plot` / `tidy` / `forrest` polish -- chunk 4.
- IPW-weighted curves -- chunk 5.
- Track B curves -- chunk 6.
- CIF under competing risks -- chunk 7.

## Seed sources

- `causatr::contrast.causatr_fit()` for the scalar contrast path and
  argument-validation patterns.
- `causatr:::apply_intervention()` for the PP-safe intervention broadcast
  (copy data, mutate the treatment column in the copy, return the copy).
- Pre-removal `causat_survival()` had no contrast implementation -- this is
  new code, but the cumulative-product math is the standard Ch. 17
  construction; cross-check against `survival::survfit` on a DGP with
  homogeneous hazard for unit-test-level sanity.

## Acceptance checklist

- [ ] `devtools::load_all()` succeeds.
- [ ] `devtools::document()` regenerates Rd for `surv_fit`, `contrast`, and
      the two S3 print methods.
- [ ] `devtools::test()` -- all tests pass (including the 79 from chunk 1).
      The NHEFS oracle smoke in `test-contrast-nhefs.R` runs locally; it is
      allowed to skip on CRAN.
- [ ] `devtools::check()` -- 0 errors, 0 warnings. NOTE about unused
      `sandwich` / `numDeriv` / `boot` imports still expected.
- [ ] `air format .` is a no-op.
- [ ] [FEATURE_COVERAGE_MATRIX.md](FEATURE_COVERAGE_MATRIX.md) updated: new
      "Contrast path" sub-section under Track A with one row per contrast
      `type` (🟢 where the closed-form or oracle pins it; 🟡 where smoke-only).
- [ ] Single commit with message `feat(chunk-2): Track A contrast -- survival
      curves + risk / RMST contrasts (no variance)`.
