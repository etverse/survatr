# Chunk 3 -- Track A sandwich variance (delta-method cross-time IF)

> **Status: ✅ done** -- landed in the chunk-3 commit.
> 173 tests passing / 1 skip. `devtools::check()`: 0 errors / 0 warnings /
> 2 NOTEs (env + unused `boot` reserved for chunk 4).

> Per-chunk implementation guide for `survatr` chunk 3. Reads alongside
> [SURVIVAL_PACKAGE_HANDOFF.md](SURVIVAL_PACKAGE_HANDOFF.md) and
> [.claude/hard-rules.md](.claude/hard-rules.md).

## Goal

Turn the `NA_real_` `se` / `ci_lower` / `ci_upper` columns from chunk 2 into
real values. Implement sandwich variance for the survival curve via the
**delta-method cross-time IF chain**. Per-individual IF becomes an
`n x |t-grid|` matrix; full cross-time covariance is `crossprod(IF) / n^2`;
pointwise SE is the diagonal; RMST SE is the quadratic form
`w' V w` with trapezoidal weights `w` from `rmst_weights()`.

## The math (for reference while coding)

Per-individual, per-period hazard under intervention a:

    h_hat_{i,k}^a = expit(X_{i,k}^a' beta_hat)

Per-individual survival at t:

    S_hat_i^a(t) = prod_{k <= t} (1 - h_hat_{i,k}^a)

Delta-method derivative along the cumulative product:

    d log S_hat_i^a(t) / d beta
      = sum_{k <= t} d log(1 - h_hat_{i,k}^a) / d beta
      = - sum_{k <= t} [mu_eta(eta_{i,k}^a) / (1 - h_hat_{i,k}^a)] X_{i,k}^a

    d S_hat_i^a(t) / d beta
      = - S_hat_i^a(t) * sum_{k <= t} [mu_eta(eta_{i,k}^a) / (1 - h_hat_{i,k}^a)] X_{i,k}^a

Call this `J_i^a(t) := d S_hat_i^a(t) / d beta` (a length-p row vector
per individual per time).

The population survival is `S_hat^a(t) = (1/n) sum_i S_hat_i^a(t)`. The
first-order IF for `S_hat^a(t)` has two pieces:

1. **Ch1** (target variability, mean-zero):
   `IF1_i^a(t) = S_hat_i^a(t) - S_hat^a(t)`.
2. **Ch2** (estimation of beta): `(1/n) sum_i J_i^a(t)` times the IF on
   beta_hat from `causatr:::prepare_model_if()`. In stacked form,
   `IF2_i^a(t) = J_bar^a(t) * A_inv * psi_i` where `psi_i` is the score
   contribution for the i'th at-risk row and `A_inv` is the bread.

Full per-individual IF on `S_hat^a(t)`:

    IF_i^a(t) = (S_hat_i^a(t) - S_hat^a(t)) + J_bar^a(t) * A_inv * psi_i

(For rows outside the fit_rows mask, the Ch2 contribution is zero.)

Cross-time covariance of `S_hat^a(.)`:

    V^a = crossprod(IF_mat^a) / n^2     ## |t| x |t|

where `IF_mat^a` is `n x |t|` with `IF_mat^a[i, j] = IF_i^a(t_j)`.

Risk: `risk_hat^a(t) = 1 - S_hat^a(t)`, IF is `-IF^a(t)` -- same variance.

Risk difference: `RD(t) = risk_hat^{a1}(t) - risk_hat^{a0}(t)`, IF is
`-(IF^{a1}(t) - IF^{a0}(t))`. Variance by the same cross-product.

Risk ratio: `RR(t) = risk^{a1}(t) / risk^{a0}(t)`. Delta method on the
ratio:

    IF_RR(t) = (1 / risk^{a0}(t)) * IF_risk^{a1}(t)
             - (risk^{a1}(t) / risk^{a0}(t)^2) * IF_risk^{a0}(t)

SE by the cross-product. We report **log-RR** SE (and exponentiate the
CI endpoints) by default -- log-scale CIs are well-behaved and symmetric.

RMST: `RMST^a(t_j) = sum_{i=1..j} w_j[i] * S_hat^a(t_i)` (trapezoidal
weights from `rmst_weights()`). IF: `w_j %*% IF_mat^a`. Variance:
`w_j' V^a w_j`. RMST difference: `w_j' (V^{a1} + V^{a0} - 2 C^{a1,a0}) w_j`
where `C` is the cross-intervention covariance of the IF matrices.

## Deliverables

### New R files

| File | Contents |
|---|---|
| `R/variance_if_survival.R` | **Internal** `compute_survival_if()` -- per-intervention, returns the `n x |t|` IF matrix for `S_hat^a(t)` and the mean `J_bar^a(t)` matrix (p x |t|). Uses `causatr:::prepare_model_if()` for the Ch2 pieces and assembles the delta chain from per-row `(mu_eta / (1 - h_hat), X)` products. |
| `R/variance_sandwich.R` | **Internal** `sandwich_result()` -- takes the list of per-intervention IF matrices and returns a filled-in `survatr_result` with `se`, `ci_lower`, `ci_upper`. Handles risk / risk_difference / risk_ratio / rmst / rmst_difference via the formulas above. |

### Updated R files

| File | Change |
|---|---|
| `R/contrast.R` | Accept `ci_method = "sandwich"`. When set, compute the IF matrices after `compute_survival_curve()` and hand off to `sandwich_result()`. Drop the `survatr_ci_not_available` branch for `"sandwich"`. |
| `R/rmst.R` | Already exposes `rmst_weights()`. No change expected but confirm the quadrature matrix lines up with what the sandwich code expects (row j = weights for `RMST(t_j)`). |
| `R/survatr-package.R` | Add `@importFrom numDeriv jacobian` and `@importFrom sandwich bread estfun` if they are used in the Tier-1 numeric fallback path. |

### Tests (`tests/testthat/`)

| File | Asserts |
|---|---|
| `test-variance_if_survival.R` | IF matrix for `S_hat^a(t)` on a constant-hazard DGP has mean 0 at every t (across the n rows), and `sqrt(crossprod(IF) / n^2)` diagonal matches the analytical asymptotic SE within 10%. |
| `test-sandwich-survival.R` | End-to-end `contrast(..., ci_method = "sandwich")` on the constant-hazard DGP: `s_hat` CI covers `(1 - h)^t` at nominal 95% across ~500 simulation replicates (marked 🟢). Single-seed test + coverage simulation. |
| `test-sandwich-vs-bootstrap.R` | **Cross-check** against a small bootstrap (B = 200) on the same fit: per-time SEs agree within 25%. Bootstrap code lives in chunk 4 but can be stood up as a local helper in this test file temporarily -- or we delay this assertion to chunk 4. **Recommendation:** include the helper here so chunk 3's sandwich is not trusted on its own. |
| `test-sandwich-rmst.R` | RMST SE via quadratic form on `V^a` matches a direct delta-method computation on a closed-form constant-hazard DGP (tolerance 1e-2). |
| `test-sandwich-risk-ratio.R` | log-RR variance on a DGP with known relative risk: CI coverage at nominal 95% in a small simulation (n = 2000, B_sim = 200). |
| `test-nhefs-sandwich.R` | **🟡 smoke**: fit on NHEFS Ch. 17 data, compute sandwich CI for `risk_difference` at t = 120; assert the CI overlaps the book's interval (-4.1%, 3.7%). Skipped on CRAN. Full replication test ships in chunk 9. |

## API contract

```r
result <- contrast(
  fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times = seq(0, 120, by = 12),
  type = "risk_difference",
  ci_method = "sandwich",
  conf_level = 0.95
)
# result$estimates:
#   intervention | time | s_hat | risk_hat | se | ci_lower | ci_upper | n
# result$contrasts:
#   contrast | time | estimate | se | ci_lower | ci_upper
```

- `conf_level` defaults to `0.95`; accepts any value in `(0, 1)`.
- `se` is the pointwise SE on `s_hat` / `risk_hat` / `rmst_hat` /
  contrast `estimate`.
- `ci_lower` / `ci_upper` are Wald CIs on the native scale for
  `s_hat` / `risk_hat` / `rmst_hat` / `risk_difference` / `rmst_difference`.
  For `risk_ratio`, CIs are computed on the log scale and exponentiated.

## Behaviour rules

- **Use causatr's IF primitives, do not reimplement.**
  `causatr:::prepare_model_if()` returns `B_inv` and `r_score`;
  `apply_model_correction(prep, gradient)` wires them through a gradient
  vector. For the cross-time chain we call `apply_model_correction()`
  once per time point with `gradient = J_bar^a(t)`, or (equivalently and
  more efficiently) compute `J %*% A_inv` once and apply to `psi` via
  a matrix product.
- **`J_bar^a(t)` is averaged across the target population at time t.**
  In Track A that's the full n (Ch1 uses `S_i - S_bar` for mean-zero
  variability). No IPW reweighting in this chunk.
- **Pointwise bands only in v1.** Simultaneous bands (Hall-Wellner,
  equal-precision) are out of scope.
- **Quasibinomial dispersion.** When `fit$family == "quasibinomial"`,
  the variance propagates through `causatr:::prepare_model_if()` which
  already accounts for dispersion via `stats::weights(model, "working")`.
  Confirm with a small DGP where we fit both families and the IPW
  weights are ≡ 1 -- SEs should differ by ≤ the expected dispersion
  factor.
- **Extrapolation outside `fit$time_grid` is still rejected** at chunk 3.
- **Numeric fallback (Tier 1 / Tier 2) is NOT wired this chunk.** If the
  analytical delta chain fails (e.g. non-GLM `model_fn`), raise
  `survatr_if_failed` with a pointer. A numeric fallback via
  `numDeriv::jacobian()` on the per-time `s_hat` can ship in a follow-up.

## Non-goals (deferred)

- Bootstrap variance -- chunk 4.
- Simultaneous confidence bands -- out of v1 scope.
- IPW-weighted curves -- chunk 5 (will extend the delta chain to
  weighted averages and density-ratio-weight IFs).
- Track B ICE variance -- chunk 6 (reuses `causatr:::variance_if_ice()`).
- CR variance (stacked EE across cause-specific hazards) -- chunk 7.

## Seed sources

- `causatr:::prepare_model_if()` -- returns the bread and score-adjusted
  residual vector. Use verbatim.
- `causatr:::apply_model_correction()` -- applies a single-gradient IF
  correction. We will likely extend this inline (the cross-time chain
  needs the same correction at each of `|t-grid|` gradient vectors; do
  the matrix product by hand rather than looping).
- `causatr:::vcov_from_if()` -- call shape: `crossprod(IF_list) / n^2`.
  Reused for the final covariance.
- `causatr:::variance_if_numeric()` -- Tier 1 / Tier 2 numeric fallback
  that uses `numDeriv::jacobian()`. Useful as a reference implementation
  if we ever wire a numeric fallback path.

## Acceptance checklist

- [ ] `devtools::load_all()` succeeds.
- [ ] `devtools::document()` regenerates Rd cleanly.
- [ ] `devtools::test()` -- all tests pass (chunk 1 + 2 unchanged, new
      sandwich tests green). Coverage simulation (`test-sandwich-survival.R`)
      runs locally; allowed to skip on CRAN.
- [ ] `devtools::check()` -- 0 errors, 0 warnings. Unused-import NOTE
      should drop to `boot` only (numDeriv + sandwich now used).
- [ ] `air format .` is a no-op.
- [ ] [FEATURE_COVERAGE_MATRIX.md](FEATURE_COVERAGE_MATRIX.md) updated:
      replace the chunk-2 🟡 placeholder rows for `se` / `ci_*` with 🟢
      rows pinned to the coverage simulation / sandwich-vs-bootstrap tests.
- [ ] Single commit with message `feat(chunk-3): Track A sandwich variance
      -- delta-method cross-time IF aggregation`.
