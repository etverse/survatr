---
paths:
  - "R/gcomp_survival.R"
  - "R/predict_hazard.R"
  - "R/survival_curve.R"
  - "R/variance_if_survival.R"
  - "R/variance_sandwich.R"
  - "R/rmst.R"
  - "R/estimand_registry.R"
  - "R/estimands_quantile.R"
  - "R/contrast*.R"
  - "R/print.R"
  - "R/plot.R"
  - "R/tidy.R"
  - "R/forrest.R"
  - "tests/testthat/test-sandwich-*.R"
  - "tests/testthat/test-estimands-*.R"
  - "tests/testthat/test-rmst*.R"
  - "tests/testthat/test-survival_curve.R"
  - "tests/testthat/test-contrast*.R"
---

# Hazard model, curves and estimands

These are deliberate. Do not report them as bugs without a numerical reproducer.

## Model and curve

- The outcome model is a pooled-logistic hazard on person-period rows, `logit h(t | A, L) = α(t) + β_A A + β_L L`, fitted on the risk set (rows before each id's first event and first censoring). The discrete-time hazard is the contract; do not substitute a Cox partial likelihood.
- The curve is a cumulative product per individual, `Sᵢᵃ(t) = ∏_{k≤t} (1 − hᵢₖᵃ)`, averaged afterwards. Averaging hazards before the product is biased by Jensen's inequality. A Monte Carlo over a stochastic intervention also averages at the cumulative-product level.
- The variance goes through the product by the delta method, `dSᵢᵃ(t) = −Sᵢᵃ(t) Σ_{k≤t} dhᵢₖᵃ / (1 − hᵢₖᵃ)`: causatr's outcome-model influence function gives per-row terms on the hazards, which are summed within individual. The scalar causatr engine cannot be used directly.
- The per-individual influence function is an `n × |t-grid|` matrix; the covariance over time is `crossprod(IF) / n²` and the pointwise SE its diagonal square root. Only pointwise intervals exist.
- A `mgcv::gam` fit under the sandwich uses the `predict(type = "lpmatrix")` basis from `causatr:::iv_design_matrix()`, not `model.matrix(terms(model))`, which turns `s(t)` into a linear `t` and is non-conformable with `B_inv = model$Vp`. `predict.gam()` returns a one-dimensional array; coerce it to numeric before the row scaling. Using `Vp` as the bread is mgcv's default and justified by Marra and Wood (2012).

## Estimands

- Estimand dispatch goes through the registry in `R/estimand_registry.R` (`estimand_field()`, `estimand_is_contrast()`, `estimand_is_curve()`, `estimand_has_cause()`, `estimand_requires_pair()`). A new estimand is one registry row, plus a new SE path only when its maths is new; do not reintroduce `switch(type)` lists.
- RMST uses `rmst_weights()`, whose diagonal entries are `dt[j] / 2`: the `S(0) = 1` term is the constant `dt[1] / 2`, which drops out of the delta method. Contract: `W %*% s_hat + dt[1]/2 == trapezoidal_rmst(times, s_hat)`.
- RMTL is `t − RMST` with the same quadratic-form SE, and `rmtl_difference = −rmst_difference`; equal SEs are correct.
- The quantile is indexed by `q`, not time, and accepts a single intervention. `τ_q = inf{t : S ≤ 1 − q}` by linear interpolation, with SE `IF_τ = −IF_S(τ_q) / slope`, and it has its own assembly (`assemble_quantile_result()`). Generic code that assumes a `time` index or two interventions must handle it: the bootstrap index column is `meta$index_col`, `forrest()` rejects it (`survatr_forrest_wrong_type`), and `setdiff(iv, ref)` feeding `paste0()` must guard the empty case, since `paste0()` turns a length-0 input into a phantom contrast. Aborts: `survatr_quantile_unreached`, `survatr_bad_q`.
- `risk_ratio` reports `se` on the natural scale (`RR · se(log RR)`) while the interval is built on the log scale, so `se` and the interval are deliberately not `point ± z·se`.
- Differences use symmetric Wald intervals and can have negative bounds. A survival or risk interval outside [0, 1] is a property of Wald intervals; a transformed interval would be an enhancement, not a fix.

## Methods

- `print.survatr_result` must not index with `show[seq_len(n)]`: its `n` argument collides with the `n` column of `estimates` under data.table evaluation.
- `forrest()` has two `t_ref` errors: an off-grid or malformed `t_ref` (`survatr_bad_t_ref`), and a grid time with no contrast rows (`survatr_forrest_no_contrasts`).
