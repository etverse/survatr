---
paths:
  - "R/**"
  - "tests/testthat/**"
---

# Combination axes

A combination audit runs over these axes. The values are the ones the code accepts; whether a combination is supported is decided by the rejection paths in `R/checks.R`, `R/contrast_validators.R` and the fitters, and recorded per combination in `FEATURE_COVERAGE_MATRIX.md`. Re-read the argument lists when an axis may have changed.

| Axis | Values | Where |
|---|---|---|
| estimator | `gcomp` (pooled-logistic hazard, point treatment), `ipw` (weighted hazard MSM), `ice` (longitudinal ICE hazards); matching is rejected | `surv_fit(estimator =)`, `R/surv_fit.R` |
| treatment | point (constant within id) or time-varying (`ice`); numeric under `ice` | `surv_fit(treatment =)` |
| hazard model | `stats::glm` (binomial or quasibinomial), `mgcv::gam`; baseline hazard through `time_formula` | `surv_fit(model_fn =, time_formula =)` |
| intervention | causatr's constructors; `ipsi()` and stochastic interventions are rejected under `ice`, and `ipsi()` is also rejected under `ipw`; `gcomp` with either is not rejected and is an open bug in the GitHub issues | `R/contrast.R`, `R/ice_survival.R` |
| estimand (`type`) | `survival`, `risk`, `risk_difference`, `risk_ratio`, `rmst`, `rmst_difference`, `rmtl`, `rmtl_difference`, `quantile`, `cif`, `cif_difference`, `cif_ratio`, `yll` | `R/estimand_registry.R` |
| event structure | single event; competing risks through `competing =` (cause-specific hazards) | `surv_fit()` |
| censoring | row filter (`censoring =`), IPCW (`ipcw =`) | `surv_fit()` |
| weights | none, external `weights`, `trim` | `surv_fit()` |
| variance | `none`, `sandwich`, `bootstrap` (`boot_ci` `percentile` or `wald`), `cluster =` for either | `contrast(ci_method =, cluster =)` |
