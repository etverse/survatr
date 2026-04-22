# survatr (development version)

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
