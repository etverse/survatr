# survatr (development version)

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
