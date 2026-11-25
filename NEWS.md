# survatr (development version)

No version has been released yet. What the package does today:

- `surv_fit()` fits a pooled-logistic discrete-time hazard on person-period data with three estimators: g-computation (`estimator = "gcomp"`, point treatment), a stabilized inverse-probability-weighted hazard marginal structural model (`estimator = "ipw"`, binary point treatment) and iterated conditional expectation hazards for a numeric, possibly time-varying treatment (`estimator = "ice"`). Matching is rejected.
- `contrast()` turns a fit and a named list of causatr interventions (`static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`) into survival and risk curves, risk differences and ratios, restricted mean survival time (RMST) and restricted mean time lost (RMTL) with their differences, and survival quantiles such as the median.
- Competing risks are fitted as cause-specific hazards with `surv_fit(competing = <cause column>)` under g-computation. `contrast()` then returns cumulative incidence functions, their differences and ratios, and per-cause years of life lost (`type = "yll"`).
- `surv_fit(ipcw = ~ L)` adds inverse-probability-of-censoring weights from a pooled-logistic censoring model to the IPW estimator, and the sandwich variance includes the censoring-model correction.
- `contrast(ci_method = "sandwich")` gives delta-method influence-function standard errors, `"bootstrap"` resamples individuals and refits the models (percentile or Wald intervals, optional parallel backend and seed), and `cluster = <column>` makes either variance cluster-robust.
- `diagnose()` returns per-period panels for hazard positivity, covariate balance, the IPW weight distribution, censoring incidence and the competing-risks cumulative incidence; `plot()`, `tidy()` and `forrest()` draw and tabulate a `survatr_result`.
