---
paths:
  - "R/surv_fit.R"
  - "R/prepare_data.R"
  - "R/checks.R"
  - "R/risk_set.R"
  - "R/contrast_validators.R"
  - "tests/testthat/test-checks.R"
  - "tests/testthat/test-prepare_data.R"
  - "tests/testthat/test-risk_set.R"
  - "tests/testthat/test-surv_fit*.R"
  - "tests/testthat/test-*rejection*.R"
---

# Input checks

Each of these is a deliberate classed rejection with a regression test. Do not report it as a bug without a numerical reproducer.

- `na.action = na.exclude` is rejected (`survatr_bad_na_action`), for causatr's reason: padded residuals against dropped `model.matrix()` rows corrupt the influence function.
- `NA` in a predictor column is rejected up front (`survatr_na_in_predictors`), for IPW's treatment-model predictors too: the contrast predicts on every person-period row, and `glm`'s `na.omit` would misalign `prep$X_fit` with `fit_idx`. `NA` in the censoring column means uncensored and is accepted.
- Person-period data must be rectangular (`survatr_ragged_pp`): every id has a row at every time. Users pad with `outcome = 0` and `censoring = 1`, and the risk-set builder drops the padded rows from the fit.
- Outcome and censoring columns are 0/1 indicators (`survatr_bad_indicator`); other values would pass silently through the cumulative-sum risk-set code.
- `censoring =` alone is a row filter: the hazard model is fitted on uncensored rows and `contrast()` predicts over everyone. IPCW is the separate `ipcw =` path.
- Matching is rejected for every treatment type (`survatr_matching_rejected`); `survival::coxph(..., weights =, cluster = subclass)` on the `MatchIt` output is the tool for that.
- A constant treatment aborts under IPW (`survatr_ipw_no_treatment_variation`) before `fit_treatment_model()`, so causatr's family detection cannot read a degenerate binary column as gaussian.
- `validate_times()` accepts numeric, `Date`, `POSIXct`, `POSIXlt` and `difftime`.
- ICE scoping: external and IPCW weights (`survatr_ice_external_weights`) and `ipsi()` or stochastic interventions (`survatr_ice_intervention_deferred`) are rejected; a factor treatment is rejected (`survatr_ice_treatment_unsupported`) because the intervention writes a numeric value; a constant-within-id treatment under `ice` informs (`survatr_ice_static_treatment`) and proceeds; a time-varying treatment under `gcomp` warns (`survatr_tv_treatment_track_a`). A numeric multi-level treatment under `ice` is modelled linearly by the user's choice.
- `competing = <col>` cannot fit a cause-deleted hazard on its own: using it outside the cause-specific path raises `survatr_competing_misuse`.
