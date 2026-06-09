# snapshot: ipcw error messages are informative

    Code
      surv_fit(dt, outcome = "Y", treatment = "A", confounders = ~L, id = "id", time = "t",
        censoring = "C", estimator = "ipw", ipcw = "L")
    Condition
      Error in `surv_fit()`:
      ! `ipcw` must be a one-sided formula (e.g. `~ L1 + L2`) or `NULL`.
      i Pass the censoring-model covariate adjustment set as a formula.

---

    Code
      surv_fit(dt, outcome = "Y", treatment = "A", confounders = ~L, id = "id", time = "t",
        censoring = "C", estimator = "gcomp", ipcw = ~L)
    Condition
      Error in `surv_fit()`:
      ! Built-in IPCW (`ipcw =`) is only supported with `estimator = "ipw"` in this release.
      i IPCW for `estimator = "gcomp"` and Track B (`"ice"`) ship in later chunks. Got `estimator = "gcomp"`.

---

    Code
      surv_fit(sim_confounded_survival(n = 100L, K = 3L, seed = 1L), outcome = "Y",
      treatment = "A", confounders = ~L, id = "id", time = "t", estimator = "ipw",
      ipcw = ~L)
    Condition
      Error in `surv_fit()`:
      ! Built-in IPCW (`ipcw =`) requires a `censoring` column (the censoring indicator is the response for the censoring hazard model).
      i Supply `censoring = "<col>"` to `surv_fit()`, where `<col>` is the 0/1 column flagging censored rows.

