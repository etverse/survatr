---
paths:
  - "R/ipcw_survival.R"
  - "tests/testthat/test-ipcw-*.R"
  - "tests/testthat/helper-ipcw-survival-oracle.R"
---

# Censoring weights

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The censoring model is fitted on rows with `prev_event == 0 & prev_cens == 0`, which includes the row where `C_k = 1`, because that is where censoring is observed. This is broader than the hazard-MSM fit rows, which also need `C_k == 0` (`is_uncensored()`).
- The per-person censoring influence function is `(n_ids / n_cens_fit) * psi_per_id %*% B_inv_cens`: the model is fitted on `n_cens_fit` person-period rows, `causatr:::prepare_model_if(..., n_total = n_cens_fit)` returns `B_inv = n_cens_fit (X'WX)⁻¹`, and the factor converts to the person level. It matches the `delicatessen` oracle within 2% (`test-ipcw-delicatessen.R`, `data-raw/delicatessen_ipcw_survival.py`).
- The IPCW running product uses all person-period rows (`compute_ipcw_running_weights()`, `ipcw_running_cumprod()`); only the MSM score is restricted to at-risk rows. The oracle uses the same convention.
- The `phi_bar_cens` closure for `numDeriv` uses trimming thresholds fixed at the point estimate (`fit$ipcw_trim_thresholds`); recomputing quantiles inside it would make the weights non-smooth in γ.
- The censoring correction can widen or narrow the SE, depending on its covariance with the rest of the influence function; it has no guaranteed sign.
- The bootstrap refits both the treatment and the censoring model in each replicate (`fit$ipcw`, `fit$censoring_model_fn`); without either refit, bootstrap and sandwich disagree by about 15%.
