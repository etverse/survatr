---
paths:
  - "R/ice_survival.R"
  - "R/variance_if_ice_survival.R"
  - "tests/testthat/test-ice-*.R"
  - "tests/testthat/helper-ice-survival-oracle.R"
---

# ICE hazards

These are deliberate. Do not report them as bugs without a numerical reproducer.

- survatr owns the ICE backward loop and its influence-function cascade and does not call `causatr:::ice_iterate()` or `variance_if_ice_one()`, which assume a single terminal outcome and have no failure carry-forward. The survival pseudo-outcome is `Ỹₖ = Dₖ + (1 − Dₖ) qₖ₊₁`, coded `ifelse(D_k == 1, 1, q_{k+1})` so it is safe with `NA` (`ice_surv_pseudo()`).
- The `(1 − Dₖ)` factor in the cascade is required. `fit_ids[[k]]` is at risk at k and includes people with the event at k, and `dỸₖ/dβₖ₊₁ = (1 − Dₖ) dqₖ₊₁/dβₖ₊₁` must zero them, so `survatr_ice_surv_chain()` multiplies each earlier contribution by `(1 − D_{k−1})` from `build_event_by_step()`. Without it the SE is inflated, more so at longer horizons. Pinned by sandwich against bootstrap and, to 1e-5, by the `delicatessen` oracle (`test-ice-survival-delicatessen.R`, `data-raw/delicatessen_ice_survival.py`).
- The family is binomial at the horizon step (a 0/1 hazard indicator) and quasibinomial at earlier steps (a pseudo-outcome in [0, 1]); swapping them changes the score equations.
- One backward pass per horizon: the pseudo-outcome is relative to a fixed final period, so `compute_ice_survival_curve()` runs a pass per requested `t`, with `R(t) = mean(pseudo_final)` and `S(t) = 1 − R(t)`.
- Lag columns hold the observed treatment; the intervention sets the current-period treatment only (`causatr:::ice_apply_intervention_long()`).
- `confounders` are baseline and never lagged; `confounders_tv` are lag-expanded; `history` is the Markov order (`Inf` for full history). `time_formula` is not part of the per-step formula. The `causatr_fit` used for the influence function is built with `causatr:::new_causatr_fit()`, not `fit_ice()`.
- ICE standardises over the people at risk at baseline. People censored at entry carry `NA` in `pseudo_final` and drop out of the means (`na.rm`), consistent under MCAR entry censoring, and the influence-function target is the at-risk-at-baseline mask. `n_ids` passed to `fill_sandwich_ses()` stays the full first-period count, because `ice_if_setup()` already scales the sampling term by `n / n_target`. Removing `na.rm` or setting `target = rep(TRUE, n)` gives an all-`NA` curve; rescaling by `n_target` counts the factor twice.
- Clustering works with ICE: `contrast_track_b()` aligns the cluster labels with the first-period id order of the influence-function matrix. Its tests need a time-varying treatment (`sim_clustered_ice()`); a static treatment makes the lag collinear with the current value, and the resulting rank-deficiency warning comes from the data, not the clustering.
