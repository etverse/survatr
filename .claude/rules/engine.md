---
paths:
  - "R/**"
---

# How survatr builds on causatr

- causatr is the engine. Influence-function primitives (`prepare_model_if()`, `apply_model_correction()`, `vcov_from_if()`, the numeric fallback), the intervention constructors and `to_person_period()` stay in causatr; survatr adds the cross-time aggregation and the curve-shaped estimands. Today survatr reaches these through `causatr:::`; the design replaces that with causatr's exported `dev_` functions, so do not reimplement a primitive in survatr to avoid the `:::` call.
- `test-causatr-integration.R` pins the signatures survatr relies on: `apply_intervention(data, treatment, iv)`, `prepare_model_if(model, fit_idx, n_total)` returning `X_fit`, `B_inv` and `r_score`, and `iv_design_matrix(model, newdata)` returning a design whose columns match the model coefficients.
- survatr inverts no matrix of its own: every bread inverse comes from causatr's `prepare_model_if()` (guarded with `MASS::ginv()`, `$Vp` for GAMs). Its variance code is matrix products plus a guarded risk-ratio division, so a singular-bread guard in survatr has nothing to guard.
- `estimator`, not `method`: `method` is reserved for `...` forwarding.
- The fit holds the hazard models; curves exist only in the `contrast()` result. A `survatr_result` has `estimates` (`intervention | time | s_hat | se_s | …`), `contrasts` (`contrast | time | estimate | se | …`) and `time_grid`, and does not fit the scalar `causatr_result` shape.
- Internal person-period columns use the reserved prefix `.survatr_` (`.survatr_prev_event`, `.survatr_prev_cens`, `.survatr_any_event`), plus `.cf_hazard` and `.cf_surv` on internal copies; all are guarded against collisions with user columns at the boundary.
