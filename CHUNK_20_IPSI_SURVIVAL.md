# Chunk 20 — IPSI (incremental propensity-score) survival

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (IPW weighted-MSM core), causatr IPSI internals.
> **Oracle:** causatr's IPSI point IPW (`helper-ipw-lmtp-oracle.R`'s Kennedy
> 2019 closed form) lifted to the survival cumulative-product functional.

## Goal

Support `causatr::ipsi(delta)` for survival. IPSI (Kennedy 2019) is an
incremental propensity-score intervention: it multiplies the treatment odds by
`delta` rather than setting `A` to a fixed value. It is therefore **not** an
MSM plug-in (you cannot "set `A = a` and predict"); the counterfactual survival
is recovered by **reweighting** with the closed-form IPSI weight

```
w_i(delta) = (delta * A_i + (1 - A_i)) / (delta * p_i + (1 - p_i))
```

and cumulating the reweighted hazards. This is why chunk 5 rejects `ipsi()` on
an IPW fit with `survatr_ipw_ipsi_deferred`.

## What changes

- A dedicated `contrast()` path for `ipsi` interventions on an IPW fit that
  routes through the IPSI weight (via `causatr:::ipsi_weight_formula()` /
  `compute_density_ratio_weights()`'s ipsi branch) rather than
  `apply_intervention_pp()`.
- The estimand is `S^{delta}(t)` (and contrasts across `delta` values); the
  sandwich stacks the IPSI weight's propensity dependence (the `p_i` term).
- Lift the `survatr_ipw_ipsi_deferred` guard in `contrast.survatr_fit()`.

## Deliverables
- `R/ipw_survival.R` (or a new `R/ipsi_survival.R`) — IPSI weight-path survival
  curve + its IF.
- Tests: point estimate vs the Kennedy closed form; sandwich vs two-stage
  bootstrap over `delta`.

## Acceptance checklist
- [ ] `ipsi(delta)` on an IPW fit returns `S^{delta}(t)` (no longer aborts).
- [ ] Curve matches the Kennedy (2019) closed-form reweighting.
- [ ] Stacked sandwich matches the two-stage bootstrap.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
