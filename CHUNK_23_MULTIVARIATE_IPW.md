# Chunk 23 — Multivariate-treatment IPW survival

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (point IPW core), causatr multivariate-IPW internals
> (`fit_treatment_models()`, `compute_density_ratio_weights_mv()`,
> `make_weight_fn_mv()`, block-diagonal propensity bread).
> **Oracle:** causatr's multivariate-IPW tests (joint density factorization);
> `lmtp` where the joint intervention is expressible.

## Goal

Support a **joint** (multivariate) point treatment `A = (A_1, ..., A_K)` for IPW
survival. The joint conditional density is factorized by the chain rule and fit
once on the original-row (baseline) data; the **product** density-ratio weight
is broadcast onto every person-period row, and the weighted marginal hazard MSM
absorbs the per-row weight (handoff §6, "Multivariate IPW + point survival").

## What changes from chunk 5

- The single `fit_treatment_model()` becomes `fit_treatment_models()` (the
  per-component chain); the stabilized weight is the product of per-component
  density ratios via `compute_density_ratio_weights_mv()`.
- The weight closure becomes `make_weight_fn_mv()` (block-diagonal in the
  per-component `alpha` blocks).
- **Variance.** The treatment-correction's cross-derivative gains one block per
  component; the propensity bread is block-diagonal, so the chunk-5 single-block
  correction sums over components (mirrors causatr's
  `compute_ipw_if_self_contained_mv_one()`).
- `surv_fit()` accepts a length-≥2 `treatment` vector for `estimator = "ipw"`
  (currently a scalar; a vector should route here rather than error generically).

## Acceptance checklist
- [ ] Joint binary treatment IPW survival curve matches an oracle.
- [ ] Stacked sandwich (block-diagonal propensity) matches the two-stage
      bootstrap.
- [ ] Reduces to chunk 5 for a single treatment.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
