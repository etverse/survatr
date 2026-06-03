# Chunk 21 — Survey / external-weight composition with IPW (transport)

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (IPW weighted-MSM core).
> **Oracle:** causatr's sampling-weighted IPW (`variance_if_ipw_sampling.R`)
> lifted to the survival functional; a survey-weighted weighted-KM cross-check.

## Goal

Allow external (survey / design) `weights` to compose with the stabilized IPW
weights so survatr can target a population other than the sample — the
target-population **transport** case. The package's supported-weights scope
already lists "survey/external (broadcast onto PP rows)"; chunk 5 rejects the
combination with `survatr_ipw_external_weights` precisely because it is a
planned-but-not-yet-built path, not because it is out of scope.

## What changes from chunk 5

- The per-id weight becomes `w_i = sw_i * f(A_i) / f(A_i | L_i)` where `sw_i`
  is the (fixed) design weight; broadcast the product onto the PP rows.
- The hazard MSM is fit with the combined weight.
- **Variance.** If the design weights are fixed (the usual case), they enter
  `phi_bar` and the score as a fixed multiplier — no extra IF block, but the
  cross-derivative and the treatment-model correction must carry the design
  weight. If the sampling fraction is itself modelled (estimated sampling
  weights / transport), add the sampling-model block, mirroring causatr's
  `variance_if_ipw_sampling.R`.
- Lift the `survatr_ipw_external_weights` guard in `surv_fit()`.

## Deliverables
- `R/ipw_survival.R` — thread the design weight into the stabilized weight, the
  broadcast, and the weight closure.
- `R/variance_if_ipw_survival.R` — design-weight-aware `phi_bar`; optional
  sampling-model block.
- Tests: survey-weighted IPW survival vs a weighted reference; sandwich vs
  two-stage bootstrap with design weights carried through the resample.

## Acceptance checklist
- [ ] `surv_fit(estimator = "ipw", weights = sw)` fits the combined-weight MSM
      (no longer aborts).
- [ ] Curve matches a survey-weighted reference on a transport DGP.
- [ ] Stacked sandwich matches the two-stage bootstrap (design weights carried).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
