# Chunk 19 — IPW extended treatment types (survival)

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (binary IPW weighted-MSM core).
> **Oracle:** `lmtp::lmtp_tmle` / `lmtp::lmtp_sdr` (continuous shift), `WeightIt`
> (categorical), causatr's count-IPW tests (Poisson/NB pushforward).

## Goal

Generalize the chunk-5 binary IPW survival path to the remaining treatment
types causatr already supports at the scalar level:

- **Continuous (gaussian).** Pushforward density-ratio weights for `shift()` /
  `scale_by()` interventions (`f_d(A|L)/f(A|L)` with the Jacobian), via
  `causatr:::compute_density_ratio_weights()`'s gaussian branch. The MSM
  becomes `logit h(t | A) = alpha(t) + beta_A A` with continuous `A`.
- **Categorical (k>2).** Multinomial propensity (`nnet::multinom`); per-level
  `static()` interventions; the sandwich treatment block routes through
  `causatr:::prepare_model_if_multinom()`.
- **Count (Poisson / NB).** `propensity_family = "poisson" | "negbin"`;
  `dpois` / `dnbinom` densities; shift / scale pushforward.

## What changes from chunk 5

- `make_survival_weight_closure()` is binary-logit only. Generalize to a
  family-dispatched closure (or reuse `causatr:::make_weight_fn()` where its
  weight form matches), so `numDeriv` perturbs the right propensity parameters.
- The binary-only gate in `fit_ipw_survival()` (`survatr_ipw_treatment_unsupported`)
  is lifted for the supported families.
- The stabilized numerator generalizes per family (marginal density model).

## Deliverables
- `R/ipw_survival.R` — family-dispatched weight + closure.
- `R/variance_if_ipw_survival.R` — multinomial treatment-block prep when the
  propensity is `nnet::multinom`.
- Tests: per-family point oracle + sandwich-vs-bootstrap pin per family.

## Acceptance checklist
- [ ] Continuous `shift` / `scale_by` survival curve matches an lmtp/SDR oracle.
- [ ] Categorical per-level `static` matches a WeightIt / lmtp oracle.
- [ ] Count Poisson/NB pushforward weights match a hand / causatr reference.
- [ ] Stacked sandwich matches the two-stage bootstrap for each family.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
