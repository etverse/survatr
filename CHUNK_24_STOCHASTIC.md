# Chunk 24 — Stochastic interventions + survival

> **Status: ⬜ Not started**
> **Depends on:** Chunk 2 (contrast spine), Chunk 3 (sandwich), causatr
> stochastic-intervention internals (causatr `PHASE_12_STOCHASTIC.md`).
> **Oracle:** `lmtp` (modified treatment policies / stochastic shifts);
> causatr's stochastic-intervention tests.

## Goal

Support **stochastic** interventions (modified treatment policies / random
shifts) for survival curves. This is the one intervention class that does NOT
compose away trivially, because of the cumulative-product nonlinearity.

## The survival-specific rule (handoff §6, non-negotiable)

Monte-Carlo draws are taken **at the individual cumulative-product level**, not
at the hazard level:

```
S^g_i(t) = (1/M) sum_m prod_{k <= t} (1 - h(k | A_{i,m}, L_{i,k}))
```

Averaging hazards across draws **before** cumulating is biased (Jensen's
inequality — the cumulative product is nonlinear). The estimand is the average
over draws of the per-draw survival curve.

## What changes

- `contrast()` accepts a stochastic intervention; for each MC draw it applies
  the drawn treatment, predicts hazards on all PP rows, cumulates per id, and
  averages the resulting **survival curves** over draws (not hazards).
- Variance: the IF must account for the MC averaging (or use bootstrap); the
  sandwich extends the chunk-3 delta chain over the draw-averaged curve.

## Acceptance checklist
- [ ] Stochastic-intervention survival curve matches an `lmtp` MTP oracle.
- [ ] MC averaging is at the cumulative-product level (regression test against
      the biased hazard-level average to lock the rule in).
- [ ] Sandwich / bootstrap variance for the stochastic estimand.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
