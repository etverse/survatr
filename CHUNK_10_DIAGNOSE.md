# Chunk 10 — Survival-aware `diagnose()`

> **Status: ⬜ Not started**
> **Depends on:** Chunk 2 (fit + curve), Chunk 7 (competing-risks decomposition);
> Chunk 5 (IPW weight diagnostics).
> **Oracle:** n/a (diagnostic summaries; assert structure + known-DGP signals).

## Goal

A survival-aware `diagnose()` returning a `survatr_diag` S3 object with
**per-period** panels — the survival analogue of causatr's `diagnose()`
(positivity, balance, weight summaries), extended across the time axis and to
competing risks.

## Panels

1. **Per-period hazard positivity.** For each period `k`, summarize the
   predicted counterfactual hazards `ĥ^a(k)` (and, under IPW, the per-id weight
   distribution + effective sample size). Flag periods where hazards approach 0
   or 1, or where late-time risk sets are sparse (the spline-baseline danger
   zone noted in the GAM/bread discussion).
2. **Cross-time covariate balance.** Balance of confounders by treatment arm at
   each time step (standardized mean differences), surfacing time-varying
   imbalance for Track B.
3. **Weight diagnostics (IPW).** Per-period (Track B) / per-id (Track A)
   density-ratio weight distribution, ESS, max weight, share of weight in the
   top quantile (the truncation decision input).
4. **Censoring-by-arm.** Censoring incidence per period per treatment arm — the
   informative-censoring smell test (motivates IPCW).
5. **Competing-risks decomposition.** Per-cause CIF summaries + the
   `Σ_j F^(j)(t) + S(t) = 1` identity check, with the **truncation-by-death
   interpretational caveat** surfaced explicitly (carried from Chunk 7).

## Deliverables

### New R files
- `R/diagnose.R` — `diagnose.survatr_fit()` building the panels above; the
  `survatr_diag` constructor.
- `R/print.R` — `print.survatr_diag` (compact per-period summary).

### Updated R files
- `R/survatr-package.R` — `@importFrom` as needed.

### Tests (`tests/testthat/`)
- `test-diagnose.R` — structure of `survatr_diag`; on a known-DGP, positivity
  flags fire when hazards are pushed near 1; balance panel reports ~0 SMD on a
  randomized DGP and non-zero on a confounded one; ESS sane under IPW; the CR
  identity check holds.

## API contract

```r
fit <- surv_fit(...)
dx <- diagnose(fit)            # survatr_diag
dx$positivity   # per-period hazard / weight summary (data.table)
dx$balance      # per-period SMD by arm
dx$weights      # ESS, max, top-quantile share (IPW only)
dx$censoring    # per-period censoring by arm
dx$competing    # per-cause CIF decomposition + identity check (CR only)
print(dx)
```

## Behaviour rules (non-negotiable — see hard-rules.md)
- Returns a `survatr_diag` S3; panels are `data.table`s indexed by `time` (and
  `cause` for the CR panel).
- Diagnostics are **non-fatal**: report and flag, never abort a fit.
- The truncation-by-death caveat is text on the CR panel, not a silent number.

## Non-goals (deferred)
- Formal positivity hypothesis tests (report distributions + flags only).
- Plot methods for `survatr_diag` (could follow; not required for v1).

## Dependencies & composition
- Chunk 2 (curve), Chunk 5 (weights), Chunk 7 (CR decomposition).

## Acceptance checklist
- [ ] `diagnose()` returns `survatr_diag` with the five panels (those
      applicable to the fit's estimator/track).
- [ ] Positivity flags fire on a near-1-hazard DGP; balance ~0 on randomized.
- [ ] CR identity check + truncation-by-death caveat present.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
