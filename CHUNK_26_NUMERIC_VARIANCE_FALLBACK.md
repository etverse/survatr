# Chunk 26 — Numeric variance fallback (Tier 1 / Tier 2)

> **Status: ⬜ Not started**
> **Source:** causatr-reuse audit (2026-06-11), finding B2.
> **Depends on:** Chunk 3 (sandwich IF spine).
> **Oracle:** on a GLM hazard the numeric fallback must reproduce the analytic
> sandwich SE (chunk 3) to numerical-Jacobian tolerance (~1e-4); a custom
> `model_fn` with no analytic bread must return a finite SE that matches a
> dual-refit bootstrap within ~15%.

## Goal

Today survatr's sandwich path assumes a GLM (analytic bread via
`prepare_model_if()`) or an mgcv GAM (`Vp` bread). A hazard model fit with any
other `model_fn` — a custom family, an exotic link, a model exposing no
`estfun` / analytic bread — has **no sandwich path** at all, even though the
§4 inheritance table claims survatr "inherits" a numeric fallback via causatr's
`variance_if_numeric()`. This chunk makes that promise real: wire causatr's
Tier-1 / Tier-2 numeric variance as the fallback the cross-time delta chain
consumes when the analytic bread is unavailable.

## The reuse (do not re-derive)

causatr already ships `variance_if_numeric()`:

- **Tier 1** — `sandwich::estfun()` for the per-row score, numeric bread.
- **Tier 2** — `numDeriv` delta on the estimating function when `estfun` is
  unavailable.

survatr must **call** this, not reimplement it. The only survatr-side work is
detecting when the analytic path is unavailable and routing the resulting
per-row IF into the existing `n × |t-grid|` cross-time delta aggregation.

## Deliverables

### Updated R files
- `R/variance_if_survival.R` — when `prepare_model_if()` cannot return an
  analytic bread for `fit$model`, fall back to `causatr:::variance_if_numeric()`
  to obtain `B_inv` / `r_score`, then feed the existing delta chain unchanged.
- `R/variance_sandwich.R` — no change to the aggregation; the fallback only
  changes how `prep` is built.

### Tests (`tests/testthat/`)
- `test-variance-numeric-fallback.R` — GLM agreement vs the analytic sandwich;
  a custom `model_fn` smoke with finite SE ≈ bootstrap.

## Behaviour rules (non-negotiable)
- **Reuse `causatr:::variance_if_numeric()`** — survatr inverts no matrix of its
  own; the numeric bread is causatr's, the cross-time delta is survatr's.
- The fallback is a **fallback**: GLM / GAM keep their analytic / `Vp` breads;
  the numeric path triggers only when neither is available.

## Acceptance checklist
- [ ] GLM numeric fallback ≈ analytic sandwich (~1e-4).
- [ ] Custom `model_fn` returns finite SE ≈ bootstrap (~15%).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
