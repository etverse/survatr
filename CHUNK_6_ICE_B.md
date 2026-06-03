# Chunk 6 — Track B: longitudinal survival via ICE hazards

> **Status: ⬜ Not started**
> **Depends on:** Chunk 3 (sandwich IF), causatr ICE engine (`ice_iterate()`,
> `variance_if_ice()` / `variance_if_ice_one()`, `expand_em_lag_terms()`).
> **Oracle:** `gfoRmula::gformula_survival()` (forward-sim reference);
> `lmtp::lmtp_tmle(outcome_type = "survival")` (longitudinal point estimate).

## Goal

Time-varying treatment, time-to-event outcome via iterated conditional
expectations (Zivich et al. 2024) extended to the **hazard link**. Reuse
causatr's ICE engine wholesale; add only the survival-aware per-step target,
per-step link forcing, and the cumulative-product survival-curve shape.

## The math: ICE backward iteration on the hazard

```
Strategy d = (a_0*, a_1*, ..., a_K*):

1. Final period K — fit the all-cause hazard among rows at risk at K:
     h_K(a, L̄_K, Ā_{K-1}) = P(Y_K = 1 | at risk, A_K = a, L̄_K, Ā_{K-1})   [binomial]

2. Backward step k < K — form the pseudo-outcome as the predicted conditional
   survival TAIL through the remaining periods under the intervention:
     Ỹ_k = 1 - Ŝ^d_{k:K}
   and fit h_k targeting Ỹ_k among rows at risk at k                        [quasibinomial]

3. Iterate to k = 0.

4. Counterfactual risk:  R^d(t) = (1/n) Σ_i [1 - Ŝ^d_{0:t, i}]
   Survival:             S^d(t) = (1/n) Σ_i Ŝ^d_{0:t, i}
```

**Why this reuses causatr's ICE variance directly.** `variance_if_ice_one()`'s
forward sensitivity recursion is *agnostic to the per-step response* as long as
the per-step model exposes `family$mu.eta` and `family$variance` — both
`binomial` and `quasibinomial` do. So the block-triangular bread, the
treatment-lag EM auto-expansion, external-weight propagation, and individual
bootstrap all carry over. Track-B survival only adds the cumulative-product
survival-curve shape on top of the ICE pseudo-outcome.

## Deliverables

### New R files
- `R/ice_survival.R` — survival-aware wrapper around causatr's ICE:
  - per-step **target construction** (hazard indicator `Y_K` at K;
    survival-tail pseudo-outcome `1 - Ŝ^d_{k:K}` at `k < K`);
  - per-step **link forcing** (`binomial()` at K, `quasibinomial()` at `k < K`);
  - cumulative-product aggregation inside the contrast path.

### Updated R files
- `R/surv_fit.R` — detect time-varying treatment (treatment varies within id);
  enable `estimator = "ice"` (currently rejected); set `track = "B"`.
- `R/contrast.R` — Track-B curve assembly + the time-indexed result shape.
- `R/variance_if_survival.R` — bridge to causatr's `variance_if_ice()` for the
  per-step IFs, then the survival-tail cumulative aggregation.

### Tests (`tests/testthat/`)
- `test-ice-survival.R` — point estimate vs `gfoRmula::gformula_survival()` on a
  treatment-confounder-feedback DGP (the Ch. 21 / Table 20.1 structure).
- `test-ice-survival-oracle.R` — vs `lmtp::lmtp_tmle(outcome_type =
  "survival")` for static + shift longitudinal strategies.
- `helper-ice-survival-oracle.R` — reusable gfoRmula + lmtp comparators.

## API contract

```r
fit <- surv_fit(
  data, outcome = "Y", treatment = "A",
  confounders = ~ L1 + L2,            # time-varying covariates, lag-expanded
  id = "id", time = "t", censoring = "C",
  estimator = "ice"                   # NEW — Track B
)
# fit$track == "B"
result <- contrast(fit, interventions = list(d = causatr::static(1)),
                   times = ..., type = "risk_difference",
                   ci_method = "sandwich")  # stacked-EE ICE sandwich
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Per-step link forcing is load-bearing:** binomial at K (0/1 indicator),
  quasibinomial at `k < K` (survival-tail pseudo-outcome in [0, 1]). Swapping
  them silently changes the score equations.
- **ICE intervention applies to current-time treatment only.** Lag columns hold
  the **observed** `A_{k-1}, A_{k-2}, ...` at every backward step (inherited
  causatr rule). Do **not** recompute lags from the intervened treatment.
- **Pseudo-outcome is the survival tail** `1 - Ŝ^d_{k:K}`, not the raw hazard.
- **Reuse causatr's ICE engine** (`ice_iterate()`, `variance_if_ice()`); do not
  reimplement the backward loop or the forward sensitivity recursion.
- Cumulative-product survival shape is applied at the contrast layer (the
  chunk-2 Jensen-safe rule still holds).

## Non-goals (deferred)
- Forward-simulation g-formula (that is `gfoRmula`'s job — an oracle, not a
  feature here).
- Longitudinal IPW survival (the IPW analogue) — separate.
- Competing risks under Track B — composes with Chunk 7 after both land.

## Dependencies & composition
- causatr ICE internals via `causatr:::` — confirm `ice_iterate()`,
  `variance_if_ice()` / `variance_if_ice_one()`, `expand_em_lag_terms()` against
  the current causatr `R/` and pin them in `test-causatr-integration.R`.

## Acceptance checklist
- [ ] Track-B risk curve matches `gfoRmula::gformula_survival()` within
      tolerance on a treatment-confounder-feedback DGP.
- [ ] Static + shift strategies match `lmtp` point estimates.
- [ ] Per-step link is binomial at K, quasibinomial at `k < K` (asserted).
- [ ] Lag columns carry observed (not intervened) treatment (asserted).
- [ ] Stacked-EE sandwich populated; bootstrap agrees.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
