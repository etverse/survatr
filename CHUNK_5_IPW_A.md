# Chunk 5 — Track A under IPW (weighted hazard MSM)

> **Status: ⬜ Not started**
> **Depends on:** Chunk 2 (contrast spine), Chunk 3 (sandwich IF), causatr IPW
> internals.
> **Oracle:** `lmtp::lmtp_tmle(outcome_type = "survival")` (point estimate);
> Cole & Hernán (2004) adjusted survival curves; `causatr::causat(estimator =
> "ipw")` on a single-period person-period table (degenerate-to-scalar check).

## Goal

Estimate `S^a(t)` for a **point** treatment by inverse-probability weighting
instead of outcome standardization. Fit a baseline treatment model on the
**original-row** data (one row per id), form stabilized density-ratio weights,
**broadcast** the per-id weight onto every person-period row, and fit a
**weighted pooled-logistic hazard MSM** `logit h(t | A) = alpha(t) + beta_A A`
(no `L` — confounding is handled by the weights). The survival curve is the
cumulative product of the marginal weighted hazards. Sandwich variance is a
stacked estimating equation: treatment-model block + weighted-hazard block.

## The math

Baseline weight for id `i` (stabilized; unstabilized drops the numerator):

```
w_i = f(A_i) / f(A_i | L_i)
```

obtained from causatr's treatment-model machinery (the same density-ratio used
by `causatr`'s point IPW). The weight is **constant within id** and broadcast to
every person-period row of `i`.

Weighted hazard MSM on uncensored at-risk rows, `quasibinomial` family (weights
=> free dispersion, drops the non-integer-successes warning):

```
logit h(t | A = a) = alpha(t) + beta_A a
S^a(t) = prod_{k <= t} (1 - h(k | A = a))
```

Because the MSM is marginal in `A`, `h^a(k)` does not depend on `i`; the
per-individual cumulative product is degenerate and `S^a(t)` reduces to the
product of the fitted marginal hazards. (If the user adds a baseline effect
modifier to the MSM, the per-individual product re-emerges and the chunk-2
average-across-ids machinery handles it unchanged.)

**Variance.** The IF of `S^a(t)` has two stacked blocks: (1) the weighted
hazard-MSM coefficients via the chunk-3 delta chain on the cumulative product,
and (2) the **weight-estimation** correction — the treatment-model score
propagated through the density ratio. Ignoring block (2) is anticonservative.
causatr's IPW IF machinery supplies the treatment-model block (density-ratio
IF); survatr stacks it under the cross-time delta chain. This is the survival
analogue of causatr's point-IPW sandwich.

## Deliverables

### New R files
- `R/ipw_survival.R` — baseline weight via causatr, broadcast onto PP rows,
  weighted MSM fit, and the weight-block IF contribution to the sandwich.

### Updated R files
- `R/surv_fit.R` — enable `estimator = "ipw"` (currently rejected with
  `survatr_bad_estimator`); accept a propensity-model spec (`confounders` feed
  the **treatment** model; the hazard MSM uses `A` only); reuse `model_fn`
  convention for the treatment model (`propensity_model_fn`, default
  `stats::glm`).
- `R/variance_if_survival.R` — add the stacked treatment-model IF block.
- `R/checks.R` — weight-truncation (`trim`) option on the id-level weight
  (winsorize at a quantile before any product; fixed cutoff for sandwich,
  per-replicate for bootstrap — inherit causatr's resolved semantics).

### Tests (`tests/testthat/`)
- `test-ipw-survival.R` — point estimate vs `lmtp::lmtp_tmle(outcome_type =
  "survival")` on a confounded DGP (🟡→🟢 once pinned); marginal-MSM sanity vs
  the unadjusted KM when weights ≡ 1.
- `test-ipw-survival-sandwich.R` — coverage simulation; **assert the stacked SE
  is wider than the naive hazard-only SE** (proves the weight block is wired).
- `helper-ipw-survival-oracle.R` — reusable lmtp comparator (mirror causatr's
  `helper-ipw-lmtp-oracle.R`).

## API contract

```r
fit <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "t", censoring = "C",
  time_formula = ~ splines::ns(t, 4),
  estimator = "ipw",                 # NEW
  propensity_model_fn = stats::glm,  # NEW — treatment model
  trim = NULL                        # NEW — optional weight winsorization
)
# fit$weights holds the per-id broadcast weight; fit$treatment_model the
# propensity fit (for the sandwich + diagnose()).
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "risk_difference",
                   ci_method = "sandwich")
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- Weights are **constant within id** and broadcast to every PP row. Never refit
  the treatment model on person-period rows (that would over-count).
- The hazard MSM is **marginal** in `A` (`alpha(t) + beta_A A`); confounders go
  into the **treatment** model, not the hazard model.
- Stabilized weights are the default. Quasibinomial family for the weighted fit
  (the T6 rule).
- The sandwich **must** include the treatment-model IF block (stacked EE).
- IPSI (`ipsi()`) is an IPW-only intervention: route it through the weight path,
  never through `apply_intervention_pp()` (which aborts on `ipsi`).

## Non-goals (deferred)
- **Longitudinal IPW** (per-period cumulative weights) — composes with Chunk 6.
- **IPCW** (censoring weights) — sibling of this chunk; **pending scope
  ratification** as its own chunk (handoff §7 calls it the motivating IPCW
  case). Do not fold silently into the treatment-weight path.
- AIPW (doubly-robust) — pending scope ratification.

## Dependencies & composition
- causatr: treatment-model + density-ratio IF (the point-IPW internals),
  `prepare_model_if()` for the treatment model. Confirm exact symbol names
  against the current causatr `R/` at implementation time and pin them in
  `test-causatr-integration.R`.

## Acceptance checklist
- [ ] `estimator = "ipw"` fits a weighted marginal hazard MSM; curve matches
      `lmtp` point estimate within tolerance on a confounded DGP.
- [ ] Stacked sandwich SE > hazard-only SE (weight block verified).
- [ ] `weights ≡ 1` reproduces the unadjusted (chunk-2-style) curve.
- [ ] `trim` winsorizes id-level weights; sandwich uses a fixed cutoff.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
