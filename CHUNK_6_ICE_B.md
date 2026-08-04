# Chunk 6 — Longitudinal ICE-hazard: survival via ICE hazards

> **Status: ✅ Done**
> **Depends on:** Chunk 3 (sandwich IF), causatr ICE per-step primitives
> (`ice_fit_step()`, `ice_predict_step()`, `ice_build_formula()`,
> `ice_apply_intervention_long()`, `create_lag_vars()`, `expand_em_lag_terms()`)
> and IF primitives (`ice_if_setup()`, `correct_model()`, `iv_design_matrix()`).
> **Oracle:** forward-simulation Monte-Carlo g-formula truth (primary point
> oracle); `delicatessen` stacked M-estimation (variance oracle, the gold
> standard); `lmtp::lmtp_tmle(outcome_type = "survival")` and
> `gfoRmula::gformula_survival()` (external cross-checks when installed).

## Goal

Time-varying treatment, time-to-event outcome via iterated conditional
expectations (Zivich et al. 2024) extended to the **hazard link**. Reuse
causatr's ICE engine wholesale; add only the survival-aware per-step target,
per-step link forcing, and the cumulative-product survival-curve shape.

### Why ICE, not forward simulation

Hernán & Robins Ch. 21 describes two longitudinal g-formula approaches.
Forward simulation (`gfoRmula`-style) needs models for every time-varying
covariate + outcome and bootstrap-only inference. ICE models the outcome only
at each time, iterates backward with a pseudo-outcome, and admits a stacked
estimating-equation sandwich. Here the pseudo-outcome is the **survival tail**
under the intervention; the per-step link is binomial at K, quasibinomial at
k < K.

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

**Variance: reuse causatr's single-model IF primitives, own the cross-step
chain.** The original plan was to feed a survatr-built `ice_result` straight to
`causatr:::variance_if_ice_one()` and reuse its forward sensitivity recursion
verbatim. That is **not** valid for the survival pseudo-outcome, and the
sandwich-vs-bootstrap test caught it: causatr's chain (built for a single
terminal `Y`) sums the cross-step gradient over `fit_ids[[k]]` = at-risk-at-`k`,
which **includes** individuals who have the event at `k` (`D_k = 1`). The
survival-tail derivative `dỸ_k/dβ_{k+1} = (1 - D_k) dq_{k+1}/dβ_{k+1}` carries a
`(1 - D_k)` factor that those failed individuals should zero out; without it the
sandwich over-covers, increasingly at later horizons (more accumulated
failures). survatr therefore reuses causatr's *single-model* IF pieces
(`ice_if_setup()` for Channel 1, `correct_model()` / `iv_design_matrix()` /
`coef_clean()` for each step's bread and score) but owns the cross-step cascade
(`survatr_ice_surv_chain()`), injecting the `(1 - D_k)` factor via per-period
event indicators. This mirrors the point-treatment g-computation architecture (causatr owns
single-model IF primitives; survatr owns the cross-time aggregation). Validated
to ~1e-5 against an independent `delicatessen` numerical sandwich.

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
  confounders = ~ V1 + V2,            # BASELINE (time-invariant) confounders
  confounders_tv = ~ L1 + L2,         # NEW — time-varying confounders, lag-expanded
  id = "id", time = "t", censoring = "C",
  estimator = "ice",                  # NEW — longitudinal ICE-hazard
  history = Inf                       # NEW — Markov lag order (default: full history)
)
# fit$track == "B"
result <- contrast(fit, interventions = list(d = causatr::static(1)),
                   times = ..., type = "risk_difference",
                   ci_method = "sandwich")  # survival-aware stacked-EE ICE sandwich
```

**Confounders API.** The longitudinal ICE path splits baseline (`confounders`, never lagged) from
time-varying (`confounders_tv`, lag-expanded at each backward step) confounders,
matching causatr's `fit_ice()`. `history` sets the Markov lag order (`Inf` =
full history, capped at `n_times - 1`). Both args are specific to the longitudinal ICE path and ignored
by point-treatment g-computation.

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Per-step link forcing is load-bearing:** binomial at K (0/1 indicator),
  quasibinomial at `k < K` (survival-tail pseudo-outcome in [0, 1]). Swapping
  them silently changes the score equations.
- **ICE intervention applies to current-time treatment only.** Lag columns hold
  the **observed** `A_{k-1}, A_{k-2}, ...` at every backward step (inherited
  causatr rule). Do **not** recompute lags from the intervened treatment.
- **Pseudo-outcome is the survival tail** `Ỹ_k = D_k + (1 - D_k) q_{k+1} =
  1 - (1 - D_k)(1 - q_{k+1})`. The `D_k` term (the failure carry-forward) is the
  piece causatr's plain ICE lacks; the NA-safe `ifelse(D_k == 1, 1, q_{k+1})`
  form avoids `0 * NA`.
- **survatr owns the backward loop** (`run_ice_survival_horizon()`) and the
  cross-step IF cascade (`survatr_ice_surv_chain()`), reusing causatr's
  *single-model* primitives (`ice_fit_step`, `ice_predict_step`,
  `ice_build_formula`, `ice_if_setup`, `correct_model`). Do **not** call
  `causatr:::ice_iterate()` / `variance_if_ice_one()` — they are built for a
  single terminal outcome and omit the `(1 - D_k)` factor.
- **`fit_ids[[k]]` is the at-risk-at-`k` fit set** (the rows model `k` is fit
  on). The cross-step gradient additionally multiplies by `(1 - D_{k-1})` so
  individuals who fail at the previous step do not propagate the next model's
  sensitivity.
- One backward pass **per horizon** (the survival-tail target is relative to a
  fixed final period); the cumulative-product survival shape is `S^d(t) =
  1 - R^d(t)` with `R^d(t) = mean(pseudo_final)` (Jensen-safe — the product is
  inside each per-individual pseudo-outcome).

## Non-goals (deferred)
- Forward-simulation g-formula (that is `gfoRmula`'s job — an oracle, not a
  feature here).
- Longitudinal IPW survival (the IPW analogue) — separate.
- Competing risks under longitudinal ICE-hazard — composes with Chunk 7 after both land.

## Dependencies & composition
- causatr per-step + IF primitives via `causatr:::` — `ice_fit_step()`,
  `ice_predict_step()`, `ice_build_formula()`, `ice_apply_intervention_long()`,
  `create_lag_vars()`, `ice_if_setup()`, `correct_model()`, `iv_design_matrix()`,
  `coef_clean()`, `new_causatr_fit()`, `parse_effect_mod()`,
  `has_stochastic_component()` — pinned in `test-causatr-integration.R`.

## Acceptance checklist
- [x] Track-B risk curve matches the forward-simulation g-formula truth on a
      treatment-confounder-feedback DGP (primary oracle); `gfoRmula` /
      `lmtp` cross-checks when installed.
- [x] Static strategies match `lmtp` / forward-sim truth point estimates.
- [x] Per-step link is binomial at K, quasibinomial at `k < K` (asserted).
- [x] Lag columns carry observed (not intervened) treatment (asserted).
- [x] Survival-aware stacked-EE sandwich populated; matches the `delicatessen`
      M-estimation oracle to ~1e-5 and agrees with the bootstrap.
- [x] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
