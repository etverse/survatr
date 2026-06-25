# Chunk 15 — Parametric doubly-robust (AIPW) survival

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (IPW: treatment model, density-ratio weights,
> stacked-EE treatment block), Chunk 11 (IPCW: censoring model + per-period
> censoring weights), Chunk 3 (sandwich IF); Chunk 7 (competing-risks CIF)
> optional.
> **Oracle:** `riskRegression::ate()` (AIPTW survival with IF-based SE — the
> closest sibling); `stdReg2` (doubly-robust standardization); `lmtp` point
> estimate (TMLE/SDR — point only, EIF-SE not directly comparable).
> **References:** Bai et al. (2013); Zhang & Schaubel (2012).

## Goal

A **parametric** doubly-robust (augmented IPW) estimator of `S^g(t)` that
composes the three nuisances survatr already fits: the outcome **hazard model**
(chunk 2), the **treatment model** (chunk 5), and the **censoring model**
(chunk 11). The estimator is consistent if **either** the outcome model **or**
the (treatment + censoring) models are correctly specified — strictly stronger
than the singly-robust gcomp (chunk 2) and IPW (chunk 5) it builds on. The
boundary is explicit and load-bearing: **parametric nuisances only.** ML /
TMLE / cross-fitted survival stays OUT of survatr — defer to `lmtp` and
`concrete`. This is a large chunk (the scale of chunk 7): it wires a new
estimator, a new stacked-EE sandwich across all nuisance blocks plus the
augmentation term, and an Invariants section pinning the double-robustness
guarantees.

## The math (Bai et al. 2013; Zhang & Schaubel 2012; handoff §7)

The augmented estimator standardizes the outcome model and adds a
weighted residual correction:

```
Ŝ^g_AIPW(t) = (1/n) Σ_i [ Ŝ^g_i(t)  +  W_i(g) · ( 𝟙{T_i > t, C_i > t} − Ŝ^{A_i}_i(t) ) ]
```

Reading the three pieces:

- **`Ŝ^g_i(t)`** — the chunk-2 g-computation prediction for id `i` under the
  intervention `g` (cumulative product of intervened hazards). This is the
  outcome-model term; on its own it is the singly-robust gcomp curve.
- **`Ŝ^{A_i}_i(t)`** — the outcome-model prediction at id `i`'s **observed**
  treatment `A_i` (cumulative product of hazards at the factual treatment). The
  augmentation centres the residual on this factual prediction.
- **`W_i(g)`** — the combined IPTW × IPCW weight for id `i` under `g`: the
  chunk-5 treatment density ratio (`𝟙{A_i = a}/f(A_i | L_i)` for a static
  intervention, the appropriate ratio for shift / ipsi) times the chunk-11
  per-period cumulative censoring weight `W^C_{i,t}` evaluated at `t`. The
  indicator `𝟙{T_i > t, C_i > t}` is the observed survival-and-uncensored event;
  the weight up-weights the observed survivors that stand in for the censored.

The augmentation term has **mean zero** when the weight model is correct
(IPW-consistent) and the residual `(𝟙{·} − Ŝ^{A_i}_i(t))` has mean zero when
the outcome model is correct (gcomp-consistent) — hence the double robustness.

**Reduction properties (these are the Invariants below):**

- If the weights are correct but the outcome model is misspecified, the
  augmentation corrects the bias → consistent (the IPW arm carries it).
- If the outcome model is correct but the weight model is misspecified,
  `E[W_i(g) · residual] = 0` because the residual is mean-zero → the estimator
  collapses to gcomp → consistent (the gcomp arm carries it).
- If the weights are set to the **correct** known values and the residual is
  exactly zero (saturated outcome model), the estimator reduces to plain
  gcomp; with `W ≡ 0` it is exactly gcomp; with a saturated outcome model the
  augmentation reproduces the IPW estimator.

**Variance (stacked EE).** The sandwich stacks every nuisance score plus the
augmentation estimating equation:

1. **outcome-hazard block** — the chunk-3 cross-time delta chain on the
   cumulative product (for both `Ŝ^g_i` and `Ŝ^{A_i}_i`);
2. **treatment-model block** — the chunk-5 density-ratio IF;
3. **censoring-model block** — the chunk-11 cross-time IF on the running-product
   IPCW weight;
4. **augmentation IF** — the influence contribution of the residual term
   itself, which couples all three nuisance blocks through `W_i(g)` and the two
   outcome predictions.

The per-individual IF is the usual `n × |t-grid|` matrix; full time-covariance
is `crossprod(IF) / n²`. Because the estimator is doubly robust, the IF is the
**efficient-style** combination of the three blocks — but it is the parametric
M-estimation IF (Bai et al. 2013; Zhang & Schaubel 2012), NOT the TMLE / SDR
EIF that `lmtp` / `riskRegression` report, so `lmtp` is a point-estimate oracle
only (the SE is not directly comparable, same rule as chunk 5/6).

## Deliverables

### New R files
- `R/aipw_survival.R` — the AIPW estimator: assemble `Ŝ^g_i(t)`, `Ŝ^{A_i}_i(t)`,
  the combined IPTW × IPCW weight `W_i(g)`, the survived-and-uncensored
  indicator, the augmented average, and the augmentation IF contribution to the
  stacked sandwich. Handles the CR variant (augment the CIF) when the fit is
  competing-risks.

### Updated R files
- `R/surv_fit.R` — `estimator = "aipw"` (currently rejected with
  `survatr_bad_estimator`): require the outcome-model spec (hazard model),
  the treatment-model spec (`propensity_model_fn`, chunk 5) and the
  censoring-model spec (`ipcw`, chunk 11); store all three nuisance fits on the
  `survatr_fit`.
- `R/contrast.R` — AIPW curve assembly (the augmented average) and the
  time-indexed result; CIF-AIPW when CR.
- `R/variance_if_survival.R` — assemble the four stacked blocks (outcome
  hazard, treatment, censoring, augmentation) into the per-individual IF.
- `R/checks.R` — reuse the chunk-5/11 `trim` plumbing on `W_i(g)`; validate that
  all three nuisance specs are present for `estimator = "aipw"`.

### Tests (`tests/testthat/`)
- `test-aipw-survival.R` — **double-robustness simulation (the headline
  test):** four cells — (outcome ✓, weights ✓), (outcome ✓, weights ✗),
  (outcome ✗, weights ✓), (outcome ✗, weights ✗); AIPW is consistent in the
  first three and biased only in the fourth, while gcomp is biased in cells 3–4
  and IPW in cells 2 + 4. Point estimate vs `riskRegression::ate()` (AIPTW) and
  `lmtp` on a correctly-specified DGP.
- `test-aipw-survival-sandwich.R` — coverage simulation under each
  double-robustness cell; the four-block stacked SE covers nominally when the
  estimator is consistent; cross-check the IF-SE structure against
  `riskRegression::ate()`'s IF-based SE (point + approximate-shape comparator,
  documented as not identical due to the M-estimation vs IF derivation).
- `test-aipw-reductions.R` — the Invariant reductions: `W ≡ 0` ⇒ exactly
  gcomp; saturated outcome model ⇒ reproduces IPW; correct weights + correct
  outcome ⇒ matches both singly-robust arms within tolerance.
- `helper-aipw-survival-oracle.R` — reusable `riskRegression::ate` + `stdReg2`
  + `lmtp` comparators.

## API contract

```r
fit <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "t", censoring = "C",
  time_formula = ~ splines::ns(t, 4),
  estimator = "aipw",                    # NEW — doubly robust
  model_fn = stats::glm,                 # outcome hazard model (gcomp arm)
  propensity_model_fn = stats::glm,      # treatment model (chunk 5)
  ipcw = ~ L1 + L2,                      # censoring model (chunk 11)
  censoring_model_fn = stats::glm,
  trim = NULL
)
# fit holds all three nuisance fits: outcome hazard, treatment, censoring.
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "risk_difference",
                   ci_method = "sandwich")
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Parametric nuisances ONLY.** No ML / Super Learner / TMLE / cross-fitting in
  survatr. The DR boundary is parametric AIPW (Bai et al. 2013; Zhang &
  Schaubel 2012); ML/TMLE survival is `lmtp` / `concrete`'s job and stays out of
  scope. Document this boundary in roxygen and the vignette.
- **All three nuisance specs are required** for `estimator = "aipw"` (outcome
  hazard, treatment model, censoring model); a missing spec aborts upfront
  (`survatr_aipw_missing_nuisance`).
- **The augmentation centres on the FACTUAL prediction** `Ŝ^{A_i}_i(t)` (outcome
  model at the observed treatment), not the counterfactual `Ŝ^g_i(t)` — the
  residual `𝟙{T_i > t, C_i > t} − Ŝ^{A_i}_i(t)` must be mean-zero under a
  correct outcome model for the gcomp arm to carry consistency.
- **The combined weight is IPTW × IPCW** (chunk-5 density ratio × chunk-11
  per-period cumulative censoring weight); reuse both, do not re-derive.
- **The sandwich MUST stack all four blocks** (outcome, treatment, censoring,
  augmentation). Treating any nuisance as fixed is anticonservative.
- **`lmtp` / `riskRegression` are point-estimate oracles** for variance: their
  EIF / IF-based SE is not identical to the M-estimation sandwich (same rule as
  chunks 5/6).
- Quasibinomial family for any weighted fit (the T6 rule).

## Invariants (double-robustness contract — regression-tested)

- **AIPW is consistent when at most one nuisance arm is misspecified.** Correct
  outcome ∨ correct (treatment + censoring) ⇒ consistent. (Tested by the
  four-cell double-robustness simulation.)
- **`W_i(g) ≡ 0` reduces AIPW exactly to gcomp** — the augmentation vanishes.
- **A saturated outcome model reduces the augmentation to the IPW estimator** —
  the residual carries all the information.
- **Correct weights + correct outcome ⇒ AIPW agrees with both singly-robust
  arms** (gcomp and IPW) within Monte-Carlo tolerance.
- **The augmentation term is mean-zero** under a correct outcome model
  (numerically: `mean(W_i · residual) → 0` as `n → ∞`).

## Non-goals (deferred)
- **ML / TMLE / SDR / cross-fitted** survival estimators — permanently out of
  scope (the explicit boundary above); defer to `lmtp` / `concrete`.
- **AIPW under longitudinal ICE-hazard** (longitudinal DR) — composes after Chunk 6; far more
  involved (sequential DR), flagged as a future major lift, not this chunk.
- **One-step / targeting iteration** (TMLE-style fluctuation) — that is the ML
  arm and stays out.

## Dependencies & composition
- Chunk 5 (treatment model + density-ratio weights + treatment IF block);
  Chunk 11 (censoring model + per-period IPCW weight + censoring IF block);
  Chunk 3 (cross-time delta chain); Chunk 7 (CIF) for the optional CR-AIPW
  variant.
- causatr: `causatr:::prepare_model_if()` for each parametric nuisance block
  (already imported by chunks 3/5). Confirm exact symbol names against the
  current causatr `R/` at implementation time and pin them in
  `test-causatr-integration.R`.

## Acceptance checklist
- [ ] `estimator = "aipw"` fits the three parametric nuisances and returns the
      augmented curve; matches `riskRegression::ate()` (AIPTW) + `lmtp` point
      estimate on a correctly-specified DGP.
- [ ] Four-cell double-robustness simulation: consistent when outcome ✓ OR
      weights ✓; biased only when both ✗.
- [ ] Reduction invariants hold: `W ≡ 0` ⇒ gcomp; saturated outcome ⇒ IPW;
      correct-both ⇒ agrees with singly-robust arms.
- [ ] Four-block stacked sandwich covers nominally when the estimator is
      consistent.
- [ ] Parametric-only boundary documented (ML/TMLE explicitly out, deferred to
      `lmtp` / `concrete`).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
