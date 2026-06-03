# Chunk 12 — Survival quantiles, RMTL, years-of-life-lost

> **Status: ⬜ Not started**
> **Depends on:** Chunk 2 (survival curve), Chunk 3 (sandwich IF matrix);
> Chunk 7 (competing-risks CIF) for the years-of-life-lost estimand only.
> **Oracle:** `adjustedCurves` (median / quantile survival, RMTL);
> `survival::survfit()` median (unadjusted sanity); closed-form on a
> constant-hazard DGP (`S(t) = exp(-λt)` ⇒ median `= log(2)/λ`).

## Goal

Add the estimands whose absence is the clearest fingerprint of the
scalar-package transplant: **survival-time quantiles / median**, **restricted
mean time LOST (RMTL)**, and **per-cause years-of-life-lost (YLL)**. All three
are smooth functionals of the survival curve `S^a(t)` (or the CIF
`F^(j),a(t)`) that survatr already computes, so they reuse the existing
`n × |t-grid|` IF matrix directly — no new estimation, only a new functional and
its delta-method propagation. These are new `type =` values on `contrast()`.

## The math

All three are functionals `θ^a = g(S^a(·))` of the already-estimated curve; the
sandwich SE is `∇g · V · ∇gᵀ` with `V = crossprod(IF) / n²` the existing
time-covariance and `∇g` the gradient of the functional in the grid values.

**Survival quantile (`type = "quantile"`, with `q =`).** The `q`-quantile of
the survival time under intervention `a` is the smallest `t` with
`S^a(t) ≤ 1 - q` (median survival is `q = 0.5`, i.e. `S^a(t) = 0.5`):

```
τ_q^a = inf { t : S^a(t) ≤ 1 - q }
```

On the discrete grid this is found by interpolation between the bracketing grid
points. The delta-method SE uses the implicit-function gradient: with `S^a`
locally differentiable in `t`, `dτ_q = - dS^a(τ_q) / s'(τ_q)`, where `s'` is the
local slope of the curve (estimated by finite difference on the grid). The
IF on `S^a(τ_q)` is read by interpolating the IF columns at `τ_q`; SE follows
from the existing time-covariance. (Bootstrap CI is the robust fallback when
the curve is near-flat at `τ_q` and `s'` is poorly estimated — offer both
`ci_method = "sandwich"` and `"bootstrap"`.)

**Restricted mean time lost (`type = "rmtl"`).** RMTL up to `t*` is the
complement of RMST within the same window:

```
RMTL^a(t*) = t* - RMST^a(t*) = ∫_0^{t*} (1 - S^a(u)) du
```

The SE reuses the **trapezoidal quadratic-form** machinery already shipped for
RMST (`rmst_weights()`): RMTL is a linear functional with the same trapezoidal
weights `a` applied to `(1 - S^a)`, so `Var(RMTL) = aᵀ V a` is the identical
quadratic form (the constant `t*` and the sign drop out of the gradient). RMTL
difference / ratio between interventions follows the chunk-3 contrast pattern.

**Years-of-life-lost (`type = "yll"`, per cause; needs Chunk 7).** Per-cause
years-of-life-lost up to `t*` is the integral of the cause-`j` CIF:

```
YLL^(j),a(t*) = ∫_0^{t*} F^(j),a(u) du
```

with the decomposition `Σ_j YLL^(j),a(t*) = RMTL^a(t*)` (all-cause time lost
equals the sum of per-cause years lost). The SE uses the chunk-7 CIF IF matrix
with the same trapezoidal weights. This estimand carries a `cause` dimension
on the result, exactly like the chunk-7 CIF estimands, and inherits the
truncation-by-death interpretational caveat.

## Deliverables

### New R files
- `R/estimands_quantile.R` — survival quantile / median solver (grid
  interpolation + implicit-function delta-method gradient), RMTL functional
  (trapezoidal complement of RMST), and the per-cause YLL integral. Each
  functional returns its point estimate and an IF/gradient row so the existing
  variance and bootstrap engines consume it unchanged.

### Updated R files
- `R/contrast.R` — register `type = "quantile"` (with a required `q =`
  argument, default `0.5` = median), `type = "rmtl"`, `type = "yll"` (with a
  `cause =` argument, gated on a competing-risks fit). Dispatch each to the new
  functionals.
- `R/contrast_types.R` — add the three estimand descriptors (point/SE shape,
  whether a `cause` dimension applies, whether `t*` is required).
- `R/contrast_validators.R` — validate `q ∈ (0, 1)`; require `t*` for `rmtl` /
  `yll`; abort `yll` on a non-competing-risks fit (`survatr_yll_needs_cr`);
  abort `quantile` when the curve never crosses `1 - q` on the grid
  (`survatr_quantile_unreached` — the grid does not extend far enough).
- `R/variance_if_survival.R` — gradient propagation for the quantile
  (implicit-function) and the integral functionals (reuse `rmst_weights()`).
- `R/rmst.R` — extract the shared trapezoidal-weight helper so RMTL/YLL reuse
  it rather than duplicating the quadrature.

### Tests (`tests/testthat/`)
- `test-estimands-quantile.R` — constant-hazard DGP: median matches
  `log(2)/λ`; quantiles match the closed form `-log(1-q)/λ`; sandwich vs
  bootstrap CI agree; `survival::survfit()` median sanity on the unadjusted
  curve. `survatr_quantile_unreached` fires on a too-short grid.
- `test-estimands-rmtl.R` — `RMTL = t* - RMST` identity holds to tolerance;
  `Var(RMTL) = Var(RMST)` (same quadratic form, sign drops out); RMTL
  difference matches the negative of the RMST difference.
- `test-estimands-yll.R` (gated on Chunk 7) — `Σ_j YLL^(j) = RMTL` identity;
  per-cause YLL matches `∫ F^(j)` on the analytic two-cause DGP; `adjustedCurves`
  comparator for the adjusted curve quantities.

## API contract

```r
# Median / quantile survival time
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "quantile", q = 0.5,
                   ci_method = "sandwich")
# result$estimates: intervention | q | tau_hat | se | ci_*
# result$contrasts: contrast    | q | estimate | se | ci_*  (difference in median time)

# Restricted mean time lost up to t*
result <- contrast(fit, ..., type = "rmtl", times = seq(0, 120, 12))

# Per-cause years-of-life-lost (competing-risks fit)
result <- contrast(cr_fit, ..., type = "yll", cause = 1, times = seq(0, 120, 12))
# result carries a `cause` dimension (chunk-7 shape)
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **No new estimation.** All three are smooth functionals of the curve +
  IF matrix already produced by chunks 2/3 (and 7 for YLL). Reuse the
  `n × |t-grid|` IF; do not refit anything.
- **RMTL reuses the RMST trapezoidal quadratic form** (`rmst_weights()`); the
  `S(0) = 1` endpoint convention and the `dt[1]/2` constant rule of the existing
  `rmst_weights()` contract carry over unchanged.
- **Quantile SE is the implicit-function delta method**; when the curve is
  near-flat at `τ_q` the sandwich slope estimate is unstable — bootstrap is the
  documented fallback, offered via `ci_method`.
- **YLL requires a competing-risks fit** (`survatr_yll_needs_cr`); it carries
  the `cause` dimension and the truncation-by-death caveat from chunk 7.
- **`Σ_j YLL^(j)(t*) = RMTL(t*)`** and **`RMTL(t*) = t* - RMST(t*)`** identities
  must hold numerically (regression-tested).
- Wald CIs on a quantile or RMTL difference may go negative — a difference is
  not a time. Same intentional symmetric-Wald-on-a-difference invariant that
  governs the risk / RMST difference CIs.

## Non-goals (deferred)
- Quantile residual life / conditional quantiles given survival to `s` —
  beyond v1.
- Simultaneous bands over the quantile process — that is Chunk 16's machinery.
- Restricted-mean **gain** under non-monotone curves — out of scope (the curve
  is a survival probability, monotone non-increasing in expectation).

## Dependencies & composition
- Chunks 2, 3 for quantile + RMTL; additionally Chunk 7 for YLL (the CIF IF
  matrix). No new causatr internals beyond those already imported by chunk 3.

## Acceptance checklist
- [ ] `type = "quantile"` solves `S^a(t) = 1 - q`; median matches
      `log(2)/λ` on the constant-hazard DGP; SE via implicit-function delta and
      bootstrap agree.
- [ ] `type = "rmtl"` satisfies `RMTL = t* - RMST` with the matching quadratic
      form SE.
- [ ] `type = "yll"` (CR fit) satisfies `Σ_j YLL^(j) = RMTL`; matches `∫ F^(j)`
      on the analytic DGP.
- [ ] Off-grid / unreached-quantile aborts fire with distinct classes.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
