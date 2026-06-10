# Chunk 11 — Built-in IPCW (per-period censoring weights)

> **Status: ✅ Done (`cdcdbae`; delicatessen oracle `15ea801`)**
> **Depends on:** Chunk 5 (IPW machinery: density-ratio weights, stacked-EE
> treatment block, `trim`), Chunk 3 (sandwich IF), Chunk 2 (contrast spine).
> **Oracle:** `lmtp::lmtp_tmle(outcome_type = "survival", cens = <col>)` (point
> estimate under informative censoring); Cole & Hernán (2004) adjusted survival
> curves; `weights ≡ 1` reproduces the chunk-2 / chunk-5 curve.

## Goal

Survival is the **motivating** use case for inverse-probability-of-censoring
weighting (IPCW), not a supported afterthought. When censoring depends on
measured covariates (and, for Track B, on time-varying history), the
row-filter `censoring =` semantics of chunks 1–4 are biased. This chunk builds
a **censoring hazard model** on the person-period grid, forms **per-period
cumulative** inverse-probability-of-censoring weights, and multiplies them into
the weighted hazard MSM. The stacked estimating equation extends with the
censoring-model block, exactly paralleling chunk 5's treatment-weight block but
keyed on `(i, k)` rather than `i`.

## The math

Let `C_{i,k} = 1` flag censoring of id `i` in period `k`. Fit a per-period
**censoring hazard** on at-risk, not-yet-censored rows:

```
logit P(C_k = 1 | uncensored & at risk through k, A, L̄_k) = gamma(k) + delta_A A + delta_L L̄_k
```

(`mgcv::gam` / `splines::ns` for `gamma(k)`, same `model_fn` contract as the
hazard model). The per-period probability of remaining uncensored is
`1 - ĝ_{i,k}`; the **cumulative** inverse-probability-of-censoring weight for id
`i` at period `k` is the running product

```
W^C_{i,k} = prod_{m <= k} ( s(m) / (1 - ĝ_{i,m}) )
```

where `s(m)` is the stabilizing numerator (a marginal / treatment-only
censoring model, or `1` for the unstabilized weight). Stabilized is the
default (Cole & Hernán 2004) — the numerator removes the structural variance
inflation from the denominator's late-time products.

The weight is **per-period (per `(i, k)` row)**, NOT constant within id: it
grows as more periods accumulate. This is the key structural difference from
chunk 5's id-level treatment weight `w_i`. Under a combined IPTW + IPCW
analysis the row weight is the product `w_i · W^C_{i,k}` (treatment weight
broadcast onto each row, times the per-period censoring weight).

Weighted hazard MSM (quasibinomial — the T6 family switch) on the **uncensored
at-risk** rows, with the cumulative IPCW weight up-weighting the survivors that
stand in for those censored:

```
logit h(t | A = a) = alpha(t) + beta_A a            [weights = w_i · W^C_{i,k}]
S^a(t) = prod_{k <= t} (1 - h(k | A = a))
```

**Variance.** The IF of `S^a(t)` stacks three blocks: (1) the weighted
hazard-MSM coefficients via the chunk-3 cross-time delta chain on the
cumulative product; (2) the **treatment-model** correction (chunk 5's
density-ratio IF); and (3) the **censoring-model** correction — the
censoring-hazard score propagated through the cumulative product
`W^C_{i,k} = prod_{m <= k} 1/(1 - ĝ_{i,m})`. Block (3) is itself a cross-time
delta chain (the weight is a running product of estimated terms), so its IF
contribution accumulates over periods within id before entering the sandwich.
Ignoring block (3) is anticonservative for the same reason ignoring the
treatment block is.

## Deliverables

### New R files
- `R/ipcw_survival.R` — censoring-hazard fit on the PP grid, per-period
  cumulative IPCW weight construction (`W^C_{i,k}` running product), stabilizing
  numerator, combination with the chunk-5 treatment weight, and the
  censoring-block IF contribution to the stacked sandwich.

### Updated R files
- `R/surv_fit.R` — accept a censoring-model spec: `ipcw = TRUE` (or a one-sided
  formula for the censoring-model covariates), `censoring_model_fn` (default
  `stats::glm`, mirrors `propensity_model_fn`). When IPCW is on, `censoring =`
  switches from a pure row filter to the modelled path; the at-risk row set is
  unchanged (fit on uncensored rows), the weight does the reweighting.
- `R/ipw_survival.R` — multiply the per-period IPCW weight into the per-row MSM
  weight (broadcast treatment weight × per-period censoring weight).
- `R/variance_if_survival.R` — add the stacked censoring-model IF block (a
  cross-time delta chain on the running-product weight).
- `R/checks.R` — extend the chunk-5 `trim` semantics to the **per-period**
  censoring weight: winsorize `W^C_{i,k}` at a quantile **per period** (fixed
  cutoff for sandwich, per-replicate for bootstrap), adopted from day one — the
  late-time products are the heaviest tails.

### Tests (`tests/testthat/`)
- `test-ipcw-survival.R` — informative-censoring DGP (censoring depends on `L`):
  IPCW curve matches `lmtp::lmtp_tmle(outcome_type = "survival", cens = <col>)`
  point estimate within tolerance (🟡→🟢 once pinned); the uncorrected
  row-filter fit is visibly biased; `weights ≡ 1` (non-informative censoring)
  reproduces the chunk-2 curve.
- `test-ipcw-survival-sandwich.R` — coverage simulation; **assert the
  three-block stacked SE exceeds the two-block (treatment-only) SE** (proves the
  censoring block is wired), which in turn exceeds the hazard-only SE.
- `helper-ipcw-survival-oracle.R` — reusable `lmtp` comparator with the `cens =`
  argument (mirror the chunk-5 oracle helper).

## API contract

```r
fit <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "t", censoring = "C",
  time_formula = ~ splines::ns(t, 4),
  estimator = "ipw",                       # IPCW rides the chunk-5 IPW path
  propensity_model_fn = stats::glm,        # treatment model (chunk 5)
  ipcw = ~ L1 + L2,                        # NEW — censoring-model covariates
  censoring_model_fn = stats::glm,         # NEW — censoring hazard model
  trim = NULL                              # winsorize per-period IPCW weights
)
# fit$censoring_model holds the censoring-hazard fit (for the sandwich +
# diagnose()); fit$ipcw_weights holds the per-(i,k) cumulative weight.
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "risk_difference",
                   ci_method = "sandwich")
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **IPCW weights are per-period (per `(i, k)` row), NOT constant within id.**
  They are a running product over periods. This is the structural distinction
  from chunk-5's id-level treatment weight — do not broadcast a single per-id
  censoring weight.
- **Stabilized weights are the default** (Cole & Hernán 2004); the unstabilized
  form drops the numerator.
- **Quasibinomial family** for the weighted hazard fit (the T6 rule).
- The censoring hazard is fit on at-risk, **not-yet-censored** rows; the hazard
  MSM is still fit on uncensored at-risk rows — the IPCW weight reweights the
  survivors, it does not change which rows are fit.
- **The sandwich MUST include the censoring-model IF block** (third stacked EE
  block). Treating estimated IPCW weights as fixed is anticonservative.
- **Per-period truncation from day one.** The `trim` winsorization is applied
  per period (the late-time cumulative products carry the heaviest tails); fixed
  cutoff for sandwich, per-replicate for bootstrap (chunk-5 semantics).
- When IPCW is active, `censoring =` is no longer a pure row filter — document
  the switch in roxygen so the chunk-1 row-filter contract is not silently
  overridden.

## Non-goals (deferred)
- **Longitudinal IPCW under Track B** (per-period censoring weights composed
  with the ICE pseudo-outcome) — composes after Chunk 6; the per-period weight
  machinery here is the prerequisite.
- **Competing-risks IPCW** (a censoring weight alongside cause-specific
  hazards) — composes after Chunk 7.
- **Doubly-robust IPCW** (augmenting the IPCW estimator with an outcome model) —
  that is the AIPW chunk; this chunk is singly-robust in the censoring nuisance.

## Dependencies & composition
- Chunk 5 supplies the IPW weighted-MSM spine, the stacked-EE treatment block,
  and the `trim` plumbing; this chunk reuses all of it and adds the censoring
  block parallel to the treatment block.
- causatr: `causatr:::prepare_model_if()` for the censoring-model IF (same
  primitive used for the treatment model in chunk 5). Confirm exact symbol
  names against the current causatr `R/` at implementation time and pin them in
  `test-causatr-integration.R`.

## Acceptance checklist
- [x] Censoring hazard model fits on the PP grid; per-period cumulative IPCW
      weight `W^C_{i,k}` is a running product (asserted to grow within id).
- [x] Under informative censoring the IPCW curve matches the `lmtp` point
      estimate within tolerance; the uncorrected row-filter fit is biased
      (direction test: naive > IPCW in >70% of 25 runs, δ=0.8).
- [x] `weights ≡ 1` (non-informative censoring, δ=0) reproduces the naive
      row-filter IPW curve (within 0.03); IPCW weights cluster around 1.
- [x] Censoring block wired: `A_beta_gamma ≠ 0`; `IF_gamma_per_id ≠ 0`;
      three-block IF matrix differs from two-block by >1e-6. *(Note: the
      original spec said three-block SE > two-block SE > hazard-only SE, but
      the censoring correction has no guaranteed sign direction — see
      hard-rules.md. The IF-matrix difference check is the correct assertion.
      Validated to <2% against the delicatessen oracle.)*
- [x] Per-period `trim` winsorizes the cumulative weights; sandwich uses a
      fixed cutoff (re-quantiling inside the numDeriv closure would break
      smoothness in γ).
- [x] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md + hard-rules.md
      updated. delicatessen Python oracle committed (`test-ipcw-delicatessen.R`,
      `data-raw/delicatessen_ipcw_survival.py`).
