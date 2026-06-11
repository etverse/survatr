# Chunk 29 — Effect modification / `by`-stratification (Track A)

> **Status: ⬜ Not started**
> **Source:** causatr-reuse audit (2026-06-11) — the §4 inheritance table claims
> survatr inherits causatr's `parse_effect_mod()` (already imported for Track B),
> but Track A exposes no `by =` surface.
> **Depends on:** Chunk 2 (survival curve), Chunk 3 (sandwich IF).
> **Oracle:** subgroup curves must equal the curves from fitting the model on
> each subgroup separately (when the hazard model is saturated in the modifier);
> the pooled-over-subgroup average must reproduce the marginal curve.

## Goal

Report subgroup-conditional survival / risk / RMST(/RMTL) curves and their
contrasts within levels of a baseline effect modifier `V`, e.g.
`S^a(t | V = v)`. causatr already exposes `parse_effect_mod()` (and survatr
imports it for Track B), so the missing piece is the Track A surface and the
**subgroup-restricted standardization** of the cumulative-product curve plus
its cross-time IF.

## The math

For a baseline modifier `V` with levels `v`, the subgroup-conditional curve is
the cumulative-product survival averaged over the **subgroup** population:

```
S^a(t | V = v) = (1 / n_v) Σ_{i : V_i = v} ∏_{k ≤ t} (1 − h^a_{i,k})
```

The IF is the chunk-3 per-individual IF restricted to the subgroup target mask
(the Channel-1 `S_i(t) − S̄_v(t)` centring uses the subgroup mean; Channel-2
`J̄_v` averages over the subgroup), with `crossprod(IF) / n_v²`. A modifier ×
treatment interaction in the hazard formula is what makes the subgroup curves
differ — the standardization itself is just a masked average.

## Deliverables

### Updated R files
- `R/contrast.R` — a `by =` argument (a baseline column, constant within id);
  loop the curve + IF construction over its levels, stacking a `subgroup`
  column into `estimates` / `contrasts`.
- `R/variance_if_survival.R` — accept a subgroup target mask so Channel-1
  centring and Channel-2 averaging use the subgroup population.
- `R/contrast_validators.R` — `validate_by()`: column exists, constant within
  id, no NA.
- S3 (`print` / `plot` / `tidy`) — carry the `subgroup` dimension (mirrors the
  competing-risks `cause` handling).

### Tests (`tests/testthat/`)
- `test-effect-modification.R` — subgroup curve equals the separately-fit
  subgroup curve on a saturated model; subgroup average reproduces the marginal;
  subgroup SEs finite, CIs well-formed.

## Behaviour rules (non-negotiable)
- **Reuse causatr's `parse_effect_mod()`** for parsing the modifier spec — do
  not re-roll it (it is already imported for Track B).
- **The modifier is baseline and constant within id**; a time-varying modifier
  is out of scope (different estimand).
- **The standardization is a masked average** on the existing IF — no refit, no
  change to the point hazard model unless the user adds the interaction term.

## Non-goals (deferred)
- Continuous effect modifiers / smooth interaction surfaces.
- Track B (ICE) effect modification — Track B already imports
  `parse_effect_mod()`; a parallel `by` surface there is a later extension.

## Acceptance checklist
- [ ] Subgroup curve equals the separately-fit subgroup curve (saturated model).
- [ ] Subgroup-averaged curve reproduces the marginal curve.
- [ ] `by` validation aborts fire (varies-within-id, NA).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
