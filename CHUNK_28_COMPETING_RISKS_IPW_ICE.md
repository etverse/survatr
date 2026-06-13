# Chunk 28 — Competing risks under IPW / ICE

> **Status: ⬜ Not started**
> **Source:** causatr-reuse audit (2026-06-11) — chunk 7 fit cause-specific
> hazards under gcomp only.
> **Depends on:** Chunk 5 (IPW), Chunk 6 (ICE), Chunk 7 (cause-specific hazards
> + CIF sandwich).
> **Oracle:** confounded competing-risks DGP — the IPW / ICE CIF must match a
> delicatessen stacked-EE oracle (point ~1e-4) and the dual-refit bootstrap;
> under no confounding it must collapse to the chunk-7 gcomp CIF.

## Goal

Chunk 7 fit the J cause-specific hazards on one shared all-cause risk set and
built the block-diagonal CIF sandwich — but under **gcomp only**. This chunk
adds the two adjusted estimators: weighted cause-specific hazard MSMs (IPW) and
per-cause ICE pseudo-outcomes (ICE), reusing the chunk-7 risk-set and CIF
decomposition.

## Scope

- **CR + IPW.** Each cause-`j` hazard model is fit with the chunk-5 stabilized
  treatment weight broadcast onto the shared all-cause risk set. The CIF
  sandwich gains the chunk-5 treatment-model correction per cause (the bread
  stays block-diagonal across causes; the treatment block is shared).
- **CR + ICE.** A per-cause survival-tail pseudo-outcome on the backward sweep;
  the cross-step cascade carries the chunk-7 all-cause `(1 − H)` sensitivity.

## Reuse (do not re-derive)
- `build_risk_set()` once (shared all-cause), `fit_competing_risks()` per cause
  (chunk 7) — only the per-cause fit step gains weights / ICE.
- The chunk-7 CIF derivative with the lagged cumulative sensitivity
  (`cumSX_lag = cumSX − SX`) and the all-cause `(1 − H)` denominator are
  unchanged — they are estimator-agnostic.
- chunk-5 `prepare_ipw_correction()` / `apply_model_correction()` for the IPW
  block; chunk-6 `survatr_ice_surv_chain()` for the ICE cascade.

## Deliverables

### Updated R files
- `R/competing_risks.R` — accept `estimator %in% c("ipw", "ice")` (chunk 7
  currently aborts with `survatr_competing_estimator`).
- `R/variance_if_competing.R` — add the treatment-model (IPW) / cross-step (ICE)
  blocks to the per-cause IF.

### Tests (`tests/testthat/`)
- `test-competing-risks-ipw.R`, `test-competing-risks-ice.R` — delicatessen +
  bootstrap agreement; collapse-to-gcomp under no confounding.

## Behaviour rules (non-negotiable)
- **One shared all-cause risk set** (chunk-7 invariant) — do not fit per-cause
  risk sets under IPW / ICE either.
- **The CIF sensitivity denominator stays all-cause `(1 − H)`**, not
  `(1 − h^(j))` — the weighting / ICE changes the fit, not the delta geometry.
- **Fine–Gray stays out of scope** — cause-specific hazards + CIF only.

## Acceptance checklist
- [ ] CR + IPW and CR + ICE match a delicatessen oracle (point) + bootstrap (SE).
- [ ] Both collapse to the chunk-7 gcomp CIF under no confounding.
- [ ] The chunk-7 rejection (`survatr_competing_estimator`) is narrowed, not
      removed (any still-unsupported estimator stays classed).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
