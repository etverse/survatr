# Chunk 27 — IPCW for gcomp and ICE

> **Status: ⬜ Not started**
> **Source:** causatr-reuse audit (2026-06-11) — gap between chunk 11 (IPW-only
> IPCW) and chunks 2/3 (gcomp) + 6 (ICE).
> **Depends on:** Chunk 6 (ICE), Chunk 11 (IPCW running product + three-block
> sandwich).
> **Oracle:** non-informative censoring (`delta_cens = 0`) must reproduce the
> uncorrected gcomp / ICE curve; informative censoring must match a delicatessen
> stacked-EE oracle (point ~1e-4) and the dual-refit bootstrap (~15% on SEs).

## Goal

Chunk 11 built the per-period cumulative IPCW weight `W^C_{i,k}` and the
three-block stacked-EE sandwich (beta + alpha_trt + gamma_cens), but coupled
them to the **IPW** path only. gcomp + IPCW and ICE + IPCW are the natural
remaining combinations: a g-computation or ICE hazard model under informative
right-censoring needs the same censoring reweighting. This chunk extends the
chunk-11 machinery to those two estimators.

## Scope

- **gcomp + IPCW.** Multiply the censoring running product onto the gcomp fit
  rows; the sandwich gains the chunk-11 gamma_cens block (no treatment-model
  block, since gcomp does not estimate a propensity). The cross-time delta and
  cumulative-product curve are unchanged.
- **ICE + IPCW.** The backward sequential regression is fit on
  censoring-reweighted rows; the cross-step IF cascade gains the censoring
  block. This is the subtler half — the running product interacts with the
  per-horizon pseudo-outcome.

## Reuse (do not re-derive)
- `compute_ipcw_running_weights()` / `ipcw_running_cumprod()` (chunk 11) — the
  running product is estimator-agnostic.
- The chunk-11 censoring-model fit-row convention (`prev_event == 0 &
  prev_cens == 0`, includes the `C_k = 1` row) and the `n_ids / n_cens_fit`
  bread scaling carry over verbatim.

## Deliverables

### Updated R files
- `R/ipcw_survival.R` — allow `estimator %in% c("gcomp", "ice")` alongside
  `"ipw"` when `ipcw =` is supplied.
- `R/variance_if_survival.R` / `R/variance_if_ice_survival.R` — add the
  censoring correction block to the gcomp and ICE IF chains.

### Tests (`tests/testthat/`)
- `test-ipcw-gcomp.R`, `test-ipcw-ice.R` — non-informative reproduction;
  informative-censoring delicatessen + bootstrap agreement.

## Behaviour rules (non-negotiable)
- **Reuse the chunk-11 running product + fit-row convention** unchanged; only
  the MSM / ICE fit step and its IF block differ by estimator.
- **gcomp + IPCW has no treatment-model block** (two-block sandwich:
  beta + gamma_cens); do not add a spurious alpha_trt correction.

## Acceptance checklist
- [ ] gcomp + IPCW and ICE + IPCW reproduce the uncorrected curve under
      non-informative censoring.
- [ ] Both match a delicatessen oracle (point) and the bootstrap (SE).
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
