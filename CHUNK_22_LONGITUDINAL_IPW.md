# Chunk 22 — Longitudinal IPW survival (time-varying treatment MSM)

> **Status: ⬜ Not started**
> **Depends on:** Chunk 5 (point IPW core), Chunk 11 (IPCW; longitudinal
> follow-up is censored), causatr longitudinal-IPW internals
> (`fit_longitudinal_ipw()`, per-period density chain).
> **Oracle:** `lmtp::lmtp_tmle` / `lmtp::lmtp_sdr` on a time-varying-treatment
> survival DGP; causatr's longitudinal-IPW tests for the per-period weights.

## Goal

Extend IPW survival from a **point** (baseline-constant) treatment to a
**time-varying** treatment with treatment-confounder feedback — the Robins
marginal structural model for survival (Hernán & Robins Ch. 17/21). This is the
home for the `survatr_ipw_time_varying_treatment` rejection that chunk 5 raises:
chunk 5 only weights a single per-id baseline weight, so a treatment that varies
within id must route here.

## The recipe (handoff architecture §)

Per-period treatment density models produce **cumulative** density-ratio weights
indexed by `(i, k)`:

```
SW_{i,k} = prod_{l <= k} f(A_{i,l}) / f(A_{i,l} | A_{i,<l}, L_{i,l})
```

(stabilized; the numerator may condition on baseline covariates only). The
weight is **no longer constant within id** — it accumulates across periods. The
weighted pooled-logistic hazard MSM is fit on uncensored person-period rows with
the per-period cumulative weight; the survival curve is the cumulative product
of the marginal weighted hazards. This is the IPW analogue of ICE-survival
(chunk 6).

## What changes from chunk 5

- The treatment model becomes a **per-period chain** (reuse causatr's
  `fit_longitudinal_ipw()` / per-period density machinery), not a single
  baseline fit. The weight is `(i, k)`-indexed, not broadcast from one per-id
  value — `broadcast_weights_to_pp()` is replaced by the cumulative product.
- Censoring is intertwined: longitudinal follow-up is censored, so this composes
  with **IPCW** (chunk 11) — the cumulative treatment weight multiplies the
  cumulative censoring weight.
- **Variance.** The stacked EE gains one alpha-block per period (block-diagonal
  per-period propensity bread) plus the within-id cumulative-weight dependence;
  the chunk-5 single-block cross-derivative generalizes to a sum over periods.

## Deliverables
- `R/ipw_longitudinal_survival.R` — per-period weight chain + cumulative product
  + weighted MSM.
- `R/variance_if_ipw_survival.R` — per-period treatment-correction blocks.
- Tests: point estimate vs `lmtp` on a feedback DGP; sandwich vs two-stage
  bootstrap; degenerate-to-chunk-5 check when treatment is baseline-constant.

## Acceptance checklist
- [ ] Time-varying treatment no longer aborts; fits the longitudinal MSM.
- [ ] Curve matches `lmtp::lmtp_tmle` on a treatment-confounder-feedback DGP.
- [ ] Reduces to the chunk-5 point estimate when treatment is constant within id.
- [ ] Stacked sandwich matches the two-stage bootstrap.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 updated.
