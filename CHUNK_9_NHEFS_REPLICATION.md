# Chunk 9 — NHEFS Ch. 17 replication + survival vignette

> **Status: ✅ Done**
> **Depends on:** Chunks 2–7 (the full point-treatment g-computation surface; IPW for the weighted
> arm; CR optional).
> **Oracle:** Hernán & Robins *Causal Inference: What If* Ch. 17 published
> numbers; `survival::survfit()` (unadjusted KM sanity).

## Goal

The headline acceptance test for point-treatment g-computation: reproduce the Hernán & Robins Ch. 17
NHEFS analysis, and ship a `survival` vignette that walks the two-step API end
to end. This is the package's "it actually works on a real, published analysis"
proof.

## Replication targets (handoff §9)

- 120-month survival: **≈ 80.7%** under treatment (`qsmk = 1`), **≈ 80.5%**
  under no treatment (`qsmk = 0`).
- Risk difference at 120 months: **≈ 0.2%** (95% CI **−4.1% to 3.7%**) —
  essentially null.

## Data note (handoff §5)

The NHEFS bundled in `causatr` is **stripped** and lacks the `death` / `yrdth`
columns this analysis needs. Chunk 9 must source the **full** NHEFS extract
(with survival time + death indicator), e.g. via a `data-raw/` script from the
`causaldata` package or the H&R course distribution, reshape to person-period
(monthly intervals to month 120) with `causatr::to_person_period()`, and either
bundle a survival-augmented `nhefs_surv` dataset or build it in the test/vignette
setup.

## Deliverables

### New files
- `data-raw/nhefs_survival.R` — source + assemble the survival-augmented NHEFS;
  document provenance + license.
- `data/nhefs_surv.rda` (or build in-test if licensing precludes bundling).
- `tests/testthat/test-nhefs-replication.R` — pin the Ch. 17 targets.
- `vignettes/survival.qmd` — the two-step API walkthrough on NHEFS (gcomp curve,
  risk difference, RMST; sandwich + bootstrap CIs; a plot + a forest at t* =
  120).

### Updated files
- `DESCRIPTION` — `survival` (KM sanity) + `quarto` in `Suggests` if not present.

## Tests

- `type = "survival"` gcomp curve at t = 120 within tolerance of 80.7% / 80.5%.
- `risk_difference` at t = 120 ≈ 0.2%, sandwich CI brackets the published
  (−4.1%, 3.7%) — assert the point within ~1% and the CI width in the right
  ballpark (the published CI is the acceptance band, not an exact pin).
- Unadjusted KM via `survival::survfit()` as a sanity cross-check on the
  marginal curve.

## Behaviour rules
- Pin the **point estimate** tightly and the CI as a **band** (sandwich SE need
  not match the book's bootstrap CI exactly — document the comparison basis).
- Monthly person-period grid to 120; censor at administrative end.

## Non-goals (deferred)
- Longitudinal ICE-hazard / longitudinal NHEFS (the Ch. 17 analysis is point-treatment).

## Dependencies & composition
- Chunks 2–4 are the minimum; Chunk 5 (IPW) adds the weighted arm Ch. 17 also
  reports (IP-weighted survival curves, Cole & Hernán 2004).

## Acceptance checklist
- [ ] Survival-augmented NHEFS sourced + documented (`data-raw/`).
- [ ] 120-month survival ≈ 80.7% / 80.5%; RD ≈ 0.2% within tolerance.
- [ ] Sandwich CI brackets the published (−4.1%, 3.7%).
- [ ] `vignettes/survival.qmd` builds and renders curve + forest.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
