# Chunk 25 — Missing data / multiple imputation for survival

> **Status: ⬜ Not started (research first)**
> **Depends on:** Chunk 1 (fit path), Chunk 2 (contrast spine). Composes with
> every estimator (gcomp / IPW / ICE / AIPW) and the variance engines.
> **Oracle:** `smcfcs` (substantive-model-compatible FCS; Bartlett et al. 2015),
> `mice` + Rubin's rules; complete-case and a fully-observed reference as
> sanity bounds.

## Motivation

survatr currently **rejects** NA in predictor columns upfront
(`survatr_na_in_predictors`) — a deliberate boundary that keeps the
influence-function row alignment intact, but it forces users to handle missing
data entirely outside the package. This chunk adds a principled, optional
missing-data path so partially-observed person-period data can be analysed
without ad-hoc preprocessing.

## Research phase (do this before scoping the implementation)

Survey and pin the method choices, then write them into this doc before coding:

1. **Imputation model for person-period / discrete-time hazard data** — how to
   impute baseline confounders, time-varying covariates, and treatment on the
   long grid without breaking the risk-set / cumulative-product structure.
2. **Congeniality with the estimator.** The imputation model must be compatible
   with the substantive g-computation / IPW / ICE model (Meng 1994; Bartlett et
   al. 2015 SMC-FCS) — an uncongenial imputation biases the causal estimate.
3. **Pooling.** Rubin's rules for a **curve-valued** estimand: pool `S^a(t)` (or
   risk / RD / RMST) across imputations per time point, and pool the
   sandwich/bootstrap variance (within- + between-imputation components) across
   the cross-time IF matrix.
4. **Interaction with the existing NA boundary.** Decide whether MI is a
   pre-step that produces `m` completed datasets fed through `surv_fit()` /
   `contrast()` and pooled, or a deeper integration. The former composes
   cleanly with the current architecture and is the likely v1.

## Likely deliverables (refine after the research phase)
- `R/impute.R` — MI driver (wrap `mice` / `smcfcs`), `m` completed PP datasets.
- `R/pool.R` — Rubin's-rules pooling of the time-indexed `survatr_result`
  (estimates + cross-time covariance).
- A `contrast()` / dedicated entry that loops imputations and pools.
- Tests: bias + coverage under MCAR/MAR vs a fully-observed reference and an
  external `smcfcs` oracle; congeniality check.

## Acceptance checklist
- [ ] Research phase written up (method choices pinned with citations).
- [ ] MI path recovers the fully-observed estimate under MCAR/MAR (low bias).
- [ ] Pooled CIs achieve nominal coverage (within- + between-imputation).
- [ ] Congenial-imputation oracle (`smcfcs`) agreement.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.

> Promoted from the handoff "Open research questions" note (2026-06-03) at the
> user's request so the missing-data work is a tracked chunk, not a loose item.
