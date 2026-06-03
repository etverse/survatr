# Chunk 14 — Left-truncation / delayed entry

> **Status: ⬜ Not started**
> **Depends on:** Chunk 1 (risk-set builder, `prepare_pp_data()` guards),
> Chunk 2 (contrast spine).
> **Oracle:** `survival::survfit(Surv(entry, time, event) ~ 1)` and
> `survival::coxph(Surv(entry, time, event) ~ A)` (delayed-entry counting-
> process form); closed-form on a constant-hazard DGP with staggered entry
> (delayed entry leaves the hazard estimate unbiased, ignoring it does not).

## Goal

Allow individuals to enter the risk set at a **delayed period** rather than at
time zero — left-truncation / delayed entry, ubiquitous in registry / EHR data
and in age-timescale designs (subjects enter the study already older than the
time origin). Pooled-logistic hazard handles this natively: an individual
simply contributes **at-risk rows starting at its entry period**, and the
hazard is estimated only over rows where the individual is observed and at
risk. The work is (1) accepting an entry column / per-id first at-risk period,
(2) relaxing the rectangular-PP requirement to permit a documented delayed
start, and (3) keying the risk-set builder on at-risk rows rather than assuming
entry at the first grid time.

## The math

Standard pooled logistic estimates the hazard from every at-risk row
`k = 1, …, T_i`. Under left-truncation, id `i` enters at period `e_i ≥ 1` and
is **not observed** for `k < e_i`; conditioning on entry, it contributes
at-risk rows only for `e_i ≤ k ≤ T_i`:

```
contribution of id i to the hazard score:  rows  { (i, k) : e_i ≤ k ≤ T_i, at risk }
```

The hazard model is unchanged (`logit h(t | A, L) = alpha(t) + beta_A A +
beta_L L`); only the **fit row set** changes — pre-entry rows are excluded, not
treated as `Y = 0`. This is the discrete-time analogue of the counting-process
risk set `Surv(entry, time, event)` in `survival`: a subject is in the
denominator only while under observation.

The counterfactual survival curve is the cumulative product of fitted hazards.
The estimand choice matters and must be documented:

- **Conditional-on-entry survival** `S^a(t | T > e)` — the curve among those
  who survived to entry, the natural object when the time origin is calendar /
  age and entry is the late-truncation point.
- A marginal curve from time zero is **not identified** for periods before the
  earliest entry — there is no at-risk information there. The contrast path
  anchors the curve at the earliest entry period (or a user-supplied landmark),
  and `S^a` is reported on the support where at-risk rows exist.

Variance is unchanged in form: the chunk-3 cross-time delta chain runs over the
**observed** at-risk rows per individual; pre-entry rows simply do not enter the
per-row IF sum. The `n × |t-grid|` IF matrix has structural zeros in an
individual's pre-entry columns.

## Design

- A new `entry =` argument on `surv_fit()`: either a column giving the per-id
  first at-risk period (delayed-entry time), or `NULL` (current behaviour,
  everyone enters at the first grid time).
- The risk-set builder keys on **at-risk rows** (entry ≤ k ≤ event/censor),
  rather than assuming the at-risk window starts at the first grid time. The
  existing `.survatr_prev_event` / `.survatr_prev_cens` cumsum machinery is
  augmented with an entry mask so pre-entry rows are dropped from the fit, not
  scored as `Y = 0`.
- The **rectangular-PP guard** (`survatr_ragged_pp`) is relaxed to allow a
  documented delayed start: ids need not have rows for `k < e_i`. The guard
  still rejects *interior* gaps and post-event/post-censor rows — only a
  leading delayed-entry gap is now permitted, and only when `entry =` is
  supplied. Without `entry =`, the rectangular requirement stands unchanged.

## Deliverables

### New R files
- `R/left_truncation.R` — entry-column validation, per-id entry period
  resolution, the at-risk mask (entry ≤ k ≤ terminal), and the relaxed
  rectangular check that distinguishes a permitted leading gap from a forbidden
  interior gap.

### Updated R files
- `R/surv_fit.R` — accept `entry =` (column name or `NULL`); thread the entry
  period into the risk-set construction.
- `R/risk_set.R` — exclude pre-entry rows from the fit row set (drop, do not
  set `Y = 0`); the cumsum risk-set bookkeeping keys on at-risk rows.
- `R/prepare_data.R` — relax `survatr_ragged_pp` to allow a leading delayed-
  entry gap **only** when `entry =` is supplied; validate entry ∈ grid, entry ≤
  terminal period, no NA (`survatr_bad_entry`).
- `R/survival_curve.R` — anchor the curve at the earliest entry period; report
  `S^a` on the at-risk support; document the conditional-on-entry estimand.

### Tests (`tests/testthat/`)
- `test-left-truncation.R` — staggered-entry constant-hazard DGP: the
  delayed-entry hazard estimate is unbiased while the naive
  pre-entry-as-`Y=0` fit is biased; the curve matches
  `survival::survfit(Surv(entry, time, event) ~ 1)` on the at-risk support;
  `entry = NULL` reproduces the chunk-2 curve exactly (degenerate-to-current);
  the interaction with `survatr_ragged_pp` — a permitted leading gap passes, an
  interior gap still aborts.

## API contract

```r
fit <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "t", censoring = "C",
  entry = "enter",                    # NEW — per-id first at-risk period
  time_formula = ~ splines::ns(t, 4),
  estimator = "gcomp"
)
# entry = NULL  -> everyone enters at the first grid time (current behaviour)
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "risk_difference")
# curve is reported on the at-risk support (>= earliest entry); the estimand is
# conditional on survival to entry.
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Pre-entry rows are EXCLUDED from the fit, never scored as `Y = 0`.**
  Treating them as event-free observation manufactures person-time that was
  never at risk and biases the hazard downward.
- **The hazard model is unchanged** — only the fit row set (at-risk window)
  changes. Discrete-time analogue of `Surv(entry, time, event)`.
- **The rectangular-PP guard is relaxed for a LEADING gap only**, and only when
  `entry =` is supplied. Interior gaps and post-event/post-censor rows still
  abort (`survatr_ragged_pp`). Without `entry =`, the chunk-1 rectangular
  contract stands.
- **The reported curve is conditional on survival to entry**; a marginal curve
  before the earliest entry is not identified. Anchor at the earliest entry
  period (or a user landmark) and document the estimand in roxygen.
- **`entry = NULL` reproduces the chunk-2 curve exactly** (degenerate-to-
  current invariant).
- Entry must lie on the time grid, be ≤ the terminal period, and be non-NA
  (`survatr_bad_entry`).

## Non-goals (deferred)
- **Interval censoring** (event known only to lie within an interval) — a
  different likelihood; out of scope for this chunk.
- **Time-dependent / dynamic entry** (re-entry after a gap) — only a single
  delayed entry per id here.
- **Landmark analysis** (re-anchoring time zero conditional on survival to a
  landmark) — related but distinct; that is Chunk 17. This chunk handles the
  data-side delayed entry; the landmark chunk handles the estimand-side
  re-anchoring.

## Dependencies & composition
- Chunks 1, 2. Composes with the variance chunks (3, 13) without change: the
  delta chain runs over observed at-risk rows; pre-entry columns of the
  `n × |t-grid|` IF matrix are structurally zero. No new causatr internals.

## Acceptance checklist
- [ ] `entry =` excludes pre-entry rows from the fit (not scored `Y = 0`);
      delayed-entry hazard estimate is unbiased while the naive fit is biased.
- [ ] Curve matches `survival::survfit(Surv(entry, time, event) ~ 1)` on the
      at-risk support.
- [ ] `entry = NULL` reproduces the chunk-2 curve exactly.
- [ ] Relaxed `survatr_ragged_pp` permits a leading gap, still aborts interior
      gaps; `survatr_bad_entry` validation fires.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
