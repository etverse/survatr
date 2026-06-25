# Chunk 17 — Target-trial alignment: landmark + immortal-time tooling

> **Status: ⬜ Not started**
> **Depends on:** Chunk 2 (fit + curve), Chunk 10 (`diagnose()` panel host);
> Chunk 14 (delayed entry) for the data-side mechanics of re-anchoring.
> **Oracle:** n/a (diagnostic + pedagogical). Validate by **structure** + a
> known-bias DGP where the immortal-time check fires (treatment defined using
> post-baseline survival ⇒ flag), and a clean DGP where it stays silent.

## Goal

Tooling to keep survatr analyses aligned with the **target-trial emulation**
framing — the discipline that prevents the most common self-inflicted survival
bias (immortal time from misaligned time zero). Three deliverables, mostly
diagnostics + docs, no new estimator:

1. **Landmark analysis** — re-anchor time zero at a user-chosen landmark,
   conditioning on survival to it, and estimate the curve from there.
2. **A `diagnose()` time-zero-alignment check** — flag when treatment is (or
   appears to be) defined using post-baseline information, the immortal-time
   signature.
3. **A target-trial-emulation vignette** — the framing that ties the eligibility
   / treatment-assignment / time-zero triad to survatr's API.

## Design

### Landmark analysis

A landmark analysis fixes a landmark time `t_L`, **restricts to individuals who
survived event-free and uncensored to `t_L`**, re-anchors time zero at `t_L`,
and estimates `S^a(t | survived to t_L)` from there. Mechanically this reuses
the chunk-14 delayed-entry plumbing in reverse: instead of variable entry, every
retained id enters at `t_L` and the pre-landmark rows are dropped from both the
fit and the curve. The estimand is explicitly **conditional on survival to the
landmark** — the same conditioning caveat documented for delayed entry. Treatment
status is the value **as of the landmark** (re-defining treatment at `t_L` is the
whole point — it is the analysis that lets a post-baseline-defined exposure be
assigned at a legitimate time origin without manufacturing immortal time).

### Immortal-time / time-zero-alignment check

A `diagnose()` panel that inspects the alignment of three events that a target
trial requires to coincide at time zero: **eligibility**, **treatment
assignment**, and **start of follow-up**. The detectable signature of
immortal-time bias is treatment status that is **determined by post-baseline
survival** — operationally, the check flags when:

- the treatment column varies **within** an id **after** the first at-risk
  period in a point-treatment g-computation fit (treatment should be fixed at
  baseline for the point-treatment path — a post-baseline change is the immortal-time smell);
- treated ids have a **structural gap** between time zero and first observed
  event/censor that untreated ids lack (a survivorship requirement built into
  the exposure definition);
- (informational) the share of person-time that is "guaranteed event-free by
  construction" — the immortal person-time.

The check is **non-fatal** (the chunk-10 rule): it reports + flags, it never
aborts a fit. It surfaces a pointer to the landmark analysis and the
target-trial vignette as the remedy.

### Target-trial vignette

A pedagogical vignette mapping the target-trial protocol (eligibility,
treatment strategies, assignment, time zero, outcome, estimand, analysis plan)
onto survatr's two-step API, with a worked immortal-time example: the naive
misaligned fit, the diagnostic flag firing, and the landmark / re-anchored fit
that fixes it.

## Deliverables

### New R files
- `R/landmark.R` — landmark re-anchoring: validate `t_L` on the grid, restrict
  to ids event-free + uncensored through `t_L`, re-anchor time zero at `t_L`
  (reuse the chunk-14 at-risk-window mechanics), take treatment as of `t_L`, and
  hand the restricted PP table to the chunk-2 fit/contrast spine. The
  conditional-on-landmark estimand is recorded on the fit.

### Updated R files
- `R/surv_fit.R` — accept `landmark = <time>` (default `NULL` → no landmark);
  thread it into the restricted-PP construction.
- `R/diagnose.R` (chunk 10) — add the **time-zero-alignment** panel:
  per-arm structural-gap summary, within-id post-baseline treatment-change flag
  (point-treatment), immortal-person-time share; non-fatal flags with a remedy pointer.
- `R/print.R` — render the alignment panel in `print.survatr_diag`.
- `vignettes/` — `target-trial.qmd` (the pedagogical vignette).

### Tests (`tests/testthat/`)
- `test-landmark.R` — landmark re-anchoring restricts to survivors of `t_L`,
  re-anchors time zero, and the curve matches a direct fit on the manually
  re-anchored data; `landmark = NULL` reproduces the chunk-2 curve; invalid
  `t_L` (off-grid, beyond last event) aborts (`survatr_bad_landmark`).
- `test-immortal-time-check.R` — **known-bias DGP** where treatment is defined
  using post-baseline survival: the alignment check **fires** (flag set);
  on a clean baseline-defined-treatment DGP the check stays **silent**;
  structure of the alignment panel; non-fatal (never aborts).

## API contract

```r
# Landmark analysis: re-anchor time zero at t_L, condition on survival to it
fit_lm <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "t", censoring = "C",
  landmark = 12,                  # NEW — re-anchor at t_L = 12
  estimator = "gcomp"
)
result <- contrast(fit_lm, interventions = list(a1 = causatr::static(1),
                                                a0 = causatr::static(0)),
                   times = ..., type = "risk_difference")
# estimand: S^a(t | survived event-free & uncensored to t_L)

# Time-zero-alignment / immortal-time diagnostic
dx <- diagnose(fit)
dx$alignment   # per-arm structural-gap summary, post-baseline-treatment flag,
               # immortal-person-time share, remedy pointer
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Landmark conditions on survival to `t_L`** and re-anchors time zero there;
  the estimand is explicitly conditional-on-landmark (same caveat as chunk-14
  delayed entry). Document it in roxygen.
- **Landmark reuses the chunk-14 at-risk-window mechanics** — do not fork a
  second risk-set builder; restrict + re-anchor on the existing one.
- **Treatment is taken as of the landmark** for a landmark fit (the legitimate
  re-definition of exposure at a valid time origin).
- **The alignment check is a `diagnose()` panel and is non-fatal** — report +
  flag, never abort (the chunk-10 rule). It is a *smell test*, not a proof of
  bias; word the flag accordingly and point to the landmark remedy + vignette.
- **`landmark = NULL` reproduces the chunk-2 curve exactly** (degenerate-to-
  current invariant).
- Invalid `t_L` (off-grid / beyond last observed event) aborts with
  `survatr_bad_landmark`.

## Non-goals (deferred)
- **Sequential / rolling landmark (super-landmark) analysis** — a single
  landmark here; the rolling extension is a future lift.
- **Cloning–censoring–weighting** for time-varying eligibility / grace periods —
  out of scope for this chunk (it is a Track-B + IPCW composition, far heavier).
- **Automated bias *correction*** — the check diagnoses; the remedy is the
  landmark / re-anchored design the user chooses, not a silent adjustment.

## Dependencies & composition
- Chunk 2 (fit + curve), Chunk 10 (the `diagnose()` host for the alignment
  panel), Chunk 14 (delayed-entry at-risk-window mechanics reused by the
  landmark restriction). No new causatr internals; the landmark fit flows
  through the existing chunk-2 spine.

## Acceptance checklist
- [ ] `landmark =` restricts to survivors of `t_L`, re-anchors time zero, and
      matches a manual re-anchored fit; `landmark = NULL` reproduces chunk 2.
- [ ] `survatr_bad_landmark` fires on off-grid / out-of-range `t_L`.
- [ ] The `diagnose()` alignment panel fires on a known immortal-time DGP and
      stays silent on a clean one; non-fatal.
- [ ] Target-trial vignette ships with the worked naive → flag → landmark
      example.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
