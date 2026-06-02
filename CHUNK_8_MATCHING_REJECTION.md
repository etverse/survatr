# Chunk 8 — Matching + survival rejection surface

> **Status: ⬜ Not started** (partially wired in Chunk 1)
> **Depends on:** — (independent).
> **Oracle:** n/a (rejection path; `MatchIt` in `Suggests` only for the test).

## Goal

Make matching + survival a clean, fully-covered **hard rejection**. `surv_fit()`
already aborts `estimator ∈ {"matching", "match"}` with
`survatr_matching_rejected` (wired in Chunk 1). Chunk 8 closes the surface: any
matching entry point (including a `MatchIt` object handed in as data, or a
`method = "matching"` style mis-call) aborts with the same classed error, and a
regression test pins the message + the redirect to the correct tool.

The rationale (handoff §1.6, §3): `MatchIt` weights + pooled logistic on
person-period data is awkward and duplicative —
`survival::coxph(..., weights = match_weights, cluster = subclass)` directly on
the `MatchIt` output is the right tool and lives outside survatr.

## Deliverables

### Updated R files
- `R/surv_fit.R` — ensure every matching alias and entry route hits the single
  upstream abort (audit: `estimator = "matching"/"match"`, a `matchit` object
  passed as `data`, `method =` mis-forwarding). One classed error, one message.

### Tests (`tests/testthat/`)
- `test-matching-rejection.R` — `skip_if_not_installed("MatchIt")` for the
  object-passing case; `expect_error(..., class = "survatr_matching_rejected")`
  for each route; snapshot the message (must name
  `survival::coxph(..., weights = , cluster = )`).

## API contract

```r
surv_fit(..., estimator = "matching")
#> Error: Matching + survival is out of scope for survatr.
#>   i Use survival::coxph(..., weights = match_weights, cluster = subclass)
#>     directly on the MatchIt output.
#  class: survatr_matching_rejected
```

## Behaviour rules (non-negotiable — see hard-rules.md)
- Single unified abort upstream, error class `survatr_matching_rejected`,
  regardless of treatment type. Message points to `survival::coxph(weights=,
  cluster=)`.
- Matching is binary-only in causatr; survival + matching is rejected
  regardless of treatment type.

## Non-goals (deferred)
- Actually supporting matching for survival (permanently out of scope).

## Acceptance checklist
- [ ] Every matching route aborts with `survatr_matching_rejected`.
- [ ] Message names the `coxph` redirect; snapshot-pinned.
- [ ] `FEATURE_COVERAGE_MATRIX.md` rejection row + handoff §10 + CLAUDE.md
      updated.
