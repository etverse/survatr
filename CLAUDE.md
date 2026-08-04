# survatr

Causal survival analysis via pooled-logistic hazard g-computation, IPW, and iterated
conditional expectation (ICE) on person-period data. Extends the
[`causatr`](https://github.com/etverse/causatr) engine to time-to-event outcomes with
delta-method sandwich variance on cumulative-product survival curves, risk / RMST
contrasts, competing-risks cause-specific hazards, and empirical bootstrap. Part of
the [etverse](https://github.com/etverse) ecosystem.

## Guide files

- `SURVIVAL_PACKAGE_HANDOFF.md` — **single source of truth** for scope, design,
  integration with causatr, and the §10 chunk roadmap. Read this first.
- `CHUNK_<N>_<name>.md` at the repo root — one design doc per chunk, carrying the
  rationale and rejected alternatives for that piece.
- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works"**;
  per-combination feature × test coverage. Every PR that changes a feature MUST
  update this file and the corresponding tests.
- `.claude/hard-rules.md` — architecture invariants review must not re-flag as bugs.

## Project structure

This is an R package: `R/` (source, layout below), `tests/testthat/` (tests,
`test-foo.R` mirrors `R/foo.R`), `man/` and `NAMESPACE` (generated — do not edit),
`DESCRIPTION` (metadata and dependencies), `vignettes/` (long-form docs).

## R/ file layout

Filenames name their contents; read the directory for the current set. The groupings:

- **Core API** — `surv_fit.R` (pooled-logistic hazard on person-period data),
  `contrast.R` (time-indexed curve / risk / RMST contrasts), `diagnose.R`
- **Estimation engines** — `gcomp_survival.R`, `ipw_survival.R`, `ice_survival.R`,
  `competing_risks.R` (cause-specific hazards + CIF)
- **Inference** — `variance_if_survival.R` (delta-method cross-time IF aggregation),
  `variance_cluster.R`, `variance_bootstrap.R`
- **Data utilities** — `risk_set.R`, `prepare_data.R`
- **S3 methods** — `print.R`, `summary.R`, `plot.R`, `tidy.R`
- **Helpers** — `utils.R`, `checks.R`, `zzz.R`, `survatr-package.R`

## S3 classes

- `survatr_fit` (`surv_fit()`) — hazard model(s), person-period data, track (A / B)
- `survatr_result` (`contrast()`) — time-indexed `estimates` / `contrasts`
  `data.table`s, `time_grid`, and the vcov as an `n × |t-grid|` IF matrix →
  `|t-grid| × |t-grid|` covariance
- `survatr_diag` (`diagnose()`) — per-period positivity and balance,
  competing-risks decomposition summaries

## Two-step API

```r
# Step 1: fit the hazard model once
fit <- surv_fit(
  data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  id = "id", time = "time", censoring = "C",
  time_formula = ~ ns(time, df = 4),  # alpha(t) baseline hazard
  estimator = "gcomp",                # or "ipw" / "ice"
  model_fn = stats::glm               # defaults to glm with binomial()
)

# Step 2: contrast many interventions, curve-shaped result
result <- contrast(
  fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times = seq(0, 120, by = 12),
  type = "risk_difference",           # or "risk_ratio" / "rmst_difference"
  ci_method = "sandwich"              # or "bootstrap"
)
```

## Development commands

```r
devtools::load_all()     # load for dev
devtools::test()         # run tests
devtools::check()        # R CMD check
devtools::document()     # regenerate roxygen
```

Shell: `air format .` (format all R files).

## Code style

- Use `pkg::fun()` for external package calls (never bare `library()` in `R/` code)
- Use `rlang::abort()` / `rlang::warn()` / `rlang::inform()` instead of `stop()` / `warning()` / `message()`
- Use `data.table` internally (`:=`, `.SD`, `.N`, `by`); return `data.table` from user-facing functions
- **Roxygen on every function, including internal/`@noRd` helpers.** Document params, return value, and a short `@description`.
- **Add inline comments generously.** Explain non-obvious logic: math derivations, design rationale, subtle invariants, tricky indexing, `data.table` semantics, and any "why this instead of the obvious approach" decisions. Obvious one-liners don't need comments; math-heavy or non-local logic does.
- Do NOT remove existing comments unless the related code is also removed; place
  exported functions at top of files, internal helpers below

## Testing rules

- Use `expect_snapshot(error = TRUE)` for testing error conditions
- NEVER delete or mock away failing tests — fix the source code; no complex mocking
- **Truth-based simulation tests are mandatory for new features.** Every supported (track × estimator × treatment type × intervention × variance) combination MUST assert against an analytical or external reference (`lmtp::lmtp_tmle(outcome_type = "survival")` for point estimates; `survival::survfit` / `survival::coxph` for unadjusted KM / Cox sanity; `gfoRmula::gformula_survival()` for the longitudinal ICE path). Smoke tests accepted only where no oracle exists and should be marked 🟡 in `FEATURE_COVERAGE_MATRIX.md`.
- **NHEFS Ch. 17 replication** (handoff §9) is the acceptance target for point-treatment g-computation: 120-month survival ≈ 80.7% under treatment / 80.5% under no treatment; risk difference ≈ 0.2% (95% CI: −4.1% to 3.7%).
- **External reference cross-checks.** When the analytical truth is hard to derive (longitudinal MTPs, competing risks under non-static interventions), validate against `lmtp::lmtp_tmle(outcome_type = "survival")` once at test design time, then pin the expected value with a comment block citing the validation.

## Constraints

- Run `devtools::test()` before committing
- Do not modify `man/` or `NAMESPACE` directly — use roxygen2 tags
- Run `devtools::document()` after changing roxygen comments

## Scope

survatr is a **causal survival analysis package** for time-to-event outcomes, built
on causatr. It owns pooled-logistic hazard g-computation, ICE hazards for
longitudinal survival, IPW weighted hazard MSM, built-in IPCW, **parametric**
doubly-robust (AIPW) survival, cause-specific hazards + CIF for competing risks,
recurrent-event / multi-state extensions, left-truncation, and a curve-valued
estimand surface with sandwich and bootstrap variance. Phased roadmap:
[SURVIVAL_PACKAGE_HANDOFF.md §10](SURVIVAL_PACKAGE_HANDOFF.md).

It is NOT: a Cox-PH modelling package (use `survival`), an **ML / TMLE** survival
package (use `lmtp` — survatr's AIPW is parametric / M-estimation only), a
forward-simulation g-formula (use `gfoRmula`), or a Fine–Gray
subdistribution-hazards package (competing risks are cause-specific only).

### Why a separate package

Six reasons, in full at `SURVIVAL_PACKAGE_HANDOFF.md` §1: person-period runtime,
a curve-shaped estimand that the `causatr_result` S3 shape cannot carry,
first-class competing risks, a discrete-time hazard outcome model, cross-time
delta-chain variance aggregation, and the matching rejection.

### Why ICE, not forward simulation

See `CHUNK_6_ICE_B.md` §Goal. In short: ICE models the outcome only at each time
and admits a stacked estimating-equation sandwich; forward simulation
(`gfoRmula`-style) needs models for every time-varying covariate and is
bootstrap-only.

## R ecosystem integration

`DESCRIPTION` is authoritative for the dependency tiers. The non-obvious relationships:

| Package | Relationship |
|---|---|
| `causatr` | Imports (GitHub remote `etverse/causatr`) — interventions, IF primitives, `to_person_period`, NHEFS |
| `lmtp` | Suggests, **point-estimate oracle only** — its EIF-SE is not comparable to survatr's |
| `gfoRmula` | Suggests, oracle for the longitudinal ICE-hazard track |
| `survival` | Suggests, **test oracles only** (unadjusted KM / Cox sanity) |
| `MatchIt` | Suggests, exercises the rejection path only |
| `mgcv` | Suggests, GAM time-spline baseline hazard via `model_fn` |

## Supported features

Combination-audit table: `.claude/hard-rules.md` §"Supported dimensions". Per-combination
pass/fail status: `FEATURE_COVERAGE_MATRIX.md`. Neither is restated here.

Boundaries worth knowing without opening either: two tracks (A point-treatment
pooled-logistic hazard, B longitudinal ICE hazards); estimators gcomp / ipw / ice / aipw
with **matching hard-rejected** and ML/TMLE out of scope; competing risks via
cause-specific hazards + CIF only, **Fine–Gray out of scope**.

## Key design decisions

The mathematical invariants behind these (cumulative-product-before-averaging /
Jensen, the delta-method cross-time IF aggregation, the cluster-robust `n²` divisor,
the GAM lpmatrix basis, the reserved `.survatr_*` prefix) are stated with their
numerical validation in `.claude/hard-rules.md`; the package-split rationale is in
`SURVIVAL_PACKAGE_HANDOFF.md` §1. What follows is only the API shape.

- **causatr is the engine.** All IF primitives (`prepare_model_if()`,
  `apply_model_correction()`, `vcov_from_if()`, numeric fallback),
  intervention constructors (`static()`, `shift()`, ...), and the
  `to_person_period()` reshape live in causatr. survatr layers the
  cross-time aggregation and survival-curve estimand shape on top.
- **`model_fn`, not hardcoded link.** User passes `stats::glm`,
  `mgcv::gam`, etc. Defaults to `stats::glm` with `family = binomial()`
  (or `quasibinomial()` for weighted fits).
- **`estimator`, not `method`.** Same rule as causatr — `method` is
  reserved for `...`-forwarding (though matching is rejected here).
- **Two-step API.** Fit the hazard model once; contrast many interventions.
  Survival curves live in the `contrast()` return, not the fit.
- **Curve-shaped `survatr_result`.** Time-indexed `data.table`s for
  `estimates` and `contrasts`. S3 methods (`print`, `plot`, `tidy`)
  dispatch on the time-indexed shape. `forrest`-style forest plots at a
  user-chosen reference time `t*`.
- **Matching + survival rejected.** `survival::coxph(..., weights =
  match_weights, cluster = subclass)` on the `MatchIt` output is the
  correct tool and lives outside survatr. Hard-abort with error class
  `survatr_matching_rejected`.
## Implementation chunks

**Canonical roadmap, per-chunk status and commit pins:**
[SURVIVAL_PACKAGE_HANDOFF.md §10](SURVIVAL_PACKAGE_HANDOFF.md) — 29 chunks, one
`CHUNK_<N>_*.md` doc each. That table is the only copy; update it there when a chunk
flips status and do not restate per-chunk status here. Done: 1–13 (v1 complete, v1.x
in progress). Next: 14 (left-truncation / delayed entry).

## Architecture notes

Hard-won "do NOT flag this as a bug" invariants live in
[`.claude/hard-rules.md`](.claude/hard-rules.md), which `implement-feature` and
`critical-review-loop` read at Step 0. Each entry records the invariant, the
numerical validation that established it, and what reviewers must not re-flag.
That file is their home — extend it rather than re-adding narrative here.
