# Chunk 1 — Skeleton + Track A fit path

> **Status: ✅ done** — landed in `6e911d3`
> (`feat(chunk-1): Track A fit skeleton -- surv_fit() pooled-logistic hazard`).
> 79 tests passing. `devtools::check()`: 0 errors / 0 warnings.

> Per-chunk implementation guide for `survatr` chunk 1. Reads alongside
> [SURVIVAL_PACKAGE_HANDOFF.md](SURVIVAL_PACKAGE_HANDOFF.md) (scope, design,
> invariants) and [.claude/hard-rules.md](.claude/hard-rules.md) (project
> conventions).

## Goal

Ship a fit-only Track A entry point: `surv_fit()` takes person-period (PP)
data + outcome/treatment/confounders/id/time/(optional censoring) +
`time_formula` + `model_fn`, builds the risk set, fits the pooled-logistic
hazard model, and returns a `survatr_fit` S3 object. **No contrast, no
variance, no curves** — those ship in later chunks.

This chunk is a direct port-and-generalize of the pre-removal
`causatr::causat_survival()` fit path, lifted into its own package with
the naming / API tightened (`estimator`, `model_fn`, `survatr_*` error
classes and reserved columns).

## Deliverables

### New R files

| File | Contents |
|---|---|
| `R/surv_fit.R` | **Exported** `surv_fit()` — main fit entry. Validates inputs, builds risk set, calls the hazard engine, returns `survatr_fit`. |
| `R/gcomp_survival.R` | Internal `fit_hazard_gcomp()` — pooled-logistic fit on at-risk PP rows. Family switch: `binomial()` unweighted, `quasibinomial()` weighted. |
| `R/prepare_data.R` | Internal `prepare_pp_data()` — coerce to `data.table`, column-presence checks, PP-shape check (every `id` has ≥2 rows, unique `(id, time)` pairs). |
| `R/risk_set.R` | Internal `build_risk_set()` — within-id lagged cumsums `.survatr_prev_event` / `.survatr_prev_cens`, returns the logical fit-row mask. |
| `R/checks.R` | Copy-adapt from causatr: `check_weights()`, `check_dots_na_action()`, `check_reserved_cols()`. Error classes `survatr_bad_weights`, `survatr_bad_na_action`, `survatr_reserved_col`. |
| `R/utils.R` | Copy from causatr: `is_uncensored()` (NA|0 ⇒ TRUE). S3 constructor `new_survatr_fit()`. Constant `SURVATR_INTERNAL_COLS = c(".survatr_prev_event", ".survatr_prev_cens")`. |
| `R/survatr-package.R` | Add `#' @importFrom` tags for `data.table`, `rlang`, `stats` as needed. |
| `R/print.R` | `print.survatr_fit()` — minimal banner (estimator, outcome, treatment, n individuals, n PP rows, time grid span). Placeholder; polished in chunk 4. |

### Tests (`tests/testthat/`)

| File | Asserts |
|---|---|
| `test-surv_fit.R` | happy path on NHEFS-like PP data: returns `survatr_fit` with `track = "A"`, `estimator = "gcomp"`, `$model` is `glm` with `family$family == "binomial"`, `$pp_data` has no `.survatr_*` columns, `$time_grid` is sorted unique times. |
| `test-prepare_data.R` | rejects wide data (`survatr_not_person_period`); rejects duplicated `(id, time)` (`survatr_duplicate_pp_row`); coerces `data.frame` to `data.table`; copies (does not mutate) input. |
| `test-risk_set.R` | drops rows at/after first event per id; drops rows at/after first censor per id; `is_uncensored()` treats `NA | 0` as at-risk. |
| `test-checks.R` | `check_weights()` rejects NA/Inf/NaN/negative/non-numeric/mis-sized; accepts zero weights. `check_dots_na_action()` rejects `na.exclude`, accepts `na.omit` / `na.fail` / absent. `check_reserved_cols()` rejects `.survatr_prev_event` / `.survatr_prev_cens` column collisions. All use `expect_snapshot(error = TRUE)`. |
| `test-surv_fit-weighted.R` | passing `weights = w` switches family to `quasibinomial()`, fits on `fit_rows` subset, coefficients match an unweighted fit when `w ≡ 1`. |
| `test-surv_fit-family-oracle.R` | **Truth-based sanity**: fit `surv_fit()` on a closed-form DGP (constant hazard, no covariates, `time_formula = ~ 1`) and check the recovered β₀ ≈ `qlogis(h_true)` within MC tolerance on n=5000. |

### Docs / infra

- Roxygen on every function (exported + internal). `@noRd` on internals.
- Run `devtools::document()`; inspect `NAMESPACE` for the expected exports
  (`surv_fit`, `print.survatr_fit`) and nothing else.
- `air format .` clean.
- Update [FEATURE_COVERAGE_MATRIX.md](FEATURE_COVERAGE_MATRIX.md) — no
  user-facing feature ships yet, so add a short "Chunk 1 (fit-only): ✓
  `surv_fit()` binomial / quasibinomial — smoke + closed-form β₀ oracle"
  line under a **Track A** section.

## API contract

```r
surv_fit(
  data,                              # person-period data.table / data.frame
  outcome,                           # event-indicator column name (0/1)
  treatment,                         # treatment column name (point-treatment: constant within id)
  confounders,                       # RHS formula, e.g. ~ L1 + L2
  id,                                # individual id column
  time,                              # time column (integer-valued period index)
  censoring = NULL,                  # optional censoring indicator column
  time_formula = ~ splines::ns(time, 4),  # alpha(t) baseline hazard RHS
  weights = NULL,                    # external weights on PP rows
  estimator = "gcomp",               # Track A only in chunk 1; "ipw"/"ice" deferred
  model_fn = stats::glm,             # fitting function
  ...                                # forwarded to model_fn; na.exclude rejected
)
# => survatr_fit: list(model, pp_data, treatment, outcome, confounders,
#                      id, time, censoring, time_grid, track = "A",
#                      estimator = "gcomp", family = "binomial"|"quasibinomial",
#                      call)
```

### Behaviour rules (non-negotiable — see hard-rules.md)

- **Pooled-logistic only** this chunk. Discrete-time hazard
  `logit h(t | A, L) = alpha(t) + beta_A A + beta_L L`.
- **Risk set = `.survatr_prev_event == 0 & (.survatr_prev_cens == 0 &
  is_uncensored(data, censoring))?`**. Lagged within-id cumsums with
  `data.table::shift(..., n = 1L, fill = 0, type = "lag")`.
- **Family switch is load-bearing.** Unweighted ⇒ `binomial()`; weighted ⇒
  `quasibinomial()` (drops "non-integer #successes" warning, preserves
  score equations).
- **`na.action = na.exclude` is a hard reject.** Error class
  `survatr_bad_na_action`.
- **Reserved columns: `.survatr_prev_event`, `.survatr_prev_cens`.** Guard
  upfront with `check_reserved_cols()`; strip from `fit$pp_data` before
  return.
- **Matching: hard reject** with class `survatr_matching_rejected` (full
  rejection path lives in chunk 8, but `estimator = "matching"` should
  abort here too rather than fall through).
- **`competing` not wired this chunk** — adding the argument but erroring
  with `survatr_competing_misuse` if passed non-`NULL` is fine, since full
  CR lives in chunk 7. Mirrors causatr's pre-removal guard.

## Non-goals (deferred)

- `contrast.survatr_fit()` — chunk 2.
- Sandwich variance / IF machinery — chunk 3.
- Bootstrap + `plot` / `tidy` / `forrest` — chunk 4.
- IPW weighted MSM — chunk 5 (family switch in place but no density-ratio
  weight construction yet).
- Track B (ICE) — chunk 6.
- Full competing risks — chunk 7.
- Full matching rejection path + tests — chunk 8.
- NHEFS Ch. 17 replication — chunk 9 (this chunk can smoke-test on a
  synthetic closed-form DGP; NHEFS replication is acceptance for Track A
  end-to-end).

## Seed sources

- Pre-removal `causatr::causat_survival()` — still present in the
  installed `causatr` (local version 0.0.0.9000). Inspect with
  `View(causatr::causat_survival)` or `deparse(causatr::causat_survival)`.
  **Rename** `causatr_*` → `survatr_*` for internal columns, error
  messages, and the S3 constructor; **retarget** `model_fn` from hardcoded
  `stats::glm`; **tighten** the `estimator` / `track` surface.
- Internals to reuse via `causatr:::` (do NOT re-implement):
  `causatr:::check_weights`, `causatr:::check_dots_na_action`,
  `causatr:::check_reserved_cols`, `causatr:::is_uncensored`. For chunk 1
  we **copy** these into `survatr` (the handoff §5 says "Copy") because
  the reserved-col list is survatr-specific and the error classes need to
  be `survatr_*`. Later chunks will reach into `causatr:::` for IF
  primitives only.

## Acceptance checklist

- [ ] `devtools::load_all()` succeeds.
- [ ] `devtools::document()` regenerates `NAMESPACE` cleanly; only
      `surv_fit` and the `print.survatr_fit` S3 method are exported.
- [ ] `devtools::test()` — all tests in the files above pass; error-path
      tests use `expect_snapshot(error = TRUE)`.
- [ ] `devtools::check()` — 0 errors, 0 warnings; the usual "new
      submission" note is acceptable.
- [ ] `air format .` is a no-op.
- [ ] `Rscript -e 'lintr::lint_package()'` — no new lints.
- [ ] [FEATURE_COVERAGE_MATRIX.md](FEATURE_COVERAGE_MATRIX.md) updated.
- [ ] Single commit (or small series) with message `feat(chunk-1): Track A
      fit skeleton — surv_fit() pooled-logistic hazard on PP data`.
