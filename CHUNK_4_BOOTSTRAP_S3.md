# Chunk 4 -- Track A bootstrap + S3 polish

> **Status: ✅ done** -- landed in the chunk-4 commit.
> `devtools::check()`: 0 errors / 0 warnings (all unused-Imports NOTEs gone
> now that `boot` is used). Full suite pending rerun confirmation.

> Per-chunk implementation guide for `survatr` chunk 4. Reads alongside
> [SURVIVAL_PACKAGE_HANDOFF.md](SURVIVAL_PACKAGE_HANDOFF.md) and
> [.claude/hard-rules.md](.claude/hard-rules.md).

## Goal

Two deliverables in one chunk, because the S3 polish wants to
demonstrate `ci_method = "bootstrap"` alongside `"sandwich"`:

1. **Empirical bootstrap** (`ci_method = "bootstrap"`): resample
   **individuals** (all of each id's PP rows together), refit the hazard
   model per replicate, recompute the per-intervention survival curves /
   contrasts, take sample quantiles or SE-based Wald bands across
   replicates.
2. **S3 polish** on `survatr_result`: `tidy()`, `plot()` (with CI
   ribbons), and `forrest()` (forest plot at a user-chosen reference
   time t*). `print()` already exists from chunk 2 but gets a minor
   polish to show CIs when available.

## Bootstrap design

### Sampling unit

Resample **individuals**, not PP rows. Each replicate draws n_ids ids
with replacement, concatenates the corresponding PP blocks, and refits
the hazard model. Per-row resampling would break the within-id
cumulative-product dependence structure and underestimate variance.

### What we compute per replicate

For each replicate b and each user time t_j:

- Per-intervention `s_hat^a_b(t_j)` (vector of length |interventions|).
- The contrast type's `estimate^a_b(t_j)` (scalar per contrast row, or
  vector if multiple non-reference interventions).

We **store only the aggregated per-time quantities**, not the per-id
IFs -- the bootstrap's role is to calibrate uncertainty via replicate
spread, not to rebuild the IF. A B x |interventions| x |t| array
(plus parallel structure for contrasts) is the data backbone.

### Re-using chunk 2 machinery

One bootstrap replicate = one `surv_fit()` call on the resampled data +
one `contrast(..., ci_method = "none")` call on the refit. Both already
exist. The replicate loop is a thin wrapper that extracts the `s_hat`
/ `estimate` / `rmst_hat` columns from the chunk-2 result.

### CI construction

Two options:

- **Percentile** (default): `ci_lower = quantile(replicates, (1-conf)/2)`,
  `ci_upper = quantile(replicates, (1+conf)/2)`. Simple, robust to
  heavy tails, but can be erratic in small B.
- **Normal / Wald** via bootstrap SE: `ci = point +/- z * sd(replicates)`.
  Matches the sandwich shape but throws away asymmetry.

Go with **percentile** as the default (`boot_ci = "percentile"`), with
`boot_ci = "wald"` as the opt-in alternative. For `risk_ratio` the
percentile path sidesteps the log-scale transform entirely.

### Parallelization

Expose `parallel = c("no", "multicore", "snow")` and `ncpus = 1L` to
match causatr's `contrast()` signature. Under the hood use
`boot::boot()` (already in Imports) for its parallel dispatch,
confidence-interval helpers, and failure tracking -- simpler than
rolling our own.

### Arguments added to `contrast()`

```r
contrast(
  fit,
  interventions,
  times,
  type        = "risk_difference",
  reference   = NULL,
  ci_method   = "none",        # "none" | "sandwich" | "bootstrap"
  conf_level  = 0.95,
  n_boot      = 500L,          # bootstrap replicates (ignored unless bootstrap)
  boot_ci     = "percentile",  # "percentile" | "wald"
  parallel    = "no",          # boot::boot parallel mode
  ncpus       = 1L,            # boot::boot cpu count
  seed        = NULL,          # optional integer for reproducibility
  ...
)
```

`n_boot`, `boot_ci`, `parallel`, `ncpus`, `seed` are all ignored when
`ci_method != "bootstrap"`. They are nonetheless validated upfront so
bad values fail fast.

### Error surface

- `survatr_bad_n_boot` -- `n_boot` not a positive integer.
- `survatr_bad_boot_ci` -- `boot_ci` not in `c("percentile", "wald")`.
- `survatr_boot_failed` -- replicate failure rate exceeds a threshold
  (say 10%); point the user at sandwich or at smaller strata.

## S3 polish

### `tidy.survatr_result()`

Shape: one long `data.table`, columns
`intervention | contrast | time | estimand | estimate | se | ci_lower | ci_upper`.
`intervention` is filled when the row comes from `estimates`, `NA`
otherwise. `contrast` is filled when the row comes from `contrasts`,
`NA` otherwise. `estimand` identifies whether the row is `s_hat` /
`risk_hat` / `rmst_hat` / `risk_difference` / `risk_ratio` /
`rmst_difference`.

Arguments:

```r
tidy(x, which = c("estimates", "contrasts", "all"), conf.int = TRUE, ...)
```

Returns a `data.frame` when called without `data.table` loaded;
structurally identical to what `broom::tidy` consumers expect.

### `plot.survatr_result()`

For `type %in% c("survival", "risk")`:

- x-axis: `time`
- y-axis: `s_hat` (or `risk_hat`)
- one line per intervention, colored
- CI ribbon when `ci_lower` / `ci_upper` present

For `type %in% c("risk_difference", "risk_ratio", "rmst_difference")`:

- x-axis: `time`
- y-axis: `estimate`
- one line per contrast; ribbon for CI; reference line at 0 (RD / RMST
  difference) or 1 (RR)

For `type == "rmst"`:

- x-axis: `time`
- y-axis: `rmst_hat`
- one line per intervention; monotone non-decreasing ribbon

Backend: **base R graphics** via `graphics::plot` / `lines` / `polygon`.
Rationale: no ggplot2 dependency; matches tidyverse-skeptic leaning of
the etverse codebase (and causatr's `plot.causatr_result` is base R
too).

Arguments:

```r
plot(x, which = c("auto", "curves", "contrasts"),
     main = NULL, xlab = NULL, ylab = NULL, col = NULL, ...)
```

`which = "auto"` picks `curves` for `survival` / `risk` / `rmst` and
`contrasts` for the three contrast-flavored types.

### `forrest.survatr_result()`

Forest plot at a single reference time `t_ref` (user-supplied). Rows
are contrasts (one per non-reference intervention); columns show
estimate, 95% CI, and a textual label. Calls causatr's `forrest`
mechanism if it is exported; otherwise roll our own.

Actually causatr does **not** export a `forrest` generic (checked the
namespace). Define our own here. Signature:

```r
forrest(x, t_ref, col = NULL, main = NULL, ...)
```

Base-R implementation using `plot.default` + `arrows`. Rejects if
`t_ref` is not in `x$time_grid`.

### `print.survatr_result()` polish

Extend to show the type-appropriate columns (`s_hat` / `risk_hat` /
`rmst_hat` / `estimate`) plus CIs when available. Currently shows the
head of `contrasts` or `estimates`; that stays, but widen the column
set so CI columns are not clipped by default.

## Deliverables

### New R files

| File | Contents |
|---|---|
| `R/variance_bootstrap.R` | **Internal** `bootstrap_survival()` -- resamples individuals, refits the hazard model, returns per-replicate `s_hat^a` and contrast estimates. Uses `boot::boot()` for parallel dispatch and `boot::boot.ci()` for percentile / normal CIs. |
| `R/tidy.R` | **Exported** `tidy.survatr_result()`. S3 method on the `tidy` generic from `generics` (causatr also imports `generics`; we do the same). |
| `R/plot.R` | **Exported** `plot.survatr_result()`. Base-R graphics. |
| `R/forrest.R` | **Exported** `forrest()` generic + `forrest.survatr_result()` method. |

### Updated R files

| File | Change |
|---|---|
| `R/contrast.R` | Accept `ci_method = "bootstrap"`. New args `n_boot`, `boot_ci`, `parallel`, `ncpus`, `seed`. New validators `validate_n_boot()`, `validate_boot_ci()`. The `survatr_ci_not_available` branch on `"bootstrap"` is removed. |
| `R/print.R` | Minor column-set polish so CI columns are shown when populated. |
| `R/survatr-package.R` | Add `@importFrom boot boot boot.ci`, `@importFrom generics tidy`. `generics` gets added to DESCRIPTION Imports (causatr already does this; mirror the pattern). |
| `DESCRIPTION` | Add `generics` to `Imports:` (runtime; pulled in by the `tidy` re-export). |

### Tests (`tests/testthat/`)

| File | Asserts |
|---|---|
| `test-bootstrap-survival.R` | End-to-end `contrast(..., ci_method = "bootstrap")` on a constant-hazard DGP: returns populated CI columns; percentile CIs cover truth at nominal 95% (200 reps of B = 200, or a single 500-replicate run). |
| `test-bootstrap-vs-sandwich.R` | Per-time SE from bootstrap (B = 500) agrees with sandwich SE within 25% on a moderate-size DGP. Skipped on CRAN. |
| `test-bootstrap-risk-ratio.R` | Percentile CI on RR is strictly positive and covers 1 on a no-effect DGP; does NOT apply the log-scale / exponentiate that sandwich uses (percentile is invariant to monotonic transforms). |
| `test-bootstrap-rejections.R` | `n_boot <= 0` / non-integer / huge, `boot_ci` not in the valid set, `parallel` mis-specified -- classed errors. |
| `test-tidy-survatr_result.R` | Happy-path: long-format output, expected column set, correct rowcount for each `which` option. Error on bad `which`. |
| `test-plot-survatr_result.R` | Calls `plot()` under each `type` and captures the output with `vdiffr::expect_doppelganger()` (if `vdiffr` available) or falls back to a `expect_silent` run. Mark 🟡 if we go the `expect_silent` route. |
| `test-forrest-survatr_result.R` | `forrest(res, t_ref)` produces a base-R plot; rejects `t_ref` not in `time_grid` with `survatr_bad_t_ref`. |

## API contract

Everything from chunks 2 / 3 plus:

```r
res_boot <- contrast(
  fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times         = seq(0, 120, by = 12),
  type          = "risk_difference",
  reference     = "a0",
  ci_method     = "bootstrap",
  n_boot        = 1000L,
  boot_ci       = "percentile",
  parallel      = "multicore",
  ncpus         = 4L,
  seed          = 2026L
)

tidy(res_boot)                              # long data.frame
plot(res_boot)                              # RD curve with percentile ribbon
forrest(res_boot, t_ref = 120)              # forest at t* = 120
```

## Behaviour rules

- **Bootstrap resamples ids, not PP rows.** Within-id dependence must be
  preserved. Implement by sampling from `unique_ids` with replacement,
  then concatenating `fit$pp_data[get(id_col) == sampled_id]` blocks.
  Re-id with a replicate-local counter so a doubled id doesn't create
  ambiguity.
- **Refit the hazard model per replicate.** Re-use `fit$model_fn`,
  `fit$confounders`, `fit$time_formula`, `fit$family` (via weights /
  unweighted). Do NOT reuse the original `fit$model`'s coefficients.
- **Percentile CI for all types** including RR, because percentile is
  transform-invariant. SE-based (`boot_ci = "wald"`) needs a log-scale
  wrapping for RR -- handle in the same branch as sandwich.
- **Replicate failures degrade gracefully.** If a single replicate fails
  (e.g. all-zero events in a resampled dataset), record NA and continue.
  If > 10% of replicates fail, abort with `survatr_boot_failed`.
- **`seed` is respected.** When non-null, `set.seed(seed)` before calling
  `boot::boot()` so the whole replicate sequence is reproducible.
- **Do NOT memoize the expensive bits across replicates** -- each
  replicate sees a different sample, so `prepare_sandwich_shared()`
  cannot be reused.

## Non-goals (deferred)

- BCa CIs -- `boot_ci = "bca"` can be a chunk-4.5 follow-up if users
  ask.
- Cluster bootstrap by a second variable -- out of scope; ids are the
  cluster.
- Simultaneous bands -- still out of v1 scope.

## Seed sources

- `causatr:::variance_bootstrap()` for the resampling boilerplate -- we
  need the individual-level resample + refit pattern, but causatr's
  version operates on scalar fits so the contrast-replicate loop is
  new code.
- `boot::boot()` + `boot::boot.ci()` for the parallel dispatch and CI
  helpers.

## Acceptance checklist

- [ ] `devtools::load_all()` succeeds.
- [ ] `devtools::document()` regenerates Rd cleanly.
- [ ] `devtools::test()` -- all chunk-1/2/3 tests unchanged; new chunk-4
      tests green. `test-bootstrap-vs-sandwich.R` skips on CRAN.
- [ ] `devtools::check()` -- 0 errors, 0 warnings. Unused-imports NOTE
      should disappear entirely now that `boot` is used.
- [ ] `air format .` is a no-op.
- [ ] [FEATURE_COVERAGE_MATRIX.md](FEATURE_COVERAGE_MATRIX.md) updated
      with a "Bootstrap variance" sub-section and S3 rows (`tidy` /
      `plot` / `forrest`).
- [ ] Single commit with message `feat(chunk-4): Track A bootstrap +
      S3 polish -- tidy / plot / forrest`.
