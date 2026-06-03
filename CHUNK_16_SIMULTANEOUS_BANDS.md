# Chunk 16 — Simultaneous (uniform) confidence bands

> **Status: ⬜ Not started · post-v1**
> **Depends on:** Chunk 3 (sandwich IF: the `n × |t-grid|` per-individual IF
> matrix).
> **Oracle:** `adjustedCurves` simultaneous bands; a coverage simulation
> (simultaneous coverage ≥ nominal across the whole grid; the pointwise band
> under-covers when read as a simultaneous statement).

## Goal

Add **simultaneous (uniform) confidence bands** over the survival curve (and any
curve-valued contrast), so a single band has the stated coverage for the
**entire curve at once** — the inferentially honest object for a curve-valued
estimand, where the shipped pointwise bands under-cover when interpreted across
all time points. The method is a **multiplier (wild) bootstrap** that reuses the
existing `n × |t-grid|` IF matrix directly: Gaussian multipliers × IF rows,
sup-statistic over the time grid, empirical critical value. No refit, no new
estimation. The handoff scopes simultaneous bands **out of v1**; this chunk is
explicitly the **first post-v1 inference chunk** that schedules them rather than
leaving a bare "out".

## The math

Chunk 3 stores the per-individual influence function `IF` as an `n × |t-grid|`
matrix; the pointwise SE at time `t_j` is `σ̂(t_j) = sqrt( (1/n²) Σ_i IF_i(t_j)² )`.
A **pointwise** `(1−α)` band is `Ŝ^a(t_j) ± z_{1−α/2} · σ̂(t_j)` at each `t_j`
independently; read across all `j` it under-covers because it ignores the
multiplicity of the grid.

The multiplier bootstrap draws `B` replicates of i.i.d. standard-normal
multipliers `ξ^(b) = (ξ_1^(b), …, ξ_n^(b))` and forms a perturbed,
**standardized** curve process:

```
G^(b)(t_j) = (1 / (n · σ̂(t_j))) Σ_i ξ_i^(b) · IF_i(t_j)
sup-statistic:  M^(b) = max_j | G^(b)(t_j) |
```

The simultaneous critical value `c_{1−α}` is the empirical `(1−α)` quantile of
`{ M^(1), …, M^(B) }`. The **simultaneous** band is

```
Ŝ^a(t_j) ± c_{1−α} · σ̂(t_j)      for all t_j on the grid
```

By construction `c_{1−α} ≥ z_{1−α/2}`, so the simultaneous band is **wider**
than the pointwise band and contains it. Standardizing by `σ̂(t_j)` inside the
sup makes the band an equal-precision (proportional-to-pointwise-SE) band — the
standard construction reused by `adjustedCurves`. The same machinery applies to
any curve-valued contrast (risk-difference curve, CIF curve): substitute that
functional's IF columns and pointwise SE.

This is a pure post-processing step on the IF matrix — it composes with every
estimator (gcomp, ipw, ice, aipw) and every curve estimand, because it operates
on the variance currency, not the functional.

## Design

- A new `band = c("pointwise", "simultaneous")` argument on `contrast()`
  (default `"pointwise"` — the shipped behaviour, unchanged). `"simultaneous"`
  triggers the multiplier bootstrap on the stored IF matrix.
- A `n_multiplier =` argument (default e.g. 2000) and a `seed =` for
  reproducibility (L'Ecuyer-CMRG under parallel, per the bootstrap RNG
  invariant in `bootstrap_survival()`).
- The result carries both the pointwise SE (always) and, when requested, the
  simultaneous critical value `c_{1−α}` and the resulting band columns; the
  band type is recorded on the `survatr_result` so `print` / `plot` render it.

## Deliverables

### New R files
- `R/bands_simultaneous.R` — the multiplier-bootstrap critical-value computation
  on the `n × |t-grid|` IF matrix (Gaussian multipliers, standardized sup
  process, empirical quantile), returning `c_{1−α}` and the band columns for the
  level curve and for curve-valued contrasts.

### Updated R files
- `R/contrast.R` — accept `band =`, `n_multiplier =`, `seed =`; route to the
  simultaneous path; attach the band type + critical value to the result.
- `R/variance_sandwich.R` — expose the standardized IF columns (IF / pointwise
  SE) the multiplier step consumes (the pointwise SE is already computed here).
- `R/plot.R` — render simultaneous bands distinctly from pointwise (e.g. a
  wider, lighter ribbon outside the pointwise one) when the result carries a
  simultaneous band.
- `R/print.R` / `R/tidy.R` — surface the band type and critical value.

### Tests (`tests/testthat/`)
- `test-bands-simultaneous.R` — coverage simulation on a known-curve DGP: the
  simultaneous band covers the **whole** true curve with probability ≥ `1−α`,
  while the pointwise band's whole-curve coverage is materially below nominal;
  `c_{1−α} ≥ z_{1−α/2}` (the band contains the pointwise band); reproducibility
  under a fixed seed; structure vs `adjustedCurves` simultaneous bands.

## API contract

```r
result <- contrast(
  fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times = seq(0, 120, 12),
  type = "survival",
  ci_method = "sandwich",
  band = "simultaneous",      # NEW — default "pointwise"
  n_multiplier = 2000,        # NEW
  seed = 1                    # NEW — reproducible multipliers
)
# result$estimates gains simultaneous band columns + the critical value c_{1-α};
# pointwise SE is always present.
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Reuse the existing `n × |t-grid|` IF matrix.** The multiplier bootstrap is
  pure post-processing — no refit, no new estimation, no change to the point
  estimate.
- **Pointwise remains the default** (`band = "pointwise"`); the shipped chunk-3
  behaviour is unchanged unless the user opts in.
- **The simultaneous band must contain the pointwise band** (`c_{1−α} ≥
  z_{1−α/2}`); a simultaneous band narrower than pointwise is a bug.
- **Standardize the sup process by the pointwise SE** (equal-precision band) —
  this is the construction the coverage test and `adjustedCurves` assume.
- **Multiplier reproducibility uses L'Ecuyer-CMRG under parallel backends** (the
  same RNG invariant `bootstrap_survival()` enforces); restore the caller's
  RNGkind on exit.
- This is **post-v1** — it does not block the v1 release; it explicitly retires
  the handoff's bare "simultaneous bands out of v1" note by scheduling them.

## Non-goals (deferred)
- **Hall–Wellner / Equal-Precision analytic bands** (closed-form transformed-
  scale bands) — the multiplier bootstrap is the chosen route; analytic bands
  are a possible later alternative, not this chunk.
- **Confidence bands on the quantile process** — composes with Chunk 12's
  quantiles once both land.
- **Simultaneous bands across multiple contrasts at once** (family-wise over
  interventions and time) — single curve / single contrast here.

## Dependencies & composition
- Chunk 3 only (the IF matrix). Composes with every estimator and curve
  estimand (chunks 5/7/11/12/15) for free, because it operates on the IF, not
  the functional. No new causatr internals.

## Acceptance checklist
- [ ] `band = "simultaneous"` produces a uniform band via multiplier bootstrap
      on the IF matrix; `c_{1−α} ≥ z_{1−α/2}`.
- [ ] Coverage simulation: simultaneous whole-curve coverage ≥ nominal;
      pointwise whole-curve coverage materially below nominal.
- [ ] Reproducible under a fixed seed (L'Ecuyer-CMRG under parallel).
- [ ] `plot()` renders the simultaneous band distinctly; `print` / `tidy`
      surface the band type + critical value.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
