# Chunk 7 — Competing risks: cause-specific hazards + CIF

> **Status: ✅ Done.** `R/competing_risks.R` + `R/variance_if_competing.R`;
> tests `test-competing-risks.R`, `test-competing-risks-sandwich.R`; oracle
> `helper-cr-oracle.R` + `data-raw/delicatessen_competing_risks.py`. gcomp /
> point-treatment g-computation only; IPW / ICE competing risks and per-cause RMST / YLL deferred.
> **Depends on:** Chunk 2 (contrast spine), Chunk 3 (sandwich IF).
> **Oracle:** closed-form two-cause constant-hazard DGP (analytic CIF);
> `riskRegression::ate()` / `survtmle` for CIF point estimates;
> `delicatessen` (Python) for the stacked-EE sandwich (no single-call R oracle).

## Goal

First-class competing risks via **cause-specific hazards + cumulative incidence
functions (CIF)**. For `J` competing event types, fit `J` parallel
cause-specific hazard models (each treating other causes as censored at their
event time), then build the CIF under intervention `a`. Contrast CIFs (and the
all-cause survival). Sandwich variance via stacked estimating equations across
the cause-specific models. **Fine–Gray / subdistribution hazards are out of
scope** — document explicitly.

## The math (Hernán & Robins Ch. 17)

Cause-specific hazard for cause `j`, fit on person-period rows treating events
of type `j' ≠ j` as censored at their time:

```
logit h^(j)(t | A, L) = alpha_j(t) + beta_{A,j} A + beta_{L,j} L
```

All-cause survival uses the **sum** of cause-specific hazards:

```
S^a(k) = prod_{m <= k} (1 - Σ_j ĥ^(j),a(m))
```

Cause-`j` cumulative incidence under intervention `a`:

```
F^(j),a(t) = Σ_{k=1}^{t} S^a(k-1) · ĥ^(j),a(k)
```

(`S^a(0) = 1`.) Per-individual versions cumulate within id and average across
ids (the Jensen-safe rule).

**Variance.** Stack the `J` cause-specific hazard-model scores. The CIF IF
combines (a) the all-cause survival IF (which depends on every cause's hazard
through the sum) and (b) the cause-`j` hazard IF, propagated through the CIF
sum by the delta method. The per-individual IF is again an `n × |t-grid|`
matrix; full time-covariance is `crossprod(IF) / n^2`. CIF-difference and
CIF-ratio contrast IFs follow the chunk-3 pattern.

## Deliverables

### New R files
- `R/competing_risks.R` — parallel cause-specific fit (one `surv_fit`-style
  model per cause, others censored), all-cause survival from the hazard sum,
  CIF computation, CIF contrasts, stacked-EE sandwich, CIF bootstrap.

### Updated R files
- `R/surv_fit.R` — `competing = <cause-col>` now **invokes** the cause-specific
  path (currently aborts with `survatr_competing_misuse`; the abort stays for
  the *misuse* case — passing a competing column without the CR entry point —
  but the proper path is enabled here).
- `R/contrast.R` — CIF / CIF-difference / CIF-ratio estimands; result carries a
  `cause` dimension.
- `R/variance_if_survival.R` — stacked cause-specific IF assembly.

### Tests (`tests/testthat/`)
- `test-competing-risks.R` — analytic two-cause constant-hazard DGP: CIF matches
  the closed form `F^(j)(t) = Σ_k S(k-1) h_j` to tolerance; `Σ_j F^(j)(t) + S(t)
  = 1` identity holds.
- `test-competing-risks-sandwich.R` — coverage simulation; CIF-difference CI
  covers 0 under a no-effect DGP. Cross-validate the stacked-EE SE against a
  `delicatessen` Python fixture (store under `tests/testthat/fixtures/python/`,
  commit the CSV).
- `helper-cr-oracle.R` — analytic CIF + `riskRegression::ate` comparator.

## API contract

```r
fit <- surv_fit(data, outcome = "event_type",   # 0 = censored, 1..J = cause
                treatment = "A", confounders = ~ L,
                id = "id", time = "t",
                competing = "event_type")        # invokes the CR path
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "cif_difference", cause = 1,
                   ci_method = "sandwich")
# result$estimates: intervention | cause | time | cif_hat | se | ci_*
# result$contrasts: contrast | cause | time | estimate | se | ci_*
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Cause-specific hazards only.** No Fine–Gray / subdistribution hazards
  (different estimand, different data structure). Document the choice in the
  vignette + roxygen.
- All-cause survival uses **Σ_j ĥ^(j)**; never a single cause's hazard.
- The `Σ_j F^(j)(t) + S(t) = 1` identity must hold numerically.
- `competing = <col>` without the CR entry point remains a
  `survatr_competing_misuse` error (no silent cause-deleted hazard).
- **Truncation-by-death caveat.** Cause-specific CIF contrasts carry a
  well-known interpretational hazard (conditioning on survival to a competing
  event). Ship an explicit caveat in `diagnose()` / the vignette — do not emit
  numbers silently.

## Non-goals (deferred)
- Fine–Gray subdistribution hazards (out of scope, documented).
- Years-of-life-lost per cause (integral of CIF) — **pending scope
  ratification**; would fold in here cheaply once approved.
- Competing risks under longitudinal ICE-hazard — composes after Chunk 6.

## Dependencies & composition
- Chunks 2, 3. Reuses `causatr:::prepare_model_if()` per cause-specific model.

## Acceptance checklist
- [x] `J` cause-specific hazards fit; CIF matches the analytic two-cause form
      (and the Aalen–Johansen estimator).
- [x] `Σ_j F^(j)(t) + S(t) = 1` to tolerance (1e-12).
- [x] Stacked-EE sandwich CI for CIF-difference covers 0 under no effect
      (≥ 88% nominal 95%); cross-checked against `delicatessen` (~1e-4) and the
      bootstrap (~2%).
- [x] `competing =` misuse still aborts (`survatr_competing_misuse` /
      `survatr_competing_estimator` / `survatr_bad_competing` /
      `survatr_competing_type`); proper path works.
- [x] Truncation-by-death caveat shipped (`print` note + one-time `inform` +
      vignette).
- [x] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md +
      `.claude/hard-rules.md` + `NEWS.md` updated.
