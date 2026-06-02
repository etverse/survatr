# survatr — Scope, Integration & Critical-Review (round 2)

**Date:** 2026-06-02. Generated during a resume-the-package session after a long
gap (active work had moved to `causatr`). Three workstreams:

1. **causatr integration health** — does survatr's reliance on the causatr engine still hold?
2. **Systematic survival-analysis scope review** — the field-driven scoping review that was *never done* (scope was transplanted wholesale from causatr's handoff at the split).
3. **Critical-review round 2** — adversarial pass over the whole `R/` tree (round 1 was 2026-04-22).

This file is the single resume artifact for that session. Per-finding fix status
is tracked in `NEWS.md` + commits; scope decisions land in `CLAUDE.md`,
`SURVIVAL_PACKAGE_HANDOFF.md` §10, and `FEATURE_COVERAGE_MATRIX.md`.

State at review time: chunks 1–4 done (Track A: fit, contrast, sandwich, bootstrap, S3). Chunks 5–10 not started.

---

## 1. causatr integration health — VERDICT: healthy

survatr has exactly **two** runtime dependencies on causatr internals:

| Symbol | Access | survatr call site | causatr def | Status |
|---|---|---|---|---|
| `apply_intervention(data, treatment, iv)` | `:::` | `R/apply_intervention.R:24` | `causatr/R/interventions.R:439` | exists, `@noRd`, signature + copy-semantics intact ✅ |
| `prepare_model_if(model, fit_idx, n_total)` | `:::` | `R/variance_if_survival.R:200` | `causatr/R/variance_if_core.R:309` | exists, `@noRd`, returns `list(X_fit, B_inv, r_score, …)` intact ✅ |

Intervention constructors (`static`/`shift`/…), `to_person_period`, `contrast`
appear in survatr's `R/` **only in roxygen + error strings** — users pass
intervention objects in; survatr never calls them directly. `hard-rules.md`
lists `apply_model_correction`/`vcov_from_if`/`variance_if_numeric` as intended
imports but **none are currently called** — survatr rolls its own cross-time IF
aggregation. They are latent future deps only.

Copied helpers in sync: `is_uncensored`, `check_weights`, `check_dots_na_action`
(all with the sanctioned `survatr_*` class rename). `check_reserved_cols` is
**ahead** of causatr (survatr added `call=` threading — causatr should backport).

**Latent risks:**
- `DESCRIPTION` floats `Remotes: github::etverse/causatr` **unpinned** while
  depending on two `causatr:::` internals → a causatr refactor could silently
  break survatr's sandwich/contrast paths with no version gate. **Pin + assert
  the internal-API shape at load.** (Top integration recommendation.)
- survatr's `contrast()` S3 generic **shadows** causatr's exported `contrast()`
  function when both are attached. Intentional + documented; `causatr::contrast()`
  is the workaround. Belongs in a user-facing vignette.

**Improvements worth porting from causatr (post-split):**
- `iv_design_matrix()` GLM/GAM split (`variance_if_core.R:222-235`) + GAM `$Vp`
  bread + missing-`$Vp` hard-abort → **this is the fix for the GAM×sandwich bug below.**
- Singular-bread `MASS::ginv()` fallback with rate-limited classed warning
  (`variance_if_core.R:167-191`). survatr inherits it inside `prepare_model_if`
  but its own delta-chain algebra has no rank-deficiency guard.
- `helper-*-oracle.R` + delicatessen Python-fixture testing pattern → needed for
  chunks 6–9 (NHEFS, Track B gfoRmula oracle, CR sandwich).
- (chunk 5+) per-period weight truncation `trim` + `replay_fit()` dots discipline.

---

## 2. Systematic survival-analysis scope review

### 2.1 What scope currently covers

Done: Track A gcomp fit, survival/risk/RD/RR/RMST/RMST-diff curves, sandwich
(delta-method cross-time IF), bootstrap, S3 (`print`/`plot`/`tidy`/`forrest`),
survey weights, `censoring=` as row filter.

Planned: IPW (5), Track B ICE (6), competing risks cause-specific+CIF (7),
matching-rejection surface (8), NHEFS Ch.17 replication (9), `diagnose()` (10).

Explicitly out (documented): Cox modelling, TMLE/ML survival, forward-sim
g-formula, Fine–Gray, matching+survival, simultaneous bands in v1.

### 2.2 The NOT-CONSIDERED list (the point of the review)

Items with **no trace in any scope doc** — silent absences, not decisions:

| Bucket | Items |
|---|---|
| **Outputs the existing curve+IF spine can already produce** | survival quantiles / median survival; RMTL; years-of-life-lost (CR); cluster-robust SE; simultaneous/uniform confidence bands |
| **Scoping limbo (mentioned, unowned)** | AIPW/doubly-robust survival (handoff §7 says "composes cleanly", no chunk); IPCW (called "the motivating use case", no chunk); stochastic+survival ("pending", no chunk) |
| **Unmentioned design/estimand topics** | left-truncation / delayed entry; recurrent events; multi-state / illness-death; semi-competing risks; landmark analysis; target-trial-emulation framing; immortal-time-bias check; interval censoring |
| **Unmentioned data/usability** | clustered/multilevel data; missing-data / MI story (NA in predictors is *rejected* with no documented workaround); large-data performance contract |
| **Unmentioned caveats** | truncation-by-death / principal-stratification (the biggest interpretation trap in CR causal inference); censoring-by-arm summary; ESS in diagnose |

### 2.3 Gaps ranked (demand × low cost on existing spine)

1. **Survival quantiles / median survival** — HIGH value, LOW cost. Smooth functional of the existing S^a(t) curve; CI via the IF columns or bootstrap. `adjustedCurves` headline feature. → new `type="quantile"`.
2. **RMTL + years-of-life-lost** — HIGH (post chunk 7), LOW. Complement/integral of computed quantities; reuses trapezoidal quadratic-form SE. → fold into chunk 2 estimands + chunk 7.
3. **Cluster-robust SE** — MEDIUM, LOW. Sum per-individual IF within cluster before `crossprod`. → `cluster=` arg in chunk-3 engine.
4. **Simultaneous bands** — MEDIUM, MEDIUM. Multiplier-bootstrap reuses the IF matrix directly. → revisit "out in v1"; first post-v1 inference chunk.
5. **AIPW / doubly-robust survival** — HIGH methodological value, MEDIUM-HIGH cost. Every sibling ships it; singly-robust gcomp+IPW is the weakest point vs the field. → resolve in-or-out (see 2.4b).
6. **IPCW → real chunk** — MEDIUM. Parallels chunk-5 IPW machinery. → assign a chunk.
7. **Left-truncation / delayed entry** — MEDIUM, LOW. Pooled-logistic handles it natively (id enters risk set at a delayed period). → small risk-set extension or explicit exclusion.
8. **MI / missing-data story** — MEDIUM. mice + per-imputation fit + Rubin pooling over the curve. → vignette + optional helper.
9. **Diagnostics gaps** (censoring-by-arm, ESS, truncation-by-death caveat) — enumerate inside chunk 10.

### 2.4 Scope decisions to RATIFY or REVISIT (maintainer's call)

- **(a) Simultaneous bands "out of v1"** — for a *curve-valued* estimand this matters more than for scalars; cheap given the IF matrix exists. Recommend: keep out of v1 but explicitly schedule, don't leave as a bare "out".
- **(b) AIPW limbo** — RESOLVE. Either promote parametric AIPW to a chunk (distinct from the TMLE/ML exclusion, which stays out), or ratify "singly-robust only, defer DR to lmtp/riskRegression" in CLAUDE.md. Currently neither — planned-but-unowned.
- **(c) IPCW unchunked** — internal inconsistency (called "motivating" but no chunk). Assign a chunk after IPW.
- **(d) Left-truncation / recurrent / multi-state / landmark / target-trial / immortal-time** — silent absences. Make an explicit ratify-or-exclude call on each in CLAUDE.md "Scope".
- **(e) Survival quantiles + RMTL** — should be IN; their absence is the clearest fingerprint of scalar-package transplant.
- **(f) Truncation-by-death caveat** — chunk 7 should ship a documented caveat, not silent numbers.
- **(g) Cluster-robust SE** — ratify as in-scope; near-free given the IF matrix.

Comparator anchors: `riskRegression::ate` (AIPTW + CR + IF SE — closest sibling),
`adjustedCurves` (g-comp/IPW/AIPW adjusted curves + quantiles + bands),
`concrete` (CT-TMLE CR), `stdReg2` (DR standardization), `lmtp` (longitudinal MTP/TMLE — already a survatr oracle).

---

## 3. Critical-review round 2 — triage

Round 1 (2026-04-22) hardened the obvious correctness surface, so round 2 is thin
on bugs (honest outcome). Findings, de-noised:

### Correctness
| ID | Finding | Verdict |
|---|---|---|
| A1 | **GAM × sandwich basis mismatch** (`variance_if_survival.R:57-58`): `stats::model.matrix(delete.response(terms(model)), …)` does not reproduce `mgcv::gam`'s penalized `lpmatrix` basis, but `B_inv` from `prepare_model_if` uses `$Vp` on that basis → mismatched bases. `hard-rules.md:34` advertises GAM `s(t)` as supported. | **REAL** — fix (guard or port `iv_design_matrix`). Verified by reading the code path. |
| A2 | **Duplicate `survatr_bad_t_ref` class** (`forrest.R:73,80`) for two distinct failures; 2nd branch near-unreachable (`t_ref ∈ time_grid` already enforced). | **REAL, minor** — distinct class. |
| A3 | hazards length-check missing before `pp[, .cf_hazard := hazards]` (`survival_curve.R`). | Defensive — add cheap guard. |
| A4 | risk_difference/RMST Wald CI can go below 0. | **Not a bug** (Wald behavior) — document, no code change. |
| A5 | NSE risk: `tbl[get(group_col) == g]` (`plot.R:128`) — hypothetical `g` column. | Cheap hardening — `.env$g`. |
| A6 | No dim assertion on `IF_mat` before `crossprod` (`variance_sandwich.R`). | Defensive — add `stopifnot`. |
| A7 | "singular-bread guard in survatr's own algebra". | survatr does no unguarded `solve()`; `B_inv` guarded upstream. **Likely retract** — verify. |

### Compliance (process gaps vs implement-feature skill)
- C1: missing `@examples` on exports (`surv_fit`, `contrast.survatr_fit`, S3 methods).
- C2: math `@details` + source on `rmst.R` (`trapezoidal_rmst`, `rmst_weights`), `variance_if_survival.R`, `variance_sandwich.R`.
- C3: split `contrast.R` (512 lines → core + validators).
- C4: split `variance_bootstrap.R` (360 lines → bootstrap + ci).
- C5: `@family` / `@description` tags on exports + S3 methods.
- C6: inline comments (mu_eta cancellation = h(1−h)/(1−h)=h for logit; weight-subset rationale; cumprod invariant).
- C7: internal validators in `contrast.R` missing `@param` for `call`.

### 3.1 Round-2 resolution (all findings closed)

| # | Type | Resolution |
|---|---|---|
| A1 | Correctness | **Fixed** `a0eee89` — lpmatrix design via `causatr:::iv_design_matrix()` + `predict.gam` 1-D-array coercion. Validated GAM sandwich ≈ GLM sandwich within 2%; new `test-sandwich-gam.R`. |
| A2 | Correctness | **Fixed** `db613de` — distinct `survatr_forrest_no_contrasts` class. |
| A3 | Correctness (defensive) | **Fixed** `d8e6569` — `survatr_hazard_misaligned` length guard. |
| A4 | Correctness | **Retracted** — symmetric Wald CI on a *difference* is valid (not a probability); survival/risk Wald-bound outside [0,1] is a known Wald property → routed to scope (cloglog option), not a bug. |
| A5 | Correctness (defensive) | **Fixed** `d8e6569` — plain-R group mask in `plot()` (no `g`-column NSE shadow). |
| A6 | Correctness (defensive) | **Fixed** `d8e6569` — IF-matrix dimension assertion before `crossprod`. |
| A7 | Correctness | **Retracted** — survatr has no `solve()`/`ginv` of its own; the only bread inversion is causatr's already-guarded `bread_inv()`. |
| B2 | Integration | **Guarded** `c72b0b6` — `test-causatr-integration.R` pins the 3 `causatr:::` internals' contract. Hard-pin of the remote deferred to release time. |
| C1–C7 | Compliance | **Fixed** `a0150d1` — file splits (`contrast.R`→`contrast_validators.R`, `variance_bootstrap.R`→`variance_bootstrap_ci.R`), `@examples` on exports, math `@details` (H&R 2020 Ch.17 + chunk docs), `@family`/`@description`, inline comments. |

Full suite after all fixes: all pass, 1 pre-existing skip (lmtp oracle declines the toy DGP). New durable invariants recorded in `.claude/hard-rules.md` (§ round-2).

### 3.2 Scope decisions — RATIFIED 2026-06-02

§2.4 was resolved with the maintainer. Every item resolved toward inclusion;
the roadmap grew by 8 chunks (11–18), folded into
`SURVIVAL_PACKAGE_HANDOFF.md` §10 + the CLAUDE.md chunk table, each with its own
`CHUNK_<n>_*.md` plan.

| Decision | Outcome | Chunk |
|---|---|---|
| Doubly-robust (AIPW) survival | **In** — parametric only (ML/TMLE stays out → lmtp/concrete) | 15 |
| IPCW (own chunk) | **In** (was "motivating" but unchunked) | 11 |
| Survival quantiles / median | **In** | 12 |
| RMTL + years-of-life-lost | **In** | 12 |
| Cluster-robust SE | **In** | 13 |
| Left-truncation / delayed entry | **In** | 14 |
| Simultaneous confidence bands | **In** (v2) | 16 |
| Target-trial / landmark / immortal-time | **In** (v2) | 17 |
| Recurrent / multi-state events | **In** (v2; large lift) | 18 |

Clear-cut items applied without a question: truncation-by-death caveat ships
with chunk 7; survival/risk CIs stay raw Wald for v1 (cloglog transform
deferred). **Phasing:** v1 = 1–10, v1.x = 11–15, v2 = 16–18. **Out of scope
(ratified):** ML/TMLE survival (→ lmtp), forward-sim g-formula (→ gfoRmula),
Fine–Gray subdistribution hazards, Cox-PH modelling, matching.
