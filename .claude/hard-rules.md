# survatr hard rules

Project-specific rules that override / extend the etverse-wide rules at
`~/Documents/personal/software/etverse/.claude/skills/*/SKILL.md`. Read by the
`implement-feature` and `critical-review-loop` skills before they do anything.

## Project conventions

- **Design-doc pattern.** `SURVIVAL_PACKAGE_HANDOFF.md` is the single source of
  truth for scope and design. Per-chunk implementation guides live at the repo
  root as `CHUNK_*.md` (e.g. `CHUNK_1_SKELETON.md`,
  `CHUNK_2_CONTRAST_A.md`). When a skill says "read the design doc", read the
  handoff first and then the relevant chunk doc.
- **Feature coverage file.** `FEATURE_COVERAGE_MATRIX.md` (same pattern as
  causatr). Every PR that adds, removes, or changes a feature MUST update
  this file.
- **Error-class prefix.** `survatr_*` for all `rlang::abort()` calls.
- **Repro-script prefix.** `/tmp/survatr_repro_<slug>.R`.
- **Test-log paths.** `/tmp/survatr-test-<file>.txt` for per-file runs,
  `/tmp/survatr-test-results.txt` for the full suite.
- **Reserved column prefix.** `.survatr_*` for internal person-period
  bookkeeping columns (e.g. `.survatr_prev_event`, `.survatr_prev_cens`).
  Guard against user-data collisions upfront.

## Supported dimensions (for combination audits)

| Dimension | Values |
|---|---|
| **Track** | A (point treatment + pooled-logistic hazard), B (longitudinal ICE-hazards) |
| **Estimator** | gcomp (pooled-logistic), ipw (weighted hazard MSM), ice (hazard pseudo-outcome). **Matching: hard-abort.** |
| **Treatment timing** | point (Track A), longitudinal (Track B) |
| **Treatment type** | binary, continuous, categorical (k>2), count (Poisson/NB, IPW only), multivariate (via causatr inheritance) |
| **Outcome family** | binomial hazard (first-step / indicator), quasibinomial (pseudo-outcome / weighted fits) |
| **Model class** | GLM (pooled logistic), GAM (via `mgcv::gam` with `s(t)` for baseline hazard) |
| **Intervention** | static, shift, scale_by, threshold (gcomp only), dynamic, ipsi (IPW only), stochastic (pending) |
| **Estimand** | survival S^a(t), risk 1 - S^a(t), risk difference, risk ratio, RMST up to t* |
| **Contrast type** | difference, ratio |
| **Variance method** | sandwich (delta-method cross-time IF), bootstrap (resample individuals), numeric Tier 1/2 fallback |
| **Weights** | none, survey/external, censoring row-filter, IPCW (cumulative, per-period) |
| **Competing risks** | cause-specific hazards + CIF (first-class); Fine-Gray / subdistribution hazards out of scope |

## Hard rules (appended to the skill's generic rules)

### Architecture invariants — DO NOT flag these as bugs without a numerical reproducer

- **Pooled-logistic hazard on person-period rows** is the outcome model for
  all of Track A. `logit h(t | A, L) = alpha(t) + beta_A A + beta_L L`. The
  risk set drops rows at/after first event per id and rows at/after first
  censor per id. Don't substitute a Cox partial likelihood — discrete-time
  hazard is the contract.
- **Survival curve is a cumulative product over predicted hazards per
  individual.** `S^a_i(t) = prod_{k <= t} (1 - h^a_{i,k})`, then average
  across individuals. Averaging hazards before cumulating is biased by
  Jensen's inequality (the cumulative product is nonlinear) — never do it.
- **Variance propagates through the cumulative product via delta method.**
  `d S^a_i(t) = -S^a_i(t) * sum_{k <= t} d h^a_{i,k} / (1 - h^a_{i,k})`. The
  outcome-model IF (via `causatr:::prepare_model_if()` /
  `apply_model_correction()`) produces per-row IFs on hazards; aggregate
  across rows within individual via the delta chain to get per-individual
  IFs on `S^a_i(t)`. Don't try to use the scalar-outcome IF engine directly.
- **Per-individual IF is an `n x |t-grid|` matrix.** Full time-covariance is
  `crossprod(IF_mat) / n^2`. Pointwise SE is the diagonal. RMST SE uses
  `a^T V a` with `a` the trapezoidal weights. **Simultaneous bands are out
  of scope for v1** — only pointwise.
- **Matching + survival is hard-aborted.** `survival::coxph(..., weights =
  match_weights, cluster = subclass)` on the `MatchIt` output is the correct
  tool and lives outside survatr. Error class `survatr_matching_rejected`.
- **Fine-Gray / subdistribution hazards are out of scope.** Competing risks
  use cause-specific hazards + CIF decomposition only. Document the choice
  explicitly — do not add Fine-Gray as a "missing feature".
- **Track B (longitudinal) per-step link forcing.** Step K: binomial (0/1
  hazard indicator). Steps k < K: quasibinomial (survival-tail
  pseudo-outcome in [0, 1]). Swapping these is a subtle bug — the IF engine
  needs `$mu.eta` and `$variance`, both families provide them, but the
  family switch is load-bearing for the score equations.
- **ICE intervention applies to current-time treatment only.** Lag columns
  hold OBSERVED `A_{k-1}, A_{k-2}, ...` at every backward step (inherited
  rule from causatr ICE). Don't recompute lag columns from the intervened
  treatment.
- **Stochastic + survival MC happens at the cumulative-product level, not
  the hazard level.** `S^g_i(t) = (1/M) sum_m prod_{k <= t} (1 -
  h(k | A_{i,m}, L_{i,k}))`. Averaging hazards across draws before
  cumulating is biased.

### Invariants enforced by code — tests must exercise these, not flag them

- **`na.action = na.exclude` is REJECTED** (error class
  `survatr_bad_na_action`). Only `na.omit` and `na.fail` are accepted.
  Inherited rationale from causatr: residuals padding vs `model.matrix`
  dropping causes silent IF corruption.
- **`censoring =` is a row filter, not IPCW.** For Track A without IPCW, the
  hazard model is fit on uncensored rows; `contrast()` predicts over all
  individuals via the cumulative product. Built-in IPCW (per-period
  cumulative weights on the person-period grid) is a separate path and is
  the motivating use case for IPCW — not an afterthought.
- **Matching is binary-only in causatr; survival + matching is rejected
  regardless of treatment type.** Single unified abort upstream with error
  class `survatr_matching_rejected`.
- **Competing risks: `competing = ` cannot silently fit a cause-deleted
  hazard.** The pre-removal `causat_survival(competing = ...)` in causatr
  aborted for this reason. In survatr, passing `competing = <col>` without
  invoking the cause-specific + CIF path is an error
  (`survatr_competing_misuse`).

### Implementation conventions

- **causatr is the engine.** Import `prepare_model_if()`,
  `apply_model_correction()`, `vcov_from_if()`,
  `variance_if_numeric()` via `causatr:::` (or re-exported internals) and
  layer the cross-time delta aggregation on top. Do NOT reimplement the
  IF primitives.
- **`to_person_period()` stays in causatr.** survatr calls
  `causatr::to_person_period()` — do not fork it.
- **Interventions are constructed via causatr.** `causatr::static()`,
  `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`. Pass them
  to survatr's fit / contrast API unchanged.
- **`model_fn` parameter** — user passes the fitting function
  (`stats::glm` default, `mgcv::gam`, ...). The hazard model defaults to
  `stats::glm` with `family = binomial()` (or `quasibinomial()` for
  weighted fits — the T6 family switch from causatr's pre-removal
  `causat_survival`).
- **`estimator`, not `method`.** Same rule as causatr — `method` is
  reserved for `...`-forwarding.
- **Two-step API.** Fit the hazard model once, then contrast many
  interventions. Curve-shape results live in the `contrast()` return, not
  the fit.
- **Bootstrap resamples individuals (all person-period rows together), not
  rows.** Refit the hazard model per replicate.

### Survival-curve result shape (S3)

- `contrast()` returns a `survatr_result` with:
  - `estimates` — `intervention | time | s_hat | se_s | ...`
  - `contrasts` — `contrast | time | estimate | se | ...`
  - `time_grid` — numeric vector of time points
- `print` / `plot` / `tidy` dispatch on the time-indexed shape. Forest plots
  via `forrest` are available at a user-chosen reference time `t*`.
- Do NOT try to coerce survival results into the scalar `causatr_result`
  shape — they don't fit.

### Review-time heuristics

- **Before flagging a variance bug**, run sandwich vs bootstrap numerically
  on a simulation with known truth. Cross-time delta derivations are subtle;
  a paper derivation without a numerical check is not a reproducer.
- **Before flagging a survival-curve bug**, cross-check against
  `lmtp::lmtp_tmle(outcome_type = "survival")` for the point estimate
  (EIF-based SE is not directly comparable to the M-estimation sandwich, so
  it is a point-estimate oracle only).
- **Pooled logistic vs Cox.** For small per-interval hazards (< 0.1) the two
  agree to high precision (Hernán & Robins TP 17.1; D'Agostino et al.
  1990). Don't flag disagreement at larger hazards as a bug without
  checking the grid spacing.
