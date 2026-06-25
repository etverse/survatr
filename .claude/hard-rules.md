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
| **Estimator** | gcomp (pooled-logistic), ipw (weighted hazard MSM), ice (hazard pseudo-outcome), aipw (parametric doubly-robust; ML/TMLE out). **Matching: hard-abort.** |
| **Treatment timing** | point (point-treatment g-computation), longitudinal (longitudinal ICE-hazard) |
| **Treatment type** | binary, continuous, categorical (k>2), count (Poisson/NB, IPW only), multivariate (via causatr inheritance) |
| **Outcome family** | binomial hazard (first-step / indicator), quasibinomial (pseudo-outcome / weighted fits) |
| **Model class** | GLM (pooled logistic), GAM (via `mgcv::gam` with `s(t)` for baseline hazard) |
| **Intervention** | static, shift, scale_by, threshold (gcomp only), dynamic, ipsi (IPW only), stochastic (pending) |
| **Estimand** | survival S^a(t), risk 1 - S^a(t), risk difference, risk ratio, RMST + RMTL up to t*, survival quantiles / median, per-cause years-of-life-lost |
| **Contrast type** | difference, ratio |
| **Variance method** | sandwich (delta-method cross-time IF; pointwise + simultaneous bands; cluster-robust), bootstrap (resample individuals), numeric Tier 1/2 fallback |
| **Weights** | none, survey/external, censoring row-filter, IPCW (cumulative, per-period) |
| **Entry / truncation** | right-censoring, left-truncation / delayed entry |
| **Event structure** | single terminal event, competing risks, recurrent events, multi-state / illness-death |
| **Competing risks** | cause-specific hazards + CIF (first-class); Fine-Gray / subdistribution hazards out of scope |

## Hard rules (appended to the skill's generic rules)

### Architecture invariants — DO NOT flag these as bugs without a numerical reproducer

- **Pooled-logistic hazard on person-period rows** is the outcome model for
  all of point-treatment g-computation. `logit h(t | A, L) = alpha(t) + beta_A A + beta_L L`. The
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
- **Longitudinal ICE-hazard per-step link forcing.** Step K: binomial (0/1
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
- **`censoring =` is a row filter, not IPCW.** For point-treatment g-computation without IPCW, the
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

### Established invariants from 2026-04-22 round-1 critical review

Do NOT flag these as bugs. Each has a classed-error boundary check
and a regression test:

- **NA in predictor columns is rejected upfront** (`survatr_na_in_predictors`).
  the point-treatment g-computation contrast path predicts on every PP row; NA would either
  silently propagate through `cumprod` or be dropped by `glm`'s
  `na.omit` (our permitted default), misaligning `prep$X_fit` with
  `fit_idx`. NA in the `censoring` column retains "uncensored"
  semantics and is still accepted.
- **Rectangular PP is required** (`survatr_ragged_pp`). Every id must
  have a row at every unique time. Ragged PP (post-event / post-censor
  rows dropped) is rejected at `prepare_pp_data()`. Users pad with
  `outcome = 0` and `censoring = 1` on the padded rows so the
  risk-set builder drops them from the fit. Single-row ids in a
  `n_times == 1` study are trivially rectangular and accepted.
- **Outcome and censoring columns are 0/1 indicators**
  (`survatr_bad_indicator`). Non-binary values walk the cumsum risk-set
  machinery without error.
- **`rmst_weights()` diagonal entries are `dt[j] / 2`, not
  `dt[j-1]/2 + dt[j]/2`.** The `S(0) = 1` contribution is the constant
  `dt[1] / 2` that drops out of the delta and is NOT in the matrix.
  Contract: `W %*% s_hat + dt[1]/2 == trapezoidal_rmst(times, s_hat)`.
- **`print.survatr_result` must not use `show[seq_len(n)]`.** The
  local `n` arg collides with the `n` column of `estimates` under
  data.table NSE. Use `utils::head()` or a non-colliding name.
- **Bootstrap reproducibility under parallel backends requires
  L'Ecuyer-CMRG.** The serial-path `set.seed(seed)` is ignored by
  `mclapply` / `parLapply` otherwise. `bootstrap_survival()` sets
  L'Ecuyer-CMRG before `set.seed()` when parallel != "no" and
  restores the caller's prior RNGkind on exit.
- **Reserved columns include `.cf_hazard` and `.cf_surv`** in addition
  to `.survatr_prev_event` / `.survatr_prev_cens`. The former live only
  on internal copies but are still guarded at the boundary.
- **`validate_times()` accepts numeric + Date + POSIXct + POSIXlt +
  difftime.** `is.numeric()` alone rejects all time-like classes.

### Established invariants from 2026-06-02 round-2 critical review

Do NOT flag these as bugs. Each has a regression test.

- **GAM × sandwich is supported and uses the lpmatrix basis + `Vp` bread.**
  The counterfactual design for an `mgcv::gam` fit MUST be the
  `predict(type = "lpmatrix")` basis, obtained via
  `causatr:::iv_design_matrix()` — NOT `model.matrix(terms(model))`, which
  silently degrades a smooth `s(t)` term to a linear `t` (fewer columns) and
  is non-conformable with `B_inv = model$Vp`. `predict.gam()` also returns a
  1-D array (a vector carrying a `dim` attribute); coerce the link / response
  predictions to plain numeric before the row-scale `X_pp * s_row` or it
  aborts with "non-conformable arrays". The `Vp`-as-bread strategy is mgcv's
  own default and is justified for frequentist coverage by Marra & Wood
  (2012, Scand. J. Stat. 39:53-74). Do not re-mark GAM × sandwich as
  unsupported.
- **survatr performs NO matrix inversion of its own.** Every bread inversion
  is causatr's `prepare_model_if()` / `bread_inv()` (already `MASS::ginv()`-
  guarded, with `$Vp` for GAMs). survatr's variance code is pure matrix
  multiplication plus a guarded risk-ratio division. Do not add a
  singular-bread guard to survatr's delta chain — there is no `solve()` there.
- **`forrest()` has two distinct `t_ref` failure classes:** off-grid /
  malformed `t_ref` → `survatr_bad_t_ref`; a valid grid time with no pairwise
  contrast rows (empty `contrasts`, e.g. a single intervention) →
  `survatr_forrest_no_contrasts`.
- **Symmetric Wald CIs are intentional.** `point ± z·se` for differences
  (risk / RMST) can legitimately produce negative bounds — a difference is
  not a probability. Risk-ratio CIs are built on the log scale and
  exponentiated. A `survival` / `risk` probability CI falling outside [0, 1]
  is a known property of Wald intervals, not a bug; a cloglog / logit
  transform is a possible future enhancement, not a correctness fix.
- **The three `causatr:::` internals are contract-pinned** by
  `test-causatr-integration.R`: `apply_intervention(data, treatment, iv)`;
  `prepare_model_if(model, fit_idx, n_total)` returning `X_fit` / `B_inv` /
  `r_score`; `iv_design_matrix(model, newdata)` returning a coef-aligned
  design. The causatr remote is unpinned `main`; hard-pinning is deferred to
  release time.

### Established invariants from 2026-06-02 round-3 critical review (IPW, chunk 5)

Do NOT flag these as bugs. Each has a regression test.

- **IPW stabilized weights are composed in survatr, not causatr.** causatr's
  point IPW is Hájek-on-`Y~1` and rejects univariate stabilization
  (`causatr_stabilize_univariate`), so it exposes no point-treatment
  `f(A)/f(A|L)` weight. survatr composes it from two
  `causatr:::evaluate_density()` calls (full `A~L` + marginal `A~1`) and
  `causatr:::truncate_weights()`. This is composition of primitives, not
  reimplementation — do not "move it into causatr".
- **The IPW sandwich treatment-correction carries exactly ONE factor of
  `n_ids`.** It is built via `numDeriv` on the weighted-MSM `phi_bar` (÷ n_fit)
  and `causatr:::apply_model_correction(prep_trt, g)` with `prep_trt` at
  `n_total = n_ids`; the `n_fit` in `h_t = n_fit·B_inv·J̄` cancels `phi_bar`'s
  `1/n_fit`. Validated two ways: stacked sandwich SE ≈ full two-stage bootstrap,
  and ≈ an independent `delicatessen` stacked-EE sandwich to ~1e-4 on shared
  data (`test-ipw-delicatessen.R`). Do not add or remove an `n_ids` factor
  without re-running those pins.
- **For stabilized weights the stacked SE is NARROWER than the naive
  weights-as-known SE.** The treatment-model correction is subtracted (it
  projects out the variance explained by the propensity score). Do not "fix"
  the sandwich to be wider than the naive hazard-only SE — that direction is
  correct (Robins 1999; Hernán, Brumback & Robins 2000; Kostouraki 2024). The
  marginal-numerator (`A~1`) estimation is intentionally omitted from the IF
  (conservative).
- **`risk_ratio` `se` is reported on the NATURAL scale** (`RR · se(log RR)`)
  in both the sandwich and bootstrap paths, while the CI is built on the log
  scale and exponentiated. The `se` column is therefore intentionally NOT
  consistent with `point ± z·se` for ratios — the log-based CI is the correct
  interval, and `se` is the natural-scale SE of the estimand. Do not flag the
  `se`-vs-CI asymmetry as a bug.
- **Missing data is rejected upfront for IPW too.** The IPW treatment-model
  predictors (`all.vars(confounders)`) flow through
  `check_no_na_in_predictors()` before dispatch, so NA is rejected
  (`survatr_na_in_predictors`) rather than silently row-dropped by causatr's
  treatment-model fit — this preserves the influence-function row alignment.
- **A constant treatment aborts with `survatr_ipw_no_treatment_variation`**
  (positivity), guarded before `fit_treatment_model()` so causatr's family
  auto-detection cannot misclassify a degenerate binary column as "gaussian".

### Established invariants from 2026-06-03 longitudinal ICE-hazard

Do NOT flag these as bugs. Each has a regression test.

- **survatr owns the ICE backward loop and the cross-step IF cascade; it does
  NOT call `causatr:::ice_iterate()` or `variance_if_ice_one()`.** Those are
  built for a single terminal `Y` (backward step regresses the raw next-step
  prediction `q_{k+1}`) and have no failure carry-forward. The survival
  pseudo-outcome is `Ỹ_k = D_k + (1−D_k)q_{k+1}`, coded NA-safe as
  `ifelse(D_k==1, 1, q_{k+1})` (`ice_surv_pseudo()`). This is reuse of causatr's
  **single-model** primitives composed into a survival recursion, not a
  reimplementation to "move into causatr".
- **The `(1−D_k)` failure carry-forward factor in the survival IF chain is
  REQUIRED and is NOT redundant with the at-risk fit set.** `fit_ids[[k]]` is
  at-risk-at-`k`, which *includes* individuals with the event at `k`
  (`D_k = 1`). The survival cross-step derivative
  `dỸ_k/dβ_{k+1} = (1−D_k) dq_{k+1}/dβ_{k+1}` must zero those out, so
  `survatr_ice_surv_chain()` multiplies each previous-step contribution by
  `(1 − D_{k−1})` from `build_event_by_step()`. Removing it inflates the
  sandwich, growing with horizon. Pinned by the sandwich-vs-bootstrap test and,
  to ~1e-5, by the independent `delicatessen` stacked-EE oracle
  (`test-ice-survival-delicatessen.R`, `data-raw/delicatessen_ice_survival.py`).
  Do not "simplify" the chain back to causatr's verbatim form.
- **Per-step link forcing: binomial at the horizon step, quasibinomial at
  earlier steps.** Asserted in `test-ice-survival.R`. Swapping changes the score
  equations.
- **One backward pass per horizon.** The survival-tail target is relative to a
  fixed final period, so `compute_ice_survival_curve()` runs an independent pass
  per requested `t`; `R^d(t) = mean(pseudo_final)`, `S^d(t) = 1 − R^d(t)`. Do
  not "optimise" to a single pass without re-deriving the intermediate-horizon
  pseudo-outcomes.
- **Lag columns hold OBSERVED treatment; the intervention sets current-period
  treatment only** (`causatr:::ice_apply_intervention_long()`). The static
  counterfactual is recovered by the backward composition, not by setting lags.
  Do not recompute lags from the intervened treatment.
- **longitudinal ICE-hazard confounders split: `confounders` (baseline, never lagged) +
  `confounders_tv` (lag-expanded).** `history` is the Markov order (`Inf` =
  full). Both longitudinal-ICE-only; `time_formula` is NOT part of the ICE per-step
  formula. The minimal `causatr_fit` for the IF is hand-built via
  `causatr:::new_causatr_fit()`, never `fit_ice()`.
- **longitudinal ICE-hazard standardises over the AT-RISK-AT-BASELINE population; the means use
  `na.rm` and the IF `target` is the at-risk-at-baseline mask.** Individuals
  censored at entry (period 1) are never in the period-1 risk set, so they carry
  `NA` in `pseudo_final`. They drop from the ICE standardisation (consistent
  under MCAR entry censoring). `n_ids` passed to `fill_sandwich_ses` stays the
  FULL first-period count — `ice_if_setup` scales Channel 1 by `n / n_target`,
  so `crossprod(IF) / n_ids^2` is the correct mean variance. Do NOT "fix" this
  to drop the `na.rm` or set `target = rep(TRUE, n)` (that re-introduces the
  all-`NA`-curve bug, Issue #1 of the 2026-06-03 review) and do NOT rescale by
  `n_target` (double-counts the `n/n_target` factor already in the IF).
- **longitudinal ICE-hazard v1 rejects external/IPCW weights (`survatr_ice_external_weights`)
  and `ipsi()` / stochastic interventions
  (`survatr_ice_intervention_deferred`).** A constant-within-id treatment under
  `estimator="ice"` informs (`survatr_ice_static_treatment`) but proceeds; a
  time-varying treatment under point-treatment g-computation warns
  (`survatr_tv_treatment_track_a`). These are deliberate scoping, not bugs.
- **longitudinal ICE-hazard requires a NUMERIC treatment** (binary, or a numeric dose modelled
  linearly). Factor / categorical (k > 2) treatments are rejected with
  `survatr_ice_treatment_unsupported` (the intervention sets a numeric value,
  which collides with a factor column; a `treatment_form` design path ships
  later). A numeric-coded multi-level treatment is modelled linearly by design
  — that is the user's modelling choice, not a bug to "fix".

### Established invariants from 2026-06-08 competing risks (chunk 7)

Do NOT flag these as bugs. Each has a regression test.

- **One shared all-cause risk set; J cause-specific models fit on it.**
  `surv_fit(competing = <col>)` derives `1{competing ≠ 0}`, calls
  `build_risk_set()` ONCE, and fits each cause-`j` model with response
  `1{competing == j}` on those rows (`fit_competing_risks()`). Sharing the risk
  set is exactly "treat other causes as censored at their event time" — do NOT
  fit per-cause risk sets, and do NOT cumsum the multi-valued cause column in
  `build_risk_set` (it expects a 0/1 indicator; `.survatr_any_event` is the
  derived all-cause one). The cause column is validated by
  `check_competing_col()` (non-negative integers, ≥ 1 positive cause), the 0/1
  `check_indicator_col(outcome)` is SKIPPED for the CR path.
- **The CIF sandwich sensitivity denominator is the all-cause `(1 − H)`, NOT
  `(1 − h^(j))`.** `s^(j')_l = h^(j')(1−h^(j'))/(1 − H_l)` with `H = Σ_j h^(j)`.
  The single-event `s = h` cancellation (chunk 3) does not hold because the
  survival uses the summed hazard; the formula reduces to chunk 3 only when
  `J = 1` (verified). Do not "simplify" it back to `s = h`.
- **The CIF derivative uses the LAGGED cumulative sensitivity.** `F^(j)_i(t) =
  Σ_{k≤t} S(k−1) h^(j)_k`, so `∂F/∂β_{j'}` weights `cumSX^(j')(k−1)` — coded as
  `cumSX_lag = cumSX − SX`. The all-cause survival IF uses the *inclusive*
  `cumSX(t)` (it reduces to chunk 3). Mixing the two is a subtle bug.
- **The bread is block-diagonal across causes; the cross-cause correlation
  enters via the summed per-individual IF.** Each cause's
  `causatr:::prepare_model_if()` is independent (a cause model's score involves
  only its own `β_j`); `crossprod(IF)/n²` on the per-individual IF (summed over
  causes) captures the within-id cross-cause correlation. Do not try to build a
  joint cross-cause bread. survatr inverts no matrix of its own here either.
- **Validated to ~1e-4 against `delicatessen` and to ~2% against the bootstrap.**
  `data-raw/delicatessen_competing_risks.py` reads the shared
  `fixtures/python/cr_survival_data.csv`; `test-competing-risks-sandwich.R`
  pins per-cause `F^(j),a(t)`, all-cause `S^a(t)`, and `RD^(j)(t)` (point + SE).
  Do not change the IF chain without re-running these pins.
- **gcomp / point-treatment g-computation only; documented rejections.** CR × ipw/ice →
  `survatr_competing_estimator`; `outcome ≠ competing` or `< 2` causes →
  `survatr_competing_misuse`; bad cause column → `survatr_bad_competing`; CIF
  estimand on a single-event fit (or single-event contrast on a CR fit) →
  `survatr_competing_type`; unknown `cause` → `survatr_bad_cause`; external
  `weights` + `competing` → `survatr_competing_weights`. These are deliberate
  scoping, not bugs. **Fine–Gray / subdistribution hazards stay out of scope.**
- **Cause-specific CIF contrasts ship a truncation-by-death caveat.** A one-time
  `rlang::inform(.frequency = "once")` at compute time + a `print` note for
  `cif_difference` / `cif_ratio` + the vignette. Do not strip the caveat to
  "clean up" the output — numbers must never be emitted silently.
- **The CR IF per-time pulls are guarded by `cr_rows_by_time()`** (one row per id
  per requested time, else `survatr_if_failed`), mirroring the single-event
  `compute_survival_if_matrix()` check. It is unreachable on point-treatment g-computation (rectangular
  PP) so it is NOT dead code to delete — it is the boundary a future ragged-PP /
  left-truncation chunk must keep. (2026-06-08 critical review, Issue #1.)

### Established invariants from 2026-06-10 IPCW (chunk 11)

Do NOT flag these as bugs. Each has a regression test.

- **Censoring-model fit rows are `prev_event==0 & prev_cens==0`, which includes
  the `C_k=1` row itself.** This is intentionally broader than the hazard-MSM
  fit rows (which also require `C_k==0` via `is_uncensored()`). The censoring
  outcome is observed at the `C_k=1` row, so the model must be fit there. Do
  not "fix" the censoring model to exclude censored rows from its fit set.
- **The `IF_gamma_per_id` scaling is `(n_ids / n_cens_fit) * psi_per_id %*%
  B_inv_cens`.** The censoring model is fit on `n_cens_fit` person-period rows
  (multiple per id), but the per-individual IF is indexed at the person level
  (n_ids). `causatr:::prepare_model_if(..., n_total = n_cens_fit)` returns
  `B_inv = n_cens_fit * (X'WX)^{-1}`; multiplying by `n_ids / n_cens_fit` gives
  `n_ids * (X'WX)^{-1} psi_per_id`, the correct per-person M-estimation IF.
  This is validated to <2% by the delicatessen oracle
  (`test-ipcw-delicatessen.R`, `data-raw/delicatessen_ipcw_survival.py`).
  Do not change the `n_ids / n_cens_fit` factor without re-running that pin.
- **The IPCW running product uses ALL person-period rows, not just censoring-
  model fit rows.** `compute_ipcw_running_weights()` computes per-row factors
  and calls `ipcw_running_cumprod()` on the full PP grid. The at-risk restriction
  applies only to which rows enter the MSM score — the product accumulates over
  all periods whether at-risk or not. Do not filter the product to cens_at_risk
  rows; the delicatessen oracle uses the same convention.
- **The numDeriv `phi_bar_cens` closure uses fixed trim thresholds.** Re-
  quantiling inside the closure would make the weight function non-smooth in γ,
  breaking the numerical Jacobian. The thresholds are computed once at the point
  estimate and stored in `fit$ipcw_trim_thresholds`. Do not move the quantile
  computation inside the closure.
- **The censoring correction can widen OR narrow the SE** (direction depends on
  the covariance between the censoring IF and the existing IF matrix). Unlike the
  treatment-model correction (which always narrows SE for stabilized weights), the
  censoring correction has no guaranteed sign. Do not flag a narrowed SE as a bug.
  Validated: IF matrices differ by up to ~0.000159 on the standard DGP; SE impact
  is ~1e-5 relative at moderate censoring (~8% per period).
- **Bootstrap refits both the treatment model AND the censoring model per
  replicate** (survatr stores `fit$ipcw` and `fit$censoring_model_fn` and passes
  them back to `surv_fit()` in the bootstrap loop). If either refit is absent,
  the bootstrap SE will disagree with the sandwich at the 15% tolerance level.

### Established invariants from 2026-06-11/12 estimands (chunk 12)

Do NOT flag these as bugs. Each has a regression test.

- **Estimand dispatch lives in the central registry (`R/estimand_registry.R`),
  not scattered `switch(type)` / `type %in% c(...)` lists.** Every dispatch site
  (`build_contrasts`, `fill_sandwich_ses`, `boot_estimand_col` + the bootstrap
  `has_contrast`, `plot` / `tidy` / `print` / `forrest`, the validators, the CR
  fillers) looks up `estimand_field()` / `estimand_is_contrast()` /
  `estimand_is_curve()` / `estimand_has_cause()` / `estimand_requires_pair()`.
  Adding an estimand is one registry row plus a new SE family only when the math
  is genuinely new. Do NOT re-introduce hardcoded type lists.
- **RMTL is `t - RMST`; the SE is the IDENTICAL `rmst_weights()` quadratic form**
  (the constant `t*` has zero gradient, so `Var(RMTL) = Var(RMST)`).
  `rmtl_difference = -(rmst_difference)` (same SE, sign cancels). Do not give
  RMTL its own SE path or flag the equal SEs as a bug.
- **The quantile is `q`-indexed, not time-indexed, and tolerates a single
  intervention.** `tau_q = inf{t : S ≤ 1 − q}` via linear interpolation; the SE
  is the implicit-function delta `IF_tau = −IF_S(tau_q)/slope` (`solve_quantile`
  / `quantile_if_vector`), reading the IF at `tau_q` by interpolating the
  bracketing IF columns (origin bracket `lo == 0` → zero IF column; `slope == 0`
  is dead defensive code, unreachable when `j ≥ 2` gives `s_lo > p ≥ s_hi`). The
  result has its OWN assembly (`assemble_quantile_result` / `finalize_quantile`),
  NOT `fill_sandwich_ses`; `estimand_requires_pair("quantile")` is FALSE so a
  lone median is accepted. Any GENERIC machinery that assumes a `time` index or
  ≥2 interventions must be made shape-aware: the bootstrap minor index is
  `meta$index_col` ∈ {`time`, `q`}; `forrest()` rejects `q`-indexed results with
  `survatr_forrest_wrong_type`; `plot(which = "auto")` falls back to `"curves"`
  when contrast-kind but the contrasts table is empty; and `setdiff(iv, ref)`
  feeding `paste0()` MUST guard the empty case (`paste0` recycles length-0 to
  length-1, creating a phantom contrast). Aborts: `survatr_quantile_unreached`
  (curve never crosses `1 − q`), `survatr_bad_q` (`q ∉ (0, 1)`).
- **YLL is the per-cause CIF integral; `Σ_j YLL^(j)(t*) = RMTL(t*)` to machine
  precision.** `YLL^(j) = ∫ F^(j)` reuses the chunk-7 CIF IF mapped through
  `rmst_weights()` — the CIF starts at `F(0) = 0`, so the quadrature is exactly
  `W %*% cif` and the IF is `CIF_IF %*% t(W)` (no `S(0) = 1` constant). The CR
  sandwich filler shares one `map_if()` keyed on the registry `level_se`
  (`"rmst"` → through `W`, for BOTH all-cause RMST/RMTL and per-cause YLL;
  `"pointwise"` / `"cif"` → direct). YLL needs a CR fit
  (`survatr_yll_needs_cr`, distinct from `survatr_competing_type`). All-cause
  `rmst` / `rmtl` now run on a CR fit (the identity needs it); the extra
  `cause = NA` column survives `add_rmst_to_estimates`'s `setcolorder` (trailing
  columns retained), so the `is.na(cause)` filter holds. Per-cause RMST and CR
  `*_difference` contrasts stay out of scope.

### Established invariants from 2026-06-24 cluster-robust SE (chunk 13)

Do NOT flag these as bugs. Each has a regression test.

- **The cluster-robust divisor is `n²` (number of INDIVIDUALS), NOT `G²`.**
  `contrast(cluster = <col>)` sums each individual's IF row into its cluster's
  row, then `crossprod(rowsum(IF, cluster)) / n²`. The estimator is the
  person-level average `(1/n) Σ_i IF_i`, so its variance is `(1/n²)
  Var(Σ_i IF_i)` and clustering only changes the *meat* (between-cluster outer
  products replace between-individual ones). This matches `sandwich::vcovCL(type
  = "HC0", cadjust = FALSE)` to 1e-10 (`test-variance-cluster.R`) and is the
  divisor `causatr:::vcov_from_if()` already uses. The chunk-13 design doc's
  original `G²` was WRONG (over-inflates by `(n/G)²`, factor 64 in the n=150
  check); the doc has been corrected. Do NOT "fix" the divisor to `G²`.
- **`cluster = id` reproduces the per-individual SE to machine tolerance** (1e-12
  for survival / risk_difference / rmst_difference / quantile / CIF). Singleton
  clusters give `G = n` and `rowsum` is the identity. This is the load-bearing
  regression invariant — if it breaks, the cluster aggregation is wrong.
- **One shared helper does the reduction everywhere.**
  `clustered_pointwise_se(if_mat, n_ids, cluster)` wraps
  `causatr:::vcov_from_if(asplit(if_mat, 2), n_ids, NULL, cluster)` and is called
  by `fill_sandwich_ses()` (gcomp / IPW / IPCW), `fill_sandwich_ses_cr()`
  (competing risks), and `assemble_quantile_result()` (the scalar `n×1` quantile
  IF). With `cluster = NULL` it is exactly the previous inline
  `sqrt(diag(crossprod(IF)) / n_ids²)`. Do NOT re-inline the reduction.
- **Cluster labels are resolved ONCE (name-keyed) then reindexed per path.**
  `validate_cluster()` returns a named (id → label) character vector;
  `cluster_for_ids(labels, unique_ids)` aligns it onto each path's IF row order
  by id NAME (not position), because the single-event / CR / quantile paths
  derive `unique_ids` separately. A missing id aborts `survatr_if_failed`.
- **The contrast SE need not widen even when the level SE does.** With treatment
  randomized WITHIN cluster, a shared cluster effect cancels in the *difference*
  IF (`IF_S_a0 − IF_S_a1`), so the difference SE may be ≈ or slightly below the
  per-individual SE — correct. Widening on the difference requires cluster-level
  treatment (the frailty does not cancel). The level (survival / risk) SE always
  widens under positive within-cluster correlation. Do not flag a non-widening
  difference SE as a bug.
- **The bootstrap resamples whole clusters when `cluster` is set.**
  `bootstrap_survival()` generalizes its resampling unit from id to cluster
  (`unit_to_ids` maps a cluster to all its member ids); a re-drawn cluster gets a
  brand-new set of bootstrap-local ids via the running counter. Cluster bootstrap
  SE ≈ clustered sandwich SE and is wider than the per-individual bootstrap.
- **longitudinal ICE-hazard + cluster is SUPPORTED.** `contrast_track_b()` aligns the
  name-keyed cluster labels onto the ICE IF matrix's FIRST-PERIOD id order
  (`cluster_for_ids()` on `base$data_lag` at `details$time_points[1]`), then
  routes through the same `fill_sandwich_ses()` / `bootstrap_survival()`.
  Entry-censored ids carry a near-zero IF row and sum into their cluster like any
  other — the at-risk-at-baseline standardisation affects the IF values, not the
  row set. Verified by `cluster = id` reproduction (1e-10) and an independent
  within-cluster `rowsum`. Tests MUST use a time-varying-treatment (feedback) DGP
  (`sim_clustered_ice()`); a static treatment makes the ICE chain rank-deficient
  (lag collinear with the current value), independent of clustering — that
  rank-deficiency warning is a DGP artifact, not a cluster bug.
- **Validation aborts:** `survatr_cluster_varies_within_id` (label not constant
  within id), `survatr_cluster_na` (NA in the column), `survatr_cluster_degenerate`
  (`G < 2`), `survatr_bad_cluster` (not a single existing column name). All
  snapshot-pinned.

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
- **`expect_equal(tolerance = t)` uses RELATIVE tolerance (waldo 0.6+).** The
  effective absolute threshold is `t × mean(|target|)`. For survival
  estimates in the 0.07–0.30 range, `tolerance = 0.04` gives an absolute
  threshold of only ~0.007 — not 0.04. For inter-estimator cross-checks
  where the known systematic offset is ~0.02–0.04, use
  `all.equal(..., scale = 1, tolerance = 0.04)` for an absolute comparison
  (verified: `waldo::compare` delegates to `all.equal` and inherits the
  `scale` arg). Emit the `all.equal` result on failure via `testthat::fail()`
  so CI shows the actual deviations. RMST values ~12–14 are less affected:
  `tolerance = 0.05` → ~0.65 absolute, which is usually appropriate for a
  cross-estimator sanity check.

### Once-per-session `rlang::inform` and test ordering

- **`rlang::inform(.frequency = "once", .frequency_id = "...")` fires exactly
  once per R session regardless of how many times the call site executes.**
  Tests that `expect_message(class = "...")` for a once-per-session inform
  will fail if any earlier test file (alphabetical order) already triggered
  that inform, because the rlang once-tracking state persists across test
  files in the same session.
- **Fix: in non-targeted tests, avoid calling code that fires the once-per-
  session inform.** If a test only needs to verify *structural* behaviour
  (e.g. "weights panel is NULL for ICE"), choose a DGP that doesn't trigger
  the inform (e.g. use a time-varying treatment DGP like `sim_ice_feedback`
  instead of a static-treatment DGP for ICE) so the slot stays unconsumed.
- **The targeted test that `expect_message(...)`s the inform should be the
  only test to trigger it.** Currently `survatr_ice_static_treatment` is
  tested in `test-ice-survival.R:289`; no other test file should call
  `surv_fit(estimator = "ice")` with a constant-within-id treatment.
