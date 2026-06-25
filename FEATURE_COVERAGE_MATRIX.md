# Feature coverage matrix

Single source of truth for **what works** in survatr. Mirrors causatr's
`FEATURE_COVERAGE_MATRIX.md` convention: every PR that adds, removes, or
changes a feature MUST update this file and the corresponding tests.

## Legend

- 🟢 — supported, truth-based test pinned against an analytical or external
  reference (`lmtp::lmtp_tmle(outcome_type = "survival")`,
  `gfoRmula::gformula_survival()`, `survival::survfit` /
  `survival::coxph`, or a closed-form DGP).
- 🟡 — smoke-tested only (runs without error; point estimate / SE not
  pinned to a truth). Acceptable for combinations where no oracle exists,
  temporary during a multi-chunk feature rollout, or where the oracle is
  too expensive to run in CI.
- 🔴 — hard-rejected by a classed error with a regression test pinning the
  rejection.

<!-- Rows are added as features ship. No speculative ⚪ entries: the matrix
reflects **current** state, not planned scope. Planned scope lives in
`SURVIVAL_PACKAGE_HANDOFF.md` §10 (implementation chunks). -->

## Track A — Point survival, pooled-logistic hazard

### Fit path (`surv_fit()`)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| `estimator = "gcomp"`, `binomial()` family, unweighted | 🟢 | `test-surv_fit-family-oracle.R` | Closed-form constant-hazard DGP: β₀ → qlogis(h) at n = 5000, K = 10, h = 0.05. Tolerance 0.1 (~1.6 × MC-SE). |
| `estimator = "gcomp"`, `quasibinomial()` family, `weights ≡ 1` | 🟢 | `test-surv_fit-family-oracle.R`, `test-surv_fit-weighted.R` | Same closed-form DGP (intercept oracle); coefficient equivalence vs unweighted `binomial()` to 1e-10. |
| Risk-set construction (drop at/after first event; drop at/after first censor) | 🟢 | `test-risk_set.R` | Hand-verified 5-id, 3-period fixture (`fixture_small_pp()`): n_at_risk = 12 without censoring, 10 with. |
| `is_uncensored()` (NA \| 0 ⇒ at-risk) | 🟢 | `test-risk_set.R` | Hand-specified truth vector. |
| Reserved-column guard (`.survatr_prev_event`, `.survatr_prev_cens`) | 🔴 | `test-checks.R`, `test-surv_fit.R` | `survatr_reserved_col`. |
| External weights validation (NA / Inf / NaN / negative / mis-sized / non-numeric) | 🔴 | `test-checks.R`, `test-surv_fit-weighted.R` | `survatr_bad_weights`. Zero weights allowed. |
| `na.action = na.exclude` | 🔴 | `test-checks.R`, `test-surv_fit.R` | `survatr_bad_na_action`. Inherited rationale from causatr (residuals padding vs `model.matrix` drop misalignment). |
| `estimator = "matching"` / `"match"` | 🔴 | `test-matching-rejection.R`, `test-surv_fit.R` | `survatr_matching_rejected`. Points to `survival::coxph(..., weights = match_weights, cluster = subclass)`. |
| `method = "matching"` / `"match"` in `...` (causatr-style mis-call) | 🔴 | `test-matching-rejection.R` | `survatr_matching_rejected`. Caught before model dispatch. |
| `data` is a `matchit` object (MatchIt output) | 🔴 | `test-matching-rejection.R` | `survatr_matching_rejected`. Detected before column lookup. |
| `estimator = "ipw"` (weighted marginal hazard MSM) | 🟢 | `test-ipw-survival.R` | See the IPW section below. |
| `estimator = "ice"` (longitudinal ICE hazards) | 🟢 | `test-ice-survival.R` | See the Track B section below. `track = "B"`, `model = NULL` (per-step models fit lazily in `contrast()`). |
| `estimator = <unknown>` | 🔴 | `test-surv_fit.R` | `survatr_bad_estimator`. |
| `competing = <non-NULL>` | 🔴 | `test-surv_fit.R` | `survatr_competing_misuse`. Cause-specific + CIF path ships in chunk 7. |
| Missing column name in `data` | 🔴 | `test-prepare_data.R` | `survatr_col_not_found`. |
| Wide input (one row per id) | 🔴 | `test-prepare_data.R` | `survatr_not_person_period`. Points to `causatr::to_person_period()`. |
| Duplicated `(id, time)` rows | 🔴 | `test-prepare_data.R` | `survatr_duplicate_pp_row`. |
| Input mutation safety (`data.frame` / `data.table` both copied) | 🟢 | `test-prepare_data.R`, `test-surv_fit.R` | Post-fit name / row-count equality with pre-fit snapshot. |

### Contrast path (`contrast.survatr_fit()`)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| `type = "survival"` on constant-hazard DGP | 🟢 | `test-contrast-survival.R`, `test-survival_curve.R` | Closed-form `S(t) = (1 - h)^t` on n = 5000, absolute tolerance 0.03 (~4 MC-SE). |
| `type = "risk"` | 🟢 | `test-contrast-survival.R` | Derived from `s_hat`; returns curve-only (empty `contrasts` stub). |
| `type = "risk_difference"` | 🟢 | `test-contrast.R` | DGP with no treatment effect: RD ≈ 0 across time, tolerance 0.02 at n = 5000. |
| `type = "risk_ratio"` | 🟢 | `test-contrast.R` | DGP with no treatment effect: RR ≈ 1, tolerance 0.15. |
| `type = "rmst"` | 🟢 | `test-rmst.R`, `test-contrast.R` | Closed-form trapezoidal integral of `(1-h)^t` matched to 1e-12; curve-only shape verified. |
| RMST cross-check vs `survRM2::rmst2()` (unadjusted KM) | 🟢 | `test-rmst-survRM2.R` | Constant-hazard DGP (n = 3000, K = 20, h = 0.05, no covariates): per-arm pooled-logistic trapezoidal RMST at t = 20 agrees with KM-based RMST from `survRM2::rmst2()` within 0.05. Validates the RMST scale (pooled-logistic ≈ KM for small h). `skip_if_not_installed("survRM2")`. |
| `type = "rmst_difference"` | 🟢 | `test-contrast.R` | DGP with no effect: RMST-diff ≈ 0, tolerance 0.1. |
| `type = "rmtl"` (restricted mean time lost) | 🟢 | `test-estimands-rmtl.R`, `test-gcomp-delicatessen.R` | `RMTL(t) = t - RMST(t)` identity to 1e-12; **direct delicatessen cross-check** (implied oracle `t - RMST_deli` with the identical SE, from the shared gcomp fixture). Wired across gcomp / IPW / IPCW / Track B (shared SE filler); sandwich vs bootstrap SE ratio in (0.7, 1.4). |
| `type = "rmtl_difference"` | 🟢 | `test-estimands-rmtl.R` | Equals `-(rmst_difference)` (point + SE) to 1e-12. |
| RMTL variance identity `Var(RMTL) = Var(RMST)` | 🟢 | `test-estimands-rmtl.R` | Same trapezoidal quadratic form: the constant `t*` drops out of the delta gradient, so the SE is identical to RMST's to 1e-12. |
| `type = "quantile"` (survival quantile / median) | 🟢 | `test-estimands-quantile.R` | Constant-hazard closed form: median `= log(2)/λ`, `τ_q = -log(1-q)/λ` matched to ~0.3% at n = 20000 (estimator returns the exact linear-interp crossing of the fitted curve; verified at n up to 1e5). `survival::survfit()` KM-median sanity. Vector `q` supported. **No delicatessen oracle** — the quantile is a non-smooth functional outside delicatessen's smooth M-estimation; validated by the closed form + KM + sandwich-vs-bootstrap instead. |
| Quantile sandwich (implicit-function delta) vs bootstrap | 🟢 | `test-estimands-quantile.R` | `dτ_q = -dS(τ_q)/S'(τ_q)` via interpolated IF columns; sandwich-vs-bootstrap SE ratio in (0.7, 1.4). |
| Median difference contrast + single-intervention quantile | 🟢 | `test-estimands-quantile.R` | `estimate = τ(a1) - τ(a0)` to 1e-10; a lone median (one intervention) is accepted (no `survatr_bad_interventions`). |
| Quantile across estimators (IPW / IPCW / ICE / CR all-cause) | 🟢 | `test-estimands-quantile.R` | Reuses the per-estimator survival IF; wired for gcomp / IPW / IPCW / Track B and all-cause on a competing-risks fit (no cause dimension). |
| `survatr_quantile_unreached` / `survatr_bad_q` | 🔴 | `test-estimands-quantile.R` | Curve never crosses `1 - q` on the grid; `q` outside `(0, 1)`. |
| `type = "yll"` (per-cause years of life lost) | 🟢 | `test-estimands-yll.R` | `∫F^(j)` matched to the analytic two-cause CIF integral (`analytic_cr()` × trapezoidal weights) within 5% at n = 6000; carries the `cause` dimension. Competing-risks fit only. Oracle: the CIF building block is delicatessen-pinned (`test-competing-risks-sandwich.R`) and YLL inherits via the `Σⱼ YLL = RMTL` identity + the analytic integral (a direct delicatessen YLL would need the full CIF curve, which the selected-time fixture does not store). |
| YLL identity `Σⱼ YLL^(j)(t*) = RMTL(t*)` | 🟢 | `test-estimands-yll.R` | Holds to 1e-10 (partition of unity `Σⱼ F^(j) = 1 - S`, integral linear). |
| YLL sandwich (CIF IF × trapezoidal) vs bootstrap | 🟢 | `test-estimands-yll.R` | CIF IF mapped through `rmst_weights()`; sandwich-vs-bootstrap SE ratio in (0.6, 1.5). |
| All-cause `rmst` / `rmtl` on a competing-risks fit | 🟢 | `test-estimands-yll.R` | Integral of the all-cause survival; needed for the `Σⱼ YLL = RMTL` identity. Per-cause RMST and `*_difference` on CR remain out of scope. |
| `survatr_yll_needs_cr` | 🔴 | `test-estimands-yll.R` | `type = "yll"` on a single-event fit (distinct class from `survatr_competing_type`). |
| Oracle cross-check vs `lmtp::lmtp_tmle(outcome_type = "survival")` | 🟢 | `test-contrast-lmtp-oracle.R` | gcomp `S^a(t)` (a1 and a0) on a confounded DGP (n = 2000) matches lmtp's TMLE survival within 0.05 at t ∈ {3, 5}. lmtp 1.5.3: one fit per horizon, `ife@x` estimate (the prior `folds`/`$theta` form silently skipped). |
| Per-individual cumulative product (Jensen-safe) | 🟢 | `test-survival_curve.R` | Cumulative product within id before averaging across ids; monotone non-increasing in t on random DGPs. |
| RMST trapezoidal quadrature (closed form) | 🟢 | `test-rmst.R` | Explicit sum of `(S(t_i) + S(t_{i+1}))/2 * (t_{i+1} - t_i)` reproduced to 1e-12. |
| `fit$pp_data` mutation safety | 🟢 | `test-contrast.R` | Names / nrow / treatment column identical pre- and post- `contrast()`. |
| Empty / unnamed / duplicate-named interventions | 🔴 | `test-contrast-rejections.R` | `survatr_bad_interventions`. |
| Non-`causatr_intervention` list elements | 🔴 | `test-contrast-rejections.R` | `survatr_bad_interventions`; error message points to `causatr::static()` et al. |
| `times` not numeric / empty / NA / outside `fit$time_grid` | 🔴 | `test-contrast-rejections.R` | `survatr_bad_times` (structural) / `survatr_time_extrapolation` (values outside grid). |
| Bad `reference` (not a name in `interventions`) | 🔴 | `test-contrast-rejections.R` | `survatr_bad_reference`. |
| `ci_method = "bootstrap"` | 🟢 | `test-bootstrap-survival.R` | Percentile and Wald CIs over B replicates; per-id resampling preserves within-id dependence. See the "Bootstrap variance" section below. |
| Unknown `ci_method` / `type` / `conf_level` | 🔴 | `test-contrast-rejections.R` | `survatr_bad_ci_method` / `match.arg` / `survatr_bad_conf_level`. |
| `se` / `ci_lower` / `ci_upper` columns when `ci_method = "none"` | 🟡 | `test-contrast.R` | All `NA_real_` by design (opt-in CI path). |

### Sandwich variance (`ci_method = "sandwich"`, delta-method cross-time IF)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| Per-individual IF matrix for `S^a(t)` (Ch1 + Ch2) | 🟢 | `test-variance_if_survival.R` | Column means ~ 0 across individuals on constant-hazard DGP; sandwich SE within 75% of hand-derived closed-form target (tolerance conservative because closed form ignores at-risk attrition). Also cross-checked: sandwich SE matches empirical SD of `s_hat` across 100 simulation replicates to ~7% at t = 5, n = 3000. |
| Sandwich CI for `s_hat` / `risk_hat` | 🟢 | `test-sandwich-survival.R` | 200-rep coverage simulation at n = 1000, K = 6, h = 0.08 achieves ≥ 88% nominal 95% coverage at t ∈ {1, 5}. Single-seed sanity test at n = 5000 covers truth. |
| Sandwich CI for `risk_difference` | 🟢 | `test-sandwich-risk-difference.R` | 200-rep coverage simulation on a no-effect DGP: RD CI covers 0 at ≥ 88% nominal 95%. |
| Sandwich CI for `risk_ratio` (log-scale) | 🟢 | `test-sandwich-risk-ratio.R` | Single-seed test: CI is strictly positive, `ci_upper > ci_lower`, and covers 1 on a no-effect DGP. |
| Sandwich CI for `rmst` (trapezoidal quadratic form) | 🟢 | `test-sandwich-rmst.R` | SE is non-negative, monotone non-decreasing in `t` (by construction of the cumulative trapezoidal integral of a positive IF), CI bounds finite. |
| Sandwich CI for `rmst_difference` | 🟢 | `test-sandwich-rmst.R` | CI at `t = 10` covers 0 on a no-effect DGP (n = 3000). |
| gcomp sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-gcomp-delicatessen.R` | `S^a(t)`, `RD(t)`, `RR(t)`, `RMST^a(t)`, and `RMST-difference(t)` point + sandwich SE match a Python `delicatessen` stacked-EE oracle to ~1e-3 at t ∈ {2,3,4,5} on a shared confounded fixture (`fixtures/python/ipw_survival_data.csv`, also used by the IPW oracle). Reference: `data-raw/delicatessen_gcomp_survival.py`. Closes the only variance chunk that previously lacked an independent M-estimation cross-check. |
| `conf_level` in (0, 1) | 🟢 | `test-contrast-rejections.R` | Rejects values outside the open interval with `survatr_bad_conf_level`. |
| `model_fn = mgcv::gam` (penalized `s(t)` baseline) × sandwich | 🟢 | `test-sandwich-gam.R` | Counterfactual design built on the gam `lpmatrix` basis via `causatr:::iv_design_matrix()` to match the `model$Vp` bread; `predict.gam` 1-D-array output coerced to plain numeric. GAM sandwich SE matches the analytically-anchored GLM sandwich SE within 2% on a constant-hazard DGP, and tracks the bootstrap SE identically to the GLM. `Vp`-as-bread justified for frequentist coverage by Marra & Wood (2012). A gam fit lacking `$Vp` still aborts in `causatr:::bread_inv()`. |

### Cluster-robust sandwich (`contrast(cluster = "<column>")`, chunk 13)

The per-individual IF rows are summed within cluster before `crossprod`; the
divisor stays `n²` (`V = crossprod(IF_g) / n²`). One shared helper
(`clustered_pointwise_se()` → `causatr:::vcov_from_if(cluster=)`) serves every
sandwich path. `cluster = "<id-column>"` reproduces the per-individual SE.

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| Aggregation primitive vs `sandwich::vcovCL` | 🟢 | `test-variance-cluster.R` | `clustered_pointwise_se()` matches `vcovCL(type = "HC0", cadjust = FALSE)` to 1e-10 on the sample mean; singleton clusters match `vcovHC(type = "HC0")`. Confirms the `n²` divisor (not `G²`). |
| `cluster = id` reproduces per-individual SE | 🟢 | `test-variance-cluster.R` | Machine-tolerance (1e-12) match for survival / risk_difference / rmst_difference / quantile (gcomp) and CIF (competing risks). The load-bearing regression invariant. |
| Clustered SE ≥ per-individual SE (positive within-cluster corr) | 🟢 | `test-variance-cluster.R` | Multi-site frailty DGP: level (risk) SE widens uniformly (~2.3×); contrast SE widens under cluster-level treatment (~3.4×). |
| Calibration to the cluster-sampling SD | 🟢 | `test-variance-cluster.R` | Skipped on CRAN. Empirical SD of risk@t* over 150 re-draws of the sites: clustered SE within 30%; per-individual SE < 0.85× the truth (under-states). |
| Cluster-resampling bootstrap | 🟢 | `test-variance-cluster.R` | Skipped on CRAN. Resamples whole sites (B = 400); SE ≈ clustered sandwich within 25%, wider than the per-individual bootstrap. |
| gcomp / IPW / IPCW / competing-risks / quantile coverage | 🟢 | `test-variance-cluster.R` | All flow through the shared helper (single-event `fill_sandwich_ses()` + CR `fill_sandwich_ses_cr()` + `assemble_quantile_result()`). IPW/IPCW non-singleton meat pinned against the cluster bootstrap. |
| Track B (ICE) cluster-robust SE | 🟢 | `test-variance-cluster.R` | IF rows aligned by first-period id; `cluster = id` reproduces the per-individual SE to 1e-10 and the clustered SE matches an independent within-cluster `rowsum`; widens under frailty. Time-varying-treatment (feedback) DGP keeps the ICE chain full-rank. |
| Validation aborts | 🔴 | `test-variance-cluster.R` | Snapshot-pinned: `survatr_cluster_varies_within_id`, `survatr_cluster_na`, `survatr_cluster_degenerate`, `survatr_bad_cluster`. |

### Bootstrap variance (`ci_method = "bootstrap"`, resample individuals)

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| Per-id resampling + per-replicate refit | 🟢 | `test-bootstrap-survival.R` | Smoke: CIs populated, point estimate in `[ci_lower, ci_upper]`. Cluster = id; each replicate draws n_ids ids with replacement, concatenates their PP blocks (renumbered), refits via `surv_fit()`, contrasts via `contrast(ci_method = "none")`. |
| Percentile CI | 🟢 | `test-bootstrap-survival.R` | Single-seed coverage of `(1-h)^t` at n = 2000, B = 300. Default `boot_ci = "percentile"` (transform-invariant; safer for ratios / RMST). |
| Wald CI | 🟢 | `test-bootstrap-survival.R` | Sample-SD × `z` bands around the observed point estimate. |
| Reproducibility with `seed` | 🟢 | `test-bootstrap-survival.R` | Two calls at the same `seed` return identical SEs and CI endpoints. |
| Bootstrap SE ≈ sandwich SE (cross-check) | 🟢 | `test-bootstrap-survival.R` | Skipped on CRAN. B = 500 at n = 1500: per-time SE agrees within 15%. |
| Empirical-SD oracle: sandwich + bootstrap ≈ sampling SD | 🟢 | `test-bootstrap-survival.R` | Skipped on CRAN. 100-replicate sim truth (n = 1000, h = 0.06, K = 6). Out-of-band 300-rep validation (2026-04-22) pinned both to within 1-2% of truth; the in-test 20% tolerance catches class-of-factor-n_ids scaling bugs like the one fixed in chunk 3 (`a3f79cb`). |
| `risk_ratio` via percentile CI | 🟢 | `test-bootstrap-survival.R` | Strictly positive; covers 1 on a no-effect DGP. |
| `risk_difference` with populated contrast CIs | 🟢 | `test-bootstrap-survival.R` | `contrasts$se` / `ci_*` non-NA, point ∈ CI. |
| Failure guard (>10% replicate failures) | 🔴 | — | `survatr_boot_failed`; point the user at sandwich or at a smaller / simpler DGP. |
| Parallel backend (`parallel` + `ncpus`) | 🟢 | `test-bootstrap-rejections.R` | Validated upfront; accepts `"no"`, `"multicore"`, `"snow"`; `ncpus` positive integer. |
| Bad `n_boot` / `boot_ci` / `parallel` | 🔴 | `test-bootstrap-rejections.R` | `survatr_bad_n_boot`, `survatr_bad_boot_ci`, `survatr_bad_parallel`. |

### S3 methods on `survatr_result`

| Surface | Status | Test file | Notes |
|---|---|---|---|
| `print()` | 🟢 | `test-surv_fit.R`, `test-contrast.R` | Snapshot-pinned. Shows type, reference, ci_method, time grid, head of contrasts (or estimates for curve-only). |
| `tidy()` | 🟢 | `test-tidy-survatr_result.R` | Long `data.frame` with `intervention`, `contrast`, `time`, `estimand`, `estimate`, `se`, `ci_lower`, `ci_upper`. `which` in `{"all", "estimates", "contrasts"}`; `conf.int = FALSE` drops CI columns. S3 method on the `generics::tidy` generic (re-exported). |
| `plot()` | 🟢 | `test-plot-survatr_result.R` | Base-R graphics: curves for `survival` / `risk` / `rmst`, contrasts with reference line at 0 / 1 for the three pairwise types. CI ribbons via `adjustcolor` when populated. Smoke-only (no `vdiffr`). |
| `forrest()` | 🟢 | `test-forrest-survatr_result.R` | Forest plot at a user-chosen `t_ref`. Aborts on curve-only types (`survatr_forrest_wrong_type`), on `t_ref` outside `time_grid` (`survatr_bad_t_ref`), and on a valid grid time with no pairwise contrast rows (`survatr_forrest_no_contrasts`). |
| causatr internal-API contract (`apply_intervention`, `prepare_model_if`, `iv_design_matrix`) | 🟢 | `test-causatr-integration.R` | Pins the formal names + return shape of the three `causatr:::` internals survatr depends on, so an upstream signature drift fails loudly in CI rather than as a cryptic runtime error. |
### IPW — weighted marginal hazard MSM (`estimator = "ipw"`, binary treatment)

Fit a baseline treatment model, form **stabilized** density-ratio weights
`w_i = f(A_i) / f(A_i | L_i)` at the observed treatment, broadcast onto the
person-period rows, and fit a weighted marginal MSM `logit h(t|A) = α(t) + β_A·A`
(quasibinomial). Confounding is handled by the weights, not the hazard model.

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| IPW survival curve `S^a(t)`, binary `static` interventions | 🟢 | `test-ipw-survival.R` | `lmtp::lmtp_tmle(outcome_type = "survival")` on a confounded DGP (n = 2000): per-arm S(t) at t ∈ {3, 5} within 0.05. Mid-test agreement ~0.01. |
| Degenerate single-period IPW ⇒ scalar point IPW | 🟢 | `test-ipw-survival.R` | `causatr::causat(estimator = "ipw")` on a one-period PP table: risk under each arm matches to 1e-2 (saturated binary MSM ⇔ Hájek point IPW). |
| Unconfounded (γ = 0) ⇒ unadjusted gcomp curve; weights ≈ 1 | 🟢 | `test-ipw-survival.R` | IPW curve matches the `confounders = ~1` gcomp curve to 0.02; mean weight ≈ 1, sd < 0.25. |
| Stabilized weight constant within id (broadcast invariant) | 🟢 | `test-ipw-survival.R` | One unique weight value per id; `fit$weights` length = `nrow(data)`. |
| `trim` winsorization (fixed cutoff for the sandwich) | 🟢 | `test-ipw-survival.R` | `trim = 0.9` caps the per-id weights at the resolved quantile; untrimmed records `NA` threshold. |
| Collinear confounder (aliased propensity coef dropped) | 🟢 | `test-ipw-survival.R` | causatr drops the aliased coef (warning); IPW fit + sandwich finite. |
| IPW stacked-EE sandwich (hazard block + treatment-model correction) | 🟢 | `test-ipw-survival-sandwich.R`, `test-ipw-delicatessen.R` | Stacked SE matches the full two-stage bootstrap SE (resamples ids, refits both models) within 15%; the treatment-model correction is strictly variance-reducing vs the naive weights-as-known SE (Robins 1999; Hernán et al. 2000). |
| IPW sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-ipw-delicatessen.R` | `S^a(t)`, `S^0(t)`, and `RD(t)` point + sandwich SE match a Python `delicatessen` stacked-EE sandwich to ~1e-4 on a shared person-period fixture. Reference generated by `data-raw/delicatessen_ipw_survival.py`; both read `fixtures/python/ipw_survival_data.csv`. |
| IPW sandwich RD CI coverage | 🟢 | `test-ipw-survival-sandwich.R` | 200-rep simulation, n = 800: RD CI covers the marginal-RD truth (∫ over L) at ≥ 88% nominal 95% at t ∈ {2, 4}. |
| IPW bootstrap (per-replicate dual refit) | 🟢 | `test-ipw-survival-bootstrap.R` | Weights re-estimated each replicate (treatment + hazard refit); CIs populated, point ∈ CI, reproducible under `seed`, SE ≈ sandwich. |
| IPW × `mgcv::gam` baseline hazard | 🟢 | `test-ipw-survival-sandwich.R` | lpmatrix-basis + `Vp` bread (as gcomp gam); IPW treatment correction composes via `prep$X_fit`. Stacked sandwich tracks the two-stage bootstrap within 15% on `s(t, k=4)`. |
| Continuous / categorical / count treatment IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_treatment_unsupported` (binary only this chunk; extended types → chunk 19). |
| Time-varying treatment under IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_time_varying_treatment` (point treatment only; longitudinal IPW out of scope). |
| External `weights` + `estimator = "ipw"` | 🔴 | `test-ipw-survival.R` | `survatr_ipw_external_weights` (survey/design-weight composition = transport; deferred → chunk 21). |
| `ipsi()` under IPW | 🔴 | `test-ipw-survival.R` | `survatr_ipw_ipsi_deferred` (weight-path estimand → chunk 20). |
| `trim` validation | 🔴 | `test-checks.R` | `survatr_bad_trim` (NULL or scalar in (0, 1]). |
| causatr IPW internal-API contract (`fit_treatment_model`, `evaluate_density`, `truncate_weights`, `apply_model_correction`) | 🟢 | `test-causatr-integration.R` | Pins formals + return shapes of the IPW `causatr:::` internals survatr reuses. |

### IPCW — per-period censoring weights (`ipcw = <formula>`, `estimator = "ipw"`)

Built-in inverse-probability-of-censoring weighting for the IPW weighted-MSM
path. Fits a per-period censoring hazard model (denominator: `C ~ time + A + L`;
numerator: `C ~ time + A`) to form stabilized running-product IPCW weights
`W^C_{i,k} = ∏_{m≤k} P(C_m=0|A)/P(C_m=0|A,L)`. The combined weight
`w_i × W^C_{i,k}` (treatment × censoring) is passed to the weighted hazard MSM.
Sandwich variance uses a **three-block stacked-EE** approach (beta + alpha + gamma).
Fit-row convention: IPCW model fits on rows where `prev_event == 0 & prev_cens == 0`
(includes the censoring-event row itself). MSM fit rows additionally require
`is_uncensored()` (excludes the C=1 row). Per-period trim stored with fixed
thresholds for the sandwich. Bootstrap refits the censoring model per replicate.

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| IPCW fit structure (censoring model + IPCW weights stored) | 🟢 | `test-ipcw-survival.R` | `fit$censoring_model` class `glm`; `fit$ipcw_weights` non-NULL, length = nrow(data); `fit$weights = fit$ipw_treatment_weights_pp × fit$ipcw_weights`. |
| IPCW weights are per-period running products (grow within id) | 🟢 | `test-ipcw-survival.R` | Within-id weight variance > 0 for most ids; treatment weights constant within id; IPCW weights vary across periods. |
| IPCW curve vs `lmtp::lmtp_tmle(outcome_type = "survival", cens = ...)` | 🟢 | `test-ipcw-survival.R` | Informative-censoring DGP (n = 2000, K = 5, δ = 0.8): per-arm S(t) at t ∈ {3, 5} within 0.05. `skip_if_not_installed("lmtp")`. |
| Naive IPW over-estimates survival under informative censoring (direction test) | 🟢 | `test-ipcw-survival.R` | R = 25 runs × n = 1000: naive IPW > IPCW at t = 5 in > 70% of runs for both arms (δ = 0.8; bias-direction property of informative censoring). `skip_on_cran()`. |
| Non-informative censoring (δ = 0): IPCW weights ≈ 1; curve ≈ naive | 🟢 | `test-ipcw-survival.R` | Mean IPCW weight ≈ 1 (tol 0.05), SD < 0.25; IPCW vs row-filter curve within 0.03. |
| Per-period trim stores fixed thresholds, caps weights | 🟢 | `test-ipcw-survival.R` | `fit$ipcw_trim_thresholds` named by period; all `ipcw_weights[t == k] ≤ threshold[k]`. |
| IPCW contrast structure (finite SE, CIs ordered, point ∈ CI) | 🟢 | `test-ipcw-survival.R` | Sandwich CI on n = 300 DGP: SE > 0 for all arms and RD. |
| Bootstrap refits censoring model per replicate | 🟢 | `test-ipcw-survival.R` | B = 200 bootstrap: SE finite + positive; censoring model re-estimated each replicate. |
| Censoring-model correction wired in the IF matrix | 🟢 | `test-ipcw-survival-sandwich.R` | `A_beta_gamma` and `IF_gamma_per_id` non-zero; three-block IF matrix differs from two-block by > 1e-6 absolute. |
| Treatment-model correction still variance-reducing vs hazard-only | 🟢 | `test-ipcw-survival-sandwich.R` | `se_two < se_hazard` at all times; > 1% relative reduction at peak. |
| Three-block sandwich ≈ full three-stage bootstrap SE | 🟢 | `test-ipcw-survival-sandwich.R` | Per-time RD SE within 15% at B = 400 (n = 700, K = 5). `skip_on_cran()`. |
| IPCW three-block sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-ipcw-delicatessen.R` | `S^1(t)`, `S^0(t)`, and `RD(t)` point + SE match a Python `delicatessen` three-block stacked-EE oracle to ~1e-4 (points) and ~2% (SEs) on a shared informative-censoring fixture. Reference: `data-raw/delicatessen_ipcw_survival.py`; both read `fixtures/python/ipcw_survival_data.csv`. Pins the `A_beta_gamma` cross-derivative and `n_ids / n_cens_fit` bread scaling. |
| IPCW sandwich CI for risk difference covers marginal truth | 🟢 | `test-ipcw-survival-sandwich.R` | 150-rep coverage simulation, n = 800: ≥ 88% nominal 95% at t ∈ {2, 5}. `skip_on_cran()`. |
| `ipcw =` with non-formula | 🔴 | `test-ipcw-survival.R` | `survatr_bad_ipcw`. |
| `ipcw =` with non-`"ipw"` estimator | 🔴 | `test-ipcw-survival.R` | `survatr_ipcw_estimator` (gcomp and ice both rejected; later chunks). |
| `ipcw =` without a `censoring` column | 🔴 | `test-ipcw-survival.R` | `survatr_ipcw_no_censoring`. |
| Snapshot: all three IPCW error messages | 🟢 | `test-ipcw-survival.R` | `_snaps/ipcw-survival.md`. |

### End-to-end acceptance test

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| NHEFS Ch. 17 replication (Track A gcomp, 120-mo survival) | 🟢 | `test-nhefs-replication.R` | H&R 2024 published targets: 120-mo S^a(t) ≈ 80.7% (qsmk=1) / 80.5% (qsmk=0) within ±0.03; RD ≈ 0.2% within ±0.01; sandwich CI spans 0; unadjusted KM in the 75–90% ballpark. Dataset: `nhefs_surv` (1629 × 120 rectangular PP, 318 events). Skipped on CRAN. |

## Track B — Longitudinal survival (ICE hazards)

Time-varying treatment + time-to-event via backward iterated conditional
expectations on the discrete-time hazard, with a survival-tail pseudo-outcome
`Ỹ_k = D_k + (1−D_k)q_{k+1}`. survatr owns the backward loop + cross-step IF
cascade (with the `(1−D_k)` failure carry-forward), reusing causatr's
single-model primitives. Confounders split into `confounders` (baseline) +
`confounders_tv` (lagged); `history` Markov order.

| Feature combination | Status | Test file | Oracle / notes |
|---|---|---|---|
| Track B risk curve `R^d(t)`, static strategies | 🟢 | `test-ice-survival.R` | Forward-simulation Monte-Carlo g-formula truth on a treatment-confounder-feedback DGP (n = 4000, K = 4): `S^a(t)` within 0.02 at every t; protective treatment ⇒ `S^1 > S^0`. |
| Per-step link forcing (binomial at horizon, quasibinomial earlier) | 🟢 | `test-ice-survival.R` | Asserts `models[[K]]$family == "binomial"`, `models[[k<K]] == "quasibinomial"`. |
| Lag columns carry observed (not intervened) treatment | 🟢 | `test-ice-survival.R` | Under `shift(1)`: current `A` shifted, `lag1_A` unchanged. |
| Survival-aware stacked-EE sandwich (`ci_method = "sandwich"`) | 🟢 | `test-ice-survival.R`, `test-ice-survival-delicatessen.R` | SEs finite, CIs bracket the point; `risk_difference` / `risk_ratio` / `rmst_difference` contrasts well-formed. |
| ICE sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-ice-survival-delicatessen.R` | `S^1`, `S^0`, `RD` point + sandwich SE match a Python `delicatessen` stacked-EE survival-tail sandwich to ~1e-5 at t ∈ {1,2,3} on a shared fixture. Reference: `data-raw/delicatessen_ice_survival.py`; both read `fixtures/python/ice_survival_data.csv`. Pins the `(1−D_k)` failure carry-forward. |
| ICE sandwich ≈ empirical bootstrap | 🟢 | `test-ice-survival.R` | Per-time sandwich vs bootstrap SE within 20% (n = 1200, 400 reps). Guards against regressing to causatr's verbatim chain (which over-covers, growing in t). |
| Bootstrap variance (per-replicate ICE refit) | 🟢 | `test-ice-survival.R` | `bootstrap_survival()` refits Track B per replicate (lags rebuilt, `confounders_tv` / `history` threaded). |
| External point-estimate oracle: `lmtp::lmtp_tmle(survival)` | 🟢 | `test-ice-survival-oracle.R` | Static strategies: lmtp and ICE both within 0.03 of the forward-sim truth (`skip_if_not_installed`). |
| External oracle: `gfoRmula::gformula_survival()` | 🟢 | `test-ice-survival-oracle.R` | Cross-check when installed; `skip_if_not_installed` + defensive `tryCatch` (API-sensitive). nsimul = 100 000; `expect_equal` tolerance 0.04 (known ~0.02–0.04 gfoRmula under-estimation vs ICE and analytic truth). Directional pin: ICE ≥ gfoRmula − 0.01. |
| `estimator = "ice"` with constant-within-id treatment | 🟢 | `test-ice-survival.R` | Informs `survatr_ice_static_treatment` (Track A cheaper) but proceeds. |
| Time-varying treatment under Track A (gcomp/ipw) | 🟢 | `test-ice-survival.R` | Warns `survatr_tv_treatment_track_a`, points to `estimator = "ice"`. |
| Non-numeric (factor / categorical k>2) treatment under Track B | 🔴 | `test-ice-survival.R` | `survatr_ice_treatment_unsupported`. Numeric (binary / linear dose) only; treatment-design-formula path → later chunk. |
| (MCAR) entry (period-1) censoring under Track B | 🟢 | `test-ice-survival.R` | Standardises over the at-risk-at-baseline population; curve + sandwich finite, matches forward-sim truth, `n` = effective count. |
| External `weights` + `estimator = "ice"` | 🔴 | `test-ice-survival.R` | `survatr_ice_external_weights` (weighted / IPCW longitudinal → later chunk). |
| `ipsi()` / stochastic interventions under Track B | 🔴 | `test-ice-survival.R` | `survatr_ice_intervention_deferred` (weight-path / Monte-Carlo → later chunks). |
| causatr ICE primitive contract pins | 🟢 | `test-causatr-integration.R` | Signatures of `ice_fit_step`, `ice_predict_step`, `ice_build_formula`, `ice_apply_intervention_long`, `create_lag_vars`, `ice_if_setup`, `correct_model`, `new_causatr_fit`, etc. |

## Competing risks (cause-specific hazards + CIF)

`surv_fit(competing = <cause-col>)` fits `J` parallel cause-specific
pooled-logistic hazards on a **shared all-cause risk set** (other causes
censored at their event time), then `contrast()` builds per-cause cumulative
incidence `F^(j),a(t) = Σ_{k≤t} S^a(k−1) ĥ^(j),a(k)` and all-cause survival
`S^a(t) = ∏(1 − Σ_j ĥ^(j))`. gcomp / Track A only this chunk. Cause-specific
hazards only — Fine–Gray / subdistribution is out of scope (documented).

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| CR fit: `J` cause-specific hazards on a shared all-cause risk set | 🟢 | `test-competing-risks.R` | Structure (`cause_models`, `causes`, `model = NULL`); both cause models fit the same `nobs`. |
| `type = "cif"` per-cause CIF, `static` interventions | 🟢 | `test-competing-risks.R` | Closed-form two-cause constant-hazard CIF `F^(j)(t) = (h_j/H)(1 − (1−H)^t)` (tol 0.02 at n = 6000); Aalen–Johansen (`survival::survfit`) on a saturated `~ factor(t)` fit (tol 0.02). |
| All-cause `type = "survival"` / `"risk"` (summed hazards) | 🟢 | `test-competing-risks.R` | Closed-form `(1 − h₁ − h₂)^t` (tol 0.02). |
| `Σ_j F^(j)(t) + S(t) = 1` identity | 🟢 | `test-competing-risks.R` | Numerical identity to 1e-12 per (intervention, time). |
| `type = "cif_difference"` | 🟢 | `test-competing-risks.R`, `test-competing-risks-sandwich.R` | No-effect DGP: CIF-diff CI covers 0 at ≥ 88% nominal 95% (150-rep sim, `skip_on_cran`). |
| `type = "cif_ratio"` (log-scale) | 🟢 | `test-competing-risks-sandwich.R` | Strictly positive, `ci_upper > ci_lower`, covers 1 on a no-effect DGP. |
| `cause = ` subsetting | 🟢 | `test-competing-risks.R` | `cause = 1` returns only cause-1 rows in estimates + contrasts. |
| Stacked-EE sandwich (per-cause hazard blocks + CIF/survival delta) | 🟢 | `test-competing-risks-sandwich.R` | IF columns mean ~ 0; sandwich SE matches the empirical bootstrap to ~2%. |
| CR sandwich vs `delicatessen` (independent analytic M-estimation) | 🟢 | `test-competing-risks-sandwich.R` | Per-cause `F^(j),a(t)`, all-cause `S^a(t)`, and `RD^(j)(t)` point + sandwich SE match a Python `delicatessen` stacked-EE oracle to ~1e-4 on a shared fixture. Reference: `data-raw/delicatessen_competing_risks.py`; both read `fixtures/python/cr_survival_data.csv`. |
| CR CIF vs `riskRegression::ate` g-formula (same-class comparator) | 🟢 | `test-competing-risks-riskregression.R` | Per-cause CIF under each arm on a confounded 2-cause DGP (n = 3000) matches `CSC` + `ate()` g-formula to ≤ 0.06 (relative; max observed ~4% for cause 1, treated arm where pooled-logistic vs Cox differ most). Oracle: `rr_ate_cif_oracle()` in `helper-cr-oracle.R`; requires `riskRegression`. |
| CR bootstrap (per-replicate dual cause-model refit) | 🟢 | `test-competing-risks.R`, `test-competing-risks-sandwich.R` | Cause models re-estimated per resample; CIs populated, point ∈ CI; SE ≈ sandwich. |
| CR × `mgcv::gam` baseline hazard | 🟢 | `test-competing-risks-sandwich.R` | lpmatrix-basis per cause; sandwich SE finite / positive on `s(t, k = 4)`. |
| Cause-aware `print` / `tidy` / `plot` / `forrest` | 🟢 | `test-competing-risks-s3.R` | `cause` column threaded; `print` shows the truncation-by-death caveat for difference / ratio; `tidy` keeps `cause` (NA for all-cause); `plot` / `forrest` render per (group, cause); single-event shapes unchanged. |
| Truncation-by-death caveat | 🟢 | `test-competing-risks.R` (via `suppressMessages`) | One-time `rlang::inform()` + `print` note + vignette. |
| CR × `estimator = "ipw"` / `"ice"` | 🔴 | `test-competing-risks.R` | `survatr_competing_estimator` (gcomp / Track A only; IPW / ICE CR → later chunks). |
| `competing != outcome`, or `< 2` distinct causes | 🔴 | `test-competing-risks.R` | `survatr_competing_misuse`. |
| Non-integer / negative / NA / all-zero cause column | 🔴 | `test-competing-risks.R` | `survatr_bad_competing`. |
| CIF estimand on a single-event fit; single-event contrast on a CR fit | 🔴 | `test-competing-risks.R` | `survatr_competing_type`. |
| Unknown `cause` label | 🔴 | `test-competing-risks.R` | `survatr_bad_cause`. |
| External `weights` + `competing =` | 🔴 | (guarded in `surv_fit()`) | `survatr_competing_weights` (weighted / IPCW competing risks → later chunk). |
| Fine–Gray / subdistribution hazards | — | — | Out of scope (cause-specific only); documented in roxygen + vignette. |
| Per-cause years-of-life-lost (`type = "yll"`) | 🟢 | `test-estimands-yll.R` | Shipped in chunk 12: `∫F^(j)`, the `Σⱼ YLL = RMTL` identity, CIF-IF-through-trapezoidal sandwich. Per-cause RMST remains out of scope. |
| Competing risks under Track B | — | — | Out of scope this chunk (composes after chunks 6 + 7). |

## `diagnose()` — survival-aware diagnostics (Chunk 10)

`diagnose.survatr_fit()` returns a `survatr_diag` with five panels, all operating on
the at-risk rows from `build_risk_set()`. For Track B (ICE), the positivity panel
uses the empirical event rate (no fitted model at `surv_fit()` time).

| Surface | Status | Test file | Oracle |
|---|---|---|---|
| `diagnose()` returns `survatr_diag` with five named panels | 🟢 | `test-diagnose.R` | Structure check; `expect_s3_class`; null panels asserted. |
| `$positivity` columns (`time`, `n_at_risk`, `h_min/h_mean/h_max`, `flag_low/high`) | 🟢 | `test-diagnose.R` | Column names; row count = number of time periods. |
| `flag_low` does NOT fire on moderate-hazard DGP | 🟢 | `test-diagnose.R` | `sim_constant_hazard(h = 0.08)`; all flags FALSE. |
| `flag_low` fires when h < 0.001 | 🟢 | `test-diagnose.R` | `sim_constant_hazard(h = 1e-4)`; at least one period flagged. |
| `$balance` SMDs ≈ 0 on randomized DGP | 🟢 | `test-diagnose.R` | `sim_confounded_survival(gamma = 0)`; all `|SMD| < 0.3`. |
| `$balance` SMDs non-trivial on confounded DGP | 🟢 | `test-diagnose.R` | `sim_confounded_survival(gamma = 1.5)`; some `|SMD| > 0.1`. |
| `$weights` non-NULL for IPW; ESS < n | 🟢 | `test-diagnose.R` | IPW fit on confounded DGP; `ess < n_ids`. |
| `$weights` NULL for gcomp | 🟢 | `test-diagnose.R` | Structural assertion. |
| `$weights` NULL for ICE | 🟢 | `test-diagnose.R` | Structural assertion. |
| `$censoring` NULL when no censoring column | 🟢 | `test-diagnose.R` | Structural assertion. |
| `$censoring` populated with correct columns when `censoring=` supplied | 🟢 | `test-diagnose.R` | Column names; proportions in [0, 1]. |
| `$competing` present with CR fit; correct causes; `Σ F^(j) + S = 1` | 🟢 | `test-diagnose.R` | Identity check to < 1e-6; caveat attribute non-empty. |
| `$competing` NULL for single-event fit | 🟢 | `test-diagnose.R` | Structural assertion. |
| `print.survatr_diag` runs without error | 🟢 | `test-diagnose.R` | Output contains `"survatr_diag"`. |
| `print` shows Weights line for IPW | 🟢 | `test-diagnose.R` | `grepl("Weights:", ...)`. |
| `print` shows Competing line for CR fit | 🟢 | `test-diagnose.R` | `grepl("Competing:", ...)`. |
