---
paths:
  - "R/competing_risks.R"
  - "R/variance_if_competing.R"
  - "R/estimand_registry.R"
  - "tests/testthat/test-competing-risks*.R"
  - "tests/testthat/test-estimands-yll.R"
  - "tests/testthat/helper-cr-oracle.R"
---

# Competing risks

These are deliberate. Do not report them as bugs without a numerical reproducer.

- Competing risks are cause-specific hazards with a cumulative incidence function. Fine-Gray subdistribution hazards are out of scope.
- One all-cause risk set, J cause-specific models on it: `surv_fit(competing = <col>)` derives `1{competing ≠ 0}`, calls `build_risk_set()` once, and fits cause j with response `1{competing == j}` on those rows (`fit_competing_risks()`). Sharing the risk set is exactly treating other causes as censored at their event time. Do not build a risk set per cause or pass the multi-valued cause column to `build_risk_set()`, which expects a 0/1 indicator. The cause column is checked by `check_competing_col()`, and the 0/1 check on `outcome` is skipped on this path.
- The CIF sandwich sensitivity divides by the all-cause `1 − H`, not `1 − h⁽ʲ⁾`: `s⁽ʲ'⁾ₗ = h⁽ʲ'⁾(1 − h⁽ʲ'⁾)/(1 − Hₗ)` with `H = Σⱼ h⁽ʲ⁾`. It reduces to the single-event form only when J = 1.
- The CIF derivative uses the lagged cumulative sensitivity, because `F⁽ʲ⁾ᵢ(t) = Σ_{k≤t} S(k−1) h⁽ʲ⁾ₖ`: it is coded `cumSX_lag = cumSX − SX`. The all-cause survival influence function uses the inclusive `cumSX(t)`. Mixing the two is a bug.
- The bread is block-diagonal across causes, since each cause model's score involves only its own coefficients; the correlation between causes enters through the per-individual influence function summed over causes. Do not build a joint bread.
- Validated to 1e-4 against `delicatessen` (`data-raw/delicatessen_competing_risks.py`, `fixtures/python/cr_survival_data.csv`) and to 2% against the bootstrap, in `test-competing-risks-sandwich.R`.
- Scope: g-computation only. Rejections: `ipw` or `ice` with `competing` (`survatr_competing_estimator`), `outcome ≠ competing` or fewer than two causes (`survatr_competing_misuse`), a bad cause column (`survatr_bad_competing`), a CIF estimand on a single-event fit or the reverse (`survatr_competing_type`), an unknown `cause` (`survatr_bad_cause`), external `weights` with `competing` (`survatr_competing_weights`).
- CIF contrasts (`cif_difference`, `cif_ratio`) carry a truncation-by-death caveat: a once-per-session `rlang::inform()`, a `print()` note and the vignette. Keep it.
- `cr_rows_by_time()` checks one row per id per requested time (else `survatr_if_failed`). It cannot fail on rectangular data today, but it is the boundary that ragged or left-truncated data will need; it is not dead code.
- Years of life lost: `YLL⁽ʲ⁾ = ∫ F⁽ʲ⁾`, through the CIF influence function and `rmst_weights()` with no `S(0) = 1` constant, so `Σⱼ YLL⁽ʲ⁾(t*) = RMTL(t*)` to machine precision. YLL needs a competing-risks fit (`survatr_yll_needs_cr`). Per-cause RMST and competing-risks `*_difference` contrasts are out of scope.
