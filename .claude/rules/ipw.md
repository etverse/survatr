---
paths:
  - "R/ipw_survival.R"
  - "R/ipw_weights_survival.R"
  - "R/variance_if_ipw_survival.R"
  - "tests/testthat/test-ipw-*.R"
  - "tests/testthat/helper-ipw-survival-oracle.R"
---

# IPW

These are deliberate. Do not report them as bugs without a numerical reproducer.

- Stabilized weights are composed in survatr from two `causatr:::evaluate_density()` calls (`A ~ L` and `A ~ 1`) and `causatr:::truncate_weights()`, because causatr's point IPW is a Hájek `Y ~ 1` model that exposes no `f(A)/f(A|L)` weight. This is composition of causatr's primitives.
- The treatment-model correction in the sandwich carries exactly one factor of `n_ids`: it comes from `numDeriv` on the weighted-MSM `phi_bar` (divided by `n_fit`) and `causatr:::apply_model_correction(prep_trt, g)` with `n_total = n_ids`, and the `n_fit` in `h_t = n_fit · B_inv · J̄` cancels `phi_bar`'s `1/n_fit`. It matches the two-stage bootstrap and an independent `delicatessen` stacked sandwich to 1e-4 (`test-ipw-delicatessen.R`); rerun those before changing the factor.
- With stabilized weights the sandwich SE is narrower than the SE that treats the weights as known, because the treatment-model correction is subtracted (Robins 1999; Hernán, Brumback and Robins 2000). The estimation of the marginal numerator `A ~ 1` is left out of the influence function, which is conservative.
- `ipsi()` is rejected for IPW survival.
