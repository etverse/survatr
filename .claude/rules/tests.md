---
paths:
  - "tests/testthat/**"
---

# Tests

- Run a file with `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-<name>.R")'`. Slow tests are behind `skip_on_cran()`; set `NOT_CRAN=true` to run them.
- Shared data-generating processes are in `helper-dgp.R`; oracle wrappers in the `helper-*-oracle.R` files; Python oracle scripts in `data-raw/` and their outputs in `fixtures/python/`.
- The NHEFS chapter 17 replication is the acceptance target for point g-computation: 120-month survival of 80.7% treated and 80.5% untreated, risk difference 0.2% (95% CI −4.1% to 3.7%).
- Before reporting a variance bug, compare sandwich with bootstrap on a simulation with known truth; a derivation on paper alone is not a reproducer.
- `lmtp::lmtp_tmle(outcome_type = "survival")` is an oracle for point estimates only; its efficient-influence-function SE is not comparable to the M-estimation sandwich.
- Pooled logistic and Cox agree closely only for per-interval hazards below about 0.1 (Hernán and Robins, technical point 17.1); check the grid spacing before reading a larger disagreement as a bug.
- `expect_equal(tolerance = t)` is relative: the absolute threshold is `t × mean(|target|)`, so `tolerance = 0.04` on survival values of 0.07 to 0.30 allows only about 0.007. For an absolute comparison between estimators use `all.equal(..., scale = 1, tolerance = )` and report its result through `testthat::fail()`.
- `rlang::inform(.frequency = "once")` fires once per R session, and the state persists across test files. Only the test that expects the message should trigger it: `survatr_ice_static_treatment` is expected in `test-ice-survival.R`, so no other file calls `surv_fit(estimator = "ice")` with a treatment constant within id; use a time-varying treatment such as `sim_ice_feedback` instead.
