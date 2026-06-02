test_that("sandwich log-RR CI is strictly positive and covers 1 under null", {
  dt <- sim_constant_hazard(n = 3000L, K = 8L, h = 0.08, seed = 181L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(1, 4, 8),
    type = "risk_ratio",
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_true(all(res$contrasts$ci_lower > 0))
  expect_true(all(res$contrasts$ci_upper > res$contrasts$ci_lower))
  ## On a no-effect DGP the CI should include 1 at every time.
  expect_true(all(res$contrasts$ci_lower <= 1 & 1 <= res$contrasts$ci_upper))
})

test_that("sandwich risk-ratio se is on the natural scale (matches bootstrap)", {
  skip_on_cran()
  ## The risk-ratio `se` column reports the natural-scale SE of RR
  ## (delta method, `se(RR) = RR * se(log RR)`), the same quantity the
  ## bootstrap reports (`sd(RR)`). The DGP is chosen so RR is well away from 1
  ## (treatment protective), so a log-scale `se` would differ from the natural
  ## one by the factor RR and the two ci_methods would disagree.
  dt <- sim_confounded_survival(n = 3000L, K = 6L, seed = 31L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "gcomp"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  sw <- contrast(
    fit,
    interventions = ivs,
    times = c(3, 6),
    type = "risk_ratio",
    reference = "a0",
    ci_method = "sandwich"
  )
  bs <- contrast(
    fit,
    interventions = ivs,
    times = c(3, 6),
    type = "risk_ratio",
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 400L,
    seed = 4L
  )
  ## RR is clearly below 1, so natural and log scales differ materially.
  expect_true(all(sw$contrasts$estimate < 0.85))
  ## Natural-scale sandwich se agrees with the bootstrap se (both SE of RR).
  expect_equal(sw$contrasts$se, bs$contrasts$se, tolerance = 0.2)
})
