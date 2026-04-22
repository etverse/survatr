## Closed-form oracle for the trapezoidal RMST on a constant-hazard curve.
## With S(t) = (1 - h)^t and S(0) = 1, the trapezoidal integral from 0 to
## integer t on unit-spaced grid is
##   RMST(t) = sum_{j=0}^{t-1} ((1-h)^j + (1-h)^(j+1)) / 2
##          = (1 + (1-h)) / 2 * sum_{j=0}^{t-1} (1-h)^j (up to boundary)
## The function itself is small enough to assert against a manual expansion
## without re-using the closed form -- we do both.

test_that("trapezoidal_rmst matches explicit trapezoid sum on (1-h)^t", {
  h <- 0.1
  times <- 1:6
  s <- (1 - h)^times
  rmst <- trapezoidal_rmst(times, s)

  ## Manual trapezoid sum with S(0) = 1 prepended.
  t_full <- c(0, times)
  s_full <- c(1, s)
  dt <- diff(t_full)
  avg <- (s_full[-length(s_full)] + s_full[-1L]) / 2
  manual <- cumsum(dt * avg)

  expect_equal(rmst, manual, tolerance = 1e-12)
  expect_equal(length(rmst), length(times))
})

test_that("trapezoidal_rmst is monotone non-decreasing for S in [0, 1]", {
  s <- c(0.99, 0.95, 0.80, 0.60, 0.30, 0.05)
  rmst <- trapezoidal_rmst(seq_along(s), s)
  expect_true(all(diff(rmst) >= 0))
})

test_that("rmst_weights gives a cumulative trapezoid quadrature matrix", {
  times <- 1:4
  W <- rmst_weights(times)
  expect_equal(dim(W), c(length(times), length(times)))
  ## W is lower triangular in the sense that row j only uses columns <= j.
  K <- length(times)
  for (j in seq_len(K - 1L)) {
    expect_true(all(W[j, (j + 1L):K] == 0))
  }
})

test_that("add_rmst_to_estimates replaces s_hat / risk_hat with rmst_hat", {
  est <- data.table::data.table(
    intervention = rep(c("a1", "a0"), each = 4L),
    time = rep(1:4, 2L),
    s_hat = c(0.9, 0.8, 0.7, 0.6, 0.95, 0.9, 0.85, 0.8),
    risk_hat = 1 - c(0.9, 0.8, 0.7, 0.6, 0.95, 0.9, 0.85, 0.8),
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    n = 100L
  )
  out <- add_rmst_to_estimates(est, times = 1:4)
  expect_true("rmst_hat" %in% names(out))
  expect_false("s_hat" %in% names(out))
  expect_false("risk_hat" %in% names(out))
  expect_equal(
    out[intervention == "a1", rmst_hat],
    trapezoidal_rmst(1:4, c(0.9, 0.8, 0.7, 0.6))
  )
})
