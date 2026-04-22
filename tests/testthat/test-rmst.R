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

## Regression test for B1 (2026-04-22 critical review, round 1):
## `rmst_weights()` previously had off-by-one indexing on the per-interval
## `dt` contributions AND double-counted `S(0) = 1` onto the first column
## of every row, inflating the sandwich SE for `rmst` and
## `rmst_difference` by up to 2x at t_1 and ~57% on irregular grids.
## Fixed in the same commit as this test. Repro script:
## `/tmp/survatr_repro_b1_rmst_weights.R`. Contract:
##   W %*% s_hat + dt[1]/2 == trapezoidal_rmst(times, s_hat)
## for any `times` (sorted, first > 0) and any `s_hat` in [0, 1].
test_that("rmst_weights matches trapezoidal_rmst on arbitrary grids (B1)", {
  for (times in list(
    1:5,
    seq(2, 10, by = 2),
    c(1, 3, 7, 12),
    c(2.5, 5, 7.5, 10),
    c(0.1, 0.5, 1.0, 2.0, 5.0)
  )) {
    W <- rmst_weights(times)
    dt1_half <- times[1L] / 2 ## dt[1] = t_1 - 0 = t_1; constant for S(0) = 1
    s_hat <- seq(1 - 1e-3, 1 - 0.9, length.out = length(times))
    rmst_from_W <- as.numeric(W %*% s_hat) + dt1_half
    rmst_direct <- trapezoidal_rmst(times, s_hat)
    expect_equal(rmst_from_W, rmst_direct, tolerance = 1e-12)
  }
})

test_that("rmst_weights has correct per-row trapezoid structure (B1)", {
  times <- c(1, 3, 4, 7)
  dt <- diff(c(0, times))
  W <- rmst_weights(times)
  ## Row j = 1: only the (0, t_1) half-interval contributes to S(t_1).
  expect_equal(W[1L, ], c(dt[1L] / 2, 0, 0, 0))
  ## Row j > 1: column 1 gets dt[1]/2 + dt[2]/2; column j gets dt[j]/2;
  ## interior columns get (dt[i] + dt[i+1]) / 2.
  expect_equal(W[2L, 1L], dt[1L] / 2 + dt[2L] / 2)
  expect_equal(W[2L, 2L], dt[2L] / 2)
  expect_equal(W[3L, 2L], (dt[2L] + dt[3L]) / 2)
  expect_equal(W[3L, 3L], dt[3L] / 2)
  expect_equal(W[4L, 3L], (dt[3L] + dt[4L]) / 2)
  expect_equal(W[4L, 4L], dt[4L] / 2)
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
