# survatr depends on three unexported causatr functions via `:::`. The remote
# is an unpinned GitHub `main`, so a refactor in causatr could change their
# signature or return shape and silently break survatr's fit / sandwich paths.
# These tests turn that into a loud failure: they pin the contract (formal
# names + return shape) that survatr relies on. If causatr drifts, these fail
# in CI instead of the breakage surfacing as a cryptic runtime error.

test_that("causatr:::apply_intervention signature is (data, treatment, iv)", {
  expect_true(exists("apply_intervention", where = asNamespace("causatr")))
  expect_identical(
    names(formals(causatr:::apply_intervention)),
    c("data", "treatment", "iv")
  )
})

test_that("causatr:::prepare_model_if keeps its signature and return shape", {
  expect_identical(
    names(formals(causatr:::prepare_model_if)),
    c("model", "fit_idx", "n_total")
  )
  d <- data.frame(
    y = rep(c(0L, 1L), 50L),
    x = as.numeric(seq_len(100L))
  )
  m <- stats::glm(y ~ x, family = stats::binomial(), data = d)
  prep <- causatr:::prepare_model_if(
    m,
    fit_idx = seq_len(100L),
    n_total = 100L
  )
  # survatr reads exactly these three components in variance_if_survival.R.
  expect_true(all(c("X_fit", "B_inv", "r_score") %in% names(prep)))
  expect_equal(ncol(prep$X_fit), length(stats::coef(m)))
  expect_equal(dim(prep$B_inv), rep(length(stats::coef(m)), 2L))
  expect_equal(length(prep$r_score), 100L)
})

test_that("causatr:::fit_treatment_model keeps its signature + return shape", {
  # survatr's IPW path fits the propensity (and marginal) model via this and
  # reads `model` / `family` / `alpha_hat` / `X_prop` / `fit_rows`.
  expect_true(exists("fit_treatment_model", where = asNamespace("causatr")))
  expect_identical(
    names(formals(causatr:::fit_treatment_model)),
    c(
      "data",
      "treatment",
      "confounders",
      "model_fn",
      "propensity_family",
      "weights",
      "..."
    )
  )
  d <- data.frame(A = rbinom(200L, 1L, 0.5), L = rnorm(200L))
  tm <- causatr:::fit_treatment_model(
    data.table::as.data.table(d),
    treatment = "A",
    confounders = ~L
  )
  expect_true(all(
    c("model", "family", "alpha_hat", "X_prop", "fit_rows") %in% names(tm)
  ))
  expect_identical(tm$family, "bernoulli")
  expect_equal(ncol(tm$X_prop), length(tm$alpha_hat))
})

test_that("causatr:::evaluate_density signature is (model, values, newdata)", {
  expect_identical(
    names(formals(causatr:::evaluate_density)),
    c("treatment_model", "treatment_values", "newdata")
  )
})

test_that("causatr:::truncate_weights + apply_model_correction signatures", {
  # `truncate_weights` does the IPW weight winsorization;
  # `apply_model_correction` supplies the treatment-model correction block.
  expect_identical(names(formals(causatr:::truncate_weights)), c("w", "trim"))
  expect_identical(
    names(formals(causatr:::apply_model_correction)),
    c("prep", "gradient")
  )
  # `apply_model_correction` returns `correction` / `d` / `h`; survatr reads
  # `correction` per time column.
  d <- data.frame(y = rep(c(0L, 1L), 50L), x = as.numeric(seq_len(100L)))
  m <- stats::glm(y ~ x, family = stats::binomial(), data = d)
  prep <- causatr:::prepare_model_if(m, fit_idx = seq_len(100L), n_total = 100L)
  res <- causatr:::apply_model_correction(prep, stats::coef(m))
  expect_true(all(c("correction", "d", "h") %in% names(res)))
  expect_equal(length(res$correction), 100L)
})

test_that("causatr:::iv_design_matrix signature + coef-aligned columns", {
  expect_identical(
    names(formals(causatr:::iv_design_matrix)),
    c("model", "newdata")
  )
  d <- data.frame(
    y = rep(c(0L, 1L), 50L),
    x = as.numeric(seq_len(100L))
  )
  m <- stats::glm(y ~ x, family = stats::binomial(), data = d)
  x_des <- causatr:::iv_design_matrix(m, d)
  # The counterfactual design must have one column per (non-aliased) coef so
  # it is conformable with the B_inv bread.
  expect_equal(ncol(x_des), length(stats::coef(m)))
  expect_equal(nrow(x_des), nrow(d))
})
