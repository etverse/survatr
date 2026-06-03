# Internal-invariant guards. These cannot be reached through the public API
# (the upstream code always supplies aligned inputs), so they are exercised by
# calling the internal helpers directly with deliberately malformed inputs.

test_that("compute_survival_curve() rejects a wrong-length hazard vector", {
  pp <- data.table::data.table(
    id = rep(1:2, each = 3L),
    t = rep(1:3, times = 2L)
  )
  expect_error(
    survatr:::compute_survival_curve(
      pp_data = pp,
      hazards = numeric(3L), # pp has 6 rows -> misaligned
      id = "id",
      time = "t",
      times = c(1, 2, 3),
      intervention_name = "a1"
    ),
    class = "survatr_hazard_misaligned"
  )
})

test_that("fill_sandwich_ses() rejects an IF matrix of the wrong dimensions", {
  # IF_mat is 3 x 2, but n_ids = 5 and length(times) = 4 are expected.
  if_list <- list(a1 = list(IF_mat = matrix(0, nrow = 3L, ncol = 2L)))
  expect_error(
    survatr:::fill_sandwich_ses(
      estimates = data.table::data.table(),
      contrasts = data.table::data.table(),
      if_list = if_list,
      type = "survival",
      reference = NULL,
      times = c(1, 2, 3, 4),
      conf_level = 0.95,
      n_ids = 5L
    ),
    class = "survatr_if_failed"
  )
})
