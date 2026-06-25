## Cluster-robust sandwich variance (chunk 13).
##
## The per-individual influence function is the variance currency in survatr;
## clustering is one aggregation step (sum the IF rows within cluster before the
## cross-product, keeping the n^2 divisor). These tests pin:
##   1. the aggregation primitive against `sandwich::vcovCL` (gold standard);
##   2. the degenerate-to-current invariant (`cluster = id` == per-individual);
##   3. calibration on a multi-site frailty DGP (naive under-covers, clustered
##      restores the cluster-sampling SD) + clustered SE >= per-individual SE;
##   4. cluster-resampling bootstrap agreement with the clustered sandwich;
##   5. the validation aborts and the Track B deferral.

## ---- 1. aggregation primitive vs sandwich::vcovCL -------------------------

test_that("clustered_pointwise_se matches sandwich::vcovCL on the mean", {
  skip_if_not_installed("sandwich")
  ## For the sample mean, the influence function is IF_i = y_i - ybar and
  ## theta_hat - theta ~ (1/n) sum IF_i, so the cluster-robust variance of the
  ## mean is exactly `crossprod(rowsum(IF, cluster)) / n^2`. `vcovCL(type =
  ## "HC0", cadjust = FALSE)` on `lm(y ~ 1)` is the canonical reference.
  set.seed(3)
  G <- 25L
  m <- 6L
  n <- G * m
  clus <- rep(seq_len(G), each = m)
  y <- 2 + rnorm(G, sd = 1.5)[clus] + rnorm(n, sd = 0.4)
  if_mat <- matrix(y - mean(y), ncol = 1L)
  fit <- stats::lm(y ~ 1)

  se_cl <- clustered_pointwise_se(if_mat, n, as.character(clus))
  v_cl <- sandwich::vcovCL(fit, cluster = clus, type = "HC0", cadjust = FALSE)
  expect_equal(se_cl, sqrt(as.numeric(v_cl)), tolerance = 1e-10)

  ## Singleton clusters (cluster = id) reduce to the ordinary HC0 sandwich.
  se_id <- clustered_pointwise_se(if_mat, n, as.character(seq_len(n)))
  se_hc0 <- sqrt(as.numeric(sandwich::vcovHC(fit, type = "HC0")))
  expect_equal(se_id, se_hc0, tolerance = 1e-10)
  ## ... and equal the no-cluster reduction (the shipped per-individual path).
  se_none <- clustered_pointwise_se(if_mat, n, NULL)
  expect_equal(se_id, se_none, tolerance = 1e-12)
})

## ---- 2. degenerate invariant: cluster = id == per-individual --------------

test_that("cluster = id reproduces the per-individual sandwich SE exactly", {
  pp <- sim_clustered_survival(n_per = 20L, G = 25L, K = 6L, seed = 31L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  for (ty in c("survival", "risk_difference", "rmst_difference")) {
    none <- contrast(fit, ivs, times = 1:6, type = ty, ci_method = "sandwich")
    idcl <- contrast(
      fit,
      ivs,
      times = 1:6,
      type = ty,
      ci_method = "sandwich",
      cluster = "id"
    )
    expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-12)
    expect_equal(none$contrasts$se, idcl$contrasts$se, tolerance = 1e-12)
  }

  ## The q-indexed quantile takes its own assembly; pin it too.
  none_q <- contrast(
    fit,
    ivs,
    times = 1:6,
    type = "quantile",
    q = 0.3,
    ci_method = "sandwich"
  )
  id_q <- contrast(
    fit,
    ivs,
    times = 1:6,
    type = "quantile",
    q = 0.3,
    ci_method = "sandwich",
    cluster = "id"
  )
  expect_equal(none_q$estimates$se, id_q$estimates$se, tolerance = 1e-12)
  expect_equal(none_q$contrasts$se, id_q$contrasts$se, tolerance = 1e-12)
})

test_that("competing-risks cluster = id reproduces the per-individual SE", {
  ## Multi-site competing-risks DGP: shared site frailty on both causes.
  set.seed(5)
  G <- 25L
  np <- 20L
  K <- 6L
  u <- rnorm(G, sd = 1.0)
  rows <- vector("list", G * np)
  idx <- 0L
  for (g in seq_len(G)) {
    for (j in seq_len(np)) {
      idx <- idx + 1L
      a_i <- rbinom(1L, 1L, 0.5)
      cause <- integer(K)
      for (k in seq_len(K)) {
        p1 <- plogis(-2.4 - 0.4 * a_i + u[g])
        p2 <- plogis(-2.8 + u[g])
        d <- sample(0:2, 1L, prob = c(1 - p1 - p2, p1, p2))
        cause[k] <- d
        if (d > 0L) {
          if (k < K) {
            cause[(k + 1L):K] <- 0L
          }
          break
        }
      }
      rows[[idx]] <- data.table::data.table(
        id = idx,
        site = g,
        t = seq_len(K),
        A = a_i,
        cause = cause
      )
    }
  }
  ppc <- data.table::rbindlist(rows)
  fc <- surv_fit(
    ppc,
    "cause",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "cause"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  none <- contrast(fc, ivs, times = 1:6, type = "cif", ci_method = "sandwich")
  idcl <- contrast(
    fc,
    ivs,
    times = 1:6,
    type = "cif",
    ci_method = "sandwich",
    cluster = "id"
  )
  expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-12)
})

## ---- 3. widening + calibration on a correlated multi-site DGP --------------

test_that("clustered SE >= per-individual SE under within-cluster corr", {
  pp <- sim_clustered_survival(n_per = 25L, G = 40L, sigma = 1.4, seed = 11L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  ## Level (risk): the shared frailty correlates within-site outcomes, so the
  ## clustered SE is uniformly wider.
  none <- contrast(fit, ivs, times = 1:8, type = "risk", ci_method = "sandwich")
  site <- contrast(
    fit,
    ivs,
    times = 1:8,
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )
  se_none <- none$estimates[get("intervention") == "a1"]$se
  se_site <- site$estimates[get("intervention") == "a1"]$se
  expect_true(all(se_site >= se_none - 1e-9))
  ## Meaningful widening (verified ~2.3x on this DGP); band guards both sides.
  ratio <- mean(se_site / se_none)
  expect_true(ratio > 1.3 && ratio < 6)
})

test_that("contrast SE widens under cluster-level treatment", {
  ## With treatment assigned at the site level the frailty does NOT cancel in
  ## the risk difference, so the contrast SE itself widens.
  pp <- sim_clustered_survival(
    n_per = 25L,
    G = 40L,
    sigma = 1.4,
    cluster_trt = TRUE,
    seed = 21L
  )
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  none <- contrast(
    fit,
    ivs,
    times = 1:8,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  site <- contrast(
    fit,
    ivs,
    times = 1:8,
    type = "risk_difference",
    ci_method = "sandwich",
    cluster = "site"
  )
  expect_true(all(site$contrasts$se >= none$contrasts$se - 1e-9))
  ratio <- mean(site$contrasts$se / none$contrasts$se)
  expect_true(ratio > 1.5 && ratio < 8)
})

test_that("clustered SE is calibrated to the cluster-sampling distribution", {
  skip_on_cran()
  ## Truth: the SD of the risk@t* estimator across independent re-draws of the
  ## G sites (fresh frailties) IS what the cluster-robust SE estimates. The
  ## per-individual SE ignores the between-site variance and under-states it.
  G <- 40L
  np <- 20L
  K <- 6L
  t_star <- 6L
  ivs <- list(a1 = causatr::static(1))
  R <- 150L
  ests <- numeric(R)
  for (r in seq_len(R)) {
    ppr <- sim_clustered_survival(
      n_per = np,
      G = G,
      K = K,
      sigma = 1.2,
      seed = 1000L + r
    )
    fr <- surv_fit(ppr, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
    rr <- contrast(fr, ivs, times = t_star, type = "risk", ci_method = "none")
    ests[r] <- rr$estimates[get("intervention") == "a1"]$risk_hat
  }
  emp_sd <- stats::sd(ests)

  pp0 <- sim_clustered_survival(
    n_per = np,
    G = G,
    K = K,
    sigma = 1.2,
    seed = 1L
  )
  fit0 <- surv_fit(pp0, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  se_site <- contrast(
    fit0,
    ivs,
    times = t_star,
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )$estimates[get("intervention") == "a1"]$se
  se_naive <- contrast(
    fit0,
    ivs,
    times = t_star,
    type = "risk",
    ci_method = "sandwich"
  )$estimates[get("intervention") == "a1"]$se

  ## Clustered SE tracks the empirical cluster-sampling SD; naive understates it.
  expect_equal(se_site, emp_sd, tolerance = 0.30)
  expect_true(se_naive < emp_sd * 0.85)
})

## ---- 4. cluster-resampling bootstrap ----------------------------------------

test_that("bootstrap resamples clusters and matches the clustered sandwich", {
  skip_on_cran()
  pp <- sim_clustered_survival(n_per = 20L, G = 35L, sigma = 1.3, seed = 41L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  se_sand <- contrast(
    fit,
    ivs,
    times = c(4L, 8L),
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )$estimates[get("intervention") == "a1"]$se
  se_boot <- contrast(
    fit,
    ivs,
    times = c(4L, 8L),
    type = "risk",
    ci_method = "bootstrap",
    cluster = "site",
    n_boot = 400L,
    seed = 7L
  )$estimates[get("intervention") == "a1"]$se
  ## Cluster bootstrap (resample sites) agrees with the clustered sandwich.
  expect_equal(se_boot, se_sand, tolerance = 0.25)

  ## ... and is wider than the naive (resample-individuals) bootstrap.
  se_boot_naive <- contrast(
    fit,
    ivs,
    times = c(4L, 8L),
    type = "risk",
    ci_method = "bootstrap",
    n_boot = 400L,
    seed = 7L
  )$estimates[get("intervention") == "a1"]$se
  expect_true(all(se_boot >= se_boot_naive * 0.9))
  expect_true(mean(se_boot / se_boot_naive) > 1.2)
})

## ---- 5. deferrals and validation aborts -------------------------------------

test_that("cluster is rejected for Track B (ICE) fits", {
  pp <- sim_clustered_survival(n_per = 15L, G = 20L, K = 5L, seed = 51L)
  fit_b <- surv_fit(
    pp,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ice"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  expect_error(
    contrast(
      fit_b,
      ivs,
      times = 1:5,
      type = "risk_difference",
      ci_method = "sandwich",
      cluster = "site"
    ),
    class = "survatr_cluster_track_b_deferred"
  )
})

test_that("cluster validation aborts fire with classed errors", {
  pp <- sim_clustered_survival(n_per = 15L, G = 20L, K = 5L, seed = 61L)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  ## Varies within id.
  pp_v <- data.table::copy(pp)
  pp_v[get("id") == 1L & get("t") == 1L, "site" := 999L]
  fit_v <- surv_fit(pp_v, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  expect_snapshot(
    contrast(
      fit_v,
      ivs,
      1:5,
      "risk_difference",
      ci_method = "sandwich",
      cluster = "site"
    ),
    error = TRUE
  )

  ## NA cluster label.
  pp_na <- data.table::copy(pp)
  pp_na[1L, "site" := NA_integer_]
  fit_na <- surv_fit(pp_na, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  expect_snapshot(
    contrast(
      fit_na,
      ivs,
      1:5,
      "risk_difference",
      ci_method = "sandwich",
      cluster = "site"
    ),
    error = TRUE
  )

  ## Single cluster (G = 1).
  pp_1 <- data.table::copy(pp)
  pp_1[, "site" := 1L]
  fit_1 <- surv_fit(pp_1, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  expect_snapshot(
    contrast(
      fit_1,
      ivs,
      1:5,
      "risk_difference",
      ci_method = "sandwich",
      cluster = "site"
    ),
    error = TRUE
  )

  ## Unknown / malformed column.
  fit_ok <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  expect_snapshot(
    contrast(
      fit_ok,
      ivs,
      1:5,
      "risk_difference",
      ci_method = "sandwich",
      cluster = "nope"
    ),
    error = TRUE
  )
  expect_error(
    contrast(
      fit_ok,
      ivs,
      1:5,
      "risk_difference",
      ci_method = "sandwich",
      cluster = c("a", "b")
    ),
    class = "survatr_bad_cluster"
  )
})

## ---- 6. every single-event estimand reproduces under cluster = id -----------

test_that("cluster = id reproduces per-individual SE for all single-event types", {
  pp <- sim_clustered_survival(n_per = 20L, G = 25L, K = 6L, seed = 33L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  types <- c(
    "survival",
    "risk",
    "risk_difference",
    "risk_ratio",
    "rmst",
    "rmst_difference",
    "rmtl",
    "rmtl_difference"
  )
  for (ty in types) {
    none <- contrast(fit, ivs, times = 1:6, type = ty, ci_method = "sandwich")
    idcl <- contrast(
      fit,
      ivs,
      times = 1:6,
      type = ty,
      ci_method = "sandwich",
      cluster = "id"
    )
    ## waldo treats NA == NA as equal, so log-scale ratio SEs with early-time
    ## NA columns compare cleanly.
    expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-10)
    if (nrow(none$contrasts) > 0L) {
      expect_equal(none$contrasts$se, idcl$contrasts$se, tolerance = 1e-10)
    }
  }
})

## ---- 7. IPW path ------------------------------------------------------------

test_that("IPW: cluster = id reproduces the per-individual stacked-EE SE", {
  pp <- sim_confounded_survival(n = 600L, K = 5L, seed = 7L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  for (ty in c("survival", "risk_difference", "risk_ratio")) {
    none <- contrast(fit, ivs, times = 1:5, type = ty, ci_method = "sandwich")
    idcl <- contrast(
      fit,
      ivs,
      times = 1:5,
      type = ty,
      ci_method = "sandwich",
      cluster = "id"
    )
    expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-9)
    if (nrow(none$contrasts) > 0L) {
      expect_equal(none$contrasts$se, idcl$contrasts$se, tolerance = 1e-9)
    }
  }
})

test_that("IPW: clustered SE >= per-individual SE on a confounded multi-site DGP", {
  pp <- sim_clustered_confounded(n_per = 25L, G = 40L, K = 6L, seed = 71L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  none <- contrast(fit, ivs, times = 1:6, type = "risk", ci_method = "sandwich")
  site <- contrast(
    fit,
    ivs,
    times = 1:6,
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )
  sn <- none$estimates[get("intervention") == "a1"]$se
  ss <- site$estimates[get("intervention") == "a1"]$se
  expect_true(all(ss >= sn - 1e-9))
  expect_true(mean(ss / sn) > 1.1)
})

## ---- 8. IPCW path -----------------------------------------------------------

test_that("IPCW: cluster = id reproduces the three-block stacked-EE SE", {
  pp <- sim_informative_censoring(n = 600L, K = 5L, seed = 9L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    censoring = "C",
    estimator = "ipw",
    ipcw = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  none <- contrast(
    fit,
    ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  idcl <- contrast(
    fit,
    ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "sandwich",
    cluster = "id"
  )
  expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-9)
  expect_equal(none$contrasts$se, idcl$contrasts$se, tolerance = 1e-9)
})

## ---- 9. competing-risks: every CR estimand reproduces under cluster = id -----

test_that("competing risks: cluster = id reproduces per-individual SE for all CR types", {
  pp <- sim_two_cause_constant_hazard(
    n = 800L,
    K = 6L,
    beta1_A = -0.3,
    beta1_L = 0.4,
    gamma = 0.3,
    seed = 4L
  )
  fit <- surv_fit(
    pp,
    "event",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "event"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  ## Per-cause CIF estimands.
  for (ty in c("cif", "cif_difference", "cif_ratio", "yll")) {
    none <- contrast(fit, ivs, times = 1:6, type = ty, ci_method = "sandwich")
    idcl <- contrast(
      fit,
      ivs,
      times = 1:6,
      type = ty,
      ci_method = "sandwich",
      cluster = "id"
    )
    expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-9)
    if (nrow(none$contrasts) > 0L) {
      expect_equal(none$contrasts$se, idcl$contrasts$se, tolerance = 1e-9)
    }
  }
  ## All-cause estimands on the CR fit.
  for (ty in c("survival", "risk", "rmst", "rmtl")) {
    none <- contrast(fit, ivs, times = 1:6, type = ty, ci_method = "sandwich")
    idcl <- contrast(
      fit,
      ivs,
      times = 1:6,
      type = ty,
      ci_method = "sandwich",
      cluster = "id"
    )
    expect_equal(none$estimates$se, idcl$estimates$se, tolerance = 1e-9)
  }
})

test_that("competing risks: clustered SE >= per-individual SE (CIF, frailty DGP)", {
  ## Multi-site competing-risks DGP with a shared site frailty on both causes.
  set.seed(13)
  G <- 35L
  np <- 22L
  K <- 6L
  u <- rnorm(G, sd = 1.1)
  rows <- vector("list", G * np)
  idx <- 0L
  for (g in seq_len(G)) {
    for (j in seq_len(np)) {
      idx <- idx + 1L
      a_i <- rbinom(1L, 1L, 0.5)
      cause <- integer(K)
      for (k in seq_len(K)) {
        p1 <- plogis(-2.4 - 0.4 * a_i + u[g])
        p2 <- plogis(-2.8 + u[g])
        d <- sample(0:2, 1L, prob = c(1 - p1 - p2, p1, p2))
        cause[k] <- d
        if (d > 0L) {
          if (k < K) {
            cause[(k + 1L):K] <- 0L
          }
          break
        }
      }
      rows[[idx]] <- data.table::data.table(
        id = idx,
        site = g,
        t = seq_len(K),
        A = a_i,
        cause = cause
      )
    }
  }
  ppc <- data.table::rbindlist(rows)
  fc <- surv_fit(
    ppc,
    "cause",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "cause"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  none <- contrast(fc, ivs, times = 1:6, type = "cif", ci_method = "sandwich")
  site <- contrast(
    fc,
    ivs,
    times = 1:6,
    type = "cif",
    ci_method = "sandwich",
    cluster = "site"
  )
  expect_true(all(site$estimates$se >= none$estimates$se - 1e-9))
  expect_true(mean(site$estimates$se / none$estimates$se) > 1.1)
})

## ---- 10. cluster-label type / ordering invariance ---------------------------

test_that("cluster SE is invariant to label type (integer / character / factor)", {
  pp <- sim_clustered_survival(n_per = 20L, G = 15L, K = 5L, seed = 81L)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  run <- function(dt) {
    fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
    contrast(
      fit,
      ivs,
      times = 1:5,
      type = "risk_difference",
      ci_method = "sandwich",
      cluster = "site"
    )$contrasts$se
  }
  se_int <- run(pp)

  ## Character labels whose sort order differs from the integer order ("S10" <
  ## "S2") exercises `cluster_for_ids()` aligning by NAME, not position.
  pp_chr <- data.table::copy(pp)
  pp_chr[, "site" := paste0("S", get("site"))]
  se_chr <- run(pp_chr)

  pp_fac <- data.table::copy(pp)
  pp_fac[, "site" := factor(paste0("S", get("site")))]
  se_fac <- run(pp_fac)

  expect_equal(se_int, se_chr, tolerance = 1e-12)
  expect_equal(se_int, se_fac, tolerance = 1e-12)
})

## ---- 11. edge cases & reproducibility ---------------------------------------

test_that("cluster is validated even when ci_method = 'none'", {
  pp <- sim_clustered_survival(n_per = 10L, G = 10L, K = 4L, seed = 91L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  ## Validation runs regardless of CI method (catches typos early).
  expect_error(
    contrast(
      fit,
      ivs,
      times = 1:4,
      type = "risk_difference",
      ci_method = "none",
      cluster = "nope"
    ),
    class = "survatr_bad_cluster"
  )
})

test_that("G = 2 clusters is accepted (degenerate boundary is G = 1)", {
  pp <- sim_clustered_survival(n_per = 30L, G = 2L, K = 5L, seed = 93L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  res <- contrast(
    fit,
    ivs,
    times = 1:5,
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )
  expect_false(any(is.na(res$estimates$se)))
  expect_true(all(res$estimates$se >= 0))
})

test_that("cluster bootstrap is reproducible under a fixed seed", {
  skip_on_cran()
  pp <- sim_clustered_survival(n_per = 15L, G = 20L, K = 5L, seed = 95L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  r1 <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "bootstrap",
    cluster = "site",
    n_boot = 50L,
    seed = 42L
  )
  r2 <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "bootstrap",
    cluster = "site",
    n_boot = 50L,
    seed = 42L
  )
  expect_equal(r1$estimates$se, r2$estimates$se)
  expect_equal(r1$estimates$ci_lower, r2$estimates$ci_lower)
})

## ---- 12. cluster_for_ids alignment unit (covers the defensive guard) --------

test_that("cluster_for_ids reindexes by name and guards missing ids", {
  expect_null(cluster_for_ids(NULL, c("a", "b")))
  labels <- c(a = "S1", b = "S2", c = "S1")
  ## Reindexed onto an arbitrary id order, by NAME not position.
  expect_equal(cluster_for_ids(labels, c("c", "a", "b")), c("S1", "S1", "S2"))
  ## An id absent from the label map is an internal row-order / id mismatch and
  ## must abort rather than form a phantom NA cluster in `rowsum()`.
  expect_error(
    cluster_for_ids(labels, c("a", "z")),
    class = "survatr_if_failed"
  )
})

## ---- 13. non-singleton clustered meat oracles (IPW / IPCW) ------------------

test_that("IPW: clustered sandwich SE agrees with the cluster bootstrap", {
  skip_on_cran()
  ## Non-singleton oracle for the two-stage clustered meat: the cluster
  ## bootstrap refits the propensity per replicate AND resamples whole sites, so
  ## agreement pins the clustered IPW stacked-EE sandwich (not just the
  ## cluster = id reduction).
  pp <- sim_clustered_confounded(n_per = 25L, G = 45L, K = 5L, seed = 72L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  s_sand <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )$estimates[get("intervention") == "a1"]$se
  s_boot <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "bootstrap",
    cluster = "site",
    n_boot = 300L,
    seed = 5L
  )$estimates[get("intervention") == "a1"]$se
  expect_equal(s_boot, s_sand, tolerance = 0.30)
})

test_that("IPCW: clustered sandwich SE agrees with the cluster bootstrap", {
  skip_on_cran()
  ## Non-singleton oracle for the three-block (beta + alpha + gamma) clustered
  ## meat. The cluster bootstrap refits BOTH the treatment and censoring models
  ## per replicate and resamples whole sites.
  pp <- sim_clustered_censoring(n_per = 25L, G = 45L, K = 5L, seed = 73L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    censoring = "C",
    estimator = "ipw",
    ipcw = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  s_sand <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "sandwich",
    cluster = "site"
  )$estimates[get("intervention") == "a1"]$se
  s_boot <- contrast(
    fit,
    ivs,
    times = c(3, 5),
    type = "risk",
    ci_method = "bootstrap",
    cluster = "site",
    n_boot = 300L,
    seed = 6L
  )$estimates[get("intervention") == "a1"]$se
  expect_equal(s_boot, s_sand, tolerance = 0.30)
})
