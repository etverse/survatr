
<!-- README.md is generated from README.Rmd. Please edit that file -->

# survatr

<!-- badges: start -->

[![R-CMD-check](https://github.com/etverse/survatr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/etverse/survatr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/etverse/survatr/graph/badge.svg)](https://app.codecov.io/gh/etverse/survatr)
<!-- badges: end -->

**survatr** extends the [causatr](https://github.com/etverse/causatr)
engine to **time-to-event outcomes**. It estimates counterfactual
survival via pooled-logistic discrete-time hazard models on
person-period data, with three complementary estimators: g-computation
(Track A), inverse-probability weighting (IPW, Track A), and iterated
conditional expectations (ICE) for longitudinal, time-varying treatments
(Track B). The **estimand is a curve** — survival $S^a(t)$, risk
$1 - S^a(t)$, risk difference, risk ratio, or restricted mean survival
time (RMST) — not a scalar. Inference is by a delta-method sandwich that
propagates per-row influence functions through the cumulative-product
survival curve, or by an individual-level nonparametric bootstrap.

The package implements the survival methods of Hernán & Robins (2025)
*Causal Inference: What If* (Chapter 17) with the same two-step API
causatr users know:

1.  **Fit** the hazard model once with `surv_fit()`
2.  **Contrast** many interventions cheaply with `contrast()`

## Installation

You can install the development version of survatr from
[GitHub](https://github.com/etverse/survatr) with:

``` r
# install.packages("pak")
pak::pak("etverse/survatr")
```

## Quick example

The data-generating process below has a constant 6% discrete-time
hazard, a single binary treatment with no causal effect, and one
confounder `L` that drives both treatment and the hazard. Under this DGP
the counterfactual survival under either intervention converges to
$(1 - 0.06)^t$.

``` r
library(causatr) # attach first so survatr's contrast() generic shadows it
library(survatr)
library(data.table)

# Simulate person-period data: constant hazard, one confounder L.
set.seed(1)
n <- 800L
K <- 10L
L <- rnorm(n)
A <- rbinom(n, 1L, plogis(0.3 * L)) # L affects treatment
pp <- rbindlist(lapply(seq_len(n), function(i) {
  Y <- rbinom(K, 1L, 0.06) # hazard independent of A (null effect)
  first <- which(Y == 1L)[1L]
  if (!is.na(first) && first < K) Y[(first + 1L):K] <- 0L
  data.table(id = i, t = seq_len(K), A = A[i], L = L[i], Y = Y)
}))

# Step 1: fit the pooled-logistic hazard model once.
fit <- surv_fit(
  pp,
  outcome      = "Y",
  treatment    = "A",
  confounders  = ~ L,
  id           = "id",
  time         = "t",
  time_formula = ~ factor(t) # non-parametric baseline hazard
)

# Step 2: contrast interventions on the fitted hazard. Risk difference with
# a delta-method sandwich CI at each of the 10 periods.
res <- contrast(
  fit,
  interventions = list(a1 = static(1), a0 = static(0)),
  times         = 1:10,
  type          = "risk_difference",
  reference     = "a0",
  ci_method     = "sandwich"
)
res$contrasts[] # `[]` forces data.table to auto-print a function-returned table
#>     contrast  time    estimate          se     ci_lower   ci_upper
#>       <char> <int>       <num>       <num>        <num>      <num>
#>  1: a1 vs a0     1 0.001650894 0.005623608 -0.009371174 0.01267296
#>  2: a1 vs a0     2 0.003298391 0.011214573 -0.018681768 0.02527855
#>  3: a1 vs a0     3 0.004732117 0.016070436 -0.026765359 0.03622959
#>  4: a1 vs a0     4 0.005742165 0.019490983 -0.032459460 0.04394379
#>  5: a1 vs a0     5 0.007329630 0.024885202 -0.041444470 0.05610373
#>  6: a1 vs a0     6 0.008164440 0.027713630 -0.046153277 0.06248216
#>  7: a1 vs a0     7 0.008834524 0.029982957 -0.049930992 0.06760004
#>  8: a1 vs a0     8 0.009416127 0.031949858 -0.053204443 0.07203670
#>  9: a1 vs a0     9 0.009897176 0.033579818 -0.055918058 0.07571241
#> 10: a1 vs a0    10 0.010389809 0.035247492 -0.058694005 0.07947362
```

The DGP has no treatment effect, so the risk-difference CI comfortably
covers 0 at every time point.

## Features

- **Three estimators**: g-computation (`estimator = "gcomp"`,
  pooled-logistic standardization), IPW (`estimator = "ipw"`, a weighted
  marginal hazard MSM with stabilized density-ratio weights), and ICE
  (`estimator = "ice"`, Track B longitudinal survival via backward
  iterated conditional expectations for a time-varying treatment).
  Matching + survival is a hard reject — use
  `survival::coxph(..., weights, cluster)` directly.
- **Curve-valued estimands**: survival $S^a(t)$, risk $1 - S^a(t)$, risk
  difference, risk ratio (log-scale CI), and RMST / RMST difference, all
  returned as time-indexed `data.table`s.
- **Flexible interventions** via causatr: `static()`, `shift()`,
  `scale_by()`, `threshold()`, and `dynamic()`. Which are available
  depends on the estimator.
- **Robust inference**: a delta-method sandwich that aggregates
  per-individual influence functions across time through the cumulative
  product (pointwise Wald bands; stacked estimating equations for IPW
  and ICE), or a nonparametric bootstrap that resamples ids and refits
  per replicate.
- **Pluggable hazard models** via `model_fn` (`stats::glm` by default;
  `mgcv::gam` with an `s(time)` baseline-hazard term, splines via
  `ns()`).
- **`tidy()` / `plot()` / `print()` / `forrest()`** S3 methods that
  dispatch on the time-indexed `survatr_result` shape (survival curves,
  and forest plots at a reference time `t*`).

## Learning more

- `vignette("survatr")` — introduction and the two-step API
- `vignette("ipw")` — Track A IPW (weighted hazard MSM)
- `vignette("ice")` — Track B longitudinal survival via ICE

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC. Chapter 17 (survival analysis via IP weighting and
standardization).

Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
sandwich variance estimator for iterated conditional expectation
g-computation. *Statistics in Medicine* 43:5562–5572.

## Acknowledgements

This package was built with the contribution of
[Claude](https://claude.ai), Anthropic’s AI assistant.
