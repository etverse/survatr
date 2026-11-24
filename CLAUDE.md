# survatr

survatr estimates causal effects on time-to-event outcomes in cohort data: survival and risk curves, restricted mean survival and time lost, quantiles and competing-risks cumulative incidence.

It estimates them by pooled-logistic hazard g-computation, a weighted hazard marginal structural model (IPW) and iterated conditional expectation (ICE) on person-period data, with sandwich, clustered and bootstrap variance. It builds on causatr's influence-function engine and intervention constructors.

It is part of the etverse, whose design (https://etverse.github.io/design.html) fixes what survatr owns and the interface every package shares. survatr does not follow that design yet; until it is brought into line, this file and the package documentation describe what the code does today.

## API

```r
fit <- surv_fit(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                id = "id", time = "time", censoring = "C",
                estimator = "gcomp")   # or "ipw", "ice"

res <- contrast(fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times = seq(0, 120, by = 12),
  type = "risk_difference",            # see R/estimand_registry.R
  ci_method = "sandwich")              # or "bootstrap", "none"
```

Competing risks use `surv_fit(competing = <cause column>)`; IPCW uses `ipcw =`; `contrast(cluster = <column>)` gives cluster-robust intervals.

S3 classes: `survatr_fit` (from `surv_fit()`), `survatr_result` (from `contrast()`, with time-indexed `estimates` and `contrasts` tables and `time_grid`) and `survatr_diag` (from `diagnose()`).

## Code layout

Files in `R/` are named by what they hold; `tests/testthat/test-foo.R` tests `R/foo.R`.

- Fit: `surv_fit.R`, `prepare_data.R`, `risk_set.R`, `checks.R`
- Estimators: `gcomp_survival.R`, `ipw_survival.R`, `ipw_weights_survival.R`, `ice_survival.R`, `ipcw_survival.R`, `competing_risks.R`
- Contrast and estimands: `contrast.R`, `contrast_types.R`, `contrast_validators.R`, `estimand_registry.R`, `estimands_quantile.R`, `rmst.R`, `survival_curve.R`, `predict_hazard.R`, `apply_intervention.R`
- Variance: `variance_if_survival.R` (the delta method across time), `variance_if_ipw_survival.R`, `variance_if_ice_survival.R`, `variance_if_competing.R`, `variance_sandwich.R`, `variance_cluster.R`, `variance_bootstrap*.R`
- S3 methods and plots: `print.R`, `plot.R`, `tidy.R`, `forrest.R`, `diagnose.R`
- Other: `utils.R`, `data.R`, `survatr-package.R`, `zzz.R`

Invariants that look like bugs and are not, and the axes a combination audit runs over, are in `.claude/rules/`; each file loads when you open the code it covers.

## Commands

The `Makefile` has `make document`, `make test`, `make check`, `make format` (Air), `make lint` and `make site` (altdoc). While iterating, run one file:

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-rmst.R")'
```

Slow tests call `skip_on_cran()`; set `NOT_CRAN=true` to run them. Python oracle scripts are in `data-raw/` and write to `tests/testthat/fixtures/python/`.

## Planning files

`SURVIVAL_PACKAGE_HANDOFF.md` holds the scope, the reasons for a separate package and the oracle list; `CHUNK_*.md` hold the design of each feature; `FEATURE_COVERAGE_MATRIX.md` records which combinations are tested.

## Dependencies

Imports: boot, causatr (GitHub remote `etverse/causatr`), data.table, generics, numDeriv, rlang, sandwich. survatr calls unexported causatr functions through `causatr:::`, so a change in causatr can break it; `test-causatr-integration.R` pins the signatures. Oracles in Suggests include lmtp (point estimates only), gfoRmula, survival, riskRegression and survRM2.
