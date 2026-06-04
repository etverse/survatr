# Python reference fixtures (delicatessen)

These CSVs pin **independent** analytic-sandwich reference values from
[`delicatessen`](https://deli.readthedocs.io) (Zivich 2022) stacked
M-estimation, used to cross-validate survatr's variance against a second
implementation. The R tests read the committed CSVs; they do **not** call
Python at test time.

## Files

- `ipw_survival_data.csv` — shared person-period fixture (`id, t, A, L, Y,
  at_risk`). Generated in R from `sim_confounded_survival(n = 1200, K = 5,
  seed = 41, gamma = 0.8)`; `at_risk` is survatr's risk set. **Both** R and
  Python read this so the comparison is on identical data.
- `ipw_survival_delicatessen.csv` — delicatessen output (`quantity, time,
  estimate, se`) for `S1`, `S0`, and `RD` at times 2–5, consumed by
  `test-ipw-delicatessen.R`.
- `ice_survival_data.csv` — shared person-period fixture (`id, t, A, L, Y`) for
  Track B (longitudinal ICE). Generated in R from `sim_ice_feedback(n = 1500,
  K = 3, seed = 41)` (treatment-confounder feedback). Both R and Python read
  this.
- `ice_survival_delicatessen.csv` — delicatessen output (`quantity, time,
  estimate, se`) for `S1`, `S0`, and `RD` at times 1–3, consumed by
  `test-ice-survival-delicatessen.R`. The Python script reimplements the
  backward survival-tail ICE recursion as a stacked M-estimator, so its
  numerical sandwich independently validates survatr's analytic survival IF
  chain — in particular the `(1 - D_k)` failure carry-forward factor.

## Environment

`delicatessen` is not an R dependency and is not installed at test time. The
generator scripts run in causatr's Python venv:

```
../../../../causatr/data-raw/zepid_venv/bin/python   # delicatessen 4.1, numpy 2.4.4, scipy 1.17.1
```

(any environment with `delicatessen`, `numpy`, `scipy`, `pandas` works).

## Regenerating

From the package root:

```sh
# (1) regenerate the shared data fixture (R), if the DGP changes:
Rscript -e 'devtools::load_all(); source("tests/testthat/helper-dgp.R"); \
  dt <- sim_confounded_survival(n=1200L, K=5L, seed=41L, gamma=0.8); \
  data.table::setkeyv(dt, c("id","t")); \
  dt[, at_risk := as.integer(data.table::shift(cumsum(Y),1L,fill=0L,"lag")==0L), by=id]; \
  data.table::fwrite(dt[,.(id,t,A,L,Y,at_risk)], \
    "tests/testthat/fixtures/python/ipw_survival_data.csv")'

# (2) regenerate the delicatessen reference (Python):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_ipw_survival.py

# Track B (ICE) fixtures:
#   (1) data fixture (R):
Rscript -e 'source("tests/testthat/helper-ice-survival-oracle.R"); \
  dt <- sim_ice_feedback(n=1500L, K=3L, seed=41L); \
  data.table::fwrite(dt[,.(id,t,A,L,Y)], \
    "tests/testthat/fixtures/python/ice_survival_data.csv")'
#   (2) delicatessen reference (Python):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_ice_survival.py
```

The R fit and the delicatessen stack agree to ~1e-4 on `S^a(t)` and the risk
difference (point and sandwich SE).
