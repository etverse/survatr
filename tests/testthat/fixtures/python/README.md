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
- `gcomp_survival_delicatessen.csv` — delicatessen output (`quantity, time,
  estimate, se`) for the g-computation sandwich: `S1` / `S0`, `RD`, `RR`,
  `RMST1` / `RMST0`, and `RMSTdiff` at times 2–5, consumed by
  `test-gcomp-delicatessen.R`. **Reuses `ipw_survival_data.csv`** — gcomp and
  IPW target the same counterfactual `S^a(t)` on the same confounded data
  (gcomp adjusts via the outcome model, with `L` in the hazard design; IPW via
  weights). Pins the foundational Track A sandwich (`estimator = "gcomp"`)
  against an independent M-estimator across the full estimand surface.
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
- `cr_survival_data.csv` — shared person-period fixture (`id, t, A, L, event,
  at_risk`) for competing risks. `event` is the multi-valued cause column
  (`0` = no event, `1` / `2` = cause); `at_risk` is the shared all-cause risk
  set. Generated in R from `sim_two_cause_constant_hazard(n = 1200, K = 5,
  seed = 41, beta1_A = -0.5, beta1_L = 0.6, beta2_L = 0.4, gamma = 0.7)`. Both
  R and Python read this.
- `cr_survival_delicatessen.csv` — delicatessen output (`quantity, time,
  estimate, se`) for per-cause CIF `F{j}_a{a}`, all-cause `S_a{a}`, and the
  cause-`j` CIF difference `RD{j}` at times 2 and 4, consumed by
  `test-competing-risks-sandwich.R`. The Python script stacks the two
  cause-specific hazard scores plus the CIF / survival functionals, so its
  sandwich independently validates survatr's cause-specific CIF IF — in
  particular the all-cause `(1 - H)` sensitivity denominator.
- `ipcw_survival_data.csv` — shared person-period fixture (`id, t, A, L, Y, C,
  at_risk, cens_at_risk`) for IPCW (chunk 11). `C` is the per-period censoring
  indicator; `at_risk` is the MSM fit mask (`prev_event==0 & prev_cens==0 &
  C==0`); `cens_at_risk` is the censoring-model fit mask (`prev_event==0 &
  prev_cens==0`, including the `C=1` rows). Generated in R from
  `sim_informative_censoring(n=2000, K=5, seed=42)` — informative censoring
  with `delta_cens = 0.8`. Both R and Python read this.
- `ipcw_survival_delicatessen.csv` — delicatessen output (`quantity, time,
  estimate, se`) for `S1`, `S0`, and `RD` at times 2–5, consumed by
  `test-ipcw-delicatessen.R`. The Python script implements the full three-block
  stacked EE (treatment propensity + censoring-model denominator + weighted
  hazard MSM), independently validating survatr's `A_beta_gamma` cross-
  derivative and the `n_ids / n_cens_fit` bread scaling.

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

# gcomp sandwich reference (Python; reuses ipw_survival_data.csv):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_gcomp_survival.py

# Track B (ICE) fixtures:
#   (1) data fixture (R):
Rscript -e 'source("tests/testthat/helper-ice-survival-oracle.R"); \
  dt <- sim_ice_feedback(n=1500L, K=3L, seed=41L); \
  data.table::fwrite(dt[,.(id,t,A,L,Y)], \
    "tests/testthat/fixtures/python/ice_survival_data.csv")'
#   (2) delicatessen reference (Python):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_ice_survival.py

# Competing-risks fixtures:
#   (1) data fixture (R):
Rscript -e 'source("tests/testthat/helper-cr-oracle.R"); library(data.table); \
  dt <- sim_two_cause_constant_hazard(n=1200L, K=5L, h1=0.10, h2=0.07, \
    seed=41L, beta1_A=-0.5, beta1_L=0.6, beta2_L=0.4, gamma=0.7); \
  setkeyv(dt, c("id","t")); \
  dt[, at_risk := as.integer(shift(cumsum(event!=0),1L,fill=0L,"lag")==0L), by=id]; \
  fwrite(dt[,.(id,t,A,L,event,at_risk)], \
    "tests/testthat/fixtures/python/cr_survival_data.csv")'
#   (2) delicatessen reference (Python):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_competing_risks.py
```

# IPCW (chunk 11) fixtures:
#   (1) data fixture (R):
Rscript -e 'devtools::load_all(); dt <- sim_informative_censoring(n=2000L, K=5L, seed=42L); \
  dt[, prev_event := data.table::shift(cumsum(Y),1L,fill=0L), by=id]; \
  dt[, prev_cens := data.table::shift(cumsum(C),1L,fill=0L), by=id]; \
  dt[, at_risk := as.integer(prev_event==0 & prev_cens==0 & C==0)]; \
  dt[, cens_at_risk := as.integer(prev_event==0 & prev_cens==0)]; \
  dt[, c("prev_event","prev_cens") := NULL]; \
  data.table::fwrite(dt, "tests/testthat/fixtures/python/ipcw_survival_data.csv")'
#   (2) delicatessen reference (Python):
../causatr/data-raw/zepid_venv/bin/python data-raw/delicatessen_ipcw_survival.py
```

The R fit and the delicatessen stack agree to ~1e-4 on `S^a(t)` and the risk
difference (point and sandwich SE).
