# survatr roadmap

survatr estimates causal effects on time-to-event outcomes in cohort data: survival and risk curves, restricted mean survival time (RMST), restricted mean time lost and competing risks, by pooled-logistic hazard g-computation, IPW and iterated conditional expectation (ICE). The [etverse design](https://etverse.github.io/design.html) fixes its scope and interface. Each item's status is in the [issues](https://github.com/etverse/survatr/issues), not in this file.

## 0.1.0

The first release comes once every etverse estimator has been checked against independent software, the bugs those checks find are fixed, and the package follows the design. Its work is the [0.1.0 milestone](https://github.com/etverse/survatr/milestone/1):

- Integrate RMST over every period, not only the requested times ([#15](https://github.com/etverse/survatr/issues/15)).
- Reject or define `stochastic()` and `ipsi()` under g-computation ([#16](https://github.com/etverse/survatr/issues/16)).
- Decide whether survival ICE stays in survatr or is built on dynamatr ([#17](https://github.com/etverse/survatr/issues/17)).
- Rename the "Track A/B" labels left in a condition class, two internal functions and the fit's `track` field ([#34](https://github.com/etverse/survatr/issues/34)).

## Estimators

- Left truncation and delayed entry ([#18](https://github.com/etverse/survatr/issues/18)).
- Doubly robust (AIPW) survival ([#19](https://github.com/etverse/survatr/issues/19)).
- IPCW for g-computation and ICE ([#29](https://github.com/etverse/survatr/issues/29)).
- Competing risks under IPW and ICE ([#30](https://github.com/etverse/survatr/issues/30)).

## IPW and interventions

- Continuous, categorical and count treatments ([#22](https://github.com/etverse/survatr/issues/22)), multivariate treatments ([#25](https://github.com/etverse/survatr/issues/25)) and external or survey weights ([#24](https://github.com/etverse/survatr/issues/24)).
- Incremental propensity score interventions ([#23](https://github.com/etverse/survatr/issues/23)) and stochastic interventions ([#26](https://github.com/etverse/survatr/issues/26)).
- Subgroup curves with `by =` ([#31](https://github.com/etverse/survatr/issues/31)).

## Variance

- Simultaneous confidence bands ([#20](https://github.com/etverse/survatr/issues/20)).
- A numeric variance fallback for models without an analytic bread ([#28](https://github.com/etverse/survatr/issues/28)).
- Multi-way clustering and small-cluster corrections ([#32](https://github.com/etverse/survatr/issues/32)).
- Multiple imputation for missing predictors, which needs research before code ([#27](https://github.com/etverse/survatr/issues/27)).

## Open scope questions

- Recurrent events and multi-state models ([#21](https://github.com/etverse/survatr/issues/21)).
- Instrumental-variable estimation by structural nested cumulative failure time models ([#33](https://github.com/etverse/survatr/issues/33)).

## Not planned here

- Target-trial emulation and landmarking: [emulatr](https://github.com/etverse/emulatr/issues/1).
- Time-varying treatment weights: [dynamatr](https://github.com/etverse/dynamatr/issues/6).
- Survival in sampled designs: matchatr. Separable effects: separatr.
- Fine-Gray subdistribution hazards.
