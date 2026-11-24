---
paths:
  - "R/variance_cluster.R"
  - "R/variance_bootstrap*.R"
  - "R/variance_sandwich.R"
  - "tests/testthat/test-variance-cluster.R"
  - "tests/testthat/test-bootstrap-*.R"
---

# Clustered and bootstrap variance

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The cluster-robust divisor is n² (individuals), not G². `contrast(cluster = <col>)` sums each person's influence-function row into its cluster and computes `crossprod(rowsum(IF, cluster)) / n²`: the estimator is a person-level mean, so clustering changes only the meat. It matches `sandwich::vcovCL(type = "HC0", cadjust = FALSE)` to 1e-10 and is the divisor of `causatr:::vcov_from_if()`. G² inflates the variance by (n/G)².
- `cluster = id` reproduces the per-individual SE to 1e-12 for every estimand; if that breaks, the aggregation is wrong.
- One helper does the reduction everywhere: `clustered_pointwise_se(if_mat, n_ids, cluster)`, called by `fill_sandwich_ses()`, `fill_sandwich_ses_cr()` and `assemble_quantile_result()`. With `cluster = NULL` it equals `sqrt(diag(crossprod(IF)) / n_ids²)`.
- Cluster labels are resolved once by `validate_cluster()` into an id-to-label vector and aligned to each path's row order by id name with `cluster_for_ids()`; a missing id raises `survatr_if_failed`. Validation errors: `survatr_cluster_varies_within_id`, `survatr_cluster_na`, `survatr_cluster_degenerate` (fewer than two clusters), `survatr_bad_cluster`.
- With treatment randomised within cluster, a shared cluster effect cancels in a difference, so the SE of a contrast can stay at or below the per-individual SE while the SE of a level widens.
- The bootstrap resamples individuals (all their person-period rows) and refits the hazard model per replicate; with `cluster` set it resamples whole clusters, giving each redrawn cluster new bootstrap ids.
- Under a parallel backend the bootstrap sets `L'Ecuyer-CMRG` before `set.seed()` and restores the caller's RNG kind on exit, since `mclapply` and `parLapply` ignore a plain seed otherwise.
