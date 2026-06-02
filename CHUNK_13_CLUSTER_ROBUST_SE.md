# Chunk 13 — Cluster-robust sandwich variance

> **Status: ⬜ Not started**
> **Depends on:** Chunk 3 (sandwich IF: the `n × |t-grid|` per-individual IF
> matrix).
> **Oracle:** `cluster = id` must reproduce the current per-individual SE
> exactly (degenerate-to-current check); a known multi-site DGP where ignoring
> clustering under-covers and the clustered SE restores nominal coverage.

## Goal

Allow a `cluster =` argument so the sandwich variance is robust to within-
cluster correlation (multi-site studies, repeated enrolment, household /
provider clustering). The per-individual influence function is already the
variance currency in survatr: clustering is a single aggregation step on the
existing `n × |t-grid|` IF matrix — **sum the IF rows within cluster before
the crossproduct**. No new estimation, no change to the fit or contrast path.
This is a short doc precisely because the IF spine makes it near-free.

## The math

Chunk 3 stores the per-individual influence function as an `n × |t-grid|`
matrix `IF`, and the time-covariance is `V = crossprod(IF) / n²` (so pointwise
SE is `sqrt(diag(V))`, RMST SE is `aᵀ V a`). The standard sandwich treats
individuals as independent. Under clustering, the unit of independence is the
**cluster**, not the individual. The cluster-robust covariance sums each
individual's IF row into its cluster's IF row, then takes the crossproduct over
the `G` cluster-level rows:

```
IF_g = Σ_{i ∈ cluster g} IF_i          (row sum within cluster → G × |t-grid|)
V_clustered = crossprod(IF_g) / G²      (G = number of clusters)
```

When every cluster is a singleton (`cluster = id`), `IF_g ≡ IF_i` and `G = n`,
so `V_clustered = crossprod(IF) / n² = V` — **identical to the current
estimate**. This is the regression invariant: the cluster-robust path is a
strict generalization that reduces to the shipped estimate at the finest
partition.

The aggregation commutes with all the downstream functionals: pointwise SE is
`sqrt(diag(V_clustered))`, RMST / RMTL SE is `aᵀ V_clustered a`, contrast
(difference / ratio) IFs are differenced row-wise **before** the cluster sum (a
contrast IF is a per-individual quantity, then aggregated to clusters exactly
like the level IF). The number of independent units drops from `n` to `G`, so
the normalization changes from `n²` to `G²` — this is the entire source of the
wider interval.

## Design

- A new `cluster =` argument on `contrast()` (a column name in the fit's
  person-period data, constant within id — validated). Default `NULL` →
  unchanged per-individual sandwich.
- The cluster id is collapsed to one value per individual (it is constant
  within id), then used as the grouping key for the row sum on the
  per-individual IF matrix.
- The change lives entirely in the variance assembly: after the per-individual
  IF matrix is built (level and contrast), `rowsum()` by cluster, then
  `crossprod() / G²`. Bootstrap gets a sibling change — resample **clusters**
  (all members together), not individuals.

## Deliverables

### New R files
- `R/variance_cluster.R` — the cluster-aggregation helper: validate the
  cluster column (constant within id, no NA), collapse to per-id cluster
  labels, `rowsum()` the `n × |t-grid|` IF matrix to `G × |t-grid|`, return the
  clustered time-covariance. One function reused by the level and contrast
  paths.

### Updated R files
- `R/contrast.R` — thread a `cluster =` argument through to the variance
  assembly; pass the per-id cluster labels alongside the IF matrix.
- `R/variance_sandwich.R` — when `cluster` is supplied, route `crossprod(IF)`
  through `variance_cluster.R` (`rowsum` then `/ G²`) instead of `/ n²`.
- `R/variance_bootstrap.R` — cluster-resampling option: when `cluster` is set,
  resample clusters (all person-period rows of all members together).
- `R/contrast_validators.R` — `validate_cluster()`: column exists, constant
  within id (`survatr_cluster_varies_within_id`), no NA
  (`survatr_cluster_na`), `G ≥ 2` (a single cluster is meaningless →
  `survatr_cluster_degenerate`).

### Tests (`tests/testthat/`)
- `test-variance-cluster.R` — `cluster = id` reproduces the current
  per-individual SE to machine tolerance (the degenerate-to-current invariant);
  a multi-site DGP with strong within-cluster correlation under-covers with the
  per-individual SE and recovers nominal coverage with the clustered SE; the
  clustered SE is wider than (or equal to) the per-individual SE on positively
  correlated clusters.

## API contract

```r
result <- contrast(
  fit,
  interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
  times = seq(0, 120, 12),
  type = "risk_difference",
  ci_method = "sandwich",
  cluster = "site"            # NEW — column name; constant within id
)
# cluster = NULL  -> current per-individual sandwich (unchanged)
# cluster = "id"  -> identical to per-individual sandwich (G = n)
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **The per-individual IF is already the variance currency.** Clustering is one
  aggregation step (`rowsum` within cluster) on the existing `n × |t-grid|`
  matrix — no refit, no new estimation, no change to the point estimate.
- **`cluster = id` must reproduce the current SE exactly** (singleton clusters,
  `G = n`). This is the load-bearing regression invariant.
- **Normalization changes from `n²` to `G²`** (independent units = clusters);
  this is the entire mechanism of the wider interval.
- **Contrast IFs are differenced per-individual first, then summed within
  cluster** — the difference is a per-individual quantity; aggregating after
  the difference is correct, aggregating before is wrong.
- **Cluster id must be constant within id** (an individual belongs to one
  cluster); validated upfront with a classed abort.
- **Bootstrap under clustering resamples clusters**, all members together — the
  same unit-of-independence change as the sandwich.

## Non-goals (deferred)
- Multi-way / nested clustering (e.g. site × period) — single clustering
  dimension in this chunk.
- Small-`G` corrections (cluster-count bias adjustments) — report `G` so the
  user can judge; a degrees-of-freedom correction is a possible later
  enhancement.
- Clustering as a modelling choice (random effects in the hazard model) — out
  of scope; this is variance-only, the point estimate is unchanged.

## Dependencies & composition
- Chunk 3 only (the IF matrix). Composes transparently with chunks 5/7/11/12:
  any estimand whose variance flows through the `n × |t-grid|` IF matrix gains
  the clustered SE for free, because the aggregation is on the IF, not the
  functional.

## Acceptance checklist
- [ ] `cluster = id` reproduces the per-individual SE to machine tolerance.
- [ ] On a correlated multi-site DGP, per-individual SE under-covers and the
      clustered SE recovers nominal coverage; clustered SE ≥ per-individual SE.
- [ ] Cluster validation aborts fire (varies-within-id, NA, single cluster).
- [ ] Bootstrap resamples clusters when `cluster` is set.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
