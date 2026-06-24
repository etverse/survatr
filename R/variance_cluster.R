#' Reindex name-keyed cluster labels onto an IF matrix's row order
#'
#' @description
#' `validate_cluster()` returns cluster labels keyed by id (as character),
#' independent of any one influence-function matrix's row order. Each variance
#' path builds its `n_ids x |t-grid|` IF matrix with rows in its own
#' `unique(id)` order; this helper looks the labels up by id name so the
#' returned vector lines up row-for-row with that IF matrix.
#'
#' @details
#' Alignment is by id **name**, not position, so it is robust to whichever
#' first-appearance ordering a given path uses (the single-event sandwich, the
#' competing-risks sandwich, and the quantile assembly each derive `unique_ids`
#' separately). Every id in `unique_ids` must be present in `labels` — a missing
#' id signals an internal row-order / id mismatch and aborts loudly rather than
#' silently dropping a row from the within-cluster sum.
#'
#' @param labels Named character vector from `validate_cluster()` (id -> cluster
#'   label), or `NULL` for the per-individual sandwich.
#' @param unique_ids The id values in the IF matrix's row order.
#'
#' @returns `NULL` when `labels` is `NULL`, otherwise an unnamed character
#'   vector of cluster labels, one per row of the IF matrix.
#' @family variance
#' @noRd
cluster_for_ids <- function(labels, unique_ids) {
  if (is.null(labels)) {
    return(NULL)
  }
  out <- labels[as.character(unique_ids)]
  ## A missing id means the cluster map and the IF rows disagree on the id set
  ## -- a bug in the calling path, not user error. Surface it as an internal
  ## variance failure rather than letting `rowsum()` form a phantom NA cluster.
  if (anyNA(out)) {
    rlang::abort(
      "Cluster labels do not cover every id in the influence-function matrix.",
      class = "survatr_if_failed"
    )
  }
  unname(out)
}

#' Pointwise SE from a survival IF matrix, optionally cluster-robust
#'
#' @description
#' Reduce an `n_ids x k` per-individual influence-function matrix to a
#' length-`k` vector of standard errors. With `cluster = NULL` this is the
#' ordinary sandwich `sqrt(diag(crossprod(if_mat)) / n_ids^2)`; with `cluster`
#' supplied it is the cluster-robust generalization that sums each individual's
#' IF row into its cluster's row before the cross-product.
#'
#' @details
#' The cluster-robust covariance keeps the **same `n_ids^2` divisor** and only
#' changes the meat: the unit of independence becomes the cluster, so the
#' between-cluster outer products replace the between-individual ones,
#'
#' \deqn{V = \frac{1}{n^2} \sum_{g} \Big(\sum_{i \in g} IF_i\Big)
#'           \Big(\sum_{i \in g} IF_i\Big)^\top.}
#'
#' This matches `sandwich::vcovCL(type = "HC0", cadjust = FALSE)` and reduces
#' exactly to the per-individual estimate when every cluster is a singleton
#' (`cluster = id`, `G = n_ids`). The inner row-sum-then-cross-product is
#' delegated to `causatr:::vcov_from_if()` (causatr-reuse audit), the one shared
#' IF-reduction primitive: survatr keeps the survival-specific gradients (the
#' RMST trapezoidal weights, the log-ratio scaling) and hands the final
#' `crossprod -> V` step to causatr so the per-individual and clustered paths
#' share one tested implementation. `k` is small here (the time grid), so
#' forming the full `k x k` covariance just to read its diagonal is cheap.
#'
#' @param if_mat The `n_ids x k` influence-function matrix (a single column for
#'   the scalar quantile IF).
#' @param n_ids Number of individuals (the divisor; unchanged by clustering).
#' @param cluster `NULL` for the per-individual sandwich, or a length-`n_ids`
#'   vector of cluster labels aligned to the rows of `if_mat` (via
#'   `cluster_for_ids()`).
#'
#' @returns Numeric vector of length `k`: the pointwise SEs.
#' @family variance
#' @seealso `cluster_for_ids()`, `validate_cluster()`
#' @noRd
clustered_pointwise_se <- function(if_mat, n_ids, cluster = NULL) {
  ## `asplit(if_mat, 2)` hands `vcov_from_if()` its expected list-of-IF-vectors
  ## (one per column / time point); it cbinds them back, applies the optional
  ## within-cluster `rowsum()`, and returns `crossprod / n_ids^2`.
  v <- causatr:::vcov_from_if(
    IF_list = asplit(if_mat, 2L),
    n = n_ids,
    int_names = NULL,
    cluster = cluster
  )
  sqrt(pmax(diag(v), 0))
}
