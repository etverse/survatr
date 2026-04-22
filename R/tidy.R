#' Tidy generic re-exported from `generics`
#'
#' @name tidy
#' @importFrom generics tidy
#' @export
NULL

#' Tidy a `survatr_result` into long format
#'
#' Stacks the per-intervention `estimates` and pairwise `contrasts` into a
#' single long `data.frame` suitable for downstream plotting or further
#' aggregation. The `estimand` column identifies which column the row
#' came from: `s_hat` / `risk_hat` / `rmst_hat` for intervention rows,
#' `risk_difference` / `risk_ratio` / `rmst_difference` for contrast rows.
#'
#' @param x A `survatr_result` from `contrast()`.
#' @param which One of `"estimates"` (per-intervention rows only),
#'   `"contrasts"` (pairwise rows only), or `"all"` (both stacked).
#'   Default `"all"`.
#' @param conf.int Logical; when `FALSE`, drop `ci_lower` / `ci_upper`
#'   columns. Default `TRUE`.
#' @param ... Unused.
#'
#' @return A `data.frame` in long format with columns
#'   `intervention | contrast | time | estimand | estimate | se |
#'   ci_lower | ci_upper` (or a subset when `which` / `conf.int`
#'   restrict it).
#' @method tidy survatr_result
#' @export
tidy.survatr_result <- function(
  x,
  which = c("all", "estimates", "contrasts"),
  conf.int = TRUE,
  ...
) {
  which <- match.arg(which)
  if (
    !is.logical(conf.int) ||
      length(conf.int) != 1L ||
      is.na(conf.int)
  ) {
    rlang::abort(
      "`conf.int` must be TRUE or FALSE.",
      class = "survatr_bad_conf_int"
    )
  }

  type <- x$type
  estimand_col <- switch(
    type,
    survival = "s_hat",
    risk = "risk_hat",
    risk_difference = "risk_hat",
    risk_ratio = "risk_hat",
    rmst = "rmst_hat",
    rmst_difference = "rmst_hat"
  )

  ## Per-intervention rows.
  est_long <- NULL
  if (which %in% c("all", "estimates")) {
    est_long <- data.table::data.table(
      intervention = x$estimates[["intervention"]],
      contrast = NA_character_,
      time = x$estimates[["time"]],
      estimand = estimand_col,
      estimate = x$estimates[[estimand_col]],
      se = x$estimates[["se"]],
      ci_lower = x$estimates[["ci_lower"]],
      ci_upper = x$estimates[["ci_upper"]]
    )
  }

  ## Pairwise contrasts.
  ctr_long <- NULL
  if (which %in% c("all", "contrasts") && nrow(x$contrasts) > 0L) {
    ctr_long <- data.table::data.table(
      intervention = NA_character_,
      contrast = x$contrasts[["contrast"]],
      time = x$contrasts[["time"]],
      estimand = type,
      estimate = x$contrasts[["estimate"]],
      se = x$contrasts[["se"]],
      ci_lower = x$contrasts[["ci_lower"]],
      ci_upper = x$contrasts[["ci_upper"]]
    )
  }

  out <- data.table::rbindlist(
    list(est_long, ctr_long),
    use.names = TRUE,
    fill = TRUE
  )
  if (!isTRUE(conf.int)) {
    out[, `:=`(ci_lower = NULL, ci_upper = NULL)]
  }
  as.data.frame(out)
}
