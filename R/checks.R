## Column names survatr reserves for internal bookkeeping. Guarded upfront so
## user-data collisions produce an informative error rather than silent overwrite.
## `.survatr_prev_event` / `.survatr_prev_cens`: lagged within-id cumsums
## written by `build_risk_set()` and stripped from `fit$pp_data` before
## return.
## `.cf_hazard` / `.cf_surv`: temporary columns written by
## `compute_survival_curve()` and `compute_survival_if_matrix()` onto a
## copy of `fit$pp_data` to carry per-row counterfactual hazards and the
## within-id cumulative product. They live only on the internal copy and
## never reach user code, but we reserve the names so a user input with
## those columns is caught upfront rather than shadowed inside a copy.
## `.survatr_any_event` / `.survatr_cause_event`: transient all-cause and
## per-cause 0/1 event indicators derived from the multi-valued competing column
## in the cause-specific hazards path. Like the lagged cumsums they live only on
## internal copies, but the names are reserved so a user input carrying them is
## caught upfront rather than shadowed inside a copy.
SURVATR_RESERVED_COLS <- c(
  ".survatr_prev_event",
  ".survatr_prev_cens",
  ".cf_hazard",
  ".cf_surv",
  ".survatr_any_event",
  ".survatr_cause_event",
  ".cf_all_haz",
  ".cf_surv_lag",
  ".cf_cif"
)

#' Validate external weights
#'
#' Copy-adapted from `causatr:::check_weights`. Accepts `NULL` (no weights).
#' Otherwise: must be numeric, length equal to `n`, finite, non-negative.
#' `NA` / `Inf` / `NaN` / negative values are rejected. Zero weights are
#' allowed (they act as row exclusions in GLM fits).
#'
#' @param weights Numeric vector, or `NULL`.
#' @param n Expected length (`nrow(data)` at the call site).
#' @param call Enclosing frame for the error signal; passed to `rlang::abort()`.
#'
#' @return Invisibly `NULL`. Called for its side effect (error or pass-through).
#' @family checks
#' @noRd
check_weights <- function(weights, n, call = rlang::caller_env()) {
  if (is.null(weights)) {
    return(invisible(NULL))
  }
  if (!is.numeric(weights)) {
    rlang::abort(
      "`weights` must be numeric.",
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (length(weights) != n) {
    rlang::abort(
      paste0(
        "`weights` must have length equal to `nrow(data)` (",
        n,
        "), got ",
        length(weights),
        "."
      ),
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (anyNA(weights)) {
    rlang::abort(
      paste0(
        "`weights` contains ",
        sum(is.na(weights)),
        " missing value(s). ",
        "Drop those rows or impute before calling `surv_fit()`."
      ),
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (any(!is.finite(weights))) {
    rlang::abort(
      "`weights` contains non-finite value(s) (Inf / NaN).",
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (any(weights < 0)) {
    rlang::abort(
      "`weights` must be non-negative.",
      class = "survatr_bad_weights",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the IPW weight-truncation argument
#'
#' Accepts `NULL` (no truncation) or a numeric scalar in `(0, 1]`. The `trim`
#' value is the quantile at which the per-id stabilized weights are winsorized
#' (Cole & Hernán 2008); `1` means no truncation. Out-of-range, non-numeric,
#' non-finite, or non-scalar values are rejected upfront.
#'
#' @param trim Numeric scalar in `(0, 1]`, or `NULL`.
#' @param call Enclosing frame for the error signal.
#'
#' @returns Invisibly `NULL`. Called for its side effect (error or pass-through).
#' @family checks
#' @noRd
check_trim <- function(trim, call = rlang::caller_env()) {
  if (is.null(trim)) {
    return(invisible(NULL))
  }
  ok <- is.numeric(trim) &&
    length(trim) == 1L &&
    is.finite(trim) &&
    trim > 0 &&
    trim <= 1
  if (!ok) {
    rlang::abort(
      c(
        "`trim` must be `NULL` or a single number in (0, 1].",
        i = paste0(
          "`trim` is the upper quantile at which to winsorize IPW ",
          "weights; 1 means no truncation."
        )
      ),
      class = "survatr_bad_trim",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the IPCW argument and its required context
#'
#' Rejects unsupported combinations of `ipcw` with other `surv_fit()`
#' arguments:
#' - IPCW is only supported with `estimator = "ipw"` (the censoring weight
#'   multiplies the treatment weight in the weighted MSM).
#' - IPCW requires a non-`NULL` `censoring` column (the censoring indicator
#'   is the response for the censoring hazard model).
#' - `ipcw` itself must be a one-sided formula or `NULL`.
#'
#' @param ipcw One-sided formula or `NULL`.
#' @param estimator Character scalar: the estimator passed to `surv_fit()`.
#' @param censoring Column name or `NULL`.
#' @param call Enclosing frame for the error signal.
#'
#' @returns Invisibly `NULL`.
#' @family checks
#' @noRd
check_ipcw <- function(
  ipcw,
  estimator,
  censoring,
  call = rlang::caller_env()
) {
  if (is.null(ipcw)) {
    return(invisible(NULL))
  }
  ## Must be a one-sided formula.
  if (!inherits(ipcw, "formula")) {
    rlang::abort(
      c(
        "`ipcw` must be a one-sided formula (e.g. `~ L1 + L2`) or `NULL`.",
        i = "Pass the censoring-model covariate adjustment set as a formula."
      ),
      class = "survatr_bad_ipcw",
      call = call
    )
  }
  ## IPCW is only supported with `estimator = "ipw"`: the censoring weight
  ## multiplies the treatment weight in the weighted MSM. Other estimators
  ## need a different integration path (deferred).
  if (!identical(estimator, "ipw")) {
    rlang::abort(
      c(
        paste0(
          "Built-in IPCW (`ipcw =`) is only supported with ",
          "`estimator = \"ipw\"` in this release."
        ),
        i = paste0(
          "IPCW for `estimator = \"gcomp\"` and Track B (`\"ice\"`) ship ",
          "in later chunks. Got `estimator = \"",
          estimator,
          "\"`."
        )
      ),
      class = "survatr_ipcw_estimator",
      call = call
    )
  }
  ## A censoring column is required: the censoring indicator is the response
  ## for the censoring hazard model.
  if (is.null(censoring)) {
    rlang::abort(
      c(
        paste0(
          "Built-in IPCW (`ipcw =`) requires a `censoring` column (the ",
          "censoring indicator is the response for the censoring hazard model)."
        ),
        i = paste0(
          "Supply `censoring = \"<col>\"` to `surv_fit()`, where `<col>` is ",
          "the 0/1 column flagging censored rows."
        )
      ),
      class = "survatr_ipcw_no_censoring",
      call = call
    )
  }
  invisible(NULL)
}

#' Reject `na.action = na.exclude` in `...`
#'
#' Copy-adapted from `causatr:::check_dots_na_action`. Under `na.exclude`,
#' `residuals(model, "working")` is padded with `NA`s while `model.matrix()`
#' silently drops them; the IF engine's `d_fit * r_score` product then
#' misaligns and silently corrupts the sandwich variance. We only accept
#' `na.omit` (default) or `na.fail`.
#'
#' @param ... Dots captured from the caller (e.g. `surv_fit()`).
#' @param call Enclosing frame for the error signal.
#'
#' @return Invisibly `NULL`. Called for its side effect.
#' @family checks
#' @noRd
check_dots_na_action <- function(..., call = rlang::caller_env()) {
  dots <- list(...)
  if (!"na.action" %in% names(dots)) {
    return(invisible(NULL))
  }
  na_action <- dots$na.action
  ok <- FALSE
  if (is.function(na_action)) {
    ok <- identical(na_action, stats::na.omit) ||
      identical(na_action, stats::na.fail)
  } else if (is.character(na_action) && length(na_action) == 1L) {
    ok <- na_action %in% c("na.omit", "na.fail")
  }
  if (!ok) {
    rlang::abort(
      c(
        "`na.action` must be `na.omit` (default) or `na.fail`.",
        i = paste0(
          "survatr builds its own row-alignment bookkeeping from ",
          "`fit_rows` and the fitted model's `na.action` attribute. ",
          "`na.exclude` pads working residuals with NA and silently ",
          "corrupts the sandwich variance."
        ),
        i = "Drop NA rows before calling `surv_fit()` or use `na.action = na.omit`."
      ),
      class = "survatr_bad_na_action",
      call = call
    )
  }
  invisible(NULL)
}

#' Assert that a column holds only `{0, 1, NA}` values
#'
#' Both the outcome and the censoring column are interpreted as
#' indicators: 0 / NA = no event / uncensored, 1 = event / censored.
#' Other values break two things:
#'
#' - `build_risk_set()` computes `cumsum(outcome)` and `cumsum(censoring)`
#'   lagged within id; non-binary values cause the cumulative sums to
#'   walk integer space, producing a risk-set filter that bears no
#'   relationship to "has the id experienced the event yet?".
#' - `is_uncensored()` treats any value that is not `NA` and not `0` as
#'   "censored". A stray `2` in a `status` column from `survival::Surv`
#'   silently marks that row censored.
#'
#' We reject non-binary indicator columns at the boundary with a classed
#' error pointing to `0/1` recoding.
#'
#' @param data Person-period `data.table`.
#' @param col Column name (character scalar).
#' @param role String describing the column role, inserted into the
#'   error message (e.g. `"outcome"`, `"censoring"`).
#' @param allow_na Logical; when `TRUE`, NA is a permitted value
#'   (censoring allows NA to mean "uncensored"); when `FALSE`, NA is
#'   rejected (outcome cannot be NA — caught separately by
#'   `check_no_na_in_predictors()`).
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @family checks
#' @noRd
check_indicator_col <- function(
  data,
  col,
  role,
  allow_na = FALSE,
  call = rlang::caller_env()
) {
  if (!col %in% names(data)) {
    return(invisible(NULL))
  }
  v <- data[[col]]
  if (is.logical(v)) {
    v <- as.integer(v)
  }
  ok <- v %in% c(0L, 1L)
  if (allow_na) {
    ok <- ok | is.na(v)
  }
  if (!all(ok)) {
    bad_vals <- utils::head(unique(v[!ok]), 5L)
    rlang::abort(
      c(
        paste0(
          "`",
          col,
          "` (",
          role,
          ") must contain only 0 / 1",
          if (allow_na) " / NA" else "",
          " values."
        ),
        i = paste0(
          "Got value(s): ",
          paste(shQuote(bad_vals), collapse = ", "),
          ". Recode to 0 / 1 before calling `surv_fit()`."
        )
      ),
      class = "survatr_bad_indicator",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the competing-risks cause column and return its cause labels
#'
#' The cause-specific hazards + CIF path encodes the event type in a single
#' integer-valued column: `0` = no event at this period (the at-risk-continuing
#' state), `j in {1, ..., J}` = a cause-`j` event at this period. Administrative
#' censoring is carried separately by the `censoring` column, exactly as for the
#' single-event path.
#'
#' Non-integer, negative, or `NA` values break the cause-specific machinery: the
#' shared all-cause risk set is built from `cumsum(competing != 0)` lagged within
#' id, and each cause-`j` model regresses `1{competing == j}`; a fractional or
#' negative code would silently produce a risk set and per-cause indicators that
#' bear no relationship to the intended event structure. We reject those at the
#' boundary rather than fit a meaningless model.
#'
#' @param data Person-period `data.table`.
#' @param competing Column name (character scalar) of the multi-valued cause
#'   column.
#' @param call Caller frame for the error signal.
#'
#' @returns Integer vector of the distinct **positive** cause labels present
#'   (sorted ascending). Called for its return value and its validation side
#'   effect (a classed `survatr_bad_competing` error on invalid input).
#' @family checks
#' @noRd
check_competing_col <- function(data, competing, call = rlang::caller_env()) {
  v <- data[[competing]]
  ## NA is rejected: a missing cause code cannot be assigned to "no event" vs a
  ## specific cause without guessing, and it would poison the lagged cumsum.
  if (anyNA(v)) {
    rlang::abort(
      paste0(
        "`",
        competing,
        "` (competing-risks cause column) contains NA. Recode each row to ",
        "0 (no event this period) or a positive integer cause label."
      ),
      class = "survatr_bad_competing",
      call = call
    )
  }
  ## Must be whole numbers >= 0. `is.integer()` is too strict (users pass
  ## doubles like `0`, `1`, `2`), so test integrality numerically.
  is_whole <- is.numeric(v) && all(abs(v - round(v)) < 1e-8)
  if (!is_whole || any(v < 0)) {
    offending <- if (!is.numeric(v)) v else v[v < 0 | abs(v - round(v)) >= 1e-8]
    bad_vals <- utils::head(unique(offending), 5L)
    rlang::abort(
      c(
        paste0(
          "`",
          competing,
          "` (competing-risks cause column) must contain only non-negative ",
          "integer cause labels (0 = no event, 1..J = cause)."
        ),
        i = paste0(
          "Got value(s): ",
          paste(shQuote(bad_vals), collapse = ", "),
          "."
        )
      ),
      class = "survatr_bad_competing",
      call = call
    )
  }
  ## The distinct positive labels are the competing causes. A column that is all
  ## zeros has no events at all -- there is nothing to model.
  causes <- sort(unique(as.integer(round(v[v > 0]))))
  if (length(causes) == 0L) {
    rlang::abort(
      paste0(
        "`",
        competing,
        "` (competing-risks cause column) has no positive cause labels -- ",
        "every row is 0 (no event). There are no events to model."
      ),
      class = "survatr_bad_competing",
      call = call
    )
  }
  causes
}

#' Reject NA values in predictor columns
#'
#' Track A's contrast path predicts the counterfactual hazard on **every**
#' person-period row (including rows that were censored / had the event
#' in reality — the cumulative product `prod_{k <= t} (1 - h^a_{i,k})`
#' needs every per-period hazard). Any NA in a predictor column would
#' either (a) silently produce NA predictions that propagate through
#' `cumprod` and corrupt `s_hat`, or (b) be dropped by `glm`'s
#' `na.action = na.omit` (the default we permit), creating a size
#' mismatch between `prep$X_fit` / `r_score` and the post-drop fit
#' rows — which crashes the sandwich IF chain downstream.
#'
#' The simplest and least surprising contract is to reject NAs in the
#' required columns at the boundary. Users should preprocess (impute or
#' drop) before calling `surv_fit()`. NA in the `censoring` column
#' retains its "uncensored" semantics (see `is_uncensored()`) and is
#' **not** rejected here.
#'
#' @param data Person-period `data.table`.
#' @param cols Character vector of predictor column names to check.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @family checks
#' @noRd
check_no_na_in_predictors <- function(data, cols, call = rlang::caller_env()) {
  cols <- intersect(cols, names(data))
  if (length(cols) == 0L) {
    return(invisible(NULL))
  }
  na_by_col <- vapply(
    cols,
    function(cn) anyNA(data[[cn]]),
    logical(1L)
  )
  if (any(na_by_col)) {
    bad <- cols[na_by_col]
    rlang::abort(
      c(
        paste0(
          "Predictor column(s) ",
          paste0("`", bad, "`", collapse = ", "),
          " contain NA values. Track A's contrast path predicts the ",
          "counterfactual hazard on every person-period row and cannot ",
          "be NA-safe."
        ),
        i = paste0(
          "Drop or impute NA rows before calling `surv_fit()`. NA in ",
          "the `censoring` column is separately treated as 'uncensored' ",
          "and is not rejected."
        )
      ),
      class = "survatr_na_in_predictors",
      call = call
    )
  }
  invisible(NULL)
}

#' Reject user-data collisions with survatr's reserved column names
#'
#' survatr adds `.survatr_prev_event` and `.survatr_prev_cens` to the
#' person-period data during risk-set construction. If the user's data
#' already contains a column with one of these names we would silently
#' overwrite it, so abort upfront.
#'
#' @param data Person-period data (anything with `names()`).
#' @param which Character vector of reserved names to check against.
#' @param call Enclosing frame for the error signal.
#'
#' @return Invisibly `NULL`. Called for its side effect.
#' @family checks
#' @noRd
check_reserved_cols <- function(
  data,
  which = SURVATR_RESERVED_COLS,
  call = rlang::caller_env()
) {
  bad <- intersect(which, names(data))
  if (length(bad) > 0L) {
    rlang::abort(
      paste0(
        "Column name(s) ",
        paste0("`", bad, "`", collapse = ", "),
        " are reserved by survatr internals. Rename the column(s) in ",
        "your input data before calling `surv_fit()`."
      ),
      class = "survatr_reserved_col",
      call = call
    )
  }
  invisible(NULL)
}
