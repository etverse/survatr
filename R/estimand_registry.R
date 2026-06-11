#' Estimand descriptor registry
#'
#' @description
#' Single source of truth for the per-estimand properties that the contrast,
#' variance, bootstrap, and S3 paths all dispatch on. Adding an estimand is a
#' new entry here plus -- only when the math is genuinely new -- a new SE family
#' in the variance code. This replaces the scattered `switch()` /
#' `type %in% c(...)` dispatch that previously had to be kept in sync by hand
#' across ~10 sites (`build_contrasts()`, `fill_sandwich_ses()`,
#' `boot_estimand_col()`, the bootstrap `has_contrast` test, `plot()`, `tidy()`,
#' `validate_reference()`, and the `contrast()` type allow-lists).
#'
#' @details
#' Each entry carries:
#' - `point_col` -- the estimand column in `estimates` (also the value S3 /
#'   bootstrap read).
#' - `kind` -- `"curve"` (per-intervention only) or `"contrast"` (also builds
#'   pairwise contrasts).
#' - `op` -- the contrast operator, `"difference"` / `"ratio"` (or `NA` for
#'   curve types).
#' - `index` -- the minor result index, `"time"` for the curve estimands or
#'   `"q"` for the quantile (which collapses the time axis).
#' - `level_se` -- the per-intervention SE family (`"pointwise"` / `"rmst"` /
#'   `"quantile"` / `"cif"`).
#' - `contrast_se` -- the pairwise-contrast SE family (`NA` for curve types).
#' - `ylab` -- default y-axis label for `plot()`.
#' - `dim_cause` -- `TRUE` when the estimand carries a competing-risks `cause`
#'   dimension.
#' - `single` / `competing` -- which fit kind(s) the estimand is valid for.
#'
#' The competing-risks estimands additionally route their SE through the
#' dedicated `fill_sandwich_ses_cr()` (per-cause, all-cause `(1 - H)`
#' sensitivity); the `level_se` / `contrast_se` fields are recorded here for
#' completeness but the single-event `fill_sandwich_ses()` never sees them.
#'
#' @format A named list of estimand descriptor lists.
#' @noRd
.SURVATR_ESTIMANDS <- list(
  survival = list(
    point_col = "s_hat",
    kind = "curve",
    op = NA_character_,
    index = "time",
    level_se = "pointwise",
    contrast_se = NA_character_,
    ylab = "S(t)",
    dim_cause = FALSE,
    single = TRUE,
    competing = TRUE
  ),
  risk = list(
    point_col = "risk_hat",
    kind = "curve",
    op = NA_character_,
    index = "time",
    level_se = "pointwise",
    contrast_se = NA_character_,
    ylab = "1 - S(t)",
    dim_cause = FALSE,
    single = TRUE,
    competing = TRUE
  ),
  risk_difference = list(
    point_col = "risk_hat",
    kind = "contrast",
    op = "difference",
    index = "time",
    level_se = "pointwise",
    contrast_se = "difference",
    ylab = "risk difference",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  risk_ratio = list(
    point_col = "risk_hat",
    kind = "contrast",
    op = "ratio",
    index = "time",
    level_se = "pointwise",
    contrast_se = "logratio",
    ylab = "risk ratio",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  rmst = list(
    point_col = "rmst_hat",
    kind = "curve",
    op = NA_character_,
    index = "time",
    level_se = "rmst",
    contrast_se = NA_character_,
    ylab = "RMST(t)",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  rmst_difference = list(
    point_col = "rmst_hat",
    kind = "contrast",
    op = "difference",
    index = "time",
    level_se = "rmst",
    contrast_se = "rmst",
    ylab = "RMST difference",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  rmtl = list(
    point_col = "rmtl_hat",
    kind = "curve",
    op = NA_character_,
    index = "time",
    level_se = "rmst",
    contrast_se = NA_character_,
    ylab = "RMTL(t)",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  rmtl_difference = list(
    point_col = "rmtl_hat",
    kind = "contrast",
    op = "difference",
    index = "time",
    level_se = "rmst",
    contrast_se = "rmst",
    ylab = "RMTL difference",
    dim_cause = FALSE,
    single = TRUE,
    competing = FALSE
  ),
  cif = list(
    point_col = "cif_hat",
    kind = "curve",
    op = NA_character_,
    index = "time",
    level_se = "cif",
    contrast_se = NA_character_,
    ylab = "F(t)",
    dim_cause = TRUE,
    single = FALSE,
    competing = TRUE
  ),
  cif_difference = list(
    point_col = "cif_hat",
    kind = "contrast",
    op = "difference",
    index = "time",
    level_se = "cif",
    contrast_se = "cif",
    ylab = "CIF difference",
    dim_cause = TRUE,
    single = FALSE,
    competing = TRUE
  ),
  cif_ratio = list(
    point_col = "cif_hat",
    kind = "contrast",
    op = "ratio",
    index = "time",
    level_se = "cif",
    contrast_se = "cif_logratio",
    ylab = "CIF ratio",
    dim_cause = TRUE,
    single = FALSE,
    competing = TRUE
  )
)

#' Look up a field of an estimand descriptor
#'
#' @param type Estimand type string (a key of `.SURVATR_ESTIMANDS`).
#' @param field Descriptor field name (e.g. `"point_col"`, `"op"`, `"index"`,
#'   `"level_se"`, `"contrast_se"`, `"ylab"`).
#'
#' @returns The field value. Errors loudly (via `[[`) if `type` is unknown,
#'   which only happens on an internal wiring bug since `contrast()` validates
#'   `type` against `estimand_types()` upfront.
#' @noRd
estimand_field <- function(type, field) {
  .SURVATR_ESTIMANDS[[type]][[field]]
}

#' Does this estimand build pairwise contrasts?
#'
#' @param type Estimand type string.
#' @returns `TRUE` for the `"contrast"` kinds (`*_difference` / `*_ratio` and
#'   the quantile), `FALSE` for the curve-only kinds.
#' @noRd
estimand_is_contrast <- function(type) {
  identical(.SURVATR_ESTIMANDS[[type]]$kind, "contrast")
}

#' Is this a curve-only estimand (no pairwise contrasts)?
#'
#' @param type Estimand type string.
#' @returns `TRUE` for `survival` / `risk` / `rmst` / `rmtl` / `cif`.
#' @noRd
estimand_is_curve <- function(type) {
  identical(.SURVATR_ESTIMANDS[[type]]$kind, "curve")
}

#' Does this estimand require at least two interventions?
#'
#' @description
#' The strict pairwise contrasts (`*_difference` / `*_ratio`) are meaningless
#' with a single intervention and abort upfront. The quantile is a contrast kind
#' too (it can report a median difference) but tolerates a single arm -- the
#' lone median is a valid request -- so it is excluded. Derived from existing
#' fields (a contrast that is not `q`-indexed) rather than a stored flag.
#'
#' @param type Estimand type string.
#' @returns `TRUE` for the strict pairwise contrasts, `FALSE` otherwise.
#' @noRd
estimand_requires_pair <- function(type) {
  estimand_is_contrast(type) &&
    !identical(.SURVATR_ESTIMANDS[[type]]$index, "q")
}

#' Does this estimand carry a competing-risks `cause` dimension?
#'
#' @param type Estimand type string.
#' @returns `TRUE` for the `cif*` estimands.
#' @noRd
estimand_has_cause <- function(type) {
  isTRUE(.SURVATR_ESTIMANDS[[type]]$dim_cause)
}

#' All registered estimand types
#'
#' @returns Character vector of every estimand key (the `contrast()` `type`
#'   allow-list).
#' @noRd
estimand_types <- function() {
  names(.SURVATR_ESTIMANDS)
}

#' Estimand types valid on a single-event fit
#'
#' @returns Character vector of the single-event estimand keys.
#' @noRd
single_event_estimands <- function() {
  names(Filter(function(s) isTRUE(s$single), .SURVATR_ESTIMANDS))
}

#' Estimand types valid on a competing-risks fit
#'
#' @returns Character vector of the competing-risks estimand keys.
#' @noRd
competing_estimands <- function() {
  names(Filter(function(s) isTRUE(s$competing), .SURVATR_ESTIMANDS))
}
