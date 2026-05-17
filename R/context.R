#' @title Internal Context Object for GWPR.light 1.0.0
#'
#' @description
#' These internal functions construct and validate the standardised
#' `gwpr_context` list that is passed between modules, eliminating
#' repetitive argument passing.
#'
#' @name context
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# new_gwpr_context
# ---------------------------------------------------------------------------

#' Construct a new gwpr_context object
#'
#' Creates the standardised internal context list used to pass state between
#' GWPR modules.  Any field not supplied defaults to `NULL` (or, for
#' `metadata` and `warnings`, to their appropriate empty types).
#'
#' @param call       The matched call from the top-level API function.
#' @param formula    A formula object.
#' @param family     Character: `"gaussian"` or `"binomial"`.
#' @param model      Character: `"pooling"`, `"within"`, or `"random"`.
#' @param effect     Character: `"individual"`, `"time"`, `"two-way"`, or
#'   `"nested"`.
#' @param id         Name of the unit ID column.
#' @param time       Name of the time column.
#' @param kernel     Kernel name.
#' @param adaptive   Logical; `TRUE` for adaptive bandwidth.
#' @param threshold  Numeric classification threshold (Logistic).
#' @param workers    Number of parallel workers.
#' @param seed       Integer random seed or `NULL`.
#' @param raw_data   The original user-supplied data frame.
#' @param raw_spatial The original user-supplied sf object.
#' @param panel_data Processed panel data frame.
#' @param spatial_data Processed sf object.
#' @param id_map     Named integer vector mapping unit IDs to row indices.
#' @param coords     Matrix of spatial coordinates.
#' @param model_frame Model frame derived from formula and panel_data.
#' @param model_matrix Design matrix.
#' @param response   Numeric response vector.
#' @param metadata   Named list of supplementary information.
#' @param warnings   Character vector of accumulated warnings.
#' @param ...        Additional named fields stored in the context list.
#'
#' @return A named list with class `"gwpr_context"`.
#' @keywords internal
new_gwpr_context <- function(
    call          = NULL,
    formula       = NULL,
    family        = NULL,
    model         = NULL,
    effect        = NULL,
    id            = NULL,
    time          = NULL,
    kernel        = NULL,
    adaptive      = NULL,
    threshold     = NULL,
    workers       = NULL,
    seed          = NULL,
    raw_data      = NULL,
    raw_spatial   = NULL,
    panel_data    = NULL,
    spatial_data  = NULL,
    id_map        = NULL,
    coords        = NULL,
    model_frame   = NULL,
    model_matrix  = NULL,
    response      = NULL,
    metadata      = list(),
    warnings      = character(),
    ...
) {
  extra <- list(...)

  ctx <- c(
    list(
      call         = call,
      formula      = formula,
      family       = family,
      model        = model,
      effect       = effect,
      id           = id,
      time         = time,
      kernel       = kernel,
      adaptive     = adaptive,
      threshold    = threshold,
      workers      = workers,
      seed         = seed,
      raw_data     = raw_data,
      raw_spatial  = raw_spatial,
      panel_data   = panel_data,
      spatial_data = spatial_data,
      id_map       = id_map,
      coords       = coords,
      model_frame  = model_frame,
      model_matrix = model_matrix,
      response     = response,
      metadata     = metadata,
      warnings     = warnings
    ),
    extra
  )

  structure(ctx, class = "gwpr_context")
}

# ---------------------------------------------------------------------------
# validate_gwpr_context
# ---------------------------------------------------------------------------

#' Validate a gwpr_context object
#'
#' Checks that all core fields required for model fitting are non-`NULL`.
#' Stops with an informative message listing every missing field.
#'
#' Core fields: `formula`, `family`, `id`, `time`, `model`, `effect`,
#' `kernel`, `adaptive`, `threshold`, `workers`.
#'
#' @param context A list (typically of class `"gwpr_context"`) to validate.
#'
#' @return Invisibly returns `TRUE` when all core fields are present.
#' @keywords internal
validate_gwpr_context <- function(context) {
  if (!is.list(context)) {
    stop("`context` must be a list.", call. = FALSE)
  }

  core_fields <- c(
    "formula", "family", "id", "time",
    "model", "effect", "kernel", "adaptive",
    "threshold", "workers"
  )

  missing_fields <- core_fields[vapply(
    core_fields,
    function(f) is.null(context[[f]]),
    logical(1L)
  )]

  if (length(missing_fields) > 0L) {
    stop(
      "gwpr_context is missing required core field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
