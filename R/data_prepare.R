#' @title Data Preparation Module for GWPR.light 1.0.0
#'
#' @description
#' Internal functions that convert user inputs into the internal data structures
#' required by the model engine: panel indices, spatial alignment, model frame,
#' and model matrix.
#'
#' @name data_prepare
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# prepare_data  (top-level entry point)
# ---------------------------------------------------------------------------

#' Prepare all data structures for GWPR fitting
#'
#' Orchestrates the full data preparation pipeline: panel data extraction,
#' spatial data extraction, ID mapping, model frame and model matrix
#' construction. Updates and returns a `gwpr_context`.
#'
#' @param context A `gwpr_context` object with at least `formula`, `family`,
#'   `id`, `time`, `model`, `raw_data`, and `raw_spatial` populated.
#'
#' @return An updated `gwpr_context` with `panel_data`, `spatial_data`,
#'   `id_map`, `coords`, `model_frame`, `model_matrix`, `response`, and
#'   `metadata` filled in.
#' @keywords internal
prepare_data <- function(context) {
  context <- prepare_panel_data(context)
  context <- prepare_spatial_data(context)

  context$id_map <- build_id_map(
    panel_data   = context$panel_data,
    spatial_data = context$spatial_data,
    id           = context$id
  )

  context <- build_model_frame(context)
  context <- build_model_matrix(context)

  context
}

# ---------------------------------------------------------------------------
# prepare_panel_data
# ---------------------------------------------------------------------------

#' Prepare panel data from the raw data in context
#'
#' Extracts the `id`, `time`, and formula variables from `raw_data`. Adds a
#' `raw_row_id` column preserving the original row positions. Sorts by id then
#' time. Records panel balance information in `metadata`. For `within` models,
#' records single-observation individuals in `metadata`.
#'
#' @param context A `gwpr_context` with `raw_data`, `formula`, `id`, `time`,
#'   and `model` populated.
#'
#' @return Updated context with `panel_data` and `metadata` populated.
#' @keywords internal
prepare_panel_data <- function(context) {
  data    <- context$raw_data
  formula <- context$formula
  id      <- context$id
  time    <- context$time
  model   <- context$model

  if (!is.data.frame(data)) {
    stop("`raw_data` must be a data.frame.", call. = FALSE)
  }

  # Identify all variables referenced in the formula
  formula_vars <- all.vars(formula)
  # Keep only those that exist in data (response + predictors + id + time)
  keep_vars <- unique(c(id, time, formula_vars))
  keep_vars <- intersect(keep_vars, names(data))

  panel_data <- data[, keep_vars, drop = FALSE]

  # Preserve original row positions
  panel_data$raw_row_id <- seq_len(nrow(data))

  # Sort by id then time for predictable ordering
  panel_data <- panel_data[order(panel_data[[id]], panel_data[[time]]), , drop = FALSE]
  rownames(panel_data) <- NULL

  # Determine panel balance
  obs_per_unit <- table(panel_data[[id]])
  n_obs_unique <- length(unique(as.integer(obs_per_unit)))
  balanced     <- n_obs_unique == 1L

  metadata <- context$metadata
  metadata$panel_balanced <- balanced
  metadata$n_units        <- length(obs_per_unit)
  metadata$obs_per_unit   <- as.integer(obs_per_unit)
  names(metadata$obs_per_unit) <- names(obs_per_unit)

  # Record single-observation individuals when using within model
  if (!is.null(model) && model == "within") {
    single_obs_ids <- names(obs_per_unit)[obs_per_unit == 1L]
    if (length(single_obs_ids) > 0L) {
      metadata$within_single_obs_ids <- single_obs_ids
      warning(
        "Within model: ", length(single_obs_ids),
        " individual(s) have only one observation and cannot contribute to ",
        "within-individual variation. Recorded in metadata$within_single_obs_ids.",
        call. = FALSE
      )
    }
  }

  context$panel_data <- panel_data
  context$metadata   <- metadata
  context
}

# ---------------------------------------------------------------------------
# prepare_spatial_data
# ---------------------------------------------------------------------------

#' Extract and align spatial data for GWPR fitting
#'
#' Extracts geometry and representative coordinates from the `raw_spatial` sf
#' object. Aligns spatial rows to the unique individual IDs present in
#' `panel_data` via `id_map`.
#'
#' @param context A `gwpr_context` with `raw_spatial`, `panel_data`, and `id`
#'   populated.
#'
#' @return Updated context with `spatial_data` and `coords` populated.
#' @keywords internal
prepare_spatial_data <- function(context) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Please install it.", call. = FALSE)
  }

  spatial    <- context$raw_spatial
  panel_data <- context$panel_data
  id         <- context$id

  if (!inherits(spatial, "sf")) {
    stop("`raw_spatial` must be an sf object.", call. = FALSE)
  }

  # Unique panel IDs (preserve order after sorting by id, then time)
  panel_ids <- unique(panel_data[[id]])

  # Filter spatial to only those IDs present in panel data
  spatial_ids     <- spatial[[id]]
  in_panel        <- spatial_ids %in% panel_ids
  spatial_aligned <- spatial[in_panel, , drop = FALSE]

  # Re-order to match the order of unique panel IDs
  row_order <- match(panel_ids, spatial_aligned[[id]])

  # Check all panel IDs are found in spatial
  missing_mask <- is.na(row_order)
  if (any(missing_mask)) {
    missing_ids <- panel_ids[missing_mask]
    stop(
      "The following panel ID(s) have no matching row in `spatial`: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  spatial_data <- spatial_aligned[row_order, , drop = FALSE]
  rownames(spatial_data) <- NULL

  # Extract representative coordinates
  coords <- extract_coords_from_sf(spatial_data)

  context$spatial_data <- spatial_data
  context$coords       <- coords
  context
}

# ---------------------------------------------------------------------------
# extract_coords_from_sf  (internal helper)
# ---------------------------------------------------------------------------

#' Extract representative XY coordinates from an sf object
#'
#' Uses point coordinates for POINT geometries and centroids for other
#' geometry types (POLYGON, MULTIPOLYGON, etc.).
#'
#' @param spatial An `sf` object.
#'
#' @return A numeric matrix with columns `X` and `Y`.
#' @keywords internal
extract_coords_from_sf <- function(spatial) {
  geom_type <- unique(as.character(sf::st_geometry_type(spatial)))

  if (length(geom_type) == 1L && geom_type == "POINT") {
    coords_sf <- sf::st_coordinates(spatial)
    coords    <- coords_sf[, c("X", "Y"), drop = FALSE]
  } else {
    # Use centroids for polygon / other types; suppress warnings about
    # geographic CRS centroids
    suppressWarnings(
      centroids <- sf::st_centroid(sf::st_geometry(spatial))
    )
    coords_sf <- sf::st_coordinates(centroids)
    coords    <- coords_sf[, c("X", "Y"), drop = FALSE]
  }

  coords
}

# ---------------------------------------------------------------------------
# build_id_map
# ---------------------------------------------------------------------------

#' Build a mapping from panel unit IDs to spatial row indices
#'
#' Returns a named integer vector where each name is a panel unit ID (as a
#' character string) and each value is the corresponding row index in
#' `spatial_data` (1-based).
#'
#' Rules:
#' - Every panel ID must have a spatial match; missing IDs cause an error.
#' - Extra spatial rows (not in panel) are silently ignored.
#'
#' @param panel_data   A data frame with a column named `id`.
#' @param spatial_data An `sf` data frame with a column named `id`.
#' @param id           Character; name of the shared ID column.
#'
#' @return Named integer vector mapping unit ID to spatial row index.
#' @keywords internal
build_id_map <- function(panel_data, spatial_data, id) {
  panel_ids   <- unique(as.character(panel_data[[id]]))
  spatial_ids <- as.character(spatial_data[[id]])

  row_indices <- match(panel_ids, spatial_ids)

  missing_mask <- is.na(row_indices)
  if (any(missing_mask)) {
    missing_ids <- panel_ids[missing_mask]
    stop(
      "The following panel ID(s) have no matching row in `spatial_data`: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  id_map        <- row_indices
  names(id_map) <- panel_ids
  id_map
}

# ---------------------------------------------------------------------------
# build_model_frame
# ---------------------------------------------------------------------------

#' Build the model frame from the formula and panel data
#'
#' @param context A `gwpr_context` with `formula` and `panel_data` populated.
#'
#' @return Updated context with `model_frame` populated.
#' @keywords internal
build_model_frame <- function(context) {
  formula    <- context$formula
  panel_data <- context$panel_data

  mf <- stats::model.frame(formula, data = panel_data, na.action = stats::na.pass)

  context$model_frame <- mf
  context
}

# ---------------------------------------------------------------------------
# build_model_matrix
# ---------------------------------------------------------------------------

#' Build the design matrix and response vector
#'
#' Extracts the response variable `y` and the design matrix `X` from the
#' model frame.  For `binomial` family, the response is standardised to 0/1
#' via `standardize_binary_response()`.
#'
#' @param context A `gwpr_context` with `model_frame`, `formula`, and `family`
#'   populated.
#'
#' @return Updated context with `model_matrix` and `response` populated.
#' @keywords internal
build_model_matrix <- function(context) {
  mf     <- context$model_frame
  family <- context$family

  # Extract response
  y_raw <- stats::model.response(mf)

  if (!is.null(family) && family == "binomial") {
    y <- standardize_binary_response(y_raw)
  } else {
    y <- as.numeric(y_raw)
  }

  # Build design matrix (removes the response column)
  X <- stats::model.matrix(context$formula, data = mf)

  context$model_matrix <- X
  context$response     <- y
  context
}

# ---------------------------------------------------------------------------
# standardize_binary_response
# ---------------------------------------------------------------------------

#' Standardise a binary response variable to 0/1 numeric
#'
#' Converts logical and two-level factor responses to 0/1 integer. Numeric
#' 0/1 vectors are returned unchanged (as numeric). Other inputs raise an
#' error.
#'
#' For factor inputs the first level is mapped to 0 and the second level is
#' mapped to 1.
#'
#' @param y A vector that is 0/1 numeric, logical, or a two-level factor.
#'
#' @return A numeric vector of 0s and 1s.
#' @keywords internal
standardize_binary_response <- function(y) {
  if (is.logical(y)) {
    return(as.integer(y))
  }

  if (is.factor(y)) {
    lvls <- levels(y)
    if (length(lvls) != 2L) {
      stop(
        "`y` must be a two-level factor for binary response standardisation. ",
        "Found ", length(lvls), " level(s).",
        call. = FALSE
      )
    }
    # First level -> 0, second level -> 1
    return(as.integer(y) - 1L)
  }

  if (is.numeric(y)) {
    unique_vals <- unique(y[!is.na(y)])
    invalid     <- setdiff(unique_vals, c(0, 1))
    if (length(invalid) > 0L) {
      stop(
        "Numeric `y` must contain only 0 and 1 for binary response. ",
        "Found: ", paste(invalid, collapse = ", "),
        call. = FALSE
      )
    }
    return(as.numeric(y))
  }

  stop(
    "`y` must be 0/1 numeric, logical, or a two-level factor.",
    call. = FALSE
  )
}
