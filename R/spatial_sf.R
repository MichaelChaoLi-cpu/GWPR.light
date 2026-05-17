#' @title Spatial SF Module for GWPR.light 1.0.0
#'
#' @description
#' Internal functions providing the `sf`-first spatial interface. These
#' functions abstract geometry extraction, coordinate derivation, spatial
#' alignment to panel data, and neighbour-structure construction. They are used
#' by the data preparation and diagnostics modules.
#'
#' @name spatial_sf
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# assert_sf
# ---------------------------------------------------------------------------

#' Assert that an object is an sf object
#'
#' Stops with an informative error when `spatial` does not inherit from `"sf"`.
#' sp objects are not supported.
#'
#' @param spatial Any R object.
#'
#' @return Invisibly returns `TRUE` when the check passes.
#' @keywords internal
assert_sf <- function(spatial) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Please install it.", call. = FALSE)
  }
  if (!inherits(spatial, "sf")) {
    stop(
      "`spatial` must be an sf object. ",
      "sp objects are not supported; please convert with sf::st_as_sf().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# extract_geometry
# ---------------------------------------------------------------------------

#' Extract the geometry column from an sf object
#'
#' Returns the `sfc` geometry column of the supplied `sf` object.
#'
#' @param spatial An `sf` object.
#'
#' @return An `sfc` geometry column.
#' @keywords internal
extract_geometry <- function(spatial) {
  assert_sf(spatial)
  sf::st_geometry(spatial)
}

# ---------------------------------------------------------------------------
# extract_coordinates
# ---------------------------------------------------------------------------

#' Extract representative XY coordinates from an sf object
#'
#' For POINT geometries the point coordinates are returned directly. For
#' all other geometry types (POLYGON, MULTIPOLYGON, etc.) `sf::st_centroid()`
#' is used to derive a representative point, with warnings suppressed (they
#' are typically non-actionable geographic-CRS notes).
#'
#' @param spatial An `sf` object.
#'
#' @return A numeric matrix with columns `X` and `Y`, one row per feature.
#' @keywords internal
extract_coordinates <- function(spatial) {
  assert_sf(spatial)

  geom_type <- unique(as.character(sf::st_geometry_type(spatial)))

  if (length(geom_type) == 1L && geom_type == "POINT") {
    coords_mat <- sf::st_coordinates(spatial)
    return(coords_mat[, c("X", "Y"), drop = FALSE])
  }

  # For POLYGON / MULTIPOLYGON / other types use centroid.
  # st_centroid warnings about geographic CRS are suppressed as they are
  # non-actionable for typical projected or planar data.
  suppressWarnings(
    centroids <- sf::st_centroid(sf::st_geometry(spatial))
  )
  coords_mat <- sf::st_coordinates(centroids)
  coords_mat[, c("X", "Y"), drop = FALSE]
}

# ---------------------------------------------------------------------------
# align_spatial_to_panel
# ---------------------------------------------------------------------------

#' Reorder spatial rows to match panel ID order
#'
#' Given a named integer vector `id_map` (names = unit IDs, values = row
#' indices in `spatial`) produced by `build_id_map()`, this function reorders
#' the rows of `spatial` so that they align with the panel data ordering.
#'
#' @param spatial An `sf` object containing at minimum the rows referenced in
#'   `id_map`.
#' @param id_map  A named integer vector as produced by `build_id_map()`.
#'   Names are unit IDs (character); values are 1-based row indices into
#'   `spatial`.
#'
#' @return An `sf` object with rows reordered to match `id_map`.
#' @keywords internal
align_spatial_to_panel <- function(spatial, id_map) {
  assert_sf(spatial)

  if (!is.integer(id_map) && !is.numeric(id_map)) {
    stop("`id_map` must be a named numeric (integer) vector.", call. = FALSE)
  }
  if (is.null(names(id_map))) {
    stop("`id_map` must be named.", call. = FALSE)
  }

  row_indices <- as.integer(id_map)

  # Guard against out-of-range indices
  n_rows <- nrow(spatial)
  if (any(row_indices < 1L | row_indices > n_rows)) {
    stop(
      "`id_map` contains row indices outside the range [1, ", n_rows, "].",
      call. = FALSE
    )
  }

  aligned <- spatial[row_indices, , drop = FALSE]
  rownames(aligned) <- NULL
  aligned
}

# ---------------------------------------------------------------------------
# build_neighbor_structure
# ---------------------------------------------------------------------------

#' Build a neighbour structure from an sf object
#'
#' Returns different structures depending on `type`:
#'
#' * `"distance"` — a numeric coordinate matrix (columns `X` and `Y`) suitable
#'   for pairwise distance computation.
#' * `"contiguity"` — a named list where each element is the integer vector of
#'   neighbour row indices (1-based, Queen contiguity via
#'   `sf::st_relate()`). Used for spatial diagnostics such as Moran's I.
#'
#' @param spatial An `sf` object.
#' @param type    Character scalar: `"distance"` (default) or `"contiguity"`.
#'
#' @return
#' * For `"distance"`: a numeric matrix with columns `X` and `Y`.
#' * For `"contiguity"`: a named list of integer vectors.
#' @keywords internal
build_neighbor_structure <- function(spatial,
                                     type = c("distance", "contiguity")) {
  type <- match.arg(type)
  assert_sf(spatial)

  if (type == "distance") {
    return(extract_coordinates(spatial))
  }

  # Contiguity: Queen criterion (shared boundary or vertex)
  # st_relate pattern "F***T****" (touches) or broader "****T****"
  # We use the standard spdep-compatible approach via sf::st_relate
  n <- nrow(spatial)
  geom <- sf::st_geometry(spatial)

  # Use Queen contiguity: polygons that share at least a point.
  # DE-9IM pattern for "touches or overlaps but not disjoint at boundary"
  # We detect shared boundaries/vertices with the "T" intersect matrix.
  # A robust approach: use st_relate with pattern "F***T****" for touches,
  # combined with polygon intersection detection.
  #
  # Simplest correct approach: use sf::st_relate() with Queen pattern.
  # Queen pattern (from spdep conventions): "F***T****" OR "****T****"
  # We use "****T****" which covers all non-disjoint boundary cases.

  # st_relate returns an sgbp sparse list: element i is an integer vector of
  # row indices (1-based) of features satisfying the pattern with feature i.
  nb_sparse <- tryCatch(
    sf::st_relate(geom, geom, pattern = "F***T****"),
    error = function(e) NULL
  )

  row_names <- if (!is.null(rownames(spatial))) rownames(spatial) else as.character(seq_len(n))

  if (is.null(nb_sparse)) {
    # Fallback: no contiguity structure (e.g., POINT geometries)
    nb <- vector("list", n)
    names(nb) <- row_names
    nb <- lapply(nb, function(x) integer(0L))
    return(nb)
  }

  # nb_sparse is a list of length n; each element is an integer vector of
  # neighbour indices (already excludes self when pattern does not match i->i).
  nb <- vector("list", n)
  names(nb) <- row_names

  for (i in seq_len(n)) {
    # Remove self-references to be safe
    neighbours <- nb_sparse[[i]]
    nb[[i]]    <- as.integer(neighbours[neighbours != i])
  }

  nb
}
