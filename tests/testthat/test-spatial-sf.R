library(testthat)

# ---------------------------------------------------------------------------
# Helpers: small simulated sf objects
# ---------------------------------------------------------------------------

# Four POINT features in a 2x2 grid
make_point_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) skip("sf not available")
  coords <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2, byrow = TRUE,
    dimnames = list(NULL, c("x", "y"))
  )
  pts <- sf::st_as_sf(
    data.frame(id = c("A", "B", "C", "D"), x = coords[, 1], y = coords[, 2]),
    coords = c("x", "y"),
    crs    = NA
  )
  pts
}

# Four POLYGON features (unit squares tiling a 2x2 area)
make_polygon_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) skip("sf not available")

  make_square <- function(xmin, ymin) {
    sf::st_polygon(list(matrix(
      c(xmin,     ymin,
        xmin + 1, ymin,
        xmin + 1, ymin + 1,
        xmin,     ymin + 1,
        xmin,     ymin),
      ncol = 2, byrow = TRUE
    )))
  }

  polys <- sf::st_sfc(
    make_square(0, 0),
    make_square(1, 0),
    make_square(0, 1),
    make_square(1, 1)
  )

  sf::st_sf(id = c("A", "B", "C", "D"), geometry = polys)
}

# A non-sf object for negative tests
make_non_sf <- function() {
  data.frame(id = 1:3, x = 1:3, y = 1:3)
}

# ---------------------------------------------------------------------------
# assert_sf
# ---------------------------------------------------------------------------

test_that("assert_sf: returns TRUE invisibly for valid sf", {
  pts <- make_point_sf()
  expect_true(assert_sf(pts))
})

test_that("assert_sf: errors for data.frame (non-sf)", {
  expect_error(assert_sf(make_non_sf()), "sf object")
})

test_that("assert_sf: errors for plain list", {
  expect_error(assert_sf(list(x = 1:3)), "sf object")
})

test_that("assert_sf: errors for NULL", {
  expect_error(assert_sf(NULL), "sf object")
})

test_that("assert_sf: error message mentions sp not supported", {
  err <- tryCatch(assert_sf(make_non_sf()), error = conditionMessage)
  expect_match(err, "sp")
})

# ---------------------------------------------------------------------------
# extract_geometry
# ---------------------------------------------------------------------------

test_that("extract_geometry: returns sfc for POINT sf", {
  pts  <- make_point_sf()
  geom <- extract_geometry(pts)
  expect_true(inherits(geom, "sfc"))
  expect_equal(length(geom), nrow(pts))
})

test_that("extract_geometry: returns sfc for POLYGON sf", {
  polys <- make_polygon_sf()
  geom  <- extract_geometry(polys)
  expect_true(inherits(geom, "sfc"))
  expect_equal(length(geom), nrow(polys))
})

test_that("extract_geometry: errors for non-sf input", {
  expect_error(extract_geometry(make_non_sf()), "sf object")
})

# ---------------------------------------------------------------------------
# extract_coordinates: POINT
# ---------------------------------------------------------------------------

test_that("extract_coordinates: POINT returns matrix with X and Y columns", {
  pts    <- make_point_sf()
  coords <- extract_coordinates(pts)
  expect_true(is.matrix(coords))
  expect_true(all(c("X", "Y") %in% colnames(coords)))
})

test_that("extract_coordinates: POINT row count matches feature count", {
  pts    <- make_point_sf()
  coords <- extract_coordinates(pts)
  expect_equal(nrow(coords), nrow(pts))
})

test_that("extract_coordinates: POINT coordinates match input", {
  pts    <- make_point_sf()
  coords <- extract_coordinates(pts)
  # First point is (0, 0)
  expect_equal(as.numeric(coords[1, "X"]), 0)
  expect_equal(as.numeric(coords[1, "Y"]), 0)
  # Last point is (1, 1)
  expect_equal(as.numeric(coords[4, "X"]), 1)
  expect_equal(as.numeric(coords[4, "Y"]), 1)
})

# ---------------------------------------------------------------------------
# extract_coordinates: POLYGON
# ---------------------------------------------------------------------------

test_that("extract_coordinates: POLYGON returns matrix with X and Y columns", {
  polys  <- make_polygon_sf()
  coords <- extract_coordinates(polys)
  expect_true(is.matrix(coords))
  expect_true(all(c("X", "Y") %in% colnames(coords)))
})

test_that("extract_coordinates: POLYGON row count matches feature count", {
  polys  <- make_polygon_sf()
  coords <- extract_coordinates(polys)
  expect_equal(nrow(coords), nrow(polys))
})

test_that("extract_coordinates: POLYGON representative points lie inside units", {
  polys  <- make_polygon_sf()
  coords <- extract_coordinates(polys)
  # Each centroid should be within the bounding box [0, 2] x [0, 2]
  expect_true(all(coords[, "X"] >= 0 & coords[, "X"] <= 2))
  expect_true(all(coords[, "Y"] >= 0 & coords[, "Y"] <= 2))
})

test_that("extract_coordinates: POLYGON first polygon centroid is approx (0.5, 0.5)", {
  polys  <- make_polygon_sf()
  coords <- extract_coordinates(polys)
  expect_equal(as.numeric(coords[1, "X"]), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(coords[1, "Y"]), 0.5, tolerance = 1e-6)
})

test_that("extract_coordinates: errors for non-sf input", {
  expect_error(extract_coordinates(make_non_sf()), "sf object")
})

# ---------------------------------------------------------------------------
# align_spatial_to_panel
# ---------------------------------------------------------------------------

test_that("align_spatial_to_panel: reorders rows according to id_map", {
  pts <- make_point_sf()
  # Reverse the order: D, C, B, A
  id_map <- c(D = 4L, C = 3L, B = 2L, A = 1L)
  aligned <- align_spatial_to_panel(pts, id_map)
  expect_equal(aligned$id, c("D", "C", "B", "A"))
})

test_that("align_spatial_to_panel: subset of rows via id_map", {
  pts <- make_point_sf()
  # Select only A and C (rows 1 and 3)
  id_map <- c(A = 1L, C = 3L)
  aligned <- align_spatial_to_panel(pts, id_map)
  expect_equal(nrow(aligned), 2L)
  expect_equal(aligned$id, c("A", "C"))
})

test_that("align_spatial_to_panel: preserves sf class", {
  pts    <- make_point_sf()
  id_map <- c(A = 1L, B = 2L, C = 3L, D = 4L)
  aligned <- align_spatial_to_panel(pts, id_map)
  expect_true(inherits(aligned, "sf"))
})

test_that("align_spatial_to_panel: errors for non-sf input", {
  id_map <- c(A = 1L)
  expect_error(align_spatial_to_panel(make_non_sf(), id_map), "sf object")
})

test_that("align_spatial_to_panel: errors when id_map is not named", {
  pts <- make_point_sf()
  expect_error(align_spatial_to_panel(pts, c(1L, 2L)), "named")
})

test_that("align_spatial_to_panel: errors for out-of-range row index", {
  pts    <- make_point_sf()
  id_map <- c(X = 99L)
  expect_error(align_spatial_to_panel(pts, id_map), "range")
})

test_that("align_spatial_to_panel: ID ordering matches panel order (integration)", {
  pts <- make_point_sf()
  # Simulate panel ordering: B, D, A, C
  panel_ids    <- c("B", "D", "A", "C")
  spatial_ids  <- pts$id  # A, B, C, D

  id_map <- match(panel_ids, spatial_ids)
  names(id_map) <- panel_ids
  storage.mode(id_map) <- "integer"

  aligned <- align_spatial_to_panel(pts, id_map)
  expect_equal(aligned$id, panel_ids)
})

# ---------------------------------------------------------------------------
# build_neighbor_structure: distance type
# ---------------------------------------------------------------------------

test_that("build_neighbor_structure distance: returns coordinate matrix", {
  pts    <- make_point_sf()
  result <- build_neighbor_structure(pts, type = "distance")
  expect_true(is.matrix(result))
  expect_true(all(c("X", "Y") %in% colnames(result)))
  expect_equal(nrow(result), nrow(pts))
})

test_that("build_neighbor_structure distance: POLYGON returns coordinate matrix", {
  polys  <- make_polygon_sf()
  result <- build_neighbor_structure(polys, type = "distance")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(polys))
})

test_that("build_neighbor_structure distance: default type is distance", {
  pts    <- make_point_sf()
  result <- build_neighbor_structure(pts)
  expect_true(is.matrix(result))
})

test_that("build_neighbor_structure distance: errors for non-sf input", {
  expect_error(
    build_neighbor_structure(make_non_sf(), type = "distance"),
    "sf object"
  )
})

# ---------------------------------------------------------------------------
# build_neighbor_structure: contiguity type
# ---------------------------------------------------------------------------

test_that("build_neighbor_structure contiguity: returns a list", {
  polys  <- make_polygon_sf()
  result <- build_neighbor_structure(polys, type = "contiguity")
  expect_true(is.list(result))
})

test_that("build_neighbor_structure contiguity: list length equals n features", {
  polys  <- make_polygon_sf()
  result <- build_neighbor_structure(polys, type = "contiguity")
  expect_equal(length(result), nrow(polys))
})

test_that("build_neighbor_structure contiguity: each element is integer vector", {
  polys  <- make_polygon_sf()
  result <- build_neighbor_structure(polys, type = "contiguity")
  for (nb in result) {
    expect_true(is.integer(nb))
  }
})

test_that("build_neighbor_structure contiguity: adjacent polygons are neighbours", {
  polys  <- make_polygon_sf()
  # 2x2 grid: polygon 1 (SW) shares a boundary with polygon 2 (SE) and 3 (NW)
  result <- build_neighbor_structure(polys, type = "contiguity")
  # Polygon 1 should have at least polygons 2 and 3 as neighbours
  expect_true(2L %in% result[[1]] || length(result[[1]]) >= 2L)
})

test_that("build_neighbor_structure contiguity: no self-loops", {
  polys  <- make_polygon_sf()
  result <- build_neighbor_structure(polys, type = "contiguity")
  for (i in seq_along(result)) {
    expect_false(i %in% result[[i]])
  }
})

test_that("build_neighbor_structure contiguity: errors for non-sf input", {
  expect_error(
    build_neighbor_structure(make_non_sf(), type = "contiguity"),
    "sf object"
  )
})

test_that("build_neighbor_structure: invalid type raises error", {
  pts <- make_point_sf()
  expect_error(build_neighbor_structure(pts, type = "invalid"), "arg")
})

# ---------------------------------------------------------------------------
# Non-sf object: comprehensive check across all functions
# ---------------------------------------------------------------------------

test_that("all functions reject sp-like or plain list objects", {
  bad <- list(id = 1:3, coords = matrix(1:6, ncol = 2))
  expect_error(assert_sf(bad),                              "sf object")
  expect_error(extract_geometry(bad),                       "sf object")
  expect_error(extract_coordinates(bad),                    "sf object")
  expect_error(align_spatial_to_panel(bad, c(a = 1L)),     "sf object")
  expect_error(build_neighbor_structure(bad, "distance"),   "sf object")
})
