library(testthat)

# ---------------------------------------------------------------------------
# Helpers: build minimal simulation data
# ---------------------------------------------------------------------------

#' Create a small balanced panel data frame
make_panel_data <- function(ids = c("A", "B", "C"),
                            times = 1:3,
                            id_col = "id",
                            time_col = "time") {
  n_ids   <- length(ids)
  n_times <- length(times)
  df <- data.frame(
    id   = rep(ids,   each = n_times),
    time = rep(times, times = n_ids),
    y    = rnorm(n_ids * n_times),
    x1   = rnorm(n_ids * n_times),
    stringsAsFactors = FALSE
  )
  names(df)[names(df) == "id"]   <- id_col
  names(df)[names(df) == "time"] <- time_col
  df
}

#' Create a small sf POINT object aligned with given IDs
make_spatial_sf <- function(ids = c("A", "B", "C"), id_col = "id",
                             extra_ids = character(0)) {
  all_ids <- c(ids, extra_ids)
  n       <- length(all_ids)
  pts     <- sf::st_sfc(lapply(seq_len(n), function(i) sf::st_point(c(i, i))))
  sf_obj  <- sf::st_sf(
    setNames(list(all_ids), id_col),
    geometry = pts,
    stringsAsFactors = FALSE
  )
  sf_obj
}

#' Build a minimal gwpr_context with raw data and spatial
make_context <- function(data, spatial,
                         formula  = y ~ x1,
                         family   = "gaussian",
                         model    = "within",
                         effect   = "individual",
                         id       = "id",
                         time     = "time",
                         kernel   = "bisquare",
                         adaptive = FALSE,
                         threshold = 0.5,
                         workers   = 1L) {
  new_gwpr_context(
    formula     = formula,
    family      = family,
    model       = model,
    effect      = effect,
    id          = id,
    time        = time,
    kernel      = kernel,
    adaptive    = adaptive,
    threshold   = threshold,
    workers     = workers,
    raw_data    = data,
    raw_spatial = spatial
  )
}

# ---------------------------------------------------------------------------
# standardize_binary_response
# ---------------------------------------------------------------------------

test_that("standardize_binary_response: numeric 0/1 returned as numeric", {
  y <- c(0, 1, 0, 1)
  result <- standardize_binary_response(y)
  expect_true(is.numeric(result))
  expect_equal(result, c(0, 1, 0, 1))
})

test_that("standardize_binary_response: logical converted to 0/1 integer", {
  y <- c(TRUE, FALSE, TRUE)
  result <- standardize_binary_response(y)
  expect_true(is.numeric(result) || is.integer(result))
  expect_equal(as.integer(result), c(1L, 0L, 1L))
})

test_that("standardize_binary_response: two-level factor first level is 0", {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  result <- standardize_binary_response(y)
  expect_equal(as.integer(result), c(0L, 1L, 0L, 1L))
})

test_that("standardize_binary_response: factor with >2 levels raises error", {
  y <- factor(c("a", "b", "c"))
  expect_error(standardize_binary_response(y), "two-level factor")
})

test_that("standardize_binary_response: numeric with invalid values raises error", {
  y <- c(0, 1, 2)
  expect_error(standardize_binary_response(y), "0 and 1")
})

test_that("standardize_binary_response: non-numeric/logical/factor raises error", {
  expect_error(standardize_binary_response("a"), "0/1 numeric")
})

# ---------------------------------------------------------------------------
# build_id_map
# ---------------------------------------------------------------------------

test_that("build_id_map: returns named integer vector with correct mappings", {
  panel_data   <- make_panel_data(ids = c("A", "B", "C"))
  spatial_data <- make_spatial_sf(ids = c("A", "B", "C"))

  id_map <- build_id_map(panel_data, spatial_data, id = "id")

  expect_true(is.integer(id_map))
  expect_equal(sort(names(id_map)), sort(c("A", "B", "C")))
  # Each index should be a valid row in spatial_data
  expect_true(all(id_map >= 1L & id_map <= nrow(spatial_data)))
})

test_that("build_id_map: extra spatial rows do not cause an error", {
  panel_data   <- make_panel_data(ids = c("A", "B"))
  # Spatial has an extra ID "C" not in panel
  spatial_data <- make_spatial_sf(ids = c("A", "B", "C"))

  expect_no_error(build_id_map(panel_data, spatial_data, id = "id"))
})

test_that("build_id_map: panel ID missing from spatial raises error", {
  panel_data   <- make_panel_data(ids = c("A", "B", "X"))
  spatial_data <- make_spatial_sf(ids = c("A", "B", "C"))

  expect_error(
    build_id_map(panel_data, spatial_data, id = "id"),
    "X"
  )
})

test_that("build_id_map: numeric IDs handled correctly", {
  panel_data   <- make_panel_data(ids = c(1L, 2L, 3L), id_col = "uid")
  panel_data$uid <- as.integer(panel_data$uid)
  spatial_data <- make_spatial_sf(ids = c(1L, 2L, 3L), id_col = "uid")
  spatial_data$uid <- as.integer(spatial_data$uid)

  id_map <- build_id_map(panel_data, spatial_data, id = "uid")
  expect_equal(length(id_map), 3L)
  expect_true(all(!is.na(id_map)))
})

# ---------------------------------------------------------------------------
# prepare_panel_data
# ---------------------------------------------------------------------------

test_that("prepare_panel_data: adds raw_row_id column", {
  data    <- make_panel_data()
  spatial <- make_spatial_sf()
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_panel_data(ctx)
  expect_true("raw_row_id" %in% names(ctx2$panel_data))
})

test_that("prepare_panel_data: sorts by id then time", {
  # Create unsorted data
  data    <- make_panel_data(ids = c("B", "A"), times = c(2, 1))
  spatial <- make_spatial_sf(ids = c("A", "B"))
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_panel_data(ctx)
  pd   <- ctx2$panel_data
  # After sorting, first id should be A
  expect_equal(pd$id[1], "A")
})

test_that("prepare_panel_data: metadata records panel balance (balanced)", {
  data    <- make_panel_data(ids = c("A", "B"), times = 1:3)
  spatial <- make_spatial_sf(ids = c("A", "B"))
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_panel_data(ctx)
  expect_true(ctx2$metadata$panel_balanced)
})

test_that("prepare_panel_data: unbalanced panel does not raise error", {
  # B has only 2 time periods
  df <- rbind(
    make_panel_data(ids = "A", times = 1:3),
    make_panel_data(ids = "B", times = 1:2)
  )
  spatial <- make_spatial_sf(ids = c("A", "B"))
  ctx     <- make_context(df, spatial)

  expect_no_error(prepare_panel_data(ctx))
  ctx2 <- prepare_panel_data(ctx)
  expect_false(ctx2$metadata$panel_balanced)
})

test_that("prepare_panel_data: within model records single-obs individuals in metadata", {
  # Create data where "C" has only 1 observation
  df <- rbind(
    make_panel_data(ids = "A", times = 1:3),
    make_panel_data(ids = "B", times = 1:3),
    make_panel_data(ids = "C", times = 1)
  )
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(df, spatial, model = "within")

  ctx2 <- suppressWarnings(prepare_panel_data(ctx))
  expect_true("within_single_obs_ids" %in% names(ctx2$metadata))
  expect_true("C" %in% ctx2$metadata$within_single_obs_ids)
})

# ---------------------------------------------------------------------------
# prepare_spatial_data
# ---------------------------------------------------------------------------

test_that("prepare_spatial_data: coords matrix has X and Y columns", {
  data    <- make_panel_data()
  spatial <- make_spatial_sf()
  ctx     <- make_context(data, spatial)
  ctx     <- prepare_panel_data(ctx)

  ctx2 <- prepare_spatial_data(ctx)
  expect_true(!is.null(ctx2$coords))
  expect_true(all(c("X", "Y") %in% colnames(ctx2$coords)))
})

test_that("prepare_spatial_data: extra spatial regions not in panel are excluded", {
  data    <- make_panel_data(ids = c("A", "B"))
  # Spatial has extra ID "C"
  spatial <- make_spatial_sf(ids = c("A", "B"), extra_ids = "C")
  ctx     <- make_context(data, spatial)
  ctx     <- prepare_panel_data(ctx)

  ctx2 <- prepare_spatial_data(ctx)
  expect_equal(nrow(ctx2$spatial_data), 2L)
  expect_false("C" %in% ctx2$spatial_data$id)
})

test_that("prepare_spatial_data: panel ID missing from spatial raises error", {
  data    <- make_panel_data(ids = c("A", "B", "X"))
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(data, spatial)
  ctx     <- prepare_panel_data(ctx)

  expect_error(prepare_spatial_data(ctx), "X")
})

# ---------------------------------------------------------------------------
# build_model_frame and build_model_matrix
# ---------------------------------------------------------------------------

test_that("build_model_frame: produces a model frame from formula", {
  data    <- make_panel_data()
  spatial <- make_spatial_sf()
  ctx     <- make_context(data, spatial)
  ctx     <- prepare_panel_data(ctx)

  ctx2 <- build_model_frame(ctx)
  expect_true(!is.null(ctx2$model_frame))
  expect_true(is.data.frame(ctx2$model_frame))
})

test_that("build_model_matrix: gaussian response is numeric", {
  data    <- make_panel_data()
  spatial <- make_spatial_sf()
  ctx     <- make_context(data, spatial, family = "gaussian")
  ctx     <- prepare_panel_data(ctx)
  ctx     <- build_model_frame(ctx)
  ctx     <- build_model_matrix(ctx)

  expect_true(is.numeric(ctx$response))
  expect_true(is.matrix(ctx$model_matrix))
})

test_that("build_model_matrix: binomial factor response standardised to 0/1", {
  data    <- make_panel_data()
  data$y  <- factor(
    ifelse(data$y > 0, "yes", "no"),
    levels = c("no", "yes")
  )
  spatial <- make_spatial_sf()
  ctx     <- make_context(data, spatial, family = "binomial",
                          formula = y ~ x1)
  ctx     <- prepare_panel_data(ctx)
  ctx     <- build_model_frame(ctx)
  ctx     <- build_model_matrix(ctx)

  expect_true(all(ctx$response %in% c(0, 1)))
})

# ---------------------------------------------------------------------------
# prepare_data  (full pipeline)
# ---------------------------------------------------------------------------

test_that("prepare_data: full pipeline with character IDs", {
  data    <- make_panel_data(ids = c("A", "B", "C"))
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_data(ctx)

  expect_false(is.null(ctx2$panel_data))
  expect_false(is.null(ctx2$spatial_data))
  expect_false(is.null(ctx2$id_map))
  expect_false(is.null(ctx2$coords))
  expect_false(is.null(ctx2$model_frame))
  expect_false(is.null(ctx2$model_matrix))
  expect_false(is.null(ctx2$response))
})

test_that("prepare_data: full pipeline with numeric IDs", {
  data    <- make_panel_data(ids = c(1, 2, 3))
  spatial <- make_spatial_sf(ids = c(1, 2, 3))
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_data(ctx)

  expect_false(is.null(ctx2$panel_data))
  expect_false(is.null(ctx2$id_map))
  expect_equal(length(ctx2$id_map), 3L)
})

test_that("prepare_data: unbalanced panel completes without error", {
  df <- rbind(
    make_panel_data(ids = "A", times = 1:4),
    make_panel_data(ids = "B", times = 1:2),
    make_panel_data(ids = "C", times = 1:3)
  )
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(df, spatial, model = "pooling")

  expect_no_error(prepare_data(ctx))
  ctx2 <- prepare_data(ctx)
  expect_false(ctx2$metadata$panel_balanced)
})

test_that("prepare_data: spatial extra regions excluded, panel IDs all matched", {
  data    <- make_panel_data(ids = c("A", "B"))
  spatial <- make_spatial_sf(ids = c("A", "B"), extra_ids = c("D", "E"))
  ctx     <- make_context(data, spatial)

  ctx2 <- prepare_data(ctx)
  # spatial_data should contain only A and B
  expect_equal(nrow(ctx2$spatial_data), 2L)
  expect_equal(sort(ctx2$spatial_data$id), sort(c("A", "B")))
})

test_that("prepare_data: panel ID missing in spatial raises error", {
  data    <- make_panel_data(ids = c("A", "B", "Z"))
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(data, spatial)

  expect_error(prepare_data(ctx), "Z")
})

test_that("prepare_data: binomial logistic factor response converted to 0/1", {
  data   <- make_panel_data(ids = c("A", "B", "C"))
  data$y <- factor(ifelse(data$y > 0, "case", "control"),
                   levels = c("control", "case"))
  spatial <- make_spatial_sf(ids = c("A", "B", "C"))
  ctx     <- make_context(data, spatial, family = "binomial", formula = y ~ x1)

  ctx2 <- prepare_data(ctx)
  expect_true(all(ctx2$response %in% c(0L, 1L)))
})
