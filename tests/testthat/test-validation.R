library(testthat)

# ---------------------------------------------------------------------------
# Helper: small simulated panel data and sf spatial object
# ---------------------------------------------------------------------------

make_test_data <- function() {
  data.frame(
    id   = rep(c("A", "B", "C"), each = 3),
    time = rep(1:3, times = 3),
    y    = c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9),
    x1   = rnorm(9, mean = 0, sd = 1),
    x2   = rnorm(9, mean = 5, sd = 2),
    stringsAsFactors = FALSE
  )
}

make_test_spatial <- function() {
  # Three simple point features aligned with ids A, B, C
  sf::st_as_sf(
    data.frame(
      id = c("A", "B", "C"),
      lon = c(0, 1, 2),
      lat = c(0, 1, 2)
    ),
    coords = c("lon", "lat"),
    crs = 4326
  )
}

# ---------------------------------------------------------------------------
# validate_formula
# ---------------------------------------------------------------------------

test_that("validate_formula: valid formula passes", {
  df <- make_test_data()
  expect_true(validate_formula(y ~ x1 + x2, df))
})

test_that("validate_formula: non-formula input raises error", {
  df <- make_test_data()
  expect_error(validate_formula("y ~ x1", df), "`formula` must be a formula object")
})

test_that("validate_formula: missing response variable raises error", {
  df <- make_test_data()
  expect_error(validate_formula(missing_y ~ x1, df), "Response variable")
})

test_that("validate_formula: missing predictor variable raises error", {
  df <- make_test_data()
  expect_error(validate_formula(y ~ x1 + x_missing, df), "Predictor variable")
})

test_that("validate_formula: non-data.frame raises error", {
  expect_error(validate_formula(y ~ x1, list(y = 1, x1 = 2)), "`data` must be a data.frame")
})

# ---------------------------------------------------------------------------
# validate_panel_index
# ---------------------------------------------------------------------------

test_that("validate_panel_index: valid columns pass", {
  df <- make_test_data()
  expect_true(validate_panel_index(df, "id", "time"))
})

test_that("validate_panel_index: missing id column raises error", {
  df <- make_test_data()
  expect_error(validate_panel_index(df, "no_such_id", "time"), "not found in `data`")
})

test_that("validate_panel_index: missing time column raises error", {
  df <- make_test_data()
  expect_error(validate_panel_index(df, "id", "no_such_time"), "not found in `data`")
})

test_that("validate_panel_index: non-character id raises error", {
  df <- make_test_data()
  expect_error(validate_panel_index(df, 1L, "time"), "`id` must be a single character string")
})

# ---------------------------------------------------------------------------
# validate_spatial
# ---------------------------------------------------------------------------

test_that("validate_spatial: valid sf object with id column passes", {
  sp <- make_test_spatial()
  expect_true(validate_spatial(sp, "id"))
})

test_that("validate_spatial: non-sf object raises error", {
  df <- make_test_data()
  expect_error(validate_spatial(df, "id"), "must be an sf object")
})

test_that("validate_spatial: missing id column in sf raises error", {
  sp <- make_test_spatial()
  expect_error(validate_spatial(sp, "no_such_id"), "not found in `spatial`")
})

test_that("validate_spatial: sp object raises informative error", {
  # A plain data.frame (not sp, not sf) should fail with the sf-required message
  fake_sp <- structure(list(), class = "SpatialPointsDataFrame")
  expect_error(validate_spatial(fake_sp, "id"), "must be an sf object")
})

# ---------------------------------------------------------------------------
# validate_family_response
# ---------------------------------------------------------------------------

test_that("validate_family_response: gaussian with numeric response passes", {
  df <- make_test_data()
  expect_true(validate_family_response(df, y ~ x1, "gaussian"))
})

test_that("validate_family_response: gaussian with non-numeric response raises error", {
  df <- make_test_data()
  df$y <- as.character(df$y)
  expect_error(validate_family_response(df, y ~ x1, "gaussian"), "must be numeric")
})

test_that("validate_family_response: binomial with 0/1 numeric passes", {
  df <- make_test_data()
  df$y <- rep(c(0L, 1L, 0L), 3)
  expect_true(validate_family_response(df, y ~ x1, "binomial"))
})

test_that("validate_family_response: binomial with logical passes", {
  df <- make_test_data()
  df$y <- rep(c(TRUE, FALSE, TRUE), 3)
  expect_true(validate_family_response(df, y ~ x1, "binomial"))
})

test_that("validate_family_response: binomial with two-level factor passes", {
  df <- make_test_data()
  df$y <- factor(rep(c("yes", "no", "yes"), 3))
  expect_true(validate_family_response(df, y ~ x1, "binomial"))
})

test_that("validate_family_response: binomial with multi-level factor raises error with 1.0.0 message", {
  df <- make_test_data()
  df$y <- factor(rep(c("a", "b", "c"), 3))
  expect_error(
    validate_family_response(df, y ~ x1, "binomial"),
    "1.0.0 does not support multinomial"
  )
})

test_that("validate_family_response: binomial with numeric not 0/1 raises error", {
  df <- make_test_data()
  df$y <- c(0, 1, 2, 0, 1, 0, 1, 0, 1)
  expect_error(validate_family_response(df, y ~ x1, "binomial"), "must contain only 0 and 1")
})

test_that("validate_family_response: unknown family raises error", {
  df <- make_test_data()
  expect_error(validate_family_response(df, y ~ x1, "poisson"))
})

# ---------------------------------------------------------------------------
# validate_workers
# ---------------------------------------------------------------------------

test_that("validate_workers: workers = 1 passes", {
  expect_true(validate_workers(1L))
})

test_that("validate_workers: workers = 4 passes", {
  expect_true(validate_workers(4L))
})

test_that("validate_workers: workers = 0 raises error", {
  expect_error(validate_workers(0), "positive integer")
})

test_that("validate_workers: workers = -1 raises error", {
  expect_error(validate_workers(-1), "positive integer")
})

test_that("validate_workers: non-integer numeric raises error", {
  expect_error(validate_workers(1.5), "positive integer")
})

test_that("validate_workers: character raises error", {
  expect_error(validate_workers("2"), "positive integer")
})

test_that("validate_workers: NA raises error", {
  expect_error(validate_workers(NA_real_), "positive integer")
})

# ---------------------------------------------------------------------------
# validate_bandwidth_control
# ---------------------------------------------------------------------------

test_that("validate_bandwidth_control: valid grid control passes", {
  ctrl <- list(lower = 10, upper = 100, step = 10)
  expect_true(validate_bandwidth_control("grid", ctrl, adaptive = FALSE))
})

test_that("validate_bandwidth_control: grid missing lower raises error", {
  ctrl <- list(upper = 100, step = 10)
  expect_error(validate_bandwidth_control("grid", ctrl, FALSE), "lower")
})

test_that("validate_bandwidth_control: grid missing upper raises error", {
  ctrl <- list(lower = 10, step = 10)
  expect_error(validate_bandwidth_control("grid", ctrl, FALSE), "upper")
})

test_that("validate_bandwidth_control: grid missing step raises error", {
  ctrl <- list(lower = 10, upper = 100)
  expect_error(validate_bandwidth_control("grid", ctrl, FALSE), "step")
})

test_that("validate_bandwidth_control: lower >= upper raises error", {
  ctrl <- list(lower = 100, upper = 10, step = 5)
  expect_error(validate_bandwidth_control("grid", ctrl, FALSE), "must be less than")
})

test_that("validate_bandwidth_control: step <= 0 raises error", {
  ctrl <- list(lower = 10, upper = 100, step = -1)
  expect_error(validate_bandwidth_control("grid", ctrl, FALSE), "must be positive")
})

test_that("validate_bandwidth_control: adaptive grid with non-integer lower raises error", {
  ctrl <- list(lower = 1.5, upper = 10, step = 1)
  expect_error(validate_bandwidth_control("grid", ctrl, adaptive = TRUE), "positive integer")
})

test_that("validate_bandwidth_control: adaptive grid with valid integer lower passes", {
  ctrl <- list(lower = 2, upper = 20, step = 2)
  expect_true(validate_bandwidth_control("grid", ctrl, adaptive = TRUE))
})

test_that("validate_bandwidth_control: sgd method with empty control passes", {
  expect_true(validate_bandwidth_control("sgd", list(), adaptive = FALSE))
})

test_that("validate_bandwidth_control: random method with empty control passes", {
  expect_true(validate_bandwidth_control("random", list(), adaptive = FALSE))
})

test_that("validate_bandwidth_control: non-list control raises error", {
  expect_error(validate_bandwidth_control("grid", "not_a_list", FALSE), "`control` must be a list")
})

# ---------------------------------------------------------------------------
# validate_inputs (top-level)
# ---------------------------------------------------------------------------

test_that("validate_inputs: all valid inputs pass", {
  df <- make_test_data()
  sp <- make_test_spatial()
  expect_true(
    validate_inputs(
      formula  = y ~ x1 + x2,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "within",
      effect   = "individual",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 1L
    )
  )
})

test_that("validate_inputs: invalid model raises error", {
  df <- make_test_data()
  sp <- make_test_spatial()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "fixed",      # not a valid model
      effect   = "individual",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 1L
    )
  )
})

test_that("validate_inputs: invalid effect raises error", {
  df <- make_test_data()
  sp <- make_test_spatial()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "within",
      effect   = "bad_effect",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 1L
    )
  )
})

test_that("validate_inputs: invalid kernel raises error", {
  df <- make_test_data()
  sp <- make_test_spatial()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "within",
      effect   = "individual",
      kernel   = "linear",     # not supported
      adaptive = FALSE,
      workers  = 1L
    )
  )
})

test_that("validate_inputs: workers = 0 raises error", {
  df <- make_test_data()
  sp <- make_test_spatial()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "within",
      effect   = "individual",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 0L
    ),
    "positive integer"
  )
})

test_that("validate_inputs: non-sf spatial raises error", {
  df <- make_test_data()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = df,           # data.frame, not sf
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "within",
      effect   = "individual",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 1L
    ),
    "must be an sf object"
  )
})

test_that("validate_inputs: binomial multi-level factor raises 1.0.0 error", {
  df <- make_test_data()
  df$y <- factor(rep(c("a", "b", "c"), 3))
  sp <- make_test_spatial()
  expect_error(
    validate_inputs(
      formula  = y ~ x1,
      data     = df,
      spatial  = sp,
      id       = "id",
      time     = "time",
      family   = "binomial",
      model    = "pooling",
      effect   = "individual",
      kernel   = "bisquare",
      adaptive = FALSE,
      workers  = 1L
    ),
    "1.0.0 does not support multinomial"
  )
})
