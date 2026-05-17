library(testthat)

# ---------------------------------------------------------------------------
# Helper: build a minimal valid context
# ---------------------------------------------------------------------------

# Defaults for all core fields; callers may override any of them.
make_minimal_context <- function(
    formula   = y ~ x1,
    family    = "gaussian",
    model     = "within",
    effect    = "individual",
    id        = "id",
    time      = "time",
    kernel    = "bisquare",
    adaptive  = FALSE,
    threshold = 0.5,
    workers   = 1L,
    ...
) {
  new_gwpr_context(
    formula   = formula,
    family    = family,
    model     = model,
    effect    = effect,
    id        = id,
    time      = time,
    kernel    = kernel,
    adaptive  = adaptive,
    threshold = threshold,
    workers   = workers,
    ...
  )
}

# ---------------------------------------------------------------------------
# new_gwpr_context: field completeness
# ---------------------------------------------------------------------------

test_that("new_gwpr_context: returns a list with class gwpr_context", {
  ctx <- new_gwpr_context()
  expect_true(is.list(ctx))
  expect_s3_class(ctx, "gwpr_context")
})

test_that("new_gwpr_context: all required fields are present and default to NULL", {
  ctx <- new_gwpr_context()
  required_null_fields <- c(
    "call", "formula", "family", "model", "effect",
    "id", "time", "kernel", "adaptive", "threshold",
    "workers", "seed",
    "raw_data", "raw_spatial",
    "panel_data", "spatial_data", "id_map",
    "coords", "model_frame", "model_matrix", "response"
  )
  for (fld in required_null_fields) {
    expect_null(ctx[[fld]], label = paste("ctx$", fld, "is NULL by default"))
  }
})

test_that("new_gwpr_context: metadata defaults to empty list", {
  ctx <- new_gwpr_context()
  expect_identical(ctx$metadata, list())
})

test_that("new_gwpr_context: warnings defaults to empty character vector", {
  ctx <- new_gwpr_context()
  expect_identical(ctx$warnings, character())
})

test_that("new_gwpr_context: supplied values are stored correctly", {
  ctx <- make_minimal_context(
    seed     = 42L,
    workers  = 1L,
    metadata = list(n_units = 3L),
    warnings = c("test warning")
  )
  expect_equal(deparse(ctx$formula), deparse(y ~ x1))
  expect_identical(ctx$family,  "gaussian")
  expect_identical(ctx$model,   "within")
  expect_identical(ctx$effect,  "individual")
  expect_identical(ctx$id,      "id")
  expect_identical(ctx$time,    "time")
  expect_identical(ctx$kernel,  "bisquare")
  expect_false(ctx$adaptive)
  expect_equal(ctx$threshold, 0.5)
  expect_identical(ctx$workers,  1L)
  expect_identical(ctx$seed,    42L)
  expect_identical(ctx$metadata, list(n_units = 3L))
  expect_identical(ctx$warnings, "test warning")
})

test_that("new_gwpr_context: extra arguments are stored in context", {
  ctx <- new_gwpr_context(custom_field = "hello")
  expect_identical(ctx$custom_field, "hello")
})

test_that("new_gwpr_context: data fields can be set and retrieved", {
  df <- data.frame(id = 1:3, time = 1:3, y = rnorm(3), x1 = rnorm(3))
  ctx <- new_gwpr_context(raw_data = df)
  expect_identical(ctx$raw_data, df)
})

# ---------------------------------------------------------------------------
# validate_gwpr_context: valid context passes
# ---------------------------------------------------------------------------

test_that("validate_gwpr_context: fully populated core fields pass", {
  ctx <- make_minimal_context()
  expect_true(validate_gwpr_context(ctx))
})

test_that("validate_gwpr_context: returns invisibly TRUE on success", {
  ctx <- make_minimal_context()
  result <- validate_gwpr_context(ctx)
  expect_true(result)
})

# ---------------------------------------------------------------------------
# validate_gwpr_context: missing core fields raise errors
# ---------------------------------------------------------------------------

test_that("validate_gwpr_context: missing formula raises error", {
  ctx <- make_minimal_context()
  ctx$formula <- NULL
  expect_error(validate_gwpr_context(ctx), "formula")
})

test_that("validate_gwpr_context: missing family raises error", {
  ctx <- make_minimal_context()
  ctx$family <- NULL
  expect_error(validate_gwpr_context(ctx), "family")
})

test_that("validate_gwpr_context: missing id raises error", {
  ctx <- make_minimal_context()
  ctx$id <- NULL
  expect_error(validate_gwpr_context(ctx), "id")
})

test_that("validate_gwpr_context: missing time raises error", {
  ctx <- make_minimal_context()
  ctx$time <- NULL
  expect_error(validate_gwpr_context(ctx), "time")
})

test_that("validate_gwpr_context: missing model raises error", {
  ctx <- make_minimal_context()
  ctx$model <- NULL
  expect_error(validate_gwpr_context(ctx), "model")
})

test_that("validate_gwpr_context: missing effect raises error", {
  ctx <- make_minimal_context()
  ctx$effect <- NULL
  expect_error(validate_gwpr_context(ctx), "effect")
})

test_that("validate_gwpr_context: missing kernel raises error", {
  ctx <- make_minimal_context()
  ctx$kernel <- NULL
  expect_error(validate_gwpr_context(ctx), "kernel")
})

test_that("validate_gwpr_context: missing adaptive raises error", {
  ctx <- make_minimal_context()
  ctx$adaptive <- NULL
  expect_error(validate_gwpr_context(ctx), "adaptive")
})

test_that("validate_gwpr_context: missing threshold raises error", {
  ctx <- make_minimal_context()
  ctx$threshold <- NULL
  expect_error(validate_gwpr_context(ctx), "threshold")
})

test_that("validate_gwpr_context: missing workers raises error", {
  ctx <- make_minimal_context()
  ctx$workers <- NULL
  expect_error(validate_gwpr_context(ctx), "workers")
})

test_that("validate_gwpr_context: multiple missing fields listed in error message", {
  ctx <- make_minimal_context()
  ctx$formula <- NULL
  ctx$family  <- NULL
  ctx$id      <- NULL
  err <- tryCatch(validate_gwpr_context(ctx), error = function(e) conditionMessage(e))
  expect_match(err, "formula")
  expect_match(err, "family")
  expect_match(err, "id")
})

test_that("validate_gwpr_context: non-list input raises error", {
  expect_error(validate_gwpr_context("not a list"), "`context` must be a list")
})

# ---------------------------------------------------------------------------
# Interface consistency: context can be reused across module boundaries
# ---------------------------------------------------------------------------

test_that("context fields cover search module requirements (bandwidth, seed, workers)", {
  ctx <- make_minimal_context(seed = 99L, workers = 1L)
  # Fields required by bandwidth search modules
  expect_false(is.null(ctx$kernel))
  expect_false(is.null(ctx$adaptive))
  expect_false(is.null(ctx$workers))
  expect_identical(ctx$seed, 99L)
})

test_that("context fields cover fit module requirements (formula, model, effect)", {
  ctx <- make_minimal_context()
  expect_false(is.null(ctx$formula))
  expect_false(is.null(ctx$model))
  expect_false(is.null(ctx$effect))
})

test_that("context fields cover diagnostics module requirements (family, threshold)", {
  ctx <- make_minimal_context()
  expect_false(is.null(ctx$family))
  expect_false(is.null(ctx$threshold))
})

test_that("context supports binomial family for logistic GWPR", {
  ctx <- make_minimal_context(family = "binomial", threshold = 0.3)
  expect_identical(ctx$family, "binomial")
  expect_equal(ctx$threshold, 0.3)
})

test_that("context supports adaptive bandwidth mode", {
  ctx <- make_minimal_context(adaptive = TRUE)
  expect_true(ctx$adaptive)
})

test_that("context metadata and warnings can be updated after construction", {
  ctx <- make_minimal_context()
  ctx$metadata$n_units       <- 5L
  ctx$metadata$unbalanced    <- FALSE
  ctx$warnings               <- c(ctx$warnings, "added warning")

  expect_identical(ctx$metadata$n_units,   5L)
  expect_false(ctx$metadata$unbalanced)
  expect_identical(ctx$warnings, "added warning")
})
