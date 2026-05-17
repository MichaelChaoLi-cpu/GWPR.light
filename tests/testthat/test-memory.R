test_that("estimate_memory returns a gwpr_memory_estimate with correct structure", {
  ctx <- list(
    n_units = 50,
    n_time  = 10,
    n_vars  = 3,
    metadata = list()
  )
  result <- estimate_memory(ctx, workers = 1, cache_distance = TRUE)

  expect_s3_class(result, "gwpr_memory_estimate")
  expect_named(
    result,
    c("n_units", "n_time", "n_vars", "n_rows", "workers",
      "cache_distance", "distance_bytes", "model_bytes",
      "total_bytes", "risk"),
    ignore.order = TRUE
  )
  expect_equal(result$n_rows, 50 * 10)
  expect_equal(result$workers, 1L)
  expect_true(result$cache_distance)
})

test_that("estimate_memory calculates distance_bytes correctly when caching", {
  ctx <- list(n_units = 100, n_time = 5, n_vars = 2, metadata = list())
  result <- estimate_memory(ctx, workers = 1, cache_distance = TRUE)

  expected_dist  <- 100^2 * 8
  expected_model <- 100 * 5 * 2 * 8 * 1
  expect_equal(result$distance_bytes, expected_dist)
  expect_equal(result$model_bytes, expected_model)
  expect_equal(result$total_bytes, expected_dist + expected_model)
})

test_that("estimate_memory sets distance_bytes to 0 when cache_distance = FALSE", {
  ctx <- list(n_units = 100, n_time = 5, n_vars = 2, metadata = list())
  result <- estimate_memory(ctx, workers = 1, cache_distance = FALSE)

  expect_equal(result$distance_bytes, 0)
  expect_gt(result$model_bytes, 0)
})

test_that("estimate_memory uses conservative default (cache_distance = TRUE) when NULL", {
  ctx <- list(n_units = 100, n_time = 5, n_vars = 2, metadata = list())
  result_null <- estimate_memory(ctx, workers = 1, cache_distance = NULL)
  result_true <- estimate_memory(ctx, workers = 1, cache_distance = TRUE)

  expect_equal(result_null$total_bytes, result_true$total_bytes)
  expect_true(result_null$cache_distance)
})

test_that("workers increase causes estimated memory to not decrease", {
  ctx <- list(n_units = 80, n_time = 8, n_vars = 4, metadata = list())

  mem1 <- estimate_memory(ctx, workers = 1)
  mem2 <- estimate_memory(ctx, workers = 2)
  mem4 <- estimate_memory(ctx, workers = 4)

  expect_gte(mem2$total_bytes, mem1$total_bytes)
  expect_gte(mem4$total_bytes, mem2$total_bytes)
})

test_that("workers increase scales model_bytes linearly", {
  ctx <- list(n_units = 50, n_time = 5, n_vars = 3, metadata = list())

  mem1 <- estimate_memory(ctx, workers = 1, cache_distance = FALSE)
  mem3 <- estimate_memory(ctx, workers = 3, cache_distance = FALSE)

  expect_equal(mem3$model_bytes, mem1$model_bytes * 3)
})

test_that("large n_units triggers high memory risk", {
  # n_units = 20000 -> distance matrix = 20000^2 * 8 = 3.2 GB -> high risk
  ctx <- list(n_units = 20000, n_time = 5, n_vars = 2, metadata = list())
  result <- estimate_memory(ctx, workers = 1, cache_distance = TRUE)

  expect_equal(result$risk, "high")
})

test_that("classify_memory_risk returns 'low' below 500 MB", {
  bytes_low <- 100 * 1024^2   # 100 MB
  expect_equal(classify_memory_risk(bytes_low), "low")
})

test_that("classify_memory_risk returns 'medium' between 500 MB and 2 GB", {
  bytes_med <- 1 * 1024^3     # 1 GB
  expect_equal(classify_memory_risk(bytes_med), "medium")
})

test_that("classify_memory_risk returns 'high' above 2 GB", {
  bytes_high <- 3 * 1024^3    # 3 GB
  expect_equal(classify_memory_risk(bytes_high), "high")
})

test_that("classify_memory_risk handles boundary values", {
  mb_500 <- 500 * 1024^2
  gb_2   <- 2   * 1024^3

  # Just below 500 MB -> low
  expect_equal(classify_memory_risk(mb_500 - 1), "low")
  # Exactly 500 MB -> medium
  expect_equal(classify_memory_risk(mb_500), "medium")
  # Exactly 2 GB -> medium
  expect_equal(classify_memory_risk(gb_2), "medium")
  # Just above 2 GB -> high
  expect_equal(classify_memory_risk(gb_2 + 1), "high")
})

test_that("classify_memory_risk rejects invalid input", {
  expect_error(classify_memory_risk("large"), regexp = "numeric")
  expect_error(classify_memory_risk(c(1, 2)), regexp = "single")
  expect_error(classify_memory_risk(-1), regexp = "non-negative")
})

test_that("format_memory_warning returns a character string", {
  ctx  <- list(n_units = 50, n_time = 5, n_vars = 3, metadata = list())
  mem  <- estimate_memory(ctx, workers = 1)
  warn <- format_memory_warning(mem)

  expect_type(warn, "character")
  expect_length(warn, 1L)
})

test_that("format_memory_warning message contains risk level", {
  ctx  <- list(n_units = 50, n_time = 5, n_vars = 3, metadata = list())
  mem  <- estimate_memory(ctx, workers = 1)
  warn <- format_memory_warning(mem)

  expect_true(grepl("LOW|MEDIUM|HIGH", warn))
})

test_that("format_memory_warning for high risk includes actionable suggestions", {
  # Build a high-risk scenario
  ctx  <- list(n_units = 20000, n_time = 5, n_vars = 2, metadata = list())
  mem  <- estimate_memory(ctx, workers = 4, cache_distance = TRUE)
  warn <- format_memory_warning(mem)

  expect_equal(mem$risk, "high")
  # Should mention workers
  expect_true(grepl("workers", warn, ignore.case = TRUE))
  # Should mention cache
  expect_true(grepl("cache_distance", warn, ignore.case = TRUE))
  # Should mention spatial units or n_units
  expect_true(grepl("spatial units|n_units", warn, ignore.case = TRUE))
})

test_that("format_memory_warning for medium risk contains a monitoring note", {
  # Build a medium-risk scenario via classify_memory_risk directly:
  # 800 MB is between 500 MB and 2 GB -> medium
  bytes_medium <- 800 * 1024^2
  ctx <- list(n_units = 10, n_time = 2, n_vars = 1, metadata = list())
  mem <- estimate_memory(ctx, workers = 1, cache_distance = FALSE)
  # Override risk to medium to test the formatting branch independently
  mem$total_bytes <- bytes_medium
  mem$risk        <- classify_memory_risk(bytes_medium)
  expect_equal(mem$risk, "medium")

  warn <- format_memory_warning(mem)
  expect_true(grepl("monitor", warn, ignore.case = TRUE))
})

test_that("format_memory_warning rejects non-estimate input", {
  expect_error(format_memory_warning(list()), regexp = "gwpr_memory_estimate")
  expect_error(format_memory_warning("text"), regexp = "gwpr_memory_estimate")
})

test_that("estimate_memory resolves dimensions from context metadata", {
  ctx <- list(
    metadata = list(n_units = 30, n_time = 4, n_vars = 2)
  )
  result <- estimate_memory(ctx, workers = 1)
  expect_equal(result$n_units, 30)
  expect_equal(result$n_time, 4)
  expect_equal(result$n_vars, 2)
})

test_that("estimate_memory works with gwpr_context object", {
  ctx <- structure(
    list(
      n_units  = 40,
      n_time   = 6,
      n_vars   = 3,
      metadata = list(),
      warnings = character()
    ),
    class = "gwpr_context"
  )
  result <- estimate_memory(ctx, workers = 1)
  expect_s3_class(result, "gwpr_memory_estimate")
  expect_equal(result$n_units, 40)
})

test_that("estimate_memory does not stop execution (only informs)", {
  ctx <- list(n_units = 20000, n_time = 10, n_vars = 5, metadata = list())
  # Must NOT throw an error even for very large (high risk) estimates
  expect_no_error(estimate_memory(ctx, workers = 1, cache_distance = TRUE))
})

test_that("estimate_memory stops with informative error for missing n_units", {
  ctx <- list(n_time = 5, n_vars = 3, metadata = list())
  expect_error(estimate_memory(ctx, workers = 1), regexp = "n_units")
})

test_that("estimate_memory stops with informative error for missing n_time", {
  ctx <- list(n_units = 50, n_vars = 3, metadata = list())
  expect_error(estimate_memory(ctx, workers = 1), regexp = "n_time")
})

test_that("estimate_memory stops with informative error for missing n_vars", {
  ctx <- list(n_units = 50, n_time = 5, metadata = list())
  expect_error(estimate_memory(ctx, workers = 1), regexp = "n_vars")
})

test_that("estimate_memory stops with error for invalid workers", {
  ctx <- list(n_units = 50, n_time = 5, n_vars = 3, metadata = list())
  expect_error(estimate_memory(ctx, workers = 0), regexp = "workers")
  expect_error(estimate_memory(ctx, workers = -1), regexp = "workers")
})
