## tests/testthat/test-model-linear.R
## Tests for the Linear GWPR Engine (R/model_linear.R)
##
## All tests use small simulated panel data.  workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helpers: small simulated panel data
# ---------------------------------------------------------------------------

# 4 spatial units, 5 time periods each => 20 rows
make_panel_data <- function(seed = 42L) {
  set.seed(seed)
  n_units <- 4L
  n_times <- 5L
  id   <- rep(as.character(1:n_units), each = n_times)
  time <- rep(1:n_times, times = n_units)
  x1   <- rnorm(n_units * n_times)
  x2   <- rnorm(n_units * n_times)
  y    <- 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n_units * n_times, sd = 0.2)
  data.frame(id = id, time = time, y = y, x1 = x1, x2 = x2,
             stringsAsFactors = FALSE)
}

# 4 units on a small 2x2 grid
make_coords_4 <- function() {
  m <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2, byrow = TRUE,
    dimnames = list(c("1", "2", "3", "4"), c("X", "Y"))
  )
  m
}

make_context <- function(model = "pooling", effect = "individual",
                          kernel = "bisquare", adaptive = FALSE,
                          seed = NULL) {
  new_gwpr_context(
    formula      = y ~ x1 + x2,
    family       = "gaussian",
    model        = model,
    effect       = effect,
    id           = "id",
    time         = "time",
    kernel       = kernel,
    adaptive     = adaptive,
    threshold    = 0.5,
    workers      = 1L,
    seed         = seed,
    panel_data   = make_panel_data(),
    coords       = make_coords_4()
  )
}

# ---------------------------------------------------------------------------
# effect_to_plm
# ---------------------------------------------------------------------------

test_that("effect_to_plm maps individual correctly", {
  expect_equal(effect_to_plm("individual"), "individual")
})

test_that("effect_to_plm maps time correctly", {
  expect_equal(effect_to_plm("time"), "time")
})

test_that("effect_to_plm maps two-way to twoways", {
  expect_equal(effect_to_plm("two-way"), "twoways")
})

test_that("effect_to_plm stops for nested without hierarchy info", {
  expect_error(effect_to_plm("nested"), regexp = "nested")
})

test_that("effect_to_plm stops for unknown effect", {
  expect_error(effect_to_plm("unknown"), regexp = "Unknown effect")
})

# ---------------------------------------------------------------------------
# fit_linear_local_model  — pooling
# ---------------------------------------------------------------------------

test_that("fit_linear_local_model: pooling returns ok status", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_false(is.null(res$fit))
  expect_true(inherits(res$fit, "lm"))
})

# ---------------------------------------------------------------------------
# fit_linear_local_model  — within
# ---------------------------------------------------------------------------

test_that("fit_linear_local_model: within individual returns ok status", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "within",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_false(is.null(res$fit))
  expect_true(inherits(res$fit, "plm"))
})

test_that("fit_linear_local_model: within records single-obs individuals in metadata", {
  # Create data with unit 99 having only 1 observation
  pd_base <- make_panel_data()
  extra   <- data.frame(id = "99", time = 1, y = 0, x1 = 0, x2 = 0,
                        stringsAsFactors = FALSE)
  pd      <- rbind(pd_base, extra)
  wts     <- rep(1, nrow(pd))

  res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "within",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  # Metadata should record the single-obs individual
  expect_true("single_obs_individuals" %in% names(res$metadata))
  expect_true("99" %in% res$metadata$single_obs_individuals)
})

# ---------------------------------------------------------------------------
# fit_linear_local_model  — random
# ---------------------------------------------------------------------------

test_that("fit_linear_local_model: random individual returns ok status", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "random",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_false(is.null(res$fit))
  expect_true(inherits(res$fit, "plm"))
})

# ---------------------------------------------------------------------------
# fit_linear_local_model  — error handling
# ---------------------------------------------------------------------------

test_that("fit_linear_local_model: bad formula gives failed status, not error", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  # Formula references a non-existent variable
  res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ nonexistent_var,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "failed")
  expect_null(res$fit)
  expect_true(is.character(res$error) && nchar(res$error) > 0L)
})

# ---------------------------------------------------------------------------
# extract_linear_local_result
# ---------------------------------------------------------------------------

test_that("extract_linear_local_result: ok result has named numeric vectors", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  ext <- extract_linear_local_result(local_res)

  expect_equal(ext$status, "ok")
  expect_true(is.numeric(ext$coefficients))
  expect_true(is.numeric(ext$se))
  expect_true(is.numeric(ext$tvalues))
  expect_true(is.numeric(ext$local_r2) || is.na(ext$local_r2))
  expect_true(is.numeric(ext$local_aic) || is.na(ext$local_aic))
})

test_that("extract_linear_local_result: failed result returns NA placeholders", {
  failed_res <- list(status = "failed", fit = NULL, error = "test error",
                     metadata = list())
  ext <- extract_linear_local_result(failed_res)

  expect_equal(ext$status, "failed")
  expect_true(is.na(ext$coefficients) || all(is.na(ext$coefficients)))
  expect_true(is.na(ext$local_r2)  || all(is.na(ext$local_r2)))
  expect_true(is.na(ext$local_aic) || all(is.na(ext$local_aic)))
  expect_equal(ext$error, "test error")
})

# ---------------------------------------------------------------------------
# predict_linear_local_model
# ---------------------------------------------------------------------------

test_that("predict_linear_local_model: returns numeric vector of correct length", {
  pd  <- make_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_linear_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  preds <- predict_linear_local_model(local_res, pd)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(pd))
})

test_that("predict_linear_local_model: failed model returns NA vector", {
  failed_res <- list(status = "failed", fit = NULL, error = "err",
                     metadata = list())
  pd   <- make_panel_data()
  preds <- predict_linear_local_model(failed_res, pd)
  expect_true(all(is.na(preds)))
  expect_equal(length(preds), nrow(pd))
})

# ---------------------------------------------------------------------------
# fit_linear_gwpr — smoke tests
# ---------------------------------------------------------------------------

test_that("fit_linear_gwpr: pooling smoke test returns correct structure", {
  ctx <- make_context(model = "pooling", effect = "individual",
                      kernel = "bisquare", adaptive = TRUE)
  res <- fit_linear_gwpr(ctx, bandwidth = 3L)

  expect_true(is.list(res))
  expect_true("local_results" %in% names(res))
  expect_true("predictions"   %in% names(res))
  expect_true("residuals"     %in% names(res))
  expect_true("metrics"       %in% names(res))
  expect_true("metadata"      %in% names(res))

  # 4 spatial units => 4 local results
  expect_equal(length(res$local_results), 4L)
  # Metrics fields
  expect_true(all(c("R2", "MSE", "RMSE", "MAE") %in% names(res$metrics)))
})

test_that("fit_linear_gwpr: within smoke test completes without error", {
  ctx <- make_context(model = "within", effect = "individual",
                      kernel = "gaussian", adaptive = FALSE)
  expect_no_error({
    res <- fit_linear_gwpr(ctx, bandwidth = 2)
  })
  expect_equal(length(res$local_results), 4L)
})

test_that("fit_linear_gwpr: random smoke test completes without error", {
  ctx <- make_context(model = "random", effect = "individual",
                      kernel = "bisquare", adaptive = FALSE)
  expect_no_error({
    res <- fit_linear_gwpr(ctx, bandwidth = 2)
  })
  expect_equal(length(res$local_results), 4L)
})

# ---------------------------------------------------------------------------
# fit_linear_gwpr — effect parameter parsing
# ---------------------------------------------------------------------------

test_that("fit_linear_gwpr: effect=individual runs without error (pooling)", {
  ctx <- make_context(model = "pooling", effect = "individual", adaptive = FALSE)
  expect_no_error(fit_linear_gwpr(ctx, bandwidth = 2))
})

test_that("fit_linear_gwpr: effect=time runs without error (within)", {
  ctx <- make_context(model = "within", effect = "time", adaptive = FALSE)
  # within + time effect uses plm effect="time"
  expect_no_error(fit_linear_gwpr(ctx, bandwidth = 2))
})

test_that("fit_linear_gwpr: effect=two-way runs without error (within)", {
  ctx <- make_context(model = "within", effect = "two-way", adaptive = FALSE)
  expect_no_error(fit_linear_gwpr(ctx, bandwidth = 2))
})

test_that("fit_linear_gwpr: effect=nested stops with informative error", {
  ctx <- make_context(model = "within", effect = "nested", adaptive = FALSE)
  expect_error(fit_linear_gwpr(ctx, bandwidth = 2), regexp = "nested")
})

# ---------------------------------------------------------------------------
# fit_linear_gwpr — failure robustness
# ---------------------------------------------------------------------------

test_that("fit_linear_gwpr: local model failure does not break overall structure", {
  # Use a formula with a non-existent predictor to force local model failures
  pd  <- make_panel_data()
  ctx <- new_gwpr_context(
    formula    = y ~ nonexistent,   # will fail every local model
    family     = "gaussian",
    model      = "pooling",
    effect     = "individual",
    id         = "id",
    time       = "time",
    kernel     = "bisquare",
    adaptive   = TRUE,
    threshold  = 0.5,
    workers    = 1L,
    panel_data = pd,
    coords     = make_coords_4()
  )

  res <- fit_linear_gwpr(ctx, bandwidth = 3L)

  # Structure intact
  expect_true("local_results" %in% names(res))
  expect_equal(length(res$local_results), 4L)
  # All failed
  statuses <- vapply(res$local_results, `[[`, character(1L), "status")
  expect_true(all(statuses == "failed"))
  # Metadata records failure count
  expect_equal(res$metadata$n_failed, 4L)
})

# ---------------------------------------------------------------------------
# fit_linear_gwpr — reproducibility with fixed weights
# ---------------------------------------------------------------------------

test_that("fit_linear_gwpr: identical weights give identical results", {
  ctx <- make_context(model = "pooling", effect = "individual",
                      kernel = "bisquare", adaptive = TRUE)

  res1 <- fit_linear_gwpr(ctx, bandwidth = 3L)
  res2 <- fit_linear_gwpr(ctx, bandwidth = 3L)

  # Coefficients of first focus should be identical
  expect_equal(
    res1$local_results[[1]]$coefficients,
    res2$local_results[[1]]$coefficients
  )
  expect_equal(res1$predictions, res2$predictions)
})

# ---------------------------------------------------------------------------
# fit_linear_gwpr — pre-supplied weights_context
# ---------------------------------------------------------------------------

test_that("fit_linear_gwpr: accepts pre-built weights_context", {
  ctx    <- make_context(model = "pooling", effect = "individual",
                         kernel = "bisquare", adaptive = TRUE)
  wctx   <- build_distance_context(
    coords  = ctx$coords,
    ids     = rownames(ctx$coords),
    longlat = FALSE,
    cache   = TRUE
  )
  expect_no_error(fit_linear_gwpr(ctx, bandwidth = 3L, weights_context = wctx))
})
