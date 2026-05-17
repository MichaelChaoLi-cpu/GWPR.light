## tests/testthat/test-bandwidth-sgd.R
## Tests for the SGD Bandwidth Search module (R/bandwidth_sgd.R)
##
## All tests use small simulated data and workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helper: minimal mock scorer
# ---------------------------------------------------------------------------
# Returns a score that is the absolute distance from a "best" bandwidth,
# so we can predict which direction gradient descent should move.

make_sgd_mock_scorer <- function(best_bw = 10, fail_bw = NULL) {
  function(context, bandwidth) {
    if (!is.null(fail_bw) && isTRUE(all.equal(bandwidth, fail_bw))) {
      stop("Intentional scorer failure for SGD testing.")
    }
    score <- abs(bandwidth - best_bw)^2 + 0.01
    list(
      score                 = score,
      criterion             = "MSE",
      status                = "ok",
      n_local_models        = 4L,
      n_failed_local_models = 0L,
      metrics               = list(
        R2   = 1 - min(score, 1),
        MSE  = score,
        RMSE = sqrt(score),
        MAE  = score * 0.9
      )
    )
  }
}

# Monotonically increasing scorer (score always increases each call)
make_worsening_scorer <- function() {
  call_count <- 0L
  function(context, bandwidth) {
    call_count <<- call_count + 1L
    score <- call_count * 1.5
    list(
      score                 = score,
      criterion             = "MSE",
      status                = "ok",
      n_local_models        = 2L,
      n_failed_local_models = 0L,
      metrics               = list(
        R2 = 0.5, MSE = score, RMSE = sqrt(score), MAE = score * 0.9
      )
    )
  }
}

# Minimal gwpr_context with an id_map for minibatch sampling
make_sgd_context <- function(adaptive = FALSE, n_units = 6L) {
  id_map <- stats::setNames(seq_len(n_units), paste0("u", seq_len(n_units)))
  new_gwpr_context(
    formula  = y ~ x1,
    family   = "gaussian",
    model    = "pooling",
    effect   = "individual",
    id       = "id",
    time     = "time",
    kernel   = "bisquare",
    adaptive = adaptive,
    threshold = 0.5,
    workers  = 1L,
    id_map   = id_map
  )
}

# ===========================================================================
# initialize_bandwidth
# ===========================================================================

test_that("initialize_bandwidth errors when lower is missing from control", {
  expect_error(
    initialize_bandwidth(list(upper = 100), adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("initialize_bandwidth errors when upper is missing from control", {
  expect_error(
    initialize_bandwidth(list(lower = 1), adaptive = FALSE),
    regexp = "upper"
  )
})

test_that("initialize_bandwidth errors when lower is NULL", {
  expect_error(
    initialize_bandwidth(list(lower = NULL, upper = 100), adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("initialize_bandwidth errors when upper is NULL", {
  expect_error(
    initialize_bandwidth(list(lower = 1, upper = NULL), adaptive = FALSE),
    regexp = "upper"
  )
})

test_that("initialize_bandwidth errors when lower >= upper", {
  expect_error(
    initialize_bandwidth(list(lower = 10, upper = 5), adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("initialize_bandwidth returns value in [lower, upper] for fixed", {
  bw <- initialize_bandwidth(list(lower = 2, upper = 50), adaptive = FALSE,
                              seed = 42L)
  expect_true(is.numeric(bw))
  expect_true(length(bw) == 1L)
  expect_true(bw >= 2)
  expect_true(bw <= 50)
})

test_that("initialize_bandwidth returns positive value for fixed bandwidth", {
  for (s in 1:5) {
    bw <- initialize_bandwidth(list(lower = 0.5, upper = 200), adaptive = FALSE,
                                seed = s)
    expect_true(bw > 0)
  }
})

test_that("initialize_bandwidth returns integer in range for adaptive", {
  bw <- initialize_bandwidth(list(lower = 3, upper = 20), adaptive = TRUE,
                              seed = 7L)
  expect_true(is.integer(bw) || bw == as.integer(bw))
  expect_true(bw >= 3L)
  expect_true(bw <= 20L)
})

test_that("initialize_bandwidth is reproducible with fixed seed", {
  bw1 <- initialize_bandwidth(list(lower = 1, upper = 100), adaptive = FALSE,
                               seed = 99L)
  bw2 <- initialize_bandwidth(list(lower = 1, upper = 100), adaptive = FALSE,
                               seed = 99L)
  expect_equal(bw1, bw2)
})

# ===========================================================================
# sample_minibatch
# ===========================================================================

test_that("sample_minibatch returns a subset of unit IDs", {
  ctx  <- make_sgd_context(n_units = 10L)
  mini <- sample_minibatch(ctx, batch_fraction = 0.5, seed = 1L)
  all_ids <- names(ctx$id_map)
  expect_true(all(mini %in% all_ids))
})

test_that("sample_minibatch returns correct approximate size", {
  ctx  <- make_sgd_context(n_units = 10L)
  mini <- sample_minibatch(ctx, batch_fraction = 0.6, seed = 2L)
  expect_equal(length(mini), as.integer(ceiling(10L * 0.6)))
})

test_that("sample_minibatch is reproducible with fixed seed", {
  ctx  <- make_sgd_context(n_units = 8L)
  m1   <- sample_minibatch(ctx, batch_fraction = 0.5, seed = 42L)
  m2   <- sample_minibatch(ctx, batch_fraction = 0.5, seed = 42L)
  expect_equal(m1, m2)
})

test_that("sample_minibatch differs across seeds", {
  ctx <- make_sgd_context(n_units = 20L)
  m1  <- sample_minibatch(ctx, batch_fraction = 0.5, seed = 1L)
  m2  <- sample_minibatch(ctx, batch_fraction = 0.5, seed = 2L)
  # With 10 out of 20 units, it is extremely unlikely they match
  expect_false(identical(sort(m1), sort(m2)))
})

test_that("sample_minibatch with batch_fraction = 1 returns all units", {
  ctx  <- make_sgd_context(n_units = 5L)
  mini <- sample_minibatch(ctx, batch_fraction = 1.0, seed = 10L)
  expect_equal(sort(mini), sort(names(ctx$id_map)))
})

test_that("sample_minibatch errors on invalid batch_fraction", {
  ctx <- make_sgd_context()
  expect_error(sample_minibatch(ctx, batch_fraction = 0),   regexp = "batch_fraction")
  expect_error(sample_minibatch(ctx, batch_fraction = 1.5), regexp = "batch_fraction")
  expect_error(sample_minibatch(ctx, batch_fraction = -0.5), regexp = "batch_fraction")
})

# ===========================================================================
# estimate_bandwidth_gradient
# ===========================================================================

test_that("estimate_bandwidth_gradient returns a numeric scalar", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer(best_bw = 10)
  grad   <- estimate_bandwidth_gradient(ctx, bandwidth = 5, delta = 0.5,
                                         minibatch = character(), scorer = scorer)
  expect_true(is.numeric(grad))
  expect_length(grad, 1L)
})

test_that("estimate_bandwidth_gradient is negative when bw < best_bw", {
  # f(b) = (b - best_bw)^2 + 0.01; gradient at b < best_bw is negative
  # → increasing b decreases score → f(b + delta) < f(b - delta) → grad < 0
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer(best_bw = 10)
  grad   <- estimate_bandwidth_gradient(ctx, bandwidth = 5, delta = 1,
                                         minibatch = character(), scorer = scorer)
  expect_true(grad < 0)
})

test_that("estimate_bandwidth_gradient is positive when bw > best_bw", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer(best_bw = 10)
  grad   <- estimate_bandwidth_gradient(ctx, bandwidth = 20, delta = 1,
                                         minibatch = character(), scorer = scorer)
  expect_true(grad > 0)
})

test_that("estimate_bandwidth_gradient returns 0 on scorer failure", {
  ctx    <- make_sgd_context()
  bad_scorer <- function(context, bw) stop("always fails")
  grad <- estimate_bandwidth_gradient(ctx, bandwidth = 5, delta = 1,
                                       minibatch = character(),
                                       scorer = bad_scorer)
  expect_equal(grad, 0)
})

# ===========================================================================
# update_bandwidth_sgd
# ===========================================================================

test_that("update_bandwidth_sgd fixed: result is strictly positive", {
  for (bw in c(0.1, 1, 10, 100)) {
    new_bw <- update_bandwidth_sgd(bw, gradient = 5, learning_rate = 0.5,
                                   adaptive = FALSE)
    expect_true(new_bw > 0, label = paste("bw =", bw))
  }
})

test_that("update_bandwidth_sgd fixed: moves in the right direction", {
  # Positive gradient → bandwidth should decrease
  bw     <- 10
  new_bw <- update_bandwidth_sgd(bw, gradient = 2, learning_rate = 0.5,
                                 adaptive = FALSE)
  expect_true(new_bw < bw)

  # Negative gradient → bandwidth should increase
  new_bw2 <- update_bandwidth_sgd(bw, gradient = -2, learning_rate = 0.5,
                                  adaptive = FALSE)
  expect_true(new_bw2 > bw)
})

test_that("update_bandwidth_sgd fixed: respects lower and upper bounds", {
  # Force a very large step that would violate bounds
  new_bw <- update_bandwidth_sgd(10, gradient = 1e8, learning_rate = 1,
                                 adaptive = FALSE, lower = 5, upper = 20)
  expect_true(new_bw >= 5)
  expect_true(new_bw <= 20)
})

test_that("update_bandwidth_sgd adaptive: result is a positive integer", {
  new_bw <- update_bandwidth_sgd(8L, gradient = -1.5, learning_rate = 2,
                                 adaptive = TRUE, lower = 3L, upper = 20L)
  expect_true(new_bw == as.integer(new_bw))
  expect_true(new_bw >= 1L)
})

test_that("update_bandwidth_sgd adaptive: clips to [lower, upper]", {
  # Push well below lower
  new_bw <- update_bandwidth_sgd(4L, gradient = -1e6, learning_rate = 1,
                                 adaptive = TRUE, lower = 3L, upper = 20L)
  expect_true(new_bw >= 3L)

  # Push well above upper
  new_bw2 <- update_bandwidth_sgd(15L, gradient = 1e6, learning_rate = 1,
                                  adaptive = TRUE, lower = 3L, upper = 20L)
  expect_true(new_bw2 <= 20L)
})

# ===========================================================================
# check_early_stopping
# ===========================================================================

test_that("check_early_stopping returns FALSE when patience = 0", {
  h <- data.frame(score = c(1, 2, 3, 4, 5))
  expect_false(check_early_stopping(h, patience = 0L))
})

test_that("check_early_stopping returns FALSE with improving scores", {
  h <- data.frame(score = c(1.0, 0.9, 0.8, 0.7, 0.6))
  expect_false(check_early_stopping(h, patience = 2L))
})

test_that("check_early_stopping returns TRUE when score stops improving", {
  # Score improved at epoch 1-3 but stagnated for 2 consecutive epochs
  h <- data.frame(score = c(1.0, 0.9, 0.85, 0.85, 0.86))
  expect_true(check_early_stopping(h, patience = 2L))
})

test_that("check_early_stopping returns FALSE when history is too short", {
  h <- data.frame(score = c(0.5, 0.6))
  expect_false(check_early_stopping(h, patience = 3L))
})

test_that("check_early_stopping counter tracks non-improving epochs", {
  # With scores that are always worsening and patience=1
  h <- data.frame(score = c(0.5, 0.6, 0.7))
  expect_true(check_early_stopping(h, patience = 1L))
})

# ===========================================================================
# search_bandwidth_sgd
# ===========================================================================

test_that("search_bandwidth_sgd returns a gwpr_bandwidth object", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer(best_bw = 10)
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 5L,
                   learning_rate = 0.01),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  expect_s3_class(result, "gwpr_bandwidth")
})

test_that("search_bandwidth_sgd method field is 'sgd'", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 4L),
    scorer  = scorer,
    workers = 1L,
    seed    = 2L
  )
  expect_equal(result$method, "sgd")
})

test_that("search_bandwidth_sgd default epoch is 10", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 100),   # epoch not specified
    scorer  = scorer,
    workers = 1L,
    seed    = 3L
  )
  expect_equal(nrow(result$history), 10L)
})

test_that("search_bandwidth_sgd respects explicit epoch parameter", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 7L),
    scorer  = scorer,
    workers = 1L,
    seed    = 4L
  )
  expect_equal(nrow(result$history), 7L)
})

test_that("search_bandwidth_sgd history epoch column is sequential", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 6L),
    scorer  = scorer,
    workers = 1L,
    seed    = 5L
  )
  expect_equal(result$history$epoch, seq_len(6L))
})

test_that("search_bandwidth_sgd history has all required columns", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 3L),
    scorer  = scorer,
    workers = 1L,
    seed    = 6L
  )
  required_cols <- c(
    "epoch", "bandwidth_before", "bandwidth_after", "score",
    "gradient", "learning_rate", "delta", "batch_size",
    "early_stop_counter", "status", "elapsed_time",
    "r2", "mse", "rmse", "mae",
    "log_loss", "accuracy", "precision", "recall", "f1_score"
  )
  expect_true(all(required_cols %in% names(result$history)))
})

test_that("search_bandwidth_sgd fixed bandwidth never goes below 0", {
  ctx    <- make_sgd_context(adaptive = FALSE)
  scorer <- make_sgd_mock_scorer(best_bw = 1)
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 0.1, upper = 100, epoch = 10L,
                   learning_rate = 0.1),
    scorer  = scorer,
    workers = 1L,
    seed    = 7L
  )
  expect_true(all(result$history$bandwidth_after > 0, na.rm = TRUE))
})

test_that("search_bandwidth_sgd adaptive bandwidth stays as integers", {
  ctx    <- make_sgd_context(adaptive = TRUE)
  scorer <- make_sgd_mock_scorer(best_bw = 8)
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 3L, upper = 20L, epoch = 8L,
                   learning_rate = 2),
    scorer  = scorer,
    workers = 1L,
    seed    = 8L
  )
  bw_after <- result$history$bandwidth_after
  expect_true(all(bw_after == as.integer(bw_after), na.rm = TRUE))
  expect_true(all(bw_after >= 3L, na.rm = TRUE))
  expect_true(all(bw_after <= 20L, na.rm = TRUE))
})

test_that("search_bandwidth_sgd patience=0 does not early-stop", {
  ctx    <- make_sgd_context()
  scorer <- make_worsening_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 100, epoch = 6L,
                   early_stopping_patience = 0L),
    scorer  = scorer,
    workers = 1L,
    seed    = 9L
  )
  # All 6 epochs should run and have status "ok"
  expect_equal(nrow(result$history), 6L)
  ok_rows <- result$history[result$history$status == "ok", , drop = FALSE]
  expect_equal(nrow(ok_rows), 6L)
})

test_that("search_bandwidth_sgd early_stop_counter increases when score worsens", {
  ctx    <- make_sgd_context()
  scorer <- make_worsening_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 100, epoch = 5L,
                   early_stopping_patience = 0L,  # don't actually stop
                   learning_rate = 0.0001),
    scorer  = scorer,
    workers = 1L,
    seed    = 10L
  )
  # From epoch 2 onward, score always increases → counter should be > 0
  counters <- result$history$early_stop_counter[result$history$epoch >= 2L]
  expect_true(any(counters > 0L, na.rm = TRUE))
})

test_that("search_bandwidth_sgd early stopping halts after patience exceeded", {
  ctx    <- make_sgd_context()
  scorer <- make_worsening_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 100, epoch = 10L,
                   early_stopping_patience = 2L,
                   learning_rate = 0.0001),
    scorer  = scorer,
    workers = 1L,
    seed    = 11L
  )
  # After patience = 2 non-improving epochs we should have some early_stopped rows
  stopped_rows <- result$history[result$history$status == "early_stopped", , drop = FALSE]
  expect_true(nrow(stopped_rows) > 0L)
})

test_that("search_bandwidth_sgd records seed in returned object", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 3L),
    scorer  = scorer,
    workers = 1L,
    seed    = 77L
  )
  expect_equal(result$seed, 77L)
})

test_that("search_bandwidth_sgd is reproducible with fixed seed", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()

  res1 <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 5L, learning_rate = 0.01),
    scorer  = scorer,
    workers = 1L,
    seed    = 123L
  )
  res2 <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 5L, learning_rate = 0.01),
    scorer  = scorer,
    workers = 1L,
    seed    = 123L
  )
  expect_equal(res1$history$bandwidth_after, res2$history$bandwidth_after)
  expect_equal(res1$best_bandwidth, res2$best_bandwidth)
})

test_that("search_bandwidth_sgd best_bandwidth is positive for fixed", {
  ctx    <- make_sgd_context(adaptive = FALSE)
  scorer <- make_sgd_mock_scorer(best_bw = 15)
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 100, epoch = 10L,
                   learning_rate = 0.001),
    scorer  = scorer,
    workers = 1L,
    seed    = 14L
  )
  expect_true(is.numeric(result$best_bandwidth))
  expect_true(result$best_bandwidth > 0)
})

test_that("search_bandwidth_sgd best_bandwidth is integer for adaptive", {
  ctx    <- make_sgd_context(adaptive = TRUE)
  scorer <- make_sgd_mock_scorer(best_bw = 10)
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 3L, upper = 25L, epoch = 8L,
                   learning_rate = 1),
    scorer  = scorer,
    workers = 1L,
    seed    = 15L
  )
  expect_true(result$best_bandwidth == as.integer(result$best_bandwidth))
})

test_that("search_bandwidth_sgd elapsed_time is non-negative", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 3L),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  expect_true(result$elapsed_time >= 0)
})

test_that("search_bandwidth_sgd convergence_info contains key fields", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  result <- search_bandwidth_sgd(
    ctx,
    control = list(lower = 1, upper = 50, epoch = 4L,
                   learning_rate = 0.005),
    scorer  = scorer,
    workers = 1L,
    seed    = 2L
  )
  expect_true("epoch"         %in% names(result$convergence_info))
  expect_true("learning_rate" %in% names(result$convergence_info))
  expect_true("patience"      %in% names(result$convergence_info))
  expect_true("delta"         %in% names(result$convergence_info))
})

test_that("search_bandwidth_sgd errors when lower is missing", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  expect_error(
    search_bandwidth_sgd(
      ctx,
      control = list(upper = 50, epoch = 3L),
      scorer  = scorer,
      workers = 1L,
      seed    = 1L
    ),
    regexp = "lower"
  )
})

test_that("search_bandwidth_sgd errors when upper is missing", {
  ctx    <- make_sgd_context()
  scorer <- make_sgd_mock_scorer()
  expect_error(
    search_bandwidth_sgd(
      ctx,
      control = list(lower = 1, epoch = 3L),
      scorer  = scorer,
      workers = 1L,
      seed    = 1L
    ),
    regexp = "upper"
  )
})
