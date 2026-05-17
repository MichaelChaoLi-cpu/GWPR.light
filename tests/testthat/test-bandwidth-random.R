## tests/testthat/test-bandwidth-random.R
## Tests for the Random Bandwidth Optimizer module (R/bandwidth_random.R)
##
## All tests use small simulated data and workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helper: minimal mock scorer
# ---------------------------------------------------------------------------
# Returns a valid scorer result.  Score is deterministic (distance from a
# best bandwidth) so we can predict which candidate wins.

make_random_mock_scorer <- function(best_bw = 5, fail_bw = NULL) {
  function(context, bandwidth) {
    if (!is.null(fail_bw) && isTRUE(all.equal(bandwidth, fail_bw))) {
      stop("Intentional scorer failure for testing.")
    }
    score <- abs(bandwidth - best_bw) + 0.01
    list(
      score                 = score,
      criterion             = "MSE",
      status                = "ok",
      n_local_models        = 4L,
      n_failed_local_models = 0L,
      metrics               = list(
        R2   = 0.8,
        MSE  = score,
        RMSE = sqrt(score),
        MAE  = score * 0.9
      )
    )
  }
}

# Minimal gwpr_context suitable for random search tests
make_random_context <- function(adaptive = FALSE) {
  new_gwpr_context(
    formula   = y ~ x1,
    family    = "gaussian",
    model     = "pooling",
    effect    = "individual",
    id        = "id",
    time      = "time",
    kernel    = "bisquare",
    adaptive  = adaptive,
    threshold = 0.5,
    workers   = 1L,
    seed      = NULL
  )
}

# ===========================================================================
# sample_random_bandwidths
# ===========================================================================

test_that("sample_random_bandwidths errors when lower is missing", {
  expect_error(
    sample_random_bandwidths(upper = 10, n_samples = 5, adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("sample_random_bandwidths errors when upper is missing", {
  expect_error(
    sample_random_bandwidths(lower = 1, n_samples = 5, adaptive = FALSE),
    regexp = "upper"
  )
})

test_that("sample_random_bandwidths errors when lower is NULL", {
  expect_error(
    sample_random_bandwidths(lower = NULL, upper = 10, n_samples = 5,
                             adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("sample_random_bandwidths errors when upper is NULL", {
  expect_error(
    sample_random_bandwidths(lower = 1, upper = NULL, n_samples = 5,
                             adaptive = FALSE),
    regexp = "upper"
  )
})

test_that("sample_random_bandwidths errors when lower >= upper", {
  expect_error(
    sample_random_bandwidths(lower = 10, upper = 5, n_samples = 5,
                             adaptive = FALSE),
    regexp = "lower"
  )
  expect_error(
    sample_random_bandwidths(lower = 5, upper = 5, n_samples = 5,
                             adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("sample_random_bandwidths returns n_samples values for adaptive = FALSE", {
  cands <- sample_random_bandwidths(lower = 1, upper = 10, n_samples = 20L,
                                    adaptive = FALSE, seed = 99L)
  expect_length(cands, 20L)
})

test_that("sample_random_bandwidths returns values within [lower, upper]", {
  cands <- sample_random_bandwidths(lower = 2, upper = 8, n_samples = 50L,
                                    adaptive = FALSE, seed = 7L)
  expect_true(all(cands >= 2))
  expect_true(all(cands <= 8))
})

test_that("sample_random_bandwidths is reproducible with fixed seed", {
  cands1 <- sample_random_bandwidths(lower = 1, upper = 100, n_samples = 10L,
                                     adaptive = FALSE, seed = 42L)
  cands2 <- sample_random_bandwidths(lower = 1, upper = 100, n_samples = 10L,
                                     adaptive = FALSE, seed = 42L)
  expect_equal(cands1, cands2)
})

test_that("sample_random_bandwidths differs across seeds", {
  cands1 <- sample_random_bandwidths(lower = 1, upper = 100, n_samples = 10L,
                                     adaptive = FALSE, seed = 1L)
  cands2 <- sample_random_bandwidths(lower = 1, upper = 100, n_samples = 10L,
                                     adaptive = FALSE, seed = 2L)
  expect_false(identical(cands1, cands2))
})

test_that("sample_random_bandwidths adaptive = TRUE returns integers", {
  cands <- sample_random_bandwidths(lower = 2, upper = 20, n_samples = 30L,
                                    adaptive = TRUE, seed = 5L)
  expect_type(cands, "integer")
  expect_true(all(cands >= 1L))
})

test_that("sample_random_bandwidths adaptive = TRUE removes duplicates", {
  cands <- sample_random_bandwidths(lower = 1, upper = 3, n_samples = 100L,
                                    adaptive = TRUE, seed = 10L)
  expect_true(!any(duplicated(cands)))
  # With range 1-3 and rounding there can be at most 3 unique integers
  expect_true(length(cands) <= 3L)
})

test_that("sample_random_bandwidths adaptive = TRUE candidates are valid positive integers", {
  cands <- sample_random_bandwidths(lower = 1, upper = 10, n_samples = 15L,
                                    adaptive = TRUE, seed = 21L)
  expect_true(all(cands == as.integer(cands)))
  expect_true(all(cands >= 1L))
})

# ===========================================================================
# search_bandwidth_random
# ===========================================================================

test_that("search_bandwidth_random errors when lower is missing from control", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  expect_error(
    search_bandwidth_random(ctx, control = list(upper = 10, n_samples = 5L),
                            scorer = scorer),
    regexp = "lower"
  )
})

test_that("search_bandwidth_random errors when upper is missing from control", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  expect_error(
    search_bandwidth_random(ctx, control = list(lower = 1, n_samples = 5L),
                            scorer = scorer),
    regexp = "upper"
  )
})

test_that("search_bandwidth_random errors when lower is NULL in control", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  expect_error(
    search_bandwidth_random(ctx,
                            control = list(lower = NULL, upper = 10,
                                           n_samples = 5L),
                            scorer = scorer),
    regexp = "lower"
  )
})

test_that("search_bandwidth_random returns a gwpr_bandwidth object", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 10L),
    scorer  = scorer,
    workers = 1L,
    seed    = 42L
  )
  expect_s3_class(result, "gwpr_bandwidth")
})

test_that("search_bandwidth_random method field is 'random'", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 8L),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  expect_equal(result$method, "random")
  expect_true(all(result$history$method == "random"))
})

test_that("search_bandwidth_random history has correct number of rows (fixed)", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 12L),
    scorer  = scorer,
    workers = 1L,
    seed    = 3L
  )
  expect_equal(nrow(result$history), 12L)
})

test_that("search_bandwidth_random best_bandwidth comes from candidate set", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer(best_bw = 5)
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 15L),
    scorer  = scorer,
    workers = 1L,
    seed    = 7L
  )
  # best_bandwidth must appear in the candidates evaluated
  expect_true(result$best_bandwidth %in% result$history$bandwidth)
})

test_that("search_bandwidth_random best_bandwidth has the minimum score", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer(best_bw = 5)
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 20L),
    scorer  = scorer,
    workers = 1L,
    seed    = 11L
  )
  ok_hist <- result$history[result$history$status == "ok", , drop = FALSE]
  expect_equal(result$best_score, min(ok_hist$score, na.rm = TRUE))
})

test_that("search_bandwidth_random is reproducible with fixed seed", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()

  res1 <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 20, n_samples = 10L),
    scorer  = scorer,
    workers = 1L,
    seed    = 99L
  )
  res2 <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 20, n_samples = 10L),
    scorer  = scorer,
    workers = 1L,
    seed    = 99L
  )
  expect_equal(res1$history$bandwidth, res2$history$bandwidth)
  expect_equal(res1$best_bandwidth, res2$best_bandwidth)
})

test_that("search_bandwidth_random records seed in returned object", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 5L),
    scorer  = scorer,
    workers = 1L,
    seed    = 55L
  )
  expect_equal(result$seed, 55L)
})

test_that("search_bandwidth_random history includes all required columns", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 5L),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  required_cols <- c(
    "candidate_id", "bandwidth", "bandwidth_type", "method", "criterion",
    "score", "rank", "status", "error_message", "warning_message",
    "elapsed_time", "n_local_models", "n_failed_local_models",
    "r2", "mse", "rmse", "mae",
    "log_loss", "accuracy", "precision", "recall", "f1_score"
  )
  expect_true(all(required_cols %in% names(result$history)))
})

test_that("search_bandwidth_random retains failed candidates in history", {
  ctx         <- make_random_context()
  # We need to know what bandwidth will be drawn to create the failing scorer.
  # Use a seed to get a fixed candidate set first, then fail one of them.
  cands <- sample_random_bandwidths(lower = 1, upper = 10, n_samples = 6L,
                                    adaptive = FALSE, seed = 77L)
  fail_bw <- cands[1L]
  scorer  <- make_random_mock_scorer(best_bw = 5, fail_bw = fail_bw)

  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 6L),
    scorer  = scorer,
    workers = 1L,
    seed    = 77L   # same seed => same candidate set
  )
  # All 6 candidates must be present
  expect_equal(nrow(result$history), 6L)
  # At least one failed
  expect_true(any(result$history$status == "failed"))
  # The failed bandwidth is recorded
  failed_rows <- result$history[result$history$status == "failed", , drop = FALSE]
  expect_true(any(abs(failed_rows$bandwidth - fail_bw) < 1e-9))
})

test_that("search_bandwidth_random best_bandwidth comes from ok candidates only", {
  ctx <- make_random_context()
  # Draw fixed candidates, make the one closest to 5 fail
  cands <- sample_random_bandwidths(lower = 1, upper = 10, n_samples = 10L,
                                    adaptive = FALSE, seed = 33L)
  # Find which candidate is closest to 5
  closest_idx <- which.min(abs(cands - 5))
  fail_bw     <- cands[closest_idx]
  scorer      <- make_random_mock_scorer(best_bw = 5, fail_bw = fail_bw)

  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 10L),
    scorer  = scorer,
    workers = 1L,
    seed    = 33L
  )
  # best_bandwidth must be among ok candidates
  ok_bws <- result$history$bandwidth[result$history$status == "ok"]
  expect_true(result$best_bandwidth %in% ok_bws)
  # It must NOT be the failed bandwidth
  expect_false(isTRUE(all.equal(result$best_bandwidth, fail_bw)))
})

test_that("search_bandwidth_random adaptive produces integer candidates", {
  ctx    <- make_random_context(adaptive = TRUE)
  scorer <- make_random_mock_scorer(best_bw = 5)
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 2, upper = 20, n_samples = 15L),
    scorer  = scorer,
    workers = 1L,
    seed    = 8L
  )
  expect_true(
    all(result$history$bandwidth == as.integer(result$history$bandwidth))
  )
})

test_that("search_bandwidth_random uses default n_samples when not provided", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 100),   # n_samples omitted
    scorer  = scorer,
    workers = 1L,
    seed    = 2L
  )
  # Default is 50; result should have 50 rows
  expect_equal(nrow(result$history), 50L)
})

test_that("search_bandwidth_random elapsed_time is non-negative", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 5L),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  expect_true(result$elapsed_time >= 0)
})

test_that("search_bandwidth_random convergence_info contains n_candidates", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer()
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 7L),
    scorer  = scorer,
    workers = 1L,
    seed    = 1L
  )
  expect_equal(result$convergence_info$n_candidates, 7L)
})

test_that("search_bandwidth_random best_bandwidth has rank 1 in history", {
  ctx    <- make_random_context()
  scorer <- make_random_mock_scorer(best_bw = 5)
  result <- search_bandwidth_random(
    ctx,
    control = list(lower = 1, upper = 10, n_samples = 15L),
    scorer  = scorer,
    workers = 1L,
    seed    = 17L
  )
  best_row <- result$history[
    abs(result$history$bandwidth - result$best_bandwidth) < 1e-9 &
      result$history$status == "ok",
    ,
    drop = FALSE
  ]
  expect_equal(best_row$rank[1L], 1L)
})
