## tests/testthat/test-bandwidth-grid.R
## Tests for the Grid Search Bandwidth module (R/bandwidth_grid.R)
##
## All tests use small simulated data and workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helper: minimal mock scorer
# ---------------------------------------------------------------------------
# Returns a valid scorer result.  Score is a simple deterministic function of
# bandwidth so we can predict which candidate wins.

make_mock_scorer <- function(best_bw = 3, fail_bw = NULL) {
  function(context, bandwidth) {
    if (!is.null(fail_bw) && bandwidth == fail_bw) {
      stop("Intentional scorer failure for testing.")
    }
    # Score: distance from best_bw; lower is better
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

# Minimal gwpr_context suitable for grid search tests
make_grid_context <- function(adaptive = FALSE) {
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
# make_grid_candidates
# ===========================================================================

test_that("make_grid_candidates errors when lower is missing", {
  expect_error(
    make_grid_candidates(upper = 5, step = 1, adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("make_grid_candidates errors when upper is missing", {
  expect_error(
    make_grid_candidates(lower = 1, step = 1, adaptive = FALSE),
    regexp = "upper"
  )
})

test_that("make_grid_candidates errors when step is missing", {
  expect_error(
    make_grid_candidates(lower = 1, upper = 5, adaptive = FALSE),
    regexp = "step"
  )
})

test_that("make_grid_candidates errors when lower >= upper", {
  expect_error(
    make_grid_candidates(lower = 5, upper = 1, step = 1, adaptive = FALSE),
    regexp = "lower"
  )
  expect_error(
    make_grid_candidates(lower = 5, upper = 5, step = 1, adaptive = FALSE),
    regexp = "lower"
  )
})

test_that("make_grid_candidates errors when step <= 0", {
  expect_error(
    make_grid_candidates(lower = 1, upper = 5, step = -1, adaptive = FALSE),
    regexp = "step"
  )
  expect_error(
    make_grid_candidates(lower = 1, upper = 5, step = 0, adaptive = FALSE),
    regexp = "step"
  )
})

test_that("make_grid_candidates returns correct number of candidates (fixed)", {
  cands <- make_grid_candidates(lower = 1, upper = 5, step = 1, adaptive = FALSE)
  expect_length(cands, 5L)
  expect_equal(cands, c(1, 2, 3, 4, 5))
})

test_that("make_grid_candidates works with fractional step", {
  cands <- make_grid_candidates(lower = 0, upper = 1, step = 0.25, adaptive = FALSE)
  expect_length(cands, 5L)
})

test_that("make_grid_candidates rounds and deduplicates for adaptive = TRUE", {
  # step = 0.4 from 1 to 3 => seq: 1.0, 1.4, 1.8, 2.2, 2.6, 3.0
  # rounded: 1, 1, 2, 2, 3, 3 => unique: 1, 2, 3
  cands <- make_grid_candidates(lower = 1, upper = 3, step = 0.4, adaptive = TRUE)
  expect_true(all(cands == as.integer(cands)))  # all integers
  expect_true(!any(duplicated(cands)))           # no duplicates
  expect_true(all(cands >= 1L))                  # all positive
})

test_that("make_grid_candidates returns integers for adaptive = TRUE", {
  cands <- make_grid_candidates(lower = 2, upper = 8, step = 2, adaptive = TRUE)
  expect_type(cands, "integer")
})

# ===========================================================================
# score_bandwidth_candidate
# ===========================================================================

test_that("score_bandwidth_candidate returns ok status for valid scorer", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 3)
  result <- score_bandwidth_candidate(ctx, bandwidth = 3, scorer = scorer)

  expect_equal(result$status, "ok")
  expect_equal(result$bandwidth, 3)
  expect_true(is.numeric(result$score))
  expect_false(is.na(result$score))
  expect_true(is.numeric(result$elapsed_time))
})

test_that("score_bandwidth_candidate returns failed status on scorer error", {
  ctx    <- make_grid_context()
  bad_scorer <- function(context, bandwidth) stop("deliberate error")
  result <- score_bandwidth_candidate(ctx, bandwidth = 2, scorer = bad_scorer)

  expect_equal(result$status, "failed")
  expect_true(!is.na(result$error_message))
  expect_true(is.na(result$score))
})

test_that("score_bandwidth_candidate populates linear metrics", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 3)
  result <- score_bandwidth_candidate(ctx, bandwidth = 3, scorer = scorer)

  expect_false(is.na(result$r2))
  expect_false(is.na(result$mse))
  expect_false(is.na(result$rmse))
  expect_false(is.na(result$mae))
})

test_that("score_bandwidth_candidate fills logistic metrics with NA for gaussian scorer", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 3)
  result <- score_bandwidth_candidate(ctx, bandwidth = 3, scorer = scorer)

  expect_true(is.na(result$log_loss))
  expect_true(is.na(result$accuracy))
  expect_true(is.na(result$precision))
  expect_true(is.na(result$recall))
  expect_true(is.na(result$f1_score))
})

test_that("score_bandwidth_candidate records elapsed_time", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- score_bandwidth_candidate(ctx, bandwidth = 2, scorer = scorer)
  expect_true(is.numeric(result$elapsed_time) && result$elapsed_time >= 0)
})

# ===========================================================================
# rank_grid_history
# ===========================================================================

test_that("rank_grid_history ranks ok candidates by ascending score", {
  h <- data.frame(
    candidate_id = 1:3,
    score        = c(0.5, 0.3, 0.7),
    status       = c("ok", "ok", "ok"),
    rank         = NA_integer_,
    stringsAsFactors = FALSE
  )
  h_ranked <- rank_grid_history(h)
  # score 0.3 should be rank 1, 0.5 rank 2, 0.7 rank 3
  expect_equal(h_ranked$rank[h$score == 0.3], 1L)
  expect_equal(h_ranked$rank[h$score == 0.5], 2L)
  expect_equal(h_ranked$rank[h$score == 0.7], 3L)
})

test_that("rank_grid_history places failed candidates after ok candidates", {
  h <- data.frame(
    candidate_id = 1:4,
    score        = c(0.5, NA, 0.3, NA),
    status       = c("ok", "failed", "ok", "failed"),
    rank         = NA_integer_,
    stringsAsFactors = FALSE
  )
  h_ranked <- rank_grid_history(h)
  ok_ranks   <- h_ranked$rank[h_ranked$status == "ok"]
  fail_ranks <- h_ranked$rank[h_ranked$status == "failed"]
  expect_true(all(ok_ranks < min(fail_ranks)))
})

test_that("rank_grid_history handles all-failed history", {
  h <- data.frame(
    candidate_id = 1:2,
    score        = c(NA_real_, NA_real_),
    status       = c("failed", "failed"),
    rank         = NA_integer_,
    stringsAsFactors = FALSE
  )
  h_ranked <- rank_grid_history(h)
  expect_equal(sort(h_ranked$rank), c(1L, 2L))
})

test_that("rank_grid_history returns data.frame unchanged if empty", {
  h <- data.frame()
  expect_equal(rank_grid_history(h), h)
})

# ===========================================================================
# search_bandwidth_grid
# ===========================================================================

test_that("search_bandwidth_grid errors when lower is missing from control", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  expect_error(
    search_bandwidth_grid(ctx, control = list(upper = 5, step = 1), scorer = scorer),
    regexp = "lower"
  )
})

test_that("search_bandwidth_grid errors when upper is missing from control", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  expect_error(
    search_bandwidth_grid(ctx, control = list(lower = 1, step = 1), scorer = scorer),
    regexp = "upper"
  )
})

test_that("search_bandwidth_grid errors when step is missing from control", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  expect_error(
    search_bandwidth_grid(ctx, control = list(lower = 1, upper = 5), scorer = scorer),
    regexp = "step"
  )
})

test_that("search_bandwidth_grid returns a gwpr_bandwidth object", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 3)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  expect_s3_class(result, "gwpr_bandwidth")
})

test_that("search_bandwidth_grid history has correct number of rows", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 3)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  # lower=1, upper=5, step=1 => 5 candidates
  expect_equal(nrow(result$history), 5L)
})

test_that("search_bandwidth_grid best_bandwidth is the candidate with lowest score", {
  ctx    <- make_grid_context()
  # best_bw=3: score for bw=3 is 0.01, for others > 0.01
  scorer <- make_mock_scorer(best_bw = 3)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  expect_equal(result$best_bandwidth, 3)
})

test_that("search_bandwidth_grid history includes all required columns", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 3, step = 1),
    scorer  = scorer,
    workers = 1L
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

test_that("search_bandwidth_grid retains failed candidates in history", {
  ctx         <- make_grid_context()
  # bw=2 will fail
  scorer <- make_mock_scorer(best_bw = 3, fail_bw = 2)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  # All 5 candidates should be present
  expect_equal(nrow(result$history), 5L)
  # The failed candidate must be recorded with status "failed"
  expect_true(any(result$history$status == "failed"))
  expect_true(any(result$history$bandwidth == 2 & result$history$status == "failed"))
})

test_that("search_bandwidth_grid best_bandwidth comes from ok candidates only", {
  ctx    <- make_grid_context()
  # Make bw=3 fail; next best should be bw=4 (score=1.01) or bw=2 (score=1.01)
  scorer <- make_mock_scorer(best_bw = 3, fail_bw = 3)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  # best_bandwidth should NOT be 3 (it failed)
  expect_false(is.na(result$best_bandwidth))
  expect_false(result$best_bandwidth == 3)
})

test_that("search_bandwidth_grid method field is 'grid'", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 3, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  expect_equal(result$method, "grid")
  expect_true(all(result$history$method == "grid"))
})

test_that("search_bandwidth_grid records seed in returned object", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 3, step = 1),
    scorer  = scorer,
    workers = 1L,
    seed    = 42L
  )
  expect_equal(result$seed, 42L)
})

test_that("search_bandwidth_grid elapsed_time is non-negative", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 3, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  expect_true(result$elapsed_time >= 0)
})

test_that("search_bandwidth_grid adaptive grid produces integer candidates", {
  ctx    <- make_grid_context(adaptive = TRUE)
  scorer <- make_mock_scorer(best_bw = 3)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 2, upper = 8, step = 2),
    scorer  = scorer,
    workers = 1L
  )
  expect_true(all(result$history$bandwidth == as.integer(result$history$bandwidth)))
})

test_that("search_bandwidth_grid: best_bandwidth has rank 1 in history", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer(best_bw = 2)
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 5, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  best_row <- result$history[result$history$bandwidth == result$best_bandwidth, , drop = FALSE]
  expect_equal(best_row$rank[1L], 1L)
})

test_that("search_bandwidth_grid convergence_info contains n_candidates", {
  ctx    <- make_grid_context()
  scorer <- make_mock_scorer()
  result <- search_bandwidth_grid(
    ctx,
    control = list(lower = 1, upper = 4, step = 1),
    scorer  = scorer,
    workers = 1L
  )
  expect_equal(result$convergence_info$n_candidates, 4L)
})
