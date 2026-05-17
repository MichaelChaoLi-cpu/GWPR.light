# tests/testthat/test-diagnostics.R
#
# Tests for the diagnostics module (R/diagnostics.R).
# All tests use small simulated data; workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helper: build a small balanced-panel matrix and weights
# ---------------------------------------------------------------------------

make_balanced_moran_inputs <- function(n = 5L, Ti = 3L, seed = 42L) {
  set.seed(seed)
  resid_mat <- matrix(rnorm(n * Ti), nrow = n, ncol = Ti)

  # Simple row-standardised weights (uniform neighbours, no self-loops)
  W <- matrix(1 / (n - 1), nrow = n, ncol = n)
  diag(W) <- 0

  list(resid_mat = resid_mat, W = W, n = n, Ti = Ti)
}

make_panel_index <- function(n, Ti) {
  data.frame(
    id   = rep(seq_len(n), each = Ti),
    time = rep(seq_len(Ti), times = n)
  )
}

# ---------------------------------------------------------------------------
# compute_logistic_pearson_residual
# ---------------------------------------------------------------------------

test_that("pearson residual: standard case produces finite values", {
  y    <- c(1L, 0L, 1L, 0L)
  prob <- c(0.9, 0.1, 0.8, 0.2)
  r    <- compute_logistic_pearson_residual(y, prob)
  expect_length(r, 4L)
  expect_true(all(is.finite(r)))
})

test_that("pearson residual: prob = 0 or 1 does not produce Inf", {
  y    <- c(1L, 0L)
  prob <- c(1.0, 0.0)   # boundary values — would divide by zero without clipping
  r    <- compute_logistic_pearson_residual(y, prob)
  expect_true(all(is.finite(r)),
              info = "Clipping must prevent Inf residuals at prob = 0/1.")
})

test_that("pearson residual: formula correctness", {
  y    <- 1L
  prob <- 0.75
  eps  <- 1e-15
  p    <- max(min(prob, 1 - eps), eps)
  expected <- (y - p) / sqrt(p * (1 - p))
  r <- compute_logistic_pearson_residual(y, prob, eps = eps)
  expect_equal(r, expected, tolerance = 1e-12)
})

test_that("pearson residual: length mismatch is an error", {
  expect_error(
    compute_logistic_pearson_residual(c(1, 0), c(0.9)),
    regexp = "same length"
  )
})

test_that("pearson residual: non-numeric y is an error", {
  expect_error(
    compute_logistic_pearson_residual("a", 0.5),
    regexp = "numeric or integer"
  )
})

test_that("pearson residual: negative eps is an error", {
  expect_error(
    compute_logistic_pearson_residual(1L, 0.5, eps = -1),
    regexp = "positive"
  )
})

# ---------------------------------------------------------------------------
# compute_panel_moran — balanced panel (matrix input)
# ---------------------------------------------------------------------------

test_that("moran balanced panel: returns list with correct names", {
  inp <- make_balanced_moran_inputs()
  res <- compute_panel_moran(inp$resid_mat, inp$W)

  expected_names <- c("statistic", "p_value", "estimated_I", "expected_I",
                      "variance", "alternative", "n_individuals", "n_periods")
  expect_true(all(expected_names %in% names(res)))
})

test_that("moran balanced panel: statistic and p_value are finite", {
  inp <- make_balanced_moran_inputs()
  res <- compute_panel_moran(inp$resid_mat, inp$W)

  expect_true(is.finite(res$statistic),
              info = "Z-statistic should be finite for well-specified input.")
  expect_true(is.finite(res$p_value))
})

test_that("moran balanced panel: p_value in [0, 1]", {
  inp <- make_balanced_moran_inputs()
  res <- compute_panel_moran(inp$resid_mat, inp$W)
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

test_that("moran balanced panel: n_individuals and n_periods are correct", {
  inp <- make_balanced_moran_inputs(n = 6L, Ti = 4L)
  res <- compute_panel_moran(inp$resid_mat, inp$W)
  expect_equal(res$n_individuals, 6L)
  expect_equal(res$n_periods, 4L)
})

test_that("moran balanced panel: expected_I equals -1/(n-1)", {
  n   <- 5L
  inp <- make_balanced_moran_inputs(n = n)
  res <- compute_panel_moran(inp$resid_mat, inp$W)
  expect_equal(res$expected_I, -1 / (n - 1))
})

test_that("moran balanced panel: alternative is 'greater'", {
  inp <- make_balanced_moran_inputs()
  res <- compute_panel_moran(inp$resid_mat, inp$W)
  expect_equal(res$alternative, "greater")
})

# ---------------------------------------------------------------------------
# compute_panel_moran — vector input with panel_index
# ---------------------------------------------------------------------------

test_that("moran: vector input + panel_index balanced panel runs without error", {
  n  <- 4L
  Ti <- 3L
  set.seed(7L)
  resid_vec <- rnorm(n * Ti)
  idx       <- make_panel_index(n, Ti)
  W         <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  res <- compute_panel_moran(resid_vec, W, panel_index = idx)
  expect_true(is.list(res))
  expect_true(is.finite(res$statistic))
})

test_that("moran: vector input requires panel_index", {
  W   <- matrix(1 / 4, 5, 5); diag(W) <- 0
  expect_error(
    compute_panel_moran(rnorm(15), W),
    regexp = "panel_index"
  )
})

# ---------------------------------------------------------------------------
# compute_panel_moran — unbalanced panel warning
# ---------------------------------------------------------------------------

test_that("moran: unbalanced panel gives a warning", {
  n  <- 4L
  Ti <- 3L
  set.seed(11L)
  resid_vec <- rnorm(n * Ti - 1L)   # one observation missing => unbalanced
  idx <- make_panel_index(n, Ti)
  idx <- idx[-nrow(idx), ]          # drop last row to make it unbalanced
  W   <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  expect_warning(
    compute_panel_moran(resid_vec, W, panel_index = idx),
    regexp = "unbalanced"
  )
})

# ---------------------------------------------------------------------------
# compute_panel_moran — dimension mismatch error
# ---------------------------------------------------------------------------

test_that("moran: spatial_weights dimension mismatch is an error", {
  n   <- 5L
  Ti  <- 3L
  inp <- make_balanced_moran_inputs(n = n, Ti = Ti)

  W_wrong <- matrix(1 / 3, nrow = 3L, ncol = 3L)   # n=3 != 5
  diag(W_wrong) <- 0

  expect_error(
    compute_panel_moran(inp$resid_mat, W_wrong),
    regexp = "n x n"
  )
})

# ---------------------------------------------------------------------------
# diagnose_moran
# ---------------------------------------------------------------------------

test_that("diagnose_moran: linear model runs and returns residual_type = 'raw'", {
  n  <- 4L
  Ti <- 3L
  set.seed(1L)
  resid_vec <- rnorm(n * Ti)
  idx       <- make_panel_index(n, Ti)
  W         <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    effect    = "individual",
    bandwidth = 1,
    residuals = resid_vec
  )

  res <- diagnose_moran(fit, W, idx)
  expect_equal(res$residual_type, "raw")
  expect_true(is.finite(res$statistic))
})

test_that("diagnose_moran: logistic model uses Pearson residuals", {
  n  <- 4L
  Ti <- 3L
  N  <- n * Ti
  set.seed(2L)
  y    <- sample(0:1, N, replace = TRUE)
  prob <- runif(N, 0.1, 0.9)
  idx  <- make_panel_index(n, Ti)
  W    <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  fit <- new_gwpr_fit(
    family    = "binomial",
    model     = "pooling",
    bandwidth = 1,
    residuals = y - prob,           # not used directly — overridden by Pearson
    metadata  = list(response = y, prob = prob)
  )

  res <- diagnose_moran(fit, W, idx)
  expect_equal(res$residual_type, "pearson")
  expect_true(is.finite(res$statistic))
})

test_that("diagnose_moran: non-gwpr_fit object is an error", {
  expect_error(
    diagnose_moran(list(), matrix(1, 2, 2), data.frame(id = 1, time = 1)),
    regexp = "gwpr_fit"
  )
})

test_that("diagnose_moran: NULL residuals is an error", {
  fit <- new_gwpr_fit(family = "gaussian", residuals = NULL)
  expect_error(
    diagnose_moran(fit, matrix(1, 2, 2), data.frame(id = 1, time = 1)),
    regexp = "residuals.*NULL"
  )
})

# ---------------------------------------------------------------------------
# diagnose_local_f
# ---------------------------------------------------------------------------

test_that("diagnose_local_f: returns data.frame with correct columns", {
  local_res <- list(
    "1" = list(within_rss = 2.0, within_df = 5L, pooling_rss = 5.0, pooling_df = 8L, status = "ok"),
    "2" = list(within_rss = 1.5, within_df = 5L, pooling_rss = 4.0, pooling_df = 8L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_local_f(fit)

  expect_true(is.data.frame(res$local_f))
  expect_true(all(c("unit_id", "statistic", "p_value", "df1", "df2", "status") %in%
                    names(res$local_f)))
  expect_equal(nrow(res$local_f), 2L)
})

test_that("diagnose_local_f: computed F-statistic is positive", {
  local_res <- list(
    "A" = list(within_rss = 2.0, within_df = 5L, pooling_rss = 5.0, pooling_df = 8L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_local_f(fit)

  expect_gt(res$local_f$statistic[1L], 0)
})

test_that("diagnose_local_f: p-value in [0, 1]", {
  local_res <- list(
    "A" = list(within_rss = 2.0, within_df = 5L, pooling_rss = 5.0, pooling_df = 8L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_local_f(fit)

  expect_gte(res$local_f$p_value[1L], 0)
  expect_lte(res$local_f$p_value[1L], 1)
})

test_that("diagnose_local_f: binomial model returns not_applicable", {
  fit <- new_gwpr_fit(family = "binomial", model = "pooling")
  res <- diagnose_local_f(fit)

  expect_equal(res$status, "not_applicable")
  expect_equal(nrow(res$local_f), 0L)
})

test_that("diagnose_local_f: pre-computed f_statistic used directly", {
  local_res <- list(
    "X" = list(f_statistic = 4.5, f_p_value = 0.02, f_df1 = 2L, f_df2 = 10L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_local_f(fit)

  expect_equal(res$local_f$statistic[1L], 4.5)
  expect_equal(res$local_f$p_value[1L], 0.02)
})

test_that("diagnose_local_f: missing data unit gets correct status", {
  local_res <- list(
    "Z" = list(status = "ok")   # no within_rss / pooling_rss / f_statistic
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_local_f(fit)

  expect_equal(res$local_f$status[1L], "missing_local_f_data")
  expect_equal(res$n_failed, 1L)
})

test_that("diagnose_local_f: empty local_results returns no_local_results", {
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = list())
  res <- diagnose_local_f(fit)
  expect_equal(res$status, "no_local_results")
})

# ---------------------------------------------------------------------------
# diagnose_hausman
# ---------------------------------------------------------------------------

test_that("diagnose_hausman: returns data.frame with correct columns", {
  local_res <- list(
    "1" = list(hausman_statistic = 3.2, hausman_p_value = 0.07, hausman_df = 2L, status = "ok"),
    "2" = list(hausman_statistic = 6.1, hausman_p_value = 0.01, hausman_df = 2L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "random", local_results = local_res)
  res <- diagnose_hausman(fit)

  expect_true(is.data.frame(res$local_hausman))
  expect_true(all(c("unit_id", "statistic", "p_value", "df", "status") %in%
                    names(res$local_hausman)))
  expect_equal(nrow(res$local_hausman), 2L)
})

test_that("diagnose_hausman: pooling model returns not_applicable", {
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling")
  res <- diagnose_hausman(fit)

  expect_equal(res$status, "not_applicable")
  expect_match(res$message, "pooling")
  expect_equal(nrow(res$local_hausman), 0L)
})

test_that("diagnose_hausman: binomial model returns not_applicable", {
  fit <- new_gwpr_fit(family = "binomial", model = "random")
  res <- diagnose_hausman(fit)

  expect_equal(res$status, "not_applicable")
  expect_equal(nrow(res$local_hausman), 0L)
})

test_that("diagnose_hausman: missing data unit is counted as failed", {
  local_res <- list(
    "A" = list(status = "ok")   # no hausman_statistic
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "random", local_results = local_res)
  res <- diagnose_hausman(fit)

  expect_equal(res$local_hausman$status[1L], "missing_hausman_data")
  expect_equal(res$n_failed, 1L)
})

test_that("diagnose_hausman: within model runs (applicable)", {
  local_res <- list(
    "1" = list(hausman_statistic = 2.5, hausman_p_value = 0.11, hausman_df = 2L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  res <- diagnose_hausman(fit)

  expect_equal(res$status, "ok")
  expect_equal(res$n_tested, 1L)
})

# ---------------------------------------------------------------------------
# diagnose_lm
# ---------------------------------------------------------------------------

test_that("diagnose_lm: returns data.frame with correct columns", {
  local_res <- list(
    "1" = list(lm_statistic = 4.1, lm_p_value = 0.04, lm_df = 1L, status = "ok"),
    "2" = list(lm_statistic = 1.2, lm_p_value = 0.27, lm_df = 1L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling", local_results = local_res)
  res <- diagnose_lm(fit)

  expect_true(is.data.frame(res$local_lm))
  expect_true(all(c("unit_id", "statistic", "p_value", "df", "status") %in%
                    names(res$local_lm)))
  expect_equal(nrow(res$local_lm), 2L)
})

test_that("diagnose_lm: within model returns not_applicable", {
  fit <- new_gwpr_fit(family = "gaussian", model = "within")
  res <- diagnose_lm(fit)

  expect_equal(res$status, "not_applicable")
  expect_match(res$message, "within")
  expect_equal(nrow(res$local_lm), 0L)
})

test_that("diagnose_lm: binomial model returns not_applicable", {
  fit <- new_gwpr_fit(family = "binomial", model = "pooling")
  res <- diagnose_lm(fit)

  expect_equal(res$status, "not_applicable")
  expect_equal(nrow(res$local_lm), 0L)
})

test_that("diagnose_lm: missing data unit is counted as failed", {
  local_res <- list(
    "A" = list(status = "ok")   # no lm_statistic
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling", local_results = local_res)
  res <- diagnose_lm(fit)

  expect_equal(res$local_lm$status[1L], "missing_lm_data")
  expect_equal(res$n_failed, 1L)
})

test_that("diagnose_lm: random model is applicable", {
  local_res <- list(
    "X" = list(lm_statistic = 3.3, lm_p_value = 0.07, lm_df = 1L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "random", local_results = local_res)
  res <- diagnose_lm(fit)
  expect_equal(res$status, "ok")
})

# ---------------------------------------------------------------------------
# diagnose_gwpr — top-level interface
# ---------------------------------------------------------------------------

test_that("diagnose_gwpr: returns gwpr_diagnostics object", {
  local_res <- list(
    "1" = list(within_rss = 2.0, within_df = 5L, pooling_rss = 5.0, pooling_df = 8L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "within", local_results = local_res)
  dg  <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman", "lm_test"))

  expect_s3_class(dg, "gwpr_diagnostics")
})

test_that("diagnose_gwpr: result object has correct structure", {
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling")
  dg  <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman", "lm_test"))

  expect_true(is.list(dg$diagnostics))
  expect_true(!is.null(dg$model_type))
  expect_true(is.character(dg$warnings))
})

test_that("diagnose_gwpr: moran skipped without spatial_weights/panel_index", {
  fit <- new_gwpr_fit(family = "gaussian", model = "within", residuals = rnorm(12))
  dg  <- diagnose_gwpr(fit, diagnostics = "moran")

  expect_equal(dg$diagnostics$moran$status, "skipped")
})

test_that("diagnose_gwpr: moran runs when spatial_weights and panel_index provided", {
  n  <- 4L
  Ti <- 3L
  set.seed(5L)
  resid_vec <- rnorm(n * Ti)
  idx       <- make_panel_index(n, Ti)
  W         <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    bandwidth = 1,
    residuals = resid_vec
  )

  dg <- diagnose_gwpr(fit,
                      diagnostics     = "moran",
                      spatial_weights = W,
                      panel_index     = idx)

  moran_res <- dg$diagnostics$moran
  expect_false(identical(moran_res$status, "skipped"),
               info = "Moran should run, not be skipped, when inputs are provided.")
  expect_true(is.finite(moran_res$statistic))
})

test_that("diagnose_gwpr: non-applicable tests return interpretable status", {
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling")
  dg  <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman", "lm_test"))

  # f_test: no local_results -> no_local_results
  expect_true(!is.null(dg$diagnostics$f_test))
  # hausman: pooling -> not_applicable
  expect_equal(dg$diagnostics$hausman$status, "not_applicable")
  # lm_test: pooling and no local results -> no_local_results
  expect_true(!is.null(dg$diagnostics$lm_test))
})

test_that("diagnose_gwpr: non-gwpr_fit input gives error", {
  expect_error(
    diagnose_gwpr(list()),
    regexp = "gwpr_fit"
  )
})

test_that("diagnose_gwpr: panel_balance is set correctly for balanced panel", {
  n  <- 3L
  Ti <- 2L
  resid_vec <- rnorm(n * Ti)
  idx       <- make_panel_index(n, Ti)
  W         <- matrix(1 / (n - 1), n, n); diag(W) <- 0

  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    bandwidth = 1,
    residuals = resid_vec
  )

  dg <- diagnose_gwpr(fit,
                      diagnostics     = "moran",
                      spatial_weights = W,
                      panel_index     = idx)

  expect_true(dg$panel_balance)
})

test_that("diagnose_gwpr: subset of diagnostics only runs requested tests", {
  local_res <- list(
    "1" = list(lm_statistic = 2.0, lm_p_value = 0.15, lm_df = 1L, status = "ok")
  )
  fit <- new_gwpr_fit(family = "gaussian", model = "pooling", local_results = local_res)
  dg  <- diagnose_gwpr(fit, diagnostics = "lm_test")

  expect_true("lm_test" %in% names(dg$diagnostics))
  expect_false("f_test" %in% names(dg$diagnostics))
  expect_false("hausman" %in% names(dg$diagnostics))
  expect_false("moran" %in% names(dg$diagnostics))
})
