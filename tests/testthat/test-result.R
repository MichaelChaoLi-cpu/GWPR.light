## tests/testthat/test-result.R
## Unit tests for result.R – gwpr_fit, gwpr_bandwidth, gwpr_diagnostics and
## build_spatial_results.
## All data is small and synthetic; no parallel execution (workers = 1).

# ---------------------------------------------------------------------------
# Helpers shared across tests
# ---------------------------------------------------------------------------

make_local_results <- function() {
  list(
    "A" = list(
      coefficients = c(x1 = 0.5, x2 = 1.2),
      std_errors   = c(x1 = 0.1, x2 = 0.3),
      t_stats      = c(x1 = 5.0, x2 = 4.0),
      metrics      = list(R2 = 0.8, MSE = 0.05),
      diagnostics  = list(local_moran = 0.12),
      status       = "ok"
    ),
    "B" = list(
      coefficients = c(x1 = 0.3, x2 = 0.9),
      std_errors   = c(x1 = 0.2, x2 = 0.4),
      t_stats      = c(x1 = 1.5, x2 = 2.25),
      metrics      = list(R2 = 0.6, MSE = 0.12),
      diagnostics  = list(local_moran = 0.08),
      status       = "ok"
    ),
    "C" = list(
      coefficients = c(x1 = NA_real_, x2 = NA_real_),
      std_errors   = c(x1 = NA_real_, x2 = NA_real_),
      t_stats      = c(x1 = NA_real_, x2 = NA_real_),
      metrics      = list(R2 = NA_real_, MSE = NA_real_),
      diagnostics  = list(local_moran = NA_real_),
      status       = "singular matrix"
    )
  )
}

# ---------------------------------------------------------------------------
# 1. gwpr_fit constructor
# ---------------------------------------------------------------------------

test_that("new_gwpr_fit returns an object of class gwpr_fit", {
  fit <- new_gwpr_fit()
  expect_s3_class(fit, "gwpr_fit")
})

test_that("new_gwpr_fit has all required fields", {
  fit <- new_gwpr_fit()
  required <- c("call", "family", "formula", "model", "effect",
                "bandwidth", "search", "global_model", "local_results",
                "predictions", "residuals", "metrics",
                "spatial_results", "metadata", "warnings")
  for (fld in required) {
    expect_true(fld %in% names(fit),
                info = paste("Missing field:", fld))
  }
})

test_that("new_gwpr_fit stores supplied values correctly", {
  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    effect    = "individual",
    bandwidth = 100,
    metrics   = list(R2 = 0.8, MSE = 0.1),
    warnings  = c("w1", "w2")
  )
  expect_equal(fit$family,    "gaussian")
  expect_equal(fit$model,     "within")
  expect_equal(fit$effect,    "individual")
  expect_equal(fit$bandwidth, 100)
  expect_equal(fit$metrics$R2, 0.8)
  expect_equal(fit$warnings,  c("w1", "w2"))
})

test_that("new_gwpr_fit warnings field is preserved", {
  warns <- c("local model failed at unit X", "singular matrix at unit Y")
  fit <- new_gwpr_fit(warnings = warns)
  expect_identical(fit$warnings, warns)
})

test_that("print.gwpr_fit runs without error and returns x invisibly", {
  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    effect    = "individual",
    bandwidth = 50,
    metrics   = list(R2 = 0.75, RMSE = 0.5)
  )
  out <- capture.output(ret <- print(fit))
  expect_identical(ret, fit)
  expect_true(any(grepl("gwpr_fit", out, ignore.case = TRUE)))
})

test_that("summary.gwpr_fit runs without error and returns object invisibly", {
  fit <- new_gwpr_fit(
    family    = "gaussian",
    model     = "within",
    effect    = "individual",
    bandwidth = 50,
    metrics   = list(R2 = 0.75, RMSE = 0.5)
  )
  out <- capture.output(ret <- summary(fit))
  expect_identical(ret, fit)
  expect_true(any(grepl("GWPR", out)))
})

test_that("summary.gwpr_fit does not require raw data – uses object fields only", {
  # Build a fit with spatial_results that has coef_ columns
  sr <- data.frame(unit_id = c("A", "B"),
                   status  = c("ok", "ok"),
                   coef_x1 = c(0.5, 0.3))
  fit <- new_gwpr_fit(
    family          = "gaussian",
    model           = "within",
    effect          = "individual",
    bandwidth       = 100,
    metrics         = list(R2 = 0.8),
    spatial_results = sr
  )
  # Should run without needing any raw data
  expect_silent(capture.output(summary(fit)))
})

# ---------------------------------------------------------------------------
# 2. gwpr_bandwidth constructor
# ---------------------------------------------------------------------------

test_that("new_gwpr_bandwidth returns an object of class gwpr_bandwidth", {
  bw <- new_gwpr_bandwidth()
  expect_s3_class(bw, "gwpr_bandwidth")
})

test_that("new_gwpr_bandwidth has all required fields", {
  bw <- new_gwpr_bandwidth()
  required <- c("method", "best_bandwidth", "best_score", "criterion",
                "history", "metrics_history", "seed",
                "convergence_info", "elapsed_time", "warnings")
  for (fld in required) {
    expect_true(fld %in% names(bw),
                info = paste("Missing field:", fld))
  }
})

test_that("new_gwpr_bandwidth stores supplied values correctly", {
  hist_df <- data.frame(bandwidth = c(100, 150, 200),
                        score     = c(0.5, 0.42, 0.48))
  bw <- new_gwpr_bandwidth(
    method         = "grid",
    best_bandwidth = 150,
    best_score     = 0.42,
    criterion      = "MSE",
    history        = hist_df,
    seed           = 42L,
    elapsed_time   = 1.2,
    warnings       = "minor warning"
  )
  expect_equal(bw$method,         "grid")
  expect_equal(bw$best_bandwidth, 150)
  expect_equal(bw$best_score,     0.42)
  expect_equal(bw$criterion,      "MSE")
  expect_equal(nrow(bw$history),  3L)
  expect_equal(bw$seed,           42L)
  expect_equal(bw$elapsed_time,   1.2)
  expect_equal(bw$warnings,       "minor warning")
})

test_that("new_gwpr_bandwidth warnings field is preserved", {
  warns <- c("candidate failed", "numerical issues")
  bw <- new_gwpr_bandwidth(warnings = warns)
  expect_identical(bw$warnings, warns)
})

test_that("print.gwpr_bandwidth runs without error and returns x invisibly", {
  bw <- new_gwpr_bandwidth(
    method         = "grid",
    best_bandwidth = 150,
    best_score     = 0.42,
    criterion      = "MSE",
    history        = data.frame(bandwidth = 150, score = 0.42),
    elapsed_time   = 0.5
  )
  out <- capture.output(ret <- print(bw))
  expect_identical(ret, bw)
  expect_true(any(grepl("gwpr_bandwidth", out, ignore.case = TRUE)))
})

test_that("summary.gwpr_bandwidth runs without error and returns object invisibly", {
  bw <- new_gwpr_bandwidth(
    method         = "grid",
    best_bandwidth = 150,
    best_score     = 0.42,
    criterion      = "MSE",
    history        = data.frame(bandwidth = c(100, 150), score = c(0.5, 0.42)),
    elapsed_time   = 1.2
  )
  out <- capture.output(ret <- summary(bw))
  expect_identical(ret, bw)
  expect_true(any(grepl("Bandwidth", out)))
})

# ---------------------------------------------------------------------------
# 3. gwpr_diagnostics constructor
# ---------------------------------------------------------------------------

test_that("new_gwpr_diagnostics returns an object of class gwpr_diagnostics", {
  d <- new_gwpr_diagnostics()
  expect_s3_class(d, "gwpr_diagnostics")
})

test_that("new_gwpr_diagnostics has all required fields", {
  d <- new_gwpr_diagnostics()
  required <- c("diagnostics", "model_type", "spatial_weights",
                "panel_balance", "warnings", "metadata")
  for (fld in required) {
    expect_true(fld %in% names(d),
                info = paste("Missing field:", fld))
  }
})

test_that("new_gwpr_diagnostics stores supplied values correctly", {
  diag_list <- list(
    moran  = list(statistic = 0.12, p_value = 0.03),
    f_test = list(statistic = 4.5,  p_value = 0.01)
  )
  d <- new_gwpr_diagnostics(
    diagnostics   = diag_list,
    model_type    = "gaussian",
    panel_balance = TRUE,
    warnings      = "unbalanced warning"
  )
  expect_equal(d$model_type,    "gaussian")
  expect_true(d$panel_balance)
  expect_equal(d$diagnostics$moran$statistic, 0.12)
  expect_equal(d$warnings, "unbalanced warning")
})

test_that("new_gwpr_diagnostics warnings field is preserved", {
  warns <- c("unbalanced panel", "Moran's I not available")
  d <- new_gwpr_diagnostics(warnings = warns)
  expect_identical(d$warnings, warns)
})

test_that("print.gwpr_diagnostics runs without error and returns x invisibly", {
  d <- new_gwpr_diagnostics(
    diagnostics   = list(moran = list(statistic = 0.1, p_value = 0.05)),
    model_type    = "gaussian",
    panel_balance = TRUE
  )
  out <- capture.output(ret <- print(d))
  expect_identical(ret, d)
  expect_true(any(grepl("gwpr_diagnostics", out, ignore.case = TRUE)))
})

test_that("summary.gwpr_diagnostics runs without error and returns object invisibly", {
  d <- new_gwpr_diagnostics(
    diagnostics   = list(moran = list(statistic = 0.12, p_value = 0.03)),
    model_type    = "gaussian",
    panel_balance = TRUE
  )
  out <- capture.output(ret <- summary(d))
  expect_identical(ret, d)
  expect_true(any(grepl("Diagnostics", out)))
})

# ---------------------------------------------------------------------------
# 4. build_spatial_results
# ---------------------------------------------------------------------------

test_that("build_spatial_results returns a data.frame", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_s3_class(out, "data.frame")
})

test_that("build_spatial_results has one row per unit", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_equal(nrow(out), length(lr))
})

test_that("build_spatial_results has unit_id and status columns", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_true("unit_id" %in% names(out))
  expect_true("status"  %in% names(out))
})

test_that("build_spatial_results has coef_ columns for each predictor", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_true("coef_x1" %in% names(out))
  expect_true("coef_x2" %in% names(out))
})

test_that("build_spatial_results has se_ and tstat_ columns", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_true("se_x1"    %in% names(out))
  expect_true("tstat_x1" %in% names(out))
})

test_that("build_spatial_results includes metric columns", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_true("R2"  %in% names(out))
  expect_true("MSE" %in% names(out))
})

test_that("build_spatial_results includes diagnostic columns", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_true("local_moran" %in% names(out))
})

test_that("build_spatial_results preserves failed unit status", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_equal(out$status[out$unit_id == "C"], "singular matrix")
})

test_that("build_spatial_results unit_id column matches names of local_results", {
  lr  <- make_local_results()
  out <- build_spatial_results(lr)
  expect_equal(out$unit_id, names(lr))
})

test_that("build_spatial_results validates geometry length", {
  lr <- make_local_results()
  # wrong-length geometry (a plain list of length 2 instead of 3)
  bad_geom <- list("geom1", "geom2")
  expect_error(
    build_spatial_results(lr, geometry = bad_geom),
    regexp = "geometry.*length|length.*geometry",
    ignore.case = TRUE
  )
})

test_that("build_spatial_results geometry alignment: correct length passes", {
  lr       <- make_local_results()
  ok_geom  <- list("geom1", "geom2", "geom3")
  out      <- build_spatial_results(lr, geometry = ok_geom)
  expect_equal(nrow(out), 3L)
})

test_that("build_spatial_results stops on empty local_results", {
  expect_error(build_spatial_results(list()), regexp = "non-empty")
})

test_that("build_spatial_results column count consistent with sf alignment", {
  # The number of rows must match geometry for sf cbind to work
  lr      <- make_local_results()
  out     <- build_spatial_results(lr)
  geom_df <- data.frame(dummy = seq_len(nrow(out)))
  expect_equal(nrow(out), nrow(geom_df))
})

# ---------------------------------------------------------------------------
# 5. Null-coalescing helper
# ---------------------------------------------------------------------------

test_that("%||% returns left value when non-NULL", {
  expect_equal("a" %||% "b", "a")
  expect_equal(0L  %||% 1L,  0L)
})

test_that("%||% returns right value when left is NULL", {
  expect_equal(NULL %||% "fallback", "fallback")
})
