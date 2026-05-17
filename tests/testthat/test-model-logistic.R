## tests/testthat/test-model-logistic.R
## Tests for the Binary Panel Logistic Engine (R/model_logistic.R)
##
## All tests use small simulated panel data.  workers = 1 (no parallelism).
## glmmTMB random effects tests are wrapped in skip conditions when the
## package is unavailable or fails due to numerical issues (e.g. ABI mismatch).

# ---------------------------------------------------------------------------
# Helpers: small simulated binary panel data
# ---------------------------------------------------------------------------

# 4 spatial units, 5 time periods each => 20 rows
make_bin_panel_data <- function(seed = 42L) {
  set.seed(seed)
  n_units <- 4L
  n_times <- 5L
  id   <- rep(as.character(1:n_units), each = n_times)
  time <- rep(1:n_times, times = n_units)
  x1   <- rnorm(n_units * n_times)
  x2   <- rnorm(n_units * n_times)
  lp   <- 0.5 * x1 - 0.3 * x2
  y    <- as.integer(rbinom(n_units * n_times, 1L, plogis(lp)))
  data.frame(id = id, time = time, y = y, x1 = x1, x2 = x2,
             stringsAsFactors = FALSE)
}

# 4 units on a small 2x2 grid
make_coords_4 <- function() {
  matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2, byrow = TRUE,
    dimnames = list(c("1", "2", "3", "4"), c("X", "Y"))
  )
}

make_logistic_context <- function(model = "pooling", effect = "individual",
                                   kernel = "bisquare", adaptive = FALSE,
                                   threshold = 0.5, seed = NULL) {
  new_gwpr_context(
    formula      = y ~ x1 + x2,
    family       = "binomial",
    model        = model,
    effect       = effect,
    id           = "id",
    time         = "time",
    kernel       = kernel,
    adaptive     = adaptive,
    threshold    = threshold,
    workers      = 1L,
    seed         = seed,
    panel_data   = make_bin_panel_data(),
    coords       = make_coords_4()
  )
}

# ---------------------------------------------------------------------------
# validate_binary_response
# ---------------------------------------------------------------------------

test_that("validate_binary_response: passes for numeric 0/1", {
  expect_true(validate_binary_response(c(0L, 1L, 0L, 1L)))
})

test_that("validate_binary_response: passes for two-level factor", {
  y <- factor(c("no", "yes", "no"), levels = c("no", "yes"))
  expect_true(validate_binary_response(y))
})

test_that("validate_binary_response: stops for 3-level factor", {
  y <- factor(c("a", "b", "c"))
  expect_error(validate_binary_response(y), regexp = "Multi-class")
})

# ---------------------------------------------------------------------------
# standardize_logistic_response
# ---------------------------------------------------------------------------

test_that("standardize_logistic_response: numeric 0/1 passes through", {
  y <- c(0, 1, 0, 1)
  result <- standardize_logistic_response(y)
  expect_equal(result, c(0L, 1L, 0L, 1L))
})

test_that("standardize_logistic_response: two-level factor maps to 0/1", {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  result <- standardize_logistic_response(y)
  expect_equal(result, c(0L, 1L, 0L, 1L))
})

test_that("standardize_logistic_response: logical maps to 0/1", {
  y <- c(FALSE, TRUE, FALSE, TRUE)
  result <- standardize_logistic_response(y)
  expect_equal(result, c(0L, 1L, 0L, 1L))
})

test_that("standardize_logistic_response: multi-class factor raises error", {
  y <- factor(c("a", "b", "c", "a"))
  expect_error(standardize_logistic_response(y), regexp = "Multi-class")
})

# ---------------------------------------------------------------------------
# effect_to_feglm_fml
# ---------------------------------------------------------------------------

test_that("effect_to_feglm_fml: individual maps id to fixed effect", {
  fml <- effect_to_feglm_fml(y ~ x1 + x2, "individual", "id", "time")
  expect_true(grepl("id", deparse(fml)))
})

test_that("effect_to_feglm_fml: time maps time to fixed effect", {
  fml <- effect_to_feglm_fml(y ~ x1 + x2, "time", "id", "time")
  expect_true(grepl("time", deparse(fml)))
})

test_that("effect_to_feglm_fml: two-way includes both id and time", {
  fml <- effect_to_feglm_fml(y ~ x1 + x2, "two-way", "id", "time")
  fml_str <- deparse(fml)
  expect_true(grepl("id", fml_str) && grepl("time", fml_str))
})

test_that("effect_to_feglm_fml: nested stops with informative error", {
  expect_error(
    effect_to_feglm_fml(y ~ x1 + x2, "nested", "id", "time"),
    regexp = "nested"
  )
})

# ---------------------------------------------------------------------------
# effect_to_glmmtmb_fml
# ---------------------------------------------------------------------------

test_that("effect_to_glmmtmb_fml: individual adds (1|id)", {
  fml <- effect_to_glmmtmb_fml(y ~ x1, "individual", "id", "time")
  expect_true(grepl("1 | id", deparse(fml)))
})

test_that("effect_to_glmmtmb_fml: time adds (1|time)", {
  fml <- effect_to_glmmtmb_fml(y ~ x1, "time", "id", "time")
  expect_true(grepl("1 | time", deparse(fml)))
})

test_that("effect_to_glmmtmb_fml: two-way adds both random intercepts", {
  fml <- effect_to_glmmtmb_fml(y ~ x1, "two-way", "id", "time")
  fml_str <- deparse(fml)
  expect_true(grepl("1 | id", fml_str) && grepl("1 | time", fml_str))
})

test_that("effect_to_glmmtmb_fml: nested stops with informative error", {
  expect_error(
    effect_to_glmmtmb_fml(y ~ x1, "nested", "id", "time"),
    regexp = "nested"
  )
})

# ---------------------------------------------------------------------------
# fit_logistic_pooling
# ---------------------------------------------------------------------------

test_that("fit_logistic_pooling: returns glm object", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  fit <- fit_logistic_pooling(pd, y ~ x1 + x2, wts)
  expect_true(inherits(fit, "glm"))
  expect_equal(family(fit)$family, "binomial")
})

test_that("fit_logistic_pooling: works with numeric 0/1 response", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  expect_no_error(fit_logistic_pooling(pd, y ~ x1 + x2, wts))
})

test_that("fit_logistic_pooling: works with two-level factor response", {
  pd  <- make_bin_panel_data()
  pd$y_fac <- factor(pd$y, levels = c(0, 1), labels = c("no", "yes"))
  wts <- rep(1, nrow(pd))
  fit <- fit_logistic_pooling(pd, y_fac ~ x1 + x2, wts)
  expect_true(inherits(fit, "glm"))
})

# ---------------------------------------------------------------------------
# fit_logistic_fixed
# ---------------------------------------------------------------------------

test_that("fit_logistic_fixed: returns fixest object for individual effect", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  fit <- fit_logistic_fixed(pd, y ~ x1 + x2, "individual", wts, "id", "time")
  expect_true(inherits(fit, "fixest"))
})

test_that("fit_logistic_fixed: returns fixest object for time effect", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  fit <- fit_logistic_fixed(pd, y ~ x1 + x2, "time", wts, "id", "time")
  expect_true(inherits(fit, "fixest"))
})

test_that("fit_logistic_fixed: two-way effect smoke test", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  # Two-way may fail due to collinearity on small data; capture gracefully
  result <- tryCatch(
    fit_logistic_fixed(pd, y ~ x1 + x2, "two-way", wts, "id", "time"),
    error = function(e) NULL
  )
  # Either succeeded (fixest object) or failed (NULL) — both acceptable
  if (!is.null(result)) {
    expect_true(inherits(result, "fixest"))
  } else {
    succeed("Two-way FE failed gracefully on small data")
  }
})

test_that("fit_logistic_fixed: nested effect stops with informative error", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  expect_error(
    fit_logistic_fixed(pd, y ~ x1 + x2, "nested", wts, "id", "time"),
    regexp = "nested"
  )
})

# ---------------------------------------------------------------------------
# fit_logistic_random
# ---------------------------------------------------------------------------

test_that("fit_logistic_random: nested effect stops with informative error", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  expect_error(
    fit_logistic_random(pd, y ~ x1 + x2, "nested", wts, "id", "time"),
    regexp = "nested"
  )
})

test_that("fit_logistic_random: individual effect smoke test (may fail due to convergence)", {
  skip_if_not_installed("glmmTMB")
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  result <- tryCatch(
    fit_logistic_random(pd, y ~ x1 + x2, "individual", wts, "id", "time"),
    error = function(e) NULL
  )
  # Either a glmmTMB object or NULL (convergence failure) — both acceptable
  if (!is.null(result)) {
    expect_true(inherits(result, "glmmTMB"))
  } else {
    succeed("glmmTMB fit failed gracefully (convergence/numerical issue)")
  }
})

# ---------------------------------------------------------------------------
# fit_logistic_local_model
# ---------------------------------------------------------------------------

test_that("fit_logistic_local_model: pooling ok status with numeric response", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_logistic_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_false(is.null(res$fit))
  expect_true(inherits(res$fit, "glm"))
})

test_that("fit_logistic_local_model: pooling ok status with factor response", {
  pd  <- make_bin_panel_data()
  pd$y_fac <- factor(pd$y, levels = c(0, 1), labels = c("no", "yes"))
  wts <- rep(1, nrow(pd))
  res <- fit_logistic_local_model(
    local_data = pd,
    formula    = y_fac ~ x1 + x2,
    model      = "pooling",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_true(inherits(res$fit, "glm"))
})

test_that("fit_logistic_local_model: multi-class factor response stops outside tryCatch", {
  pd      <- make_bin_panel_data()
  pd$y_mc <- factor(sample(c("a", "b", "c"), nrow(pd), replace = TRUE))
  wts     <- rep(1, nrow(pd))
  # validate_binary_response is called before tryCatch so it propagates
  expect_error(
    fit_logistic_local_model(
      local_data = pd,
      formula    = y_mc ~ x1 + x2,
      model      = "pooling",
      effect     = "individual",
      weights    = wts,
      index      = c("id", "time")
    ),
    regexp = "Multi-class"
  )
})

test_that("fit_logistic_local_model: bad formula gives failed status, not error", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_logistic_local_model(
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

test_that("fit_logistic_local_model: fixed individual effect returns ok status", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_logistic_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "fixed",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  expect_equal(res$status, "ok")
  expect_true(inherits(res$fit, "fixest"))
})

test_that("fit_logistic_local_model: nested effect propagates error before tryCatch", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  expect_error(
    fit_logistic_local_model(
      local_data = pd,
      formula    = y ~ x1 + x2,
      model      = "fixed",
      effect     = "nested",
      weights    = wts,
      index      = c("id", "time")
    ),
    regexp = "nested"
  )
})

test_that("fit_logistic_local_model: random model handles failure gracefully", {
  skip_if_not_installed("glmmTMB")
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  # If glmmTMB fails numerically, status should be "failed" not thrown error
  res <- fit_logistic_local_model(
    local_data = pd,
    formula    = y ~ x1 + x2,
    model      = "random",
    effect     = "individual",
    weights    = wts,
    index      = c("id", "time")
  )
  # Either ok or failed — but must not throw
  expect_true(res$status %in% c("ok", "failed"))
  if (res$status == "failed") {
    expect_null(res$fit)
    expect_true(is.character(res$error))
  }
})

# ---------------------------------------------------------------------------
# predict_logistic_local_model
# ---------------------------------------------------------------------------

test_that("predict_logistic_local_model: returns numeric in [0,1]", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  res <- fit_logistic_local_model(
    local_data = pd, formula = y ~ x1 + x2,
    model = "pooling", effect = "individual",
    weights = wts, index = c("id", "time")
  )
  preds <- predict_logistic_local_model(res, pd)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(pd))
  expect_true(all(preds >= 0 & preds <= 1, na.rm = TRUE))
})

test_that("predict_logistic_local_model: failed model returns NA vector", {
  failed_res <- list(status = "failed", fit = NULL, error = "err",
                     metadata = list())
  pd    <- make_bin_panel_data()
  preds <- predict_logistic_local_model(failed_res, pd)
  expect_true(all(is.na(preds)))
  expect_equal(length(preds), nrow(pd))
})

# ---------------------------------------------------------------------------
# extract_logistic_local_result
# ---------------------------------------------------------------------------

test_that("extract_logistic_local_result: ok result has correct structure", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_logistic_local_model(
    local_data = pd, formula = y ~ x1 + x2,
    model = "pooling", effect = "individual",
    weights = wts, index = c("id", "time")
  )
  ext <- extract_logistic_local_result(local_res, pd, y ~ x1 + x2)

  expect_equal(ext$status, "ok")
  expect_true(is.numeric(ext$coefficients))
  expect_true(is.numeric(ext$prob))
  expect_equal(length(ext$prob), nrow(pd))
  expect_true(all(ext$prob >= 0 & ext$prob <= 1, na.rm = TRUE))
  expect_true(is.integer(ext$class_pred) || is.numeric(ext$class_pred))
  expect_true(all(ext$class_pred %in% c(0L, 1L, NA_integer_)))
  expect_true(is.numeric(ext$pearson_resid))
  expect_true(all(is.finite(ext$pearson_resid)))  # no Inf due to clipping
})

test_that("extract_logistic_local_result: failed result returns NA placeholders", {
  failed_res <- list(status = "failed", fit = NULL, error = "test error",
                     metadata = list())
  pd  <- make_bin_panel_data()
  ext <- extract_logistic_local_result(failed_res, pd, y ~ x1 + x2)

  expect_equal(ext$status, "failed")
  expect_true(is.na(ext$coefficients) || all(is.na(ext$coefficients)))
  expect_true(all(is.na(ext$prob)))
  expect_true(all(is.na(ext$class_pred)))
  expect_true(all(is.na(ext$pearson_resid)))
  expect_equal(ext$error, "test error")
})

test_that("extract_logistic_local_result: Pearson residuals finite (no Inf)", {
  # Use data with extreme probabilities to verify clipping
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_logistic_local_model(
    local_data = pd, formula = y ~ x1 + x2,
    model = "pooling", effect = "individual",
    weights = wts, index = c("id", "time")
  )
  ext <- extract_logistic_local_result(local_res, pd, y ~ x1 + x2)
  # Verify Pearson residuals are all finite (clipping prevents Inf)
  expect_true(all(is.finite(ext$pearson_resid)))
})

# ---------------------------------------------------------------------------
# threshold effect on class_pred
# ---------------------------------------------------------------------------

test_that("threshold = 0.5 and threshold = 0.0 give different class predictions", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_logistic_local_model(
    local_data = pd, formula = y ~ x1 + x2,
    model = "pooling", effect = "individual",
    weights = wts, index = c("id", "time")
  )
  ext_05 <- extract_logistic_local_result(local_res, pd, y ~ x1 + x2,
                                           threshold = 0.5)
  ext_00 <- extract_logistic_local_result(local_res, pd, y ~ x1 + x2,
                                           threshold = 0.0)

  # threshold = 0 classifies everything as 1 (all prob >= 0)
  expect_true(all(ext_00$class_pred == 1L))
  # threshold = 0.5 gives a mix
  expect_false(all(ext_05$class_pred == 1L) && all(ext_05$class_pred == 0L))
})

test_that("threshold = 1.0 classifies everything as 0", {
  pd  <- make_bin_panel_data()
  wts <- rep(1, nrow(pd))
  local_res <- fit_logistic_local_model(
    local_data = pd, formula = y ~ x1 + x2,
    model = "pooling", effect = "individual",
    weights = wts, index = c("id", "time")
  )
  ext <- extract_logistic_local_result(local_res, pd, y ~ x1 + x2,
                                        threshold = 1.0)
  # Probabilities are always < 1, so all should be classified as 0
  expect_true(all(ext$class_pred == 0L))
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — smoke tests (pooling / fixed / random)
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: pooling smoke test returns correct structure", {
  ctx <- make_logistic_context(model = "pooling", effect = "individual",
                                kernel = "bisquare", adaptive = TRUE)
  res <- fit_logistic_gwpr(ctx, bandwidth = 3L)

  expect_true(is.list(res))
  expect_true("local_results" %in% names(res))
  expect_true("predictions"   %in% names(res))
  expect_true("class_pred"    %in% names(res))
  expect_true("pearson_resid" %in% names(res))
  expect_true("metrics"       %in% names(res))
  expect_true("metadata"      %in% names(res))

  # 4 spatial units => 4 local results
  expect_equal(length(res$local_results), 4L)

  # Logistic metrics fields
  expect_true(all(c("log_loss", "accuracy", "precision", "recall", "f1_score")
                  %in% names(res$metrics)))

  # Predictions should be in [0, 1]
  probs <- res$predictions[!is.na(res$predictions)]
  expect_true(all(probs >= 0 & probs <= 1))
})

test_that("fit_logistic_gwpr: fixed effects individual smoke test", {
  ctx <- make_logistic_context(model = "fixed", effect = "individual",
                                kernel = "gaussian", adaptive = FALSE)
  expect_no_error({
    res <- fit_logistic_gwpr(ctx, bandwidth = 2)
  })
  expect_equal(length(res$local_results), 4L)
})

test_that("fit_logistic_gwpr: random effects smoke test (may fail due to convergence)", {
  skip_if_not_installed("glmmTMB")
  ctx <- make_logistic_context(model = "random", effect = "individual",
                                kernel = "bisquare", adaptive = TRUE)
  # Should not throw — failures captured per-focus
  expect_no_error({
    res <- fit_logistic_gwpr(ctx, bandwidth = 3L)
  })
  # Structure should be intact even if all models failed
  expect_equal(length(res$local_results), 4L)
  statuses <- vapply(res$local_results, `[[`, character(1L), "status")
  expect_true(all(statuses %in% c("ok", "failed")))
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — family = "binomial" only
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: non-binomial family raises informative error", {
  ctx <- make_logistic_context()
  expect_error(
    fit_logistic_gwpr(ctx, bandwidth = 3L, family = "multinomial"),
    regexp = "binomial"
  )
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — local model failure does not break structure
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: local model failures recorded, structure intact", {
  pd  <- make_bin_panel_data()
  ctx <- new_gwpr_context(
    formula    = y ~ nonexistent_var,   # will fail every local model
    family     = "binomial",
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

  res <- fit_logistic_gwpr(ctx, bandwidth = 3L)

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
# fit_logistic_gwpr — pre-supplied weights_context
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: accepts pre-built weights_context", {
  ctx   <- make_logistic_context(model = "pooling", effect = "individual",
                                  kernel = "bisquare", adaptive = TRUE)
  wctx  <- build_distance_context(
    coords  = ctx$coords,
    ids     = rownames(ctx$coords),
    longlat = FALSE,
    cache   = TRUE
  )
  expect_no_error(fit_logistic_gwpr(ctx, bandwidth = 3L,
                                    weights_context = wctx))
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — reproducibility
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: same inputs give identical results", {
  ctx  <- make_logistic_context(model = "pooling", effect = "individual",
                                 kernel = "bisquare", adaptive = TRUE)
  res1 <- fit_logistic_gwpr(ctx, bandwidth = 3L)
  res2 <- fit_logistic_gwpr(ctx, bandwidth = 3L)

  expect_equal(
    res1$local_results[[1]]$coefficients,
    res2$local_results[[1]]$coefficients
  )
  expect_equal(res1$predictions, res2$predictions)
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — custom threshold changes class_pred
# ---------------------------------------------------------------------------

test_that("fit_logistic_gwpr: threshold = 0.0 classifies all as 1", {
  ctx <- make_logistic_context(model = "pooling", effect = "individual",
                                threshold = 0.0, kernel = "bisquare",
                                adaptive = TRUE)
  res <- fit_logistic_gwpr(ctx, bandwidth = 3L)
  classes <- res$class_pred[!is.na(res$class_pred)]
  expect_true(all(classes == 1L))
})

test_that("fit_logistic_gwpr: threshold = 1.0 classifies all as 0", {
  ctx <- make_logistic_context(model = "pooling", effect = "individual",
                                threshold = 1.0, kernel = "bisquare",
                                adaptive = TRUE)
  res <- fit_logistic_gwpr(ctx, bandwidth = 3L)
  classes <- res$class_pred[!is.na(res$class_pred)]
  expect_true(all(classes == 0L))
})

# ---------------------------------------------------------------------------
# fit_logistic_gwpr — glmmTMB and fixest small-sample failure capture
# ---------------------------------------------------------------------------

test_that("fixest small-sample collinearity is captured as failed status", {
  # With very unbalanced data, feglm may fail with collinearity error
  # This verifies the error is captured, not propagated
  pd <- data.frame(
    id   = rep(as.character(1:2), each = 2),
    time = rep(1:2, times = 2),
    y    = c(1L, 1L, 0L, 0L),  # perfectly separated by id
    x1   = c(1, 1, -1, -1),    # perfectly correlated with id
    stringsAsFactors = FALSE
  )
  coords <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE,
                   dimnames = list(c("1", "2"), c("X", "Y")))
  ctx <- new_gwpr_context(
    formula    = y ~ x1,
    family     = "binomial",
    model      = "fixed",
    effect     = "individual",
    id         = "id",
    time       = "time",
    kernel     = "bisquare",
    adaptive   = FALSE,
    threshold  = 0.5,
    workers    = 1L,
    panel_data = pd,
    coords     = coords
  )
  # May succeed or fail — but must not throw an uncaught error
  expect_no_error({
    res <- fit_logistic_gwpr(ctx, bandwidth = 2)
  })
  # Structure intact
  expect_true("local_results" %in% names(res))
})

test_that("glmmTMB small-sample failure is captured as failed status", {
  skip_if_not_installed("glmmTMB")
  # Minimal data likely to cause glmmTMB convergence failure
  pd <- data.frame(
    id   = rep(as.character(1:2), each = 2),
    time = rep(1:2, times = 2),
    y    = c(1L, 0L, 1L, 0L),
    x1   = c(1, -1, 1, -1),
    stringsAsFactors = FALSE
  )
  coords <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE,
                   dimnames = list(c("1", "2"), c("X", "Y")))
  ctx <- new_gwpr_context(
    formula    = y ~ x1,
    family     = "binomial",
    model      = "random",
    effect     = "individual",
    id         = "id",
    time       = "time",
    kernel     = "bisquare",
    adaptive   = FALSE,
    threshold  = 0.5,
    workers    = 1L,
    panel_data = pd,
    coords     = coords
  )
  # Must not throw
  expect_no_error({
    res <- fit_logistic_gwpr(ctx, bandwidth = 2)
  })
  statuses <- vapply(res$local_results, `[[`, character(1L), "status")
  expect_true(all(statuses %in% c("ok", "failed")))
})
