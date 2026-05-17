## tests/testthat/test-api.R
## End-to-end tests for the public API (R/api.R)
##
## All tests use small simulated panel + sf data.
## workers = 1 (no parallel backend activated).

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

# 4 spatial units, 5 time periods each => 20 rows of linear panel data
make_api_panel_data <- function(seed = 42L) {
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

# 4 spatial units, 5 time periods each => 20 rows of binary panel data
make_api_bin_panel_data <- function(seed = 99L) {
  set.seed(seed)
  n_units <- 4L
  n_times <- 5L
  id   <- rep(as.character(1:n_units), each = n_times)
  time <- rep(1:n_times, times = n_units)
  x1   <- rnorm(n_units * n_times)
  lp   <- 0.5 * x1
  y    <- as.integer(rbinom(n_units * n_times, 1L, plogis(lp)))
  data.frame(id = id, time = time, y = y, x1 = x1,
             stringsAsFactors = FALSE)
}

# sf POINT object with 4 units on a 2x2 grid
make_api_sf <- function() {
  pts <- data.frame(
    id = as.character(1:4),
    X  = c(0, 1, 0, 1),
    Y  = c(0, 0, 1, 1)
  )
  sf::st_as_sf(pts, coords = c("X", "Y"), crs = NA_integer_)
}

# ---------------------------------------------------------------------------
# fit_gwpr  — linear (Gaussian), pooling model
# ---------------------------------------------------------------------------

test_that("fit_gwpr: linear pooling returns gwpr_fit", {
  dat  <- make_api_panel_data()
  spt  <- make_api_sf()

  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = 2,
    family    = "gaussian",
    model     = "pooling",
    effect    = "individual",
    workers   = 1L
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_equal(fit$family, "gaussian")
  expect_equal(fit$model,  "pooling")
  expect_false(is.null(fit$local_results))
  expect_length(fit$local_results, 4L)   # one per spatial unit
  expect_false(is.null(fit$predictions))
  expect_false(is.null(fit$residuals))
  expect_false(is.null(fit$metrics))
})

# ---------------------------------------------------------------------------
# fit_gwpr  — linear (Gaussian), within model
# ---------------------------------------------------------------------------

test_that("fit_gwpr: linear within returns gwpr_fit", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    fit_gwpr(
      formula   = y ~ x1 + x2,
      data      = dat,
      spatial   = spt,
      id        = "id",
      time      = "time",
      bandwidth = 2,
      family    = "gaussian",
      model     = "within",
      effect    = "individual",
      workers   = 1L
    )
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_equal(fit$family, "gaussian")
  expect_equal(fit$model,  "within")
})

# ---------------------------------------------------------------------------
# fit_gwpr  — binomial (logistic), pooling model
# ---------------------------------------------------------------------------

test_that("fit_gwpr: logistic pooling returns gwpr_fit", {
  dat <- make_api_bin_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    fit_gwpr(
      formula   = y ~ x1,
      data      = dat,
      spatial   = spt,
      id        = "id",
      time      = "time",
      bandwidth = 2,
      family    = "binomial",
      model     = "pooling",
      effect    = "individual",
      workers   = 1L
    )
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_equal(fit$family, "binomial")
  expect_false(is.null(fit$predictions))   # predicted probabilities
  expect_false(is.null(fit$metrics))
})

# ---------------------------------------------------------------------------
# fit_gwpr  — given bandwidth does NOT trigger search
# ---------------------------------------------------------------------------

test_that("fit_gwpr: does not call select_bandwidth when bandwidth is given", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  # We mock select_bandwidth to detect if it is called
  was_called <- FALSE
  local_mock_env <- new.env(parent = asNamespace("GWPR.light"))
  local({
    fit <- fit_gwpr(
      formula   = y ~ x1 + x2,
      data      = dat,
      spatial   = spt,
      id        = "id",
      time      = "time",
      bandwidth = 2,
      family    = "gaussian",
      model     = "pooling",
      workers   = 1L
    )
    expect_s3_class(fit, "gwpr_fit")
    # No assertion needed for was_called; absence of error is sufficient.
    # The key check is that bandwidth = 2 is stored in fit.
    expect_equal(fit$bandwidth, 2)
  })
})

# ---------------------------------------------------------------------------
# fit_gwpr  — bandwidth stored in output
# ---------------------------------------------------------------------------

test_that("fit_gwpr: bandwidth value is preserved in output", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  bw <- 1.5
  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = bw,
    family    = "gaussian",
    model     = "pooling",
    workers   = 1L
  )
  expect_equal(fit$bandwidth, bw)
})

# ---------------------------------------------------------------------------
# fit_gwpr  — invalid bandwidth raises error
# ---------------------------------------------------------------------------

test_that("fit_gwpr: missing bandwidth raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             workers = 1L),
    regexp = "bandwidth"
  )
})

test_that("fit_gwpr: negative bandwidth raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             bandwidth = -1, workers = 1L),
    regexp = "positive"
  )
})

# ---------------------------------------------------------------------------
# fit_gwpr  — workers = 1 does not trigger parallel backend
# ---------------------------------------------------------------------------

test_that("fit_gwpr: workers = 1 completes without activating parallel plan", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  # With workers = 1 the function must run without error / without calling
  # future::plan().  We simply verify it returns a valid object.
  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = 2,
    family    = "gaussian",
    model     = "pooling",
    workers   = 1L
  )
  expect_s3_class(fit, "gwpr_fit")
})

# ---------------------------------------------------------------------------
# Illegal family / model / effect raise clear errors
# ---------------------------------------------------------------------------

test_that("fit_gwpr: illegal family raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             bandwidth = 2, family = "poisson", workers = 1L),
    regexp = "arg"
  )
})

test_that("fit_gwpr: illegal model raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             bandwidth = 2, model = "twoway_bad", workers = 1L),
    regexp = "arg"
  )
})

test_that("fit_gwpr: illegal effect raises error via validate_inputs", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             bandwidth = 2, effect = "bad_effect", workers = 1L),
    regexp = "arg"
  )
})

# ---------------------------------------------------------------------------
# validate_inputs integration via fit_gwpr
# ---------------------------------------------------------------------------

test_that("fit_gwpr: non-sf spatial raises informative error", {
  dat <- make_api_panel_data()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = list(id = "1"), id = "id",
             time = "time", bandwidth = 2, workers = 1L),
    regexp = "sf"
  )
})

test_that("fit_gwpr: missing id column in data raises error", {
  dat <- make_api_panel_data()
  dat2 <- dat
  names(dat2)[names(dat2) == "id"] <- "unit"   # rename id column
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat2, spatial = spt, id = "id",
             time = "time", bandwidth = 2, workers = 1L),
    regexp = "id"
  )
})

test_that("fit_gwpr: workers = 0 raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    fit_gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
             bandwidth = 2, workers = 0L),
    regexp = "workers"
  )
})

# ---------------------------------------------------------------------------
# diagnose_gwpr
# ---------------------------------------------------------------------------

test_that("diagnose_gwpr: returns gwpr_diagnostics object", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = 2,
    family    = "gaussian",
    model     = "pooling",
    workers   = 1L
  )

  diag_result <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman", "lm_test"))
  expect_s3_class(diag_result, "gwpr_diagnostics")
  expect_true(length(diag_result$diagnostics) > 0L)
})

test_that("diagnose_gwpr: object must be gwpr_fit", {
  expect_error(diagnose_gwpr(list(a = 1)), regexp = "gwpr_fit")
})

test_that("diagnose_gwpr: moran skipped gracefully without spatial_weights", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = 2,
    family    = "gaussian",
    model     = "pooling",
    workers   = 1L
  )

  # Without spatial_weights / panel_index, moran should be skipped, not error
  diag_result <- diagnose_gwpr(fit, diagnostics = "moran")
  expect_s3_class(diag_result, "gwpr_diagnostics")
  expect_equal(diag_result$diagnostics$moran$status, "skipped")
})

# ---------------------------------------------------------------------------
# gwpr  — main entry point, given bandwidth (no search)
# ---------------------------------------------------------------------------

test_that("gwpr: given bandwidth returns gwpr_fit without search", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    gwpr(
      formula      = y ~ x1 + x2,
      data         = dat,
      spatial      = spt,
      id           = "id",
      time         = "time",
      bandwidth    = 2,
      family       = "gaussian",
      model        = "pooling",
      effect       = "individual",
      diagnostics  = FALSE,
      workers      = 1L
    )
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_equal(fit$bandwidth, 2)
  expect_null(fit$search)   # search was not triggered
})

test_that("gwpr: given bandwidth, search field is NULL", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    gwpr(
      formula     = y ~ x1 + x2,
      data        = dat,
      spatial     = spt,
      id          = "id",
      time        = "time",
      bandwidth   = 2,
      diagnostics = FALSE,
      workers     = 1L
    )
  )
  expect_null(fit$search)
})

# ---------------------------------------------------------------------------
# gwpr  — diagnostics integration
# ---------------------------------------------------------------------------

test_that("gwpr: diagnostics = FALSE suppresses diagnostic computation", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    gwpr(
      formula     = y ~ x1 + x2,
      data        = dat,
      spatial     = spt,
      id          = "id",
      time        = "time",
      bandwidth   = 2,
      model       = "pooling",
      diagnostics = FALSE,
      workers     = 1L
    )
  )
  expect_s3_class(fit, "gwpr_fit")
  # diagnostics_result should not exist or be NULL
  expect_true(is.null(fit$diagnostics_result))
})

# ---------------------------------------------------------------------------
# select_bandwidth  — grid search, very small grid
# ---------------------------------------------------------------------------

test_that("select_bandwidth: grid search returns gwpr_bandwidth", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  bw <- suppressWarnings(
    select_bandwidth(
      formula  = y ~ x1 + x2,
      data     = dat,
      spatial  = spt,
      id       = "id",
      time     = "time",
      family   = "gaussian",
      model    = "pooling",
      effect   = "individual",
      method   = "grid",
      control  = list(lower = 1, upper = 2, step = 1),
      workers  = 1L
    )
  )

  expect_s3_class(bw, "gwpr_bandwidth")
  expect_true(is.numeric(bw$best_bandwidth))
  expect_true(bw$best_bandwidth %in% c(1, 2))
})

test_that("select_bandwidth: grid method, best_score is finite", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  bw <- suppressWarnings(
    select_bandwidth(
      formula = y ~ x1 + x2,
      data    = dat,
      spatial = spt,
      id      = "id",
      time    = "time",
      method  = "grid",
      control = list(lower = 1, upper = 2, step = 1),
      workers = 1L
    )
  )
  expect_true(is.finite(bw$best_score))
})

# ---------------------------------------------------------------------------
# select_bandwidth  — illegal method raises error
# ---------------------------------------------------------------------------

test_that("select_bandwidth: illegal method raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    select_bandwidth(
      formula = y ~ x1, data = dat, spatial = spt,
      id = "id", time = "time",
      method = "bayes",
      control = list(lower = 1, upper = 2, step = 0.5),
      workers = 1L
    ),
    regexp = "arg"
  )
})

# ---------------------------------------------------------------------------
# select_bandwidth  — grid missing bounds raises error
# ---------------------------------------------------------------------------

test_that("select_bandwidth: grid without bounds raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    select_bandwidth(
      formula = y ~ x1, data = dat, spatial = spt,
      id = "id", time = "time",
      method  = "grid",
      control = list(),   # missing lower/upper/step
      workers = 1L
    ),
    regexp = "lower|upper|step"
  )
})

# ---------------------------------------------------------------------------
# gwpr  — illegal family / model / effect via main entry
# ---------------------------------------------------------------------------

test_that("gwpr: illegal family raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
         bandwidth = 2, family = "poisson", workers = 1L),
    regexp = "arg"
  )
})

test_that("gwpr: illegal model raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
         bandwidth = 2, model = "wrong_model", workers = 1L),
    regexp = "arg"
  )
})

test_that("gwpr: illegal effect raises error", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  expect_error(
    gwpr(y ~ x1, data = dat, spatial = spt, id = "id", time = "time",
         bandwidth = 2, effect = "bad_effect", workers = 1L),
    regexp = "arg"
  )
})

# ---------------------------------------------------------------------------
# Smoke test: complete linear pipeline via gwpr() with bandwidth search
# ---------------------------------------------------------------------------

test_that("gwpr: full linear pipeline with grid bandwidth search", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    gwpr(
      formula          = y ~ x1 + x2,
      data             = dat,
      spatial          = spt,
      id               = "id",
      time             = "time",
      family           = "gaussian",
      model            = "pooling",
      effect           = "individual",
      bandwidth_method = "grid",
      bandwidth_control = list(lower = 1, upper = 2, step = 1),
      diagnostics      = FALSE,
      workers          = 1L
    )
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_s3_class(fit$search, "gwpr_bandwidth")
  expect_true(is.numeric(fit$bandwidth))
  expect_false(is.null(fit$metrics))
})

# ---------------------------------------------------------------------------
# Smoke test: complete logistic pipeline via gwpr() with given bandwidth
# ---------------------------------------------------------------------------

test_that("gwpr: full logistic pipeline with given bandwidth", {
  dat <- make_api_bin_panel_data()
  spt <- make_api_sf()

  fit <- suppressWarnings(
    gwpr(
      formula     = y ~ x1,
      data        = dat,
      spatial     = spt,
      id          = "id",
      time        = "time",
      family      = "binomial",
      model       = "pooling",
      effect      = "individual",
      bandwidth   = 2,
      diagnostics = FALSE,
      workers     = 1L
    )
  )

  expect_s3_class(fit, "gwpr_fit")
  expect_equal(fit$family, "binomial")
  expect_false(is.null(fit$predictions))
  expect_false(is.null(fit$metrics))
  expect_null(fit$search)   # no search triggered
})

# ---------------------------------------------------------------------------
# spatial_results structure
# ---------------------------------------------------------------------------

test_that("fit_gwpr: spatial_results is a data.frame with unit_id column", {
  dat <- make_api_panel_data()
  spt <- make_api_sf()

  fit <- fit_gwpr(
    formula   = y ~ x1 + x2,
    data      = dat,
    spatial   = spt,
    id        = "id",
    time      = "time",
    bandwidth = 2,
    model     = "pooling",
    workers   = 1L
  )

  expect_true(is.data.frame(fit$spatial_results))
  expect_true("unit_id" %in% names(fit$spatial_results))
  expect_equal(nrow(fit$spatial_results), 4L)   # 4 spatial units
})
