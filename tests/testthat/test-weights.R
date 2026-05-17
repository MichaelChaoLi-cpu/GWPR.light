## tests/testthat/test-weights.R
## Tests for the Distance and Kernel Weights module (R/weights.R)
##
## All tests use small simulated data.  workers = 1 (no parallelism).

# ---------------------------------------------------------------------------
# Helper: small coordinate set (5 points on a unit grid)
# ---------------------------------------------------------------------------
make_coords <- function() {
  matrix(
    c(0, 0,
      1, 0,
      2, 0,
      0, 1,
      1, 1),
    ncol = 2, byrow = TRUE,
    dimnames = list(NULL, c("X", "Y"))
  )
}

# ===========================================================================
# validate_bandwidth
# ===========================================================================

test_that("validate_bandwidth: fixed bandwidth must be positive", {
  expect_true(validate_bandwidth(1.5, adaptive = FALSE))
  expect_error(validate_bandwidth(0,   adaptive = FALSE), regexp = "strictly positive")
  expect_error(validate_bandwidth(-1,  adaptive = FALSE), regexp = "strictly positive")
})

test_that("validate_bandwidth: adaptive bandwidth must be a positive integer", {
  expect_true(validate_bandwidth(3,   adaptive = TRUE))
  expect_error(validate_bandwidth(0,   adaptive = TRUE), regexp = "positive integer")
  expect_error(validate_bandwidth(-2,  adaptive = TRUE), regexp = "positive integer")
  expect_error(validate_bandwidth(2.5, adaptive = TRUE), regexp = "positive integer")
})

test_that("validate_bandwidth: non-numeric input triggers error", {
  expect_error(validate_bandwidth("big",  adaptive = FALSE), regexp = "single numeric")
  expect_error(validate_bandwidth(c(1,2), adaptive = FALSE), regexp = "single numeric")
})

test_that("validate_bandwidth: non-finite input triggers error", {
  expect_error(validate_bandwidth(Inf, adaptive = FALSE), regexp = "finite")
  expect_error(validate_bandwidth(NA_real_, adaptive = FALSE), regexp = "finite")
})

# ===========================================================================
# compute_distance
# ===========================================================================

test_that("compute_distance: Euclidean distance matrix has correct dimensions", {
  coords <- make_coords()
  d <- compute_distance(coords, longlat = FALSE)
  expect_true(is.matrix(d))
  expect_equal(dim(d), c(5L, 5L))
})

test_that("compute_distance: diagonal is zero", {
  coords <- make_coords()
  d <- compute_distance(coords, longlat = FALSE)
  expect_equal(unname(diag(d)), rep(0, 5))
})

test_that("compute_distance: matrix is symmetric", {
  coords <- make_coords()
  d <- compute_distance(coords, longlat = FALSE)
  expect_equal(d, t(d))
})

test_that("compute_distance: known Euclidean distances are correct", {
  # Point (0,0) and (1,0) -> distance 1
  # Point (0,0) and (1,1) -> distance sqrt(2)
  coords <- matrix(c(0,0, 1,0, 1,1), ncol = 2, byrow = TRUE)
  d <- compute_distance(coords, longlat = FALSE)
  expect_equal(d[1, 2], 1, tolerance = 1e-10)
  expect_equal(d[1, 3], sqrt(2), tolerance = 1e-10)
})

test_that("compute_distance: great-circle distance is non-negative and symmetric", {
  # Two known cities (approximate lon/lat)
  # Tokyo ~(139.7, 35.7), New York ~(-74.0, 40.7)
  coords <- matrix(c(139.7, 35.7, -74.0, 40.7), ncol = 2, byrow = TRUE)
  d <- compute_distance(coords, longlat = TRUE)
  expect_equal(dim(d), c(2L, 2L))
  expect_equal(d[1, 1], 0, tolerance = 1e-10)
  expect_true(d[1, 2] > 0)
  expect_equal(d[1, 2], d[2, 1], tolerance = 1e-6)
  # Rough sanity: Tokyo-NY is ~10800 km
  expect_true(d[1, 2] > 10000 && d[1, 2] < 12000)
})

test_that("compute_distance: accepts data.frame and coerces to matrix", {
  df <- as.data.frame(make_coords())
  expect_no_error(compute_distance(df, longlat = FALSE))
})

# ===========================================================================
# get_local_distances
# ===========================================================================

test_that("get_local_distances: returns correct vector for integer focus_id", {
  coords <- make_coords()
  d <- compute_distance(coords)
  v <- get_local_distances(d, focus_id = 1L)
  expect_equal(length(v), 5L)
  expect_equal(v[1], 0)           # focal unit distance to itself
  expect_equal(v[2], d[1, 2])
})

test_that("get_local_distances: works with character focus_id via row names", {
  coords <- make_coords()
  d <- compute_distance(coords)
  rownames(d) <- colnames(d) <- c("A", "B", "C", "D", "E")
  v <- get_local_distances(d, focus_id = "C")
  expect_equal(length(v), 5L)
  expect_equal(v[3], 0)           # C to itself
  expect_equal(v[1], d["C", "A"])
})

test_that("get_local_distances: works with list-wrapped distance context", {
  coords <- make_coords()
  ctx <- build_distance_context(coords, ids = letters[1:5])
  v <- get_local_distances(ctx, focus_id = 1L)
  expect_equal(length(v), 5L)
  expect_equal(v[1], 0)
})

test_that("get_local_distances: out-of-range focus_id triggers error", {
  coords <- make_coords()
  d <- compute_distance(coords)
  expect_error(get_local_distances(d, focus_id = 0L), regexp = "out of range")
  expect_error(get_local_distances(d, focus_id = 6L), regexp = "out of range")
})

test_that("get_local_distances: unknown character focus_id triggers error", {
  coords <- make_coords()
  d <- compute_distance(coords)
  rownames(d) <- colnames(d) <- c("A", "B", "C", "D", "E")
  expect_error(get_local_distances(d, focus_id = "Z"), regexp = "not found")
})

# ===========================================================================
# build_distance_context
# ===========================================================================

test_that("build_distance_context: returns correct class and fields", {
  coords <- make_coords()
  ctx <- build_distance_context(coords, ids = as.character(1:5))
  expect_s3_class(ctx, "gwpr_distance_context")
  expect_true(!is.null(ctx$distance_matrix))
  expect_equal(dim(ctx$distance_matrix), c(5L, 5L))
  expect_equal(ctx$ids, as.character(1:5))
})

test_that("build_distance_context: cache=FALSE leaves distance_matrix NULL", {
  coords <- make_coords()
  ctx <- build_distance_context(coords, ids = as.character(1:5), cache = FALSE)
  expect_null(ctx$distance_matrix)
})

test_that("build_distance_context: row/column names match ids", {
  coords <- make_coords()
  ids <- c("p1", "p2", "p3", "p4", "p5")
  ctx <- build_distance_context(coords, ids = ids)
  expect_equal(rownames(ctx$distance_matrix), ids)
  expect_equal(colnames(ctx$distance_matrix), ids)
})

test_that("build_distance_context: id length mismatch triggers error", {
  coords <- make_coords()
  expect_error(build_distance_context(coords, ids = 1:3), regexp = "length")
})

# ===========================================================================
# compute_kernel_weights — bisquare
# ===========================================================================

test_that("bisquare kernel: weight at d=0 equals 1", {
  w <- compute_kernel_weights(c(0, 0.5, 1, 1.5), bandwidth = 1, kernel = "bisquare", adaptive = FALSE)
  expect_equal(w[1], 1)
})

test_that("bisquare kernel: weight is 0 beyond bandwidth", {
  w <- compute_kernel_weights(c(0, 0.5, 1.01), bandwidth = 1, kernel = "bisquare", adaptive = FALSE)
  expect_equal(w[3], 0)
})

test_that("bisquare kernel: weight at boundary (d==bw) equals 0", {
  w <- compute_kernel_weights(c(1), bandwidth = 1, kernel = "bisquare", adaptive = FALSE)
  expect_equal(w[1], 0)
})

test_that("bisquare kernel: known mid-point value", {
  # d = 0.5, bw = 1: (1 - 0.25)^2 = 0.5625
  w <- compute_kernel_weights(0.5, bandwidth = 1, kernel = "bisquare", adaptive = FALSE)
  expect_equal(w, 0.5625, tolerance = 1e-10)
})

# ===========================================================================
# compute_kernel_weights — gaussian
# ===========================================================================

test_that("gaussian kernel: weight at d=0 equals 1", {
  w <- compute_kernel_weights(0, bandwidth = 1, kernel = "gaussian", adaptive = FALSE)
  expect_equal(w, 1)
})

test_that("gaussian kernel: weight is positive for moderate distances", {
  # Use distances that don't underflow to 0 (exp(-0.5*(5/2)^2) is still positive)
  d <- c(0, 1, 2, 3)
  w <- compute_kernel_weights(d, bandwidth = 5, kernel = "gaussian", adaptive = FALSE)
  expect_true(all(w > 0))
})

test_that("gaussian kernel: weight decreases with distance", {
  d <- c(0, 1, 2, 3)
  w <- compute_kernel_weights(d, bandwidth = 2, kernel = "gaussian", adaptive = FALSE)
  expect_true(all(diff(w) < 0))
})

test_that("gaussian kernel: known value at d=bw", {
  # exp(-0.5 * 1^2) = exp(-0.5)
  w <- compute_kernel_weights(1, bandwidth = 1, kernel = "gaussian", adaptive = FALSE)
  expect_equal(w, exp(-0.5), tolerance = 1e-10)
})

# ===========================================================================
# compute_kernel_weights — exponential
# ===========================================================================

test_that("exponential kernel: weight at d=0 equals 1", {
  w <- compute_kernel_weights(0, bandwidth = 2, kernel = "exponential", adaptive = FALSE)
  expect_equal(w, 1)
})

test_that("exponential kernel: weight is positive everywhere", {
  w <- compute_kernel_weights(c(0, 1, 10, 100), bandwidth = 2,
                              kernel = "exponential", adaptive = FALSE)
  expect_true(all(w > 0))
})

test_that("exponential kernel: known value at d=bw", {
  # exp(-1) when d = bw
  w <- compute_kernel_weights(3, bandwidth = 3, kernel = "exponential", adaptive = FALSE)
  expect_equal(w, exp(-1), tolerance = 1e-10)
})

# ===========================================================================
# compute_kernel_weights — tricube
# ===========================================================================

test_that("tricube kernel: weight at d=0 equals 1", {
  w <- compute_kernel_weights(0, bandwidth = 1, kernel = "tricube", adaptive = FALSE)
  expect_equal(w, 1)
})

test_that("tricube kernel: weight is 0 beyond bandwidth", {
  w <- compute_kernel_weights(1.01, bandwidth = 1, kernel = "tricube", adaptive = FALSE)
  expect_equal(w, 0)
})

test_that("tricube kernel: known mid-point value", {
  # d = 0.5, bw = 1: (1 - 0.5^3)^3 = (1 - 0.125)^3 = 0.875^3
  w <- compute_kernel_weights(0.5, bandwidth = 1, kernel = "tricube", adaptive = FALSE)
  expect_equal(w, 0.875^3, tolerance = 1e-10)
})

# ===========================================================================
# compute_kernel_weights — boxcar
# ===========================================================================

test_that("boxcar kernel: weight is 1 inside bandwidth", {
  w <- compute_kernel_weights(c(0, 0.5, 1), bandwidth = 1, kernel = "boxcar", adaptive = FALSE)
  expect_equal(w, c(1, 1, 1))
})

test_that("boxcar kernel: weight is 0 outside bandwidth", {
  w <- compute_kernel_weights(1.01, bandwidth = 1, kernel = "boxcar", adaptive = FALSE)
  expect_equal(w, 0)
})

# ===========================================================================
# fixed bandwidth — general correctness
# ===========================================================================

test_that("fixed bandwidth: all weights are non-negative", {
  d <- c(0, 0.5, 1, 2, 3)
  for (ker in c("bisquare", "gaussian", "exponential", "tricube", "boxcar")) {
    w <- compute_kernel_weights(d, bandwidth = 2, kernel = ker, adaptive = FALSE)
    expect_true(all(w >= 0), label = paste("kernel:", ker))
  }
})

test_that("fixed bandwidth: focal unit (d=0) always has maximum weight", {
  d <- c(0, 0.5, 1, 2)
  for (ker in c("bisquare", "gaussian", "exponential", "tricube", "boxcar")) {
    w <- compute_kernel_weights(d, bandwidth = 2, kernel = ker, adaptive = FALSE)
    expect_equal(which.max(w), 1L, label = paste("kernel:", ker))
  }
})

# ===========================================================================
# adaptive bandwidth
# ===========================================================================

test_that("adaptive bandwidth: uses k-th distance as effective bw", {
  # Distances from focal unit to 4 others: 1, 2, 3, 4
  # With k=2 the effective bw should be 2 (the 2nd smallest distance incl. self=0)
  # Sorted: 0, 1, 2, 3, 4  => sorted[2] = 1
  d <- c(0, 1, 2, 3, 4)
  # k=2: effective bw = sorted_d[2] = 1
  w_adaptive <- compute_kernel_weights(d, bandwidth = 2L, kernel = "bisquare", adaptive = TRUE)
  # Manually: effective bw = 1, bisquare at d=0: 1, d=1: 0, d>1: 0
  expect_equal(w_adaptive[1], 1)
  expect_equal(w_adaptive[2], 0)  # d=1 == bw -> weight 0
  expect_equal(w_adaptive[3], 0)  # d=2 > bw
})

test_that("adaptive bandwidth: k exceeding n issues warning and uses all units", {
  d <- c(0, 1, 2)
  expect_warning(
    w <- compute_kernel_weights(d, bandwidth = 10L, kernel = "gaussian", adaptive = TRUE),
    regexp = "exceeds number"
  )
  expect_equal(length(w), 3L)
})

test_that("adaptive bandwidth: result is consistent for bisquare", {
  # 5 units at distances 0,1,2,3,4; k=3 => effective bw = sorted[3] = 2
  d <- c(0, 1, 2, 3, 4)
  w <- compute_kernel_weights(d, bandwidth = 3L, kernel = "bisquare", adaptive = TRUE)
  # bw = 2; d=0: (1-0)^2 = 1; d=1: (1-0.25)^2=0.5625; d=2: 0; d=3: 0; d=4: 0
  expect_equal(w[1], 1,      tolerance = 1e-10)
  expect_equal(w[2], 0.5625, tolerance = 1e-10)
  expect_equal(w[3], 0,      tolerance = 1e-10)
  expect_equal(w[4], 0,      tolerance = 1e-10)
})

test_that("adaptive bandwidth: non-integer bandwidth triggers error", {
  d <- c(0, 1, 2, 3)
  expect_error(
    compute_kernel_weights(d, bandwidth = 2.5, kernel = "bisquare", adaptive = TRUE),
    regexp = "positive integer"
  )
})

# ===========================================================================
# distance matrix order / ID mapping consistency
# ===========================================================================

test_that("distance matrix row/column order matches input coordinate order", {
  coords <- matrix(c(0, 0,
                     10, 0,
                     0, 10), ncol = 2, byrow = TRUE)
  d <- compute_distance(coords, longlat = FALSE)
  # d[1,2] should be 10 (horizontal), d[1,3] = 10 (vertical), d[2,3] = sqrt(200)
  expect_equal(d[1, 2], 10, tolerance = 1e-10)
  expect_equal(d[1, 3], 10, tolerance = 1e-10)
  expect_equal(d[2, 3], sqrt(200), tolerance = 1e-10)
})

test_that("build_distance_context: distance_matrix rows match ids ordering", {
  coords <- matrix(c(0, 0,
                     5, 0,
                     0, 5), ncol = 2, byrow = TRUE)
  ids <- c("alpha", "beta", "gamma")
  ctx <- build_distance_context(coords, ids = ids)
  # Distance from alpha to beta: 5
  expect_equal(ctx$distance_matrix["alpha", "beta"], 5, tolerance = 1e-10)
  # Distance from beta to gamma: sqrt(50)
  expect_equal(ctx$distance_matrix["beta", "gamma"], sqrt(50), tolerance = 1e-10)
})

test_that("get_local_distances: integer and character access return same values", {
  coords <- make_coords()
  ids <- paste0("u", 1:5)
  ctx <- build_distance_context(coords, ids = ids)
  v_int  <- get_local_distances(ctx, focus_id = 3L)
  v_char <- get_local_distances(ctx$distance_matrix, focus_id = "u3")
  expect_equal(v_int, v_char)
})

# ===========================================================================
# illegal kernel name
# ===========================================================================

test_that("compute_kernel_weights: unknown kernel name triggers error", {
  expect_error(
    compute_kernel_weights(c(0, 1), bandwidth = 1, kernel = "parabola", adaptive = FALSE),
    regexp = "should be one of"
  )
})
