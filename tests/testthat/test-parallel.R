# Tests for parallel.R
# All tests use workers = 1 (serial) to stay CRAN-friendly.

# ---------------------------------------------------------------------------
# parallel_map — basic serial behaviour
# ---------------------------------------------------------------------------

test_that("parallel_map with workers=1 returns correct results", {
  result <- parallel_map(1:5, function(x) x^2, workers = 1)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_equal(result[[1]], 1)
  expect_equal(result[[3]], 9)
  expect_equal(result[[5]], 25)
})

test_that("parallel_map with workers=1 passes extra arguments via ...", {
  result <- parallel_map(1:3, function(x, add) x + add, workers = 1, add = 10)
  expect_equal(result[[1]], 11)
  expect_equal(result[[2]], 12)
  expect_equal(result[[3]], 13)
})

test_that("parallel_map with workers=1 returns a list of same length as input", {
  input  <- list(a = 1, b = 2, c = 3)
  result <- parallel_map(input, function(x) x * 2, workers = 1)
  expect_length(result, 3)
})

# ---------------------------------------------------------------------------
# parallel_map — worker error capture
# ---------------------------------------------------------------------------

test_that("parallel_map captures worker errors without crashing (workers=1)", {
  fn <- function(x) {
    if (x == 2) stop("deliberate error")
    x * 10
  }
  result <- parallel_map(1:3, fn, workers = 1)
  expect_equal(result[[1]], 10)
  expect_match(result[[2]], "^ERROR: deliberate error")
  expect_equal(result[[3]], 30)
})

test_that("parallel_map worker error returns a character string", {
  fn <- function(x) stop("boom")
  result <- parallel_map(list(1), fn, workers = 1)
  expect_type(result[[1]], "character")
  expect_match(result[[1]], "^ERROR: ")
})

# ---------------------------------------------------------------------------
# parallel_map — workers=1 does not touch future plan
# ---------------------------------------------------------------------------

test_that("parallel_map with workers=1 does not load or change future plan", {
  skip_if_not_installed("future")
  plan_before <- class(future::plan())
  parallel_map(1:3, function(x) x, workers = 1)
  plan_after  <- class(future::plan())
  expect_identical(plan_before, plan_after)
})

# ---------------------------------------------------------------------------
# parallel_map — reproducibility with fixed seed (serial)
# ---------------------------------------------------------------------------

test_that("parallel_map with fixed seed gives reproducible results (workers=1)", {
  fn <- function(x) runif(1)
  r1 <- parallel_map(1:5, fn, workers = 1, seed = 123)
  r2 <- parallel_map(1:5, fn, workers = 1, seed = 123)
  expect_identical(r1, r2)
})

test_that("parallel_map with different seeds gives different results (workers=1)", {
  fn <- function(x) runif(1)
  r1 <- parallel_map(1:5, fn, workers = 1, seed = 1)
  r2 <- parallel_map(1:5, fn, workers = 1, seed = 2)
  # It is extremely unlikely all 5 values coincide for different seeds
  expect_false(identical(r1, r2))
})

# ---------------------------------------------------------------------------
# with_reproducible_seed
# ---------------------------------------------------------------------------

test_that("with_reproducible_seed gives reproducible output for same seed", {
  r1 <- with_reproducible_seed(42, runif(5))
  r2 <- with_reproducible_seed(42, runif(5))
  expect_equal(r1, r2)
})

test_that("with_reproducible_seed gives different output for different seeds", {
  r1 <- with_reproducible_seed(1, runif(5))
  r2 <- with_reproducible_seed(2, runif(5))
  expect_false(identical(r1, r2))
})

test_that("with_reproducible_seed restores RNG state afterwards", {
  set.seed(99)
  before <- runif(3)           # consume 3 draws from seed 99

  set.seed(99)
  with_reproducible_seed(42, runif(3))   # should not affect subsequent state
  after <- runif(3)            # next 3 draws from seed 99

  expect_equal(before, after)
})

test_that("with_reproducible_seed returns the value of expr", {
  val <- with_reproducible_seed(7, 1 + 1)
  expect_equal(val, 2)
})

# ---------------------------------------------------------------------------
# validate_parallel_result
# ---------------------------------------------------------------------------

test_that("validate_parallel_result returns TRUE invisibly for matching lists", {
  sr <- parallel_map(1:4, function(x) x^2, workers = 1)
  result <- validate_parallel_result(sr, sr)
  expect_true(result)
})

test_that("validate_parallel_result errors on length mismatch", {
  sr <- list(1, 2, 3)
  pr <- list(1, 2)
  expect_error(validate_parallel_result(sr, pr), "Length mismatch")
})

test_that("validate_parallel_result errors on class mismatch", {
  sr <- list(1L, 2L)
  pr <- list("a", "b")
  expect_error(validate_parallel_result(sr, pr), "class mismatch")
})

test_that("validate_parallel_result errors when inputs are not lists", {
  expect_error(validate_parallel_result(1:3, list(1, 2, 3)), "must be a list")
  expect_error(validate_parallel_result(list(1, 2, 3), 1:3), "must be a list")
})

# ---------------------------------------------------------------------------
# Structure consistency: serial result matches expected output shape
# ---------------------------------------------------------------------------

test_that("parallel_map serial result structure is consistent across calls", {
  fn <- function(x) list(value = x, squared = x^2)
  r1 <- parallel_map(1:3, fn, workers = 1)
  r2 <- parallel_map(1:3, fn, workers = 1)
  # Validate structure
  expect_true(validate_parallel_result(r1, r2))
  # Check content
  for (i in 1:3) {
    expect_named(r1[[i]], c("value", "squared"))
    expect_equal(r1[[i]]$squared, i^2)
  }
})
