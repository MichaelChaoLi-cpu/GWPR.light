#' @title Parallel Execution Module
#' @description Unified parallel execution interface for bandwidth search,
#'   local model fitting, and diagnostics. Shields backend differences and
#'   ensures CRAN-friendly behaviour.
#' @name parallel_module
NULL

#' Apply a function over a list with optional parallel execution
#'
#' A thin wrapper that uses plain \code{lapply} when \code{workers = 1} and
#' switches to \code{future.apply::future_lapply} with a \code{multisession}
#' plan for \code{workers > 1}. The global \code{future} plan is always
#' restored to \code{sequential} after the call, preventing side-effects.
#'
#' Worker-level errors are caught and returned as character strings (prefixed
#' with \code{"ERROR: "}) rather than aborting the entire call.
#'
#' @param x A list (or vector) of inputs to iterate over.
#' @param fn A function to apply to each element of \code{x}. The function
#'   receives the element as its first argument; additional arguments are
#'   passed via \code{...}.
#' @param workers Integer scalar. Number of parallel workers. \code{1}
#'   (default) uses serial \code{lapply} and never touches \code{future}.
#' @param seed Integer scalar or \code{NULL}. Random seed for reproducibility.
#'   In serial mode the R RNG is seeded with \code{set.seed(seed)}. In
#'   parallel mode the seed is forwarded to \code{future_lapply} via the
#'   \code{future.seed} argument using L'Ecuyer-CMRG streams.
#' @param ... Additional arguments forwarded to \code{fn}.
#' @param packages Character vector of package names that workers need to load,
#'   or \code{NULL} (default). Ignored in serial mode.
#'
#' @return A list of the same length as \code{x}. Elements where \code{fn}
#'   threw an error are replaced with a character string \code{"ERROR: <msg>"}.
#'
#' @examples
#' result <- parallel_map(1:3, function(x) x^2, workers = 1)
#' stopifnot(identical(result, list(1, 4, 9)))
#'
#' @export
parallel_map <- function(x, fn, workers = 1, seed = NULL, ..., packages = NULL) {
  stopifnot(is.numeric(workers), length(workers) == 1L, workers >= 1L)
  workers <- as.integer(workers)

  safe_fn <- function(elem) {
    tryCatch(fn(elem, ...), error = function(e) paste0("ERROR: ", conditionMessage(e)))
  }

  if (workers == 1L) {
    # Serial path — never touch future plan
    if (!is.null(seed)) set.seed(seed)
    return(lapply(x, safe_fn))
  }

  # Parallel path — requires future and future.apply
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for parallel_map() with workers > 1.")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for parallel_map() with workers > 1.")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = workers)

  future_seed <- if (is.null(seed)) FALSE else seed

  result <- future.apply::future_lapply(
    X            = x,
    FUN          = safe_fn,
    future.seed  = future_seed,
    future.packages = packages
  )

  result
}

#' Execute an expression with a reproducible seed
#'
#' Sets \code{set.seed(seed)} before evaluating \code{expr} and restores the
#' previous RNG state (kind, seed, normal.kind) afterwards.
#'
#' @param seed Integer scalar. Seed value passed to \code{set.seed}.
#' @param expr An R expression to evaluate (passed unevaluated; use
#'   \code{quote()} or wrap in \code{with_reproducible_seed(seed, { ... })}).
#'
#' @return The value of \code{expr}.
#'
#' @examples
#' r1 <- with_reproducible_seed(42, runif(3))
#' r2 <- with_reproducible_seed(42, runif(3))
#' stopifnot(identical(r1, r2))
#'
#' @export
with_reproducible_seed <- function(seed, expr) {
  # Save current RNG state
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }

  set.seed(seed)
  expr
}

#' Validate that parallel and serial results have the same structure
#'
#' Checks that \code{parallel_result} and \code{serial_result} have the same
#' length and the same element types. This function is intended for use in
#' tests to verify that switching to a parallel backend does not alter the
#' output structure.
#'
#' @param serial_result A list produced by \code{parallel_map} with
#'   \code{workers = 1}.
#' @param parallel_result A list produced by \code{parallel_map} with
#'   \code{workers > 1}.
#'
#' @return \code{TRUE} (invisibly) if the structures match; otherwise an
#'   informative error is thrown.
#'
#' @examples
#' sr <- parallel_map(1:3, function(x) x * 2, workers = 1)
#' # In real usage parallel_result would use workers > 1; here we reuse sr.
#' validate_parallel_result(sr, sr)
#'
#' @export
validate_parallel_result <- function(serial_result, parallel_result) {
  if (!is.list(serial_result))   stop("`serial_result` must be a list.")
  if (!is.list(parallel_result)) stop("`parallel_result` must be a list.")

  if (length(serial_result) != length(parallel_result)) {
    stop(sprintf(
      "Length mismatch: serial has %d elements, parallel has %d.",
      length(serial_result), length(parallel_result)
    ))
  }

  for (i in seq_along(serial_result)) {
    s_class <- class(serial_result[[i]])
    p_class <- class(parallel_result[[i]])
    if (!identical(s_class, p_class)) {
      stop(sprintf(
        "Element %d class mismatch: serial = '%s', parallel = '%s'.",
        i, paste(s_class, collapse = "/"), paste(p_class, collapse = "/")
      ))
    }
  }

  invisible(TRUE)
}
