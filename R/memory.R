#' @title Memory Estimation Module for GWPR.light 1.0.0
#'
#' @description
#' Functions for estimating memory usage before running GWPR models.
#' The module provides warnings about memory risk levels to help users
#' avoid out-of-memory errors. These functions only warn; they never
#' stop execution.
#'
#' @name memory
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# estimate_memory
# ---------------------------------------------------------------------------

#' Estimate memory usage for a GWPR run
#'
#' Calculates an approximate memory requirement (in bytes) based on the
#' data dimensions stored in a `gwpr_context` object, the number of
#' parallel workers, and whether the full distance matrix will be cached.
#'
#' The two main cost components are:
#' \itemize{
#'   \item Distance matrix (when \code{cache_distance = TRUE}):
#'         \eqn{n\_units^2 \times 8} bytes.
#'   \item Local model working copies (one per worker):
#'         \eqn{n\_rows \times n\_vars \times 8 \times workers} bytes.
#' }
#'
#' @param context A `gwpr_context` list containing at least
#'   \code{metadata$n_units}, \code{metadata$n_time}, and
#'   \code{metadata$n_vars}.  If these keys are absent, the function falls
#'   back to direct arguments \code{n_units}, \code{n_time}, and
#'   \code{n_vars} when supplied via \code{...}.  Alternatively, pass a
#'   plain list with the required scalar fields directly.
#' @param workers Positive integer.  Number of parallel workers.
#' @param cache_distance Logical or \code{NULL}.  When \code{TRUE} the full
#'   \eqn{n\_units \times n\_units} distance matrix is assumed to be kept in
#'   memory.  When \code{NULL} (default) the function assumes caching is
#'   enabled for a conservative estimate.
#'
#' @return A named list with class \code{"gwpr_memory_estimate"}:
#'   \describe{
#'     \item{n_units}{Number of spatial units.}
#'     \item{n_time}{Number of time periods.}
#'     \item{n_vars}{Number of explanatory variables.}
#'     \item{n_rows}{Total panel rows (\code{n_units * n_time}).}
#'     \item{workers}{Workers used for the estimate.}
#'     \item{cache_distance}{Whether distance caching was assumed.}
#'     \item{distance_bytes}{Bytes for the distance matrix (0 if not cached).}
#'     \item{model_bytes}{Bytes for local model copies across workers.}
#'     \item{total_bytes}{Total estimated bytes.}
#'     \item{risk}{Character risk level: \code{"low"}, \code{"medium"}, or
#'       \code{"high"}.}
#'   }
#' @keywords internal
estimate_memory <- function(context, workers = 1, cache_distance = NULL) {
  # ---- resolve dimensions ------------------------------------------------
  if (inherits(context, "gwpr_context") || is.list(context)) {
    # Try metadata sub-list first (populated by data_prepare)
    meta <- context[["metadata"]]
    n_units <- meta[["n_units"]] %||% context[["n_units"]]
    n_time  <- meta[["n_time"]]  %||% context[["n_time"]]
    n_vars  <- meta[["n_vars"]]  %||% context[["n_vars"]]

    # Fall back to model_matrix dimensions when available
    if (is.null(n_vars) && !is.null(context[["model_matrix"]])) {
      n_vars <- ncol(context[["model_matrix"]])
    }
    # Fall back to counting unique IDs / times from panel_data
    if (!is.null(context[["panel_data"]])) {
      pd <- context[["panel_data"]]
      id_col   <- context[["id"]]
      time_col <- context[["time"]]
      if (is.null(n_units) && !is.null(id_col) && id_col %in% names(pd)) {
        n_units <- length(unique(pd[[id_col]]))
      }
      if (is.null(n_time) && !is.null(time_col) && time_col %in% names(pd)) {
        n_time <- length(unique(pd[[time_col]]))
      }
    }
  } else {
    stop("`context` must be a list or gwpr_context object.", call. = FALSE)
  }

  # ---- validate resolved values ------------------------------------------
  if (is.null(n_units) || !is.numeric(n_units) || n_units <= 0) {
    stop(
      "Cannot resolve `n_units` from `context`. ",
      "Set `context$metadata$n_units` or `context$n_units`.",
      call. = FALSE
    )
  }
  if (is.null(n_time) || !is.numeric(n_time) || n_time <= 0) {
    stop(
      "Cannot resolve `n_time` from `context`. ",
      "Set `context$metadata$n_time` or `context$n_time`.",
      call. = FALSE
    )
  }
  if (is.null(n_vars) || !is.numeric(n_vars) || n_vars <= 0) {
    stop(
      "Cannot resolve `n_vars` from `context`. ",
      "Set `context$metadata$n_vars` or `context$n_vars`.",
      call. = FALSE
    )
  }

  if (!is.numeric(workers) || length(workers) != 1L || workers < 1) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)

  # Default: assume caching enabled (conservative)
  if (is.null(cache_distance)) cache_distance <- TRUE

  n_units <- as.numeric(n_units)
  n_time  <- as.numeric(n_time)
  n_vars  <- as.numeric(n_vars)
  n_rows  <- n_units * n_time

  bytes_per_double <- 8

  # Distance matrix: n_units x n_units doubles
  distance_bytes <- if (isTRUE(cache_distance)) {
    n_units^2 * bytes_per_double
  } else {
    0
  }

  # Local model working copies: one per worker
  model_bytes <- n_rows * n_vars * bytes_per_double * workers

  total_bytes <- distance_bytes + model_bytes

  result <- list(
    n_units        = n_units,
    n_time         = n_time,
    n_vars         = n_vars,
    n_rows         = n_rows,
    workers        = workers,
    cache_distance = cache_distance,
    distance_bytes = distance_bytes,
    model_bytes    = model_bytes,
    total_bytes    = total_bytes,
    risk           = classify_memory_risk(total_bytes)
  )
  structure(result, class = "gwpr_memory_estimate")
}

# ---------------------------------------------------------------------------
# classify_memory_risk
# ---------------------------------------------------------------------------

#' Classify memory risk level
#'
#' Maps an estimated byte count to a human-readable risk category.
#'
#' Thresholds:
#' \itemize{
#'   \item \strong{low}:    < 500 MB
#'   \item \strong{medium}: 500 MB -- 2 GB
#'   \item \strong{high}:   > 2 GB
#' }
#'
#' @param estimated_bytes Non-negative numeric; total estimated bytes.
#'
#' @return A character string: \code{"low"}, \code{"medium"}, or
#'   \code{"high"}.
#' @keywords internal
classify_memory_risk <- function(estimated_bytes) {
  if (!is.numeric(estimated_bytes) || length(estimated_bytes) != 1L) {
    stop("`estimated_bytes` must be a single numeric value.", call. = FALSE)
  }
  if (estimated_bytes < 0) {
    stop("`estimated_bytes` must be non-negative.", call. = FALSE)
  }

  mb_500 <- 500  * 1024^2   # 500 MB in bytes
  gb_2   <- 2    * 1024^3   #   2 GB in bytes

  if (estimated_bytes < mb_500) {
    "low"
  } else if (estimated_bytes <= gb_2) {
    "medium"
  } else {
    "high"
  }
}

# ---------------------------------------------------------------------------
# format_memory_warning
# ---------------------------------------------------------------------------

#' Format a human-readable memory warning message
#'
#' Converts a `gwpr_memory_estimate` object (produced by
#' \code{\link{estimate_memory}}) into a readable character string.  For
#' high-risk estimates the message also includes actionable suggestions.
#'
#' @param memory_estimate A `gwpr_memory_estimate` list, typically the
#'   return value of \code{\link{estimate_memory}}.
#'
#' @return A character string containing the warning text.  The string is
#'   suitable for passing to \code{message()} or \code{warning()}.
#' @keywords internal
format_memory_warning <- function(memory_estimate) {
  if (!inherits(memory_estimate, "gwpr_memory_estimate")) {
    stop(
      "`memory_estimate` must be a `gwpr_memory_estimate` object ",
      "(returned by `estimate_memory()`).",
      call. = FALSE
    )
  }

  total_gb  <- memory_estimate$total_bytes / 1024^3
  risk      <- memory_estimate$risk
  workers   <- memory_estimate$workers
  cache     <- memory_estimate$cache_distance
  n_units   <- memory_estimate$n_units

  risk_label <- switch(
    risk,
    low    = "LOW",
    medium = "MEDIUM",
    high   = "HIGH"
  )

  header <- sprintf(
    "[GWPR Memory Estimate] Risk: %s | Estimated: %.2f GB",
    risk_label, total_gb
  )

  detail <- sprintf(
    "  n_units=%g, n_time=%g, n_vars=%g, n_rows=%g, workers=%d, cache_distance=%s",
    memory_estimate$n_units,
    memory_estimate$n_time,
    memory_estimate$n_vars,
    memory_estimate$n_rows,
    workers,
    if (isTRUE(cache)) "TRUE" else "FALSE"
  )

  msg <- paste(header, detail, sep = "\n")

  if (risk == "high") {
    suggestions <- paste(
      "  Suggestions to reduce memory usage:",
      "    - Reduce `workers` (currently %d). Try workers = 1.",
      "    - Disable distance caching: set cache_distance = FALSE.",
      "    - Reduce the number of spatial units (currently %g).",
      "    - Subset the data to fewer time periods.",
      sep = "\n"
    )
    suggestions <- sprintf(suggestions, workers, n_units)
    msg <- paste(msg, suggestions, sep = "\n")
  } else if (risk == "medium") {
    msg <- paste(
      msg,
      "  Consider monitoring system memory during the run.",
      sep = "\n"
    )
  }

  msg
}

# ---------------------------------------------------------------------------
# Helper: null-coalescing operator (local, not exported)
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x
