#' @title Random Bandwidth Optimizer (`bandwidth_random.R`)
#'
#' @description
#' Implements a bounded random search for bandwidth selection.  A user-
#' specified number of candidate bandwidths are drawn uniformly at random from
#' `[lower, upper]`, scored by a user-supplied scorer function, and the
#' candidate with the lowest score is returned as the best bandwidth.
#'
#' The search boundaries (`lower`, `upper`) **must** be set explicitly by the
#' caller; automatic inference is intentionally not supported.
#'
#' @name bandwidth_random
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# sample_random_bandwidths
# ---------------------------------------------------------------------------

#' Sample random bandwidth candidates
#'
#' Draws `n_samples` values uniformly from `[lower, upper]`.  When
#' `adaptive = TRUE` the draws are rounded to the nearest positive integer and
#' duplicates are removed (the resulting vector may therefore be shorter than
#' `n_samples`; this is intentional and no error is raised).
#'
#' Both `lower` and `upper` **must** be supplied explicitly — they are never
#' inferred automatically.
#'
#' @param lower     Numeric scalar.  Lower bound of the search range
#'   (inclusive).
#' @param upper     Numeric scalar.  Upper bound of the search range
#'   (inclusive).
#' @param n_samples Positive integer.  Number of random draws before any
#'   deduplication.
#' @param adaptive  Logical scalar.  When `TRUE`, draws are rounded to the
#'   nearest positive integer and duplicates are removed.
#' @param seed      Integer random seed, or `NULL`.  When supplied, `set.seed()`
#'   is called before sampling so that the candidate set is reproducible.
#'
#' @return A numeric vector of candidate bandwidth values.  The vector has
#'   length `n_samples` when `adaptive = FALSE`, and length at most
#'   `n_samples` when `adaptive = TRUE`.
#'
#' @examples
#' sample_random_bandwidths(lower = 1, upper = 10, n_samples = 5,
#'                          adaptive = FALSE, seed = 42L)
#' sample_random_bandwidths(lower = 2, upper = 10, n_samples = 20,
#'                          adaptive = TRUE, seed = 1L)
#'
#' @noRd
sample_random_bandwidths <- function(lower, upper, n_samples,
                                     adaptive = FALSE, seed = NULL) {

  if (missing(lower) || is.null(lower)) {
    stop(
      "Random search requires `lower` to be set explicitly. ",
      "Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (missing(upper) || is.null(upper)) {
    stop(
      "Random search requires `upper` to be set explicitly. ",
      "Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower)) {
    stop("`lower` must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("`upper` must be a single finite numeric value.", call. = FALSE)
  }
  if (lower >= upper) {
    stop("`lower` must be strictly less than `upper`.", call. = FALSE)
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples < 1L) {
    stop("`n_samples` must be a single positive integer.", call. = FALSE)
  }
  n_samples <- as.integer(n_samples)

  if (!is.null(seed)) set.seed(seed)

  candidates <- stats::runif(n_samples, min = lower, max = upper)

  if (adaptive) {
    candidates <- unique(pmax(1L, round(candidates)))
    candidates <- as.integer(candidates)
  }

  candidates
}

# ---------------------------------------------------------------------------
# search_bandwidth_random  (main entry point)
# ---------------------------------------------------------------------------

#' Random bandwidth search
#'
#' Samples `control$n_samples` candidate bandwidths uniformly at random from
#' `[control$lower, control$upper]`, scores each with the user-supplied
#' `scorer` function, and returns the candidate with the lowest score.  The
#' complete evaluation history is included in the returned `gwpr_bandwidth`
#' object.
#'
#' The `method` field is fixed to `"random"`.  The history structure is
#' identical to that produced by `search_bandwidth_grid()`.
#'
#' @param context  A `gwpr_context` list.  Must contain at minimum `adaptive`.
#' @param control  A list with random search parameters.  **Required** fields:
#'   \describe{
#'     \item{`lower`}{Numeric.  Lower bound of the search range (inclusive).}
#'     \item{`upper`}{Numeric.  Upper bound of the search range (inclusive).}
#'   }
#'   Optional fields:
#'   \describe{
#'     \item{`n_samples`}{Positive integer.  Number of random candidates to
#'       draw.  Defaults to `50L`.}
#'   }
#' @param scorer   A function with signature `scorer(context, bandwidth)`
#'   that returns a named list as described in `score_bandwidth_candidate()`.
#' @param workers  Integer.  Number of parallel workers.  `1` (default) uses
#'   serial `lapply`; `> 1` uses `parallel_map()`.
#' @param seed     Integer random seed, or `NULL`.  Applied before sampling and
#'   before parallel execution.  Stored in the returned object.
#'
#' @return A `gwpr_bandwidth` object (see `new_gwpr_bandwidth()`).
#'
#' @examples
#' \donttest{
#' result <- search_bandwidth_random(
#'   context = ctx,
#'   control = list(lower = 1, upper = 10, n_samples = 20L),
#'   scorer  = my_scorer_fn,
#'   workers = 1L,
#'   seed    = 42L
#' )
#' result$best_bandwidth
#' }
#'
#' @noRd
search_bandwidth_random <- function(context, control, scorer,
                                    workers = 1L, seed = NULL) {

  t_total_start <- proc.time()[["elapsed"]]

  adaptive  <- context$adaptive %||% FALSE
  n_samples <- control$n_samples %||% 50L

  # Sample candidates — will stop with informative errors if bounds missing
  candidates <- sample_random_bandwidths(
    lower    = control$lower,
    upper    = control$upper,
    n_samples = n_samples,
    adaptive  = adaptive,
    seed      = seed
  )

  n_candidates <- length(candidates)

  # Evaluate candidates
  if (workers <= 1L) {
    results_list <- lapply(
      candidates,
      function(bw) score_bandwidth_candidate(context, bw, scorer)
    )
  } else {
    results_list <- parallel_map(
      x       = as.list(candidates),
      fn      = function(bw) score_bandwidth_candidate(context, bw, scorer),
      workers = workers,
      seed    = seed
    )
  }

  # Determine bandwidth_type and criterion from first successful result
  bandwidth_type <- if (adaptive) "adaptive" else "fixed"
  first_ok <- Find(function(r) r$status == "ok", results_list)
  criterion <- if (!is.null(first_ok)) {
    first_ok$criterion %||% NA_character_
  } else {
    NA_character_
  }

  # Assemble history data.frame (same structure as grid search)
  history <- data.frame(
    candidate_id          = seq_len(n_candidates),
    bandwidth             = vapply(results_list, `[[`, numeric(1L), "bandwidth"),
    bandwidth_type        = bandwidth_type,
    method                = "random",
    criterion             = vapply(results_list, function(r) {
      r$criterion %||% NA_character_
    }, character(1L)),
    score                 = vapply(results_list, function(r) {
      r$score %||% NA_real_
    }, numeric(1L)),
    rank                  = NA_integer_,
    status                = vapply(results_list, `[[`, character(1L), "status"),
    error_message         = vapply(results_list, function(r) {
      r$error_message %||% NA_character_
    }, character(1L)),
    warning_message       = vapply(results_list, function(r) {
      r$warning_message %||% NA_character_
    }, character(1L)),
    elapsed_time          = vapply(results_list, `[[`, numeric(1L), "elapsed_time"),
    n_local_models        = vapply(results_list, function(r) {
      v <- r$n_local_models
      if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
    }, integer(1L)),
    n_failed_local_models = vapply(results_list, function(r) {
      v <- r$n_failed_local_models
      if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
    }, integer(1L)),
    r2                    = vapply(results_list, function(r) {
      r$r2        %||% NA_real_
    }, numeric(1L)),
    mse                   = vapply(results_list, function(r) {
      r$mse       %||% NA_real_
    }, numeric(1L)),
    rmse                  = vapply(results_list, function(r) {
      r$rmse      %||% NA_real_
    }, numeric(1L)),
    mae                   = vapply(results_list, function(r) {
      r$mae       %||% NA_real_
    }, numeric(1L)),
    log_loss              = vapply(results_list, function(r) {
      r$log_loss  %||% NA_real_
    }, numeric(1L)),
    accuracy              = vapply(results_list, function(r) {
      r$accuracy  %||% NA_real_
    }, numeric(1L)),
    precision             = vapply(results_list, function(r) {
      r$precision %||% NA_real_
    }, numeric(1L)),
    recall                = vapply(results_list, function(r) {
      r$recall    %||% NA_real_
    }, numeric(1L)),
    f1_score              = vapply(results_list, function(r) {
      r$f1_score  %||% NA_real_
    }, numeric(1L)),
    stringsAsFactors      = FALSE
  )

  # Rank candidates (reuse grid ranking logic)
  history <- rank_grid_history(history)

  # Identify best candidate (rank == 1 among ok)
  ok_rows <- history[history$status == "ok", , drop = FALSE]
  if (nrow(ok_rows) > 0L) {
    best_row       <- ok_rows[which.min(ok_rows$rank), , drop = FALSE]
    best_bandwidth <- best_row$bandwidth
    best_score     <- best_row$score
  } else {
    best_bandwidth <- NA_real_
    best_score     <- NA_real_
  }

  total_elapsed <- proc.time()[["elapsed"]] - t_total_start

  new_gwpr_bandwidth(
    method           = "random",
    best_bandwidth   = best_bandwidth,
    best_score       = best_score,
    criterion        = criterion,
    history          = history,
    metrics_history  = NULL,
    seed             = seed,
    convergence_info = list(n_candidates = n_candidates,
                            n_samples    = n_samples),
    elapsed_time     = total_elapsed,
    warnings         = character()
  )
}

# ---------------------------------------------------------------------------
# Internal helpers (local copies to avoid cross-file ordering issues)
# ---------------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y
