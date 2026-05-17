#' @title Grid Search for Bandwidth Selection
#'
#' @description
#' Functions implementing an exhaustive grid search over a user-specified
#' range of bandwidth candidates.  Each candidate is evaluated by a user-
#' supplied scorer function; the full search history and the best bandwidth
#' are returned as a `gwpr_bandwidth` object.
#'
#' @name bandwidth_grid
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# make_grid_candidates
# ---------------------------------------------------------------------------

#' Generate grid search bandwidth candidates
#'
#' Enumerates candidate bandwidths from `lower` to (but not necessarily
#' reaching) `upper` with increment `step`.  All three boundary arguments
#' are mandatory — they are never inferred automatically.
#'
#' @param lower   Numeric scalar.  Lower boundary (inclusive).
#' @param upper   Numeric scalar.  Upper boundary (inclusive if reachable by
#'   `step`).
#' @param step    Numeric scalar > 0.  Increment between successive candidates.
#' @param adaptive Logical scalar.  When `TRUE`, candidate values are rounded
#'   to the nearest positive integer and duplicates are removed.
#'
#' @return A numeric vector of candidate bandwidth values.
#'
#' @examples
#' make_grid_candidates(lower = 1, upper = 5, step = 1, adaptive = FALSE)
#' make_grid_candidates(lower = 2, upper = 10, step = 2.5, adaptive = TRUE)
#'
#' @export
make_grid_candidates <- function(lower, upper, step, adaptive) {
  if (missing(lower) || is.null(lower)) {
    stop(
      "Grid search requires `lower` to be set explicitly. ",
      "Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (missing(upper) || is.null(upper)) {
    stop(
      "Grid search requires `upper` to be set explicitly. ",
      "Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (missing(step) || is.null(step)) {
    stop(
      "Grid search requires `step` to be set explicitly. ",
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
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0) {
    stop("`step` must be a single finite positive numeric value.", call. = FALSE)
  }
  if (lower >= upper) {
    stop("`lower` must be strictly less than `upper`.", call. = FALSE)
  }

  candidates <- seq(from = lower, to = upper, by = step)

  if (adaptive) {
    candidates <- unique(pmax(1L, round(candidates)))
    candidates <- as.integer(candidates)
  }

  candidates
}

# ---------------------------------------------------------------------------
# score_bandwidth_candidate
# ---------------------------------------------------------------------------

#' Score a single bandwidth candidate
#'
#' Calls the `scorer` function for a given bandwidth and records the result
#' together with timing information, model counts, and metric values.
#'
#' @param context   A `gwpr_context` list.
#' @param bandwidth Numeric scalar; the candidate bandwidth to evaluate.
#' @param scorer    A function with signature `scorer(context, bandwidth)`
#'   returning a named list with at minimum:
#'   \describe{
#'     \item{`score`}{Numeric scalar.  The criterion value (lower is better).}
#'     \item{`criterion`}{Character scalar.  Name of the criterion.}
#'     \item{`n_local_models`}{Integer.  Number of local models attempted.}
#'     \item{`n_failed_local_models`}{Integer.  Number of local models that
#'       failed.}
#'     \item{`metrics`}{Named list of aggregate metrics (e.g. R2, MSE, etc.).}
#'   }
#'
#' @return A named list describing the candidate result:
#' \describe{
#'   \item{`bandwidth`}{The evaluated bandwidth.}
#'   \item{`score`}{Numeric criterion score, or `NA_real_` on failure.}
#'   \item{`criterion`}{Name of the scoring criterion.}
#'   \item{`status`}{`"ok"` or `"failed"`.}
#'   \item{`error_message`}{`NA_character_` or error text.}
#'   \item{`warning_message`}{`NA_character_` or warning text.}
#'   \item{`elapsed_time`}{Elapsed wall-clock time in seconds.}
#'   \item{`n_local_models`}{Number of local models attempted.}
#'   \item{`n_failed_local_models`}{Number of local models that failed.}
#'   \item{`r2`, `mse`, `rmse`, `mae`}{Linear metrics, or `NA_real_`.}
#'   \item{`log_loss`, `accuracy`, `precision`, `recall`, `f1_score`}{Logistic
#'     metrics, or `NA_real_`.}
#' }
#' @keywords internal
score_bandwidth_candidate <- function(context, bandwidth, scorer) {
  t_start <- proc.time()[["elapsed"]]

  warn_msgs <- character()

  raw <- tryCatch(
    withCallingHandlers(
      scorer(context, bandwidth),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(
        status          = "failed",
        error_message   = conditionMessage(e),
        score           = NA_real_,
        criterion       = NA_character_,
        n_local_models        = NA_integer_,
        n_failed_local_models = NA_integer_,
        metrics         = list()
      )
    }
  )

  elapsed <- proc.time()[["elapsed"]] - t_start

  status        <- raw$status        %||% "ok"
  error_message <- raw$error_message %||% NA_character_
  warning_message <- if (length(warn_msgs) > 0L) {
    paste(warn_msgs, collapse = "; ")
  } else {
    NA_character_
  }

  metrics <- raw$metrics %||% list()

  # Linear metrics
  r2   <- metrics$R2   %||% NA_real_
  mse  <- metrics$MSE  %||% NA_real_
  rmse <- metrics$RMSE %||% NA_real_
  mae  <- metrics$MAE  %||% NA_real_

  # Logistic metrics
  log_loss  <- metrics$log_loss  %||% NA_real_
  accuracy  <- metrics$accuracy  %||% NA_real_
  precision <- metrics$precision %||% NA_real_
  recall    <- metrics$recall    %||% NA_real_
  f1_score  <- metrics$f1_score  %||% NA_real_

  list(
    bandwidth             = bandwidth,
    score                 = raw$score             %||% NA_real_,
    criterion             = raw$criterion         %||% NA_character_,
    status                = status,
    error_message         = error_message,
    warning_message       = warning_message,
    elapsed_time          = elapsed,
    n_local_models        = raw$n_local_models        %||% NA_integer_,
    n_failed_local_models = raw$n_failed_local_models %||% NA_integer_,
    r2                    = r2,
    mse                   = mse,
    rmse                  = rmse,
    mae                   = mae,
    log_loss              = log_loss,
    accuracy              = accuracy,
    precision             = precision,
    recall                = recall,
    f1_score              = f1_score
  )
}

# ---------------------------------------------------------------------------
# rank_grid_history
# ---------------------------------------------------------------------------

#' Rank grid search history
#'
#' Assigns integer ranks to candidates: successful candidates are ranked by
#' ascending `score` (rank 1 = best / lowest score); failed candidates are
#' ranked after all successful ones.
#'
#' @param history A data.frame as produced by `search_bandwidth_grid()`.
#'
#' @return The same data.frame with a `rank` column added (or replaced).
#'
#' @examples
#' h <- data.frame(
#'   candidate_id = 1:3,
#'   score  = c(0.5, 0.3, NA),
#'   status = c("ok", "ok", "failed"),
#'   stringsAsFactors = FALSE
#' )
#' rank_grid_history(h)
#'
#' @export
rank_grid_history <- function(history) {
  if (!is.data.frame(history) || nrow(history) == 0L) {
    return(history)
  }

  n <- nrow(history)
  ranks <- integer(n)

  ok_idx   <- which(history$status == "ok")
  fail_idx <- which(history$status != "ok")

  if (length(ok_idx) > 0L) {
    ok_scores <- history$score[ok_idx]
    ok_ranks  <- rank(ok_scores, ties.method = "first")
    ranks[ok_idx] <- ok_ranks
  }

  if (length(fail_idx) > 0L) {
    fail_start <- length(ok_idx) + 1L
    ranks[fail_idx] <- seq(fail_start, length.out = length(fail_idx))
  }

  history$rank <- ranks
  history
}

# ---------------------------------------------------------------------------
# search_bandwidth_grid  (main entry point)
# ---------------------------------------------------------------------------

#' Bandwidth grid search
#'
#' Exhaustively evaluates every bandwidth candidate in the grid defined by
#' `control$lower`, `control$upper`, and `control$step`.  Candidates are
#' scored by the user-supplied `scorer` function.  The complete history and
#' the best bandwidth are returned as a `gwpr_bandwidth` object.
#'
#' @param context  A `gwpr_context` list.  Must contain at minimum `adaptive`.
#' @param control  A list with grid search parameters.  **Required** fields:
#'   \describe{
#'     \item{`lower`}{Numeric.  Lower bound of the search grid (inclusive).}
#'     \item{`upper`}{Numeric.  Upper bound (included when reachable by
#'       `step`).}
#'     \item{`step`}{Numeric > 0.  Increment between candidates.}
#'   }
#' @param scorer   A function with signature `scorer(context, bandwidth)`
#'   that returns a named list as described in `score_bandwidth_candidate()`.
#' @param workers  Integer.  Number of parallel workers.  `1` (default) uses
#'   serial `lapply`; `> 1` uses `parallel_map()`.
#' @param seed     Integer random seed, or `NULL`.  Applied before parallel
#'   execution.
#'
#' @return A `gwpr_bandwidth` object (see `new_gwpr_bandwidth()`).
#'
#' @examples
#' \dontrun{
#' result <- search_bandwidth_grid(
#'   context = ctx,
#'   control = list(lower = 1, upper = 5, step = 1),
#'   scorer  = my_scorer_fn,
#'   workers = 1
#' )
#' result$best_bandwidth
#' }
#'
#' @export
search_bandwidth_grid <- function(context, control, scorer,
                                  workers = 1L, seed = NULL) {

  t_total_start <- proc.time()[["elapsed"]]

  adaptive <- context$adaptive %||% FALSE

  # Generate candidates — will stop with informative errors if bounds missing
  candidates <- make_grid_candidates(
    lower    = control$lower,
    upper    = control$upper,
    step     = control$step,
    adaptive = adaptive
  )

  n_candidates <- length(candidates)

  # Set seed before evaluation
  if (!is.null(seed)) set.seed(seed)

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

  # Determine bandwidth_type and method/criterion from first successful result
  bandwidth_type <- if (adaptive) "adaptive" else "fixed"
  first_ok <- Find(function(r) r$status == "ok", results_list)
  criterion <- if (!is.null(first_ok)) first_ok$criterion %||% NA_character_ else NA_character_

  # Assemble history data.frame
  history <- data.frame(
    candidate_id          = seq_len(n_candidates),
    bandwidth             = vapply(results_list, `[[`, numeric(1L), "bandwidth"),
    bandwidth_type        = bandwidth_type,
    method                = "grid",
    criterion             = vapply(results_list, function(r) r$criterion %||% NA_character_, character(1L)),
    score                 = vapply(results_list, function(r) r$score %||% NA_real_, numeric(1L)),
    rank                  = NA_integer_,
    status                = vapply(results_list, `[[`, character(1L), "status"),
    error_message         = vapply(results_list, function(r) r$error_message   %||% NA_character_, character(1L)),
    warning_message       = vapply(results_list, function(r) r$warning_message %||% NA_character_, character(1L)),
    elapsed_time          = vapply(results_list, `[[`, numeric(1L), "elapsed_time"),
    n_local_models        = vapply(results_list, function(r) {
      v <- r$n_local_models; if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
    }, integer(1L)),
    n_failed_local_models = vapply(results_list, function(r) {
      v <- r$n_failed_local_models; if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
    }, integer(1L)),
    r2                    = vapply(results_list, function(r) r$r2        %||% NA_real_, numeric(1L)),
    mse                   = vapply(results_list, function(r) r$mse       %||% NA_real_, numeric(1L)),
    rmse                  = vapply(results_list, function(r) r$rmse      %||% NA_real_, numeric(1L)),
    mae                   = vapply(results_list, function(r) r$mae       %||% NA_real_, numeric(1L)),
    log_loss              = vapply(results_list, function(r) r$log_loss  %||% NA_real_, numeric(1L)),
    accuracy              = vapply(results_list, function(r) r$accuracy  %||% NA_real_, numeric(1L)),
    precision             = vapply(results_list, function(r) r$precision %||% NA_real_, numeric(1L)),
    recall                = vapply(results_list, function(r) r$recall    %||% NA_real_, numeric(1L)),
    f1_score              = vapply(results_list, function(r) r$f1_score  %||% NA_real_, numeric(1L)),
    stringsAsFactors      = FALSE
  )

  # Rank candidates
  history <- rank_grid_history(history)

  # Identify best candidate (rank == 1 among ok)
  ok_rows <- history[history$status == "ok", , drop = FALSE]
  if (nrow(ok_rows) > 0L) {
    best_row      <- ok_rows[which.min(ok_rows$rank), , drop = FALSE]
    best_bandwidth <- best_row$bandwidth
    best_score     <- best_row$score
  } else {
    best_bandwidth <- NA_real_
    best_score     <- NA_real_
  }

  total_elapsed <- proc.time()[["elapsed"]] - t_total_start

  new_gwpr_bandwidth(
    method           = "grid",
    best_bandwidth   = best_bandwidth,
    best_score       = best_score,
    criterion        = criterion,
    history          = history,
    metrics_history  = NULL,
    seed             = seed,
    convergence_info = list(n_candidates = n_candidates),
    elapsed_time     = total_elapsed,
    warnings         = character()
  )
}

# ---------------------------------------------------------------------------
# Internal helper (local copy to avoid cross-file ordering issues)
# ---------------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# parallel_map stub (used when parallel.R is not yet loaded)
# ---------------------------------------------------------------------------

#' Minimal parallel_map implementation
#'
#' Wraps `lapply` when `workers = 1` and `parallel::mclapply` when
#' `workers > 1`.  This stub is superseded once `parallel.R` is available.
#'
#' @param x       List of inputs.
#' @param fn      Function to apply to each element.
#' @param workers Integer number of workers.
#' @param seed    Random seed (passed to `mc.set.seed`).
#' @param ...     Additional arguments passed to `fn`.
#' @keywords internal
parallel_map <- function(x, fn, workers = 1L, seed = NULL, ...) {
  if (workers <= 1L) {
    lapply(x, fn, ...)
  } else {
    if (!is.null(seed)) set.seed(seed)
    parallel::mclapply(x, fn, mc.cores = workers, mc.set.seed = TRUE, ...)
  }
}
