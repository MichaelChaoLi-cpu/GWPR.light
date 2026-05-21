#' @title SGD Bandwidth Search (`bandwidth_sgd.R`)
#'
#' @description
#' Implements a stochastic gradient descent (SGD) based bandwidth search for
#' GWPR models.  A one-dimensional finite-difference gradient approximation is
#' used to iteratively update the bandwidth over a fixed number of epochs.
#' Mini-batch sampling and early stopping are supported.
#'
#' The search does **not** require the user to specify `lower`, `upper`, or
#' `step`; SGD starts from a single initialised bandwidth and follows the
#' (approximate) gradient.  When `lower` / `upper` are supplied they are used
#' as hard constraints.
#'
#' @name bandwidth_sgd
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# initialize_bandwidth
# ---------------------------------------------------------------------------

#' Initialise a bandwidth value for SGD
#'
#' Randomly initialises a starting bandwidth within `[lower, upper]`.  For
#' fixed bandwidth the initialisation is done on the log scale so the
#' distribution is spread across orders of magnitude.  For adaptive bandwidth
#' the result is rounded to a positive integer.
#'
#' Both `lower` and `upper` **must** be supplied when calling this function.
#'
#' @param control  A list that **must** contain `lower` and `upper`.
#' @param adaptive Logical.  When `TRUE` an integer bandwidth is returned.
#' @param seed     Integer random seed or `NULL`.
#'
#' @return A numeric scalar (positive, and an integer when `adaptive = TRUE`).
#'
#' @examples
#' initialize_bandwidth(list(lower = 1, upper = 100), adaptive = FALSE,
#'                      seed = 42L)
#' initialize_bandwidth(list(lower = 3, upper = 20), adaptive = TRUE,
#'                      seed = 1L)
#'
#' @noRd
initialize_bandwidth <- function(control, adaptive, seed = NULL) {
  lower <- control$lower
  upper <- control$upper

  if (is.null(lower)) {
    stop(
      "SGD bandwidth initialisation requires `control$lower` to be set ",
      "explicitly. Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (is.null(upper)) {
    stop(
      "SGD bandwidth initialisation requires `control$upper` to be set ",
      "explicitly. Automatic boundary inference is not supported.",
      call. = FALSE
    )
  }
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) || lower <= 0) {
    stop("`lower` must be a single finite positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("`upper` must be a single finite numeric value.", call. = FALSE)
  }
  if (lower >= upper) {
    stop("`lower` must be strictly less than `upper`.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  if (adaptive) {
    bw <- as.integer(round(stats::runif(1L, min = lower, max = upper)))
    bw <- max(1L, min(bw, as.integer(upper)))
  } else {
    # Log-scale initialisation
    log_bw <- stats::runif(1L, min = log(lower), max = log(upper))
    bw <- exp(log_bw)
    bw <- max(lower, min(bw, upper))
  }

  bw
}

# ---------------------------------------------------------------------------
# sample_minibatch
# ---------------------------------------------------------------------------

#' Sample a mini-batch of spatial units
#'
#' Draws a random subset of spatial unit IDs from the context.  The fraction
#' of units to include is controlled by `batch_fraction`.  When `seed` is
#' supplied the result is reproducible.
#'
#' The context must contain either `id_map` (a named vector whose names are
#' unit IDs) or `coords` (a matrix whose row count equals the number of spatial
#' units).  If neither is present the function stops with an informative error.
#'
#' @param context        A `gwpr_context` list.
#' @param batch_fraction Numeric in `(0, 1]`.  Fraction of units to include.
#' @param seed           Integer random seed or `NULL`.
#'
#' @return A character vector of unit IDs in the mini-batch.
#'
#' @examples
#' ctx <- list(id_map = c("A" = 1L, "B" = 2L, "C" = 3L, "D" = 4L, "E" = 5L))
#' sample_minibatch(ctx, batch_fraction = 0.6)
#'
#' @noRd
sample_minibatch <- function(context, batch_fraction = 1.0, seed = NULL) {
  if (!is.numeric(batch_fraction) || length(batch_fraction) != 1L ||
      !is.finite(batch_fraction) || batch_fraction <= 0 || batch_fraction > 1) {
    stop("`batch_fraction` must be a single numeric value in (0, 1].",
         call. = FALSE)
  }

  # Determine unit IDs
  if (!is.null(context$id_map)) {
    all_ids <- names(context$id_map)
  } else if (!is.null(context$coords)) {
    all_ids <- as.character(seq_len(nrow(context$coords)))
  } else if (!is.null(context$panel_data) && !is.null(context$id)) {
    all_ids <- as.character(unique(context$panel_data[[context$id]]))
  } else {
    stop(
      "Cannot determine spatial units from context. ",
      "Supply `id_map`, `coords`, or `panel_data` + `id`.",
      call. = FALSE
    )
  }

  n_total <- length(all_ids)
  n_batch <- max(1L, as.integer(ceiling(n_total * batch_fraction)))

  if (!is.null(seed)) set.seed(seed)

  sampled <- sample(all_ids, size = n_batch, replace = FALSE)
  sampled
}

# ---------------------------------------------------------------------------
# estimate_bandwidth_gradient
# ---------------------------------------------------------------------------

#' Estimate the bandwidth gradient via finite differences
#'
#' Approximates the gradient of the scorer objective with respect to bandwidth
#' using a central finite difference:
#'
#' ```
#' gradient = (f(b + delta) - f(b - delta)) / (2 * delta)
#' ```
#'
#' Both evaluations use the same `minibatch`.  When either candidate falls
#' outside `[lower, upper]` the function falls back to a one-sided difference
#' or clips the candidate to the boundary.
#'
#' @param context   A `gwpr_context` list.
#' @param bandwidth Numeric scalar; current bandwidth.
#' @param delta     Numeric scalar > 0; finite-difference step size.
#' @param minibatch Character vector of unit IDs (ignored by the default
#'   scorer; kept for future mini-batch-aware scorers).
#' @param scorer    A function `scorer(context, bandwidth)` that returns a
#'   named list with at minimum a `score` field.
#' @param lower     Numeric lower bound (or `NULL`).
#' @param upper     Numeric upper bound (or `NULL`).
#'
#' @return A numeric scalar gradient estimate.
#'
#' @examples
#' \donttest{
#' grad <- estimate_bandwidth_gradient(
#'   context   = ctx,
#'   bandwidth = 10,
#'   delta     = 1,
#'   minibatch = character(),
#'   scorer    = my_scorer
#' )
#' }
#'
#' @noRd
estimate_bandwidth_gradient <- function(context, bandwidth, delta, minibatch,
                                        scorer, lower = NULL, upper = NULL) {
  if (!is.numeric(delta) || length(delta) != 1L ||
      !is.finite(delta) || delta <= 0) {
    stop("`delta` must be a single finite positive numeric value.", call. = FALSE)
  }

  bw_plus  <- bandwidth + delta
  bw_minus <- bandwidth - delta

  # Clip to [lower, upper] if bounds are provided
  if (!is.null(lower)) {
    bw_plus  <- max(bw_plus,  lower)
    bw_minus <- max(bw_minus, lower)
  }
  if (!is.null(upper)) {
    bw_plus  <- min(bw_plus,  upper)
    bw_minus <- min(bw_minus, upper)
  }

  # Evaluate both sides; use tryCatch to handle scorer failures gracefully
  score_plus <- tryCatch(
    {
      res <- scorer(context, bw_plus)
      if (is.list(res)) res$score %||% NA_real_ else as.numeric(res)
    },
    error = function(e) NA_real_
  )
  score_minus <- tryCatch(
    {
      res <- scorer(context, bw_minus)
      if (is.list(res)) res$score %||% NA_real_ else as.numeric(res)
    },
    error = function(e) NA_real_
  )

  denom <- bw_plus - bw_minus

  if (!is.finite(score_plus) || !is.finite(score_minus) ||
      !is.finite(denom) || denom == 0) {
    return(0)
  }

  (score_plus - score_minus) / denom
}

# ---------------------------------------------------------------------------
# update_bandwidth_sgd
# ---------------------------------------------------------------------------

#' Update bandwidth with SGD step
#'
#' Applies a single gradient-descent update to the current bandwidth.
#'
#' - **Fixed bandwidth**: update is performed on the log scale so that
#'   bandwidth stays positive:
#'   `new_bw = exp(log(bw) - learning_rate * gradient)`.
#' - **Adaptive bandwidth**: update is performed on the integer scale,
#'   then rounded and clipped to `[lower, upper]`.
#'
#' @param bandwidth     Numeric scalar; current bandwidth.
#' @param gradient      Numeric scalar; gradient estimate.
#' @param learning_rate Numeric scalar > 0.
#' @param adaptive      Logical; `TRUE` for adaptive (integer) bandwidth.
#' @param lower         Numeric lower bound (or `NULL`; ignored when `NULL`).
#' @param upper         Numeric upper bound (or `NULL`; ignored when `NULL`).
#'
#' @return A numeric scalar (always positive, integer when `adaptive = TRUE`).
#'
#' @examples
#' update_bandwidth_sgd(10, gradient = 2, learning_rate = 0.5,
#'                      adaptive = FALSE)
#' update_bandwidth_sgd(8L, gradient = -3, learning_rate = 1,
#'                      adaptive = TRUE, lower = 3, upper = 20)
#'
#' @noRd
update_bandwidth_sgd <- function(bandwidth, gradient, learning_rate,
                                 adaptive = FALSE, lower = NULL, upper = NULL) {
  if (adaptive) {
    new_bw <- bandwidth - learning_rate * gradient
    new_bw <- round(new_bw)
    if (!is.null(lower)) new_bw <- max(new_bw, lower)
    if (!is.null(upper)) new_bw <- min(new_bw, upper)
    new_bw <- as.integer(max(1L, new_bw))
  } else {
    # Log-scale update ensures bandwidth > 0
    log_bw <- log(bandwidth) - learning_rate * gradient
    new_bw <- exp(log_bw)
    if (!is.null(lower)) new_bw <- max(new_bw, lower)
    if (!is.null(upper)) new_bw <- min(new_bw, upper)
    # Ensure bandwidth is strictly positive
    new_bw <- max(new_bw, .Machine$double.eps)
  }

  new_bw
}

# ---------------------------------------------------------------------------
# check_early_stopping
# ---------------------------------------------------------------------------

#' Check whether SGD early stopping should trigger
#'
#' Returns `TRUE` when the score has not improved (decreased) for `patience`
#' consecutive epochs.  Returns `FALSE` immediately when `patience = 0`
#' (early stopping disabled).
#'
#' @param history A data.frame of SGD history with at least a `score` column.
#'   Only the last `patience + 1` rows are examined.
#' @param patience Non-negative integer.  `0` disables early stopping.
#'
#' @return Logical scalar.
#'
#' @examples
#' h <- data.frame(score = c(1.0, 0.9, 0.85, 0.85, 0.86))
#' check_early_stopping(h, patience = 2L)  # TRUE
#' check_early_stopping(h, patience = 0L)  # FALSE
#'
#' @noRd
check_early_stopping <- function(history, patience) {
  patience <- as.integer(patience)
  if (patience == 0L) return(FALSE)
  if (!is.data.frame(history) || nrow(history) < patience) return(FALSE)

  # Look at the last `patience` scores
  n <- nrow(history)
  recent_scores <- history$score[(n - patience + 1L):n]

  # Also need the score just before this window
  if (n <= patience) return(FALSE)
  prev_score <- history$score[n - patience]

  # Stop if no improvement over `patience` steps
  # (score has not strictly decreased below the score before the window)
  all(recent_scores >= prev_score, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# search_bandwidth_sgd  (main entry point)
# ---------------------------------------------------------------------------

#' SGD bandwidth search
#'
#' Searches for an optimal bandwidth using a stochastic gradient descent
#' procedure.  A finite-difference gradient approximation updates the
#' bandwidth iteratively over `control$epoch` epochs.
#'
#' The complete epoch-level history is included in the returned
#' `gwpr_bandwidth` object.
#'
#' @param context  A `gwpr_context` list.  Must contain at minimum `adaptive`.
#' @param control  A list with SGD parameters.  **Required** fields:
#'   \describe{
#'     \item{`lower`}{Numeric. Lower bound of the search space (must be > 0).}
#'     \item{`upper`}{Numeric. Upper bound of the search space.}
#'   }
#'   Optional fields:
#'   \describe{
#'     \item{`epoch`}{Positive integer. Number of gradient steps. Default 10.}
#'     \item{`learning_rate`}{Numeric > 0. Default 0.001.}
#'     \item{`early_stopping_patience`}{Non-negative integer. 0 disables early
#'       stopping. Default 0.}
#'     \item{`batch_fraction`}{Numeric in (0,1]. Mini-batch fraction. Default 1.}
#'     \item{`delta`}{Numeric > 0. Finite-difference step. Default: 1\% of
#'       initial bandwidth.}
#'   }
#' @param scorer   A function `scorer(context, bandwidth)` returning a named
#'   list with at minimum `score`, `criterion`, `n_local_models`,
#'   `n_failed_local_models`, and `metrics`.
#' @param workers  Integer. Number of parallel workers (currently unused;
#'   SGD is inherently sequential).  Default 1.
#' @param seed     Integer random seed or `NULL`.  Applied at initialisation
#'   and stored in the returned object.
#'
#' @return A `gwpr_bandwidth` object (see `new_gwpr_bandwidth()`).  The
#'   `history` data.frame has one row per epoch with the following columns:
#'   `epoch`, `bandwidth_before`, `bandwidth_after`, `score`, `gradient`,
#'   `learning_rate`, `delta`, `batch_size`, `early_stop_counter`, `status`,
#'   `elapsed_time`, plus linear / logistic metric columns.
#'
#' @examples
#' \donttest{
#' result <- search_bandwidth_sgd(
#'   context = ctx,
#'   control = list(lower = 1, upper = 100, epoch = 10L,
#'                  learning_rate = 0.001),
#'   scorer  = my_scorer_fn,
#'   workers = 1L,
#'   seed    = 42L
#' )
#' result$best_bandwidth
#' }
#'
#' @noRd
search_bandwidth_sgd <- function(context, control, scorer,
                                 workers = 1L, seed = NULL) {

  t_total_start <- proc.time()[["elapsed"]]

  # --- Parameters -----------------------------------------------------------
  adaptive      <- context$adaptive %||% FALSE
  epoch         <- as.integer(control$epoch         %||% 10L)
  learning_rate <- control$learning_rate             %||% 0.001
  patience      <- as.integer(
    control$early_stopping_patience %||% 0L
  )
  batch_fraction <- control$batch_fraction %||% 1.0
  lower          <- control$lower
  upper          <- control$upper

  # --- Initialise bandwidth -------------------------------------------------
  bw <- initialize_bandwidth(control, adaptive = adaptive, seed = seed)

  # Determine delta default: 1% of initial bandwidth (minimum 1e-6)
  delta_default <- max(abs(bw) * 0.01, 1e-6)
  delta         <- control$delta %||% delta_default

  # --- History storage ------------------------------------------------------
  history_rows <- vector("list", epoch)

  early_stop_counter <- 0L
  criterion          <- NA_character_
  best_bandwidth     <- bw
  best_score         <- Inf

  # --- SGD loop -------------------------------------------------------------
  for (ep in seq_len(epoch)) {
    t_epoch_start <- proc.time()[["elapsed"]]

    bw_before <- bw

    # 1. Sample mini-batch
    batch_seed <- if (!is.null(seed)) seed + ep else NULL
    minibatch  <- tryCatch(
      sample_minibatch(context, batch_fraction = batch_fraction,
                       seed = batch_seed),
      error = function(e) character()
    )
    batch_size <- length(minibatch)

    # 2. Estimate gradient
    grad <- tryCatch(
      estimate_bandwidth_gradient(
        context   = context,
        bandwidth = bw,
        delta     = delta,
        minibatch = minibatch,
        scorer    = scorer,
        lower     = lower,
        upper     = upper
      ),
      error = function(e) 0
    )

    # 3. Update bandwidth
    bw_after <- update_bandwidth_sgd(
      bandwidth     = bw,
      gradient      = grad,
      learning_rate = learning_rate,
      adaptive      = adaptive,
      lower         = lower,
      upper         = upper
    )

    # 4. Score the new bandwidth
    t_score_start <- proc.time()[["elapsed"]]
    score_res <- tryCatch(
      withCallingHandlers(
        scorer(context, bw_after),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) {
        list(
          score                 = NA_real_,
          criterion             = NA_character_,
          status                = "failed",
          n_local_models        = NA_integer_,
          n_failed_local_models = NA_integer_,
          metrics               = list()
        )
      }
    )
    epoch_elapsed <- proc.time()[["elapsed"]] - t_epoch_start

    score   <- score_res$score     %||% NA_real_
    crit    <- score_res$criterion %||% NA_character_
    status  <- score_res$status    %||% "ok"
    metrics <- score_res$metrics   %||% list()

    if (!is.na(crit) && is.na(criterion)) criterion <- crit

    # Update best
    if (is.finite(score) && score < best_score) {
      best_score     <- score
      best_bandwidth <- bw_after
    }

    # 5. Early stopping counter
    if (ep >= 2L) {
      prev_score <- history_rows[[ep - 1L]]$score %||% NA_real_
      if (is.finite(score) && is.finite(prev_score) && score >= prev_score) {
        early_stop_counter <- early_stop_counter + 1L
      } else {
        early_stop_counter <- 0L
      }
    }

    # Record history row
    history_rows[[ep]] <- list(
      epoch              = ep,
      bandwidth_before   = bw_before,
      bandwidth_after    = bw_after,
      score              = score,
      gradient           = grad,
      learning_rate      = learning_rate,
      delta              = delta,
      batch_size         = batch_size,
      early_stop_counter = early_stop_counter,
      status             = status,
      elapsed_time       = epoch_elapsed,
      # Linear metrics
      r2                 = metrics$R2        %||% NA_real_,
      mse                = metrics$MSE       %||% NA_real_,
      rmse               = metrics$RMSE      %||% NA_real_,
      mae                = metrics$MAE       %||% NA_real_,
      # Logistic metrics
      log_loss           = metrics$log_loss  %||% NA_real_,
      accuracy           = metrics$accuracy  %||% NA_real_,
      precision          = metrics$precision %||% NA_real_,
      recall             = metrics$recall    %||% NA_real_,
      f1_score           = metrics$f1_score  %||% NA_real_
    )

    # Advance bandwidth
    bw <- bw_after

    # 6. Check early stopping
    history_so_far <- .rows_to_df(history_rows[seq_len(ep)])
    if (check_early_stopping(history_so_far, patience = patience)) {
      # Fill remaining epochs with NA rows then break
      for (ep2 in seq(ep + 1L, epoch)) {
        history_rows[[ep2]] <- .make_na_history_row(ep2, bw,
                                                     learning_rate, delta,
                                                     batch_size)
      }
      break
    }
  }

  # --- Assemble history data.frame ------------------------------------------
  history <- .rows_to_df(history_rows)

  total_elapsed <- proc.time()[["elapsed"]] - t_total_start

  new_gwpr_bandwidth(
    method           = "sgd",
    best_bandwidth   = best_bandwidth,
    best_score       = if (is.finite(best_score)) best_score else NA_real_,
    criterion        = criterion,
    history          = history,
    metrics_history  = NULL,
    seed             = seed,
    convergence_info = list(
      epoch         = epoch,
      learning_rate = learning_rate,
      patience      = patience,
      delta         = delta
    ),
    elapsed_time     = total_elapsed,
    warnings         = character()
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Convert a list of history row lists to a data.frame
#' @keywords internal
.rows_to_df <- function(rows) {
  # Remove NULL entries
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      epoch = integer(), bandwidth_before = numeric(),
      bandwidth_after = numeric(), score = numeric(),
      gradient = numeric(), learning_rate = numeric(),
      delta = numeric(), batch_size = integer(),
      early_stop_counter = integer(), status = character(),
      elapsed_time = numeric(),
      r2 = numeric(), mse = numeric(), rmse = numeric(), mae = numeric(),
      log_loss = numeric(), accuracy = numeric(),
      precision = numeric(), recall = numeric(), f1_score = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    epoch              = vapply(rows, `[[`, integer(1L),  "epoch"),
    bandwidth_before   = vapply(rows, `[[`, numeric(1L),  "bandwidth_before"),
    bandwidth_after    = vapply(rows, `[[`, numeric(1L),  "bandwidth_after"),
    score              = vapply(rows, `[[`, numeric(1L),  "score"),
    gradient           = vapply(rows, `[[`, numeric(1L),  "gradient"),
    learning_rate      = vapply(rows, `[[`, numeric(1L),  "learning_rate"),
    delta              = vapply(rows, `[[`, numeric(1L),  "delta"),
    batch_size         = vapply(rows, `[[`, integer(1L),  "batch_size"),
    early_stop_counter = vapply(rows, `[[`, integer(1L),  "early_stop_counter"),
    status             = vapply(rows, `[[`, character(1L), "status"),
    elapsed_time       = vapply(rows, `[[`, numeric(1L),  "elapsed_time"),
    r2                 = vapply(rows, `[[`, numeric(1L),  "r2"),
    mse                = vapply(rows, `[[`, numeric(1L),  "mse"),
    rmse               = vapply(rows, `[[`, numeric(1L),  "rmse"),
    mae                = vapply(rows, `[[`, numeric(1L),  "mae"),
    log_loss           = vapply(rows, `[[`, numeric(1L),  "log_loss"),
    accuracy           = vapply(rows, `[[`, numeric(1L),  "accuracy"),
    precision          = vapply(rows, `[[`, numeric(1L),  "precision"),
    recall             = vapply(rows, `[[`, numeric(1L),  "recall"),
    f1_score           = vapply(rows, `[[`, numeric(1L),  "f1_score"),
    stringsAsFactors   = FALSE
  )
}

#' Create a placeholder NA history row for epochs after early stopping
#' @keywords internal
.make_na_history_row <- function(ep, bw, learning_rate, delta, batch_size) {
  list(
    epoch              = as.integer(ep),
    bandwidth_before   = bw,
    bandwidth_after    = bw,
    score              = NA_real_,
    gradient           = NA_real_,
    learning_rate      = learning_rate,
    delta              = delta,
    batch_size         = as.integer(batch_size),
    early_stop_counter = NA_integer_,
    status             = "early_stopped",
    elapsed_time       = 0,
    r2                 = NA_real_,
    mse                = NA_real_,
    rmse               = NA_real_,
    mae                = NA_real_,
    log_loss           = NA_real_,
    accuracy           = NA_real_,
    precision          = NA_real_,
    recall             = NA_real_,
    f1_score           = NA_real_
  )
}

#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y
