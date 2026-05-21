#' @title Result Object Module for GWPR.light 1.0.0
#'
#' @description
#' S3 classes and constructor functions for the three result objects returned
#' by the GWPR.light public API: `gwpr_fit`, `gwpr_bandwidth`, and
#' `gwpr_diagnostics`.  Also provides `build_spatial_results()` for assembling
#' a data.frame that can be aligned with an `sf` geometry column.
#'
#' @name result
#' @importFrom stats quantile
NULL

# ---------------------------------------------------------------------------
# gwpr_fit
# ---------------------------------------------------------------------------

#' Construct a new gwpr_fit object
#'
#' Creates a standardised result object for a fitted geographically weighted
#' panel regression model.
#'
#' @param call       The matched call from the top-level API function, or
#'   `NULL`.
#' @param family     Character; model family: `"gaussian"` or `"binomial"`.
#' @param formula    A formula object, or `NULL`.
#' @param model      Character; panel model type: `"pooling"`, `"within"`, or
#'   `"random"`.
#' @param effect     Character; panel effect: `"individual"`, `"time"`,
#'   `"two-way"`, or `"nested"`.
#' @param bandwidth  Numeric scalar; the bandwidth used for fitting.
#' @param search     A `gwpr_bandwidth` object, or `NULL` when bandwidth was
#'   supplied directly.
#' @param global_model A summary / model object for the global (non-spatial)
#'   reference model, or `NULL`.
#' @param local_results A list of per-unit local model results.
#' @param predictions Numeric vector of in-sample predicted values, or `NULL`.
#' @param residuals  Numeric vector of residuals, or `NULL`.
#' @param metrics    Named list of overall goodness-of-fit metrics.
#' @param spatial_results A data.frame produced by `build_spatial_results()`,
#'   or `NULL`.
#' @param metadata   Named list of supplementary metadata.
#' @param warnings   Character vector of warnings accumulated during fitting.
#' @param ...        Additional named fields stored in the object.
#'
#' @return A named list of class `"gwpr_fit"`.
#'
#' @examples
#' fit <- new_gwpr_fit(
#'   family  = "gaussian",
#'   model   = "within",
#'   effect  = "individual",
#'   bandwidth = 100,
#'   metrics = list(R2 = 0.8, MSE = 0.1, RMSE = 0.316, MAE = 0.25),
#'   warnings = character()
#' )
#' print(fit)
#'
#' @noRd
new_gwpr_fit <- function(
    call           = NULL,
    family         = NULL,
    formula        = NULL,
    model          = NULL,
    effect         = NULL,
    bandwidth      = NULL,
    search         = NULL,
    global_model   = NULL,
    local_results  = list(),
    predictions    = NULL,
    residuals      = NULL,
    metrics        = list(),
    spatial_results = NULL,
    metadata       = list(),
    warnings       = character(),
    ...
) {
  extra <- list(...)

  obj <- c(
    list(
      call            = call,
      family          = family,
      formula         = formula,
      model           = model,
      effect          = effect,
      bandwidth       = bandwidth,
      search          = search,
      global_model    = global_model,
      local_results   = local_results,
      predictions     = predictions,
      residuals       = residuals,
      metrics         = metrics,
      spatial_results = spatial_results,
      metadata        = metadata,
      warnings        = warnings
    ),
    extra
  )

  structure(obj, class = "gwpr_fit")
}

#' Print a gwpr_fit object
#'
#' Displays a concise summary of a fitted GWPR model: family, panel model
#' type, effect, bandwidth, and top-level goodness-of-fit metrics.
#'
#' @param x   A `gwpr_fit` object.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' print(fit)
#' }
#'
#' @export
print.gwpr_fit <- function(x, ...) {
  cat("Geographically Weighted Panel Regression (gwpr_fit)\n")
  cat("----------------------------------------------------\n")
  cat("Family   :", x$family  %||% "<unknown>", "\n")
  cat("Model    :", x$model   %||% "<unknown>", "\n")
  cat("Effect   :", x$effect  %||% "<unknown>", "\n")
  cat("Bandwidth:", x$bandwidth %||% "<unknown>", "\n")

  if (length(x$metrics) > 0L) {
    cat("Metrics  :\n")
    for (nm in names(x$metrics)) {
      val <- x$metrics[[nm]]
      cat(sprintf("  %-12s %s\n", nm,
                  if (is.numeric(val)) formatC(val, digits = 4, format = "g")
                  else as.character(val)))
    }
  }

  if (length(x$warnings) > 0L) {
    cat("Warnings : (", length(x$warnings), ")\n", sep = "")
  }

  invisible(x)
}

#' Summarise a gwpr_fit object
#'
#' Prints the global model overview, quantile summary of local coefficients
#' (when available), and goodness-of-fit metrics.
#'
#' @param object A `gwpr_fit` object.
#' @param ...    Currently ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' summary(fit)
#' }
#'
#' @exportS3Method summary gwpr_fit
summary.gwpr_fit <- function(object, ...) {
  cat("=== GWPR Fit Summary ===\n\n")

  cat("Family   :", object$family  %||% "<unknown>", "\n")
  cat("Model    :", object$model   %||% "<unknown>", "\n")
  cat("Effect   :", object$effect  %||% "<unknown>", "\n")
  cat("Bandwidth:", object$bandwidth %||% "<unknown>", "\n\n")

  # Global model
  if (!is.null(object$global_model)) {
    cat("--- Global Model ---\n")
    print(object$global_model)
    cat("\n")
  }

  # Local coefficients quantile summary
  if (!is.null(object$spatial_results) && nrow(object$spatial_results) > 0L) {
    coef_cols <- grep("^coef_", names(object$spatial_results), value = TRUE)
    if (length(coef_cols) > 0L) {
      cat("--- Local Coefficient Quantiles ---\n")
      for (col in coef_cols) {
        vals <- object$spatial_results[[col]]
        if (is.numeric(vals)) {
          q <- quantile(vals, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
          cat(sprintf("  %-20s Min=%.4g  Q1=%.4g  Med=%.4g  Q3=%.4g  Max=%.4g\n",
                      col, q[1], q[2], q[3], q[4], q[5]))
        }
      }
      cat("\n")
    }
  }

  # Metrics
  if (length(object$metrics) > 0L) {
    cat("--- Metrics ---\n")
    for (nm in names(object$metrics)) {
      val <- object$metrics[[nm]]
      cat(sprintf("  %-12s %s\n", nm,
                  if (is.numeric(val)) formatC(val, digits = 4, format = "g")
                  else as.character(val)))
    }
    cat("\n")
  }

  if (length(object$warnings) > 0L) {
    cat("Warnings: (", length(object$warnings), " recorded)\n", sep = "")
  }

  invisible(object)
}

# ---------------------------------------------------------------------------
# gwpr_bandwidth
# ---------------------------------------------------------------------------

#' Construct a new gwpr_bandwidth object
#'
#' Creates a standardised result object for a bandwidth search.
#'
#' @param method           Character; search method: `"grid"`, `"sgd"`, or
#'   `"random"`.
#' @param best_bandwidth   Numeric scalar; the selected bandwidth.
#' @param best_score       Numeric scalar; the criterion score at the best
#'   bandwidth.
#' @param criterion        Character; name of the scoring criterion used.
#' @param history          A data.frame with one row per candidate / epoch.
#' @param metrics_history  A data.frame with metric values per candidate /
#'   epoch, or `NULL`.
#' @param seed             Integer random seed that was used, or `NULL`.
#' @param convergence_info Named list with convergence diagnostics, or `NULL`.
#' @param elapsed_time     Numeric scalar; total elapsed time in seconds.
#' @param warnings         Character vector of warnings.
#' @param ...              Additional named fields.
#'
#' @return A named list of class `"gwpr_bandwidth"`.
#'
#' @examples
#' bw <- new_gwpr_bandwidth(
#'   method         = "grid",
#'   best_bandwidth = 150,
#'   best_score     = 0.42,
#'   criterion      = "MSE",
#'   history        = data.frame(bandwidth = c(100, 150, 200),
#'                               score     = c(0.5, 0.42, 0.48)),
#'   elapsed_time   = 1.2
#' )
#' print(bw)
#'
#' @noRd
new_gwpr_bandwidth <- function(
    method           = NULL,
    best_bandwidth   = NULL,
    best_score       = NULL,
    criterion        = NULL,
    history          = NULL,
    metrics_history  = NULL,
    seed             = NULL,
    convergence_info = NULL,
    elapsed_time     = NULL,
    warnings         = character(),
    ...
) {
  extra <- list(...)

  obj <- c(
    list(
      method           = method,
      best_bandwidth   = best_bandwidth,
      best_score       = best_score,
      criterion        = criterion,
      history          = history,
      metrics_history  = metrics_history,
      seed             = seed,
      convergence_info = convergence_info,
      elapsed_time     = elapsed_time,
      warnings         = warnings
    ),
    extra
  )

  structure(obj, class = "gwpr_bandwidth")
}

#' Print a gwpr_bandwidth object
#'
#' Displays the search method, best bandwidth, criterion score, and number of
#' iterations explored.
#'
#' @param x   A `gwpr_bandwidth` object.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' bw <- select_bandwidth(y ~ x1, data = dat, spatial = pts, id = "id",
#'   time = "time", method = "grid",
#'   control = list(lower = 1, upper = 3, step = 1), workers = 1)
#' print(bw)
#' }
#'
#' @export
print.gwpr_bandwidth <- function(x, ...) {
  cat("Bandwidth Search Result (gwpr_bandwidth)\n")
  cat("----------------------------------------\n")
  cat("Method        :", x$method        %||% "<unknown>", "\n")
  cat("Criterion     :", x$criterion     %||% "<unknown>", "\n")
  cat("Best bandwidth:", x$best_bandwidth %||% "<unknown>", "\n")
  cat("Best score    :", x$best_score    %||% "<unknown>", "\n")

  if (!is.null(x$history) && is.data.frame(x$history)) {
    cat("Iterations    :", nrow(x$history), "\n")
  }

  if (!is.null(x$elapsed_time)) {
    cat("Elapsed (s)   :", round(x$elapsed_time, 3), "\n")
  }

  if (length(x$warnings) > 0L) {
    cat("Warnings      : (", length(x$warnings), ")\n", sep = "")
  }

  invisible(x)
}

#' Summarise a gwpr_bandwidth object
#'
#' Prints the search method, best bandwidth, and a brief history overview.
#'
#' @param object A `gwpr_bandwidth` object.
#' @param ...    Currently ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' bw <- select_bandwidth(y ~ x1, data = dat, spatial = pts, id = "id",
#'   time = "time", method = "grid",
#'   control = list(lower = 1, upper = 3, step = 1), workers = 1)
#' summary(bw)
#' }
#'
#' @exportS3Method summary gwpr_bandwidth
summary.gwpr_bandwidth <- function(object, ...) {
  cat("=== Bandwidth Search Summary ===\n\n")
  cat("Method        :", object$method        %||% "<unknown>", "\n")
  cat("Criterion     :", object$criterion     %||% "<unknown>", "\n")
  cat("Best bandwidth:", object$best_bandwidth %||% "<unknown>", "\n")
  cat("Best score    :", object$best_score    %||% "<unknown>", "\n")

  if (!is.null(object$seed)) {
    cat("Seed          :", object$seed, "\n")
  }

  if (!is.null(object$elapsed_time)) {
    cat("Elapsed (s)   :", round(object$elapsed_time, 3), "\n")
  }

  if (!is.null(object$history) && is.data.frame(object$history) &&
      nrow(object$history) > 0L) {
    cat("\n--- Search History (", nrow(object$history), " rows) ---\n", sep = "")
    print(utils::head(object$history, 6L))
    if (nrow(object$history) > 6L) {
      cat("  ... (", nrow(object$history) - 6L, " more rows)\n", sep = "")
    }
  }

  if (!is.null(object$convergence_info)) {
    cat("\n--- Convergence ---\n")
    for (nm in names(object$convergence_info)) {
      cat(sprintf("  %-20s %s\n", nm,
                  as.character(object$convergence_info[[nm]])))
    }
  }

  if (length(object$warnings) > 0L) {
    cat("\nWarnings: (", length(object$warnings), " recorded)\n", sep = "")
  }

  invisible(object)
}

# ---------------------------------------------------------------------------
# gwpr_diagnostics
# ---------------------------------------------------------------------------

#' Construct a new gwpr_diagnostics object
#'
#' Creates a standardised result object for model diagnostic tests.
#'
#' @param diagnostics     Named list of individual diagnostic test results.
#' @param model_type      Character; `"gaussian"` or `"binomial"`.
#' @param spatial_weights Spatial weights matrix or object, or `NULL`.
#' @param panel_balance   Logical; `TRUE` when the panel is balanced.
#' @param warnings        Character vector of warnings.
#' @param metadata        Named list of supplementary metadata.
#' @param ...             Additional named fields.
#'
#' @return A named list of class `"gwpr_diagnostics"`.
#'
#' @examples
#' diag_obj <- new_gwpr_diagnostics(
#'   diagnostics   = list(moran = list(statistic = 0.12, p_value = 0.03)),
#'   model_type    = "gaussian",
#'   panel_balance = TRUE,
#'   warnings      = character()
#' )
#' print(diag_obj)
#'
#' @noRd
new_gwpr_diagnostics <- function(
    diagnostics     = list(),
    model_type      = NULL,
    spatial_weights = NULL,
    panel_balance   = NULL,
    warnings        = character(),
    metadata        = list(),
    ...
) {
  extra <- list(...)

  obj <- c(
    list(
      diagnostics     = diagnostics,
      model_type      = model_type,
      spatial_weights = spatial_weights,
      panel_balance   = panel_balance,
      warnings        = warnings,
      metadata        = metadata
    ),
    extra
  )

  structure(obj, class = "gwpr_diagnostics")
}

#' Print a gwpr_diagnostics object
#'
#' Displays each diagnostic test name and, where available, its statistic and
#' p-value.
#'
#' @param x   A `gwpr_diagnostics` object.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' diag_obj <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman"))
#' print(diag_obj)
#' }
#'
#' @export
print.gwpr_diagnostics <- function(x, ...) {
  cat("GWPR Diagnostics (gwpr_diagnostics)\n")
  cat("------------------------------------\n")
  cat("Model type   :", x$model_type    %||% "<unknown>", "\n")
  cat("Panel balance:", x$panel_balance %||% "<unknown>", "\n")

  if (length(x$diagnostics) > 0L) {
    cat("Tests run    :", paste(names(x$diagnostics), collapse = ", "), "\n")
  } else {
    cat("Tests run    : (none)\n")
  }

  if (length(x$warnings) > 0L) {
    cat("Warnings     : (", length(x$warnings), ")\n", sep = "")
  }

  invisible(x)
}

#' Summarise a gwpr_diagnostics object
#'
#' Prints each diagnostic test result with statistic and p-value where
#' available.
#'
#' @param object A `gwpr_diagnostics` object.
#' @param ...    Currently ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(id = rep(1:4, each = 5), time = rep(1:5, 4),
#'                   y = rnorm(20), x1 = rnorm(20))
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' diag_obj <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman"))
#' summary(diag_obj)
#' }
#'
#' @exportS3Method summary gwpr_diagnostics
summary.gwpr_diagnostics <- function(object, ...) {
  cat("=== GWPR Diagnostics Summary ===\n\n")
  cat("Model type   :", object$model_type    %||% "<unknown>", "\n")
  cat("Panel balance:", object$panel_balance %||% "<unknown>", "\n\n")

  if (length(object$diagnostics) == 0L) {
    cat("No diagnostic tests recorded.\n")
  } else {
    for (test_name in names(object$diagnostics)) {
      res <- object$diagnostics[[test_name]]
      cat("---", test_name, "---\n")
      if (is.list(res)) {
        for (nm in names(res)) {
          cat(sprintf("  %-16s %s\n", nm, as.character(res[[nm]])))
        }
      } else {
        cat(" ", as.character(res), "\n")
      }
    }
  }

  if (length(object$warnings) > 0L) {
    cat("\nWarnings: (", length(object$warnings), " recorded)\n", sep = "")
    for (w in object$warnings) cat("  -", w, "\n")
  }

  invisible(object)
}

# ---------------------------------------------------------------------------
# build_spatial_results
# ---------------------------------------------------------------------------

#' Build a data.frame of spatial results
#'
#' Merges per-unit local coefficients, standard errors, t-statistics, local
#' fit quality metrics, and diagnostic statistics into a single data.frame
#' that can be aligned row-by-row with an `sf` geometry column.
#'
#' @param local_results A named list (one element per spatial unit) where each
#'   element is itself a named list with optional sub-lists:
#'   \describe{
#'     \item{`coefficients`}{Named numeric vector of local coefficients.}
#'     \item{`std_errors`}{Named numeric vector of local standard errors.}
#'     \item{`t_stats`}{Named numeric vector of local t-statistics.}
#'     \item{`metrics`}{Named list / numeric vector of local fit metrics.}
#'     \item{`diagnostics`}{Named list / numeric vector of local diagnostics.}
#'     \item{`status`}{Character scalar; `"ok"` or an error description.}
#'   }
#' @param id_map A named character or integer vector mapping unit IDs to
#'   positions, or `NULL`.  The names of `local_results` are used as unit IDs
#'   when `id_map` is `NULL`.
#' @param geometry An `sf` geometry column (`sfc`), a data.frame with the same
#'   number of rows as unique units, or `NULL`.  Used only for length
#'   validation; the geometry itself is not attached to the returned data.frame.
#'
#' @return A data.frame with one row per spatial unit.  Column names:
#'   \describe{
#'     \item{`unit_id`}{Unit identifier.}
#'     \item{`status`}{Fit status (`"ok"` or error description).}
#'     \item{`coef_<name>`}{Local coefficient for each predictor.}
#'     \item{`se_<name>`}{Local standard error for each predictor.}
#'     \item{`tstat_<name>`}{Local t-statistic for each predictor.}
#'     \item{`<metric name>`}{Local fit metric (e.g. `R2`, `MSE`, ...).}
#'     \item{`<diagnostic name>`}{Local diagnostic value.}
#'   }
#'
#' @examples
#' local_res <- list(
#'   "A" = list(
#'     coefficients = c(x1 = 0.5, x2 = 1.2),
#'     std_errors   = c(x1 = 0.1, x2 = 0.3),
#'     t_stats      = c(x1 = 5.0, x2 = 4.0),
#'     metrics      = list(R2 = 0.8, MSE = 0.05),
#'     status       = "ok"
#'   ),
#'   "B" = list(
#'     coefficients = c(x1 = 0.3, x2 = 0.9),
#'     std_errors   = c(x1 = 0.2, x2 = 0.4),
#'     t_stats      = c(x1 = 1.5, x2 = 2.25),
#'     metrics      = list(R2 = 0.6, MSE = 0.12),
#'     status       = "ok"
#'   )
#' )
#' build_spatial_results(local_res)
#'
#' @noRd
build_spatial_results <- function(local_results, id_map = NULL,
                                  geometry = NULL) {
  if (!is.list(local_results) || length(local_results) == 0L) {
    stop("`local_results` must be a non-empty list.", call. = FALSE)
  }

  unit_ids <- names(local_results)
  if (is.null(unit_ids)) {
    unit_ids <- as.character(seq_along(local_results))
  }

  # Validate geometry length if provided
  if (!is.null(geometry)) {
    geom_len <- if (is.data.frame(geometry)) nrow(geometry) else length(geometry)
    if (geom_len != length(unit_ids)) {
      stop(
        "`geometry` length (", geom_len, ") does not match the number of ",
        "units in `local_results` (", length(unit_ids), ").",
        call. = FALSE
      )
    }
  }

  rows <- vector("list", length(unit_ids))

  for (i in seq_along(unit_ids)) {
    uid <- unit_ids[[i]]
    res <- local_results[[i]]

    row <- list(unit_id = uid,
                status  = res$status %||% "ok")

    # Coefficients
    if (!is.null(res$coefficients) && length(res$coefficients) > 0L) {
      coef_names <- names(res$coefficients)
      if (is.null(coef_names)) {
        coef_names <- paste0("V", seq_along(res$coefficients))
      }
      for (j in seq_along(res$coefficients)) {
        row[[paste0("coef_", coef_names[j])]] <- res$coefficients[[j]]
      }
    }

    # Standard errors
    if (!is.null(res$std_errors) && length(res$std_errors) > 0L) {
      se_names <- names(res$std_errors)
      if (is.null(se_names)) {
        se_names <- paste0("V", seq_along(res$std_errors))
      }
      for (j in seq_along(res$std_errors)) {
        row[[paste0("se_", se_names[j])]] <- res$std_errors[[j]]
      }
    }

    # t-statistics
    if (!is.null(res$t_stats) && length(res$t_stats) > 0L) {
      ts_names <- names(res$t_stats)
      if (is.null(ts_names)) {
        ts_names <- paste0("V", seq_along(res$t_stats))
      }
      for (j in seq_along(res$t_stats)) {
        row[[paste0("tstat_", ts_names[j])]] <- res$t_stats[[j]]
      }
    }

    # Metrics
    if (!is.null(res$metrics) && length(res$metrics) > 0L) {
      for (nm in names(res$metrics)) {
        row[[nm]] <- res$metrics[[nm]]
      }
    }

    # Diagnostics
    if (!is.null(res$diagnostics) && length(res$diagnostics) > 0L) {
      for (nm in names(res$diagnostics)) {
        row[[nm]] <- res$diagnostics[[nm]]
      }
    }

    rows[[i]] <- row
  }

  # Separate character columns (unit_id, status) from numeric columns
  char_cols    <- c("unit_id", "status")
  all_row_cols <- unique(unlist(lapply(rows, names)))
  num_cols     <- setdiff(all_row_cols, char_cols)

  # Build numeric columns
  num_list <- lapply(
    stats::setNames(num_cols, num_cols),
    function(col) {
      vapply(rows, function(r) {
        v <- r[[col]]
        if (is.null(v)) NA_real_ else as.numeric(v)
      }, numeric(1L))
    }
  )

  # Build character columns
  unit_id_col <- unit_ids
  status_col  <- vapply(rows, function(r) r$status %||% "ok", character(1L))

  # Assemble final data.frame: unit_id, status first, then numeric columns
  out <- data.frame(
    unit_id = unit_id_col,
    status  = status_col,
    stringsAsFactors = FALSE
  )

  for (col in num_cols) {
    out[[col]] <- num_list[[col]]
  }

  out
}

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y
