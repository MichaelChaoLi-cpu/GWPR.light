#' @title Public API for GWPR.light 1.0.0
#'
#' @description
#' High-level user-facing functions for Geographically Weighted Panel
#' Regression.  These four functions form the complete public interface;
#' all internal complexity is hidden behind them.
#'
#' \itemize{
#'   \item \code{\link{gwpr}} — full pipeline (bandwidth search + fitting +
#'     optional diagnostics).
#'   \item \code{\link{select_bandwidth}} — standalone bandwidth search.
#'   \item \code{\link{fit_gwpr}} — fit with a known bandwidth.
#'   \item \code{\link{diagnose_gwpr}} — run diagnostics on a fitted model.
#' }
#'
#' @name api
NULL

# ---------------------------------------------------------------------------
# Internal helper: null-coalescing operator
# ---------------------------------------------------------------------------

`.api_null_coalesce` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# gwpr
# ---------------------------------------------------------------------------

#' Fit a Geographically Weighted Panel Regression (main entry point)
#'
#' Orchestrates the complete GWPR pipeline: input validation, data preparation,
#' optional memory estimation, optional bandwidth search, model fitting, and
#' optional diagnostics.
#'
#' @param formula A \code{formula} object specifying the model (e.g.
#'   \code{y ~ x1 + x2}).
#' @param data A \code{data.frame} containing the panel data.  Must include
#'   the columns referenced by \code{id} and \code{time}.
#' @param spatial An \code{sf} object with one row per spatial unit.  Must
#'   include the column referenced by \code{id}.
#' @param id Character scalar.  Name of the unit (individual) ID column shared
#'   by \code{data} and \code{spatial}.
#' @param time Character scalar.  Name of the time-period column in \code{data}.
#' @param family Character scalar.  Model family: \code{"gaussian"} (default,
#'   linear GWPR) or \code{"binomial"} (binary panel logistic GWPR).
#' @param model Character scalar.  Panel model type: \code{"within"} (default),
#'   \code{"pooling"}, or \code{"random"}.
#' @param effect Character scalar.  Panel effect: \code{"individual"} (default),
#'   \code{"time"}, \code{"two-way"}, or \code{"nested"}.
#' @param bandwidth Numeric scalar or \code{NULL} (default).  When \code{NULL}
#'   the bandwidth is selected automatically via \code{select_bandwidth()}.
#' @param bandwidth_method Character scalar.  Method for automatic bandwidth
#'   search: \code{"sgd"} (default), \code{"grid"}, or \code{"random"}.
#'   Ignored when \code{bandwidth} is supplied.
#' @param bandwidth_control Named list of control parameters passed to the
#'   bandwidth search function.  See \code{search_bandwidth_grid()},
#'   \code{search_bandwidth_sgd()}, or \code{search_bandwidth_random()} for
#'   accepted fields.
#' @param kernel Character scalar.  Kernel function: \code{"bisquare"}
#'   (default), \code{"gaussian"}, \code{"exponential"}, \code{"tricube"}, or
#'   \code{"boxcar"}.
#' @param adaptive Logical scalar.  \code{FALSE} (default) uses a fixed
#'   distance bandwidth; \code{TRUE} uses an adaptive (k-nearest-neighbour)
#'   bandwidth.
#' @param threshold Numeric scalar.  Classification threshold for
#'   \code{family = "binomial"} (default \code{0.5}).
#' @param workers Positive integer.  Number of parallel workers.  \code{1}
#'   (default) uses serial execution; values \code{> 1} enable explicit
#'   parallelism.
#' @param seed Integer or \code{NULL}.  Random seed for reproducibility of any
#'   stochastic steps (bandwidth search, parallel RNG).
#' @param diagnostics Logical scalar.  When \code{TRUE} (default),
#'   \code{diagnose_gwpr()} is called after fitting and its results are stored
#'   in the returned object.
#' @param ... Additional arguments passed to the bandwidth search or fitting
#'   functions.
#'
#' @return A \code{gwpr_fit} object.  Key fields:
#' \describe{
#'   \item{\code{local_results}}{Per-unit local model results.}
#'   \item{\code{predictions}}{In-sample predicted values / probabilities.}
#'   \item{\code{residuals}}{Residuals or Pearson residuals.}
#'   \item{\code{metrics}}{Overall goodness-of-fit metrics.}
#'   \item{\code{spatial_results}}{Data frame of per-unit coefficients.}
#'   \item{\code{search}}{Bandwidth search result (\code{gwpr_bandwidth}), or
#'     \code{NULL} when bandwidth was supplied directly.}
#'   \item{\code{diagnostics}}{A \code{gwpr_diagnostics} object, or
#'     \code{NULL}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal linear GWPR with a fixed bandwidth
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(
#'   id   = rep(1:4, each = 5),
#'   time = rep(1:5, 4),
#'   y    = rnorm(20),
#'   x1   = rnorm(20)
#' )
#' fit <- gwpr(y ~ x1, data = dat, spatial = pts, id = "id", time = "time",
#'             bandwidth = 2, diagnostics = FALSE, workers = 1)
#' print(fit)
#' }
#'
#' @export
gwpr <- function(formula,
                 data,
                 spatial,
                 id,
                 time,
                 family             = c("gaussian", "binomial"),
                 model              = c("within", "pooling", "random"),
                 effect             = c("individual", "time", "two-way", "nested"),
                 bandwidth          = NULL,
                 bandwidth_method   = c("sgd", "grid", "random"),
                 bandwidth_control  = list(),
                 kernel             = c("bisquare", "gaussian", "exponential",
                                        "tricube", "boxcar"),
                 adaptive           = FALSE,
                 threshold          = 0.5,
                 workers            = 1L,
                 seed               = NULL,
                 diagnostics        = TRUE,
                 ...) {
  cl <- match.call()

  # Match arguments
  family           <- match.arg(family)
  model            <- match.arg(model)
  effect           <- match.arg(effect)
  kernel           <- match.arg(kernel)
  bandwidth_method <- match.arg(bandwidth_method)

  # 1. Validate inputs
  validate_inputs(
    formula  = formula,
    data     = data,
    spatial  = spatial,
    id       = id,
    time     = time,
    family   = family,
    model    = model,
    effect   = effect,
    kernel   = kernel,
    adaptive = adaptive,
    workers  = workers
  )

  # 2. Build context
  ctx <- new_gwpr_context(
    call        = cl,
    formula     = formula,
    family      = family,
    model       = model,
    effect      = effect,
    id          = id,
    time        = time,
    kernel      = kernel,
    adaptive    = adaptive,
    threshold   = threshold,
    workers     = workers,
    seed        = seed,
    raw_data    = data,
    raw_spatial = spatial
  )

  # 3. Prepare data
  ctx <- prepare_data(ctx)

  # Enrich metadata with n_time and n_vars for memory estimation
  ctx$metadata$n_time <- length(unique(ctx$panel_data[[time]]))
  ctx$metadata$n_vars <- ncol(ctx$model_matrix)

  # 4. Memory estimate; warn on high risk
  mem <- tryCatch(
    estimate_memory(ctx, workers = workers, cache_distance = TRUE),
    error = function(e) NULL
  )
  if (!is.null(mem) && mem$risk == "high") {
    warning(format_memory_warning(mem), call. = FALSE)
  }

  # 5. Bandwidth selection (if not provided)
  bw_result <- NULL
  if (is.null(bandwidth)) {
    bw_result <- select_bandwidth(
      formula           = formula,
      data              = data,
      spatial           = spatial,
      id                = id,
      time              = time,
      family            = family,
      model             = model,
      effect            = effect,
      method            = bandwidth_method,
      control           = bandwidth_control,
      kernel            = kernel,
      adaptive          = adaptive,
      threshold         = threshold,
      workers           = workers,
      seed              = seed,
      ...
    )
    bandwidth <- bw_result$best_bandwidth
  }

  # 6. Fit the model
  fit_result <- fit_gwpr(
    formula   = formula,
    data      = data,
    spatial   = spatial,
    id        = id,
    time      = time,
    bandwidth = bandwidth,
    family    = family,
    model     = model,
    effect    = effect,
    kernel    = kernel,
    adaptive  = adaptive,
    threshold = threshold,
    workers   = workers,
    seed      = seed,
    ...
  )

  # Attach search result
  fit_result$search <- bw_result
  fit_result$call   <- cl

  # 7. Diagnostics
  if (isTRUE(diagnostics)) {
    diag_result <- tryCatch(
      diagnose_gwpr(fit_result),
      error   = function(e) NULL,
      warning = function(w) {
        withCallingHandlers(
          diagnose_gwpr(fit_result),
          warning = function(w2) invokeRestart("muffleWarning")
        )
      }
    )
    fit_result$diagnostics_result <- diag_result
  }

  fit_result
}

# ---------------------------------------------------------------------------
# select_bandwidth
# ---------------------------------------------------------------------------

#' Select an optimal bandwidth for GWPR
#'
#' Validates inputs, prepares data, and dispatches to the appropriate bandwidth
#' search algorithm: \code{\link{search_bandwidth_grid}},
#' \code{\link{search_bandwidth_sgd}}, or \code{\link{search_bandwidth_random}}.
#'
#' @param formula  A \code{formula} object.
#' @param data     A \code{data.frame} with panel data.
#' @param spatial  An \code{sf} object.
#' @param id       Character scalar; unit ID column name.
#' @param time     Character scalar; time column name.
#' @param family   \code{"gaussian"} or \code{"binomial"}.
#' @param model    \code{"within"}, \code{"pooling"}, or \code{"random"}.
#' @param effect   \code{"individual"}, \code{"time"}, \code{"two-way"}, or
#'   \code{"nested"}.
#' @param method   Bandwidth search method: \code{"sgd"} (default),
#'   \code{"grid"}, or \code{"random"}.
#' @param control  Named list of search control parameters.
#' @param kernel   Kernel function name.
#' @param adaptive Logical; \code{TRUE} for adaptive (k-NN) bandwidth.
#' @param threshold Numeric; classification threshold (binomial only).
#' @param workers  Positive integer; number of parallel workers.
#' @param seed     Integer random seed, or \code{NULL}.
#' @param ...      Additional arguments (currently unused).
#'
#' @return A \code{gwpr_bandwidth} object with fields:
#' \describe{
#'   \item{\code{best_bandwidth}}{The selected bandwidth value.}
#'   \item{\code{best_score}}{The criterion value at the best bandwidth.}
#'   \item{\code{method}}{The search method used.}
#'   \item{\code{history}}{Search history data frame.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(
#'   id   = rep(1:4, each = 5),
#'   time = rep(1:5, 4),
#'   y    = rnorm(20),
#'   x1   = rnorm(20)
#' )
#' bw <- select_bandwidth(
#'   y ~ x1, data = dat, spatial = pts, id = "id", time = "time",
#'   method  = "grid",
#'   control = list(lower = 0.5, upper = 2, step = 0.5),
#'   workers = 1
#' )
#' bw$best_bandwidth
#' }
#'
#' @export
select_bandwidth <- function(formula,
                              data,
                              spatial,
                              id,
                              time,
                              family    = c("gaussian", "binomial"),
                              model     = c("within", "pooling", "random"),
                              effect    = c("individual", "time", "two-way", "nested"),
                              method    = c("sgd", "grid", "random"),
                              control   = list(),
                              kernel    = c("bisquare", "gaussian", "exponential",
                                            "tricube", "boxcar"),
                              adaptive  = FALSE,
                              threshold = 0.5,
                              workers   = 1L,
                              seed      = NULL,
                              ...) {
  family <- match.arg(family)
  model  <- match.arg(model)
  effect <- match.arg(effect)
  kernel <- match.arg(kernel)
  method <- match.arg(method)

  # Validate
  validate_inputs(
    formula  = formula,
    data     = data,
    spatial  = spatial,
    id       = id,
    time     = time,
    family   = family,
    model    = model,
    effect   = effect,
    kernel   = kernel,
    adaptive = adaptive,
    workers  = workers
  )
  validate_bandwidth_control(method = method, control = control,
                              adaptive = adaptive)

  # Build context
  ctx <- new_gwpr_context(
    formula     = formula,
    family      = family,
    model       = model,
    effect      = effect,
    id          = id,
    time        = time,
    kernel      = kernel,
    adaptive    = adaptive,
    threshold   = threshold,
    workers     = workers,
    seed        = seed,
    raw_data    = data,
    raw_spatial = spatial
  )
  ctx <- prepare_data(ctx)

  # Scorer: fits GWPR for a given bandwidth and returns score + metrics
  scorer <- .make_scorer(family = family, threshold = threshold)

  # Dispatch to search method
  bw_result <- switch(
    method,
    grid   = search_bandwidth_grid(ctx, control = control, scorer = scorer,
                                    workers = workers, seed = seed),
    sgd    = search_bandwidth_sgd(ctx,  control = control, scorer = scorer,
                                   workers = workers, seed = seed),
    random = search_bandwidth_random(ctx, control = control, scorer = scorer,
                                      workers = workers, seed = seed),
    stop("Unknown bandwidth method: ", method, call. = FALSE)
  )

  bw_result
}

# ---------------------------------------------------------------------------
# fit_gwpr
# ---------------------------------------------------------------------------

#' Fit GWPR with a given bandwidth
#'
#' Validates inputs, prepares data, builds spatial weights, and fits the
#' Geographically Weighted Panel Regression for the specified bandwidth.
#' Returns a \code{gwpr_fit} object.
#'
#' @param formula   A \code{formula} object.
#' @param data      A \code{data.frame} with panel data.
#' @param spatial   An \code{sf} object.
#' @param id        Character scalar; unit ID column name.
#' @param time      Character scalar; time column name.
#' @param bandwidth Numeric scalar.  The bandwidth to use (fixed distance or
#'   number of neighbours when \code{adaptive = TRUE}).
#' @param family    \code{"gaussian"} (default) or \code{"binomial"}.
#' @param model     \code{"within"} (default), \code{"pooling"}, or
#'   \code{"random"}.
#' @param effect    \code{"individual"} (default), \code{"time"},
#'   \code{"two-way"}, or \code{"nested"}.
#' @param kernel    Kernel function name (default \code{"bisquare"}).
#' @param adaptive  Logical; \code{FALSE} (default) for fixed bandwidth.
#' @param threshold Numeric; classification threshold (binomial only,
#'   default \code{0.5}).
#' @param workers   Positive integer; number of parallel workers (default 1).
#' @param seed      Integer random seed, or \code{NULL}.
#' @param ...       Currently unused.
#'
#' @return A \code{gwpr_fit} object.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- sf::st_as_sf(
#'   data.frame(id = 1:4, X = c(0,1,0,1), Y = c(0,0,1,1)),
#'   coords = c("X", "Y"), crs = NA_integer_
#' )
#' dat <- data.frame(
#'   id   = rep(1:4, each = 5),
#'   time = rep(1:5, 4),
#'   y    = rnorm(20),
#'   x1   = rnorm(20)
#' )
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' print(fit)
#' }
#'
#' @export
fit_gwpr <- function(formula,
                     data,
                     spatial,
                     id,
                     time,
                     bandwidth,
                     family    = c("gaussian", "binomial"),
                     model     = c("within", "pooling", "random"),
                     effect    = c("individual", "time", "two-way", "nested"),
                     kernel    = c("bisquare", "gaussian", "exponential",
                                   "tricube", "boxcar"),
                     adaptive  = FALSE,
                     threshold = 0.5,
                     workers   = 1L,
                     seed      = NULL,
                     ...) {
  family <- match.arg(family)
  model  <- match.arg(model)
  effect <- match.arg(effect)
  kernel <- match.arg(kernel)

  if (missing(bandwidth) || is.null(bandwidth)) {
    stop("`bandwidth` must be supplied to `fit_gwpr()`. ",
         "To search automatically use `gwpr()` or `select_bandwidth()`.",
         call. = FALSE)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth)) {
    stop("`bandwidth` must be a single finite numeric value.", call. = FALSE)
  }
  if (bandwidth <= 0) {
    stop("`bandwidth` must be positive.", call. = FALSE)
  }

  # Validate
  validate_inputs(
    formula  = formula,
    data     = data,
    spatial  = spatial,
    id       = id,
    time     = time,
    family   = family,
    model    = model,
    effect   = effect,
    kernel   = kernel,
    adaptive = adaptive,
    workers  = workers
  )

  # Build context and prepare data
  ctx <- new_gwpr_context(
    formula     = formula,
    family      = family,
    model       = model,
    effect      = effect,
    id          = id,
    time        = time,
    kernel      = kernel,
    adaptive    = adaptive,
    threshold   = threshold,
    workers     = workers,
    seed        = seed,
    raw_data    = data,
    raw_spatial = spatial
  )
  ctx <- prepare_data(ctx)

  # Fit by family
  if (family == "gaussian") {
    engine_result <- fit_linear_gwpr(ctx, bandwidth = bandwidth)

    spatial_res <- tryCatch(
      build_spatial_results(engine_result$local_results),
      error = function(e) NULL
    )

    fit <- new_gwpr_fit(
      call            = NULL,
      family          = family,
      formula         = formula,
      model           = model,
      effect          = effect,
      bandwidth       = bandwidth,
      search          = NULL,
      local_results   = engine_result$local_results,
      predictions     = engine_result$predictions,
      residuals       = engine_result$residuals,
      metrics         = engine_result$metrics,
      spatial_results = spatial_res,
      metadata        = c(ctx$metadata, engine_result$metadata),
      warnings        = ctx$warnings
    )

  } else {
    # binomial — map model names: "within" -> "fixed" for logistic engine
    logistic_model <- switch(
      model,
      within  = "fixed",
      pooling = "pooling",
      random  = "random"
    )
    ctx_logistic        <- ctx
    ctx_logistic$model  <- logistic_model

    engine_result <- fit_logistic_gwpr(ctx_logistic, bandwidth = bandwidth)

    spatial_res <- tryCatch(
      build_spatial_results(engine_result$local_results),
      error = function(e) NULL
    )

    fit <- new_gwpr_fit(
      call            = NULL,
      family          = family,
      formula         = formula,
      model           = model,
      effect          = effect,
      bandwidth       = bandwidth,
      search          = NULL,
      local_results   = engine_result$local_results,
      predictions     = engine_result$predictions,
      residuals       = engine_result$pearson_resid,
      metrics         = engine_result$metrics,
      spatial_results = spatial_res,
      metadata        = c(
        ctx$metadata,
        engine_result$metadata,
        list(
          response = as.numeric(ctx$response),
          prob     = as.numeric(engine_result$predictions)
        )
      ),
      warnings        = ctx$warnings
    )
  }

  fit
}

# ---------------------------------------------------------------------------
# diagnose_gwpr
# ---------------------------------------------------------------------------

#' Run diagnostics on a fitted GWPR model
#'
#' Wraps the individual diagnostic functions
#' (\code{\link{diagnose_moran}}, \code{\link{diagnose_local_f}},
#' \code{\link{diagnose_hausman}}, \code{\link{diagnose_lm}}) and collects
#' their results into a \code{gwpr_diagnostics} object.
#'
#' @param object      A \code{gwpr_fit} object returned by \code{\link{fit_gwpr}}
#'   or \code{\link{gwpr}}.
#' @param diagnostics Character vector of diagnostic tests to run.  One or
#'   more of \code{"moran"}, \code{"f_test"}, \code{"hausman"},
#'   \code{"lm_test"}.  Defaults to all four.
#' @param ...         Additional arguments forwarded to individual diagnostic
#'   functions (e.g. \code{spatial_weights}, \code{panel_index} for Moran's I).
#'
#' @return A \code{gwpr_diagnostics} object.
#'
#' @examples
#' \dontrun{
#' fit <- fit_gwpr(y ~ x1, data = dat, spatial = pts, id = "id",
#'                 time = "time", bandwidth = 2, workers = 1)
#' diag_result <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman"))
#' print(diag_result)
#' }
#'
#' @export
diagnose_gwpr <- function(object,
                           diagnostics = c("moran", "f_test", "hausman",
                                           "lm_test"),
                           ...) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }

  diagnostics <- match.arg(diagnostics, several.ok = TRUE)

  dots <- list(...)
  results  <- list()
  warnings_acc <- character()

  family        <- `.api_null_coalesce`(object$family, "gaussian")
  panel_balance <- `.api_null_coalesce`(object$metadata$panel_balanced, TRUE)

  for (test in diagnostics) {
    res <- switch(
      test,

      moran = tryCatch({
        sw <- dots$spatial_weights
        pi <- dots$panel_index
        if (is.null(sw) || is.null(pi)) {
          list(status  = "skipped",
               message = paste(
                 "Moran's I requires `spatial_weights` and `panel_index`.",
                 "Pass them via `...` to `diagnose_gwpr()`."
               ))
        } else {
          diagnose_moran(object, spatial_weights = sw, panel_index = pi, ...)
        }
      }, error = function(e) {
        list(status = "error", message = conditionMessage(e))
      }),

      f_test = tryCatch(
        diagnose_local_f(object, ...),
        error = function(e) list(status = "error", message = conditionMessage(e))
      ),

      hausman = tryCatch(
        diagnose_hausman(object, ...),
        error = function(e) list(status = "error", message = conditionMessage(e))
      ),

      lm_test = tryCatch(
        diagnose_lm(object, ...),
        error = function(e) list(status = "error", message = conditionMessage(e))
      ),

      # default: unknown test name
      list(status = "unknown_test", message = paste("Unknown diagnostic:", test))
    )

    results[[test]] <- res
  }

  new_gwpr_diagnostics(
    diagnostics   = results,
    model_type    = object$family,
    panel_balance = panel_balance,
    warnings      = warnings_acc
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build a scorer function for bandwidth search
#'
#' Returns a closure that fits the full GWPR engine for a given bandwidth and
#' returns the score (MSE for linear, log_loss for logistic), together with
#' aggregate metrics.
#'
#' @param family    Character; `"gaussian"` or `"binomial"`.
#' @param threshold Numeric; classification threshold (binomial only).
#'
#' @return A function with signature `scorer(context, bandwidth)`.
#' @keywords internal
.make_scorer <- function(family = "gaussian", threshold = 0.5) {
  function(context, bandwidth) {
    if (family == "gaussian") {
      engine_result <- fit_linear_gwpr(context, bandwidth = bandwidth)
      metrics  <- engine_result$metrics
      score    <- `.api_null_coalesce`(metrics$MSE, NA_real_)
      n_local  <- length(engine_result$local_results)
      n_failed <- sum(vapply(engine_result$local_results,
                             function(r) !is.null(r$status) && r$status == "failed",
                             logical(1L)))
      list(
        score                 = score,
        criterion             = "MSE",
        n_local_models        = n_local,
        n_failed_local_models = n_failed,
        metrics               = metrics,
        status                = "ok"
      )
    } else {
      # binomial — remap model
      ctx_logistic        <- context
      ctx_logistic$model  <- switch(
        context$model,
        within  = "fixed",
        pooling = "pooling",
        random  = "random"
      )
      engine_result <- fit_logistic_gwpr(ctx_logistic, bandwidth = bandwidth)
      metrics  <- engine_result$metrics
      score    <- `.api_null_coalesce`(metrics$log_loss, NA_real_)
      n_local  <- length(engine_result$local_results)
      n_failed <- sum(vapply(engine_result$local_results,
                             function(r) !is.null(r$status) && r$status == "failed",
                             logical(1L)))
      list(
        score                 = score,
        criterion             = "log_loss",
        n_local_models        = n_local,
        n_failed_local_models = n_failed,
        metrics               = metrics,
        status                = "ok"
      )
    }
  }
}
