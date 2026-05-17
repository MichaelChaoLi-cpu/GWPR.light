#' @title Linear GWPR Engine for GWPR.light 1.0.0
#'
#' @description
#' Internal functions for fitting Geographically Weighted Panel Regression
#' with a Gaussian (linear) response.  Supports pooling, within, and random
#' panel models, plus individual, time, two-way, and nested effects.
#'
#' @name model_linear
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# effect_to_plm
# ---------------------------------------------------------------------------

#' Map user-facing effect string to plm effect parameter
#'
#' @param effect Character scalar: one of `"individual"`, `"time"`,
#'   `"two-way"`, `"nested"`.
#'
#' @return Character scalar accepted by `plm::plm()` for its `effect`
#'   argument: `"individual"`, `"time"`, or `"twoways"`.  For `"nested"`,
#'   the function stops with an informative message because `plm` does not
#'   support nested effects without additional data conventions.
#' @keywords internal
effect_to_plm <- function(effect) {
  switch(
    effect,
    individual = "individual",
    time       = "time",
    `two-way`  = "twoways",
    nested     = stop(
      '`effect = "nested"` requires additional hierarchical level information ',
      "that is not available from `id` and `time` alone. ",
      "Please provide a nested-level column and use a custom formula.",
      call. = FALSE
    ),
    stop(
      sprintf('Unknown effect: "%s". Must be one of: individual, time, two-way, nested.', effect),
      call. = FALSE
    )
  )
}

# ---------------------------------------------------------------------------
# fit_linear_local_model
# ---------------------------------------------------------------------------

#' Fit a single local panel linear model
#'
#' Fits a geographically weighted panel linear model for one focal spatial
#' unit.  Errors are caught and returned as a structured failure result rather
#' than propagating to the caller.
#'
#' @param local_data  A `pdata.frame` or plain `data.frame` with columns for
#'   all formula variables plus the panel indices.
#' @param formula     A formula object.
#' @param model       Character scalar: `"pooling"`, `"within"`, or
#'   `"random"`.
#' @param effect      Character scalar: `"individual"`, `"time"`, `"two-way"`,
#'   or `"nested"`.
#' @param weights     Numeric vector of kernel weights aligned with the rows of
#'   `local_data`.
#' @param index       Character vector of length 2 giving the panel index
#'   column names: `c(id_col, time_col)`.
#' @param random_method Character scalar; estimation method for variance
#'   components when `model = "random"` (default `"swar"`).
#'
#' @return A list with elements:
#' \describe{
#'   \item{`fit`}{The fitted model object, or `NULL` on failure.}
#'   \item{`status`}{`"ok"` or `"failed"`.}
#'   \item{`error`}{`NULL` or character string with the error message.}
#'   \item{`metadata`}{Named list with additional model-fitting metadata,
#'     e.g. flagging single-observation individuals for within models.}
#' }
#' @keywords internal
fit_linear_local_model <- function(local_data, formula, model, effect,
                                   weights, index,
                                   random_method = "swar") {
  model  <- match.arg(model,  choices = c("pooling", "within", "random"))
  effect <- match.arg(effect, choices = c("individual", "time", "two-way", "nested"))

  # Validate effect BEFORE entering tryCatch so that nested/unknown effects
  # propagate to the caller rather than being silently caught.
  if (model != "pooling") {
    plm_effect <- effect_to_plm(effect)   # stops for "nested"
  } else {
    plm_effect <- NULL
  }

  metadata <- list()

  # Embed kernel weights into the data frame as a dedicated column.
  # This avoids lazy-evaluation scoping issues that arise when a function
  # parameter named `weights` is passed to lm()/plm(), whose own `weights`
  # argument triggers name-shadowing in the calling environment.
  .local <- local_data
  .local[[".gwpr_wt"]] <- weights
  .wts <- weights   # plain numeric copy in this scope

  result <- tryCatch({

    if (model == "pooling") {
      # Weighted OLS via lm() — use do.call to keep data visible
      fit <- do.call(
        stats::lm,
        list(formula = formula, data = quote(.local), weights = .wts),
        envir = environment()
      )

    } else {
      # plm_effect was already computed and validated above

      # Warn / record single-observation individuals for within models
      if (model == "within") {
        id_col  <- index[1]
        obs_per_id <- table(.local[[id_col]])
        singles <- names(obs_per_id[obs_per_id == 1L])
        if (length(singles) > 0L) {
          metadata$single_obs_individuals <- singles
        }
      }

      # Build pdata.frame
      .pdata <- plm::pdata.frame(
        .local, index = index,
        drop.index = FALSE, row.names = FALSE,
        stringsAsFactors = FALSE
      )
      .pwts <- .pdata[[".gwpr_wt"]]

      fit <- do.call(
        plm::plm,
        list(formula = formula, data = quote(.pdata), model = model,
             effect = plm_effect, index = index, weights = .pwts,
             random.method = random_method),
        envir = environment()
      )
    }

    list(fit = fit, status = "ok", error = NULL, metadata = metadata)

  }, error = function(e) {
    list(fit = NULL, status = "failed", error = conditionMessage(e),
         metadata = metadata)
  })

  result
}

# ---------------------------------------------------------------------------
# extract_linear_local_result
# ---------------------------------------------------------------------------

#' Extract coefficients and diagnostics from a local linear model
#'
#' @param local_result  A list as returned by `fit_linear_local_model()`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`coefficients`}{Named numeric vector of local coefficient estimates.}
#'   \item{`se`}{Named numeric vector of standard errors (same names).}
#'   \item{`tvalues`}{Named numeric vector of t-statistics.}
#'   \item{`local_r2`}{Numeric scalar or `NA_real_`.}
#'   \item{`local_aic`}{Numeric scalar or `NA_real_`.}
#'   \item{`status`}{`"ok"` or `"failed"`.}
#'   \item{`error`}{`NULL` or character error message.}
#' }
#' @keywords internal
extract_linear_local_result <- function(local_result) {
  status <- local_result$status
  err    <- local_result$error

  if (status != "ok" || is.null(local_result$fit)) {
    return(list(
      coefficients = NA_real_,
      se           = NA_real_,
      tvalues      = NA_real_,
      local_r2     = NA_real_,
      local_aic    = NA_real_,
      status       = "failed",
      error        = err
    ))
  }

  fit <- local_result$fit

  # lmtest::coeftest works for both lm and plm objects
  ct <- tryCatch(
    lmtest::coeftest(fit),
    error = function(e) NULL
  )

  if (is.null(ct)) {
    coefs  <- tryCatch(stats::coef(fit),    error = function(e) NA_real_)
    ses    <- NA_real_
    tvals  <- NA_real_
  } else {
    coefs <- ct[, 1L]
    ses   <- ct[, 2L]
    tvals <- ct[, 3L]
  }

  # Local R^2
  local_r2 <- tryCatch({
    if (inherits(fit, "plm")) {
      plm::r.squared(fit)
    } else {
      summary(fit)$r.squared
    }
  }, error = function(e) NA_real_)

  # Local AIC
  local_aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)

  list(
    coefficients = coefs,
    se           = ses,
    tvalues      = tvals,
    local_r2     = if (is.null(local_r2) || length(local_r2) == 0) NA_real_ else as.numeric(local_r2),
    local_aic    = if (is.null(local_aic) || length(local_aic) == 0) NA_real_ else as.numeric(local_aic),
    status       = "ok",
    error        = NULL
  )
}

# ---------------------------------------------------------------------------
# predict_linear_local_model
# ---------------------------------------------------------------------------

#' Predict response values for a local linear model
#'
#' Returns predicted values for the rows of `local_data`, using the fitted
#' model object stored in `local_result`.
#'
#' @param local_result A list as returned by `fit_linear_local_model()`.
#' @param local_data   The data frame used for prediction (same structure as
#'   the training data).
#'
#' @return Numeric vector of predicted values (same length as `nrow(local_data)`).
#'   Returns a vector of `NA_real_` on failure.
#' @keywords internal
predict_linear_local_model <- function(local_result, local_data) {
  if (local_result$status != "ok" || is.null(local_result$fit)) {
    return(rep(NA_real_, nrow(local_data)))
  }

  fit <- local_result$fit

  # For plm models, predict() issues a warning when newdata is not a
  # pdata.frame. Since the local model is fit on the full panel data,
  # extract fitted values and match to focus rows via the model frame rownames.
  if (inherits(fit, "plm")) {
    fv <- tryCatch(as.numeric(stats::fitted(fit)), error = function(e) NULL)
    if (!is.null(fv) && length(fv) == nrow(local_data)) {
      return(fv)
    }
    # Fallback: predict with pdata.frame conversion
    tryCatch({
      pd_new <- plm::pdata.frame(local_data,
                                 index = attr(fit$model, "index") |> names(),
                                 drop.index = FALSE, row.names = FALSE)
      suppressWarnings(as.numeric(stats::predict(fit, newdata = pd_new)))
    }, error = function(e) rep(NA_real_, nrow(local_data)))
  } else {
    tryCatch(
      as.numeric(stats::predict(fit, newdata = local_data)),
      error = function(e) rep(NA_real_, nrow(local_data))
    )
  }
}

# ---------------------------------------------------------------------------
# fit_linear_gwpr  (main entry point)
# ---------------------------------------------------------------------------

#' Fit a Geographically Weighted Panel Linear Regression
#'
#' Iterates over each spatial unit (focus), computes kernel weights,
#' fits a local panel linear model, and collects coefficients, predictions,
#' residuals, and summary metrics.
#'
#' @param context         A `gwpr_context` list that must contain at minimum:
#'   `formula`, `model`, `effect`, `panel_data`, `coords`, `id`, `time`,
#'   `kernel`, `adaptive`, and `workers`.  Optionally `seed`.
#' @param bandwidth       Numeric scalar.  Fixed distance or adaptive
#'   (k-nearest-neighbour) bandwidth.
#' @param weights_context Optional pre-built distance context as returned by
#'   `build_distance_context()`.  If `NULL`, one is built from
#'   `context$coords`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`local_results`}{A list of per-focus extracted results (one per
#'     unique spatial unit).}
#'   \item{`predictions`}{Named numeric vector of focus-unit predictions
#'     (indexed by the user ID).}
#'   \item{`residuals`}{Named numeric vector of residuals (y - y_hat).}
#'   \item{`metrics`}{A list from `compute_linear_metrics()` computed over
#'     all focus observations.}
#'   \item{`metadata`}{A list with fitting metadata and any per-focus
#'     warnings.}
#' }
#' @export
fit_linear_gwpr <- function(context, bandwidth, weights_context = NULL) {

  formula  <- context$formula
  model    <- context$model
  effect   <- context$effect
  panel_data <- context$panel_data
  coords   <- context$coords
  id_col   <- context$id
  time_col <- context$time
  kernel   <- context$kernel
  adaptive <- context$adaptive
  workers  <- context$workers %||% 1L
  seed     <- context$seed

  index <- c(id_col, time_col)

  # Build distance context if not supplied
  if (is.null(weights_context)) {
    # coords rows correspond to unique spatial units
    unit_ids <- rownames(coords)
    if (is.null(unit_ids)) {
      unit_ids <- as.character(seq_len(nrow(coords)))
    }
    weights_context <- build_distance_context(
      coords  = coords,
      ids     = unit_ids,
      longlat = FALSE,
      cache   = TRUE
    )
  }

  unique_ids <- weights_context$ids

  # Helper: fit one focus unit
  fit_one_focus <- function(focus_id) {
    # Kernel weights (one per spatial unit, not per panel row)
    distances <- get_local_distances(weights_context, focus_id = focus_id)
    kw        <- compute_kernel_weights(
      distance  = distances,
      bandwidth = bandwidth,
      kernel    = kernel,
      adaptive  = adaptive
    )
    names(kw) <- unique_ids

    # Join kernel weights to panel rows based on id column
    row_weights <- kw[as.character(panel_data[[id_col]])]
    row_weights[is.na(row_weights)] <- 0

    # Fit local model
    local_res <- fit_linear_local_model(
      local_data = panel_data,
      formula    = formula,
      model      = model,
      effect     = effect,
      weights    = row_weights,
      index      = index
    )

    # Extract coefficients etc.
    extracted <- extract_linear_local_result(local_res)

    # Predict for focus rows only
    focus_rows  <- which(as.character(panel_data[[id_col]]) == focus_id)
    focus_data  <- panel_data[focus_rows, , drop = FALSE]
    preds       <- predict_linear_local_model(local_res, focus_data)
    y_focus     <- focus_data[[all.vars(formula)[1L]]]

    list(
      focus_id   = focus_id,
      extracted  = extracted,
      y          = y_focus,
      y_hat      = preds,
      resid      = y_focus - preds,
      metadata   = local_res$metadata
    )
  }

  # Run (serial for workers == 1, parallel otherwise)
  if (workers == 1L) {
    all_results <- lapply(unique_ids, fit_one_focus)
  } else {
    if (!is.null(seed)) set.seed(seed)
    all_results <- parallel::mclapply(
      unique_ids, fit_one_focus,
      mc.cores = workers, mc.set.seed = TRUE
    )
  }

  names(all_results) <- unique_ids

  # Assemble global predictions and residuals
  all_y    <- unlist(lapply(all_results, `[[`, "y"),    use.names = FALSE)
  all_yhat <- unlist(lapply(all_results, `[[`, "y_hat"), use.names = FALSE)
  all_resid <- unlist(lapply(all_results, `[[`, "resid"), use.names = FALSE)

  pred_named  <- stats::setNames(all_yhat,  rep(unique_ids, vapply(all_results, function(r) length(r$y), integer(1L))))
  resid_named <- stats::setNames(all_resid, names(pred_named))

  metrics <- compute_linear_metrics(all_y, all_yhat)

  local_results <- lapply(all_results, `[[`, "extracted")
  names(local_results) <- unique_ids

  meta <- list(
    n_focuses           = length(unique_ids),
    n_failed            = sum(vapply(local_results, function(r) r$status == "failed", logical(1L))),
    per_focus_metadata  = lapply(all_results, `[[`, "metadata")
  )

  list(
    local_results = local_results,
    predictions   = pred_named,
    residuals     = resid_named,
    metrics       = metrics,
    metadata      = meta
  )
}

# ---------------------------------------------------------------------------
# Utility: %||%  (null coalescing, local to avoid conflicts)
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x
