#' @title Binary Panel Logistic Engine for GWPR.light 1.0.0
#'
#' @description
#' Internal functions for fitting Geographically Weighted Panel Regression
#' with a binary (binomial) response.  Supports pooling (stats::glm),
#' fixed effects (fixest::feglm), and random effects (glmmTMB::glmmTMB)
#' panel models, plus individual, time, two-way, and nested effects.
#'
#' The \code{family} parameter is reserved for future multi-class extension;
#' in version 1.0.0 only \code{"binomial"} is supported.
#'
#' @name model_logistic
#' @keywords internal
#' @importFrom stats binomial
NULL

# ---------------------------------------------------------------------------
# validate_binary_response
# ---------------------------------------------------------------------------

#' Validate that a response vector is suitable for binary logistic regression
#'
#' Stops with an informative error when \code{y} is a factor with more than two
#' levels.  Passes through numeric 0/1 vectors and two-level factors unchanged.
#'
#' @param y A numeric, logical, or factor vector.
#'
#' @return Invisibly \code{TRUE} when validation passes.
#' @keywords internal
validate_binary_response <- function(y) {
  if (is.factor(y)) {
    lvls <- levels(y)
    if (length(lvls) > 2L) {
      stop(
        "Multi-class factor response is not supported in GWPR.light 1.0.0. ",
        "The response must be a two-level factor or a numeric 0/1 vector. ",
        "Found ", length(lvls), " factor levels.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# standardize_logistic_response
# ---------------------------------------------------------------------------

#' Coerce a binary response to a numeric 0/1 integer vector
#'
#' Converts logical and two-level factor responses to 0/1 integer.
#' Numeric 0/1 vectors are returned unchanged.
#' Factors with more than two levels raise an error via
#' \code{validate_binary_response()}.
#'
#' @param y A numeric, logical, or factor vector.
#'
#' @return An integer vector of 0s and 1s.
#' @keywords internal
standardize_logistic_response <- function(y) {
  validate_binary_response(y)

  if (is.logical(y)) return(as.integer(y))

  if (is.factor(y)) {
    lvls <- levels(y)
    return(as.integer(y) - 1L)
  }

  if (is.numeric(y)) {
    unique_vals <- unique(y[!is.na(y)])
    invalid     <- setdiff(unique_vals, c(0, 1))
    if (length(invalid) > 0L) {
      stop(
        "Numeric response must contain only 0 and 1 for binary logistic regression. ",
        "Found: ", paste(invalid, collapse = ", "),
        call. = FALSE
      )
    }
    return(as.integer(y))
  }

  stop(
    "Response must be a numeric 0/1 vector, a logical vector, ",
    "or a two-level factor.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# effect_to_feglm_fml
# ---------------------------------------------------------------------------

#' Build a fixest::feglm formula string from the user-facing effect string
#'
#' @param base_formula A formula object (without fixed-effect terms).
#' @param effect       Character scalar: one of \code{"individual"},
#'   \code{"time"}, \code{"two-way"}, \code{"nested"}.
#' @param id_col       Name of the individual ID column.
#' @param time_col     Name of the time column.
#'
#' @return A formula suitable for \code{fixest::feglm()}.
#' @keywords internal
effect_to_feglm_fml <- function(base_formula, effect, id_col, time_col) {
  fe_part <- switch(
    effect,
    individual = id_col,
    time       = time_col,
    `two-way`  = paste(id_col, time_col, sep = " + "),
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

  # Build feglm formula: lhs ~ rhs | fe_part
  lhs  <- deparse(base_formula[[2L]])
  rhs  <- deparse(base_formula[[3L]])
  stats::as.formula(paste0(lhs, " ~ ", rhs, " | ", fe_part))
}

# ---------------------------------------------------------------------------
# effect_to_glmmtmb_fml
# ---------------------------------------------------------------------------

#' Append random-effect terms to a base formula for glmmTMB
#'
#' @param base_formula A formula object (without random-effect terms).
#' @param effect       Character scalar: one of \code{"individual"},
#'   \code{"time"}, \code{"two-way"}, \code{"nested"}.
#' @param id_col       Name of the individual ID column.
#' @param time_col     Name of the time column.
#'
#' @return A formula suitable for \code{glmmTMB::glmmTMB()}.
#' @keywords internal
effect_to_glmmtmb_fml <- function(base_formula, effect, id_col, time_col) {
  re_part <- switch(
    effect,
    individual = paste0("(1 | ", id_col, ")"),
    time       = paste0("(1 | ", time_col, ")"),
    `two-way`  = paste0("(1 | ", id_col, ") + (1 | ", time_col, ")"),
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

  lhs  <- deparse(base_formula[[2L]])
  rhs  <- deparse(base_formula[[3L]])
  stats::as.formula(paste0(lhs, " ~ ", rhs, " + ", re_part))
}

# ---------------------------------------------------------------------------
# fit_logistic_pooling
# ---------------------------------------------------------------------------

#' Fit a pooled logistic model via stats::glm
#'
#' @param local_data A \code{data.frame} containing all formula variables.
#' @param formula    A formula object with a binary response.
#' @param weights    Numeric vector of kernel weights (same length as
#'   \code{nrow(local_data)}).
#' @param family     Character scalar reserved for future extension;
#'   currently only \code{"binomial"} is supported.
#'
#' @return A \code{glm} object.
#' @keywords internal
fit_logistic_pooling <- function(local_data, formula, weights,
                                 family = "binomial") {
  .local <- local_data
  .local[[".gwpr_wt"]] <- weights
  .wts   <- weights

  do.call(
    stats::glm,
    list(
      formula = formula,
      data    = quote(.local),
      family  = binomial(link = "logit"),
      weights = .wts
    ),
    envir = environment()
  )
}

# ---------------------------------------------------------------------------
# fit_logistic_fixed
# ---------------------------------------------------------------------------

#' Fit a fixed-effects logistic model via fixest::feglm
#'
#' @param local_data A \code{data.frame} containing all formula variables plus
#'   the panel index columns.
#' @param formula    A formula object with a binary response (no fixed-effect
#'   terms; those are added automatically from \code{effect}).
#' @param effect     Character scalar: one of \code{"individual"},
#'   \code{"time"}, \code{"two-way"}, \code{"nested"}.
#' @param weights    Numeric vector of kernel weights.
#' @param id_col     Name of the individual ID column.
#' @param time_col   Name of the time column.
#' @param family     Character scalar reserved for future extension.
#'
#' @return A \code{fixest} object.
#' @keywords internal
fit_logistic_fixed <- function(local_data, formula, effect, weights,
                               id_col, time_col, family = "binomial") {
  fe_fml <- effect_to_feglm_fml(formula, effect, id_col, time_col)
  .local <- local_data
  .local[[".gwpr_wt"]] <- weights

  fixest::feglm(
    fml     = fe_fml,
    data    = .local,
    family  = binomial(link = "logit"),
    weights = ~.gwpr_wt
  )
}

# ---------------------------------------------------------------------------
# fit_logistic_random
# ---------------------------------------------------------------------------

#' Fit a random-effects logistic model via glmmTMB::glmmTMB
#'
#' @param local_data A \code{data.frame} containing all formula variables plus
#'   the panel index columns.
#' @param formula    A formula object with a binary response (no random-effect
#'   terms; those are added automatically from \code{effect}).
#' @param effect     Character scalar: one of \code{"individual"},
#'   \code{"time"}, \code{"two-way"}, \code{"nested"}.
#' @param weights    Numeric vector of kernel weights.
#' @param id_col     Name of the individual ID column.
#' @param time_col   Name of the time column.
#' @param family     Character scalar reserved for future extension.
#'
#' @return A \code{glmmTMB} object.
#' @keywords internal
fit_logistic_random <- function(local_data, formula, effect, weights,
                                id_col, time_col, family = "binomial") {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for random-effects logistic models. ",
      "Install it with: install.packages('glmmTMB')",
      call. = FALSE
    )
  }

  re_fml <- effect_to_glmmtmb_fml(formula, effect, id_col, time_col)
  .local <- local_data
  # Store weights as a plain column; glmmTMB with binomial(link="logit") and
  # a 0/1 response interprets `weights` as prior observation weights when the
  # response is a single-column 0/1 integer vector.  To avoid the
  # "non-integer #successes" warning we round to integers (proportional to
  # the kernel weights, scaled so the largest weight = 100).
  max_wt <- max(weights, na.rm = TRUE)
  if (is.finite(max_wt) && max_wt > 0) {
    int_wts <- as.integer(round(weights / max_wt * 100))
    int_wts[int_wts < 1L] <- 1L
  } else {
    int_wts <- rep(1L, length(weights))
  }
  .local[[".gwpr_wt"]] <- int_wts

  glmmTMB::glmmTMB(
    formula = re_fml,
    data    = .local,
    family  = binomial(link = "logit"),
    weights = .local[[".gwpr_wt"]]
  )
}

# ---------------------------------------------------------------------------
# fit_logistic_local_model
# ---------------------------------------------------------------------------

#' Fit a single local binary logistic panel model
#'
#' Dispatches to the correct backend (\code{pooling}, \code{fixed}, or
#' \code{random}) based on \code{model}, wrapping execution in a
#' \code{tryCatch} so that convergence failures or complete-separation errors
#' are captured rather than propagated.
#'
#' @param local_data  A \code{data.frame} with all formula variables and panel
#'   index columns.
#' @param formula     A formula object.
#' @param model       Character scalar: \code{"pooling"}, \code{"fixed"}, or
#'   \code{"random"}.
#' @param effect      Character scalar: \code{"individual"}, \code{"time"},
#'   \code{"two-way"}, or \code{"nested"}.
#' @param weights     Numeric vector of kernel weights aligned with
#'   \code{nrow(local_data)}.
#' @param index       Character vector of length 2: \code{c(id_col, time_col)}.
#' @param threshold   Numeric scalar; classification threshold (default 0.5).
#' @param family      Character scalar reserved for future extension (currently
#'   only \code{"binomial"} is implemented).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{fit}}{The fitted model object, or \code{NULL} on failure.}
#'   \item{\code{status}}{\code{"ok"} or \code{"failed"}.}
#'   \item{\code{error}}{\code{NULL} or character error message.}
#'   \item{\code{metadata}}{Named list with additional fitting metadata.}
#' }
#' @keywords internal
fit_logistic_local_model <- function(local_data, formula, model, effect,
                                     weights, index,
                                     threshold = 0.5,
                                     family = "binomial") {
  model  <- match.arg(model,  choices = c("pooling", "fixed", "random"))
  effect <- match.arg(effect, choices = c("individual", "time", "two-way", "nested"))

  # Validate effect BEFORE tryCatch for nested/unknown — propagate to caller
  if (model != "pooling") {
    # This will stop() for nested or unknown effect
    if (model == "fixed") {
      effect_to_feglm_fml(formula, effect, index[1L], index[2L])
    } else {
      effect_to_glmmtmb_fml(formula, effect, index[1L], index[2L])
    }
  }

  # Validate response
  resp_name <- all.vars(formula)[1L]
  if (resp_name %in% names(local_data)) {
    validate_binary_response(local_data[[resp_name]])
  }

  metadata <- list()
  id_col   <- index[1L]
  time_col <- index[2L]

  result <- tryCatch({
    fit <- switch(
      model,
      pooling = fit_logistic_pooling(local_data, formula, weights, family),
      fixed   = fit_logistic_fixed(local_data, formula, effect, weights,
                                   id_col, time_col, family),
      random  = fit_logistic_random(local_data, formula, effect, weights,
                                    id_col, time_col, family)
    )
    list(fit = fit, status = "ok", error = NULL, metadata = metadata)
  }, error = function(e) {
    list(fit = NULL, status = "failed", error = conditionMessage(e),
         metadata = metadata)
  })

  result
}

# ---------------------------------------------------------------------------
# predict_logistic_local_model
# ---------------------------------------------------------------------------

#' Predict probabilities for a local logistic model
#'
#' Returns predicted probabilities (type = \code{"response"}) for the rows of
#' \code{local_data}, using the fitted model stored in \code{local_result}.
#'
#' @param local_result A list as returned by \code{fit_logistic_local_model()}.
#' @param local_data   A \code{data.frame} used for prediction.
#'
#' @return Numeric vector of probabilities in \code{[0, 1]}.  Returns a vector
#'   of \code{NA_real_} on failure.
#' @keywords internal
predict_logistic_local_model <- function(local_result, local_data) {
  if (local_result$status != "ok" || is.null(local_result$fit)) {
    return(rep(NA_real_, nrow(local_data)))
  }

  fit <- local_result$fit

  tryCatch({
    if (inherits(fit, "fixest")) {
      # fixest::predict returns vector of length equal to the training data
      pv <- stats::predict(fit, newdata = local_data, type = "response")
    } else if (inherits(fit, "glmmTMB")) {
      pv <- stats::predict(fit, newdata = local_data, type = "response",
                           allow.new.levels = TRUE)
    } else {
      # glm
      pv <- stats::predict(fit, newdata = local_data, type = "response")
    }
    as.numeric(pv)
  }, error = function(e) {
    rep(NA_real_, nrow(local_data))
  })
}

# ---------------------------------------------------------------------------
# extract_logistic_local_result
# ---------------------------------------------------------------------------

#' Extract coefficients and diagnostics from a local logistic model
#'
#' Computes predicted probabilities, predicted classes, Pearson residuals,
#' and extracts coefficient estimates.  Pearson residuals are defined as
#' \code{(y - p) / sqrt(p * (1 - p))} with \code{p} clipped to
#' \code{[eps, 1 - eps]} to avoid division by zero.
#'
#' @param local_result  A list as returned by \code{fit_logistic_local_model()}.
#' @param local_data    The \code{data.frame} used to fit the model (needed to
#'   extract \code{y} and to predict).
#' @param formula       The model formula (used to extract the response name).
#' @param threshold     Numeric scalar; classification threshold (default 0.5).
#' @param eps           Numeric scalar; clipping bound for probability (default
#'   \code{1e-15}).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{coefficients}}{Named numeric vector of local coefficient
#'     estimates, or \code{NA_real_} on failure.}
#'   \item{\code{prob}}{Numeric vector of predicted probabilities for the local
#'     data rows, or \code{NA_real_} on failure.}
#'   \item{\code{class_pred}}{Integer vector (0/1) of predicted classes, or
#'     \code{NA_real_} on failure.}
#'   \item{\code{pearson_resid}}{Numeric vector of Pearson residuals, or
#'     \code{NA_real_} on failure.}
#'   \item{\code{status}}{\code{"ok"} or \code{"failed"}.}
#'   \item{\code{error}}{\code{NULL} or character error message.}
#' }
#' @keywords internal
extract_logistic_local_result <- function(local_result, local_data, formula,
                                          threshold = 0.5, eps = 1e-15) {
  status <- local_result$status
  err    <- local_result$error
  n      <- nrow(local_data)

  if (status != "ok" || is.null(local_result$fit)) {
    return(list(
      coefficients  = NA_real_,
      prob          = rep(NA_real_, n),
      class_pred    = rep(NA_integer_, n),
      pearson_resid = rep(NA_real_, n),
      status        = "failed",
      error         = err
    ))
  }

  fit <- local_result$fit

  # Coefficients
  coefs <- tryCatch(stats::coef(fit), error = function(e) NA_real_)
  # For glmmTMB, coef() returns a list; extract the fixed part
  if (is.list(coefs)) {
    coefs <- tryCatch(
      {
        fe <- glmmTMB::fixef(fit)
        if (is.list(fe)) unlist(fe$cond) else unlist(fe)
      },
      error = function(e) NA_real_
    )
  }

  # Predicted probabilities
  prob <- predict_logistic_local_model(local_result, local_data)

  # Predicted class
  class_pred <- tryCatch(
    as.integer(prob >= threshold),
    error = function(e) rep(NA_integer_, n)
  )

  # Response y
  resp_name <- all.vars(formula)[1L]
  y_raw <- local_data[[resp_name]]
  y <- tryCatch(
    standardize_logistic_response(y_raw),
    error = function(e) rep(NA_integer_, n)
  )

  # Pearson residuals: (y - p) / sqrt(p * (1-p))
  pearson_resid <- tryCatch({
    p_clipped <- pmax(pmin(prob, 1 - eps), eps)
    (y - p_clipped) / sqrt(p_clipped * (1 - p_clipped))
  }, error = function(e) rep(NA_real_, n))

  list(
    coefficients  = coefs,
    prob          = prob,
    class_pred    = class_pred,
    pearson_resid = pearson_resid,
    status        = "ok",
    error         = NULL
  )
}

# ---------------------------------------------------------------------------
# fit_logistic_gwpr  (main entry point)
# ---------------------------------------------------------------------------

#' Fit a Geographically Weighted Panel Binary Logistic Regression
#'
#' Iterates over each spatial unit (focus), computes kernel weights, fits a
#' local panel logistic model, and collects coefficients, probabilities,
#' predicted classes, Pearson residuals, and summary metrics.
#'
#' Supports three panel model backends:
#' \itemize{
#'   \item \code{"pooling"} — \code{stats::glm(family = binomial)}.
#'   \item \code{"fixed"}   — \code{fixest::feglm()}.
#'   \item \code{"random"}  — \code{glmmTMB::glmmTMB(family = binomial)}.
#' }
#'
#' @param context         A \code{gwpr_context} list containing at minimum:
#'   \code{formula}, \code{model}, \code{effect}, \code{panel_data},
#'   \code{coords}, \code{id}, \code{time}, \code{kernel}, \code{adaptive},
#'   \code{threshold}, and \code{workers}.
#' @param bandwidth       Numeric scalar.  Fixed distance or adaptive
#'   (k-nearest-neighbour) bandwidth.
#' @param weights_context Optional pre-built distance context as returned by
#'   \code{build_distance_context()}.  If \code{NULL}, one is built from
#'   \code{context$coords}.
#' @param family          Character scalar reserved for future multi-class
#'   extension; currently only \code{"binomial"} is accepted (default).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{local_results}}{A list of per-focus extracted results (one per
#'     unique spatial unit).}
#'   \item{\code{predictions}}{Named numeric vector of predicted probabilities
#'     for each focus unit row (indexed by user ID).}
#'   \item{\code{class_pred}}{Named integer vector of predicted classes.}
#'   \item{\code{pearson_resid}}{Named numeric vector of Pearson residuals.}
#'   \item{\code{metrics}}{A list from \code{compute_logistic_metrics()} over
#'     all focus observations.}
#'   \item{\code{metadata}}{A list with fitting metadata and any per-focus
#'     warnings.}
#' }
#' @export
fit_logistic_gwpr <- function(context, bandwidth, weights_context = NULL,
                               family = "binomial") {
  if (!identical(family, "binomial")) {
    stop(
      "Only family = 'binomial' is supported in GWPR.light 1.0.0. ",
      "Multi-class logistic regression is reserved for a future release.",
      call. = FALSE
    )
  }

  formula    <- context$formula
  model      <- context$model
  effect     <- context$effect
  panel_data <- context$panel_data
  coords     <- context$coords
  id_col     <- context$id
  time_col   <- context$time
  kernel     <- context$kernel
  adaptive   <- context$adaptive
  threshold  <- .null_coalesce(context$threshold, 0.5)
  workers    <- .null_coalesce(context$workers, 1L)

  index <- c(id_col, time_col)

  # Build distance context if not supplied
  if (is.null(weights_context)) {
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
    distances <- get_local_distances(weights_context, focus_id = focus_id)
    kw        <- compute_kernel_weights(
      distance  = distances,
      bandwidth = bandwidth,
      kernel    = kernel,
      adaptive  = adaptive
    )
    names(kw) <- unique_ids

    row_weights <- kw[as.character(panel_data[[id_col]])]
    row_weights[is.na(row_weights)] <- 0

    local_res <- fit_logistic_local_model(
      local_data = panel_data,
      formula    = formula,
      model      = model,
      effect     = effect,
      weights    = row_weights,
      index      = index,
      threshold  = threshold,
      family     = family
    )

    extracted <- extract_logistic_local_result(
      local_result = local_res,
      local_data   = panel_data,
      formula      = formula,
      threshold    = threshold
    )

    # Subset to focus rows only
    focus_rows <- which(as.character(panel_data[[id_col]]) == focus_id)
    focus_data <- panel_data[focus_rows, , drop = FALSE]

    focus_prob  <- extracted$prob[focus_rows]
    focus_class <- extracted$class_pred[focus_rows]
    focus_resid <- extracted$pearson_resid[focus_rows]
    focus_y     <- tryCatch(
      standardize_logistic_response(focus_data[[all.vars(formula)[1L]]]),
      error = function(e) rep(NA_integer_, length(focus_rows))
    )

    list(
      focus_id   = focus_id,
      extracted  = extracted,
      y          = focus_y,
      prob       = focus_prob,
      class_pred = focus_class,
      resid      = focus_resid,
      metadata   = local_res$metadata
    )
  }

  # Run (serial or parallel)
  if (workers == 1L) {
    all_results <- lapply(unique_ids, fit_one_focus)
  } else {
    all_results <- parallel::mclapply(
      unique_ids, fit_one_focus,
      mc.cores = workers, mc.set.seed = TRUE
    )
  }

  names(all_results) <- unique_ids

  # Assemble global vectors
  all_y     <- unlist(lapply(all_results, `[[`, "y"),          use.names = FALSE)
  all_prob  <- unlist(lapply(all_results, `[[`, "prob"),       use.names = FALSE)
  all_class <- unlist(lapply(all_results, `[[`, "class_pred"), use.names = FALSE)
  all_resid <- unlist(lapply(all_results, `[[`, "resid"),      use.names = FALSE)

  focus_lengths <- vapply(all_results, function(r) length(r$y), integer(1L))
  id_names <- rep(unique_ids, focus_lengths)

  pred_named  <- stats::setNames(all_prob,  id_names)
  class_named <- stats::setNames(all_class, id_names)
  resid_named <- stats::setNames(all_resid, id_names)

  # Remove NAs (from failed local models) before computing global metrics
  keep <- !is.na(all_y) & !is.na(all_prob)
  metrics <- if (sum(keep) > 0L) {
    suppressWarnings(
      compute_logistic_metrics(all_y[keep], all_prob[keep], threshold = threshold)
    )
  } else {
    list(log_loss = NA_real_, accuracy = NA_real_, precision = NA_real_,
         recall = NA_real_, f1_score = NA_real_)
  }

  local_results <- lapply(all_results, `[[`, "extracted")
  names(local_results) <- unique_ids

  meta <- list(
    n_focuses          = length(unique_ids),
    n_failed           = sum(vapply(local_results,
                                   function(r) r$status == "failed",
                                   logical(1L))),
    per_focus_metadata = lapply(all_results, `[[`, "metadata")
  )

  list(
    local_results = local_results,
    predictions   = pred_named,
    class_pred    = class_named,
    pearson_resid = resid_named,
    metrics       = metrics,
    metadata      = meta
  )
}

# ---------------------------------------------------------------------------
# Utility: .null_coalesce  (avoids re-definition of %||% across files)
# ---------------------------------------------------------------------------

.null_coalesce <- function(x, y) if (is.null(x)) y else x
