#' @title Diagnostics Module for GWPR.light 1.0.0
#'
#' @description
#' Unified diagnostic interface for Geographically Weighted Panel Regression
#' models.  Provides Moran's I (spatial autocorrelation), local F test
#' (fixed vs. pooling), local Hausman test (fixed vs. random), and local
#' Breusch-Pagan LM test.
#'
#' @details
#' **Model applicability**
#'
#' | Diagnostic   | Linear | Logistic | Notes                                      |
#' |--------------|--------|----------|--------------------------------------------|
#' | moran        | yes    | yes      | Logistic uses Pearson residual              |
#' | f_test       | yes    | no       | Requires within and pooling models          |
#' | hausman      | yes    | no       | Only meaningful for random-effect models    |
#' | lm_test      | yes    | no       | Pooling or random-effect models             |
#'
#' **Panel balance**
#'
#' `compute_panel_moran()` is fully supported for balanced panels.  For
#' unbalanced panels a `warning()` is issued and the function attempts
#' computation using only the time periods present in every individual;
#' results may be unreliable.
#'
#' **Logistic interpretation**
#'
#' Moran's I computed from Pearson residuals of a Logistic model does not
#' follow the same asymptotic distribution as for linear models.  Treat the
#' test result as an exploratory heuristic, not a formal test.
#'
#' @name diagnostics
NULL

# ---------------------------------------------------------------------------
# compute_logistic_pearson_residual
# ---------------------------------------------------------------------------

#' Compute Pearson residuals for a binary logistic model
#'
#' Calculates the Pearson (standardised) residual for each observation:
#'
#' \deqn{r_i = \frac{y_i - p_i}{\sqrt{p_i (1 - p_i)}}}
#'
#' where \eqn{p_i} is the fitted probability clipped to
#' \eqn{[\texttt{eps},\; 1 - \texttt{eps}]} to prevent division by zero.
#'
#' @param y   Integer or numeric vector of true binary labels (0 or 1).
#' @param prob Numeric vector of fitted probabilities (same length as `y`).
#'   Values are clipped to `[eps, 1 - eps]`.
#' @param eps Numeric scalar; clipping bound (default `1e-15`).
#'
#' @return A numeric vector of Pearson residuals, the same length as `y`.
#'   No `Inf` or `NaN` values are produced as long as `eps > 0`.
#'
#' @details
#' **Applicable model**: binomial (logistic).
#'
#' **Logistic interpretation limit**: The Pearson residual is used here
#' only to provide a residual-like quantity for spatial diagnostics.  It is
#' not a formal goodness-of-fit statistic in the chi-squared sense unless the
#' model is correctly specified and observations are independent.
#'
#' @examples
#' y    <- c(1L, 0L, 1L, 0L)
#' prob <- c(0.9, 0.1, 0.8, 0.2)
#' compute_logistic_pearson_residual(y, prob)
#'
#' @export
compute_logistic_pearson_residual <- function(y, prob, eps = 1e-15) {
  if (!is.numeric(y) && !is.integer(y)) {
    stop("`y` must be a numeric or integer vector of 0/1 values.", call. = FALSE)
  }
  if (!is.numeric(prob)) {
    stop("`prob` must be a numeric vector of fitted probabilities.", call. = FALSE)
  }
  if (length(y) != length(prob)) {
    stop("`y` and `prob` must have the same length.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("`eps` must be a single positive numeric value.", call. = FALSE)
  }

  p <- pmax(pmin(as.numeric(prob), 1 - eps), eps)
  (as.numeric(y) - p) / sqrt(p * (1 - p))
}

# ---------------------------------------------------------------------------
# compute_panel_moran
# ---------------------------------------------------------------------------

#' Compute Moran's I for panel residuals
#'
#' Computes a panel-level Moran's I statistic by averaging Moran's I across
#' all time periods.  The variance and Z-statistic follow Beenstock &
#' Felsenstein (2019).
#'
#' @param residuals      Numeric matrix of residuals with individuals in rows
#'   and time periods in columns, **or** a numeric vector.  When a vector is
#'   supplied `panel_index` must be provided so that the vector can be reshaped
#'   into the required matrix.
#' @param spatial_weights Numeric n x n spatial weights matrix (rows should
#'   sum to 1 after row-standardisation).  Diagonal is set to 0 internally.
#'   `n` must equal the number of individuals (rows in the residual matrix).
#' @param panel_index    Optional data.frame or list with two columns / elements
#'   named `id` and `time` (in that order) that identify each element of
#'   `residuals` when it is a vector.  Must have the same length as the vector.
#'
#' @return A named list:
#' \describe{
#'   \item{`statistic`}{Z-score of the average Moran's I.}
#'   \item{`p_value`}{One-sided p-value (upper tail, alternative "greater").}
#'   \item{`estimated_I`}{Observed average Moran's I across time periods.}
#'   \item{`expected_I`}{Expected value of Moran's I under the null: \eqn{-1/(n-1)}.}
#'   \item{`variance`}{Variance of the average Moran's I.}
#'   \item{`alternative`}{`"greater"`.}
#'   \item{`n_individuals`}{Number of spatial units \eqn{n}.}
#'   \item{`n_periods`}{Number of time periods \eqn{T}.}
#' }
#'
#' @details
#' **Panel balance requirement**: Balanced panels are fully supported.  For
#' unbalanced panels a `warning()` is issued.  The function then uses only the
#' intersection of time periods present in every individual, which may discard
#' observations.
#'
#' **Applicable models**: both linear (model residuals) and logistic (Pearson
#' residuals via `compute_logistic_pearson_residual()`).
#'
#' **Logistic interpretation limit**: When used with Pearson residuals from a
#' logistic model the test is exploratory only.
#'
#' @references
#' Beenstock, M., Felsenstein, D. (2019). *The Econometric Analysis of
#' Non-Stationary Spatial Panel Data*. Springer.
#'
#' @examples
#' set.seed(1L)
#' n  <- 5L   # individuals
#' Ti <- 3L   # time periods
#' resid_mat <- matrix(rnorm(n * Ti), nrow = n, ncol = Ti)
#' W <- matrix(1 / (n - 1), nrow = n, ncol = n)
#' diag(W) <- 0
#' compute_panel_moran(resid_mat, W)
#'
#' @export
compute_panel_moran <- function(residuals, spatial_weights, panel_index = NULL) {

  # ---- reshape vector -> matrix if panel_index provided -------------------
  if (is.numeric(residuals) && !is.matrix(residuals)) {
    if (is.null(panel_index)) {
      stop(
        "When `residuals` is a vector, `panel_index` must be provided.",
        call. = FALSE
      )
    }
    idx <- if (is.data.frame(panel_index)) panel_index else as.data.frame(panel_index)
    if (ncol(idx) < 2L) {
      stop("`panel_index` must have at least two columns: id and time.", call. = FALSE)
    }
    colnames(idx)[1:2] <- c("pid", "tid")
    if (nrow(idx) != length(residuals)) {
      stop(
        "`panel_index` must have the same number of rows as `residuals` has elements.",
        call. = FALSE
      )
    }
    idx$resid <- residuals
    individuals <- unique(idx$pid)
    times       <- unique(idx$tid)

    # Check balance
    counts <- tapply(idx$resid, list(idx$pid, idx$tid), length)
    is_balanced <- all(counts == 1L, na.rm = TRUE) &&
      all(!is.na(counts))

    if (!is_balanced) {
      warning(
        "The panel is unbalanced. Moran's I for unbalanced panels is not ",
        "formally supported in version 1.0.0. Only the intersection of ",
        "time periods present in all individuals will be used. ",
        "Results may be unreliable.",
        call. = FALSE
      )
      # Keep only time periods present for every individual
      times_per_id <- tapply(idx$tid, idx$pid, unique)
      common_times <- Reduce(intersect, times_per_id)
      if (length(common_times) == 0L) {
        stop(
          "No common time periods found across all individuals. ",
          "Cannot compute panel Moran's I.",
          call. = FALSE
        )
      }
      idx <- idx[idx$tid %in% common_times, ]
      times <- common_times
    }

    # Build matrix: rows = individuals, cols = time periods
    n  <- length(individuals)
    Ti <- length(times)
    resid_mat <- matrix(NA_real_, nrow = n, ncol = Ti,
                        dimnames = list(as.character(individuals),
                                        as.character(times)))
    for (i in seq_along(individuals)) {
      id_rows <- idx[idx$pid == individuals[i], ]
      for (j in seq_along(times)) {
        hit <- id_rows[id_rows$tid == times[j], "resid"]
        if (length(hit) == 1L) resid_mat[i, j] <- hit
      }
    }
    residuals <- resid_mat
  }

  if (!is.matrix(residuals)) {
    stop("`residuals` must be a numeric matrix (n_individuals x n_periods).", call. = FALSE)
  }

  n  <- nrow(residuals)
  Ti <- ncol(residuals)

  if (!is.matrix(spatial_weights) || nrow(spatial_weights) != n || ncol(spatial_weights) != n) {
    stop(
      "`spatial_weights` must be an n x n numeric matrix where n equals ",
      "the number of rows in `residuals`.",
      call. = FALSE
    )
  }

  # Row-standardise and zero diagonal
  W <- spatial_weights
  diag(W) <- 0
  row_sums <- rowSums(W)
  zero_rows <- row_sums == 0
  if (any(zero_rows)) {
    warning(
      sum(zero_rows), " row(s) in `spatial_weights` sum to 0 after removing the diagonal. ",
      "Those units receive uniform weights.",
      call. = FALSE
    )
    row_sums[zero_rows] <- 1
  }
  W <- W / row_sums

  # Compute Moran's I for each time period
  I_vec <- numeric(Ti)
  for (t in seq_len(Ti)) {
    e      <- residuals[, t]
    ss     <- sum(e^2)
    if (ss == 0) {
      I_vec[t] <- 0
      next
    }
    # sum_ij w_ij * e_i * e_j  (vectorised)
    I_vec[t] <- as.numeric(t(e) %*% W %*% e) / ss
  }

  I_mean <- mean(I_vec)

  # Variance of average I (Beenstock & Felsenstein 2019 formula)
  V2_num   <- n * sum(W^2) + 3 * sum(W)^2 - n * sum(colSums(W)^2)
  V2_denom <- Ti * (n^2 - 1) * sum(W)^2
  V2       <- if (V2_denom == 0) NA_real_ else V2_num / V2_denom

  # Expected I
  E_I <- -1 / (n - 1)

  # Z-statistic
  ZI <- if (is.na(V2) || V2 <= 0) NA_real_ else (I_mean - E_I) / sqrt(V2)

  # One-sided p-value (greater)
  p_val <- if (is.na(ZI)) NA_real_ else stats::pnorm(ZI, lower.tail = FALSE)

  list(
    statistic    = ZI,
    p_value      = p_val,
    estimated_I  = I_mean,
    expected_I   = E_I,
    variance     = V2,
    alternative  = "greater",
    n_individuals = n,
    n_periods    = Ti
  )
}

# ---------------------------------------------------------------------------
# diagnose_moran
# ---------------------------------------------------------------------------

#' Run Moran's I diagnostic on a gwpr_fit object
#'
#' Extracts residuals from a fitted GWPR model (Pearson residuals for logistic
#' models, raw residuals for linear models) and computes the panel Moran's I
#' statistic.
#'
#' @param object A `gwpr_fit` object returned by `fit_gwpr()` or similar.
#' @param spatial_weights A row-standardised n x n spatial weights matrix.
#'   `n` must equal the number of spatial individuals in the fitted model.
#' @param panel_index A data.frame or list with columns/elements `id` and
#'   `time` that identify each element of `object$residuals`.
#' @param ...  Currently ignored.
#'
#' @return A named list compatible with the `diagnostics` slot of a
#'   `gwpr_diagnostics` object.  Contains the elements returned by
#'   `compute_panel_moran()` plus `residual_type`.
#'
#' @details
#' **Applicable models**: gaussian (linear residuals) and binomial (Pearson
#' residuals).
#'
#' **Panel balance**: See `compute_panel_moran()`.
#'
#' **Failure conditions**: Fails if `object` is not a `gwpr_fit`, if
#' `object$residuals` is `NULL`, or if `spatial_weights` dimensions do not
#' match the number of individuals.
#'
#' **Logistic interpretation limit**: Moran's I computed on Pearson residuals
#' is exploratory; the asymptotic distribution differs from the linear case.
#'
#' @examples
#' fit <- new_gwpr_fit(
#'   family    = "gaussian",
#'   model     = "within",
#'   effect    = "individual",
#'   bandwidth = 1,
#'   residuals = c(0.1, -0.2, 0.3, -0.1, 0.2, -0.3),
#'   metadata  = list(response = NULL, prob = NULL)
#' )
#' W <- matrix(c(0, 0.5, 0.5,
#'               0.5, 0, 0.5,
#'               0.5, 0.5, 0), nrow = 3, byrow = TRUE)
#' idx <- data.frame(id = c(1,1,2,2,3,3), time = c(1,2,1,2,1,2))
#' diagnose_moran(fit, W, idx)
#'
#' @export
diagnose_moran <- function(object, spatial_weights, panel_index, ...) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }
  if (is.null(object$residuals)) {
    stop("`object$residuals` is NULL; cannot compute Moran's I.", call. = FALSE)
  }

  family <- object$family %||% "gaussian"

  if (identical(family, "binomial")) {
    # Extract predicted probabilities and response from metadata
    prob <- object$metadata$prob
    y    <- object$metadata$response

    if (is.null(prob) || is.null(y)) {
      stop(
        "For logistic models, `object$metadata` must contain `prob` and `response`. ",
        "Cannot compute Pearson residuals.",
        call. = FALSE
      )
    }
    resid_vec     <- compute_logistic_pearson_residual(y, prob)
    residual_type <- "pearson"
  } else {
    resid_vec     <- as.numeric(object$residuals)
    residual_type <- "raw"
  }

  moran_res <- compute_panel_moran(resid_vec, spatial_weights, panel_index)
  moran_res$residual_type <- residual_type
  moran_res
}

# ---------------------------------------------------------------------------
# diagnose_local_f
# ---------------------------------------------------------------------------

#' Local F test diagnostic on a gwpr_fit object
#'
#' Performs a local F test (fixed effects vs. pooling) using per-unit local
#' residuals stored in the fitted model object.
#'
#' @param object A `gwpr_fit` object.
#' @param ...  Currently ignored.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`local_f`}{Data frame with columns `unit_id`, `statistic`, `p_value`,
#'     `df1`, `df2`, `status`.}
#'   \item{`n_tested`}{Number of units tested.}
#'   \item{`n_failed`}{Number of units where the test could not be computed.}
#' }
#'
#' @details
#' **Applicable models**: gaussian (linear).  Not applicable to logistic models;
#' returns a `status = "not_applicable"` result when `family = "binomial"`.
#'
#' **Panel balance requirement**: No constraint; the test uses per-unit local
#' residuals already computed during fitting.
#'
#' **Failure conditions**: If `local_results` is empty or missing, all units
#' are reported as failed.  If a unit's local result does not contain the
#' information needed (within and pooling residuals), that unit is reported as
#' failed with an informative `status`.
#'
#' **Logistic interpretation limit**: Not applicable; see above.
#'
#' @examples
#' local_res <- list(
#'   "1" = list(
#'     within_rss  = 2.0, within_df  = 5L,
#'     pooling_rss = 5.0, pooling_df = 8L,
#'     status = "ok"
#'   ),
#'   "2" = list(
#'     within_rss  = 1.5, within_df  = 5L,
#'     pooling_rss = 4.0, pooling_df = 8L,
#'     status = "ok"
#'   )
#' )
#' fit <- new_gwpr_fit(
#'   family       = "gaussian",
#'   model        = "within",
#'   local_results = local_res
#' )
#' diagnose_local_f(fit)
#'
#' @export
diagnose_local_f <- function(object, ...) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }

  family <- object$family %||% "gaussian"
  if (identical(family, "binomial")) {
    return(list(
      local_f    = data.frame(
        unit_id   = character(0),
        statistic = numeric(0),
        p_value   = numeric(0),
        df1       = integer(0),
        df2       = integer(0),
        status    = character(0),
        stringsAsFactors = FALSE
      ),
      n_tested   = 0L,
      n_failed   = 0L,
      status     = "not_applicable",
      message    = "Local F test is not applicable to logistic (binomial) models."
    ))
  }

  local_results <- object$local_results
  if (length(local_results) == 0L) {
    return(list(
      local_f  = data.frame(
        unit_id   = character(0),
        statistic = numeric(0),
        p_value   = numeric(0),
        df1       = integer(0),
        df2       = integer(0),
        status    = character(0),
        stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "no_local_results"
    ))
  }

  unit_ids  <- names(local_results)
  if (is.null(unit_ids)) unit_ids <- as.character(seq_along(local_results))

  rows <- vector("list", length(unit_ids))

  for (i in seq_along(unit_ids)) {
    uid <- unit_ids[[i]]
    res <- local_results[[i]]

    row <- list(unit_id = uid, statistic = NA_real_, p_value = NA_real_,
                df1 = NA_integer_, df2 = NA_integer_, status = "ok")

    if (!is.null(res$f_statistic) && !is.null(res$f_p_value)) {
      # Pre-computed values stored by model engine
      row$statistic <- as.numeric(res$f_statistic)
      row$p_value   <- as.numeric(res$f_p_value)
      row$df1       <- as.integer(res$f_df1 %||% NA_integer_)
      row$df2       <- as.integer(res$f_df2 %||% NA_integer_)
      row$status    <- res$status %||% "ok"
    } else if (!is.null(res$within_rss) && !is.null(res$pooling_rss)) {
      # Compute from RSS values
      rss_w  <- as.numeric(res$within_rss)
      rss_p  <- as.numeric(res$pooling_rss)
      df_w   <- as.integer(res$within_df)
      df_p   <- as.integer(res$pooling_df)
      df1    <- df_p - df_w
      df2    <- df_w

      if (df1 <= 0L || df2 <= 0L || rss_w <= 0) {
        row$status <- "insufficient_df"
      } else {
        f_stat       <- ((rss_p - rss_w) / df1) / (rss_w / df2)
        row$statistic <- f_stat
        row$p_value   <- stats::pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
        row$df1       <- df1
        row$df2       <- df2
        row$status    <- res$status %||% "ok"
      }
    } else {
      row$status <- "missing_local_f_data"
    }

    rows[[i]] <- row
  }

  out_df <- data.frame(
    unit_id   = vapply(rows, `[[`, character(1L), "unit_id"),
    statistic = vapply(rows, `[[`, numeric(1L),   "statistic"),
    p_value   = vapply(rows, `[[`, numeric(1L),   "p_value"),
    df1       = vapply(rows, `[[`, integer(1L),   "df1"),
    df2       = vapply(rows, `[[`, integer(1L),   "df2"),
    status    = vapply(rows, `[[`, character(1L), "status"),
    stringsAsFactors = FALSE
  )

  n_failed <- sum(out_df$status != "ok")

  list(local_f  = out_df,
       n_tested = nrow(out_df),
       n_failed = n_failed,
       status   = "ok")
}

# ---------------------------------------------------------------------------
# diagnose_hausman
# ---------------------------------------------------------------------------

#' Local Hausman test diagnostic on a gwpr_fit object
#'
#' Performs a local Hausman test (within vs. random) for each spatial unit
#' using test statistics pre-computed during model fitting or stored in
#' `local_results`.
#'
#' @param object A `gwpr_fit` object.
#' @param ...  Currently ignored.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`local_hausman`}{Data frame with columns `unit_id`, `statistic`,
#'     `p_value`, `df`, `status`.}
#'   \item{`n_tested`}{Number of units tested.}
#'   \item{`n_failed`}{Number of units where the test could not be computed.}
#'   \item{`status`}{Overall status: `"ok"`, `"not_applicable"`, or
#'     `"no_local_results"`.}
#' }
#'
#' @details
#' **Applicable models**: gaussian with `model = "random"`.  For pooling
#' models the Hausman test is not meaningful; the function returns a
#' `status = "not_applicable"` result.  For logistic models the function also
#' returns `status = "not_applicable"`.
#'
#' **Panel balance requirement**: No constraint at the unit level.
#'
#' **Failure conditions**: Returns `status = "missing_hausman_data"` for any
#' unit where the required statistics are absent from `local_results`.
#'
#' **Logistic interpretation limit**: Not applicable.
#'
#' @examples
#' local_res <- list(
#'   "1" = list(hausman_statistic = 3.2, hausman_p_value = 0.07,
#'              hausman_df = 2L, status = "ok"),
#'   "2" = list(hausman_statistic = 6.1, hausman_p_value = 0.01,
#'              hausman_df = 2L, status = "ok")
#' )
#' fit <- new_gwpr_fit(
#'   family        = "gaussian",
#'   model         = "random",
#'   local_results = local_res
#' )
#' diagnose_hausman(fit)
#'
#' @export
diagnose_hausman <- function(object, ...) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }

  family <- object$family %||% "gaussian"
  model  <- object$model  %||% "within"

  if (identical(family, "binomial")) {
    return(list(
      local_hausman = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "not_applicable",
      message  = "Hausman test is not applicable to logistic (binomial) models."
    ))
  }

  if (identical(model, "pooling")) {
    return(list(
      local_hausman = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "not_applicable",
      message  = paste0(
        "Hausman test compares fixed vs. random effects. ",
        "It is not meaningful for model = \"pooling\"."
      )
    ))
  }

  local_results <- object$local_results
  if (length(local_results) == 0L) {
    return(list(
      local_hausman = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "no_local_results"
    ))
  }

  unit_ids <- names(local_results)
  if (is.null(unit_ids)) unit_ids <- as.character(seq_along(local_results))

  rows <- vector("list", length(unit_ids))

  for (i in seq_along(unit_ids)) {
    uid <- unit_ids[[i]]
    res <- local_results[[i]]

    row <- list(unit_id = uid, statistic = NA_real_, p_value = NA_real_,
                df = NA_integer_, status = "ok")

    if (!is.null(res$hausman_statistic)) {
      row$statistic <- as.numeric(res$hausman_statistic)
      row$p_value   <- as.numeric(res$hausman_p_value %||% NA_real_)
      row$df        <- as.integer(res$hausman_df       %||% NA_integer_)
      row$status    <- res$status %||% "ok"
    } else {
      row$status <- "missing_hausman_data"
    }

    rows[[i]] <- row
  }

  out_df <- data.frame(
    unit_id   = vapply(rows, `[[`, character(1L), "unit_id"),
    statistic = vapply(rows, `[[`, numeric(1L),   "statistic"),
    p_value   = vapply(rows, `[[`, numeric(1L),   "p_value"),
    df        = vapply(rows, `[[`, integer(1L),   "df"),
    status    = vapply(rows, `[[`, character(1L), "status"),
    stringsAsFactors = FALSE
  )

  n_failed <- sum(out_df$status != "ok")

  list(local_hausman = out_df,
       n_tested      = nrow(out_df),
       n_failed      = n_failed,
       status        = "ok")
}

# ---------------------------------------------------------------------------
# diagnose_lm
# ---------------------------------------------------------------------------

#' Local Breusch-Pagan LM test diagnostic on a gwpr_fit object
#'
#' Performs a local Breusch-Pagan Lagrange Multiplier test for random effects
#' for each spatial unit using test statistics stored in `local_results`.
#'
#' @param object A `gwpr_fit` object.
#' @param ...  Currently ignored.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`local_lm`}{Data frame with columns `unit_id`, `statistic`,
#'     `p_value`, `df`, `status`.}
#'   \item{`n_tested`}{Number of units tested.}
#'   \item{`n_failed`}{Number of units where the test could not be computed.}
#'   \item{`status`}{Overall status.}
#' }
#'
#' @details
#' **Applicable models**: gaussian with `model = "pooling"` or
#' `model = "random"`.  For `within` models the test is not directly
#' applicable (it tests for random effects vs. OLS); the function returns
#' `status = "not_applicable"`.  For logistic models also not applicable.
#'
#' **Panel balance requirement**: No constraint at the unit level.
#'
#' **Failure conditions**: Returns `status = "missing_lm_data"` for units
#' missing the required statistics.
#'
#' **Logistic interpretation limit**: Not applicable.
#'
#' @examples
#' local_res <- list(
#'   "1" = list(lm_statistic = 4.1, lm_p_value = 0.04,
#'              lm_df = 1L, status = "ok"),
#'   "2" = list(lm_statistic = 1.2, lm_p_value = 0.27,
#'              lm_df = 1L, status = "ok")
#' )
#' fit <- new_gwpr_fit(
#'   family        = "gaussian",
#'   model         = "pooling",
#'   local_results = local_res
#' )
#' diagnose_lm(fit)
#'
#' @export
diagnose_lm <- function(object, ...) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }

  family <- object$family %||% "gaussian"
  model  <- object$model  %||% "pooling"

  if (identical(family, "binomial")) {
    return(list(
      local_lm = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "not_applicable",
      message  = "LM test is not applicable to logistic (binomial) models."
    ))
  }

  if (identical(model, "within")) {
    return(list(
      local_lm = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "not_applicable",
      message  = paste0(
        "The Breusch-Pagan LM test is most appropriate for pooling or random ",
        "effect models. It is not applicable to within (fixed-effects) models."
      )
    ))
  }

  local_results <- object$local_results
  if (length(local_results) == 0L) {
    return(list(
      local_lm = data.frame(
        unit_id   = character(0), statistic = numeric(0),
        p_value   = numeric(0),   df        = integer(0),
        status    = character(0), stringsAsFactors = FALSE
      ),
      n_tested = 0L,
      n_failed = 0L,
      status   = "no_local_results"
    ))
  }

  unit_ids <- names(local_results)
  if (is.null(unit_ids)) unit_ids <- as.character(seq_along(local_results))

  rows <- vector("list", length(unit_ids))

  for (i in seq_along(unit_ids)) {
    uid <- unit_ids[[i]]
    res <- local_results[[i]]

    row <- list(unit_id = uid, statistic = NA_real_, p_value = NA_real_,
                df = NA_integer_, status = "ok")

    if (!is.null(res$lm_statistic)) {
      row$statistic <- as.numeric(res$lm_statistic)
      row$p_value   <- as.numeric(res$lm_p_value %||% NA_real_)
      row$df        <- as.integer(res$lm_df       %||% NA_integer_)
      row$status    <- res$status %||% "ok"
    } else {
      row$status <- "missing_lm_data"
    }

    rows[[i]] <- row
  }

  out_df <- data.frame(
    unit_id   = vapply(rows, `[[`, character(1L), "unit_id"),
    statistic = vapply(rows, `[[`, numeric(1L),   "statistic"),
    p_value   = vapply(rows, `[[`, numeric(1L),   "p_value"),
    df        = vapply(rows, `[[`, integer(1L),   "df"),
    status    = vapply(rows, `[[`, character(1L), "status"),
    stringsAsFactors = FALSE
  )

  n_failed <- sum(out_df$status != "ok")

  list(local_lm = out_df,
       n_tested = nrow(out_df),
       n_failed = n_failed,
       status   = "ok")
}

# ---------------------------------------------------------------------------
# diagnose_gwpr  (top-level interface)
# ---------------------------------------------------------------------------

#' Run diagnostic tests on a fitted GWPR model
#'
#' Top-level interface that dispatches to individual diagnostic sub-functions.
#' Returns a `gwpr_diagnostics` object containing all requested test results.
#'
#' @param object     A `gwpr_fit` object returned by `fit_gwpr()` or similar.
#' @param diagnostics Character vector naming the tests to run.  Any subset of
#'   `c("moran", "f_test", "hausman", "lm_test")`.  Default is all four.
#' @param spatial_weights Required when `"moran"` is in `diagnostics`.
#'   A row-standardised n x n spatial weights matrix.
#' @param panel_index Required when `"moran"` is in `diagnostics`.
#'   A data.frame with columns `id` and `time` identifying each element of
#'   `object$residuals`.
#' @param ...  Additional arguments passed to individual diagnostic functions.
#'
#' @return A `gwpr_diagnostics` object (list with class `"gwpr_diagnostics"`)
#'   whose `diagnostics` slot contains the result of each requested test.
#'
#' @details
#' Tests that are not applicable to the fitted model (e.g., Hausman test on a
#' pooling model) return a list with `status = "not_applicable"` and an
#' explanatory `message`, rather than an error.
#'
#' @examples
#' fit <- new_gwpr_fit(
#'   family    = "gaussian",
#'   model     = "pooling",
#'   effect    = "individual",
#'   bandwidth = 1,
#'   residuals = c(0.1, -0.1, 0.2, -0.2, 0.3, -0.3)
#' )
#' diagnose_gwpr(fit, diagnostics = c("f_test", "hausman", "lm_test"))
#'
#' @export
diagnose_gwpr <- function(
    object,
    diagnostics    = c("moran", "f_test", "hausman", "lm_test"),
    spatial_weights = NULL,
    panel_index    = NULL,
    ...
) {
  if (!inherits(object, "gwpr_fit")) {
    stop("`object` must be a `gwpr_fit` object.", call. = FALSE)
  }

  diagnostics <- match.arg(
    diagnostics,
    choices  = c("moran", "f_test", "hausman", "lm_test"),
    several.ok = TRUE
  )

  results  <- list()
  warnings_acc <- character()

  # ---- Moran's I -----------------------------------------------------------
  if ("moran" %in% diagnostics) {
    if (is.null(spatial_weights) || is.null(panel_index)) {
      results[["moran"]] <- list(
        status  = "skipped",
        message = paste0(
          "Moran's I requires `spatial_weights` and `panel_index`. ",
          "Provide both arguments to run this test."
        )
      )
    } else {
      results[["moran"]] <- tryCatch(
        withCallingHandlers(
          diagnose_moran(object, spatial_weights, panel_index, ...),
          warning = function(w) {
            warnings_acc <<- c(warnings_acc, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          list(status = "error", message = conditionMessage(e))
        }
      )
    }
  }

  # ---- Local F test --------------------------------------------------------
  if ("f_test" %in% diagnostics) {
    results[["f_test"]] <- tryCatch(
      diagnose_local_f(object, ...),
      error = function(e) list(status = "error", message = conditionMessage(e))
    )
  }

  # ---- Hausman test --------------------------------------------------------
  if ("hausman" %in% diagnostics) {
    results[["hausman"]] <- tryCatch(
      diagnose_hausman(object, ...),
      error = function(e) list(status = "error", message = conditionMessage(e))
    )
  }

  # ---- LM test -------------------------------------------------------------
  if ("lm_test" %in% diagnostics) {
    results[["lm_test"]] <- tryCatch(
      diagnose_lm(object, ...),
      error = function(e) list(status = "error", message = conditionMessage(e))
    )
  }

  # ---- Panel balance check ------------------------------------------------
  panel_balance <- TRUE   # default; update when index available
  if (!is.null(panel_index)) {
    idx <- if (is.data.frame(panel_index)) panel_index else as.data.frame(panel_index)
    if (ncol(idx) >= 2L) {
      colnames(idx)[1:2] <- c("pid", "tid")
      counts    <- table(idx$pid)
      ref_count <- counts[[1L]]
      panel_balance <- all(counts == ref_count)
    }
  }

  new_gwpr_diagnostics(
    diagnostics   = results,
    model_type    = object$family  %||% "gaussian",
    panel_balance = panel_balance,
    warnings      = warnings_acc,
    metadata      = list(diagnostics_run = diagnostics)
  )
}
