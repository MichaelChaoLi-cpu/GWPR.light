#' @title Input Validation Functions for GWPR.light 1.0.0
#'
#' @description
#' These internal functions validate user inputs before any expensive computation,
#' providing early and consistent error messages.
#'
#' @name validation
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# validate_formula
# ---------------------------------------------------------------------------

#' Validate a model formula against a data frame
#'
#' @param formula A formula object.
#' @param data    A data frame containing the model variables.
#'
#' @return Invisibly returns `TRUE` when valid; stops with an informative message
#'   when a problem is detected.
#' @keywords internal
validate_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula object.", call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  # Extract response variable (left-hand side)
  lhs <- all.vars(formula[[2]])
  missing_lhs <- setdiff(lhs, names(data))
  if (length(missing_lhs) > 0) {
    stop(
      "Response variable(s) not found in `data`: ",
      paste(missing_lhs, collapse = ", "),
      call. = FALSE
    )
  }

  # Extract right-hand side variables
  rhs_vars <- all.vars(formula[[3]])
  # Remove bare "." which means all columns
  rhs_vars <- rhs_vars[rhs_vars != "."]
  missing_rhs <- setdiff(rhs_vars, names(data))
  if (length(missing_rhs) > 0) {
    stop(
      "Predictor variable(s) not found in `data`: ",
      paste(missing_rhs, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_panel_index
# ---------------------------------------------------------------------------

#' Validate panel index columns in a data frame
#'
#' @param data A data frame.
#' @param id   Name of the unit (individual) index column.
#' @param time Name of the time index column.
#'
#' @return Invisibly returns `TRUE` when valid; stops otherwise.
#' @keywords internal
validate_panel_index <- function(data, id, time) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(id) || length(id) != 1L) {
    stop("`id` must be a single character string.", call. = FALSE)
  }

  if (!is.character(time) || length(time) != 1L) {
    stop("`time` must be a single character string.", call. = FALSE)
  }

  if (!id %in% names(data)) {
    stop(
      "ID column '", id, "' not found in `data`.",
      call. = FALSE
    )
  }

  if (!time %in% names(data)) {
    stop(
      "Time column '", time, "' not found in `data`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_spatial
# ---------------------------------------------------------------------------

#' Validate a spatial sf object
#'
#' @param spatial An `sf` object representing spatial units.
#' @param id      Name of the ID column that must be present in `spatial`.
#'
#' @return Invisibly returns `TRUE` when valid; stops otherwise.
#' @keywords internal
validate_spatial <- function(spatial, id) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop(
      "Package 'sf' is required for spatial operations. ",
      "Please install it with: install.packages('sf')",
      call. = FALSE
    )
  }

  if (!inherits(spatial, "sf")) {
    stop(
      "`spatial` must be an sf object. ",
      "sp objects are not supported in GWPR.light 1.0.0. ",
      "Please convert using sf::st_as_sf().",
      call. = FALSE
    )
  }

  if (!is.character(id) || length(id) != 1L) {
    stop("`id` must be a single character string.", call. = FALSE)
  }

  if (!id %in% names(spatial)) {
    stop(
      "ID column '", id, "' not found in `spatial`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_family_response
# ---------------------------------------------------------------------------

#' Validate the response variable against the specified family
#'
#' @param data    A data frame.
#' @param formula A formula; the left-hand side is the response variable.
#' @param family  Character string: `"gaussian"` or `"binomial"`.
#'
#' @return Invisibly returns `TRUE` when valid; stops otherwise.
#' @keywords internal
validate_family_response <- function(data, formula, family) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula object.", call. = FALSE)
  }

  family <- match.arg(family, c("gaussian", "binomial"))

  lhs_name <- all.vars(formula[[2]])[[1]]

  if (!lhs_name %in% names(data)) {
    stop(
      "Response variable '", lhs_name, "' not found in `data`.",
      call. = FALSE
    )
  }

  y <- data[[lhs_name]]

  if (family == "gaussian") {
    if (!is.numeric(y)) {
      stop(
        "For family = 'gaussian', the response variable '", lhs_name,
        "' must be numeric.",
        call. = FALSE
      )
    }
  }

  if (family == "binomial") {
    if (is.factor(y)) {
      n_levels <- nlevels(y)
      if (n_levels > 2L) {
        stop(
          "For family = 'binomial', the response variable '", lhs_name,
          "' is a factor with ", n_levels, " levels. ",
          "GWPR.light 1.0.0 does not support multinomial response variables. ",
          "Please provide a binary (two-level) factor or 0/1 numeric.",
          call. = FALSE
        )
      }
      if (n_levels < 2L) {
        stop(
          "For family = 'binomial', the response variable '", lhs_name,
          "' must have exactly two levels.",
          call. = FALSE
        )
      }
    } else if (is.logical(y)) {
      # logical is fine (coercible to 0/1)
    } else if (is.numeric(y)) {
      unique_vals <- unique(y[!is.na(y)])
      invalid_vals <- setdiff(unique_vals, c(0, 1))
      if (length(invalid_vals) > 0) {
        stop(
          "For family = 'binomial', numeric response '", lhs_name,
          "' must contain only 0 and 1. ",
          "Found unexpected value(s): ",
          paste(invalid_vals, collapse = ", "),
          call. = FALSE
        )
      }
    } else {
      stop(
        "For family = 'binomial', the response variable '", lhs_name,
        "' must be 0/1 numeric, logical, or a two-level factor.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_workers
# ---------------------------------------------------------------------------

#' Validate the workers argument
#'
#' @param workers Number of parallel workers. Must be a positive integer.
#'
#' @return Invisibly returns `TRUE` when valid; stops otherwise.
#' @keywords internal
validate_workers <- function(workers) {
  if (!is.numeric(workers) || length(workers) != 1L) {
    stop("`workers` must be a single positive integer.", call. = FALSE)
  }

  if (is.na(workers) || workers != floor(workers) || workers < 1L) {
    stop(
      "`workers` must be a positive integer (>= 1). ",
      "Received: ", workers,
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_bandwidth_control
# ---------------------------------------------------------------------------

#' Validate bandwidth search control parameters
#'
#' @param method   Character: `"grid"`, `"sgd"`, or `"random"`.
#' @param control  Named list of control parameters.
#' @param adaptive Logical; if `TRUE`, bandwidth is a neighbour count.
#'
#' @return Invisibly returns `TRUE` when valid; stops otherwise.
#' @keywords internal
validate_bandwidth_control <- function(method, control, adaptive) {
  method <- match.arg(method, c("grid", "sgd", "random"))

  if (!is.list(control)) {
    stop("`control` must be a list.", call. = FALSE)
  }

  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("`adaptive` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  if (method == "grid") {
    required <- c("lower", "upper", "step")
    missing_keys <- setdiff(required, names(control))
    if (length(missing_keys) > 0) {
      stop(
        "Grid bandwidth search requires control list entries: ",
        paste(missing_keys, collapse = ", "),
        ". Please supply all of `lower`, `upper`, and `step`.",
        call. = FALSE
      )
    }

    lower <- control[["lower"]]
    upper <- control[["upper"]]
    step  <- control[["step"]]

    if (!is.numeric(lower) || length(lower) != 1L || is.na(lower)) {
      stop("`control$lower` must be a single numeric value.", call. = FALSE)
    }
    if (!is.numeric(upper) || length(upper) != 1L || is.na(upper)) {
      stop("`control$upper` must be a single numeric value.", call. = FALSE)
    }
    if (!is.numeric(step) || length(step) != 1L || is.na(step)) {
      stop("`control$step` must be a single numeric value.", call. = FALSE)
    }

    if (lower >= upper) {
      stop(
        "`control$lower` (", lower, ") must be less than `control$upper` (", upper, ").",
        call. = FALSE
      )
    }

    if (step <= 0) {
      stop(
        "`control$step` must be positive. Received: ", step,
        call. = FALSE
      )
    }

    if (adaptive) {
      if (lower != floor(lower) || lower < 1L) {
        stop(
          "For adaptive bandwidth grid search, `control$lower` must be a positive integer. ",
          "Received: ", lower,
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# validate_inputs  (top-level orchestrator)
# ---------------------------------------------------------------------------

#' Validate all inputs to the main GWPR functions
#'
#' This is the top-level validation entry point.  It calls the individual
#' `validate_*` helpers and also checks `model`, `effect`, and `kernel`.
#'
#' @param formula  A formula object.
#' @param data     A data frame.
#' @param spatial  An `sf` object.
#' @param id       Name of the unit ID column (character).
#' @param time     Name of the time column (character).
#' @param family   `"gaussian"` or `"binomial"`.
#' @param model    `"pooling"`, `"within"`, or `"random"`.
#' @param effect   `"individual"`, `"time"`, `"two-way"`, or `"nested"`.
#' @param kernel   One of `"bisquare"`, `"gaussian"`, `"exponential"`,
#'   `"tricube"`, `"boxcar"`.
#' @param adaptive Logical.
#' @param workers  Positive integer.
#'
#' @return Invisibly returns `TRUE` when all checks pass; stops otherwise.
#' @keywords internal
validate_inputs <- function(formula, data, spatial, id, time,
                            family  = c("gaussian", "binomial"),
                            model   = c("within", "pooling", "random"),
                            effect  = c("individual", "time", "two-way", "nested"),
                            kernel  = c("bisquare", "gaussian", "exponential",
                                        "tricube", "boxcar"),
                            adaptive = FALSE,
                            workers  = 1L) {
  # --- sub-validators ---
  validate_formula(formula, data)
  validate_panel_index(data, id, time)
  validate_spatial(spatial, id)
  validate_family_response(data, formula, family)
  validate_workers(workers)

  # --- model ---
  valid_models <- c("pooling", "within", "random")
  family_arg  <- match.arg(family,  c("gaussian", "binomial"))
  model_arg   <- match.arg(model,   valid_models)
  effect_arg  <- match.arg(effect,  c("individual", "time", "two-way", "nested"))
  kernel_arg  <- match.arg(kernel,  c("bisquare", "gaussian", "exponential",
                                       "tricube", "boxcar"))

  invisible(TRUE)
}
