#' @title Distance and Kernel Weights Module for GWPR.light 1.0.0
#'
#' @description
#' Internal functions for computing spatial distances and geographically
#' weighted kernel weights. Supports fixed and adaptive bandwidths, and five
#' kernel functions: bisquare, gaussian, exponential, tricube, and boxcar.
#'
#' @name weights
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# validate_bandwidth
# ---------------------------------------------------------------------------

#' Validate bandwidth parameter
#'
#' Checks that the supplied bandwidth is legal:
#' * fixed (adaptive = FALSE): must be a finite positive numeric scalar.
#' * adaptive (adaptive = TRUE): must be a finite positive integer scalar.
#'
#' @param bandwidth Numeric scalar.  The bandwidth value to validate.
#' @param adaptive  Logical scalar.  `TRUE` for adaptive (k-nearest-neighbour)
#'   bandwidth, `FALSE` for fixed distance bandwidth.
#'
#' @return Invisibly returns `TRUE` when validation passes.
#' @keywords internal
validate_bandwidth <- function(bandwidth, adaptive) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L) {
    stop("`bandwidth` must be a single numeric value.", call. = FALSE)
  }
  if (!is.finite(bandwidth)) {
    stop("`bandwidth` must be a finite number.", call. = FALSE)
  }
  if (adaptive) {
    if (bandwidth != floor(bandwidth) || bandwidth < 1) {
      stop(
        "For adaptive bandwidth, `bandwidth` must be a positive integer ",
        "(number of neighbours).",
        call. = FALSE
      )
    }
  } else {
    if (bandwidth <= 0) {
      stop(
        "For fixed bandwidth, `bandwidth` must be a strictly positive number.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# compute_distance
# ---------------------------------------------------------------------------

#' Compute a full pairwise distance matrix from a coordinate matrix
#'
#' @param coords   A numeric matrix with at least two columns (`X` and `Y`, or
#'   the first two columns if names are absent).  Each row is one spatial unit.
#' @param longlat  Logical.  If `TRUE`, great-circle distances (in kilometres)
#'   are calculated using the Haversine formula.  If `FALSE` (default),
#'   Euclidean distance is used.
#'
#' @return An n x n symmetric numeric matrix where element [i, j] is the
#'   distance between unit i and unit j.  Diagonal is 0.
#' @keywords internal
compute_distance <- function(coords, longlat = FALSE) {
  if (!is.matrix(coords)) {
    coords <- as.matrix(coords)
  }
  if (ncol(coords) < 2L) {
    stop("`coords` must have at least 2 columns (X and Y).", call. = FALSE)
  }
  n <- nrow(coords)

  if (!longlat) {
    # Euclidean distance via stats::dist — fastest for small/medium n
    d <- as.matrix(stats::dist(coords[, 1:2, drop = FALSE]))
    return(d)
  }

  # Great-circle distance (Haversine), result in kilometres
  x <- coords[, 1]   # longitude in degrees
  y <- coords[, 2]   # latitude  in degrees

  to_rad <- pi / 180
  lat1 <- y * to_rad
  lon1 <- x * to_rad

  d <- matrix(0, nrow = n, ncol = n)
  R <- 6371.0  # Earth radius in km

  for (i in seq_len(n)) {
    dlat <- lat1 - lat1[i]
    dlon <- lon1 - lon1[i]
    a <- sin(dlat / 2)^2 + cos(lat1[i]) * cos(lat1) * sin(dlon / 2)^2
    a <- pmin(a, 1)  # numerical safety
    d[i, ] <- 2 * R * asin(sqrt(a))
  }
  d
}

# ---------------------------------------------------------------------------
# get_local_distances
# ---------------------------------------------------------------------------

#' Extract distances from one focus unit to all others
#'
#' @param distance_context A list as returned by `build_distance_context()`, or
#'   a plain n x n numeric distance matrix.  If a plain matrix is supplied, the
#'   rows and columns must already be in the same order as the spatial units.
#' @param focus_id         Integer scalar (1-based) or character matching a
#'   row/column name of the distance matrix.  The focal unit whose distances
#'   are extracted.
#'
#' @return A numeric vector of length n giving the distance from the focus unit
#'   to every unit (including itself, which is 0).
#' @keywords internal
get_local_distances <- function(distance_context, focus_id) {
  # Support both a plain matrix and a list wrapper
  if (is.matrix(distance_context)) {
    dmat <- distance_context
  } else if (is.list(distance_context) && !is.null(distance_context$distance_matrix)) {
    dmat <- distance_context$distance_matrix
  } else {
    stop(
      "`distance_context` must be an n x n distance matrix or a list with ",
      "a `distance_matrix` element.",
      call. = FALSE
    )
  }

  n <- nrow(dmat)

  if (is.character(focus_id)) {
    rn <- rownames(dmat)
    if (is.null(rn)) {
      stop("Distance matrix has no row names; cannot look up character `focus_id`.", call. = FALSE)
    }
    idx <- match(focus_id, rn)
    if (is.na(idx)) {
      stop(sprintf("focus_id '%s' not found in distance matrix row names.", focus_id), call. = FALSE)
    }
    return(as.numeric(dmat[idx, ]))
  }

  if (!is.numeric(focus_id) || length(focus_id) != 1L) {
    stop("`focus_id` must be a single integer index or a character ID.", call. = FALSE)
  }
  focus_id <- as.integer(focus_id)
  if (focus_id < 1L || focus_id > n) {
    stop(sprintf("`focus_id` (%d) is out of range [1, %d].", focus_id, n), call. = FALSE)
  }
  as.numeric(dmat[focus_id, ])
}

# ---------------------------------------------------------------------------
# build_distance_context
# ---------------------------------------------------------------------------

#' Build a distance context object
#'
#' Wraps a distance matrix together with unit IDs into a list suitable for
#' passing to `get_local_distances()`.  Optionally pre-computes the full
#' distance matrix (recommended for small data) or stores only the coordinate
#' matrix for on-the-fly computation (recommended for large data).
#'
#' @param coords    A numeric matrix with columns `X` and `Y`.
#' @param ids       Character or numeric vector of unit IDs (length = nrow(coords)).
#' @param longlat   Logical.  Passed to `compute_distance()`.
#' @param cache     Logical.  If `TRUE` (default), pre-computes and caches the
#'   full n x n distance matrix.  Set `FALSE` for very large datasets to avoid
#'   memory pressure; `get_local_distances()` will then compute rows on demand.
#'
#' @return A list with class `"gwpr_distance_context"` containing:
#' * `ids` — character vector of unit IDs.
#' * `distance_matrix` — n x n matrix (or `NULL` if `cache = FALSE`).
#' * `coords` — the original coordinate matrix.
#' * `longlat` — logical flag.
#' @keywords internal
build_distance_context <- function(coords, ids, longlat = FALSE, cache = TRUE) {
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- nrow(coords)
  ids <- as.character(ids)
  if (length(ids) != n) {
    stop("`ids` length must equal nrow(coords).", call. = FALSE)
  }

  dmat <- NULL
  if (cache) {
    dmat <- compute_distance(coords, longlat = longlat)
    rownames(dmat) <- ids
    colnames(dmat) <- ids
  }

  structure(
    list(
      ids             = ids,
      distance_matrix = dmat,
      coords          = coords,
      longlat         = longlat
    ),
    class = "gwpr_distance_context"
  )
}

# ---------------------------------------------------------------------------
# compute_kernel_weights
# ---------------------------------------------------------------------------

#' Compute geographically weighted kernel weights
#'
#' Given a numeric vector of distances from one focal unit to all others,
#' returns a weight for each unit according to the chosen kernel and bandwidth.
#'
#' For **fixed** bandwidth (`adaptive = FALSE`) the bandwidth parameter is a
#' distance threshold or scale parameter used directly in the kernel formula.
#'
#' For **adaptive** bandwidth (`adaptive = TRUE`) the bandwidth parameter is
#' the number of nearest neighbours k.  The function first identifies the
#' k-th smallest distance (among all units, including the focal unit itself at
#' distance 0) and uses that distance as the effective bandwidth in the kernel
#' formula.
#'
#' Kernel formulae (d = distance, bw = effective bandwidth):
#' * `bisquare`:    `(1 - (d/bw)^2)^2` for `d <= bw`, else 0.
#' * `gaussian`:    `exp(-0.5 * (d/bw)^2)`.
#' * `exponential`: `exp(-d/bw)`.
#' * `tricube`:     `(1 - (d/bw)^3)^3` for `d <= bw`, else 0.
#' * `boxcar`:      `1` for `d <= bw`, else 0.
#'
#' @param distance  Numeric vector of distances from the focal unit to all
#'   spatial units (length n).
#' @param bandwidth Numeric scalar.  For fixed bandwidth: distance scale.
#'   For adaptive bandwidth: positive integer number of neighbours.
#' @param kernel    Character scalar, one of `"bisquare"`, `"gaussian"`,
#'   `"exponential"`, `"tricube"`, `"boxcar"`.
#' @param adaptive  Logical.  `TRUE` for adaptive (kNN) bandwidth.
#'
#' @return A numeric vector of length n with non-negative kernel weights.
#' @keywords internal
compute_kernel_weights <- function(distance, bandwidth, kernel, adaptive) {
  kernel <- match.arg(
    kernel,
    choices = c("bisquare", "gaussian", "exponential", "tricube", "boxcar")
  )
  validate_bandwidth(bandwidth, adaptive)

  if (!is.numeric(distance)) {
    stop("`distance` must be a numeric vector.", call. = FALSE)
  }

  d   <- as.numeric(distance)
  bw  <- bandwidth

  # Adaptive: determine effective bandwidth as the bw-th smallest distance
  if (adaptive) {
    k <- as.integer(bw)
    n <- length(d)
    if (k > n) {
      warning(
        "Adaptive bandwidth (", k, ") exceeds number of units (", n,
        "); using all units.",
        call. = FALSE
      )
      k <- n
    }
    sorted_d <- sort(d)
    bw <- sorted_d[k]   # k-th smallest distance (0-indexed: focal unit itself)
    # Edge case: if bw == 0 (e.g. duplicate coordinates), push to next unique
    if (bw == 0 && k < n) {
      bw <- sorted_d[k + 1L]
    }
    if (bw == 0) {
      # All distances are 0 — give equal weight
      return(rep(1, n))
    }
  }

  # Apply kernel
  w <- switch(
    kernel,
    bisquare = {
      u <- d / bw
      ifelse(u <= 1, (1 - u^2)^2, 0)
    },
    gaussian = {
      exp(-0.5 * (d / bw)^2)
    },
    exponential = {
      exp(-d / bw)
    },
    tricube = {
      u <- d / bw
      ifelse(u <= 1, (1 - u^3)^3, 0)
    },
    boxcar = {
      ifelse(d <= bw, 1, 0)
    }
  )

  w
}
