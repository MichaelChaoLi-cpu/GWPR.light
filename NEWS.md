# GWPR.light 1.0.0

## Breaking changes

- **New public API**: the primary interface is now four `snake_case` functions:
  `gwpr()`, `select_bandwidth()`, `fit_gwpr()`, `diagnose_gwpr()`.
- **`sf` first**: spatial input must be an `sf` object. The legacy `sp` /
  `SpatialPolygonsDataFrame` interface of the 0.x API is not supported by the
  new functions.
- **Old functions deprecated**: `GWPR()`, `bw.GWPR()`, `GWPR.moran.test()`,
  `GWPR.pFtest()`, `GWPR.phtest()`, and `GWPR.plmtest()` remain in the package
  for backwards compatibility but will be removed in a future version.

## New features

- **`gwpr()`**: one-call pipeline that validates inputs, optionally searches for
  the optimal bandwidth, fits the model, and runs diagnostics.
- **`select_bandwidth()`**: standalone bandwidth optimisation with three search
  strategies:
  - `"grid"` — exhaustive grid search over a user-defined range.
  - `"sgd"` — stochastic gradient descent (fast, scales to large grids).
  - `"random"` — random sampling over a candidate set.
- **`fit_gwpr()`**: fit a GWPR model with a known bandwidth. Supports
  `family = "gaussian"` (linear) and `family = "binomial"` (panel logistic).
- **`diagnose_gwpr()`**: collects results from Moran's I (`"moran"`), local
  F-test (`"f_test"`), Hausman test (`"hausman"`), and LM test (`"lm_test"`).
- **Five kernel functions**: `"bisquare"` (default), `"gaussian"`,
  `"exponential"`, `"tricube"`, `"boxcar"`.
- **Adaptive bandwidth** (`adaptive = TRUE`) using k-nearest-neighbour approach.
- **Optional parallelism** via `workers` argument (requires `future` and
  `future.apply`; default `workers = 1` is always safe).
- **Memory estimation**: `gwpr()` warns when the estimated memory footprint is
  high before starting long computations.
- **Reproducibility**: `seed` argument propagates to all stochastic steps.

## Dependency changes

- Added: `fixest`, `glmmTMB` (for binomial panel logistic models).
- `future` and `future.apply` moved to **Suggests** (only needed for
  `workers > 1`).
- `R (>= 3.5.0)` minimum version (was `R (>= 2.10)`).

---

# GWPR.light 0.2.1

- Fixed compatibility with development version of tidyselect.

# GWPR.light 0.2.0

- Revised examples to reduce run time.
- Removed manual documents for internal helper functions.

# GWPR.light 0.1.0

- Initial CRAN release.
