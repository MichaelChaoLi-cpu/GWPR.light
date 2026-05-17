
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GWPR.light

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/MichaelChaoLi-cpu/GWPR.light/workflows/R-CMD-check/badge.svg)](https://github.com/MichaelChaoLi-cpu/GWPR.light/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/GWPR.light)](https://CRAN.R-project.org/package=GWPR.light)
[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

A modern, `sf`-first implementation of **Geographically Weighted Panel
Regression (GWPR)** for spatial panel data. Version 1.0.0 provides a
clean public API, three bandwidth search strategies (grid, SGD, random),
Gaussian and binomial model families, and optional parallel execution
via the `future` framework.

## Authors

Chao Li <chaoli0394@gmail.com>  
Shunsuke Managi <managi@doc.kyushu-u.ac.jp>

## Installation

Install the released version from [CRAN](https://CRAN.R-project.org):

``` r
install.packages("GWPR.light")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("MichaelChaoLi-cpu/GWPR.light")
```

## Quick start (1.0.0 API)

The four public functions form the complete interface:

| Function | Purpose |
|---|---|
| `gwpr()` | Full pipeline: bandwidth search + fit + diagnostics |
| `select_bandwidth()` | Standalone bandwidth optimisation |
| `fit_gwpr()` | Fit with a known bandwidth |
| `diagnose_gwpr()` | Diagnostic tests on a fitted model |

``` r
library(GWPR.light)
library(sf)

# --- Simulate a small spatial panel ---
set.seed(1)
pts <- sf::st_as_sf(
  data.frame(id = 1:6, X = c(0,1,2,0,1,2), Y = c(0,0,0,1,1,1)),
  coords = c("X", "Y"), crs = NA_integer_
)
dat <- data.frame(
  id   = rep(1:6, each = 4),
  time = rep(1:4, 6),
  x1   = rnorm(24),
  x2   = rnorm(24)
)
dat$y <- 1.5 * dat$x1 - 0.8 * dat$x2 + rnorm(24, sd = 0.3)

# --- Fit with a known bandwidth ---
fit <- fit_gwpr(
  y ~ x1 + x2, data = dat, spatial = pts,
  id = "id", time = "time",
  bandwidth = 2, model = "pooling", workers = 1
)
print(fit)
```

### Automatic bandwidth search

``` r
bw <- select_bandwidth(
  y ~ x1 + x2, data = dat, spatial = pts,
  id = "id", time = "time",
  method  = "grid",
  control = list(lower = 0.5, upper = 3, step = 0.5),
  workers = 1
)
cat("Optimal bandwidth:", bw$best_bandwidth, "\n")
```

### Diagnostics

``` r
diag_result <- diagnose_gwpr(fit, diagnostics = c("f_test", "hausman"))
print(diag_result)
```

## Vignettes

- `vignette("gwpr-introduction")` — Getting started with the 1.0.0 API
- `vignette("introduction_of_GWPR")` — Legacy 0.2.x API reference

## References

- Fotheringham, A., Brunsdon, C., Charlton, M. (2002). *Geographically
  Weighted Regression*. John Wiley & Sons. ISBN: 978-0-470-85525-6.
- Beenstock, M., Felsenstein, D. (2019). *The Econometric Analysis of
  Non-Stationary Spatial Panel Data*. Springer.
  ISBN: 978-3-030-03614-0.

## Bug reports

<https://github.com/MichaelChaoLi-cpu/GWPR.light/issues>
