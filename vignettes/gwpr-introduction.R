## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----minimal-linear-----------------------------------------------------------
library(GWPR.light)
library(sf)

set.seed(42)

# Simulate a tiny spatial panel: 6 units, 4 time periods
n_units <- 6
n_time  <- 4

pts <- sf::st_as_sf(
  data.frame(
    id = 1:n_units,
    X  = c(0, 1, 2, 0, 1, 2),
    Y  = c(0, 0, 0, 1, 1, 1)
  ),
  coords = c("X", "Y"),
  crs    = NA_integer_
)

dat <- data.frame(
  id   = rep(1:n_units, each = n_time),
  time = rep(1:n_time,  n_units),
  x1   = rnorm(n_units * n_time),
  x2   = rnorm(n_units * n_time)
)
dat$y <- 1.5 * dat$x1 - 0.8 * dat$x2 + rnorm(n_units * n_time, sd = 0.3)

# Fit with a known bandwidth (skip automatic search for speed)
fit <- fit_gwpr(
  formula   = y ~ x1 + x2,
  data      = dat,
  spatial   = pts,
  id        = "id",
  time      = "time",
  bandwidth = 2,
  family    = "gaussian",
  model     = "pooling",
  workers   = 1
)

print(fit)

## ----results------------------------------------------------------------------
# Overall goodness-of-fit metrics
str(fit$metrics)

# Per-unit spatial coefficients (one row per spatial unit)
if (!is.null(fit$spatial_results)) {
  head(fit$spatial_results)
}

## ----bandwidth-search---------------------------------------------------------
bw <- select_bandwidth(
  formula  = y ~ x1 + x2,
  data     = dat,
  spatial  = pts,
  id       = "id",
  time     = "time",
  family   = "gaussian",
  model    = "pooling",
  method   = "grid",
  control  = list(lower = 0.5, upper = 3, step = 0.5),
  workers  = 1
)

print(bw)
cat("Best bandwidth:", bw$best_bandwidth, "\n")

## ----full-pipeline------------------------------------------------------------
# Use the best bandwidth found above to avoid re-running search
full_fit <- gwpr(
  formula     = y ~ x1 + x2,
  data        = dat,
  spatial     = pts,
  id          = "id",
  time        = "time",
  bandwidth   = bw$best_bandwidth,
  family      = "gaussian",
  model       = "pooling",
  diagnostics = FALSE,   # skip diagnostics for speed
  workers     = 1
)

print(full_fit)

## ----diagnostics--------------------------------------------------------------
diag_result <- diagnose_gwpr(
  full_fit,
  diagnostics = c("f_test", "hausman", "lm_test")
)

print(diag_result)

## ----long-examples, eval=FALSE------------------------------------------------
#  # Automatic SGD bandwidth search + fit (may take several seconds)
#  fit_auto <- gwpr(
#    formula          = y ~ x1 + x2,
#    data             = dat,
#    spatial          = pts,
#    id               = "id",
#    time             = "time",
#    bandwidth_method = "sgd",
#    bandwidth_control = list(n_iter = 20, step_size = 0.1),
#    workers          = 1,
#    seed             = 123
#  )
#  
#  # Binomial GWPR
#  dat$y_bin <- as.integer(dat$y > 0)
#  fit_logit <- fit_gwpr(
#    formula   = y_bin ~ x1 + x2,
#    data      = dat,
#    spatial   = pts,
#    id        = "id",
#    time      = "time",
#    bandwidth = 2,
#    family    = "binomial",
#    model     = "pooling",
#    workers   = 1
#  )

