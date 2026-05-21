## Resubmission

This is a resubmission of GWPR.light. The package was archived on 2024-02-10.
This is a complete rewrite as version 1.0.0 with a new, modern API.

### Changes in this version

- Complete API rewrite: new `sf`-first public interface with four main functions:
  `gwpr()`, `select_bandwidth()`, `fit_gwpr()`, `diagnose_gwpr()`
- Removed legacy `sp`/`SpatialPolygonsDataFrame` interface
- Added Gaussian and binomial family support, three bandwidth search strategies,
  five kernel functions, and optional parallel execution via the `future` framework
- Internal helper functions are no longer exported from the NAMESPACE

---

## Test environments

- macOS 26.3 (aarch64-apple-darwin20), R 4.3.2

## R CMD check results

```
Status: 0 ERRORs, 0 WARNINGs, 2 NOTEs (on CRAN servers)
```

**Local check shows 1 ERROR + 1 WARNING** because `pdflatex` is not installed
on the local test machine. The errors are:

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
* checking PDF version of manual without index ... ERROR
pdflatex is not available
```

These are local environment issues only and will not occur on CRAN servers.

**NOTEs (expected on CRAN):**

1. `Maintainer: 'Chao Li <chaoli0394@gmail.com>'` — informational
2. `New submission / Package was archived on CRAN` — expected for re-submission
3. `unable to verify current time` — network issue in check environment
