#' @title California (sf)
#'
#' @description The counties' boundary in California
#'
#' @usage data(California)
#'
#' @format An \code{sf} object with 58 rows (one per county) and two columns:
#' \describe{
#' \item{GEOID}{a numeric vector, fips IDs of the counties}
#' \item{geometry}{sfc_MULTIPOLYGON, county boundary polygons (CRS: NAD83 longlat)}
#' }
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi <managi.s@gmail.com>
#'
#' @examples
#' data(California)
#' class(California)
"California"
