#' @title Panel Dataset for Testing GWPR
#'
#' @description Panel dataset to estimate the relationship between county-level PM2.5 concentration and on-road transporation in California.
#'
#' @usage data(TransAirPolCalif)
#'
#' @format A \code{data.frame} with 23 variables, and 928 observations, which are:
#' \describe{
#' \item{GEOID}{a numeric vector, fips IDs of the counties}
#' \item{year}{a numeric vector, year}
#' \item{pm25}{a numeric vector, annually average PM2.5 concentration in the counties}
#' \item{co2_mean}{a numeric vector, geographically average CO2 emission from on-road transportation in each year, million tons/km2}
#' \item{Developed_Open_Space_perc}{a numeric vector, percentage of developed open space of total area in each county}
#' \item{Developed_Low_Intensity_perc}{a numeric vector, percentage of low-intensity developed area of total area in each county}
#' \item{Developed_Medium_Intensity_perc}{a numeric vector, percentage of medium-intensity developed area of total area in each county}
#' \item{Developed_High_Intensity_perc}{a numeric vector, percentage of high-intensity develope area of total area in each county}
#' \item{Open_Water_perc}{a numeric vector, percentage of open water of total area in each county}
#' \item{Woody_Wetlands_perc}{a numeric vector, percentage of woody wetland of total area in each county}
#' \item{Emergent_Herbaceous_Wetlands_perc}{a numeric vector, percentage of emergent herbaceous wetland of total area in each county}
#' \item{Deciduous_Forest_perc}{a numeric vector, percentage of deciduous forest of total area in each county}
#' \item{Evergreen_Forest_perc}{a numeric vector, percentage of evergreen forest of total area in each county}
#' \item{Mixed_Forest_perc}{a numeric vector, percentage of mixed forest of total area in each county}
#' \item{Shrub_perc}{a numeric vector, percentage of shrub of total area in each county}
#' \item{Grassland_perc}{a numeric vector, percentage of grassland of total area in each county}
#' \item{Pasture_perc}{a numeric vector, percentage of pasture of total area in each county}
#' \item{Cultivated_Crops_perc}{a numeric vector, percentage of cultivated crops of total area in each county}
#' \item{pop_density}{a numeric vector, average population density in each county}
#' \item{summer_tmmx}{a numeric vector, average temperature in summer}
#' \item{winter_tmmx}{a numeric vector, average temperature in winter}
#' \item{summer_rmax}{a numeric vector, average humidity in summer}
#' \item{winter_rmax}{a numeric vector, average humidity in winter}
#' }
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi <managi.s@gmail.com>
#'
#' @examples \dontrun{
#' data(TransAirPolCalif)
#' head(TransAirPolCalif)
#' }
"TransAirPolCalif"


