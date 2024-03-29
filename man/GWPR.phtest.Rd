% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWPR.phtest.R
\name{GWPR.phtest}
\alias{GWPR.phtest}
\title{Locally Hausman Test Based on GWPR}
\usage{
GWPR.phtest(formula, data, index, SDF, bw = NULL, adaptive = FALSE,
                   p = 2, effect = "individual", random.method = "swar",
                   kernel = "bisquare", longlat = FALSE)
}
\arguments{
\item{formula}{The regression formula: : Y ~ X1 + ... + Xk}

\item{data}{A data.frame for the Panel data.}

\item{index}{A vector of the two indexes: (c("ID", "Time")).}

\item{SDF}{Spatial*DataFrame on which is based the data, with the "ID" in the index.}

\item{bw}{The optimal bandwidth, either adaptive or fixed distance.}

\item{adaptive}{If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.}

\item{p}{The power of the Minkowski distance, default is 2, i.e. the Euclidean distance}

\item{effect}{The effects introduced in the fixed effects model, one of "individual" (default) , "time", "twoways"}

\item{random.method}{Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"}

\item{kernel}{bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise}

\item{longlat}{If TRUE, great circle distances will be calculated}
}
\value{
A list of result:
\describe{
\item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
\item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points, test value, p value, df}
}
}
\description{
Locally Hausman Test Based on GWPR
}
\note{
If the random method is "swar", to perform this test, bandwidth selection must guarantee that enough individuals in the subsamples. Using bw.GWPR function can avoid mistake.
}
\examples{
\donttest{
data(TransAirPolCalif)
data(California)
formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
   Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
   Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
   Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
   Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
   pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax

#precomputed bandwidth
bw.AIC.Fix <- 7.508404

GWPR.phtest.resu.F <- GWPR.phtest(formula = formula.GWPR, data = TransAirPolCalif,
                                  index = c("GEOID", "year"), SDF = California,
                                  bw = bw.AIC.Fix, adaptive = FALSE, p = 2,
                                  effect = "individual", kernel = "bisquare",
                                  longlat = FALSE)
library(tmap)
tm_shape(GWPR.phtest.resu.F$SDF) +
     tm_polygons(col = "p.value", breaks = c(0, 0.05, 1))
}
}
\author{
Chao Li <chaoli0394@gmail.com> Shunsuke Managi
}
