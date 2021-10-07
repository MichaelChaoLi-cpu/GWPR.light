## ---- include = FALSE---------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----example_GWmodel, message=FALSE-------------------------------------------
library(GWPR.light)
library(sp)
library(tmap)
data("California")
data("TransAirPolCalif")
formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
   Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
   Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
   Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
   Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
   pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
TransAirPolCalif <- dplyr::filter(TransAirPolCalif, year == 2001)
lm.2001 <- lm(formula = formula.GWPR, TransAirPolCalif)

## -----------------------------------------------------------------------------
summary(lm.2001)
resid.lm <- cbind(lm.2001$residuals, TransAirPolCalif$GEOID)
colnames(resid.lm) <- c("resid", "GEOID")
head(resid.lm)

## ----cross_sectional_residual, fig.align='center', out.height= '80%', out.width='80%'----
to_show <- sp::merge(California, resid.lm, by = "GEOID")
tm_shape(to_show) + tm_polygons(col = 'resid')

## ----setup, message=FALSE-----------------------------------------------------
library(GWPR.light)

## ---- eval = FALSE------------------------------------------------------------
#  bw.GWPR(formula = formula, data = data, index = index, SDF = SDF,
#          adaptive = F, p = 2,bigdata = F, upperratio = 0.25,
#          effect = "individual", model = c("pooling", "within", "random"),
#          random.method = "swar", approach = c("CV","AIC"), kernel = "bisquare",
#          longlat = F, doParallel = F, cluster.number = 2, human.set.range = F,
#          h.upper = NULL, h.lower = NULL)

## ---- out.height="60%", out.width="60%", fig.align="center"-------------------
h = 1:10 # just an example, not truth.
SS = h^2
plot(h, SS, col = "red")

## ---- fig.align='center', out.height= '80%', out.width='80%'------------------
tm_shape(to_show) +
  tm_polygons(col = 'resid')

## -----------------------------------------------------------------------------
data(TransAirPolCalif)
pdata <- plm::pdata.frame(TransAirPolCalif, index = c("GEOID", "year"))
moran.plm.model <- plm::plm(formula = formula.GWPR, data = pdata, model = "within")
summary(moran.plm.model)

## ---- eval = FALSE------------------------------------------------------------
#  bw.AIC.F <- bw.GWPR(formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
#                      SDF = California, adaptive = F, p = 2, bigdata = F, effect = "individual",
#                      model = "within", approach = "AIC", kernel = "bisquare", longlat = F,
#                      doParallel = T, cluster.number = 2)

## -----------------------------------------------------------------------------
bw.AIC.F <- 2.010529 #precomputed results from #150:155
GWPR.moran.test(moran.plm.model, SDF = California, bw = bw.AIC.F, kernel = "bisquare",
                 adaptive = F, p = 2, longlat=F, alternative = "greater")

## ----eg.GWPR.pFtest,  fig.align='center', out.height= '80%', out.width='80%'----
GWPR.pFtest.resu.F <- GWPR.pFtest(formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
                                  SDF = California, bw = bw.AIC.F, adaptive = F, p = 2, effect = "individual",
                                  kernel = "bisquare", longlat = F)
tm_shape(GWPR.pFtest.resu.F$SDF) +
     tm_polygons(col = "p.value", breaks = c(0, 0.05, 1))

## ----eg.GWPR.plmtest, fig.align='center', out.height= '80%', out.width='80%'----
GWPR.plmtest.resu.F <- GWPR.plmtest(formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
                                    SDF = California, bw = bw.AIC.F, adaptive = F, p = 2,
                                    kernel = "bisquare", longlat = F)
tm_shape(GWPR.plmtest.resu.F$SDF) +
     tm_polygons(col = "p.value", breaks = c(0, 0.05, 1))

## ----eg.GWPR.phtest, fig.align='center', out.height= '80%', out.width='80%'----
GWPR.phtest.resu.F <- GWPR.phtest(formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
                                  SDF = California, bw = bw.AIC.F, adaptive = F, p = 2, effect = "individual",
                                  kernel = "bisquare", longlat = F, random.method = "amemiya")
tm_shape(GWPR.phtest.resu.F$SDF) +
     tm_polygons(col = "p.value", breaks = c(0, 0.05, 1))

## ---- eval=FALSE--------------------------------------------------------------
#  GWPR(formula = formula, data = data, index = index, SDF = SDF, bw = NULL, adaptive = F, p = 2,
#              effect = "individual", model = c("pooling", "within", "random"),
#              random.method = "swar", kernel = "bisquare", longlat = F)

## ----eg.GWPR, fig.align='center', out.height= '80%', out.width='80%'----------
result.F.AIC <- GWPR(bw = bw.AIC.F, formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
                     SDF = California, adaptive = F, p = 2, effect = "individual", model = "within",
                     kernel = "bisquare", longlat = F)
summary(result.F.AIC$SDF$Local_R2)
tm_shape(result.F.AIC$SDF) +
  tm_polygons(col = "Local_R2", pal = "Reds",auto.palette.mapping = F,
              style = 'cont')

