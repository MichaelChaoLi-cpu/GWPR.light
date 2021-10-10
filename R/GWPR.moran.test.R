#' Moran's I Test for Panel Regression
#'
#' @description Moran's I test for spatial autocorrelation in residuals from
#'              an estimated panel linear model (plm).
#'
#' @usage GWPR.moran.test(plm_model, SDF, bw, adaptive = FALSE, p = 2,
#'                        kernel = "bisquare", longlat = FALSE, alternative = "greater")
#'
#' @param plm_model     An object of class inheriting from "plm", see plm
#' @param SDF           Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param bw            The optimal bandwidth, either adaptive or fixed distance
#' @param adaptive      If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param p             The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param kernel        bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param longlat       If TRUE, great circle distances will be calculated
#' @param alternative   A character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom sp merge coordinates
#' @importFrom plm pdim index
#' @importFrom methods is
#' @importFrom stats pnorm
#'
#' @return A list of result:
#' \describe{
#' \item{statistic}{the value of the standard deviate of Moran's I.}
#' \item{p.value}{the p-value of the test.}
#' \item{Estimated.I}{the value of the observed Moran's I.}
#' \item{Excepted.I}{the value of the expectation of Moran's I.}
#' \item{V2}{the value of the variance of Moran's I.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' }
#' @export
#'
#' @note: Current version of panel Moran's I test can only chech the balanced panel data.
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi <managi.s@gmail.com>
#'
#' @references Beenstock, M., Felsenstein, D., 2019. The econometric analysis of non-stationary spatial panel data. Springer.
#'
#' @examples
#' data(TransAirPolCalif)
#' data(California)
#' formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
#'    Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
#'    Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
#'    Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
#'    Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
#'    pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
#'
#' pdata <- plm::pdata.frame(TransAirPolCalif, index = c("GEOID", "year"))
#' moran.plm.model <- plm::plm(formula = formula.GWPR, data = pdata, model = "within")
#' summary(moran.plm.model)
#'
#' #precomputed bandwidth
#' bw.AIC.Fix <- 2.010529
#'
#' # moran's I test
#' GWPR.moran.test(moran.plm.model, SDF = California, bw = bw.AIC.Fix, kernel = "bisquare",
#'                  adaptive = FALSE, p = 2, longlat = FALSE, alternative = "greater")
GWPR.moran.test <- function(plm_model, SDF, bw, adaptive = FALSE, p = 2, kernel = "bisquare",
                             longlat = FALSE, alternative = "greater")
{
  if(class(plm_model)[1] != "plm")
  {
    stop("This test only accepts the object \"plm\".")
  }
  if(!plm::pdim(plm_model$model)$balanced)
  {
    stop("Current version only accepts balanced panel regression.")
  }
  if(!methods::is(SDF, "Spatial"))
  {
    stop("SDF should be the spatial data frame based on \"sp\".")
  }
  if(!(colnames(plm::index(plm_model$model))[1] %in% colnames(SDF@data)))
  {
    stop("Indexes in plm and SDP are not consistent.")
  }

  plm.resid <- as.data.frame(as.matrix(plm_model$residuals))
  n <- nrow(plm.resid)
  Ti <- ncol(plm.resid)
  N <- Ti*n
  id <- as.data.frame(as.numeric(row.names(plm.resid)))
  colnames(id) <- "id"
  SDF <- sp::merge(SDF, id, by.x = dplyr::all_of((colnames(plm::index(plm_model$model))[1])), by.y = "id" )
  row.id <- as.vector(as.matrix(dplyr::select(SDF@data, dplyr::all_of((colnames(plm::index(plm_model$model))[1])))))
  dp.locat <- as.matrix(sp::coordinates(SDF))
  coord <- cbind(as.data.frame(dp.locat), dplyr::select(SDF@data, dplyr::all_of((colnames(plm::index(plm_model$model))[1]))) )
  colnames(coord) <- c("X", "Y", "id")
  coord <- dplyr::arrange(coord, "id")
  dp.locat <- as.matrix(dplyr::select(coord, "X", "Y"))
  dMat <- GWmodel::gw.dist(dp.locat = dp.locat, rp.locat = dp.locat,
                           focus = 0, p = p, longlat=longlat)
  if(adaptive)
  {
    bw <- bw * Ti
  } # if adaptive bandwidth, the input number is the numbers of individuals rather than records.
  weight <- GWmodel::gw.weight(dMat, bw=bw, kernel=kernel, adaptive=adaptive)
  diag(weight) <- 0
  weight <- weight/rowSums(weight)

  I.vector <- c()
  loop_time <- 1
  while(loop_time < (ncol(plm.resid)+1) )
  {
    plm.resid.single <- plm.resid[,loop_time]
    sum.weight.residuals <- 0
    i <- 1
    j <- 1
    while(i < (length(plm.resid.single)+1))
    {
      while(j < (length(plm.resid.single)+1))
      {
        sum.weight.residuals <- sum.weight.residuals + weight[i,j]*plm.resid.single[i]*plm.resid.single[j]
        j <- j + 1
      }
      i <- i + 1
    }
    sum.residuals <- sum(plm.resid.single^2)
    I <- sum.weight.residuals/sum.residuals
    I.vector <- append(I.vector, I)
    loop_time <- loop_time + 1
  }
  I.mean <- mean(I.vector)

  # V2 of average I
  V2.upper <-  n*sum(weight^2) + 3 * (sum(weight))^2 - n* sum((colSums(weight))^2)
  V2.lower <- Ti * (n^2 - 1) * (sum(weight))^2
  V2 <- V2.upper/V2.lower

  # Expected value of I
  E <- -1/(n - 1)

  ZI <- (I.mean - E) / sqrt(V2)

  if (alternative == "two.sided")
  {
    PrI <- 2 * stats::pnorm(abs(ZI), lower.tail=FALSE)
  }
  else
  {
    if (alternative == "greater")
    {
      PrI <- stats::pnorm(ZI, lower.tail=FALSE)
    }
    else
    {
      PrI <- stats::pnorm(ZI)
    }
  }

  res <- list(statistic = ZI, p.value=PrI, Estimated.I = I.mean, Excepted.I = E, V2 = V2,
              alternative=alternative)
  return(res)
}
