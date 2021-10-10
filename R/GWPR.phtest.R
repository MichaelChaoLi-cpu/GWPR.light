#' Locally Hausman Test based on GWPR
#'
#' @usage GWPR.phtest(formula, data, index, SDF, bw = NULL, adaptive = FALSE,
#'                    p = 2, effect = "individual", random.method = "swar",
#'                    kernel = "bisquare", longlat = FALSE)
#'
#' @param formula        The regression formula: : Y ~ X1 + ... + Xk
#' @param data           A data.frame for the Panel data.
#' @param index          A vector for the indexes : (c("ID", "Time")).
#' @param SDF            Spatial*DataFrame on which is based the data, with the "ID" in the index.
#' @param bw             The optimal bandwidth, either adaptive or fixed distance.
#' @param adaptive       If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param p              The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param effect         The effects introduced in the fixed effects model, one of "individual" (default) , "time", "twoways"
#' @param random.method  Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param kernel         bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param longlat        If TRUE, great circle distances will be calculated
#'
#' @import dplyr
#' @importFrom data.table setDT
#' @importFrom sp coordinates bbox
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points, test value, p value, df}
#' }
#' @export
#'
#' @note If the random method is "swar", to perform this test, bandwidth selection must guarantee that enough individuals in the subsamples. Using bw.GWPR function can avoid mistake.
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi <managi.s@gmail.com>
#'
#' @examples
#' \donttest{
#' data(TransAirPolCalif)
#' data(California)
#' formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
#'    Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
#'    Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
#'    Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
#'    Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
#'    pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
#'
#' #precomputed bandwidth
#' bw.AIC.Fix <- 7.508404
#'
#' GWPR.phtest.resu.F <- GWPR.phtest(formula = formula.GWPR, data = TransAirPolCalif,
#'                                   index = c("GEOID", "year"), SDF = California,
#'                                   bw = bw.AIC.Fix, adaptive = FALSE, p = 2,
#'                                   effect = "individual", kernel = "bisquare",
#'                                   longlat = FALSE)
#' library(tmap)
#' tm_shape(GWPR.phtest.resu.F$SDF) +
#'      tm_polygons(col = "p.value", breaks = c(0, 0.05, 1))
#' }
GWPR.phtest <- function(formula, data, index, SDF, bw = NULL, adaptive = FALSE, p = 2, effect = "individual",
                        random.method = "swar", kernel = "bisquare", longlat = FALSE)
{
  if(length(index) != 2)
  {
    stop("The \"index\" have included \"ID\" or/and \"time\" index.")
  }
  if(!((index[1] %in% colnames(data)) & (index[2] %in% colnames(data))))
  {
    stop("The data.frame(data) does not have the index columns.")
  }
  if(!(index[1] %in% colnames(SDF@data)))
  {
    stop("The SDF does not have the \"ID\" columns.")
  }
  if(is.null(bw))
  {
    stop("The bw must be set.")
  }

  if(random.method == "swar")
  {
    message("***************************************** Caution! **********************************************\n",
            "Your random method is", random.method,". According to the \"plm\" package requirement, data sizes of each subsample should be larger\n",
            "than the number of estimated parameters.\nOtherwise, errors. Therefore, we strongly recommend using bw.GWPR to get bandwidth.\n")
  }
  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  data <- dplyr::select(data, index, varibale_name_in_equation)
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"
  model <- "within"

  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  ID <- dplyr::select(data, index[1])
  .N <- 0
  ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  if(model == "within")
  {
    data <- drop_ID_with_single_observation(data, ID_num)
    ID <- dplyr::select(data, index[1])
    ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  }


  # Judge the data size of calculation
  if (nrow(ID_num) > 1000)
  {
    message("Dear my friend, thanks for your patience!. We pass the bandwidth\n",
            "selection part. Now, regression! This should be faster. Thanks.\n",
            "................................................................\n")
    huge_data_size <- TRUE
  }
  else
  {
    huge_data_size <- FALSE
  }

  # Panel SDF preparation
  SDF@data <- dplyr::select(SDF@data, dplyr::all_of(raw_id))
  colnames(SDF@data)[1] <- "id"
  dp.locat <- sp::coordinates(SDF)
  coord <- cbind(as.data.frame(dp.locat), SDF@data$id)
  colnames(coord) <- c("X", "Y", "id")
  data <- dplyr::left_join(data, coord, by = "id")

  lvl1_data <- data

  if(huge_data_size)
  {
    message("Data Prepared! Go!............................................\n")
  }

  # GWPR test
  if (adaptive)
  {
    result <- gwpr_A_phtest(bw = bw, data = lvl1_data, SDF = SDF, index = index, ID_list = ID_num,
                            random.method = random.method, formula = formula, p = p, longlat = longlat,
                            adaptive = adaptive, kernel = kernel, effect = effect, huge_data_size = huge_data_size)

  }
  else
  {
    result <- gwpr_F_phtest(bw = bw, data = lvl1_data, SDF = SDF, index = index, ID_list = ID_num,
                            random.method = random.method, formula = formula, p = p, longlat = longlat,
                            adaptive = adaptive, kernel = kernel, effect = effect, huge_data_size = huge_data_size)
  }
  return(result)
}
