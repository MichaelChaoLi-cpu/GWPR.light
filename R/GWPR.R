#' Geographically Weighted Panel Regression Model
#'
#' @description This function implements GWPR
#'
#' @usage GWPR(formula, data, index, SDF, bw = NULL, adaptive = F, p = 2,
#'             effect = "individual", model = c("pooling", "within", "random"),
#'             random.method = "swar", kernel = "bisquare", longlat = F)
#'
#' @param formula        The regression formula: : Y ~ X1 + ... + Xk
#' @param data           A data.frame for the Panel data
#' @param index          A vector for the indexes : (c("ID", "Time"))
#' @param SDF            Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param bw             The optimal bandwidth, either adaptive or fixed distance
#' @param adaptive       If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param p              The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param effect         The effects introduced in the model, one of "individual" (default) , "time", "twoways", or "nested"
#' @param model          Panel model transformation: (c("within", "random", "pooling"))
#' @param random.method  Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param kernel         bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param longlat        If TRUE, great circle distances will be calculated
#'
#' @import dplyr
#' @import sp
#' @importFrom data.table setDT
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{R2}{global r2}
#' \item{index}{the index used in the result, Note: in order to avoid mistakes, we forced a rename of the individuals'ID as id.}
#' \item{plm.result}{an object of class inheriting from plm, see plm}
#' \item{raw.data}{the data.frame used in the regression}
#' \item{GWPR.residuals}{the data.frame includes Y, Y hat, and residuals from GWPR}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points,GWPR coefficient estimates,coefficient standard errors and t-values in its data slot.}
#' }
#' @export
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi <managi.s@gmail.com>
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#'
#' @examples
#' \dontrun{
#' data(TransAirPolCalif)
#' data(California)
#' formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
#'    Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
#'    Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
#'    Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
#'    Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
#'    pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
#'
#' bw.AIC.F <- bw.GWPR(formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"), SDF = California,
#'                     adaptive = F, p = 2, bigdata = F, effect = "individual",
#'                     model = "within", approach = "AIC", kernel = "bisquare", longlat = F,
#'                     doParallel = T, cluster.number = 4)
#'
#' result.F.AIC <- GWPR(bw = bw.AIC.F, formula = formula.GWPR, data = TransAirPolCalif, index = c("GEOID", "year"),
#'                      SDF = California, adaptive = F, p = 2, effect = "individual", model = "within",
#'                      kernel = "bisquare", longlat = F)
#' summary(result.F.AIC$SDF$Local_R2)
#' library(tmap)
#' tm_shape(result.F.AIC$SDF) +
#' tm_polygons(col = "Local_R2", pal = "Reds",auto.palette.mapping = F,
#'             style = 'cont')
#' }
GWPR <- function(formula, data, index, SDF, bw = NULL, adaptive = F, p = 2,
                 effect = "individual", model = c("pooling", "within", "random"), random.method = "swar",
                 kernel = "bisquare", longlat = F)
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
  if(!(model %in% c("pooling", "within", "random")))
  {
    stop("This version GWPR only accept \"pooling\", \"within\", or \"random\"")
  }

  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  data <- dplyr::select(data, index, varibale_name_in_equation)
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"

  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  ID <- dplyr::select(data, index[1])
  ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  if(model == "within")
  {
    data <- drop_ID_with_single_observation(data, ID_num)
    ID <- dplyr::select(data, index[1])
    ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  }

  # Judge the data size of calculation
  if (nrow(ID_num) > 1000)
  { # for test 40, real number should be 1000
    cat("Dear my friend, thanks for your patience!. We pass the bandwidth\n",
        "selection part. Now, regression! This should be faster. Thanks.\n",
        ".............................................................\n")
    huge_data_size <- T
  }
  else
  {
    huge_data_size <- F
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
    cat("Data Prepared! Go!............................................\n")
  }

  # GWPRegression
  if (adaptive)
  {
    result <- gwpr_A(bw = bw, data = lvl1_data, SDF, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
  }
  else
  {
    result <- gwpr_F(bw = bw, data = lvl1_data, SDF, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
  }
  cat("The R2 is :", result$R2,"\n")
  cat("Note: in order to avoid mistakes, we forced a rename of the individuals'ID as \"id\". \n")
  return(result)
}
