#' Locally Breusch-Pagan Lagrange Multiplier Test Based on Fixed Distance Bandwidth
#'
#' @param bw               The optimal bandwidth, fixed distance
#' @param data             A data.frame for the Panel data
#' @param SDF              Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param index            A vector for the indexes : (c("ID", "Time"))
#' @param ID_list          The data.frame with individuals' ID
#' @param formula          The regression formula: : Y ~ X1 + ... + Xk
#' @param p                The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param longlat          If TRUE, great circle distances will be calculated
#' @param adaptive         If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param kernel           bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param huge_data_size   If TRUE, the "progress_bar" function will be launched
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom  plm plm pdata.frame plmtest
#' @importFrom  sp merge
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points, test value, p value, df1, df2}
#' }
#' @noRd
gwpr_F_plmtest <- function(bw = bw, data , SDF, index, ID_list,
                           formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                           kernel = kernel, huge_data_size = huge_data_size)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat,
                       test.name = "Breusch-Pagan Lagrange Multiplier Test",
                       result.explain = "If the p-value is lower than the specific level (0.01, 0.05, etc.), significant effects exist.")
  message("**************** Breusch-Pagan Lagrange Multiplier test in each subsample *******************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: Pooling"," ---- ", "\n",
          "If the p-value is lower than the specific level (0.01, 0.05, etc.), significant effects exist.\n")
  ID_list_single <- as.vector(ID_list[[1]])
  output_result <- data.frame(Doubles = double(), Characters = character())
  loop_times <- 1
  wgt <- 0
  for (ID_individual in ID_list_single)
  {
    data$aim[data$id == ID_individual] <- 1
    data$aim[data$id != ID_individual] <- 0
    subsample <- data
    subsample <- subsample[order(-subsample$aim),]
    dp_locat_subsample <- dplyr::select(subsample, 'X', 'Y')
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- weight[,1]
    subsample <- subsample[(subsample$wgt > 0),]
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = FALSE)
    plm_subsample_ols <- plm::plm(formula=formula, model="pooling", data=Psubsample,
                                  index=index, weights = wgt)
    test <- plm::plmtest(plm_subsample_ols, type=c("bp"))
    result_line <- c(ID_individual, test$statistic, test$p.value, test$parameter)
    output_result <- rbind(output_result, result_line)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  colnames(output_result) <- c("id", "statistic", "p.value", "df")
  SDF <- sp::merge(SDF, output_result, by = "id")
  result_list <- list(GW.arguments = GW.arguments, SDF = SDF)
  return(result_list)
}
