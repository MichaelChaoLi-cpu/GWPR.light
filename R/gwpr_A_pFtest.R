#' Locally F Test Based on Adaptive Distance Bandwidth
#'
#' @description This function perform F test in each regression based on different subsamples divided by the adaptive distance bandwidth.
#'
#' @param bw               The optimal bandwidth, adaptive distance
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
#' @param effect           The effects introduced in the model, one of "individual" (default) , "time", "twoways", or "nested"
#' @param huge_data_size   If TRUE, the "progress_bar" function will be launched
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom  plm plm pdata.frame pFtest
#' @importFrom  sp merge
#'
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points, test value, p value, df1, df2}
#' }
#' @noRd
gwpr_A_pFtest <- function(bw = bw, data, SDF, index, ID_list, formula = formula, p = p, longlat = longlat,
                          adaptive = adaptive, kernel = kernel, effect = effect,
                          huge_data_size = huge_data_size)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat,
                       test.name = "F test for individual effects",
                       result.explain = "If the p-value is lower than the specific level (0.01, 0.05, etc.), significant effects exist.")
  cat("**************************** F test in each subsample *********************************\n")
  cat("Formula:", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals:", nrow(ID_list), "\n")
  cat("Bandwidth:", bw, " ---- ", "Adaptive:", adaptive, "\n")
  cat("Model: Fixed Effects vs Pooling"," ---- ", "Effect:", effect, "\n")
  cat("If the p-value is lower than the specific level (0.01, 0.05, etc.), significant effects exist.\n")
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
    subsample$distance <- dMat[,1]
    subsample <- subsample[order(subsample$distance),]
    in_subsample_id <- dplyr::select(subsample, index[1])
    in_subsample_id <- dplyr::distinct(in_subsample_id)
    in_subsample_id$yes <- 1
    in_subsample_id <- in_subsample_id[1:bw, ]
    bw_panel <- dplyr::left_join(ID_list, in_subsample_id, by = index[1])
    bw_panel$usingCount <- bw_panel$Count * bw_panel$yes
    bw_panel <- sum(bw_panel$usingCount, na.rm = T)
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw_panel, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- weight[,1]
    subsample <- subsample[(subsample$wgt > 0),]
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample_fem <- plm::plm(formula=formula, model="within", data=Psubsample,
                                  effect = effect, index=index, weights = wgt)
    plm_subsample_ols <- plm::plm(formula=formula, model="pooling", data=Psubsample,
                                  index=index, weights = wgt)
    test <- plm::pFtest(plm_subsample_fem, plm_subsample_ols)
    result_line <- c(ID_individual, test$statistic, test$p.value, test$parameter)
    output_result <- rbind(output_result, result_line)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  colnames(output_result) <- c("id", "statistic", "p.value", "df1", "df2")
  SDF <- sp::merge(SDF, output_result, by = "id")
  result_list <- list(GW.arguments = GW.arguments, SDF = SDF)
  return(result_list)
}
