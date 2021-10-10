#' Geographically Weighted Panel Regression Based on the Optimal Adaptive Distance Bandwidth
#'
#' @param bw                The optimal adaptive bandwidth
#' @param data              The data.frame has been washed
#' @param SDF               Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param ID_list           The data.frame with individuals' ID
#' @param formula           The regression formula: : Y ~ X1 + ... + Xk
#' @param p                 The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat           If TRUE, great circle distances will be calculated
#' @param adaptive          If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param model             Panel model transformation: (c("within", "random", "pooling"))
#' @param index             A vector for the indexes : (c("ID", "Time"))
#' @param kernel            bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                           gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                           exponential: wgt = exp(-vdist/bw);
#'                           tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                           boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param effect            The effects introduced in the model, one of "individual" (default) , "time", "twoways", or "nested"
#' @param random.method     Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param huge_data_size    If TRUE, the "progress_bar" function will be launched
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom  plm plm pdata.frame r.squared
#' @importFrom lmtest coeftest
#' @import sp
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{R2}{global r2}
#' \item{index}{the index used in the result, Note: in order to avoid mistakes, we forced a rename of the individuals'ID as id.}
#' \item{plm.result}{an object of class inheriting from plm, see plm}
#' \item{raw.data}{the data.frame used in the regression}
#' \item{GWPR.residuals}{the data.frame includes Y, Y hat, and residuals from GWPR}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points,GWPR coefficient estimates,coefficient standard errors and t-values in its "data" slot.}
#' }
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#' @noRd
gwpr_A <- function(bw, data, SDF, ID_list, formula, p, longlat, adaptive,
                   model, index, kernel = "bisquare", effect = "individual",
                   random.method = "swar", huge_data_size = huge_data_size)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat)
  message("************************ GWPR Begin *************************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: ", model, " ---- ", "Effect: ", effect, "\n")
  global_plm <- plm::plm(formula=formula, model=model, data=data,
                         effect = effect, index=index, random.method = random.method)
  ID_list_single <- as.vector(ID_list[[1]])
  output_result <- data.frame(Doubles = double())
  y_yhat_resid <- data.frame(Doubles = double())
  loop_times <- 1
  wgt = 0
  for (ID_individual in ID_list_single)
  {
    data_input$aim[data_input$id == ID_individual] <- 1
    data_input$aim[data_input$id != ID_individual] <- 0
    subsample <- data_input
    subsample <- subsample[order(-subsample$aim),]
    dp_locat_subsample <- dplyr::select(subsample, 'X', 'Y')
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    subsample$dist <- as.vector(dMat)
    subsample <- subsample[order(subsample$dist),]
    id_subsample <- dplyr::select(subsample, "id")
    id_subsample <- id_subsample[!duplicated(id_subsample$id),]
    id_subsample <- as.data.frame(id_subsample)
    id_subsample <- id_subsample[1:bw,]
    colnames(id_subsample) <- "id"
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    bw_of_total <- nrow(subsample)
    weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_of_total, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    coefMat <- lmtest::coeftest(plm_subsample)
    local_r2 <- plm::r.squared(plm_subsample)
    result_line <- c(ID_individual, coefMat[,1], coefMat[,2], coefMat[,3], local_r2)
    output_result <- rbind(output_result, result_line)
    dataset_add_resid <- cbind(Psubsample, plm_subsample$residuals)
    dataset_add_resid <- as.data.frame(dataset_add_resid)
    varibale_name_in_equation <- all.vars(formula)
    dataset_add_resid <- dplyr::select(dataset_add_resid, dplyr::all_of(index), dplyr::all_of(varibale_name_in_equation)[1],
                                       "plm_subsample$residuals")
    colnames(dataset_add_resid) <- c(index, "y", "resid")
    dataset_add_resid$yhat <- dataset_add_resid$y - dataset_add_resid$resid
    dataset_add_resid <- dplyr::filter(dataset_add_resid, id == ID_individual)
    y_yhat_resid <- rbind(y_yhat_resid, dataset_add_resid)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  varibale_name_in_equation <- all.vars(formula)
  if (model == "within")
  {
    varibale_name_in_equation_out <- varibale_name_in_equation[2:length(varibale_name_in_equation)]
  }
  else
  {
    varibale_name_in_equation_out <- varibale_name_in_equation
    varibale_name_in_equation_out[1] <- "Intercept"
  }
  colnames(output_result) <- c("id", varibale_name_in_equation_out, paste0(varibale_name_in_equation_out,"_SE"),
                               paste0(varibale_name_in_equation_out,"_TVa"), "Local_R2")
  SDF <- sp::merge(SDF, output_result, by = "id")
  y_yhat_resid[,1] <- as.numeric(as.character(y_yhat_resid[,1]))
  y_yhat_resid[,2] <- as.numeric(as.character(y_yhat_resid[,2]))
  r2 <- 1 - sum(y_yhat_resid$resid^2)/(sum((y_yhat_resid$y - mean(y_yhat_resid$y))^2))
  result_list <- list(GW.arguments = GW.arguments, R2 = r2, index = index, plm.result = global_plm,
                      raw.data = data, GWPR.residuals = y_yhat_resid, SDF = SDF)
  return(result_list)
}
