#' CV Score Calculator with Adaptive Distance Bandwidth
#'
#' @description Get CV score with an adaptive distance bandwidth
#'
#' @param bw                Current potential bandwidth put into the calculation.
#' @param data              The data.frame has been washed
#' @param ID_list           The data.frame with individuals' ID
#' @param formula           The regression formula: : Y ~ X1 + ... + Xk
#' @param p                 The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat           If TRUE, great circle distances will be calculated (see GWmodel::bw.gwr)
#' @param adaptive          If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param kernel            Kernel, default "bisquare". gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param model             Panel models transformation : (c("within", "random", "pooling"))
#' @param index             The index C("id", "time"), here "id" is always "id", but "time" is set by user
#' @param effect            The effects introduced in the model, one of "individual", "time", "twoways", or "nested"
#' @param random.method     Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param huge_data_size    If TRUE, the "progress_bar" function will be launched
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom  plm plm pdata.frame
#'
#' @return A CV score
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#' @noRd
CV_A <- function(bw, data, ID_list, formula, p, longlat, adaptive, kernel,
                 model, index, effect,
                 random.method, huge_data_size)
{
  CVscore_vector <- c()
  ID_list_single <- as.vector(ID_list[[1]])
  loop_times <- 1
  varibale_name_in_equation <- all.vars(formula)
  for (ID_individual in ID_list_single)
  {
    subsample <- dplyr::mutate(data, aim = ifelse(id == ID_individual, 1, 0))
    subsample <- dplyr::arrange(subsample, desc(aim))
    dp_locat_subsample <- dplyr::select(subsample, X, Y)
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    subsample$dist <- as.vector(dMat)
    subsample <- dplyr::arrange(subsample, dist)
    id_subsample <- dplyr::select(subsample, "id")
    id_subsample <- dplyr::distinct(id_subsample, id)
    id_subsample <- id_subsample[1:bw,]
    id_subsample <- as.data.frame(id_subsample)
    colnames(id_subsample) <- "id"
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    bw_to_total <- nrow(subsample)
    weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_to_total, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
      (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
    CVscore_vector <- append(CVscore_vector, CVscore)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  mean_CVscore <- mean(CVscore_vector)
  cat("Adaptive Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  return(mean_CVscore)
}
