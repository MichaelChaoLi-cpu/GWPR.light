#' CV Score Calculator with Adaptive Distance Bandwidth with Parallel Process
#'
#' @description Get CV score with an adaptive distance bandwidth
#'
#' @param bw              Current potential bandwidth put into the calculation.
#' @param data            The data.frame has been washed
#' @param ID_list         The data.frame with individuals' ID
#' @param formula         The regression formula: : Y ~ X1 + ... + Xk
#' @param p               The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat         If TRUE, great circle distances will be calculated (see GWmodel::bw.gwr)
#' @param adaptive        If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param kernel          Kernel, default "bisquare". gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param model           Panel models transformation : (c("within", "random", "pooling"))
#' @param index           The index C("id", "time"), here "id" is always "id", but "time" is set by user
#' @param effect          The effects introduced in the model, one of "individual", "time", "twoways", or "nested"
#' @param random.method   Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param cluster.number  Cluster number used in calculation
#'
#' @import dplyr
#' @import GWmodel
#' @import parallel
#' @import foreach
#' @import iterators
#' @import doParallel
#' @importFrom  plm plm pdata.frame
#'
#' @return A CV score
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#' @noRd
CV_A_para <- function(bw, data, ID_list = ID_num, formula, p, longlat, adaptive, kernel,
                      model, index, effect,
                      random.method,  cluster.number = cluster.number)
{
  ID_list_single <- as.vector(ID_list[[1]])
  varibale_name_in_equation <- all.vars(formula)
  cl <- parallel::makeCluster(cluster.number)
  doParallel::registerDoParallel(cl)
  CVscore_vector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
  {
    subsample <- dplyr::mutate(data, aim = ifelse(id == ID_individual, 1, 0))
    subsample <- dplyr::arrange(subsample, desc(aim))
    dp_locat_subsample <- dplyr::select(subsample, X, Y)
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    subsample$distance <- dMat[,1]
    subsample <- dplyr::arrange(subsample, distance)
    in_subsample_id <- dplyr::select(subsample, index[1])
    in_subsample_id <- dplyr::distinct(in_subsample_id)
    in_subsample_id$yes <- 1
    in_subsample_id <- in_subsample_id[1:bw, ]
    bw_panel <- dplyr::left_join(ID_list, in_subsample_id, by = index[1])
    bw_panel$usingCount <- bw_panel$Count * bw_panel$yes
    bw_panel <- sum(bw_panel$usingCount, na.rm = T)
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw_panel, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    subsample <- dplyr::filter(subsample, wgt > 0)
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
      (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
  }
  parallel::stopCluster(cl)
  mean_CVscore <- mean(CVscore_vector)
  cat("Adaptive Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  return(mean_CVscore)
}
