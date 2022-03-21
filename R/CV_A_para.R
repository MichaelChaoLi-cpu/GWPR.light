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

CV_A_para.step <- function(bw, data, ID_list, formula, p, longlat, adaptive, kernel,
                           model = model, index = index, effect = effect,
                           random.method = random.method, cluster.number = cluster.number)
{
  ID_list_single <- as.vector(ID_list[[1]])
  wgt <- 0
  ID_individual <- 0
  varibale_name_in_equation <- all.vars(formula)
  cl <- parallel::makeCluster(cluster.number)
  doParallel::registerDoParallel(cl)
  #  v0.1.1 the loss function is based on local r2
  # CVscore_vector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
  # v0.1.2
  residualsVector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
    {
      data$aim[data$id == ID_individual] <- 1
      data$aim[data$id != ID_individual] <- 0
      subsample <- data
      #v0.1.2
      numberOfAim <- nrow(subsample[subsample$aim == 1,])
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
      id_subsample <- as.data.frame(id_subsample)
      colnames(id_subsample) <- "id"
      id_subsample$flag <- 1
      subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
      bw_to_total <- nrow(subsample)
      weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_to_total, kernel=kernel, adaptive=adaptive)
      subsample$wgt <- as.vector(weight)
      Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                     stringsAsFactors = FALSE)
      plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                                effect = effect, index=index, weights = wgt,
                                random.method = random.method)
      # v0.1.1
      #    if(!inherits(plm_subsample, "try-error"))
      #    {
      #      CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
      #        (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
      #    }
      #    else
      #    {
      #      CVscore <- Inf
      #    }
      #0.2.0
      if(!inherits(plm_subsample, "try-error"))
      {
        residualsLocalAim <-  plm_subsample$residuals[1:numberOfAim]
      }
      else
      {
        residualsLocalAim <- Inf
      }
    }
  parallel::stopCluster(cl)
  #  v0.1.1 the loss function is based on local r2
  #  mean_CVscore <- mean(CVscore_vector)
  #  cat("Fixed Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  #  return(mean_CVscore)
  #v0.2.0
  CVscore <- nrow(data) * sum(residualsVector^2) /
    (nrow(data) - length(varibale_name_in_equation) + 1)^2
  cat("Fixed Bandwidth:", bw, "CV score:", CVscore, "\n")
  return(CVscore)
}
