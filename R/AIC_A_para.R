#' AIC Score Calculator with Adaptive Distance Bandwidth with Parallel Process
#'
#' @description Get AIC score with an adaptive distance bandwidth
#'
#' @param bw              Current potential bandwidth put into the calculation.
#' @param data_input      The data.frame has been washed
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
#' @importFrom  stats sd aggregate
#'
#' @return A AIC score
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#' @noRd
AIC_A_para <- function(bw, data_input, ID_list, formula, p, longlat, adaptive, kernel,
                       model = model, index = index, effect = effect,
                       random.method = random.method, cluster.number = cluster.number)
{
  ID_list_single <- as.vector(ID_list[[1]])
  ID_individual <- 0
  wgt <- 0
  cl <- parallel::makeCluster(cluster.number)
  doParallel::registerDoParallel(cl)
  AICscore_vector <- foreach::foreach(ID_individual = ID_list_single, .combine = c) %dopar%
  {
    subsample <- data_input
    subsample$aim[subsample$id == ID_individual] <- 1
    subsample$aim[subsample$id != ID_individual] <- 0
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
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    bw_to_total <- nrow(subsample)
    weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_to_total, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = FALSE)
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    if(model == "within")
    {
      theta = 1
    }
    if(model == "pooling")
    {
      theta = 0
    }
    if(model == "random")
    {
      theta = as.numeric(plm_subsample$ercomp$theta)
    }
    varibale_name_in_equation <- all.vars(formula)
    indep_varibale_name_in_equation <- varibale_name_in_equation[2:length(varibale_name_in_equation)]
    X <- as.data.frame(Psubsample[,c("id", indep_varibale_name_in_equation)])
    X$id <- as.character(X$id)
    if((model == "random")|(model == "pooling"))
    {
      X$intercept <- 1
    }
    if((model == "random")|(model == "pooling"))
    {
      X$intercept <- 1
    }
    if (model == "pooling")
    {
      X_trans <- (dplyr::select(X, -"id"))
    }
    else
    {
      X_mean <- stats::aggregate(X[,indep_varibale_name_in_equation], by = list(X[,'id']), mean)
      colnames(X_mean)[1] <- "id"
      X_mean <- dplyr::left_join(dplyr::select(X, "id"), X_mean, by = "id")
      X_trans <- (dplyr::select(X, -"id")) - (dplyr::select(X_mean, -"id")) * theta
    }
    X_trans <- as.matrix(X_trans)
    W <- as.vector(Psubsample$wgt)
    P <- try(X_trans %*%  solve(t(X_trans) %*% (W * X_trans)) %*% t(X_trans) * W, silent=TRUE)
    if(!inherits(P, "try-error"))
    {
      tr_hatmat <- sum(diag(P))
      n <- nrow(Psubsample)
      AICscore <- 2*n*log(sd(plm_subsample$residuals)) + n*log(2*pi) +  n * (tr_hatmat + n) / (n - 2 - tr_hatmat)
    }
    else
    {
        AICscore <- Inf
    }
  }
  parallel::stopCluster(cl)
  mean_AICscore <- mean(AICscore_vector)
  cat("Adaptive Bandwidth:", bw, "AIC score:", mean_AICscore, "\n")
  return(mean_AICscore)
}
