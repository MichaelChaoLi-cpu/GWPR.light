#' Protect each Subsample with Least Individuals to Run the Regression
#'
#' @description When we use random effect model and "swar" method, the model required number of individuals the same as the estimated parameters in the subsamples. Moreover, if the number of individuals in a subsample is one, then there is no difference between panel regression and time-serious regression. At least 2. Therefore, we use this function to protect the model.
#'
#' @param data      The data.frame has been washed
#' @param ID_list   The data.frame with individuals' ID
#' @param index     The index C("id", "time"), here "id" is always "id", but "time" is set by user
#' @param kernel    Kernel, default "bisquare". gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param p         The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat   If TRUE, great circle distances will be calculated (see GWmodel::bw.gwr)
#' @param bw_panel  The  required numbers of individuals
#'
#' @import dplyr
#' @import GWmodel
#'
#' @return A distance to gaurantee there is enought individuals in every subsample
#' @noRd
protect_model_with_least_individuals <- function(data, ID_list, index,
                                                 kernel, p, longlat, bw_panel)
{
  ID_list_single <- as.vector(ID_list[[1]])
  max_dist <- c()
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
    subsample$dist <- as.vector(dMat)
    subsample <- subsample[order(subsample$dist),]
    id_subsample <- dplyr::select(subsample, "id")
    id_subsample <- id_subsample[!duplicated(id_subsample$id),]
    id_subsample <- as.data.frame(id_subsample)
    id_subsample <- id_subsample[1:bw_panel,]
    id_subsample <- as.data.frame(id_subsample)
    colnames(id_subsample) <- "id"
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    max_dist <- append(max_dist ,max(subsample$dist))
  }
  lower <- max(max_dist) * 1.011 # because individuals with the weight lower than 0.01 would be ignored,
  # to guarantee all the individuals used in local panel model, we use 1.011 here.
  return(lower)
}
