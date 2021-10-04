#' Protect Local Models with Enough Freedom
#'
#' @description Freedom is the requirement. We need to guarantee the model has enough freedom to estimate the parameters.
#'
#' @param formula   The regression formular: : Y ~ X1 + ... + Xk
#' @param data      The data.frame has been washed
#' @param ID_list   The data.frame with individuals' ID
#' @param index     The index C("id", "time"), here "id" is always "id", but "time" is set by user
#' @param p         The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat   If TRUE, great circle distances will be calculated (see GWmodel::bw.gwr)
#'
#' @import dplyr
#' @import GWmodel
#'
#' @return A minimum of the number of observation in each subsample
#' @noRd
protect_model_with_enough_freedom <- function(formula, data, ID_list, index,
                                              p, longlat)
{
  ID_list_single <- as.vector(ID_list[[1]])
  step_increase_lower <- 1
  lower <- 1
  go_out <- T
  required_freedom <- length(all.vars(formula))
  while((step_increase_lower < 1001)&(go_out == T))
  {
    step_increase_lower <- step_increase_lower + 1 #required at least two individuals
    lower <- lower + 1
    for (ID_individual in ID_list_single)
    {
      go_out <- F
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
      id_subsample <- id_subsample[1:lower,]
      id_subsample <- as.data.frame(id_subsample)
      colnames(id_subsample) <- "id"
      id_subsample <- dplyr::mutate(id_subsample, flag = 1)
      subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
      if(nrow(subsample) < required_freedom)
      {
        go_out <- T
        break
      }
    }
  }
  return(lower)
}
