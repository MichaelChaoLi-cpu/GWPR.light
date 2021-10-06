#' @title Drop the Individuals with only one observation
#'
#' @description Drop the individual with only one observation when using "within" to avoid error in merging data frame and SDF
#'
#' @param data        The data.frame of the raw dataset
#' @param ID_num      The data.frame includes individuals' ID and the numbers of records (always named: Count)
#'
#' @return A data.frame that all the individuals have more than one observations
#' @import dplyr
#' @noRd
drop_ID_with_single_observation <- function(data, ID_num)
{
  data <- dplyr::left_join(data, ID_num, by = "id")
  data <- data[(data$Count != 1),]
  data <- dplyr::select(data, -"Count")
  return(data)
}
