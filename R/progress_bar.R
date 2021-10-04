#' Progress Bar Function
#'
#' @description if big data size, this function show the progress bar
#'
#' @param loop_times     Current loop times
#' @param total_times    Total loop times
#' @noRd
progress_bar <- function(loop_times, total_times)
{
  flag_dash <- floor(loop_times/total_times*100/2.5) - floor((loop_times-1)/total_times*100/2.5)
  if (flag_dash == 1)
  {
    cat("-")
  }
  flag_vert <- floor(loop_times/total_times*100/25) - floor((loop_times-1)/total_times*100/25)
  if (flag_vert ==1)
  {
    cat("|")
  }
  if (loop_times == total_times)
  {
    cat("\n")
  }
}
