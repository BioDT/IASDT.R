# |---------------------------------------------------| #
# GetMode ----
# |---------------------------------------------------| #

#' Calculate the mode of numbers
#'
#' Calculate the mode of numbers
#' @param v vector
#' @name GetMode
#' @references https://www.tutorialspoint.com/r/r_mean_median_mode.htm
#' @return NULL
#' @examples
#' GetMode(c(1:10,1,1,3,3,3,3))
#' @export

GetMode <- function(v) {
  unique_vals <- unique(v)
  unique_vals[which.max(tabulate(match(v, unique_vals)))]
}
