## |------------------------------------------------------------------------| #
# NUnique ----
## |------------------------------------------------------------------------| #

#' Number of unique values for all columns of a data frame
#'
#' Number of unique values for all columns of a data frame
#'
#' @name NUnique
#' @param Data input data
#' @references https://stackoverflow.com/questions/22196078/count-unique-values-for-every-column
#' @export
#' @examples
#' NUnique(mtcars)
#'
#' NUnique(iris)

NUnique <- function(Data) {
  Data %>%
    dplyr::summarise(
      dplyr::across(tidyselect::everything(), dplyr::n_distinct)) %>%
    tidyr::pivot_longer(tidyselect::everything()) %>%
    stats::setNames(c("Variable", "NUnique")) %>%
    dplyr::arrange(NUnique)
}
