## |------------------------------------------------------------------------| #
# NUnique ----
## |------------------------------------------------------------------------| #

#' Number of unique values for all columns of a data frame
#'
#' This function calculates the number of unique values for each column in a given data frame and returns a data frame with two columns: Variable and NUnique. The Variable column lists the names of the original columns, and the NUnique column lists the corresponding number of unique values in each column. The result is sorted by the number of unique values in ascending order.
#'
#' @name NUnique
#' @param Data A data frame for which the number of unique values per column will be calculated.
#' @references [Click here](https://stackoverflow.com/questions/22196078/count-unique-values-for-every-column)
#' @export
#' @return A data frame with two columns: Variable and NUnique. The Variable column lists the names of the original columns, and the NUnique column lists the number of unique values in each column. The result is sorted by NUnique in ascending order.
#' @examples
#' NUnique(mtcars)
#'
#' NUnique(iris)

NUnique <- function(Data) {

  if (is.null(Data)) {
    stop("Data cannot be NULL")
  }

  Data %>%
    dplyr::summarise(
      dplyr::across(tidyselect::everything(), dplyr::n_distinct)) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(),
        names_to = "Variable", values_to = "NUnique") %>%
    dplyr::arrange(NUnique) %>%
    return()
}
