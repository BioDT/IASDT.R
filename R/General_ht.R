## |------------------------------------------------------------------------| #
# ht ----
## |------------------------------------------------------------------------| #

#' Print head and tail of data frame
#'
#' Print head and tail of data frame
#' @name ht
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param DF data frame to print
#' @param NRows N umber of rows to print at the top and bottom of data frame
#' @examples
#' ht(mtcars)
#'
#' # -------------------------------------------
#'
#' ht(mtcars, 2)
#'
#' # -------------------------------------------
#'
#' ht(mtcars, 6)
#' @export

ht <- function(DF, NRows = 5) {
  DF %>%
    data.table::data.table() %>%
    print(topn = NRows)
}
