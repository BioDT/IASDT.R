## |------------------------------------------------------------------------| #
# ht ----
## |------------------------------------------------------------------------| #

#' Print head and tail of data frame
#'
#' This function takes a data frame and an optional number of rows to print from
#' both the top (head) and bottom (tail) of the data frame. It is useful for
#' quickly inspecting the first few and last few rows of a large data frame to
#' understand its structure and contents.
#' @name ht
#' @author Ahmed El-Gabbas
#' @return The function is used for its side effect (printing) and does not
#'   return any value.
#' @param DF `data.frame. A data frame to print. This parameter cannot be
#'   `NULL`.
#' @param NRows Integer. Number of rows to print from both the head and tail of
#'   the data frame. Defaults to 5.
#' @examples
#' ht(mtcars)
#'
#' # -------------------------------------------
#'
#' ht(DF = mtcars, 2)
#'
#' # -------------------------------------------
#'
#' ht(DF = mtcars, 6)
#' @export

ht <- function(DF = NULL, NRows = 5L) {

  if (is.null(DF)) {
    stop("DF cannot be NULL", call. = FALSE)
  }

  DF %>%
    data.table::data.table() %>%
    print(topn = NRows)

  return(invisible(NULL))
}
