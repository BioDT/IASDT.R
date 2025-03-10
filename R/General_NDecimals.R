## |------------------------------------------------------------------------| #
# NDecimals ----
## |------------------------------------------------------------------------| #

#' Number of decimal places in a numeric value
#'
#' This function calculates the number of decimal places in a numeric value. It
#' is designed to work with numeric inputs that can be coerced to character
#' format.
#'
#' @param x Numeric (or character) numeric value.
#' @name NDecimals
#' @author Ahmed El-Gabbas
#' @return An integer representing the number of decimal places in the input
#'   value. If the input value does not have any decimal places, the function
#'   returns 0.
#' @examples
#' NDecimals(x = "13.45554545")
#'
#' # -------------------------------------------
#'
#' NDecimals(x = 15.01500)
#'
#' NDecimals(x = '15.01500')
#'
#' # -------------------------------------------
#'
#' NDecimals(x = 13.45554545)
#' @export

NDecimals <- function(x) {

  if (is.null(x)) {
    stop("x cannot be NULL", call. = FALSE)
  }

  Split <- x %>%
    as.character() %>%
    format(scientific = FALSE) %>%
    stringr::str_split(pattern = "\\.", n = Inf, simplify = TRUE)

  if (ncol(Split) == 2) {
    Out <- Split %>%
      as.vector() %>%
      utils::tail(1) %>%
      nchar() %>%
      as.integer()
    return(Out)
  } else {
    return(0L)
  }
}
