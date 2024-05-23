## |------------------------------------------------------------------------| #
# NDecimals ----
## |------------------------------------------------------------------------| #

#' Number of decimal places in a vector
#'
#' Number of decimal places in a vector
#' @param x vector
#' @name NDecimals
#' @author Ahmed El-Gabbas
#' @return integer; number of decimal places
#' @examples
#' NDecimals("13.45554545")
#'
#' # -------------------------------------------
#'
#' NDecimals(15.01500)
#'
#' NDecimals('15.01500')
#'
#' # -------------------------------------------
#'
#' NDecimals(13.45554545)
#' @export

NDecimals <- function(x) {
  Split <- x %>%
    format(scientific = FALSE) %>%
    stringr::str_split(pattern = "\\.", n = Inf, simplify = TRUE)

  if (length(Split) == 2) {
    Split[, 2] %>%
      nchar() %>%
      as.integer() %>%
      return()
  } else {
    return(0)
  }
}
