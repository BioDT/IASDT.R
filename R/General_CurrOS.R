## |------------------------------------------------------------------------| #
# CurrOS ----
## |------------------------------------------------------------------------| #

#' Current operating system
#'
#' This function returns the name of the current operating system the R session
#' is running on.
#' @name CurrOS
#' @author Ahmed El-Gabbas
#' @return A character string representing the name of the operating system.
#' @examples
#' CurrOS()
#'
#' @export

CurrOS <- function() {
  as.character(Sys.info()["sysname"])
}
