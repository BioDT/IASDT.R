## |------------------------------------------------------------------------| #
# CurrOS ----
## |------------------------------------------------------------------------| #

#' Current operating system
#'
#' Current operating system
#'
#' @name CurrOS
#' @author Ahmed El-Gabbas
#' @return String for the current operating system
#' @examples
#' CurrOS()
#'
#' @export

CurrOS <- function() {
  as.character(Sys.info()["sysname"])
}
