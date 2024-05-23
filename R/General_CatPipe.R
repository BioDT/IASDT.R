## |------------------------------------------------------------------------| #
# CatPipe ----
## |------------------------------------------------------------------------| #

#' print a message with time in the middle of the pipe
#'
#' print a message with time in the middle of the pipe
#' @param x input object to pipe
#' @param message printed message
#' @name CatPipe
#' @author Ahmed El-Gabbas
#' @return NULL
#' @keywords internal
#' @noRd

CatPipe <- function(x, message = NULL) {
  CatTime(message)
  return(x)
}
