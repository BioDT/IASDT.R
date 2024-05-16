# |---------------------------------------------------| #
# ClearConsole ----
# |---------------------------------------------------| #
#
#' Clear the console
#'
#' Clear the console
#' This function is a lazy equivalent of `cat("\014")`
#'
#' @name ClearConsole
#' @export
#' @examples
#' \dontrun{
#' ClearConsole()
#' }

ClearConsole <- function() {
  cat("\014")
}
