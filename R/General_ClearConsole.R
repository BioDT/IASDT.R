## |------------------------------------------------------------------------| #
# ClearConsole ----
## |------------------------------------------------------------------------| #
#
#' Clear the console
#'
#' Clear the console. This function is equivalent to `cat("\014")`
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
