## |------------------------------------------------------------------------| #
# System ----
## |------------------------------------------------------------------------| #

#' Run a system command in a cross-platform manner
#'
#' This function executes a system command, using either `shell` on Windows or `system` on Linux. It allows the output of the command to be captured into an R object.
#'
#' @param command A string representing the bash command to be executed.
#' @param RObj A logical indicating whether to capture the output of the command as an R object. If `TRUE`, the output is captured; if `FALSE`, the output is printed to the console. Defaults to `TRUE`.
#' @param ... Additional arguments passed to either `shell` or `system` function, depending on the operating system.
#' @name System
#' @author Ahmed El-Gabbas
#' @return Depending on the value of `RObj`, either the output of the executed command as an R object or `NULL` if `RObj` is `FALSE` and the output is printed to the console.
#' @examples
#' # print working directory
#' System("pwd")
#'
#' # first 5 files on the working directory
#' (A <- System("ls | head -n 5"))
#'
#' (A <- System("ls | head -n 5", RObj = FALSE))
#' @export

System <- function(command, RObj = TRUE, ...) {
  # Use shell() in windows OS; system() for linux OS
  # The running operating system (make also available in the global environment outside of the function)

  if (IASDT.R::CurrOS() == "Windows") {
    Out <- shell(cmd = command, intern = RObj, ...)
  }
  if (IASDT.R::CurrOS() == "Linux") {
    Out <- system(command = command, intern = RObj, ...)
  }
  return(Out)
}
