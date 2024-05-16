# |---------------------------------------------------| #
# System ----
# |---------------------------------------------------| #

#' Run bash script depending on the operating system
#'
#' Run bash script depending on the operating system
#' @param command bash command to implement
#' @param RObj whether to make the output of the command an R object
#' @param ... additional arguments to `shell` pr `system` functions
#' @name System
#' @author Ahmed El-Gabbas
#' @return NULL
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
