## |------------------------------------------------------------------------| #
# loaded_packages ----
## |------------------------------------------------------------------------| #
#
#' List of loaded packages
#'
#' This function returns a character vector listing all the packages that are
#' currently loaded in the R session.
#'
#' @return A character vector containing the names of all loaded packages.
#' @name loaded_packages
#' @examples
#' loaded_packages()
#'
#' require(Hmsc)
#'
#' loaded_packages()
#'
#' @export

loaded_packages <- function() {
  Packages <- .packages()
  return(Packages)
}
