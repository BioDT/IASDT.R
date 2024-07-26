## |------------------------------------------------------------------------| #
# LoadedPackages ----
## |------------------------------------------------------------------------| #
#
#' List of loaded packages
#'
#' This function returns a character vector listing all the packages that are currently loaded in the R session.
#'
#' @return A character vector containing the names of all loaded packages.
#' @name LoadedPackages
#' @examples
#' LoadedPackages()
#'
#' require(Hmsc)
#'
#' LoadedPackages()
#'
#' @export

LoadedPackages <- function() {
  Packages <- .packages()
  return(Packages)
}
