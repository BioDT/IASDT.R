# |---------------------------------------------------| #
# FunctionsInPackage ----
# |---------------------------------------------------| #
#
#' List of functions in a package
#'
#' List of functions in a package
#' @name FunctionsInPackage
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param Package name of the package
#' @export
#' @examples
#' FunctionsInPackage("raster")

FunctionsInPackage <- function(Package) {
  library(package = Package, character.only = TRUE)
  ls(paste0("package:", Package))
}
