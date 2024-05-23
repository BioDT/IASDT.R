## |------------------------------------------------------------------------| #
# LoadedPackages ----
## |------------------------------------------------------------------------| #
#
#' List of loaded packages
#'
#' List of loaded packages
#' @name LoadedPackages
#' @examples
#' LoadedPackages()
#'
#' LoadPackages(sf)
#'
#' LoadedPackages()
#'
#' @export

LoadedPackages <- function() {
  (.packages())
}
