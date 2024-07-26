## |------------------------------------------------------------------------| #
# ReloadPackage ----
## |------------------------------------------------------------------------| #
#
#' Reload an R package
#'
#' This function reloads a specified R package. If the package is not already loaded, it attempts to load it. If the package is already loaded, it reloads the package from its library location. This can be useful during package development when changes are made to the package code. This function depends on the [devtools::reload] function.
#'
#' @name ReloadPackage
#' @param Package A character string specifying the name of the package to reload.
#' @return The function is used for its side effect of reloading a package rather than for its return value.
#' @examples
#' LoadPackages(sf)
#'
#' ReloadPackage("sf")
#'
#' @export

ReloadPackage <- function(Package) {

  if (is.null(Package)) {
    stop("Package name cannot be NULL")
  }

  if (!requireNamespace(Package, quietly = TRUE)) {
    library(package = Package, character.only = TRUE)
  } else {

    PackagesFolders <- paste0(.libPaths(), "/", Package)
    PackagesFolders <- PackagesFolders[file.exists(PackagesFolders)]
    if (length(PackagesFolders) > 0) {
      purrr::walk(
      .x = PackagesFolders, .f = ~ devtools::reload(pkg = .x, quiet = FALSE))
  }
}
  return(invisible(NULL))
}
