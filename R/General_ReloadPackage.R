## |------------------------------------------------------------------------| #
# ReloadPackage ----
## |------------------------------------------------------------------------| #
#
#' Reload an R package
#'
#' Reload an R package
#' @name ReloadPackage
#' @param Package The name of packages to reload
#' @examples
#' LoadPackages(sf)
#'
#' ReloadPackage("sf")
#'
#' @export

ReloadPackage <- function(Package = NULL) {
  if (is.null(Package)) {
    stop()
  }

  if (!Package %in% LoadedPackages()) {
    library(package = Package, character.only = TRUE)
  } else {
    PackagesFolders <- paste0(.libPaths(), "/")
    PackagesFolders <- c(outer(PackagesFolders, Package, FUN = "paste0"))
    PackagesFolders <- PackagesFolders[file.exists(PackagesFolders)]
    PackagesFolders <- gsub(pattern = "/", replacement = "//", x = PackagesFolders)
    lapply(PackagesFolders, function(x) {
      devtools::reload(x, quiet = FALSE)
    })
  }

  return(invisible(NULL))
}
