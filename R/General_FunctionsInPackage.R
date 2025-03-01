## |------------------------------------------------------------------------| #
# FunctionsInPackage ----
## |------------------------------------------------------------------------| #
#
#' List of functions in a package
#'
#' This function returns a character vector listing all the functions available
#' in the specified R package. It first checks if the package is installed and
#' can be loaded; if not, it raises an error.
#' @author Ahmed El-Gabbas
#' @return A character vector containing the names of all functions in the
#'   specified package.
#' @param Package Character. Package name.
#' @name FunctionsInPackage
#' @export
#' @examples
#' str(FunctionsInPackage(Package = "raster"))
#'
#' str(FunctionsInPackage(Package = "sf"))

FunctionsInPackage <- function(Package) {
  if (is.null(Package)) {
    stop("Package cannot be NULL or empty", call. = FALSE)
  }

  if (!requireNamespace(Package, quietly = TRUE)) {
    stop("Package", Package, "not found", call. = FALSE)
  }
  library(
    package = eval(Package), character.only = TRUE, quietly = TRUE,
    verbose  = FALSE)
  return(ls(paste0("package:", Package)))
}
