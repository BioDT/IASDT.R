## |------------------------------------------------------------------------| #
# package_functions ----
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
#' @param package Character. Package name.
#' @name package_functions
#' @export
#' @examples
#' str(package_functions(package = "raster"))
#'
#' str(package_functions(package = "sf"))

package_functions <- function(package) {

  if (is.null(package)) {
    stop("`package` cannot be NULL or empty", call. = FALSE)
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop("package", package, "not found", call. = FALSE)
  }
  library(
    package = eval(package), character.only = TRUE, quietly = TRUE,
    verbose  = FALSE)
  return(ls(paste0("package:", package)))

}
