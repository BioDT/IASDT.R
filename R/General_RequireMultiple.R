## |------------------------------------------------------------------------| #
# RequireMultiple ----
## |------------------------------------------------------------------------| #

#' Load (or install) multiple packages at once
#'
#' This function attempts to load multiple R packages specified by the user. If a package is not installed, the function can optionally install it before loading. It also provides an option to reload already loaded packages.
#'
#' @param ... character. Names of the packages to be loaded or installed.
#' @param Reload logical. Determines whether already loaded packages should be reloaded. Defaults to `FALSE`.
#' @param InstallMissing logical. Determines whether missing packages should be automatically installed and then loaded. Defaults to `TRUE`.
#' @return Invisible NULL. The function is used for its side effects (loading/installing packages) and does not return any value.
#' @author Ahmed El-Gabbas
#' @export
#' @name RequireMultiple
#' @examples
#' # Currently loaded packages
#' (P1 <- LoadedPackages())
#'
#' RequireMultiple(dismo, MASS, mgcv)
#'
#' # Loaded packages after implementing the function
#' (P2 <- LoadedPackages())
#'
#' # Which packages were loaded?
#' setdiff(P2, P1)

RequireMultiple <- function(..., Reload = FALSE, InstallMissing = TRUE) {

  options(warnings = -1)

  varnames <- sapply(substitute(list(...))[-1], deparse) %>%
    stringr::str_replace_all(pattern = '\"', replacement = "")

  sapply(
    X = varnames,
    FUN = function(x) {

      if (x %in% rownames(utils::installed.packages())) {

        if (!x %in% LoadedPackages()) {
          suppressPackageStartupMessages(suppressMessages(suppressWarnings(
            require(x, character.only = TRUE))))
          cat(crayon::red(">>"), "Library", crayon::blue(x), "was loaded successfully.\n")

        } else {

          if (Reload == TRUE) {
            suppressPackageStartupMessages(suppressMessages(suppressWarnings({
              ReloadPackage(Package = x)
              cat(crayon::red(">>"), "Library", crayon::blue(x), "was already loaded (re-loaded).\n")

            })))
          } else {
            cat(crayon::red(">>"), "Library", crayon::blue(x), "was already loaded (not re-loaded).\n")

          }
        }
      } else {

        if (InstallMissing) {
          suppressPackageStartupMessages(suppressMessages(suppressWarnings({
            utils::install.packages(
              x, quiet = TRUE, verbose = FALSE,
              repos = "https://cloud.r-project.org")
          })))

          if (x %in% rownames(utils::installed.packages())) {

            suppressPackageStartupMessages(suppressMessages(suppressWarnings({
              require(x, character.only = TRUE)
            })))
            cat(crayon::yellow(">>>>>>"), "Library", crayon::blue(x), "was installed and loaded. \n")
          } else {
            cat(crayon::yellow(">>>>>>"), "Library", crayon::blue(x), "cannot be installed. \n")
          }
        } else {
          cat(crayon::red(">>"), "Library", crayon::blue(x), "is not installed. ")
        }
      }
    })
  return(invisible(NULL))
}
