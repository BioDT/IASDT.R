## |------------------------------------------------------------------------| #
# RequireMultiple ----
## |------------------------------------------------------------------------| #

#' Load (or install) multiple packages at once
#'
#' Load (or install) multiple packages at once
#'
#' @param ... Packages to load or install
#' @param Reload Should the already loaded packages re-loaded? Default: `FALSE`
#' @param InstallMissing Should the missing packages be installed and loaded? Default: `TRUE`
#' @author Ahmed El-Gabbas
#' @export
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
    stringr::str_replace_all(pattern = '\"', "")

  sapply(varnames, function(x) {

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
          cat(crayon::yellow(">>>>>>"), "Library", crayon::blue(x), "can not be installed. \n")
        }
      } else {
        cat(crayon::red(">>"), "Library", crayon::blue(x), "is not installed. ")
      }
    }
  })
  return(invisible(NULL))
}
