## |------------------------------------------------------------------------| #
# LoadPackages ----
## |------------------------------------------------------------------------| #

#' Load or Install Multiple R Packages
#'
#' This function attempts to load multiple R packages specified by the user. If a package is not installed, the function can optionally install it before loading. It also provides an option to print the names and versions of the loaded packages.
#'
#' @param ... character. Names of the packages to be loaded or installed.
#' @param List character vector. An alternative or additional way to specify package names as a vector.
#' @param Verbose logical. If `TRUE`, prints the names and versions of the loaded packages. Defaults to `FALSE`.
#' @param InstallMissing logical. If `TRUE` (default), missing packages are automatically installed and then loaded.
#' @return This function is used for its side effects (loading/installing packages) and does not return any value.
#' @author Ahmed El-Gabbas
#' @export
#' @name LoadPackages
#' @examples
#' # Currently loaded packages
#' (P1 <- IASDT.R::LoadedPackages())
#'
#' LoadPackages(Hmsc)
#' # Loaded packages after implementing the function
#' (P2 <- IASDT.R::LoadedPackages())
#'
#' # Which packages were loaded?
#' setdiff(P2, P1)
#'
#' LoadPackages(terra, Verbose = TRUE)

LoadPackages <- function(..., List = NULL, Verbose = FALSE, InstallMissing = TRUE) {

  # Packages to load
  PG <- rlang::ensyms(...) %>%
    as.character() %>%
    c(List) %>%
    unique() %>%
    sort()

  # List of installed packages
  Pkg_Installed <- rownames(utils::installed.packages())

  # packages to install
  Pkg_ToInst <- setdiff(PG, Pkg_Installed)

  if (InstallMissing && length(Pkg_ToInst) > 0) {
    message("The following packages will be installed:\n", paste("  >>>>>  ", Pkg_ToInst, collapse = "\n"))

    # Installing missing packages
    purrr::walk(
      .x = Pkg_ToInst, .f = utils::install.packages,
      repos = "http://cran.us.r-project.org",
      dependencies = TRUE, quiet = TRUE) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()

  }  else if (length(Pkg_ToInst) > 0) {
    message("The following packages are neither available nor installed as InstallMissing = FALSE:\n", paste("  >>>>>  ", Pkg_ToInst, collapse = "\n"))
  }

  # Packages to load
  Pkg_ToLoad <- setdiff(PG, IASDT.R::LoadedPackages())

  purrr::walk(.x = Pkg_ToLoad, .f = library, character.only = TRUE,
              quietly = TRUE, warn.conflicts = FALSE) %>%
    invisible() %>%
    suppressWarnings() %>%
    suppressMessages()

  if (Verbose && length(Pkg_ToLoad) > 0) {

    loaded_packages <- purrr::map_chr(
      .x = Pkg_ToLoad,
      .f = function(pkg) {
        version <- utils::packageDescription(pkg)$Version
        paste("  >>>>>  ", pkg, ": ", as.character(version))
      })
    message("\nThe following packages were loaded:\n", paste(loaded_packages, collapse = "\n"))
  }
  return(invisible(NULL))
}
