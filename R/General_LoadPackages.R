## |------------------------------------------------------------------------| #
# LoadPackages ----
## |------------------------------------------------------------------------| #

#' Load package silently and install missing packages
#'
#' This function attempts to load specified R packages. If any package is not installed, it installs them first. It can operate silently or verbosely, indicating which packages were loaded and their versions.
#'
#' @param ... character. Names of packages to load/install, passed as individual arguments.
#' @param List character vector. An alternative or additional way to specify package names as a vector.
#' @param Verbose logical. If `TRUE`, prints the names and versions of loaded packages. Defaults to `FALSE`.
#' @name LoadPackages
#' @author Ahmed El-Gabbas
#' @return NULL This function does not return anything.
#' @export
#' @examples
#' LoadPackages(raster)
#'
#' LoadPackages(terra, Verbose = TRUE)

LoadPackages <- function(..., List = NULL, Verbose = FALSE) {

  # Packages to load
  PG <- rlang::ensyms(...) %>%
    as.character() %>%
    c(List) %>%
    unique() %>%
    sort() %>%
    setdiff(as.character(.packages(all.available = TRUE)))

  # packages to install
  InstPack <- setdiff(PG, rownames(utils::installed.packages()))

  if (length(InstPack) > 0) {
    cat("The following packages will be installed\n")
    InstPack %>%
      stringr::str_c("  >>>>>  ", .) %>%
      stringr::str_c("", collapse = "\n") %>%
      cat()

    # Installing missing packages
    purrr::map(
      .x = InstPack, .f = utils::install.packages,
      repos = "http://cran.us.r-project.org",
      dependencies = TRUE, quiet = TRUE) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()

    cat("\n")
  }

  # load packages
  sapply(PG, library, character.only = TRUE, quietly = TRUE, verbose = FALSE) %>%
    invisible() %>%
    suppressWarnings() %>%
    suppressMessages()

  if (Verbose && length(PG) > 0) {
    cat("\nThe following packages were loaded:\n")

    purrr::map_chr(
      .x = PG,
      .f = ~{
        .x %>%
          utils::packageDescription() %>%
          magrittr::extract2("Version") %>%
          as.character() %>%
          paste0("  >>>>>  ", .x, ": ", .)
      }) %>%
      stringr::str_c(collapse = "\n") %>%
      cat()
    cat("\n")
  }
  return(invisible(NULL))
}
