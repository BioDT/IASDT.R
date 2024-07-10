## |------------------------------------------------------------------------| #
# LoadPackages ----
## |------------------------------------------------------------------------| #

#' Load package silently (+ install missing packages)
#'
#' Load package silently (+ install missing packages)
#'
#' @param ... packages to load / install
#' @param List vector for the name of packages to be loaded
#' @param Verbose print a message of the package name and version
#' @name LoadPackages
#' @author Ahmed El-Gabbas
#' @return NULL
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
    sort() %>%
    setdiff(as.character(.packages()))

  # packages to install
  InstPack <- utils::installed.packages() %>%
    as.data.frame() %>%
    "["("Package") %>%
    unlist() %>%
    setdiff(x = PG)

  if (length(InstPack) > 0) {
    cat("The following packages will be installed\n")
    InstPack %>%
      stringr::str_c("  >>>>>  ", .) %>%
      stringr::str_c("", collapse = "\n") %>%
      cat()

    InstPack %>%
      purrr::map(
        utils::install.packages, repos = "http://cran.us.r-project.org",
        dependencies = TRUE, quiet = TRUE) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()

    cat("\n")
  }

  # load packages
  PG %>%
    sapply(library, character.only = TRUE, quietly = TRUE, verbose = FALSE) %>%
    invisible() %>%
    suppressWarnings() %>%
    suppressMessages()

  if (Verbose && length(PG) > 0) {
    cat("\nThe following packages were loaded:\n")
    PG %>%
      purrr::map_chr(~{
        .x %>%
          utils::packageDescription() %>%
          "$"("Version") %>%
          as.character() %>%
          paste0("  >>>>>  ", .x, ": ", .)
      }) %>%
      stringr::str_c(collapse = "\n") %>%
      cat()
    cat("\n")
  }
  return(invisible(NULL))
}
