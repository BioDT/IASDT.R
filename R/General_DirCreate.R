## |------------------------------------------------------------------------| #
# DirCreate ----
## |------------------------------------------------------------------------| #
#
#' Create directory if not existed (recursively)
#'
#' Create directory if not existed (recursively)
#' @param Path character; folder path
#' @param Verbose logical; print a message of whether the folder was created or already available. Default: `TRUE`
#' @name DirCreate
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' # create a new folder (random name) in the temporary folder
#' Path2Create <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' file.exists(Path2Create)
#'
#' DirCreate(Path2Create)
#' DirCreate(Path2Create)
#' DirCreate(Path2Create, Verbose = FALSE)
#' file.exists(Path2Create)
#'
#' @noRd

DirCreate <- function(Path, Verbose = TRUE) {
  Path2 <- gsub("\\\\", "/", Path)
  if (dir.exists(Path) && Verbose) {
    CatTime(stringr::str_glue("Path: {crayon::bold(Path2)} - already exists"), Date = TRUE)
  } else {
    dir.create(Path, recursive = TRUE, showWarnings = FALSE)
    if (Verbose) {
      "Path: {crayon::bold(Path2)} created" %>%
        stringr::str_glue() %>%
        CatTime(Date = TRUE)
    }
  }
}
