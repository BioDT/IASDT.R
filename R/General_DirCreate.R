## |------------------------------------------------------------------------| #
# DirCreate ----
## |------------------------------------------------------------------------| #
#
#' Create or verify the existence of a directory
#'
#' This function creates a directory at the specified path. If the directory already exists, it can optionally print a message indicating so. The creation is recursive, meaning it will create any necessary parent directories.
#'
#' @param Path character; the path where the directory is to be created. Cannot be NULL.
#' @param Verbose logical; if `TRUE`, the function prints messages about the operation's outcome. Defaults to `TRUE`.
#' @name DirCreate
#' @author Ahmed El-Gabbas
#' @return NULL invisibly. The function is used for its side effect (creating a directory or printing a message) rather than any return value.
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
#' @references This function is currently not exported. I use `fs::dir_create()` instead.


DirCreate <- function(Path, Verbose = TRUE) {
  
  if (is.null(Path)) {
    stop("Path cannot be NULL")
  }

  Path2 <- gsub("\\\\", "/", Path)

  if (dir.exists(Path) && Verbose) {
    IASDT.R::CatTime(
      stringr::str_glue("Path: {crayon::bold(Path2)} - already exists"),
      Date = TRUE)
  } else {
    dir.create(Path, recursive = TRUE, showWarnings = FALSE)
    if (Verbose) {
      "Path: {crayon::bold(Path2)} created" %>%
        stringr::str_glue() %>% 
        IASDT.R::CatTime(Date = TRUE)
    }
  }
  return(invisible(NULL))
}
