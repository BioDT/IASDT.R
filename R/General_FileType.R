## |------------------------------------------------------------------------| #
# FileType ----
## |------------------------------------------------------------------------| #
#
#' Determine the file type of a given file path
#'
#' This function uses the system's `file` command to determine the type of the
#' file specified by the `Path` parameter. It returns a character string
#' describing the file type.
#' @param Path Character. The path to the file whose type is to be determined.
#'   The path must not be `NULL`, and the file must exist.
#' @return A character string describing the file type.
#' @name FileType
#' @export
#' @note This function relies on the system's `file` command and therefore
#'    might produce different outputs on different platforms.
#' @author Ahmed El-Gabbas
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' FileType(Path = f)

FileType <- function(Path) {

  # Check `file` system command
  if (isFALSE(IASDT.R::CheckCommands("file"))) {
    stop("The system command 'file' is not available", call. = FALSE)
  }

  if (is.null(Path)) {
    stop("Path cannot be NULL", call. = FALSE)
  }

  # Ensure Path is a character string
  if (!is.character(Path)) {
    stop("Path must be a character string", call. = FALSE)
  }

  # Ensure file exists
  if (!file.exists(Path)) {
    stop("File does not exist", call. = FALSE)
  }

  Out <- paste0("file ", IASDT.R::NormalizePath(Path)) %>%
    system(intern = TRUE) %>%
    stringr::str_extract_all(": .+", simplify = TRUE) %>%
    as.vector() %>%
    stringr::str_remove("^: ")

  return(Out)
}
