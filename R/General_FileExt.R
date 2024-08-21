## |------------------------------------------------------------------------| #
# FileExt ----
## |------------------------------------------------------------------------| #
#
#' Get the file extension from a file path
#'
#' This function extracts the file extension from a given file path. It first
#' checks if the input is not NULL and is a character string. If these
#' conditions are met, it then uses the [tools::file_ext] function to extract
#' and return the file extension. The function does not check the existence of
#' the file or explicitly get the file type from its content. It merely guess
#' file extension from the file name.
#' @name FileExt
#' @author Ahmed El-Gabbas
#' @return A character string representing the file extension of the given file
#'   path. If the path does not have an extension, an empty string is returned.
#' @param Path A character string representing the file path from which the file
#'   extension is to be extracted. It must not be `NULL` and has to be a
#'   character string.
#' @seealso [IASDT.R::FileType()]
#' @examples
#' FileExt(Path = "File.doc")
#'
#' FileExt(Path = "D:/File.doc")
#'
#' FileExt(Path = "File.1.doc")
#'
#' FileExt(Path = "D:/Files.All")
#'
#' FileExt(Path = "D:/Files.All/")
#'
#' FileExt("example.txt") # Returns "txt"
#'
#' FileExt("archive.tar.gz") # Returns "gz"
#' @export

FileExt <- function(Path) {
  if (is.null(Path)) {
    stop("Path cannot be NULL", .call = FALSE)
  }

  # Ensure Path is a character string
  if (!is.character(Path)) {
    stop("Path must be a character string", .call = FALSE)
  }

  return(tools::file_ext(Path))
}
