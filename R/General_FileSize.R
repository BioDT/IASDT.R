## |------------------------------------------------------------------------| #
# FileSize ----
## |------------------------------------------------------------------------| #

#' File size in a Human Readable format
#'
#' This function takes the path to a file and returns its size in a format that
#' is easy to read (e.g., KB, MB, GB), using the [gdata::humanReadable]
#' function.
#' @param File character; the path to the file whose size you want to check.
#' @param ... additional arguments passed to the [gdata::humanReadable]
#'   function, allowing for further customization of the output format.
#' @return A character string representing the size of the file in a
#'   human-readable format.
#' @name FileSize
#' @export
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' FileSize(File = f)

FileSize <- function(File, ...) {

  if (is.null(File)) {
    stop("File cannot be NULL", call. = FALSE)
  }

  return(gdata::humanReadable(fs::file_size(File), ...))
}
