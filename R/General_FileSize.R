# |---------------------------------------------------| #
# FileSize ----
# |---------------------------------------------------| #

#' File size in a Human Readable format
#'
#' File size in a Human Readable format
#'
#' @param File character; file path
#' @param ... additional arguments for `gdata::humanReadable` function
#' @name FileSize
#' @export

FileSize <- function(File, ...) {
  file.size(File) %>%
    gdata::humanReadable()
}
