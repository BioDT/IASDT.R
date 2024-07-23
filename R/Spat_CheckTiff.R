## |------------------------------------------------------------------------| #
# CheckTiff ----
## |------------------------------------------------------------------------| #

#' Check if a tiff file corrupted
#'
#' This function checks if the provided tiff file is corrupted by attempting to describe it using the `terra` package and searching for the presence of a "Driver" string in the description, which indicates a valid tiff file. If the string is found, the function returns `TRUE` and `FALSE` otherwise. 
#'
#' @param x Character; the file path of the tiff file to be checked. The function will stop with an error if `x` is `NULL` or if the file does not exist.
#' @name CheckTiff
#' @author Ahmed El-Gabbas
#' @return Logical; returns `TRUE` if the tiff file is not corrupted (i.e., it can be described and contains "Driver" in its description),
#' and `FALSE` otherwise.
#' @export
#' @note This function depends on the `terra`, `magrittr`, and `stringr` packages.
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' CheckTiff(x = f)

CheckTiff <- function(x = NULL) {
  
  # Check input argument
  if (is.null(x)) {
    stop("Input file cannot be NULL")
  }

  # Check if file exists
  if (magrittr::not(file.exists(x))) {
    stop("Input file does not exist")
  }

  x %>%
    terra::describe() %>%
    as.character() %>% 
    stringr::str_detect("Driver") %>%
    any() %>% 
    return()
}
