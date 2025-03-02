## |------------------------------------------------------------------------| #
# CheckTiff ----
## |------------------------------------------------------------------------| #

#' Check if a tiff file corrupted
#'
#' This function checks if the provided tiff file is corrupted by attempting to
#' describe it using the `terra` package and searching for the presence of a
#' "Driver" string in the description, which indicates a valid tiff file. If the
#' string is found, the function returns `TRUE` and `FALSE` otherwise. The
#' function works also for reading netCDF files with the `terra` package.
#' @param x Character; the file path of the tiff file to be checked. The
#'   function will stop with an error if `x` is `NULL` or if the file does not
#'   exist.
#' @name CheckTiff
#' @author Ahmed El-Gabbas
#' @return Logical; returns `TRUE` if the tiff file is not corrupted (i.e., it
#'   can be described and contains "Driver" in its description), and `FALSE`
#'   otherwise.
#' @export
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' CheckTiff(x = f)

CheckTiff <- function(x = NULL) {

  # Check input argument
  if (is.null(x)) {
    stop("Input file cannot be NULL", call. = FALSE)
  }

  # # ..................................................................... ###

  # Check if file exists
  if (!file.exists(x)) {
    warning("Input file does not exist", call. = FALSE)
    return(FALSE)
  }

  # # ..................................................................... ###

  # Check file metadata using terra's describe
  MetadataOkay <- as.character(terra::describe(x = x)) %>%
    stringr::str_detect("Driver") %>%
    any()

  if (isFALSE(MetadataOkay)) {
    return(FALSE)
  }

  # # ..................................................................... ###

  return(terra::hasValues(terra::rast(x)))
}
