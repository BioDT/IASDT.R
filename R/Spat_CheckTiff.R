## |------------------------------------------------------------------------| #
# CheckTiff ----
## |------------------------------------------------------------------------| #

#' Check if tiff file corrupted
#'
#' Check if tiff file corrupted
#' @param x tiff file path
#' @name CheckTiff
#' @author Ahmed El-Gabbas
#' @return logical: `TRUE` when the tiff file is okay; `FALSE` when corrupted
#' @export
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' CheckTiff(f)

CheckTiff <- function(x) {
  x %>%
    terra::describe() %>%
    stringr::str_detect("Driver") %>%
    any()
}
