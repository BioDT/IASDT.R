## |------------------------------------------------------------------------| #
# boundary_to_WKT ----
## |------------------------------------------------------------------------| #

#' Determine the boundaries of the requested GBIF data
#'
#' This function constructs a Well-Known Text (WKT) string representing a
#' polygon that outlines the specified boundaries. It is used to define the area
#' of interest for downloading GBIF data through the [rgbif::pred_within()]
#' function.
#' @param left,right,bottom,top Numeric, the left, right, bottom, and top
#'   boundary of the area.
#' @name boundary_to_WKT
#' @author Ahmed El-Gabbas
#' @return A character string representing the WKT of the polygon that outlines
#'   the specified boundaries.
#' @export
#' @examples
#' IASDT.R::boundary_to_WKT(left = 20, right = 30, bottom = 40, top = 50)

boundary_to_WKT <- function(
  left = NULL, right = NULL, bottom = NULL, top = NULL) {

  if (any(c(is.null(left), is.null(right), is.null(bottom), is.null(top)))) {
    stop("none of left, right, bottom, or top can be NULL", call. = FALSE)
  }

  Bound <- stringr::str_glue(
    "POLYGON(({left} {bottom},{right} {bottom},{right} {top},",
    "{left} {top},{left} {bottom}))") %>%
    as.character()

  return(Bound)
}
