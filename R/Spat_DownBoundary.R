## |------------------------------------------------------------------------| #
# DownBoundary ----
## |------------------------------------------------------------------------| #

#' Determine the boundaries of the requested GBIF data
#'
#' This function constructs a Well-Known Text (WKT) string representing a polygon that outlines the specified boundaries. It is used to define the area of interest for downloading GBIF data through the `rgbif::pred_within()` function.
#' @param L,R,B,T Numeric, the left, right, bottom, and top boundary of the area.
#' @name DownBoundary
#' @author Ahmed El-Gabbas
#' @return A character string representing the WKT of the polygon that outlines the specified boundaries.
#' @export
#' @description `rgbif::pred_within()` function used to download GBIF data only accepts a WKT string. This function takes the values of the boundary and converts it to a WKT string.
#' @examples
#' IASDT.R::DownBoundary(20, 30, 40, 50)

DownBoundary <- function(L = NULL, R = NULL, B = NULL, T = NULL) {

  if (any(c(is.null(L), is.null(R), is.null(B), is.null(T)))) {
    stop("none of L, R, B, or T can be NULL")
  }

  "POLYGON(({L} {B},{R} {B},{R} {T},{L} {T},{L} {B}))" %>%
    stringr::str_glue() %>%
    as.character() %>% 
    return()
}
