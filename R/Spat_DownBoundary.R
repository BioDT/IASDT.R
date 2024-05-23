## |------------------------------------------------------------------------| #
# DownBoundary ----
## |------------------------------------------------------------------------| #

#' Determine the boundaries of the requested GBIF data
#'
#' Determine the boundaries of the requested GBIF data
#' @param L left boundary
#' @param R right boundary
#' @param B bottom boundary
#' @param T top boundary
#' @name DownBoundary
#' @author Ahmed El-Gabbas
#' @return WKT string
#' @export
#' @description `rgbif::pred_within()` function used to download GBIF data only accepts a WKT string. This function takes the values of the boundary and converts it to a WKT string. Default values are determined by the variables: Bound_L, R = Bound_R, Bound_B, Bound_T...
#' @examples
#' IASDT.R::DownBoundary(20, 30, 40, 50)

DownBoundary <- function(L, R, B, T) {
  "POLYGON(({L} {B},{R} {B},{R} {T},{L} {T},{L} {B}))" %>%
    stringr::str_glue() %>%
    as.character()
}
