## |------------------------------------------------------------------------| #
# DownBoundary ----
## |------------------------------------------------------------------------| #

#' Determine the boundaries of the requested GBIF data
#'
#' This function constructs a Well-Known Text (WKT) string representing a
#' polygon that outlines the specified boundaries. It is used to define the area
#' of interest for downloading GBIF data through the [rgbif::pred_within()]
#' function.
#' @param Left,Right,Bottom,Top Numeric, the left, right, bottom, and top
#'   boundary of the area.
#' @name DownBoundary
#' @author Ahmed El-Gabbas
#' @return A character string representing the WKT of the polygon that outlines
#'   the specified boundaries.
#' @export
#' @examples
#' IASDT.R::DownBoundary(Left = 20, Right = 30, Bottom = 40, Top = 50)

DownBoundary <- function(Left = NULL, Right = NULL, Bottom = NULL, Top = NULL) {

  if (any(c(is.null(Left), is.null(Right), is.null(Bottom), is.null(Top)))) {
    stop("none of Left, Right, Bottom, or Top can be NULL", call. = FALSE)
  }

  "POLYGON(({Left} {Bottom},{Right} {Bottom},{Right} {Top},{Left} {Top},{Left} {Bottom}))" %>%
    stringr::str_glue() %>%
    as.character() %>%
    return()
}
