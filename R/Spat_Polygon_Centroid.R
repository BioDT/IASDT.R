## |------------------------------------------------------------------------| #
# Polygon_Centroid ----
## |------------------------------------------------------------------------| #

#' Replace the geometry of a polygon with its centroid point
#'
#' This function replaces the geometry of a simple feature (`sf`) polygon object
#' with the geometry of its centroid point. It can optionally rename the
#' geometry column of the modified `sf` object.
#' @param x A simple feature (`sf`) object; the polygon whose geometry is to be
#'   replaced with its centroid. Cannot be `NULL`.
#' @param Rename Logical. Whether to rename the geometry column of the sf
#'   object. Defaults to `FALSE`.
#' @param NewName Character. New name for the geometry column. Only valid if
#'   `Rename = TRUE`.
#' @name Polygon_Centroid
#' @author Ahmed El-Gabbas
#' @references [Click here](https://github.com/r-spatial/sf/issues/480)
#' @return The modified sf object with its geometry replaced by the centroid of
#'   the original polygon geometry. If Rename is `TRUE`, the geometry column
#'   will also be renamed as specified by NewName.
#' @export

Polygon_Centroid <- function(x = NULL, Rename = FALSE, NewName = NULL) {

  if (is.null(x)) {
    stop("Input sf object cannot be NULL", call. = FALSE)
  }

  suppressWarnings(sf::st_geometry(x) <- sf::st_geometry(sf::st_centroid(x)))

  if (Rename) {
    if (is.null(NewName)) {
      stop("NewName cannot be NULL", call. = FALSE)
    }

    x <- Rename_geometry(g = x, name = NewName)
  }

  return(x)
}
