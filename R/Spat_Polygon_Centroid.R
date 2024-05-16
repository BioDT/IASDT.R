# |---------------------------------------------------| #
# Polygon_Centroid ----
# |---------------------------------------------------| #

#' Replace the geometry of a polygon with its centroid point
#'
#' Replace the geometry of a polygon with its centroid point
#' @param x simple feature object
#' @param Rename should the geometry field renamed
#' @param NewName the new name of geometry
#' @name Polygon_Centroid
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Polygon_Centroid <- function(x, Rename = FALSE, NewName = "") {
  # https://github.com/r-spatial/sf/issues/480
  suppressWarnings(sf::st_geometry(x) <- sf::st_geometry(sf::st_centroid(x)))
  if (Rename) {
    x %>%
      Rename_geometry(name = NewName) %>%
      return()
  } else {
    return(x)
  }
}
