## |------------------------------------------------------------------------| #
# Rename_geometry ----
## |------------------------------------------------------------------------| #

#' Rename active geometry column of an sf object
#'
#' Rename active geometry column of an sf object
#'
#' @param g simple feature object
#' @param name the new name of geometry
#' @name Rename_geometry
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Rename_geometry <- function(g, name) {
  # Source: https://gis.stackexchange.com/a/386589/30390
  # https://gis.stackexchange.com/questions/386584/sf-geometry-column-naming-differences-r
  current <- attr(g, "sf_column")
  names(g)[names(g) == current] <- name
  sf::st_geometry(g) <- name
  return(g)
}
