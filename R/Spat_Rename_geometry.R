## |------------------------------------------------------------------------| #
# Rename_geometry ----
## |------------------------------------------------------------------------| #

#' Rename active geometry column of an `sf` object.
#'
#' This function renames the active geometry column of a simple feature (`sf`) object to a new name provided by the user.
#'
#' @param g `sf` object; the simple feature object whose geometry column name is to be changed.
#' @param name `character`; the new name for the geometry column.
#' @name Rename_geometry
#' @author Ahmed El-Gabbas
#' @return The modified `sf` object with the renamed geometry column.
#' @references [Click here](https://gis.stackexchange.com/a/386589/30390)
#' @export

Rename_geometry <- function(g = NULL, name = NULL) {

  if (any(is.null(g) | is.null(name))) {
    stop("The input sf object or name cannot be 'NULL'.")
  }

  current <- attr(g, "sf_column")
  names(g)[names(g) == current] <- name
  sf::st_geometry(g) <- name
  return(g)
}
