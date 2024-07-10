## |------------------------------------------------------------------------| #
# Set_geometry ----
## |------------------------------------------------------------------------| #

#' Set geometry of an sf object in the pipe line
#'
#' Set geometry of an sf object in the pipe line
#'
#' @param x simple feature data frame
#' @param Name the name of the geometry column to be used
#' @name Set_geometry
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Set_geometry <- function(x, Name) {
  sf::st_geometry(x) <- Name
  return(x)
}
