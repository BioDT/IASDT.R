# |---------------------------------------------------| #
# sf_add_coords ------
# |---------------------------------------------------| #

#' Add coordinates of the current sf object as columns
#'
#' Add coordinates of the current sf object as columns
#' @name sf_add_coords
#' @param Sf_Obj input sf object
#' @param NameX Name of the longitude column to be added
#' @param NameY Name of the latitude column to be added
#' @param Overwrite Should columns `NameX` or `NameY` be overwritten if these columns exist in the original data
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' pt1 = sf::st_point(c(0,1))
#' pt2 = sf::st_point(c(1,1))
#' d = data.frame(a = 1:2)
#' d$geom = sf::st_sfc(pt1, pt2)
#' df = sf::st_as_sf(d)
#' df
#' (df <- sf_add_coords(df))
#'
#' (sf_add_coords(df))
#'
#' (sf_add_coords(df, Overwrite = TRUE))

sf_add_coords <- function(Sf_Obj, NameX = "Long", NameY = "Lat", Overwrite = FALSE) {
  ColNames <- names(Sf_Obj)
  Coords <- Sf_Obj %>%
    sf::st_coordinates() %>%
    tibble::as_tibble() %>%
    stats::setNames(c(NameX, NameY))

  if (any(NameX %in% ColNames || NameY %in% ColNames)) {
    if (Overwrite) {
      warning("Provided column names for longitude and Latitude already exist in the data; these columns were overwritten")
      Sf_Obj <- dplyr::select(Sf_Obj, -dplyr::all_of(c(NameX, NameY)))
    } else {
      warning('Provided column names for longitude and Latitude already exist in the data; "_NEW" is used as suffix')
      Coords <- Coords %>%
        stats::setNames(c(paste0(NameX, "_NEW"), paste0(NameY, "_NEW")))
    }
  } else {
    Coords <- stats::setNames(Coords, c(NameX, NameY))
  }
  return(dplyr::bind_cols(Sf_Obj, Coords))
}
