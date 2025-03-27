## |------------------------------------------------------------------------| #
# sf_add_coords ------
## |------------------------------------------------------------------------| #

#' Add longitude and latitude coordinates to an sf object
#'
#' Add longitude and latitude coordinates as new columns to an sf object
#' (`Sf_Obj`). It extracts the coordinates from the sf object, converts them
#' into a tibble, and appends them to the original sf object as new columns. If
#' `NameX` or `NameY`, provided as arguments respectively, already exist in the
#' sf object, the function either  1) overwrites these columns if `Overwrite` is
#' set to `TRUE` or 2) appends "_NEW" to the new column names to avoid overwrite
#' if `Overwrite` is set to `FALSE`.
#' @name sf_add_coords
#' @param Sf_Obj An `sf` object to which longitude and latitude columns will be
#'   added.
#' @param NameX,NameY Character. Name of the longitude column to be added.
#'   Defaults to `Long` and `Lat`.
#' @param Overwrite Logical. Whether to overwrite existing columns with names
#'   specified by `NameX` and `NameY`. If `FALSE` and columns with these names
#'   exist, new columns are appended with "_NEW" suffix. Defaults to `FALSE`.
#' @return An sf object with added longitude and latitude columns.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' pt1 = sf::st_point(c(0,1))
#' pt2 = sf::st_point(c(1,1))
#' d = data.frame(a = c(1, 2))
#' d$geom = sf::st_sfc(pt1, pt2)
#' df = sf::st_as_sf(d)
#' df
#' (df <- sf_add_coords(df))
#'
#' (sf_add_coords(df))
#'
#' (sf_add_coords(df, Overwrite = TRUE))
#' @note If the Overwrite parameter is `FALSE` (default) and columns with the
#'   specified names already exist, the function will issue a warning and append
#'   "_NEW" to the names of the new columns to avoid overwriting.

sf_add_coords <- function(
    Sf_Obj, NameX = "Long", NameY = "Lat", Overwrite = FALSE) {

  if (!inherits(Sf_Obj, "sf")) {
    stop("`Sf_Obj` must be an sf object", call. = FALSE)
  }

  if (!is.character(NameX) || !is.character(NameY) ||
      !nzchar(NameX) || !nzchar(NameY)) {
    stop("`NameX` and `NameY` must be non-empty strings", call. = FALSE)
  }

  ColNames <- names(Sf_Obj)

  # Coordinate Extraction
  # extract the coordinates from the sf object and converts them into a tibble,
  # naming the columns according to NameX and NameY
  Coords <- sf::st_coordinates(Sf_Obj) %>%
    tibble::as_tibble() %>%
    stats::setNames(c(NameX, NameY))

  # Column Name Check
  # Before adding the new columns, check if columns with the
  # names NameX and NameY already exist. If they do, it either 1) Overwrites
  # these columns if Overwrite is TRUE, after issuing a warning. 2) Appends
  # "_NEW" to the new column names to avoid overwriting, if Overwrite is FALSE.

  if (any(c(NameX, NameY) %in% ColNames)) {
    if (Overwrite) {
      warning(
        "Provided column names for longitude and Latitude ",
        "already exist in the data; these columns were overwritten",
        call. = FALSE)
      Sf_Obj <- dplyr::select(Sf_Obj, -dplyr::all_of(c(NameX, NameY)))
    } else {
      warning(
        "Provided column names for longitude and Latitude already exist ",
        "in the data; `_NEW` is used as suffix", call. = FALSE)
      Coords <- Coords %>%
        stats::setNames(c(paste0(NameX, "_NEW"), paste0(NameY, "_NEW")))
    }
  } else {
    Coords <- stats::setNames(Coords, c(NameX, NameY))
  }
  return(dplyr::bind_cols(Sf_Obj, Coords))
}
