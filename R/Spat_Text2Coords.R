## |------------------------------------------------------------------------| #
# Text2Coords ----
## |------------------------------------------------------------------------| #

#' Extract longitude and latitude from string
#'
#' Extract longitude and latitude from string representing a geographical point
#' in the format `"POINT (longitude latitude)"` and converts it into a
#' two-column tibble containing the longitude and latitude as numeric values.
#' The names of the columns in the resulting tibble can be customized. The
#' default names for the longitude and latitude columns are "Longitude" and
#' "Latitude", respectively.
#' @param String Character. Coordinates in the format `"POINT (longitude
#'   latitude)"`. This parameter is required and cannot be `NULL`.
#' @param Long_Name,Lat_Name Character. Name to be used for the longitude and
#'   Longitude columns in the output tibble. Defaults to "Longitude" and
#'   "Latitude".
#' @return A tibble with two columns containing the longitude and latitude
#'   values extracted from the input string. The names of these columns are
#'   determined by the `Long_Name` and `Lat_Name` parameters. If no names are
#'   provided, the default names ("Longitude" and "Latitude") are used.
#' @name Text2Coords
#' @author Ahmed El-Gabbas
#' @return two column tibble for Longitude & Latitude
#' @export
#' @examples
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)",
#'   "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords)
#'
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)",
#'   "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords, Long_Name = "Long", Lat_Name = "Lat")

Text2Coords <- function(
    String = NULL, Long_Name = "Longitude", Lat_Name = "Latitude") {

  if (is.null(String)) {
    stop("Input string cannot be NULL", call. = FALSE)
  }

  String %>%
    # convert string to 2-columns data frame
    stringr::str_remove_all("POINT \\(|\\)") %>%
    stringr::str_split_fixed(" ", 2) %>%
    as.numeric() %>%
    matrix(ncol = 2, byrow = FALSE) %>%
    as.data.frame() %>%
    stats::setNames(c(Long_Name, Lat_Name)) %>%
    tibble::as_tibble()
}
