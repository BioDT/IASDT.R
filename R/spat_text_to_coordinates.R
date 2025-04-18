## |------------------------------------------------------------------------| #
# text_to_coordinates ----
## |------------------------------------------------------------------------| #

#' Extract longitude and latitude from string
#'
#' Extract longitude and latitude from string representing a geographical point
#' in the format `"POINT (longitude latitude)"` and converts it into a
#' two-column tibble containing the longitude and latitude as numeric values.
#' The names of the columns in the resulting tibble can be customized. The
#' default names for the longitude and latitude columns are "Longitude" and
#' "Latitude", respectively.
#' @param text Character. Coordinates in the format `"POINT (longitude
#'   latitude)"`. This parameter is required and cannot be `NULL`.
#' @param name_x,name_y Character. Name to be used for the longitude and
#'   Longitude columns in the output tibble. Defaults to "Longitude" and
#'   "Latitude".
#' @return A tibble with two columns containing the longitude and latitude
#'   values extracted from the input string. The names of these columns are
#'   determined by the `name_x` and `name_y` parameters. If no names are
#'   provided, the default names ("Longitude" and "Latitude") are used.
#' @name text_to_coordinates
#' @author Ahmed El-Gabbas
#' @return two column tibble for Longitude & Latitude
#' @export
#' @examples
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)",
#'   "POINT (16.179999 38.427214)") %>%
#'  lapply(text_to_coordinates)
#'
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)",
#'   "POINT (16.179999 38.427214)") %>%
#'  lapply(text_to_coordinates, name_x = "Long", name_y = "Lat")

text_to_coordinates <- function(
    text = NULL, name_x = "Longitude", name_y = "Latitude") {

  if (is.null(text)) {
    stop("Input string cannot be NULL", call. = FALSE)
  }

  text %>%
    # convert string to 2-columns data frame
    stringr::str_remove_all("POINT \\(|POINT\\(|\\)") %>%
    stringr::str_split_fixed(" ", 2) %>%
    as.numeric() %>%
    matrix(ncol = 2, byrow = FALSE) %>%
    as.data.frame() %>%
    stats::setNames(c(name_x, name_y)) %>%
    tibble::as_tibble()
}
