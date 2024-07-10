## |------------------------------------------------------------------------| #
# Text2Coords ----
## |------------------------------------------------------------------------| #

#' Extract longitude / latitude from text
#'
#' Extract longitude / latitude from text
#'
#' @param String string of coordinate
#' @param Long_Name Column name for the longitude information; default: "Longitude"
#' @param Lat_Name Column name for the latitude information; default: "Latitude"
#' @name Text2Coords
#' @author Ahmed El-Gabbas
#' @return two column tibble for Longitude & Latitude
#' @export
#' @examples
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)", "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords)
#'
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)", "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords, Long_Name = "Long", Lat_Name = "Lat")

Text2Coords <- function(String, Long_Name = "Longitude", Lat_Name = "Latitude") {
  String %>%
    # convert string to 2-columns data frame
    stringr::str_remove_all("POINT \\(|\\)") %>%
    stringr::str_split(" ", simplify = TRUE) %>%
    as.numeric() %>%
    matrix(nrow = 1, ncol = 2) %>%
    as.data.frame() %>%
    stats::setNames(c(Long_Name, Lat_Name)) %>%
    tibble::tibble()
}
