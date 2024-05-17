# |---------------------------------------------------| #
# GetCV ----
# |---------------------------------------------------| #

#' Prepare cross-validation data
#'
#' Prepare cross-validation data
#' @param DT Input data set
#' @param NR Integer. Number of rows
#' @param NC Integer. Number of columns
#' @param Path_Grid String. Path of the reference grid.
#' @name GetCV
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetCV <- function(DT, NR = 4, NC = 3, Path_Grid) {
  Grid10 <- Path_Grid %>%
    file.path("Grid_10_Land_Crop.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap()

  DT_R <- DT %>%
    dplyr::select(x, y) %>%
    as.matrix() %>%
    terra::rasterize(Grid10) %>%
    terra::trim()

  DT_SF <- sf::st_as_sf(DT, coords = c("x", "y"), crs = 3035, remove = FALSE)
  CV <- blockCV::cv_spatial(
    x = DT_SF, r = DT_R, hexagon = FALSE, rows_cols = c(NR, NC),
    plot = FALSE, progress = FALSE, report = FALSE)

  sf::st_join(DT_SF, dplyr::select(CV$blocks, CV = folds), largest = TRUE) %>%
    sf::st_drop_geometry() %>%
    suppressWarnings()
}
