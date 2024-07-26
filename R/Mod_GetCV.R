## |------------------------------------------------------------------------| #
# GetCV ----
## |------------------------------------------------------------------------| #

#' Prepare cross-validation data for spatial analysis
#'
#' This function prepares cross-validation data by dividing the input dataset into spatial blocks. The function uses [blockCV::cv_spatial] function to assign grid cells to cross-validation folds. It is designed to work with spatial data, specifically for scenarios where spatial autocorrelation might influence the cross-validation process.
#' The function may need some improvement to implement more than cross-validation strategy.
#'
#' @param DT A data frame or tibble containing the input dataset with at least two columns for x and y coordinates.
#' @param NR,NC Integer, the number of rows and columns to divide the spatial area into. Defaults to 4 row and 3 columns.
#' @param Path_Grid String, the file path to the reference grid used for rasterizing the input data.
#' @name GetCV
#' @author Ahmed El-Gabbas
#' @return The function returns a modified version of the input dataset `DT` with an additional column indicating the cross-validation fold each record belongs to. This allows for spatially-aware cross-validation processes.
#' @export

GetCV <- function(DT, NR = 4, NC = 3, Path_Grid) {

  if (is.null(DT) || is.null(Path_Grid)) {
    stop("DT and Path_Grid can not be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  folds <- NULL

  Grid10 <- file.path(Path_Grid, "Grid_10_Land_Crop.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap()

  DT_R <- DT %>%
    dplyr::select("x", "y") %>%
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
