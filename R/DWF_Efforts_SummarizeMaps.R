## |------------------------------------------------------------------------| #
# Efforts_SummarizeMaps ----
## |------------------------------------------------------------------------| #

#' Summarize maps for efforts data
#'
#' This function processes spatial data (as an `sf` object), summarizes it based
#' on the number of observations or distinct species, and generates a raster
#' layer.
#' @param Data An `sf` object containing spatial data, with a column named
#'   `CellCode`.
#' @param NSp Logical. Whether to generate distinct species counts (`TRUE`) or 
#'   total observation counts (`FALSE`).
#' @param Name Character. Name of the count field and the prefix for the final 
#'   raster layer's name.
#' @param ClassOrder Character. The class and order combination (separated by 
#'   an underscore) represented in the `Data`.
#' @param Grid_SF,Grid_R Reference grid in the form of simple feature and
#'   raster.
#' @return A processed `terra` raster object representing the summarized data.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Efforts_Process] and [Efforts_Summarize]
#'   functions.
#' @author Ahmed El-Gabbas
#' @name Efforts_SummarizeMaps
#' @export

Efforts_SummarizeMaps <- function(
    Data, NSp, Name, ClassOrder, Grid_SF, Grid_R) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CellCode <- speciesKey <- NULL

  # # ..................................................................... ###

  # Validate if Data is an sf object
  if (!inherits(Data, "sf")) {
    stop(
      paste0(
        "Input data must be a simple feature (sf) object. ",
        "Provided data is of type: ",
        paste0(class(Data), collapse = "+")),
      call. = FALSE)
  }

  # Validate if NSp is logical
  if (!is.logical(NSp) || length(NSp) != 1) {
    stop(
      paste0(
        "The parameter `NSp` must be a single logical value (TRUE or FALSE). ",
        "Provided value is of type: ", paste0(class(NSp), collapse = "+")),
      call. = FALSE)
  }

  # Validate the Name parameter
  if (is.null(Name)) {
    stop("The parameter `Name` can not be empty", call. = FALSE)
  }

  # # ..................................................................... ###

  # Drop geometry from Data
  Data <- sf::st_drop_geometry(Data)

  # Generate distinct species counts if NSp is TRUE
  if (NSp) {
    Data <- dplyr::distinct(Data, CellCode, speciesKey)
  }

  # Count observations or species, join with the grid, and rasterize
  Data <- Data %>%
    dplyr::count(CellCode, name = Name) %>%
    dplyr::left_join(Grid_SF, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(Grid_R, field = Name) %>%
    terra::classify(cbind(NA, 0)) %>%
    terra::mask(Grid_R) %>%
    IASDT.R::setRastCRS() %>%
    IASDT.R::setRastVals() %>%
    stats::setNames(paste0(Name, "_", ClassOrder)) %>%
    terra::wrap()

  return(Data)
}
