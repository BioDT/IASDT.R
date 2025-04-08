## |------------------------------------------------------------------------| #
# raster_to_PA ------
## |------------------------------------------------------------------------| #

#' Convert raster map into binary (1/0)
#'
#' This function converts raster values into a binary format where positive
#' values are set to 1 (presence) and zeros remain 0 (absence). Additionally, it
#' allows for the conversion of NA values to 0, and/or 0 values to NA, based on
#' the user's choice.
#' @name raster_to_PA
#' @param raster The input raster map. It must be of class `PackedSpatRaster`,
#'   `RasterLayer`, or `SpatRaster`. This parameter cannot be NULL.
#' @param NA_to_0 A logical value indicating whether NA values should be
#'   converted to 0. Defaults to `TRUE`.
#' @param zero_to_NA A logical value indicating whether 0 values should be
#'   converted to NA. Defaults to `FALSE`.
#' @return A raster map where values have been converted according to the
#'   specified parameters. This object is of the same class as the input object.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::load_packages(
#'    package_list = c("dplyr", "raster", "ggplot2", "tidyterra"))
#'
#' r <- raster::raster(system.file("external/test.grd", package = "raster"))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::rast(r), maxcell = Inf) +
#'   ggplot2::theme_minimal()
#'
#' R2 <- raster::stack(
#'   # NA replaced with 0
#'   raster_to_PA(raster = r),
#'   # NA is kept as NA
#'   raster_to_PA(raster = r, NA_to_0 = FALSE),
#'   # 0 replaced with NA in the second map
#'   raster_to_PA(raster = raster_to_PA(r), zero_to_NA = TRUE))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(
#'     data = terra::as.factor(terra::rast(R2)), maxcell = Inf) +
#'   ggplot2::facet_wrap(~lyr) +
#'   ggplot2::scale_fill_manual(values = c("grey30", "red"),
#'     na.value = "transparent") +
#'   ggplot2::theme_minimal()

raster_to_PA <- function(raster = NULL, NA_to_0 = TRUE, zero_to_NA = FALSE) {

  if (is.null(raster)) {
    stop("raster can not be NULL", call. = FALSE)
  }

  if (inherits(raster, "PackedSpatRaster")) {
    raster <- terra::unwrap(raster)
  }

  if (inherits(raster, "RasterLayer")) {
    MaxVal <- raster::cellStats(raster, max)
    if (MaxVal > 0) raster[raster > 0] <- 1
    if (NA_to_0) raster <- raster::reclassify(raster, cbind(NA, 0))
    if (zero_to_NA) raster <- raster::reclassify(raster, cbind(0, NA))
  } else {

    if (!inherits(raster, "SpatRaster")) {
      stop(
        "Input map should be either PackedSpatRaster, ",
        "RasterLayer, or SpatRaster", call. = FALSE)
    }

    raster <- terra::classify(raster, cbind(0, Inf, 1))
    if (NA_to_0) raster <- terra::classify(raster, cbind(NA, 0))
    if (zero_to_NA) raster <- terra::classify(raster, cbind(0, NA))
  }
  return(raster)
}
