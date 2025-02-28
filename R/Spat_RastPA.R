## |------------------------------------------------------------------------| #
# RastPA ------
## |------------------------------------------------------------------------| #

#' Convert raster map into binary (1/0)
#'
#' This function converts raster values into a binary format where positive
#' values are set to 1 (presence) and zeros remain 0 (absence). Additionally, it
#' allows for the conversion of NA values to 0, and/or 0 values to NA, based on
#' the user's choice.
#' @name RastPA
#' @param x The input raster map. It must be of class `PackedSpatRaster`,
#'   `RasterLayer`, or `SpatRaster`. This parameter cannot be NULL.
#' @param NA_to_0 A logical value indicating whether NA values should be
#'   converted to 0. Defaults to `TRUE`.
#' @param Zero_to_NA A logical value indicating whether 0 values should be
#'   converted to NA. Defaults to `FALSE`.
#' @return A raster map where values have been converted according to the
#'   specified parameters. This object is of the same class as the input object.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(List = c("dplyr", "raster", "ggplot2", "tidyterra"))
#'
#' r <- raster::raster(system.file("external/test.grd", package = "raster"))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::rast(r), maxcell = Inf) +
#'   ggplot2::theme_minimal()
#'
#' R2 <- raster::stack(
#'   # NA replaced with 0
#'   RastPA(x = r),
#'   # NA is kept as NA
#'   RastPA(x = r, NA_to_0 = FALSE),
#'   # 0 replaced with NA in the second map
#'   RastPA(x = RastPA(r), Zero_to_NA = TRUE))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(
#'     data = terra::as.factor(terra::rast(R2)), maxcell = Inf) +
#'   ggplot2::facet_wrap(~lyr) +
#'   ggplot2::scale_fill_manual(values = c("grey30", "red"),
#'     na.value = "transparent") +
#'   ggplot2::theme_minimal()

RastPA <- function(x = NULL, NA_to_0 = TRUE, Zero_to_NA = FALSE) {

  if (is.null(x)) {
    stop("x can not be NULL", call. = FALSE)
  }

  if (inherits(x, "PackedSpatRaster")) {
    x <- terra::unwrap(x)
  }

  if (inherits(x, "RasterLayer")) {
    MaxVal <- raster::cellStats(x, max)
    if (MaxVal > 0) x[x > 0] <- 1
    if (NA_to_0) x <- raster::reclassify(x, cbind(NA, 0))
    if (Zero_to_NA) x <- raster::reclassify(x, cbind(0, NA))
  } else {

    if (!inherits(x, "SpatRaster")) {
      stop(
        "Input map should be either PackedSpatRaster, ",
        "RasterLayer, or SpatRaster", call. = FALSE)
    }

    x <- terra::classify(x, cbind(0, Inf, 1))
    if (NA_to_0) x <- terra::classify(x, cbind(NA, 0))
    if (Zero_to_NA) x <- terra::classify(x, cbind(0, NA))
  }
  return(x)
}
