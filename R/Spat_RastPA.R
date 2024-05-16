# |---------------------------------------------------| #
# RastPA ------
# |---------------------------------------------------| #

#' convert raster map into binary (1/0)
#'
#' convert raster map into binary (1/0)
#'
#' @name RastPA
#' @param x input map; either `PackedSpatRaster`, `RasterLayer`, or `SpatRaster`
#' @param NA_to_0 logical (default: `TRUE`); should NA be replaced with 0
#' @param Zero_to_NA logical (default: `FALSE`); should 0 be replaced with NA
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(dplyr, raster, ggplot2, tidyterra)
#'
#' r <- raster::raster(system.file("external/test.grd", package = "raster"))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::rast(r), maxcell = Inf) +
#'   ggplot2::theme_minimal()
#'
#' R2 <- raster::stack(
#'   RastPA(r),                            # NA replaced with 0
#'   RastPA(r, NA_to_0 = FALSE),           # NA is kept as NA
#'   RastPA(RastPA(r), Zero_to_NA = TRUE)) # 0 replaced with NA in the second map
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::as.factor(terra::rast(R2)), maxcell = Inf) +
#'   ggplot2::facet_wrap(~lyr) +
#'   ggplot2::scale_fill_manual(values = c("grey30", "red"), na.value = "transparent") +
#'   ggplot2::theme_minimal()

RastPA <- function(x, NA_to_0 = TRUE, Zero_to_NA = FALSE) {
  if (inherits(x, "PackedSpatRaster")) x <- terra::unwrap(x)

  if (inherits(x, "RasterLayer")) {
    MaxVal <- raster::cellStats(x, max)
    if (MaxVal > 0) x[x > 0] <- 1
    if (NA_to_0) x <- raster::reclassify(x, cbind(NA, 0))
    if (Zero_to_NA) x <- raster::reclassify(x, cbind(0, NA))
  } else {
    if (inherits(x, "SpatRaster")) {
      x <- terra::classify(x, cbind(0, Inf, 1))
      if (NA_to_0) x <- terra::classify(x, cbind(NA, 0))
      if (Zero_to_NA) x <- terra::classify(x, cbind(0, NA))
    } else {
      stop("Input map should be either PackedSpatRaster, RasterLayer, or SpatRaster")
    }
  }
  x
}
