## |------------------------------------------------------------------------| #
# ClipRasterByPolygon ------
## |------------------------------------------------------------------------| #

#' Clip raster by a spatial polygon
#'
#' Clip raster by a spatial polygon
#' @param raster raster layer
#' @param shape Polygon
#' @export
#' @examples
#' LoadPackages(sp)
#' LoadPackages(raster)
#' LoadPackages(rworldmap)
#'
#' # Example Polygon
#' SPDF <- getMap(resolution = "low") %>%
#'    subset(NAME == "Germany")
#'
#' # Example RasterLayer
#' r <- raster::raster(nrow = 1e3, ncol = 1e3, crs = proj4string(SPDF))
#' r[] <- 1:length(r)
#' plot(r)
#' plot(SPDF, add = TRUE)
#'
#' # ----------------------------------
#'
#' SPDF_DE <- ClipRasterByPolygon(r, SPDF)
#' plot(raster::extent(SPDF_DE), axes = FALSE, xlab = "", ylab = "")
#' plot(SPDF_DE, add = TRUE)
#' plot(SPDF, add = TRUE)

ClipRasterByPolygon <- function(raster = NULL, shape = NULL) {
  a1_crop <- raster::crop(raster, shape)
  step1 <- raster::rasterize(shape, a1_crop)
  ClippedRaster <- a1_crop * step1
  names(ClippedRaster) <- names(raster)
  return(ClippedRaster)
}
