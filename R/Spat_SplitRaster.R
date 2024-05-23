## |------------------------------------------------------------------------| #
# SplitRaster ----
## |------------------------------------------------------------------------| #

#' Split a raster object into a list of smaller rasters
#'
#' Split a raster object into a list of smaller rasters
#' @param raster raster object to split
#' @param Ncol number of columns
#' @param Nrow number of rows
#' @param save save output to disk?
#' @param SplitPath file path
#' @param plot plot the results
#' @param Extent return only the extent of split raster
#' @name SplitRaster
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @details
#' #' References:
#' https://stackoverflow.com/questions/29784829/
#' https://stackoverflow.com/questions/22109774/r-raster-mosaic-from-list-of-rasters
#' @examples
#' LoadPackages(raster)
#' logo <- raster(system.file("external/rlogo.grd", package = "raster"))
#' plot(logo, axes = FALSE, legend = FALSE, bty = "n",
#'      box = FALSE, main = "Original raster layer")
#' # --------------------------------------------------
#'
#' # Split into 3 rows and 3 columns
#' logoSplit <- SplitRaster(logo, Ncol = 3, Nrow = 3, plot = TRUE)
#'
#' print(logoSplit) # a list object of 9 items
#'
#' # --------------------------------------------------
#'
#' # Merging split maps again
#' logoSplit$fun <- mean
#' logoSplit$na.rm <- TRUE
#' logoSplit2 <- do.call(mosaic, logoSplit)
#' par(mfrow = c(1, 1))
#' plot(logoSplit2, axes = FALSE, legend = FALSE, bty = "n",
#'      box = FALSE, main = "Merged raster layers")
#'
#' print({logoSplit2 - logo}) # No value difference!
#'
#' # --------------------------------------------------
#'
#' logoSplit <- SplitRaster(logo, Ncol = 3, Nrow = 3, Extent = TRUE)
#' print(logoSplit)
#'

SplitRaster <- function(
    raster, Ncol = 4, Nrow = 4, save = FALSE,
    SplitPath = "", plot = FALSE, Extent = FALSE) {

  h <- ceiling(ncol(raster) / Ncol)
  v <- ceiling(nrow(raster) / Nrow)
  agg <- raster::aggregate(raster, fact = c(h, v))
  agg[] <- 1:raster::ncell(agg)
  agg_poly <- raster::rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for (i in 1:raster::ncell(agg)) {
    e1 <- raster::extent(agg_poly[agg_poly$polis == i, ])
    if (Extent) {
      r_list[[i]] <- e1
    } else {
      r_list[[i]] <- raster::crop(raster, e1)
    }
  }

  if (save == TRUE) {
    for (i in seq_along(r_list)) {
      raster::writeRaster(
        x = r_list[[i]], filename = paste(SplitPath, "/SplitRas", i, sep = ""),
        format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
    }
  }

  if (plot == TRUE) {
    graphics::par(mfrow = c(Nrow, Ncol))
    for (i in seq_along(r_list)) {
      raster::plot(r_list[[i]], axes = FALSE, legend = FALSE, bty = "n", box = FALSE)
    }
  }
  return(r_list)
}
