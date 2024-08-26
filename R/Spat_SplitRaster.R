## |------------------------------------------------------------------------| #
# SplitRaster ----
## |------------------------------------------------------------------------| #

#' Split a raster object into a list of smaller rasters
#'
#' Split a raster object into a list of smaller rasters based on specified
#' numbers of rows and columns. It can optionally save the resulting rasters to
#' disk, plot them, or return just their extents.
#' @param raster A raster object to be split. If NULL (the default), the
#'   function will not execute.
#' @param Ncol,Nrow Integer. The desired number of columns and rows to split the
#'   raster into. Default is 4 columns and 4 rows.
#' @param save Logical. Indicates whether to save the split rasters to disk.
#'   Default is `FALSE`.
#' @param SplitPath Character string. Specifies the directory path where the
#'   split rasters should be saved if `save` is `TRUE`. If the directory does
#'   not exist, it will be created.
#' @param plot Logical. Indicates whether to plot the split rasters. Default is
#'   `FALSE`.
#' @param Extent Logical. If `TRUE`, the function returns only the extents of
#'   the split rasters instead of the raster data. Default is `FALSE`.
#' @name SplitRaster
#' @author Ahmed El-Gabbas
#' @return A list of raster objects or extents of the split rasters, depending
#'   on the value of the `Extent` parameter.
#' @export
#' @references Click [here](https://stackoverflow.com/questions/29784829/) and
#'   [here](https://stackoverflow.com/q/22109774)
#' @examples
#' library(raster)
#' logo <- raster(system.file("external/rlogo.grd", package = "raster"))
#' plot(logo, axes = FALSE, legend = FALSE, bty = "n",
#'      box = FALSE, main = "Original raster layer")
#' # --------------------------------------------------
#'
#' # Split into 3 rows and 3 columns
#' logoSplit <- SplitRaster(raster = logo, Ncol = 3, Nrow = 3, plot = TRUE)
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
    raster = NULL, Ncol = 4, Nrow = 4, save = FALSE,
    SplitPath = "", plot = FALSE, Extent = FALSE) {

  # Directory Check
  if (save && !dir.exists(SplitPath)) {
    fs::dir_create(SplitPath)
  }

  # Check input raster
  if (is.null(raster)) {
    stop("raster cannot be NULL.", call. = FALSE)
  }

  h <- ceiling((ncol(raster) / Ncol))
  v <- ceiling((nrow(raster) / Nrow))
  agg <- raster::aggregate(raster, fact = c(h, v))
  agg[] <- seq_len(raster::ncell(agg))
  agg_poly <- raster::rasterToPolygons(agg)
  names(agg_poly@data) <- "polis"

  r_list <- list()

   for (i in seq_len(raster::ncell(agg))) {
    e1 <- raster::extent(agg_poly[agg_poly@data$polis == i, ])

    if (Extent) {
      r_list[[i]] <- e1
    } else {
      r_list[[i]] <- raster::crop(raster, e1)
    }
  }

  if (save) {
    for (i in seq_along(r_list)) {
      raster::writeRaster(
        x = r_list[[i]],
        filename = paste0(SplitPath, "/SplitRas", i),
        format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
    }
  }

  if (plot) {
    graphics::par(mfrow = c(Nrow, Ncol))
    for (i in seq_along(r_list)) {
      raster::plot(
        r_list[[i]], axes = FALSE, legend = FALSE, bty = "n",  box = FALSE)
    }
  }
  return(r_list)
}
