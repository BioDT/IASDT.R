## |------------------------------------------------------------------------| #
# set_raster_values ----
## |------------------------------------------------------------------------| #

#' Set values for `SpatRaster` Objects
#'
#' This function processes a `SpatRaster` object or a file path to it, ensuring
#' the raster is loaded and unpacked if necessary. It converts raster objects
#' from the `raster` package to `SpatRaster` objects. The function ensures that
#' the `SpatRaster` object is read from memory, not from file.
#' @param raster A `SpatRaster` object, a file path to a raster file, or an
#'   object from the `raster` package (e.g., `RasterLayer`, `RasterStack`,
#'   `RasterBrick`).
#' @return A `SpatRaster` object with values loaded into memory.
#' @name set_raster_values
#' @author Ahmed El-Gabbas
#' @details The function handles various types of input:
#'   - If a file path is provided, it attempts to load the raster using
#'   [terra::rast()].
#'   - If the input is a packed `SpatRaster`, it unpacks the raster using
#'   [terra::unwrap()].
#'   - If the input is a raster object from the `raster` package, it is
#'   converted to `SpatRaster`.
#' @export

set_raster_values <- function(raster) {

  # If raster is character object, try to read it as SpatRast object
  if (inherits(raster, "character")) {
    if (!file.exists(raster)) {
      stop("Input file path does not exist", call. = FALSE)
    }
    raster <- tryCatch(
      terra::rast(raster),
      error = function(e) {
        stop(
          "Failed to load raster from the provided file path: ", e$message,
          call. = FALSE)
      }
    )
  }

  # Convert raster package objects to SpatRaster
  if (!inherits(raster, "SpatRaster")) {
    if (!inherits(raster, c("RasterLayer", "RasterStack", "RasterBrick"))) {
      stop(
        "Input object must be a `SpatRaster`, `RasterLayer`, ",
        "`RasterStack`, or `RasterBrick` object", call. = FALSE)
    }
    raster <- terra::rast(raster)
  }


  # Unwrap when necessary
  if (inherits(raster, "PackedSpatRaster")) {
    raster <- tryCatch(
      terra::unwrap(raster),
      error = function(e) {
        stop("Failed to unwrap PackedSpatRaster: ", e$message, call. = FALSE)
      })
  }

  # Set values
  if (!all(terra::inMemory(raster))) {
    terra::values(raster) <- terra::values(raster)
  }

  return(raster)
}
