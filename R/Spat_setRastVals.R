## |------------------------------------------------------------------------| #
# setRastVals ----
## |------------------------------------------------------------------------| #

#' Set values for `SpatRaster` Objects
#'
#' This function processes a `SpatRaster` object or a file path to it, ensuring
#' the raster is loaded and unpacked if necessary. It converts raster objects
#' from the `raster` package to `SpatRaster` objects. The function ensures that
#' the `SpatRaster` object is read from memory, not from file.
#' @param R A `SpatRaster` object, a file path to a raster file, or an object
#'   from the `raster` package (e.g., `RasterLayer`, `RasterStack`,
#'   `RasterBrick`).
#' @return A `SpatRaster` object with values loaded into memory.
#' @name setRastVals
#' @author Ahmed El-Gabbas
#' @details The function handles various types of input:
#'   - If a file path is provided, it attempts to load the raster using
#'   [terra::rast()].
#'   - If the input is a packed `SpatRaster`, it unpacks the raster using
#'   [terra::unwrap()].
#'   - If the input is a raster object from the `raster` package, it is
#'   converted to `SpatRaster`.
#' @export

setRastVals <- function(R) {

  # If R is character object, try to read it as SpatRast object
  if (inherits(R, "character")) {
    if (!file.exists(R)) {
      stop("Input file path does not exist", call. = FALSE)
    }
    R <- tryCatch(
      terra::rast(R),
      error = function(e) {
        stop("Failed to load raster from the provided file path: ",
             e$message, call. = FALSE)
      }
    )
  }

  # Convert raster package objects to SpatRaster
  if (!inherits(R, "SpatRaster")) {
    if (inherits(R, "RasterLayer") || inherits(R, "RasterStack") ||
        inherits(R, "RasterBrick")) {
      R <- terra::rast(R)
    } else {
      stop(
        paste0(
          "Input object must be a `SpatRaster`, `RasterLayer`, ",
          "`RasterStack`, or `RasterBrick` object"),
        call. = FALSE)
    }
  }

  # Unwrap when necessary
  if (inherits(R, "PackedSpatRaster")) {
    R <- tryCatch(
      terra::unwrap(R),
      error = function(e) {
        stop("Failed to unwrap PackedSpatRaster: ", e$message, call. = FALSE)
      }
    )
  }

  # Set values
  if (any(!terra::inMemory(R))) {
    terra::values(R) <- terra::values(R)
  }

  return(R)
}
