## |------------------------------------------------------------------------| #
# setRastCRS  ------
## |------------------------------------------------------------------------| #

#' sets CRS for a SpatRaster
#'
#' This function sets the coordinate reference system (CRS) for a SpatRaster
#' object using the specified EPSG code. This is a wrapper function for
#' `terra::crs(R) <- CRS` but allowing to set the CRS in the pipe.
#' @name setRastCRS
#' @param R A SpatRaster object whose CRS needs to be set.
#' @param CRS A string specifying the CRS to be set, default is "epsg:3035".
#' @return The SpatRaster object with the updated CRS.
#' @author Ahmed El-Gabbas
#' @export

setRastCRS <- function(R, CRS = "epsg:3035") {
  terra::crs(R) <- CRS
  return(R)
}
