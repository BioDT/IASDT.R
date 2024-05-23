## |------------------------------------------------------------------------| #
# Scale_0_1 ----
## |------------------------------------------------------------------------| #

#' Scale a spatRaster object between 0 and 1
#'
#' Scale a spatRaster object between 0 and 1
#' 
#' @param x x
#' @name Scale_0_1
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Scale_0_1 <- function(x) {
  stopifnot("x should be a SpatRaster object" = inherits(x, "SpatRaster"))
  Num <- x - unlist(terra::global(x, "min", na.rm = TRUE))
  Denom <- unlist(terra::global(x, "max", na.rm = TRUE)) - unlist(terra::global(x, "min", na.rm = TRUE))
  (Num / Denom)
}
