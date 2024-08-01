## |------------------------------------------------------------------------| #
# CheckStackInMemory ------
## |------------------------------------------------------------------------| #

#' Check if a raster stack reads from disk or memory
#'
#' This function checks whether the layers of a RasterStack object are stored in
#' memory or read from disk.  It prints messages indicating whether all layers
#' are in memory, all layers are on disk, or a mix of both. If there's a mix, it
#' specifies which layers are on disk.
#' @author Ahmed El-Gabbas
#' @export
#' @name CheckStackInMemory
#' @param Stack A RasterStack object. If `NULL` or not a RasterStack, the
#'   function will stop with an error.
#' @return No return value, but prints messages to the console.
#' @examples
#' library(raster)
#' logo <- raster(system.file("external/rlogo.grd", package = "raster"))
#' logo@data@inmemory
#' logo@data@fromdisk
#' logo@file@name
#'
#' # -------------------------------------------
#'
#' # A raster stack reading from files
#' ST2 <- raster::stack(logo, logo)
#' CheckStackInMemory(ST2)
#' c(ST2[[1]]@data@inmemory, ST2[[2]]@data@inmemory)
#' c(ST2[[1]]@data@fromdisk, ST2[[2]]@data@fromdisk)
#' c(ST2[[1]]@file@name, ST2[[2]]@file@name)
#'
#' # -------------------------------------------
#'
#' logo2 <- raster::readAll(logo)
#' ST3 <- raster::stack(logo, logo2)
#' CheckStackInMemory(ST3)
#' c(ST3[[1]]@data@inmemory, ST3[[2]]@data@inmemory)
#' c(ST3[[1]]@data@fromdisk, ST3[[2]]@data@fromdisk)
#' c(ST3[[1]]@file@name, ST3[[2]]@file@name)

CheckStackInMemory <- function(Stack = NULL) {

  # Check input argument
  if (is.null(Stack)) {
    stop("Input Stack cannot be NULL")
  }

  if (magrittr::not(inherits(Stack, "RasterStack"))) {
    stop("The object should be a RasterStack object")
  }

  InMem <- sapply(raster::unstack(Stack), raster::inMemory)

  if (all(InMem)) {
    message(paste0("All stack layers reads from ", crayon::bold("disk")))
  }
  if (all(magrittr::not(InMem))) {
    message(paste0("All stack layers reads from ", crayon::bold("memory")))
  }

  if (sum(InMem) > 0 && (sum(InMem) < raster::nlayers(Stack))) {
    paste0("Layers numbered (",
      paste0(which(!InMem), collapse = "-"), ") reads from disk")
  }
}
