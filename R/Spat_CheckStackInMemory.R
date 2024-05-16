# |---------------------------------------------------| #
# CheckStackInMemory ------
# |---------------------------------------------------| #

#' Check if the raster stack reads from disk or memory
#'
#' Check if the raster stack reads from disk or memory
#' @author Ahmed El-Gabbas
#' @export
#' @param Stack Stack
#' @examples
#' LoadPackages(raster)
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

  if (!inherits(Stack, "RasterStack")) {
    message("The object should be a Raster stack object")
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    invisible(stop())
  }

  InMem <- sapply(raster::unstack(Stack), raster::inMemory)
  if (all(InMem)) {
    message(paste0("All stack layers reads from ", crayon::bold("disk")))
  }
  if (all(!InMem)) {
    message(paste0("All stack layers reads from ", crayon::bold("memory")))
  }

  if (sum(InMem) > 0 && sum(InMem) < raster::nlayers(Stack)) {
    paste0("Layers numbered (", paste0(which(!InMem), collapse = "-"),
           ") reads from disk")
  }
}
