
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# CheckTiff ----
# |---------------------------------------------------| #

#' Check if tiff file corrupted
#'
#' Check if tiff file corrupted
#' @param x tiff file path
#' @name CheckTiff
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

CheckTiff <- function(x) {
  x %>%
    terra::describe() %>%
    stringr::str_detect("Driver") %>%
    any() %>%
    magrittr::not()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# SplitRaster ----
# |---------------------------------------------------| #

#' Split a raster object into a list of smaller rasters
#'
#' Split a raster object into a list of smaller rasters
#' @param raster raster object to split
#' @param Ncol number of columns
#' @param NRow number of rows
#' @param save save output to disk?
#' @param SplitPath file path
#' @param plot plot the results
#' @param Extent return only the extent of split raster
#' @name SplitRaster
#' @author Ahmed El-Gabbas
#' @return NULL
#' @details
#' #' References:
#' https://stackoverflow.com/questions/29784829/
#' https://stackoverflow.com/questions/22109774/r-raster-mosaic-from-list-of-rasters
#'
#' @export

SplitRaster <- function(
    raster, Ncol = 4, NRow = 4, save = FALSE,
    SplitPath = "", plot = FALSE, Extent = FALSE) {

  h <- ceiling(ncol(raster) / Ncol)
  v <- ceiling(nrow(raster) / NRow)
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
    for (i in seq_len(r_list)) {
      raster::writeRaster(
        x = r_list[[i]], filename = paste(SplitPath, "/SplitRas", i, sep = ""),
        format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
    }
  }

  if (plot == TRUE) {
    graphics::par(mfrow = c(NRow, Ncol))
    for (i in seq_len(r_list)) {
      raster::plot(r_list[[i]], axes = FALSE, legend = FALSE, bty = "n", box = FALSE)
    }
  }
  return(r_list)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# DownBoundary ----
# |---------------------------------------------------| #

#' Determine the boundaries of the requested GBIF data
#'
#' Determine the boundaries of the requested GBIF data
#' @param L left boundary
#' @param R right boundary
#' @param B bottom boundary
#' @param T top boundary
#' @name DownBoundary
#' @author Ahmed El-Gabbas
#' @return NULL
#' @description `rgbif::pred_within()` function used to download GBIF data only accepts a WKT string. This function takes the values of the boundary and converts it to a WKT string. Default values are determined by the variables: Bound_L, R = Bound_R, Bound_B, Bound_T...
#' @export

DownBoundary <- function(L = Bound_L, R = Bound_R, B = Bound_B, T = Bound_T) {
  "POLYGON(({L} {B},{R} {B},{R} {T},{L} {T},{L} {B}))" %>%
    stringr::str_glue() %>%
    as.character()
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# rename_geometry ----
# |---------------------------------------------------| #

#' Rename active geometry column of an sf object
#'
#' Rename active geometry column of an sf object
#' @param g simple feature object
#' @param name the new name of geometry
#' @name rename_geometry
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

rename_geometry <- function(g, name) {
  # Source: https://gis.stackexchange.com/a/386589/30390
  # https://gis.stackexchange.com/questions/386584/sf-geometry-column-naming-differences-r
  current <- attr(g, "sf_column")
  names(g)[names(g) == current] <- name
  sf::st_geometry(g) <- name
  return(g)
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Polygon_Centroid ----
# |---------------------------------------------------| #

#' Replace the geometry of a polygon with its centroid point
#'
#' Replace the geometry of a polygon with its centroid point
#' @param x simple feature object
#' @param Rename should the geometry field renamed
#' @param NewName the new name of geometry
#' @name Polygon_Centroid
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Polygon_Centroid <- function(x, Rename = FALSE, NewName = "") {
  # https://github.com/r-spatial/sf/issues/480
  suppressWarnings(sf::st_geometry(x) <- sf::st_geometry(sf::st_centroid(x)))
  if (Rename) {
    x %>%
      rename_geometry(name = NewName) %>%
      return()
  } else {
    return(x)
  }
}



# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Set_geometry ----
# |---------------------------------------------------| #

#' Set geometry of an sf object in the pipe line
#'
#' Set geometry of an sf object in the pipe line
#' @param x simple feature data frame
#' @param Name the name of the geometry column to be used
#' @name Set_geometry
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Set_geometry <- function(x, Name) {
  sf::st_geometry(x) <- Name
  return(x)
}
