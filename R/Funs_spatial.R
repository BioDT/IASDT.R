
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
#' @return logical: `TRUE` when the tiff file is okay; `FALSE` when corrupted
#' @export
#' @examples
#' (f <- system.file("ex/elev.tif", package="terra"))
#'
#' CheckTiff(f)

CheckTiff <- function(x) {
  x %>%
    terra::describe() %>%
    stringr::str_detect("Driver") %>%
    any()
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
#' @return WKT string
#' @export
#' @description `rgbif::pred_within()` function used to download GBIF data only accepts a WKT string. This function takes the values of the boundary and converts it to a WKT string. Default values are determined by the variables: Bound_L, R = Bound_R, Bound_B, Bound_T...
#' @examples
#' IASDT.R::DownBoundary(20, 30, 40, 50)

DownBoundary <- function(L, R, B, T) {
  "POLYGON(({L} {B},{R} {B},{R} {T},{L} {T},{L} {B}))" %>%
    stringr::str_glue() %>%
    as.character()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Rename_geometry ----
# |---------------------------------------------------| #

#' Rename active geometry column of an sf object
#'
#' Rename active geometry column of an sf object
#' @param g simple feature object
#' @param name the new name of geometry
#' @name Rename_geometry
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Rename_geometry <- function(g, name) {
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
      Rename_geometry(name = NewName) %>%
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


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Text2Coords ----
# |---------------------------------------------------| #

#' Extract longitude / latitude from text
#'
#' Extract longitude / latitude from text
#' @param String string of coordinate
#' @param Long_Name Column name for the longitude information; default: "Longitude"
#' @param Lat_Name Column name for the latitude information; default: "Latitude"
#' @name Text2Coords
#' @author Ahmed El-Gabbas
#' @return two column tibble for Longitude & Latitude
#' @export
#' @examples
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)", "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords)
#'
#' c("POINT (11.761 46.286)", "POINT (14.8336 42.0422)", "POINT (16.179999 38.427214)") %>%
#'  lapply(Text2Coords, Long_Name = "Long", Lat_Name = "Lat")

Text2Coords <- function(String, Long_Name = "Longitude", Lat_Name = "Latitude") {
  String %>%
    # convert string to 2-columns data frame
    stringr::str_remove_all("POINT \\(|\\)") %>%
    stringr::str_split(" ", simplify = TRUE) %>%
    as.numeric() %>%
    matrix(nrow = 1, ncol = 2) %>%
    as.data.frame() %>%
    stats::setNames(c(Long_Name, Lat_Name)) %>%
    tibble::tibble()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# ClipRasterByPolygon ------
# |---------------------------------------------------| #

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


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# sf_add_coords ------
# |---------------------------------------------------| #

#' Add coordinates of the current sf object as columns
#'
#' Add coordinates of the current sf object as columns
#' @name sf_add_coords
#' @param Sf_Obj input sf object
#' @param NameX Name of the longitude column to be added
#' @param NameY Name of the latitude column to be added
#' @param Overwrite Should columns `NameX` or `NameY` be overwritten if these columns exist in the original data
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' pt1 = sf::st_point(c(0,1))
#' pt2 = sf::st_point(c(1,1))
#' d = data.frame(a = 1:2)
#' d$geom = sf::st_sfc(pt1, pt2)
#' df = sf::st_as_sf(d)
#' df
#' (df <- sf_add_coords(df))
#'
#' (sf_add_coords(df))
#'
#' (sf_add_coords(df, Overwrite = TRUE))

sf_add_coords <- function(Sf_Obj, NameX = "Long", NameY = "Lat", Overwrite = FALSE) {
  ColNames <- names(Sf_Obj)
  Coords <- Sf_Obj %>%
    sf::st_coordinates() %>%
    tibble::as_tibble() %>%
    stats::setNames(c(NameX, NameY))

  if (any(NameX %in% ColNames || NameY %in% ColNames)) {
    if (Overwrite) {
      warning("Provided column names for longitude and Latitude already exist in the data; these columns were overwritten")
      Sf_Obj <- dplyr::select(Sf_Obj, -dplyr::all_of(c(NameX, NameY)))
    } else {
      warning('Provided column names for longitude and Latitude already exist in the data; "_NEW" is used as suffix')
      Coords <- Coords %>%
        stats::setNames(c(paste0(NameX, "_NEW"), paste0(NameY, "_NEW")))
    }
  } else {
    Coords <- stats::setNames(Coords, c(NameX, NameY))
  }
  return(dplyr::bind_cols(Sf_Obj, Coords))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# GridDiagOff ------
# |---------------------------------------------------| #

#' Create a `multilinestring` (diagonal and off-diagonal lines) sf object from each grid cell
#'
#' Create a `multilinestring` (diagonal and off-diagonal lines) sf object from each grid cell
#'
#' @name GridDiagOff
#' @param DT input sf tibble
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(dplyr, sf, raster, ggplot2)
#'
#' Grid <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10,
#'                        ymn = 0, ymx = 10, crs = 4326) %>%
#'   setNames("Grid") %>%
#'   raster::setValues(1) %>%
#'   raster::rasterToPolygons() %>%
#'   sf::st_as_sf()
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()
#'
#' Grid_X <- GridDiagOff(Grid)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
#'                    linewidth = 0.5, inherit.aes = TRUE) +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()

GridDiagOff <- function(DT) {

  InputCRS <- sf::st_crs(DT)

  DT %>%
    dplyr::pull("geometry") %>%
    purrr::map(
      .f = ~{
        SS2 <- .x %>%
          sf::st_coordinates() %>%
          tibble::as_tibble() %>%
          dplyr::select(X, Y)

        OffDiag <- dplyr::bind_rows(
          dplyr::filter(SS2, X == min(X), Y == max(Y)),
          dplyr::filter(SS2, X == max(X), Y == min(Y))) %>%
          as.matrix() %>%
          sf::st_linestring()

        Diag <- dplyr::bind_rows(
          dplyr::filter(SS2, X == min(X), Y == min(Y)),
          dplyr::filter(SS2, X == max(X), Y == max(Y))) %>%
          as.matrix() %>%
          sf::st_linestring()

        list(OffDiag, Diag) %>%
          sf::st_multilinestring() %>%
          sf::st_geometry()  %>%
          sf::st_set_crs(InputCRS) %>%
          sf::st_sfc() %>%
          return()

      }) %>%
    tibble::tibble(geometry = .) %>%
    tidyr::unnest("geometry") %>%
    sf::st_as_sf()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# GridCross ------
# |---------------------------------------------------| #

#' Create a `multilinestring` (cross in the middle of the grid) sf object from each grid cell
#'
#' Create a `multilinestring` (cross in the middle of the grid) sf object from each grid cell
#'
#' @name GridCross
#' @param DT input sf tibble
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(dplyr, sf, raster, ggplot2)
#'
#' Grid <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10, crs = 4326) %>%
#'   setNames("Grid") %>%
#'   raster::setValues(1) %>%
#'   raster::rasterToPolygons() %>%
#'   sf::st_as_sf()
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()
#'
#' Grid_X <- GridCross(Grid)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
#'                    linewidth = 0.5, inherit.aes = TRUE) +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()

GridCross <- function(DT) {

  InputCRS <- sf::st_crs(DT)

  DT %>%
    dplyr::pull("geometry") %>%
    purrr::map(
      .f = ~{
        SS2 <- .x %>%
          sf::st_coordinates() %>%
          tibble::as_tibble() %>%
          dplyr::select(X, Y)

        Horiz <- SS2 %>%
          dplyr::reframe(
            X = X,
            Y = (min(Y) + (max(Y) - min(Y)) / 2)) %>%
          dplyr::distinct() %>%
          as.matrix() %>%
          sf::st_linestring()
        Vert <- SS2 %>%
          dplyr::reframe(
            X = (min(X) + (max(X) - min(X)) / 2),
            Y = Y) %>%
          dplyr::distinct() %>%
          as.matrix() %>%
          sf::st_linestring()

        list(Horiz, Vert) %>%
          sf::st_multilinestring() %>%
          sf::st_geometry()  %>%
          sf::st_set_crs(InputCRS) %>%
          sf::st_sfc() %>%
          return()
      }) %>%
    tibble::tibble(geometry = .) %>%
    tidyr::unnest("geometry") %>%
    sf::st_as_sf()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# RastPA ------
# |---------------------------------------------------| #

#' convert raster map into binary (1/0)
#'
#' convert raster map into binary (1/0)
#'
#' @name RastPA
#' @param x input raster map
#' @param NA_to_0 logical (default: `TRUE`); should NA be replaced with 0
#' @param Zero_to_NA logical (default: `FALSE`); should 0 be replaced with NA
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(dplyr, raster, ggplot2, tidyterra)
#'
#' r <- raster::raster(system.file("external/test.grd", package = "raster"))
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::rast(r), maxcell = Inf) +
#'   ggplot2::theme_minimal()
#'
#' R2 <- raster::stack(
#'   RastPA(r),                            # NA replaced with 0
#'   RastPA(r, NA_to_0 = FALSE),           # NA is kept as NA
#'   RastPA(RastPA(r), Zero_to_NA = TRUE)) # 0 replaced with NA in the second map
#'
#' ggplot2::ggplot() +
#'   tidyterra::geom_spatraster(data = terra::as.factor(terra::rast(R2)), maxcell = Inf) +
#'   ggplot2::facet_wrap(~lyr) +
#'   ggplot2::scale_fill_manual(values = c("grey30", "red"), na.value = "transparent") +
#'   ggplot2::theme_minimal()

RastPA <- function(x, NA_to_0 = TRUE, Zero_to_NA = FALSE) {
  MaxVal <- raster::cellStats(x, max)
  if (MaxVal > 0) x[x > 0] <- 1
  if (NA_to_0) x <- raster::reclassify(x, cbind(NA, 0))
  if (Zero_to_NA) x <- raster::reclassify(x, cbind(0, NA))
  x
}
