## |------------------------------------------------------------------------| #
# Chelsa_Project ----
## |------------------------------------------------------------------------| #

#' Project Chelsa data to the study area
#'
#' Project Chelsa data to the study area
#'
#' @name Chelsa_Project
#' @param InputFile character; file path or URL for input tif file
#' @param OutFile character; file path for output nc/tif file. File extension will be replaced when needed.
#' @param GridFile reference grid; either as `RasterLayer`, `SpatRaster`, or `PackedSpatRaster`. The reference grid will be used for projection and the resulted maps will be masked to it
#' @param ReturnMap logical; should the processed map be returned by the end of the function? Default: `FALSE`
#' @param DownPath character; Where to save downloaded files (if the `InputFile` is a URL)
#' @param KeepDownloaded logical; if URL is provided as input file, the file will be downloaded to disk first before processing. Should the downloaded file be kept in disk. Default: `FALSE`
#' @param SaveTiff logical; also save output map as *.tif file. Default: `FALSE`
#' @param CompressLevel integer; compression level of the exported NetCDF file (see: `terra::writeCDF` for more information). Can be set to an integer between 1 (least compression) and 9 (most compression). Default: `5`.
#' @returns If `ReturnMap = TRUE`, a wrapped SpatRaster object (`PackedSpatRaster`) is returned; otherwise nothing is returned. By default the function exports a NetCDF file. If `SaveTiff` is set as `TRUE`, additional tiff file will be saved to disk.
#' @author Ahmed El-Gabbas
#' @export

Chelsa_Project <- function(
    InputFile = NULL, OutFile = NULL, GridFile = NULL, ReturnMap = FALSE,
    DownPath = NULL, KeepDownloaded = TRUE, SaveTiff = FALSE,
    CompressLevel = 5) {

  # Ensure that the reference grid is not null
  if (is.null(GridFile)) stop("GridFile can not be empty")

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use official parameters (overriding the ones from GeoTIFF keys)
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # Input file name
  InputName <- basename(InputFile)

  # Variable description
  VarDesc <- IASDT.R::Chelsa_Info(InputFile)
  # Scale and offset information
  VarScale <- VarDesc$scale
  VarOffset <- VarDesc$offset

  # ||||||||||||||||||||||||||||||||||||||||
  # Ensure that some packages are loaded
  # ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::LoadPackages("dplyr", "raster", "terra") %>%
    suppressMessages() %>%
    suppressWarnings()

  # ||||||||||||||||||||||||||||||||||||||||
  # Loading reference grid
  # ||||||||||||||||||||||||||||||||||||||||

  if (inherits(GridFile, "RasterLayer")) {
    GridR <- terra::rast(GridFile)
  } else {
    if (inherits(GridFile, "PackedSpatRaster")) {
      GridR <- terra::unwrap(GridFile)
    } else {
      GridR <- GridFile
    }
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Remote or local
  # ||||||||||||||||||||||||||||||||||||||||

  Remote <- dplyr::if_else(stringr::str_detect(InputFile, "^http"), TRUE, FALSE)

  if (Remote) {

    if (is.null(OutFile)) stop("OutFile can not be empty if the input file is an URL")

    # Folder to download the file
    if (is.null(DownPath)) {
      # download as temporary file
      DownPath <- tempfile(fileext = ".tif")
      # delete the temporary file after processing
      KeepDownloaded <- FALSE
    } else {
      DownPath <- file.path(DownPath, InputName)
    }

    # Ensure that output dir exists
    fs::dir_create(dirname(DownPath))

    # increase time out to allow more time to download input file
    # This may not be necessary for LUMI HPC
    # EVE has a download limit. This allows more time for download
    options(timeout = max(300, getOption("timeout")))

    # Download file to disk
    utils::download.file(
      url = InputFile, destfile = DownPath, quiet = TRUE, mode = "wb")

  } else {
    # if input file is located locally,
    DownPath <- InputFile
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Land mask
  # ||||||||||||||||||||||||||||||||||||||||

  # Extent to crop the maps prior to processing.
  # This ensures that the object reads from the memory. See below
  CropExtent <- terra::ext(-26, 37.5, 34, 72)

  # source: CHELSA-W5E5 v1.0: W5E5 v1.0 downscaled with CHELSA v2.0
  # https://data.isimip.org/10.48364/ISIMIP.836809.3
  # Version: 1.0.3
  # This file was copied to the package data to make it easier to use it
  LandMaskL <- system.file(
    "extdata", "LandMask.nc", package = "IASDT.R", mustWork = TRUE) %>%
    terra::rast() %>%
    terra::crop(CropExtent) %>%
    suppressWarnings() %>%  # suppress warning on LUMI while cropping
    terra::classify(cbind(0, NA))

  # ||||||||||||||||||||||||||||||||||||||||
  # read tif file as terra SpatRaster object
  # ||||||||||||||||||||||||||||||||||||||||

  # terra package by default considers the scale and offset information stored in the tiff files. Here I disable this to read the raw values as it is and later consider the scale and offset information manually. This is more safe as I found that some of the future projections do not include such information in the tiff files.
  Rstr <- terra::rast(DownPath, raw = TRUE) %>%
    stats::setNames(paste0(tools::file_path_sans_ext(InputName), ".tif")) %>%
    # crop to European boundaries
    # although it is not necessary to crop the input maps into the European boundaries, we will crop the data prior to projection. Cropping will make the values of the raster read from memory not from the file. This is a workaround to avoid wrong extreme values in the output file because of a bug in terra package (see this issue: https://github.com/rspatial/terra/issues/1356) [18.02.2023]
    terra::crop(CropExtent) %>%
    # mask by land mask
    terra::mask(LandMaskL) %>%
    # `gsp` maps contains extremely high values instead of NA; the following replace extreme values with NA
    terra::classify(cbind(420000000, Inf, NA))

  rm(LandMaskL, CropExtent)

  # ||||||||||||||||||||||||||||||||||||||||
  # Manually considering offset and scale
  # ||||||||||||||||||||||||||||||||||||||||

  # For `npp` layers, all tiff maps except for current climate does have a scaling factor
  # all scale and offset information were set manually
  if (VarScale != 1) Rstr <- Rstr * VarScale
  if (VarOffset != 0) Rstr <- Rstr + VarOffset

  # ||||||||||||||||||||||||||||||||||||||||
  # Projecting to reference grid EPSG 3035
  # ||||||||||||||||||||||||||||||||||||||||

  # projection method
  # Use `mode` for `kg` variables (categorical variables)
  Method <- dplyr::if_else(
    stringr::str_detect(names(Rstr), "kg[0-5]"), "mode", "average")

  Rstr <- Rstr %>%
    # project to reference grid
    terra::project(GridR, method = Method, threads = TRUE) %>%
    # mask to the reference grid
    terra::mask(GridR)

  # Ensure that the object is located in memory, not reading from temporary file
  # This may not be necessary as we save the file as .tif file not .RData
  if (magrittr::not(terra::inMemory(Rstr))) {
    terra::values(Rstr) <- terra::values(Rstr)
  }

  terra::crs(Rstr) <- "epsg:3035"

  # ||||||||||||||||||||||||||||||||||||||||
  # Write file to disk --- tiff
  # ||||||||||||||||||||||||||||||||||||||||

  if (SaveTiff) {

    OutFileTif <- dplyr::if_else(
      stringr::str_detect(OutFile, ".tif$"),
      OutFile,
      file.path(dirname(OutFile), paste0(
        tools::file_path_sans_ext(InputName), ".tif")))

    terra::writeRaster(x = Rstr, filename = OutFileTif, overwrite = TRUE)

    if (Remote && magrittr::not(KeepDownloaded)) file.remove(DownPath)
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Write file to disk --- nc
  # ||||||||||||||||||||||||||||||||||||||||

  OutFileNC <- dplyr::if_else(
    stringr::str_detect(OutFile, ".nc$"),
    OutFile,
    file.path(dirname(OutFile),
              paste0(tools::file_path_sans_ext(InputName), ".nc")))

  # Variable name of the output *.nc file
  VarName4NC <- c(
    VarDesc$TimePeriod, VarDesc$ClimModel, VarDesc$ClimScenario) %>%
    unique() %>%
    paste0(collapse = "__") %>%
    paste0(VarDesc$Variable, "__", .)

  # global attributes to be added to the *.nc file
  Attrs <- c(
    paste0("Var=", VarDesc$Variable),
    paste0("TimePeriod=", VarDesc$TimePeriod),
    paste0("ClimModel=", VarDesc$ClimModel),
    paste0("ClimScenario=", VarDesc$ClimScenario),
    paste0("Long_name=", VarDesc$Long_name),
    paste0("unit=", VarDesc$unit),
    paste0("explanation=", VarDesc$explanation))

  # save as *.nc file
  terra::writeCDF(
    Rstr, filename = OutFileNC, varname = VarName4NC,
    unit = VarDesc$unit, zname = VarDesc$TimePeriod,
    atts = Attrs, overwrite = TRUE, compression = CompressLevel)

  # ||||||||||||||||||||||||||||||||||||||||
  # Return map?
  # ||||||||||||||||||||||||||||||||||||||||

  invisible(gc())

  if (ReturnMap) {
    return(terra::wrap(Rstr))
  } else {
    return(invisible(NULL))
  }
}
