## |------------------------------------------------------------------------| #
# CHELSA_Project ----
## |------------------------------------------------------------------------| #

#' Project CHELSA data to the study area
#'
#' This function processes CHELSA climate data, projects it to a specified
#' reference grid, and optionally saves the output in NetCDF or TIFF format. It
#' supports downloading data from a URL, applying a land mask, and adjusting
#' data with scale and offset values.
#' @name CHELSA_Project
#' @param Metadata Single row tibble for the metadata of the input file. This
#'   should be prepared in the [CHELSA_Prepare] function and provided via the
#'   the [CHELSA_Process] function.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param CompressLevel integer; specifies the compression level for the
#'   exported NetCDF file, ranging from 1 (least compression) to 9 (most
#'   compression). Defaults to 5.
#' @param ReturnMap logical; if `TRUE`, the processed map (as a
#'   `PackedSpatRaster` object) is returned. Defaults to `FALSE`.
#' @return Depending on the `ReturnMap` parameter, either a `PackedSpatRaster`
#'   object is returned, or nothing is returned. The function always writes
#'   NetCDF and TIFF files to disk.
#' @author Ahmed El-Gabbas
#' @export

CHELSA_Project <- function(
    Metadata = NULL, EnvFile = ".env", FromHPC = TRUE,
    CompressLevel = 5, ReturnMap = FALSE) {


  # # ..................................................................... ###

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "EnvFile", Type = "character")

  LogicArgs <- c("FromHPC", "ReturnMap")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "CompressLevel", Type = "numeric")

  rm(AllArgs, LogicArgs)

  if (!inherits(Metadata, "tbl_df")) {
    stop("Input metadata has to be a tibble")
  }

  if (nrow(Metadata) != 1) {
    stop("Input metadata has to be a single-row tibble")
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- NULL

  # # ..................................................................... ###

  if (is.null(Metadata)) {
    stop("Input metadata can not be `NULL`", call. = FALSE)
  }

  if (!file.exists(Metadata$Path_Down)) {
    stop(
      paste0("Input file does not exist: ", Metadata$Path_Down),
      call. = FALSE)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use official
  # parameters (overriding the ones from GeoTIFF keys)
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # Environment variables -----

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # Loading reference grid -----

  GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop(
      paste0("Path for the Europe boundaries does not exist: ", GridR),
      call. = FALSE)
  }
  GridR <- terra::unwrap(IASDT.R::LoadAs(GridR))

  # # ..................................................................... ###

  # Land mask ------

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

  # # ..................................................................... ###

  # Read tif file as terra SpatRaster object ----

  # terra package by default considers the scale and offset information stored
  # in the tiff files. Here I disable this to read the raw values as it is and
  # later consider the scale and offset information manually. This is more safe
  # as I found that some of the future projections do not include such
  # information in the tiff files.
  Rstr <- terra::rast(Metadata$Path_Down, raw = TRUE) %>%
    stats::setNames(basename(Metadata$Path_Down)) %>%
    # crop to European boundaries although it is not necessary to crop the input
    # maps into the European boundaries, we will crop the data prior to
    # projection. Cropping will make the values of the raster read from memory
    # not from the file. This is a workaround to avoid wrong extreme values in
    # the output file because of a bug in terra package (see this issue:
    # https://github.com/rspatial/terra/issues/1356) [18.02.2023]
    terra::crop(CropExtent) %>%
    # mask by land mask
    terra::mask(LandMaskL) %>%
    # `gsp` maps contains extremely high values instead of NA; the following
    # replace extreme values with NA
    terra::classify(cbind(420000000, Inf, NA))

  rm(LandMaskL, CropExtent)

  # # ..................................................................... ###

  # Manually considering offset and scale -----

  # For `npp` layers, all tiff maps except for current climate does have a
  # scaling factor all scale and offset information were set manually
  if (Metadata$scale != 1) Rstr <- Rstr * Metadata$scale
  if (Metadata$offset != 0) Rstr <- Rstr + Metadata$offset

  # # ..................................................................... ###

  # Projecting to reference grid EPSG 3035 ------

  Rstr <- Rstr %>%
    # project to reference grid
    terra::project(GridR, method = "average", threads = TRUE) %>%
    # mask to the reference grid
    terra::mask(GridR) %>%
    # Ensure that values are read from memory
    IASDT.R::setRastVals() %>%
    IASDT.R::setRastCRS()

  # # ..................................................................... ###

  # Write file to disk --- tiff -----

  terra::writeRaster(
    x = Rstr, filename = Metadata$Path_Out_tif, overwrite = TRUE)

  # # ..................................................................... ###

  # Write file to disk --- nc -----

  # Variable name of the output *.nc file
  VarName4NC <- c(
    Metadata$TimePeriod, Metadata$ClimateModel, Metadata$ClimScenario) %>%
    unique() %>%
    paste0(collapse = "__") %>%
    paste0(Metadata$Variable, "__", .)

  # global attributes to be added to the *.nc file
  Attrs <- c(
    paste0("URL=", Metadata$URL),
    paste0("OriginalFile=", Metadata$Path_Down),
    paste0("Variable=", Metadata$Variable),
    paste0("TimePeriod=", Metadata$TimePeriod),
    paste0("ClimModel=", Metadata$ClimateModel),
    paste0("ClimScenario=", Metadata$ClimScenario),
    paste0("Long_name=", Metadata$Long_name),
    paste0("unit=", Metadata$unit),
    paste0("explanation=", Metadata$explanation))

  # save as *.nc file
  terra::writeCDF(
    Rstr, filename = Metadata$Path_Out_NC, varname = VarName4NC,
    unit = Metadata$unit, zname = Metadata$TimePeriod,
    atts = Attrs, overwrite = TRUE, compression = CompressLevel)

  # # ..................................................................... ###

  # Return map? -----

  if (ReturnMap) {
    return(terra::wrap(Rstr))
  } else {
    return(invisible(NULL))
  }
}
