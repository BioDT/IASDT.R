## |------------------------------------------------------------------------| #
# CHELSA_Project ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 3

CHELSA_Project <- function(
    Metadata = NULL, EnvFile = ".env", CompressLevel = 5) {

  # Checking input arguments -----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "EnvFile", Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = "CompressLevel", Type = "numeric")
  rm(AllArgs, envir = environment())

  if (!inherits(Metadata, "tbl_df")) {
    stop("Input metadata has to be a tibble", call. = FALSE)
  }

  if (nrow(Metadata) != 1) {
    stop("Input metadata has to be a single-row tibble", call. = FALSE)
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
    stop("Input file does not exist: ", Metadata$Path_Down, call. = FALSE)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # Environment variables -----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Loading reference grid -----

  GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop(
      "Path for the Europe boundaries does not exist: ", GridR, call. = FALSE)
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
    suppressWarnings() %>% # suppress warning on LUMI while cropping
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

  rm(LandMaskL, CropExtent, envir = environment())

  # # ..................................................................... ###

  # Manually considering offset and scale -----

  # For `npp` layers, all tiff maps except for current climate does have a
  # scaling factor all scale and offset information were set manually
  if (Metadata$scale != 1) {
    Rstr <- Rstr * Metadata$scale
  }
  if (Metadata$offset != 0) {
    Rstr <- Rstr + Metadata$offset
  }

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
    Metadata$TimePeriod, Metadata$ClimateModel, Metadata$ClimateScenario) %>%
    unique() %>%
    paste(collapse = "__") %>%
    paste0(Metadata$Variable, "__", .)

  # global attributes to be added to the *.nc file
  Attrs <- c(
    paste0("URL=", Metadata$URL),
    paste0("OriginalFile=", Metadata$Path_Down),
    paste0("Variable=", Metadata$Variable),
    paste0("TimePeriod=", Metadata$TimePeriod),
    paste0("ClimateModel=", Metadata$ClimateModel),
    paste0("ClimateScenario=", Metadata$ClimateScenario),
    paste0("Long_name=", Metadata$Long_name),
    paste0("unit=", Metadata$unit),
    paste0("explanation=", Metadata$explanation))

  # save as *.nc file
  terra::writeCDF(
    Rstr,
    filename = Metadata$Path_Out_NC, varname = VarName4NC,
    unit = Metadata$unit, zname = Metadata$TimePeriod,
    atts = Attrs, overwrite = TRUE, compression = CompressLevel)

  # # ..................................................................... ###

  return(invisible(NULL))
}
