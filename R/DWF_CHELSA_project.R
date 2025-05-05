## |------------------------------------------------------------------------| #
# CHELSA_project ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 3

CHELSA_project <- function(
    metadata = NULL, env_file = ".env", compression_level = 5L) {

  # Checking input arguments -----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_to_check = "env_file", args_type = "character")
  ecokit::check_args(
    args_all = AllArgs, args_to_check = "compression_level",
    args_type = "numeric")
  rm(AllArgs, envir = environment())

  if (!inherits(metadata, "tbl_df")) {
    ecokit::stop_ctx(
      "Input metadata has to be a tibble", class_metadata = class(metadata))
  }

  if (nrow(metadata) != 1) {
    ecokit::stop_ctx(
      "Input metadata has to be a single-row tibble",
      metadata = metadata, class_metadata = class(metadata),
      nrow_metadata = nrow(metadata))
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- NULL

  # # ..................................................................... ###

  if (is.null(metadata)) {
    ecokit::stop_ctx("Input metadata can not be `NULL`", metadata = metadata)
  }

  if (!file.exists(metadata$Path_Down)) {
    ecokit::stop_ctx("Input file does not exist", path = metadata$Path_Down)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # Environment variables -----
  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Loading reference grid -----

  GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", GridR = GridR)
  }
  GridR <- terra::unwrap(ecokit::load_as(GridR))

  # # ..................................................................... ###

  # Land mask ------

  # The land mask layer was sourced from the CHELSA-W5E5 v1.0 dataset, provided
  # as a NetCDF file. It was cropped to an extent encompassing the study area to
  # reduce file size. The resulting land mask layer is exported as LandMask.nc
  # and included as a data file within the package.
  #
  # source: CHELSA-W5E5 v1.0: W5E5 v1.0 downscaled with CHELSA v2.0
  # https://data.isimip.org/10.48364/ISIMIP.836809.3 Version: 1.0.3

  # Extent to crop the maps prior to processing. This ensures that the object
  # reads from the memory. See below
  CropExtent <- terra::ext(-26, 37.5, 34, 72)

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
  Rstr <- terra::rast(metadata$Path_Down, raw = TRUE) %>%
    stats::setNames(basename(metadata$Path_Down)) %>%
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
  if (metadata$scale != 1) {
    Rstr <- Rstr * metadata$scale
  }
  if (metadata$offset != 0) {
    Rstr <- Rstr + metadata$offset
  }

  # # ..................................................................... ###

  # Projecting to reference grid EPSG 3035 ------

  Rstr <- Rstr %>%
    # project to reference grid
    terra::project(GridR, method = "average", threads = TRUE) %>%
    # mask to the reference grid
    terra::mask(GridR) %>%
    # Ensure that values are read from memory
    ecokit::set_raster_values() %>%
    ecokit::set_raster_CRS()

  # # ..................................................................... ###

  # Write file to disk --- tiff -----

  terra::writeRaster(
    x = Rstr, filename = metadata$Path_Out_tif, overwrite = TRUE)

  # # ..................................................................... ###

  # Write file to disk --- nc -----

  # Variable name of the output *.nc file
  VarName4NC <- c(
    metadata$TimePeriod, metadata$ClimateModel, metadata$ClimateScenario) %>%
    unique() %>%
    paste(collapse = "__") %>%
    paste0(metadata$Variable, "__", .)

  # global attributes to be added to the *.nc file
  Attrs <- c(
    paste0("URL=", metadata$URL),
    paste0("OriginalFile=", metadata$Path_Down),
    paste0("Variable=", metadata$Variable),
    paste0("TimePeriod=", metadata$TimePeriod),
    paste0("ClimateModel=", metadata$ClimateModel),
    paste0("ClimateScenario=", metadata$ClimateScenario),
    paste0("Long_name=", metadata$Long_name),
    paste0("unit=", metadata$unit),
    paste0("explanation=", metadata$explanation))

  # save as *.nc file
  terra::writeCDF(
    Rstr,
    filename = metadata$Path_Out_NC, varname = VarName4NC,
    unit = metadata$unit, zname = metadata$TimePeriod,
    atts = Attrs, overwrite = TRUE, compression = compression_level)

  # # ..................................................................... ###

  return(invisible(NULL))
}
