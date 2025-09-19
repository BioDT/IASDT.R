## |------------------------------------------------------------------------| #
# chelsa_project ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name chelsa_data
#' @rdname chelsa_data
#' @order 3

chelsa_project <- function(
    metadata = NULL, env_file = ".env", compression_level = 5L) {

  # Checking input arguments -----

  ecokit::check_args(args_to_check = "compression_level", args_type = "numeric")

  if (is.null(metadata)) {
    ecokit::stop_ctx(
      "Input metadata can not be `NULL`", metadata = metadata,
      include_backtrace = TRUE)
  }

  if (!inherits(metadata, "tbl_df")) {
    ecokit::stop_ctx(
      "Input metadata has to be a tibble", class_metadata = class(metadata),
      include_backtrace = TRUE)
  }

  if (nrow(metadata) != 1) {
    ecokit::stop_ctx(
      "Input metadata has to be a single-row tibble",
      metadata = metadata, class_metadata = class(metadata),
      nrow_metadata = nrow(metadata), include_backtrace = TRUE)
  }

  # Check metadata columns
  needed_columns <- c(
    "scale", "offset", "path_out_tif", "climate_model", "climate_scenario",
    "variable", "url", "path_download", "long_name", "explanation",
    "path_out_netcdf", "unit", "time_period")
  missing_columns <- setdiff(needed_columns, names(metadata))
  if (length(missing_columns) > 0) {
    ecokit::stop_ctx(
      "Input metadata is missing columns", missing_columns = missing_columns,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid <- NULL

  # # ..................................................................... ###

  if (!file.exists(metadata$path_download)) {
    ecokit::stop_ctx(
      "Input file does not exist", path = metadata$path_download,
      include_backtrace = TRUE)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # Environment variables -----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # Loading reference grid -----

  grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", grid_r = grid_r,
      include_backtrace = TRUE)
  }
  grid_r <- ecokit::load_as(grid_r, unwrap_r = TRUE)

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
  crop_extent <- terra::ext(-26, 37.5, 34, 72)

  land_mask_l <- system.file(
    "extdata", "LandMask.nc", package = "IASDT.R", mustWork = TRUE) %>%
    terra::rast() %>%
    terra::crop(crop_extent) %>%
    suppressWarnings() %>% # suppress warning on LUMI while cropping
    terra::classify(cbind(0, NA))

  # # ..................................................................... ###

  # Read tif file as terra SpatRaster object ----

  # terra package by default considers the scale and offset information stored
  # in the tiff files. Here I disable this to read the raw values as it is and
  # later consider the scale and offset information manually. This is more safe
  # as I found that some of the future projections do not include such
  # information in the tiff files.
  r_map <- terra::rast(metadata$path_download, raw = TRUE) %>%
    stats::setNames(basename(metadata$path_download)) %>%
    # crop to European boundaries although it is not necessary to crop the input
    # maps into the European boundaries, we will crop the data prior to
    # projection. Cropping will make the values of the raster read from memory
    # not from the file. This is a workaround to avoid wrong extreme values in
    # the output file because of a bug in terra package (see this issue:
    # https://github.com/rspatial/terra/issues/1356) [18.02.2023]
    terra::crop(crop_extent) %>%
    # mask by land mask
    terra::mask(land_mask_l) %>%
    # `gsp` maps contains extremely high values instead of NA; the following
    # replace extreme values with NA
    terra::classify(cbind(420000000, Inf, NA))

  rm(land_mask_l, crop_extent, envir = environment())

  # # ..................................................................... ###

  # Manually considering offset and scale -----

  # For `npp` layers, all tiff maps except for current climate does have a
  # scaling factor all scale and offset information were set manually
  if (metadata$scale != 1) {
    r_map <- r_map * metadata$scale
  }
  if (metadata$offset != 0) {
    r_map <- r_map + metadata$offset
  }

  # # ..................................................................... ###

  # Projecting to reference grid EPSG 3035 ------

  r_map <- r_map %>%
    # project to reference grid
    terra::project(grid_r, method = "average", threads = TRUE) %>%
    # mask to the reference grid
    terra::mask(grid_r) %>%
    # Ensure that values are read from memory
    terra::toMemory() %>%
    ecokit::set_raster_crs(crs = "epsg:3035")

  # # ..................................................................... ###

  # Write file to disk --- tiff -----

  terra::writeRaster(
    x = r_map, filename = metadata$path_out_tif, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

  # # ..................................................................... ###

  # Write file to disk --- nc -----

  # variable name of the output *.nc file
  var_name_for_nc <- c(
    metadata$time_period, metadata$climate_model, metadata$climate_scenario) %>%
    unique() %>%
    paste(collapse = "__") %>%
    paste0(metadata$variable, "__", .)

  # global attributes to be added to the *.nc file
  attributes <- c(
    paste0("url=", metadata$url),
    paste0("original_file=", metadata$path_download),
    paste0("variable=", metadata$variable),
    paste0("time_period=", metadata$time_period),
    paste0("climate_model=", metadata$climate_model),
    paste0("climate_scenario=", metadata$climate_scenario),
    paste0("long_name=", metadata$long_name),
    paste0("unit=", metadata$unit),
    paste0("explanation=", metadata$explanation))

  # save as *.nc file
  terra::writeCDF(
    r_map,
    filename = metadata$path_out_netcdf, varname = var_name_for_nc,
    unit = metadata$unit, zname = metadata$time_period,
    atts = attributes, overwrite = TRUE, compression = compression_level)

  # # ..................................................................... ###

  return(invisible(NULL))
}
