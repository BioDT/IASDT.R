# # |------------------------------------------------------------------------| #
# soil_density_process ----
## |------------------------------------------------------------------------| #

#' Retrieve and project soil bulk density data
#'
#' Downloads and projects [SoilGrids](https://soilgrids.org)' bulk density of
#' the fine earth fraction (bdod) layers for specified depth intervals. See
#' [here](https://www.isric.org/explore/soilgrids/faq-soilgrids) for more
#' details. This function project the original data at different depths
#' intervals to the modelling reference grid.
#'
#' @param depths Character vector of depth intervals (e.g. `c("0-5","5-15")`).
#'   If `NULL` (default) all valid depths are processed. Valid `depths` are the
#'   SoilGrids BDOD standard horizons (cm): "0-5", "5-15", "15-30", "30-60",
#'   "60-100", and "100-200". Each depth string should omit the trailing `cm`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6.
#' @return Returns invisibly the path to `RData` file containing the processed
#'   soil bulk density data.
#' @references
#' - SoilGrids: <https://soilgrids.org>
#' - SoilGrids250m 2.0 - [Bulk density](https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/713396f9-1687-11ea-a7c0-a0481ca9e724)
#' - Poggio *et al.* (2021): SoilGrids 2.0: producing soil information for the
#' globe with quantified spatial uncertainty. Soil.
#' [10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021)
#' @author Ahmed El-Gabbas
#' @export

soil_density_process <- function(
    depths = NULL, env_file = ".env", n_cores = 6) {

  dir_grid <- dir_soil <- NULL

  # # ..................................................................... ###

  # Validate depths ----
  valid_depths <- c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200")

  if (is.null(depths)) {
    depths <- valid_depths
  } else if (!all(depths %in% valid_depths)) {
    ecokit::stop_ctx(
      "Input depths do not match valid depths",
      depths = depths, valid_depths = valid_depths)
  }

  # # ..................................................................... ###

  # Validate n_cores ----
  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      n_cores < 1L || is.na(n_cores)) {
    ecokit::stop_ctx(
      "n_cores must be a positive integer of length 1",
      n_cores = n_cores, class_n_cores = class(n_cores))
  }
  n_cores <- as.integer(n_cores)
  max_cores <- parallelly::availableCores()
  if (n_cores > max_cores) {
    warning(
      stringr::str_glue(
        "`n_cores` exceeds available cores: {n_cores}. Using all available",
        " cores: {max_cores}"),
      call. = FALSE)
    n_cores <- max_cores
  }

  # # ..................................................................... ###

  # Environment variables -----
  if (!fs::file_exists(env_file)) {
    ecokit::stop_ctx("Environment file does not exist", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "dir_grid", "DP_R_Grid_processed", TRUE, FALSE,
    "dir_soil", "DP_R_soil_density", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  soil_density_tif <- fs::path(dir_soil, "soil_density.tif")
  soil_density_RData <- fs::path(dir_soil, "soil_density.RData")

  # Check if topographic wetness index data already exists
  if (ecokit::check_data(soil_density_RData, warning = FALSE) &&
      ecokit::check_tiff(soil_density_tif, warning = FALSE)) {
    ecokit::cat_time("Soil bulk density data already exists")
    return(invisible(soil_density_RData))
  }

  fs::dir_create(dir_soil)

  # # ..................................................................... ###

  # Load the reference grid
  reference_grid <- fs::path(dir_grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(reference_grid)) {
    ecokit::stop_ctx(
      "Reference grid file does not exist",
      reference_grid = reference_grid)
  }
  reference_grid <- ecokit::load_as(reference_grid)

  # # ..................................................................... ###

  # Project soil bulk density rasters to the reference grid

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, length(depths)), show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  soil_bulk_density <- ecokit::quietly(
    future.apply::future_lapply(
      X = depths,
      FUN = function(depth) {

        soil_vrt <- paste0(
          "/vsicurl/",
          "https://files.isric.org/soilgrids/latest/data/bdod/bdod_",
          depth, "cm_mean.vrt") %>%
          terra::rast()

        if (!inherits(soil_vrt, "SpatRaster") || terra::nlyr(soil_vrt) != 1L ||
            terra::ncell(soil_vrt) == 0L) {
          ecokit::stop_ctx(
            "Failed to read soil bulk density VRT file",
            class_soil_vrt = class(soil_vrt), soil_vrt = soil_vrt)
        }

        terra::project(
          x = soil_vrt, y = terra::unwrap(reference_grid),
          method = "average") %>%
          terra::mask(terra::unwrap(reference_grid)) %>%
          stats::setNames(
            paste0(
              "bdod_", stringr::str_replace(depth, "-", "_"), "_mean")) %>%
          terra::wrap()
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = c("terra", "stringr", "dplyr"),
      future.globals = "reference_grid")) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast()

  ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
  future::plan("sequential", gc = TRUE)

  # # ..................................................................... ###

  terra::writeRaster(
    x = soil_bulk_density, filename = soil_density_tif, overwrite = TRUE)
  ecokit::save_as(
    object = terra::wrap(soil_bulk_density), object_name = "soil_bulk_density",
    out_path = soil_density_RData)

  return(invisible(soil_density_RData))
}
