## |------------------------------------------------------------------------| #
# bioreg_process ----
## |------------------------------------------------------------------------| #

#' Process biogeographical regions dataset
#'
#' Downloads and processes the Biogeographical Regions dataset (Europe 2016, v1)
#' from the [European Environment
#' Agency](https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618).
#' This function extracts biogeographical region names per reference grid cell
#' for use in counting species presence across biogeographical regions.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @details
#'   - *Temporal coverage*: 2011-2015
#'   - *Spatial coverage*: 28°E to 81°E, 31.27°W to 62°E
#'   - *CRS*: EPSG:3035
#'   - *file format*: shapefile (compressed in zip file)
#'   - *Requirements*: `curl` (download) and `unzip` (extraction)
#' @author Ahmed El-Gabbas
#' @name bioreg_process
#' @return Invisible `NULL`. Processed data is saved to disk as raster, vector,
#'   and RData files.
#' @export

bioreg_process <- function(env_file = ".env") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_bioreg <- path_raw <- path_interim <- bioreg_URL <- path_grid <-
    rast <- path_grid_raw <- geometry <- ID <- NULL

  # # ..................................................................... ###

  # Check input arguments ------
  ecokit::cat_time("Check input arguments")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  # # ..................................................................... ###

  # Check system commands ------

  if (isFALSE(ecokit::check_system_command("curl"))) {
    ecokit::stop_ctx(
      "`curl` is required for downloading data but was not found.",
      include_backtrace = TRUE)
  }
  if (isFALSE(ecokit::check_system_command("unzip"))) {
    ecokit::stop_ctx(
      "`unzip` is required for extracting data but was not found.",
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Environment variables
  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_grid_raw", "DP_R_grid_raw", TRUE, FALSE,
    "path_raw", "DP_R_bioreg_raw", FALSE, FALSE,
    "path_interim", "DP_R_bioreg_interim", FALSE, FALSE,
    "path_bioreg", "DP_R_bioreg_processed", FALSE, FALSE,
    "bioreg_URL", "DP_R_bioreg_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # Ensure necessary directories exist
  fs::dir_create(c(path_bioreg, path_raw, path_interim))

  # # ..................................................................... ###

  # Download biogeographical regions dataset ------

  ecokit::cat_time("Downloading biogeographical regions dataset")

  # name of downloaded zip file
  zip_file_name <- "Biogeog_regions_original.zip"

  # Extract download link
  ecokit::cat_time("Extract download link", level = 1L)
  bioreg_URL2 <- tryCatch({
    bioreg_URL %>%
      # extract download link
      httr::GET(config = httr::timeout(100L)) %>%
      rvest::read_html() %>%
      rvest::html_elements(css = "#334349 .list:nth-child(2) .content") %>%
      rvest::html_nodes("a") %>%
      rvest::html_attr("href") %>%
      # extract direct download link
      httr::GET(config = httr::timeout(100L)) %>%
      rvest::read_html() %>%
      rvest::html_elements(css = "#header-primary-action .button") %>%
      rvest::html_attr("href")
  }, error = function(e) {
    ecokit::stop_ctx(
      paste0("Failed to extract download URL: ", e$message),
      include_backtrace = TRUE)
  })

  if (length(bioreg_URL2) != 1L) {
    ecokit::stop_ctx(
      paste0("Download link extraction failed. Found: ", length(bioreg_URL2)),
      bioreg_URL2 = bioreg_URL2, length_bioreg_URL2 = length(bioreg_URL2),
      include_backtrace = TRUE)
  }
  ecokit::cat_time(bioreg_URL2, level = 2L, cat_timestamp = FALSE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Downloading using `curl`
  ecokit::cat_time("Download using `curl`", level = 1L)
  zip_file <- fs::path(path_raw, zip_file_name)
  DownCommand <- stringr::str_glue(
    'curl -J --create-dirs --output-dir {path_raw} -o\\
    "{zip_file_name}" -L {bioreg_URL2} --silent --max-time 300')

  attempt <- 1L
  repeat {
    ecokit::cat_time(
      paste0("Attempt ", attempt), level = 2L, cat_timestamp = FALSE)

    invisible(ecokit::system_command(DownCommand))

    Sys.sleep(5L)

    if (ecokit::check_zip(zip_file)) {
      break
    }

    if (attempt >= 5L) {
      ecokit::stop_ctx(
        "Error: Maximum download attempts reached. Zip file check failed.",
        include_backtrace = TRUE)
    }
    attempt <- attempt + 1L
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # unzip to interim directory
  ecokit::cat_time("unzip to interim directory", level = 1L)
  stringr::str_glue("unzip -o -qq -j {zip_file} -d {path_interim}") %>%
    ecokit::system_command() %>%
    invisible()

  # # ..................................................................... ###

  # Processing biogeographical regions data ------
  ecokit::cat_time("Processing biogeographical regions data")

  # Reading data from original shapefile
  ecokit::cat_time("Read data from original shapefile", level = 1L)
  bioreg_data <- fs::dir_ls(path = path_interim, type = "file", glob = "*.shp$")
  if (length(bioreg_data) != 1L) {
    ecokit::stop_ctx(
      paste0("Expected one .shp file, found: ", length(bioreg_data)),
      bioreg_data = bioreg_data, length_bioreg_data = length(bioreg_data),
      include_backtrace = TRUE)
  }
  bioreg_data <- sf::st_read(bioreg_data, quiet = TRUE) %>%
    # project to EPSG:3035
    sf::st_transform(3035L)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Extract metadata
  ecokit::cat_time("Extract metadata", level = 1L)
  bioreg_metadata <- sf::st_drop_geometry(bioreg_data) %>%
    tibble::as_tibble() %>%
    dplyr::rename(ID = "PK_UID") %>%
    list()

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Rasterize/masking
  ecokit::cat_time("Rasterize & masking", level = 1L)

  grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", grid_r = grid_r,
      include_backtrace = TRUE)
  }

  bioreg_r <- bioreg_data %>%
    dplyr::mutate(
      rast = purrr::map(
        .x = geometry,
        .f = ~{
          .x %>%
            sf::st_geometry() %>%
            sf::st_as_sf() %>%
            terra::rasterize(
              y = ecokit::load_as(grid_r, unwrap_r = TRUE), cover = TRUE) %>%
            terra::classify(cbind(NA, 0L))
        })) %>%
    dplyr::pull(rast) %>%
    terra::rast() %>%
    terra::which.max() %>%
    stats::setNames("ID") %>%
    terra::mask(ecokit::load_as(grid_r, unwrap_r = TRUE))

  rm(bioreg_data, grid_r, envir = environment())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Remove unused levels and adjust ID column
  ecokit::cat_time("Remove unused levels", level = 1L)
  levels(bioreg_r) <- bioreg_metadata
  bioreg_r <- terra::droplevels(bioreg_r)

  map_levels <- terra::levels(bioreg_r)[[1L]]
  map_levels_new <- dplyr::mutate(map_levels, ID = seq_len(dplyr::n()))
  map_levels_m <- map_levels %>%
    dplyr::left_join(map_levels_new, by = "short_name") %>%
    dplyr::select("short_name", tidyselect::everything())
  bioreg_r <- terra::classify(bioreg_r, map_levels_m[, -1L])
  levels(bioreg_r) <- list(map_levels_new)
  terra::crs(bioreg_r) <- "epsg:3035"

  ecokit::cat_time("Convert to sf object", level = 1L)
  grid_sf <- fs::path(path_grid_raw, "Grid_10_sf.RData") %>%
    ecokit::load_as() %>%
    magrittr::extract2("Grid_10_sf_s")
  bioreg_sf <- terra::as.polygons(
    x = bioreg_r, aggregate = FALSE, na.rm = TRUE) %>%
    sf::st_as_sf() %>%
    tibble::tibble() %>%
    sf::st_as_sf() %>%
    dplyr::left_join(bioreg_metadata[[1L]], by = "short_name") %>%
    sf::st_join(grid_sf) %>%
    dplyr::relocate(geometry, .after = tidyselect::everything())

  # # ..................................................................... ###

  # Saving processed data ----
  ecokit::cat_time("Saving processed data")

  ecokit::cat_time("tiff", level = 1L, cat_timestamp = FALSE)
  terra::writeRaster(
    x = bioreg_r, overwrite = TRUE,
    filename = file.path(path_bioreg, "bioreg_r.tif"))
  # Write attributes to file
  terra::levels(bioreg_r)[[1L]] %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = file.path(path_bioreg, "bioreg_r.tif.vat.dbf"),
      factor2char = TRUE, max_nchar = 254L)

  ecokit::cat_time("RData - raster object", level = 1L, cat_timestamp = FALSE)
  ecokit::save_as(
    object = terra::wrap(bioreg_r), object_name = "bioreg_r",
    out_path = file.path(path_bioreg, "bioreg_r.RData"))

  ecokit::cat_time("RData - sf object", level = 1L, cat_timestamp = FALSE)
  save(bioreg_sf, file = file.path(path_bioreg, "bioreg_sf.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Processing biogeographical regions data took ")

  return(invisible(NULL))
}
