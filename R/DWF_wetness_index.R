# # |------------------------------------------------------------------------| #
# wetness_index_process ----
## |------------------------------------------------------------------------| #

#' Download and Process Topographic Wetness Index Data
#'
#' Downloads, extracts, and processes the global topographic wetness index data
#' at 30 arc-second resolution (Title & Bemmels, 2018). The function checks for
#' existing processed data, downloads the required dataset if necessary,
#' extracts the relevant TIFF file, re-projects and masks it to match a
#' reference grid, and saves the processed raster in both TIFF and `RData`
#' formats.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @return (Invisibly) The path to the processed wetness index `RData` file.
#' @references
#' - Title & Bemmels (2018): ENVIREM: an expanded set of bioclimatic and
#' topographic variables increases flexibility and improves performance of
#' ecological niche modeling. Ecography.
#' [10.1111/ecog.02880](https://doi.org/10.1111/ecog.02880)
#' - [ENVIREM](https://deepblue.lib.umich.edu/data/concern/generic_works/gt54kn05f): ENVIronmental Rasters for Ecological Modeling version 1.0
#'
#' @export
#' @author Ahmed El-Gabbas

wetness_index_process <- function(env_file = ".env") {

  dir_grid <- dir_wetness_raw <- dir_wetness <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  if (!fs::file_exists(env_file)) {
    ecokit::stop_ctx("Environment file does not exist", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "dir_grid", "DP_R_Grid_processed", TRUE, FALSE,
    "dir_wetness_raw", "DP_R_wetness_raw", FALSE, FALSE,
    "dir_wetness", "DP_R_wetness_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Load the reference grid
  reference_grid <- fs::path(dir_grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(reference_grid)) {
    ecokit::stop_ctx(
      "Reference grid file does not exist", reference_grid = reference_grid)
  }
  reference_grid <- ecokit::load_as(reference_grid)

  # # ..................................................................... ###

  fs::dir_create(c(dir_wetness_raw, dir_wetness))
  wetness_file_name <- "elev_global_current_30arcsec_geotiff.zip"
  wetness_index_tif <- fs::path(dir_wetness, "wetness_index.tif")
  wetness_index_RData <- fs::path(dir_wetness, "wetness_index.RData")
  topo_wet_file_name <- "current_30arcsec_topoWet.tif"
  wetness_raw_file <- fs::path(dir_wetness_raw, wetness_file_name)

  # Check if topographic wetness index data already exists
  if (ecokit::check_data(wetness_index_RData, warning = FALSE) &&
      ecokit::check_tiff(wetness_index_tif, warning = FALSE)) {
    ecokit::cat_time("Topographic wetness index data already exists")
    return(invisible(wetness_index_RData))
  }

  # # ..................................................................... ###

  # Download data -----

  if (!ecokit::check_zip(wetness_raw_file, warning = FALSE)) {

    if (fs::file_exists(wetness_raw_file)) fs::file_delete(wetness_raw_file)

    page_url <- paste0(
      "https://deepblue.lib.umich.edu/data/", "concern/generic_works/gt54kn05f")

    if (!ecokit::check_url(page_url, timeout = 10)) {
      ecokit::stop_ctx(
        "The source link does not exist or is not accessible",
        page_url = page_url)
    }

    x_path <- paste0(
      "//a[contains(text(), '",
      wetness_file_name, "') or contains(@href, '", wetness_file_name, "')]")

    # extract zip link for topographic wetness index
    zip_link <- rvest::read_html(page_url) %>%
      rvest::html_elements("a") %>%
      rvest::html_elements(xpath = x_path) %>%
      rvest::html_attr("href")

    if (length(zip_link) != 1L) {
      ecokit::stop_ctx(
        "Expected exactly one download link for the topographic wetness index",
        length = length(zip_link))
    }

    zip_link <- dplyr::if_else(
      startsWith(zip_link, "http"),
      zip_link, paste0("https://deepblue.lib.umich.edu", zip_link))
    if (!ecokit::check_url(zip_link, timeout = 10)) {
      ecokit::stop_ctx(
        "The extracted zip link does not exist or is not accessible",
        zip_link = zip_link)
    }

    wetness_download <- httr::GET(
      url = zip_link,
      config = httr::write_disk(wetness_raw_file, overwrite = TRUE),
      httr::timeout(1200))

    if (httr::status_code(wetness_download) != 200) {
      ecokit::stop_ctx(
        "Failed to download the topographic wetness index data",
        status_code = httr::status_code(wetness_download),
        zip_link = zip_link, wetness_raw_file = wetness_raw_file)
    }

    if (!ecokit::check_zip(wetness_raw_file, warning = FALSE)) {
      ecokit::stop_ctx(
        "The downloaded file is not a valid zip file or is corrupted",
        wetness_raw_file = wetness_raw_file)
    }
  }

  # # ..................................................................... ###

  # Extract topographic wetness index tiff file -----
  suppressMessages(
    suppressWarnings(
      archive::archive_extract(
        archive = wetness_raw_file, dir = dir_wetness_raw,
        files = topo_wet_file_name)))

  topo_wet_file <- fs::path(dir_wetness_raw, topo_wet_file_name)
  if (!ecokit::check_tiff(topo_wet_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "The extracted file is not a valid TIFF file or is corrupted",
      topo_wet_file = topo_wet_file)
  }

  # # ..................................................................... ###

  # Process and save data -----

  wetness_index <- terra::rast(topo_wet_file) %>%
    terra::project(terra::unwrap(reference_grid), method = "average") %>%
    terra::mask(terra::unwrap(reference_grid)) %>%
    stats::setNames("wetness_index")

  terra::writeRaster(
    x = wetness_index, filename = wetness_index_tif, overwrite = TRUE)
  ecokit::save_as(
    object = terra::wrap(wetness_index), object_name = "wetness_index",
    out_path = wetness_index_RData)

  return(invisible(wetness_index_RData))

}
