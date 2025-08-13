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
  Path_BioReg <- Path_Raw <- Path_Interim <- BioReg_URL <- Path_Grid <-
    Rast <- Path_Grid_Ref <- geometry <- ID <- NULL

  # # ..................................................................... ###

  # Check input arguments ------
  ecokit::cat_time("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  ecokit::check_args(
    args_all = AllArgs, args_type = "character", args_to_check = "env_file")
  rm(AllArgs, envir = environment())

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
  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_Raw", "DP_R_BioReg_raw", FALSE, FALSE,
    "Path_Interim", "DP_R_BioReg_interim", FALSE, FALSE,
    "Path_BioReg", "DP_R_BioReg_processed", FALSE, FALSE,
    "BioReg_URL", "DP_R_BioReg_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Ensure necessary directories exist
  fs::dir_create(c(Path_BioReg, Path_Raw, Path_Interim))

  # # ..................................................................... ###

  # Download biogeographical regions dataset ------

  ecokit::cat_time("Downloading biogeographical regions dataset")

  # name of downloaded zip file
  zip_file_name <- "Biogeog_regions_original.zip"

  # Extract download link
  ecokit::cat_time("Extract download link", level = 1L)
  BioReg_URL2 <- tryCatch({
    BioReg_URL %>%
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

  if (length(BioReg_URL2) != 1L) {
    ecokit::stop_ctx(
      paste0("Download link extraction failed. Found: ", length(BioReg_URL2)),
      BioReg_URL2 = BioReg_URL2, length_BioReg_URL2 = length(BioReg_URL2),
      include_backtrace = TRUE)
  }
  ecokit::cat_time(BioReg_URL2, level = 2L, cat_timestamp = FALSE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Downloading using `curl`
  ecokit::cat_time("Download using `curl`", level = 1L)
  Zip_file <- fs::path(Path_Raw, zip_file_name)
  DownCommand <- stringr::str_glue(
    'curl -J --create-dirs --output-dir {Path_Raw} -o\\
    "{zip_file_name}" -L {BioReg_URL2} --silent --max-time 300')

  attempt <- 1L
  repeat {
    ecokit::cat_time(
      paste0("Attempt ", attempt), level = 2L, cat_timestamp = FALSE)

    invisible(ecokit::system_command(DownCommand))

    Sys.sleep(5L)

    if (ecokit::check_zip(Zip_file)) {
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
  stringr::str_glue("unzip -o -qq -j {Zip_file} -d {Path_Interim}") %>%
    ecokit::system_command() %>%
    invisible()

  # # ..................................................................... ###

  # Processing biogeographical regions data ------
  ecokit::cat_time("Processing biogeographical regions data")

  # Reading data from original shapefile
  ecokit::cat_time("Read data from original shapefile", level = 1L)
  BioReg_DT <- fs::dir_ls(path = Path_Interim, type = "file", glob = "*.shp$")
  if (length(BioReg_DT) != 1L) {
    ecokit::stop_ctx(
      paste0("Expected one .shp file, found: ", length(BioReg_DT)),
      BioReg_DT = BioReg_DT, length_BioReg_DT = length(BioReg_DT),
      include_backtrace = TRUE)
  }
  BioReg_DT <- sf::st_read(BioReg_DT, quiet = TRUE) %>%
    # project to EPSG:3035
    sf::st_transform(3035L)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Extract metadata
  ecokit::cat_time("Extract metadata", level = 1L)
  BioReg_Metadata <- sf::st_drop_geometry(BioReg_DT) %>%
    tibble::as_tibble() %>%
    dplyr::rename(ID = "PK_UID") %>%
    list()

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Rasterize/masking
  ecokit::cat_time("Rasterize & masking", level = 1L)

  GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", GridR = GridR,
      include_backtrace = TRUE)
  }

  BioReg_R <- BioReg_DT %>%
    dplyr::mutate(
      Rast = purrr::map(
        .x = geometry,
        .f = ~{
          .x %>%
            sf::st_geometry() %>%
            sf::st_as_sf() %>%
            terra::rasterize(
              y = ecokit::load_as(GridR, unwrap_r = TRUE), cover = TRUE) %>%
            terra::classify(cbind(NA, 0L))
        })) %>%
    dplyr::pull(Rast) %>%
    terra::rast() %>%
    terra::which.max() %>%
    stats::setNames("ID") %>%
    terra::mask(ecokit::load_as(GridR, unwrap_r = TRUE))

  rm(BioReg_DT, GridR, envir = environment())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Remove unused levels and adjust ID column
  ecokit::cat_time("Remove unused levels", level = 1L)
  levels(BioReg_R) <- BioReg_Metadata
  BioReg_R <- terra::droplevels(BioReg_R)

  MapLevels <- terra::levels(BioReg_R)[[1L]]
  MapLevelsNew <- dplyr::mutate(MapLevels, ID = seq_len(dplyr::n()))
  MapLevelsM <- MapLevels %>%
    dplyr::left_join(MapLevelsNew, by = "short_name") %>%
    dplyr::select("short_name", tidyselect::everything())
  BioReg_R <- terra::classify(BioReg_R, MapLevelsM[, -1L])
  levels(BioReg_R) <- list(MapLevelsNew)
  terra::crs(BioReg_R) <- "epsg:3035"

  ecokit::cat_time("Convert to sf object", level = 1L)
  Grid_sf <- fs::path(Path_Grid_Ref, "Grid_10_sf.RData") %>%
    ecokit::load_as() %>%
    magrittr::extract2("Grid_10_sf_s")
  BioReg_sf <- terra::as.polygons(
    x = BioReg_R, aggregate = FALSE, na.rm = TRUE) %>%
    sf::st_as_sf() %>%
    tibble::tibble() %>%
    sf::st_as_sf() %>%
    dplyr::left_join(BioReg_Metadata[[1L]], by = "short_name") %>%
    sf::st_join(Grid_sf) %>%
    dplyr::relocate(geometry, .after = tidyselect::everything())

  # # ..................................................................... ###

  # Saving processed data ----
  ecokit::cat_time("Saving processed data")

  ecokit::cat_time("tiff", level = 1L, cat_timestamp = FALSE)
  terra::writeRaster(
    x = BioReg_R, overwrite = TRUE,
    filename = file.path(Path_BioReg, "BioReg_R.tif"))
  # Write attributes to file
  terra::levels(BioReg_R)[[1L]] %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = file.path(Path_BioReg, "BioReg_R.tif.vat.dbf"),
      factor2char = TRUE, max_nchar = 254L)

  ecokit::cat_time("RData - raster object", level = 1L, cat_timestamp = FALSE)
  ecokit::save_as(
    object = terra::wrap(BioReg_R), object_name = "BioReg_R",
    out_path = file.path(Path_BioReg, "BioReg_R.RData"))

  ecokit::cat_time("RData - sf object", level = 1L, cat_timestamp = FALSE)
  save(BioReg_sf, file = file.path(Path_BioReg, "BioReg_sf.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Processing biogeographical regions data took ")

  return(invisible(NULL))
}
