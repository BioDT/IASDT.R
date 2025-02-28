## |------------------------------------------------------------------------| #
# BioReg_Process ----
## |------------------------------------------------------------------------| #

#' Process biogeographical regions dataset
#'
#' Downloads and processes the Biogeographical Regions dataset (Europe 2016, v1)
#' from the [European Environment
#' Agency](https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618).
#' This function extracts biogeographical region names per reference grid cell
#' for use in counting species presence across biogeographical regions.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical. Whether the processing is being done on an
#'   High-Performance Computing (HPC) environment, to adjust file paths
#'   accordingly. Default: `TRUE`.
#' @details
#'   - *Temporal coverage*: 2011-2015
#'   - *Spatial coverage*: 28°E to 81°E, 31.27°W to 62°E
#'   - *CRS*: EPSG:3035
#'   - *file format*: shapefile (compressed in zip file)
#'   - *Requirements*: `curl` (download) and `unzip` (extraction)
#' @author Ahmed El-Gabbas
#' @name BioReg_Process
#' @return Invisible `NULL`. Processed data is saved to disk as raster, vector,
#'   and RData files.
#' @export

BioReg_Process <- function(FromHPC = TRUE, EnvFile = ".env") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_BioReg <- Path_Raw <- Path_Interim <- BioReg_URL <- Path_Grid <-
    Rast <- Path_Grid_Ref <- geometry <- ID <- NULL

  # # ..................................................................... ###

  # Check input arguments ------
  IASDT.R::CatTime("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  rm(AllArgs, envir = environment())

  if (isFALSE(IASDT.R::CheckCommands("curl"))) {
    stop(
      "`curl` is required for downloading data but was not found.",
      call. = FALSE)
  }
  if (isFALSE(IASDT.R::CheckCommands("unzip"))) {
    stop(
      "`unzip` is required for extracting data but was not found.",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Environment variables
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_Raw", "DP_R_BioReg_Raw", FALSE, FALSE,
      "Path_Interim", "DP_R_BioReg_Interim", FALSE, FALSE,
      "Path_BioReg", "DP_R_BioReg", FALSE, FALSE,
      "BioReg_URL", "DP_R_BioReg_URL", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_Raw", "DP_R_BioReg_Raw_Local", FALSE, FALSE,
      "Path_Interim", "DP_R_BioReg_Interim_Local", FALSE, FALSE,
      "Path_BioReg", "DP_R_BioReg_Local", FALSE, FALSE,
      "BioReg_URL", "DP_R_BioReg_URL", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # Ensure necessary directories exist
  fs::dir_create(c(Path_BioReg, Path_Raw, Path_Interim))

  # # ..................................................................... ###

  # Download biogeographical regions dataset ------

  IASDT.R::CatTime("Downloading biogeographical regions dataset")

  # name of downloaded zip file
  ZipFileName <- "Biogeog_regions_original.zip"

  # Extract download link
  IASDT.R::CatTime("Extract download link", Level = 1)
  BioReg_URL2 <- tryCatch({
    BioReg_URL %>%
      # extract download link
      httr::GET(config = httr::timeout(100)) %>%
      rvest::read_html() %>%
      rvest::html_elements(css = "#334349 .list:nth-child(2) .content") %>%
      rvest::html_nodes("a") %>%
      rvest::html_attr("href") %>%
      # extract direct download link
      httr::GET(config = httr::timeout(100)) %>%
      rvest::read_html() %>%
      rvest::html_elements(css = "#header-primary-action .button") %>%
      rvest::html_attr("href")
  }, error = function(e) {
    stop("Failed to extract download URL: ", e$message, call. = FALSE)
  })

  if (length(BioReg_URL2) != 1) {
    stop(
      "Download link extraction failed. Found: ", length(BioReg_URL2),
      call. = FALSE)
  }
  IASDT.R::CatTime(BioReg_URL2, Level = 2, Time = FALSE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Downloading using `curl`
  IASDT.R::CatTime("Download using `curl`", Level = 1)
  Zip_file <- IASDT.R::Path(Path_Raw, ZipFileName)
  DownCommand <- stringr::str_glue(
    'curl -J --create-dirs --output-dir {Path_Raw} -o\\
    "{ZipFileName}" -L {BioReg_URL2} --silent --max-time 300')

  attempt <- 1
  repeat {
    IASDT.R::CatTime(paste0("Attempt ", attempt), Level = 2, Time = FALSE)

    invisible(IASDT.R::System(DownCommand))

    if (IASDT.R::CheckZip(Zip_file)) {
      break
    }

    if (attempt >= 5) {
      stop(
        "Error: Maximum download attempts reached. Zip file check failed.",
        call. = FALSE)
    }
    attempt <- attempt + 1
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # unzip to interim directory
  IASDT.R::CatTime("unzip to interim directory", Level = 1)
  stringr::str_glue("unzip -o -qq -j {Zip_file} -d {Path_Interim}") %>%
    IASDT.R::System() %>%
    invisible()

  # # ..................................................................... ###

  # Processing biogeographical regions data ------
  IASDT.R::CatTime("Processing biogeographical regions data")

  # Reading data from original shapefile
  IASDT.R::CatTime("Read data from original shapefile", Level = 1)
  BioReg_DT <- fs::dir_ls(path = Path_Interim, type = "file", glob = "*.shp$")
  if (length(BioReg_DT) != 1) {
    stop("Expected one .shp file, found: ", length(BioReg_DT), call. = FALSE)
  }
  BioReg_DT <- sf::st_read(BioReg_DT, quiet = TRUE) %>%
    # project to EPSG:3035
    sf::st_transform(3035)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Extract metadata
  IASDT.R::CatTime("Extract metadata", Level = 1)
  BioReg_Metadata <- sf::st_drop_geometry(BioReg_DT) %>%
    tibble::as_tibble() %>%
    dplyr::rename(ID = "PK_UID") %>%
    list()

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Rasterize/masking
  IASDT.R::CatTime("Rasterize/masking", Level = 1)

  GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop(
      "Path for the Europe boundaries does not exist: ", GridR, call. = FALSE)
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
              y = terra::unwrap(IASDT.R::LoadAs(GridR)), cover = TRUE) %>%
            terra::classify(cbind(NA, 0))
        })) %>%
    dplyr::pull(Rast) %>%
    terra::rast() %>%
    terra::which.max() %>%
    stats::setNames("ID") %>%
    terra::mask(terra::unwrap(IASDT.R::LoadAs(GridR)))

  rm(BioReg_DT, GridR, envir = environment())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Remove unused levels and adjust ID column
  IASDT.R::CatTime("Remove unused levels", Level = 1)
  levels(BioReg_R) <- BioReg_Metadata
  BioReg_R <- terra::droplevels(BioReg_R)

  MapLevels <- terra::levels(BioReg_R)[[1]]
  MapLevelsNew <- dplyr::mutate(MapLevels, ID = seq_len(dplyr::n()))
  MapLevelsM <- MapLevels %>%
    dplyr::left_join(MapLevelsNew, by = "short_name") %>%
    dplyr::select("short_name", tidyselect::everything())
  BioReg_R <- terra::classify(BioReg_R, MapLevelsM[, -1])
  levels(BioReg_R) <- list(MapLevelsNew)
  terra::crs(BioReg_R) <- "epsg:3035"

  IASDT.R::CatTime("Convert to sf object", Level = 1)
  Grid_sf <- IASDT.R::Path(Path_Grid_Ref, "Grid_10_sf.RData") %>%
    IASDT.R::LoadAs() %>%
    magrittr::extract2("Grid_10_sf_s")
  BioReg_sf <- terra::as.polygons(
    x = BioReg_R, aggregate = FALSE, na.rm = TRUE) %>%
    sf::st_as_sf() %>%
    tibble::tibble() %>%
    sf::st_as_sf() %>%
    dplyr::left_join(BioReg_Metadata[[1]], by = "short_name") %>%
    sf::st_join(Grid_sf) %>%
    dplyr::relocate(geometry, .after = tidyselect::everything())

  # # ..................................................................... ###

  # Saving processed data ----
  IASDT.R::CatTime("Saving processed data")

  IASDT.R::CatTime("tiff", Level = 1, Time = FALSE)
  terra::writeRaster(
    x = BioReg_R, overwrite = TRUE,
    filename = file.path(Path_BioReg, "BioReg_R.tif"))
  # Write attributes to file
  terra::levels(BioReg_R)[[1]] %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = file.path(Path_BioReg, "BioReg_R.tif.vat.dbf"),
      factor2char = TRUE, max_nchar = 254)

  IASDT.R::CatTime("RData - raster object", Level = 1, Time = FALSE)
  IASDT.R::SaveAs(
    InObj = terra::wrap(BioReg_R), OutObj = "BioReg_R",
    OutPath = file.path(Path_BioReg, "BioReg_R.RData"))

  IASDT.R::CatTime("RData - sf object", Level = 1, Time = FALSE)
  save(BioReg_sf, file = file.path(Path_BioReg, "BioReg_sf.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "Processing biogeographical regions data took ")

  return(invisible(NULL))
}
