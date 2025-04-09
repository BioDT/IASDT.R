#' Process EASIN data for the `IAS-pDT`
#'
#' Extracts, processes, and visualizes data from the [European Alien Species
#' Information Network (EASIN)](https://easin.jrc.ec.europa.eu/) for the
#' Invasive Alien Species prototype Digital Twin (`IAS-pDT`). Manages taxonomy,
#' occurrence data, and plots, handling API pagination and server limits.
#' Orchestrated by `EASIN_process()` with helpers `EASIN_taxonomy()`,
#' `EASIN_download()`, and `EASIN_plot()`.
#'
#' @param extract_taxa Logical. If `TRUE`, extracts taxonomy using
#'   `EASIN_taxonomy()`. Default: `TRUE`.
#' @param extract_data Logical.If `TRUE`, downloads occurrence data with
#'   `EASIN_download()`. Default: `TRUE`.
#' @param n_download_attempts Integer. Retry attempts for downloads. Default:
#'   `10`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6. The maximum number of allowed cores are 8.
#' @param sleep_time Numeric. Seconds to wait between download attempts/chunks.
#'   Default: `10`.
#' @param n_search Integer. Records per taxonomy or data request (max 1000).
#'   Default: `1000`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param delete_chunks Logical. If `TRUE`, removes intermediate files. Default:
#'   `FALSE`.
#' @param start_year Integer. Earliest year for occurrence data (excludes
#'   earlier records). Default: `1981` (aligned with CHELSA climate data).
#' @param plot Logical. If `TRUE`, generates plots via `EASIN_plot()`. Default:
#'   `TRUE`.
#' @param species_key Character. EASIN taxon ID for which data is to be
#'   retrieved. This parameter cannot be `NULL`.
#' @param timeout Integer. Download timeout in seconds. Default: `200`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default:
#'   `FALSE`.
#' @param n_search Integer. Number of records to attempt to retrieve per
#'   request. Default: 1000, which is the current maximum allowed by the API.
#' @param n_attempts Integer. Max download attempts per chunk. Default: `10`.
#' @param sleep_time Integer. Number of seconds to pause between each data
#'   retrieval request to prevent overloading the server. Default: 5 second.
#' @param delete_chunks Logical. Whether to delete temporary files for data
#'   chunks from the `FileParts` subdirectory. Defaults to `TRUE`.
#' @param return_data Logical. If `TRUE`, returns data as a dataframe;
#'   otherwise, saves to disk and returns `invisible(NULL)`. Default: `FALSE`.
#' @param kingdom Character. Taxonomic kingdom to query. Default: `"Plantae"`.
#' @param phylum Character. Taxonomic phylum within kingdom. Default:
#'   `"Tracheophyta"`
#' @note Uses a static RDS file with EASIN-GBIF taxonomic standardization,
#'   prepared by Marina Golivets (Feb 2024).
#' @section Functions details:
#' - **`EASIN_process()`**: Orchestrates taxonomy extraction, data downloads,
#'   and plotting for EASIN species data.
#' - **`EASIN_taxonomy()`**: Fetches taxonomy data in chunks via the EASIN API,
#'   filtered by kingdom and phylum. Returns a tibble.
#' - **`EASIN_download()`**: Downloads occurrence data for a given EASIN ID,
#'   handling pagination and pauses. Returns a dataframe if `return_data =
#'   TRUE`, else `invisible(NULL)`.
#' - **`EASIN_plot()`**: Creates summary plots (observations count, species
#'   count, distribution by partner) as JPEGs. Returns `invisible(NULL)`.

# # |------------------------------------------------------------------------| #
# EASIN_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name EASIN_data
#' @rdname EASIN_data
#' @export
#' @order 1

EASIN_process <- function(
    extract_taxa = TRUE, extract_data = TRUE, n_download_attempts = 10L,
    n_cores = 6L, sleep_time = 10L, n_search = 1000L, env_file = ".env",
    delete_chunks = TRUE, start_year = 1981L, plot = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::cat_time("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character", args_to_check = "env_file")
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "extract_taxa", "extract_data", "verbose", "plot", "delete_chunks"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c(
      "DownTries", "n_cores", "sleep_time", "n_search", "start_year"))

  if (n_cores > 8) {
    message("Number of cores were reset from ", n_cores, " to 8")
    n_cores <- 8
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- TaxaInfoFile <- Name <- speciesKey <- EASINID <- taxon_name <-
    SpeciesId <- CellCode <- DataPartnerName <- Species_name <- Species_File <-
    EASIN_Ref <- Year <- WKT <- Path_EASIN <- Path_EASIN_Interim <- n <-
    Path_Grid_Ref <- Points <- Longitude <- Latitude <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_EASIN", "DP_R_EASIN_processed", FALSE, FALSE,
    "Path_EASIN_Interim", "DP_R_EASIN_interim", FALSE, FALSE,
    "TaxaInfoFile", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "EASIN_Ref", "DP_R_Taxa_easin", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # General input data + Paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Checking input and output paths")

  Path_EASIN_DT <- IASDT.R::path(Path_EASIN, "Sp_DT")
  Path_EASIN_Grid <- IASDT.R::path(Path_EASIN, "Sp_Grid")
  Path_EASIN_R <- IASDT.R::path(Path_EASIN, "Sp_R")
  Path_EASIN_Summary <- IASDT.R::path(Path_EASIN, "Summary")

  fs::dir_create(
    c(
      Path_EASIN_DT, Path_EASIN_Grid, Path_EASIN_R, Path_EASIN_Summary,
      Path_EASIN_Interim))

  # # ||||||||||||||||||||||||||||||||||||||||||||||

  ## Grid - raster ----
  GridR <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop("Path for the reference grid does not exist: ", GridR, call. = FALSE)
  }

  ## Grid - sf ----
  GridSf <- IASDT.R::path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Path for the reference grid does not exist: ", GridSf, call. = FALSE)
  }

  ## Grid - sf - study area ----
  # Grid ID overlapping with study area
  LandGrids <- IASDT.R::path(Path_Grid, "Grid_10_Land_sf.RData")
  if (!file.exists(LandGrids)) {
    stop(
      "Path for the reference grid does not exist: ", LandGrids, call. = FALSE)
  }

  ## Species list ----
  IASDT.R::cat_time("Loading species list", level = 1)
  TaxaList <- IASDT.R::load_as(TaxaInfoFile)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Download EASIN taxonomy ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Extract EASIN taxonomy list")

  Path_EASIN_Taxa <- IASDT.R::path(Path_EASIN, "EASIN_Taxa.RData")
  Taxa_Okay <- IASDT.R::check_data(Path_EASIN_Taxa, warning = FALSE)

  if (extract_taxa || isFALSE(Taxa_Okay)) {
    IASDT.R::cat_time("Download EASIN taxa", level = 1)

    # Download EASIN taxa
    EASIN_Taxa_Orig <- IASDT.R::EASIN_taxonomy(
      kingdom = "Plantae", phylum = "Tracheophyta", n_search = n_search)
    save(
      EASIN_Taxa_Orig,
      file = IASDT.R::path(Path_EASIN, "EASIN_Taxa_Orig.RData"))

    IASDT.R::cat_time("Loading pre-standardized EASIN taxonomy", level = 1)
    # Update 08.2024: While writing this function, I noticed that there 4 taxa
    # in EASIN taxonomy list that is not matched at the most recent
    # pre-standardization (2024-02-07). These four names are temporarily
    # hard-coded in this function until the next version of the taxa
    # standardization.

    # tibble::tribble(
    #   ~EASIN_New, ~EASIN_Matched, ~EASIN_ID, ~speciesKey,
    #   "Cenchrus setaceus", "Pennisetum setaceum", "R03000", 5828232,
    #   "Neltuma juliflora", "Prosopis juliflora", "R12278", 5358460,
    #   "Persicaria perfoliata", "Polygonum perfoliatum", "R19287", 4033648,
    #   "Pueraria montana (Lour.) Merr. var. lobata",
    #    "Pueraria montana var. lobata", "R12644", 2977636)

    EASIN_Ref <- readRDS(EASIN_Ref) %>%
      dplyr::select(Name, speciesKey, EASINID) %>%
      # Add missing species
      dplyr::bind_rows(
        tibble::tribble(
          ~Name, ~speciesKey, ~EASINID,
          "Cenchrus setaceus", 5828232, "R03000",
          "Neltuma juliflora", 5358460, "R12278",
          "Persicaria perfoliata", 4033648, "R19287",
          "Pueraria montana (Lour.) Merr. var. lobata", 2977636, "R12644"))

    EASIN_Taxa <- EASIN_Taxa_Orig %>%
      # Merge with EASIN reference list of taxonomy standardization
      dplyr::left_join(EASIN_Ref, by = c("Name", "EASINID")) %>%
      dplyr::select(dplyr::all_of(c("EASINID", "Name", "speciesKey"))) %>%
      # Merge with final standardized taxa list
      dplyr::left_join(TaxaList, by = "speciesKey") %>%
      # Only keep species of interest
      dplyr::filter(!is.na(taxon_name)) %>%
      dplyr::select(-speciesKey) %>%
      dplyr::distinct() %>%
      # Extract species name (exclude authors)
      dplyr::mutate(Species_name = stringr::word(taxon_name, 1, 2)) %>%
      dplyr::rename(EASIN_Name = Name)

    IASDT.R::cat_time("Check EASIN taxa not in the list", level = 1)
    # EASIN taxa not in the reference list
    New_EASIN_Taxa <- setdiff(EASIN_Taxa$EASIN_Name, EASIN_Ref$Name)
    if (length(New_EASIN_Taxa) > 0) {
      tibble::tibble(NewTaxa = New_EASIN_Taxa) %>%
        readr::write_tsv(
          file = IASDT.R::path(Path_EASIN, "New_EASIN_Taxa.txt"),
          progress = FALSE)
    }
    ## Save EASIN taxa - RData ----
    IASDT.R::cat_time("Save EASIN taxa - RData", level = 1)
    save(EASIN_Taxa, file = Path_EASIN_Taxa)
  } else {
    IASDT.R::cat_time("Loading EASIN taxa list")
    EASIN_Taxa <- IASDT.R::load_as(Path_EASIN_Taxa)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Download EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Download EASIN data")

  if (extract_data) {

    TimeStartData <- lubridate::now(tzone = "CET")

    ## Prepare working on parallel ----

    IASDT.R::cat_time(
      paste0("Prepare working on parallel using ", n_cores, " cores"),
      level = 1)

    if (n_cores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(n_cores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      withr::defer(future::plan("future::sequential", gc = TRUE))
    }

    IASDT.R::cat_time("Processing EASIN data", level = 1)

    # Start downloading, allow for a maximum of `NumDownTries` trials
    Try <- 0

    # Start downloading ----
    repeat {
      Try <- Try + 1
      NotProcessed <- list.files(Path_EASIN_Interim, pattern = ".RData$") %>%
        stringr::str_remove(".RData") %>%
        setdiff(EASIN_Taxa$EASINID, .)

      if (length(NotProcessed) == 0) {
        IASDT.R::cat_time("Data for all EASIN taxa were downloaded", level = 1)
        break
      }

      IASDT.R::cat_time(paste0("Try number: ", Try), level = 1)
      IASDT.R::cat_time(
        paste0(
          "There are ", length(NotProcessed), " EASIN taxa to be downloaded"),
        level = 2)

      if (Try > n_download_attempts) {
        IASDT.R::cat_time(
          paste0(
            "Download failed for ", length(NotProcessed),
            " EASIN taxa after ", n_download_attempts, " download attempts"),
          level = 2)
        break
      }

      Down <- try(
        future.apply::future_lapply(
          X = NotProcessed, FUN = IASDT.R::EASIN_download, env_file = env_file,
          delete_chunks = delete_chunks, n_search = n_search,
          sleep_time = sleep_time,
          future.scheduling = Inf, future.seed = TRUE,
          future.globals = c(
            "Path_EASIN_Interim", "n_search", "delete_chunks", "sleep_time",
            "env_file"),
          future.packages = c(
            "dplyr", "jsonlite", "purrr", "IASDT.R", "withr", "fs",
            "stringr", "RCurl", "tibble")),
        silent = TRUE)

      if (inherits(Down, "try-error")) {
        next
      } else {
        break
      }

      rm(Down, envir = environment())
      Sys.sleep(sleep_time)
    }

    IASDT.R::cat_diff(
      TimeStartData,
      prefix = "Downloading EASIN data was finished in ", level = 1)

    # Stopping cluster ----
    IASDT.R::cat_time("Stopping cluster", level = 1)

    if (n_cores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Merging EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Merging EASIN data")

  ## Checking taxa with no data -----
  IASDT.R::cat_time("Checking taxa with no data", level = 1)

  EASIN_Files <- list.files(
    Path_EASIN_Interim, full.names = TRUE, pattern = ".RData")

  NotProcessed <- setdiff(
    paste0(EASIN_Taxa$EASINID, ".RData"), basename(EASIN_Files))

  if (length(NotProcessed) > 0) {
    Path_NotProcessed <- IASDT.R::path(Path_EASIN, "NotProcessedID.txt")
    stringr::str_remove_all(NotProcessed, ".RData") %>%
      cat(file = Path_NotProcessed, sep = "\n")
    IASDT.R::cat_time(
      paste0(
        "There are ", length(NotProcessed), " not processed EASIN ID(s)."),
      level = 2)

    IASDT.R::cat_time(
      paste0("EASIN IDs are saved to: ", Path_NotProcessed),
      level = 2)
  }

  ## Loading/merging EASIN data -----
  IASDT.R::cat_time("Loading/merging EASIN data", level = 1)
  EASIN_Data_Orig <- purrr::map_dfr(.x = EASIN_Files, .f = IASDT.R::load_as)

  ## Save merged EASIN data - RData -----
  IASDT.R::cat_time("Save merged EASIN data - RData", level = 1)
  save(
    EASIN_Data_Orig, file = IASDT.R::path(Path_EASIN, "EASIN_Data_Orig.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Cleaning EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Cleaning EASIN data")
  EASIN_Data <- EASIN_Data_Orig %>%
    dplyr::mutate(Year = as.integer(Year), SpeciesName = NULL) %>%
    dplyr::rename(EASINID = SpeciesId) %>%
    # exclude observations with no spatial information or < start_year
    dplyr::filter(!is.na(WKT), Year >= start_year) %>%
    # Join with EASIN Taxa information
    dplyr::left_join(EASIN_Taxa, by = "EASINID")

  ## Extract coordinates from WKT string ----
  IASDT.R::cat_time("Extract coordinates from WKT string", level = 1)

  WKTs <- dplyr::distinct(EASIN_Data, WKT) %>%
    dplyr::mutate(
      Points = purrr::map(
        .x = WKT,
        .f = ~ {
          # Extract POINT coordinates from WKT string
          Points <- stringr::str_extract_all(
            .x, "POINT\\s*\\(\\s*-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+\\s*\\)")[[1]]

          if (length(Points) > 0) {
            purrr::map(Points, IASDT.R::text_to_coordinates) %>%
              dplyr::bind_rows() %>%
              dplyr::mutate(
                dplyr::across(
                  .cols = c("Longitude", "Latitude"), .fns = ~ round(.x, 5)))
          } else {
            tibble::tibble(Longitude = NA_real_, Latitude = NA_real_)
          }
        }))

  ## Add coordinates to data and convert to sf ----
  IASDT.R::cat_time("Add coordinates to data and convert to sf", level = 1)
  GridSf <- IASDT.R::load_as(GridSf) %>%
    magrittr::extract2("Grid_10_sf_s")
  LandGrids <- IASDT.R::load_as(LandGrids) %>%
    sf::st_drop_geometry() %>%
    dplyr::pull("CellCode")

  EASIN_Data <- dplyr::left_join(EASIN_Data, WKTs, by = "WKT") %>%
    tidyr::unnest(Points) %>%
    dplyr::filter(!is.na(Longitude) & !is.na(Latitude)) %>%
    # convert to sf object, while keeping original coordinates as columns
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    # spatial transform the data to EPSG:3035
    sf::st_transform(3035) %>%
    # Add grid ID
    sf::st_join(GridSf) %>%
    # exclude observations not overlapping with the study area
    dplyr::filter(CellCode %in% LandGrids)

  rm(LandGrids, WKTs, envir = environment())

  ## Save cleaned EASIN Data - RData ----
  IASDT.R::cat_time("Save cleaned EASIN Data - RData", level = 1)
  save(EASIN_Data, file = IASDT.R::path(Path_EASIN, "EASIN_Data.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Summarizing EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Summarizing EASIN data")

  GridR <- terra::unwrap(IASDT.R::load_as(GridR))

  ## NObs ----
  IASDT.R::cat_time("Number of observations per grid cell", level = 1)

  IASDT.R::cat_time("All data", level = 2)
  EASIN_NObs <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(NObs = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(GridR, field = "NObs") %>%
    terra::mask(GridR) %>%
    terra::wrap()

  Path_NObs <- IASDT.R::path(Path_EASIN_Summary, "EASIN_NObs.RData")
  save(EASIN_NObs, file = Path_NObs)

  ### NObs per partner ----
  IASDT.R::cat_time("Number of observations per partner", level = 2)

  EASIN_NObs_PerPartner <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(NObs = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::replace_space(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = GridR, field = "NObs") %>%
          terra::mask(GridR)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NObs_PerPartner <- IASDT.R::path(
    Path_EASIN_Summary, "EASIN_NObs_PerPartner.RData")
  save(EASIN_NObs_PerPartner, file = Path_NObs_PerPartner)

  ## NSpecies ----
  IASDT.R::cat_time("Number of species per grid cell", level = 1)

  IASDT.R::cat_time("All data", level = 2)
  EASIN_NSp <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::distinct(CellCode, Species_name) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(NSp = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(GridR, field = "NSp") %>%
    terra::mask(GridR) %>%
    terra::wrap()

  Path_NSp <- IASDT.R::path(Path_EASIN_Summary, "EASIN_NSp.RData")
  save(EASIN_NSp, file = Path_NSp)

  ### NSp per partner ----
  IASDT.R::cat_time("Number of species per per partner", level = 2)

  EASIN_NSp_PerPartner <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::distinct(CellCode, Species_name, DataPartnerName) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(NSp = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::replace_space(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = GridR, field = "NSp") %>%
          terra::mask(GridR)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NSp_PerPartner <- IASDT.R::path(
    Path_EASIN_Summary, "EASIN_NSp_PerPartner.RData")
  save(EASIN_NSp_PerPartner, file = Path_NSp_PerPartner)

  if (delete_chunks) {
    fs::dir_delete(IASDT.R::path(Path_EASIN_Interim, "FileParts"))
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  ## Species-specific data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Species-specific data")

  TimeStartData <- lubridate::now(tzone = "CET")

  EASIN_Data2 <- dplyr::group_by(EASIN_Data, Species_File) %>%
    dplyr::group_split() %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::replace_space(.x$Species_File[1]))
    )

  seq_len(length(EASIN_Data2)) %>%
    purrr::walk(
      .f = ~ {
        # # ||||||||||||||||||||||||||||||||||
        # Data as sf object
        # # ||||||||||||||||||||||||||||||||||
        IASDT.R::save_as(
          object = EASIN_Data2[[.x]],
          object_name = names(EASIN_Data2)[.x],
          out_path = IASDT.R::path(
            Path_EASIN_DT, paste0(names(EASIN_Data2)[.x], "_DT.RData"))
        )

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - sf
        # # ||||||||||||||||||||||||||||||||||
        DT_Grid <- sf::st_drop_geometry(EASIN_Data2[[.x]]) %>%
          dplyr::count(
            EASINID, taxon_name, Species_name, Species_name2,
            Species_File, CellCode) %>%
          dplyr::rename(NObs = n) %>%
          dplyr::left_join(GridSf, by = "CellCode") %>%
          sf::st_as_sf()

        IASDT.R::save_as(
          object = DT_Grid,
          object_name = paste0(names(EASIN_Data2)[.x], "_Grid"),
          out_path = IASDT.R::path(
            Path_EASIN_Grid, paste0(names(EASIN_Data2)[.x], "_Grid.RData")))

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - raster
        # # ||||||||||||||||||||||||||||||||||
        DT_R <- terra::rasterize(DT_Grid, GridR, field = "NObs")
        DT_R$PA <- terra::as.int(DT_R > 0)
        IASDT.R::save_as(
          object = terra::wrap(DT_R),
          object_name = paste0(names(EASIN_Data2)[.x], "_R"),
          out_path = IASDT.R::path(
            Path_EASIN_R, paste0(names(EASIN_Data2)[.x], "_R.RData")))

        return(invisible(NULL))
      }, .progress = FALSE)

  rm(EASIN_Data2, envir = environment())

  IASDT.R::cat_diff(
    TimeStartData,
    prefix = "Preparing Species-specific data was finished in ", level = 1)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  if (plot) {
    IASDT.R::cat_time("Plotting")
    IASDT.R::EASIN_plot(env_file = env_file)
  }

  # # ..................................................................... ###

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nProcessing EASIN data was finished in ")

  return(invisible(NULL))
}
