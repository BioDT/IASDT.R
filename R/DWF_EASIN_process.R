#' Process EASIN data for the `IASDT`
#'
#' Extracts, processes, and visualises data from the [European Alien Species
#' Information Network (EASIN)](https://easin.jrc.ec.europa.eu/) for the
#' Invasive Alien Species Digital Twin (`IASDT`). Manages taxonomy, occurrence
#' data, and plots, handling API pagination and server limits. Orchestrated by
#' `EASIN_process()` with helpers `EASIN_taxonomy()`, `EASIN_download()`, and
#' `EASIN_plot()`.
#'
#' @param extract_taxa Logical. If `TRUE`, extracts taxonomy using
#'   `EASIN_taxonomy()`. Default: `TRUE`.
#' @param extract_data Logical.If `TRUE`, downloads occurrence data with
#'   `EASIN_download()`. Default: `TRUE`.
#' @param n_download_attempts Integer. Retry attempts for downloads. Default:
#'   `10`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6. The maximum number of allowed cores are 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param sleep_time Numeric. Seconds to wait between download attempts/chunks.
#'   Default: `10`.
#' @param n_search Integer. Records per taxonomy or data request (max 1000).
#'   Default: `1000`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
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
    n_cores = 6L, strategy = "multisession", sleep_time = 10L,
    n_search = 1000L, env_file = ".env", delete_chunks = TRUE,
    start_year = 1981L, plot = TRUE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")
  ecokit::check_args(
    args_to_check = c(
      "extract_taxa", "extract_data", "plot", "delete_chunks"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c(
      "n_download_attempts", "sleep_time", "n_search", "start_year"),
    args_type = "numeric")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  if (n_cores > 8) {
    message("Number of cores were reset from ", n_cores, " to 8")
    n_cores <- 8
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid <- taxa_info_file <- Name <- speciesKey <- EASINID <- taxon_name <-
    SpeciesId <- CellCode <- DataPartnerName <- Species_name <- Species_File <-
    EASIN_ref <- Year <- WKT <- path_EASIN <- path_EASIN_interim <- n <-
    path_grid_ref <- Points <- Longitude <- Latitude <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_Grid_processed", TRUE, FALSE,
    "path_grid_ref", "DP_R_Grid_raw", TRUE, FALSE,
    "path_EASIN", "DP_R_EASIN_processed", FALSE, FALSE,
    "path_EASIN_interim", "DP_R_EASIN_interim", FALSE, FALSE,
    "taxa_info_file", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "EASIN_ref", "DP_R_Taxa_easin", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "jsonlite", "purrr", "IASDT.R", "withr", "fs", "magrittr",
      "stringr", "RCurl", "tibble", "ecokit"),
    strategy = strategy)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # General input data + Paths ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Checking input and output paths")

  path_EASIN_data <- fs::path(path_EASIN, "species_data")
  path_EASIN_grid <- fs::path(path_EASIN, "species_grid")
  path_EASIN_r <- fs::path(path_EASIN, "species_r")
  path_EASIN_summary <- fs::path(path_EASIN, "EASIN_summary")

  fs::dir_create(
    c(
      path_EASIN_data, path_EASIN_grid, path_EASIN_r, path_EASIN_summary,
      path_EASIN_interim))

  # # ||||||||||||||||||||||||||||||||||||||||||||||

  ## Grid - raster ----
  grid_r <- fs::path(path_grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Path for the reference grid does not exist", grid_r = grid_r,
      include_backtrace = TRUE)
  }

  ## Grid - sf ----
  grid_sf <- fs::path(path_grid_ref, "Grid_10_sf.RData")
  if (!file.exists(grid_sf)) {
    ecokit::stop_ctx(
      "Path for the reference grid does not exist", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }

  ## Grid - sf - study area ----
  # Grid ID overlapping with study area
  land_grids <- fs::path(path_grid, "Grid_10_Land_sf.RData")
  if (!file.exists(land_grids)) {
    ecokit::stop_ctx(
      "Path for the reference grid does not exist", land_grids = land_grids,
      include_backtrace = TRUE)
  }

  ## Species list ----
  ecokit::cat_time("Loading species list", level = 1L)
  taxa_list <- ecokit::load_as(taxa_info_file)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Download EASIN taxonomy ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Extract EASIN taxonomy list")

  path_EASIN_taxa <- fs::path(path_EASIN, "EASIN_taxa.RData")
  taxa_okay <- ecokit::check_data(path_EASIN_taxa, warning = FALSE)

  if (extract_taxa || isFALSE(taxa_okay)) {
    ecokit::cat_time("Download EASIN taxa", level = 1L)

    # Download EASIN taxa
    EASIN_taxa_original <- IASDT.R::EASIN_taxonomy(
      kingdom = "Plantae", phylum = "Tracheophyta", n_search = n_search)
    save(
      EASIN_taxa_original,
      file = fs::path(path_EASIN, "EASIN_taxa_original.RData"))

    ecokit::cat_time("Loading pre-standardized EASIN taxonomy", level = 1L)
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

    EASIN_ref <- readRDS(EASIN_ref) %>%
      dplyr::select(Name, speciesKey, EASINID) %>%
      # Add missing species
      dplyr::bind_rows(
        tibble::tribble(
          ~Name, ~speciesKey, ~EASINID,
          "Cenchrus setaceus", 5828232, "R03000",
          "Neltuma juliflora", 5358460, "R12278",
          "Persicaria perfoliata", 4033648, "R19287",
          "Pueraria montana (Lour.) Merr. var. lobata", 2977636, "R12644"))

    EASIN_taxa <- EASIN_taxa_original %>%
      # Merge with EASIN reference list of taxonomy standardization
      dplyr::left_join(EASIN_ref, by = c("Name", "EASINID")) %>%
      dplyr::select(dplyr::all_of(c("EASINID", "Name", "speciesKey"))) %>%
      # Merge with final standardized taxa list
      dplyr::left_join(taxa_list, by = "speciesKey") %>%
      # Only keep species of interest
      dplyr::filter(!is.na(taxon_name)) %>%
      dplyr::select(-speciesKey) %>%
      dplyr::distinct() %>%
      # Extract species name (exclude authors)
      dplyr::mutate(Species_name = stringr::word(taxon_name, 1, 2)) %>%
      dplyr::rename(EASIN_Name = Name)

    ecokit::cat_time("Check EASIN taxa not in the list", level = 1L)
    # EASIN taxa not in the reference list
    new_EASIN_taxa <- setdiff(EASIN_taxa$EASIN_Name, EASIN_ref$Name)
    if (length(new_EASIN_taxa) > 0) {
      tibble::tibble(NewTaxa = new_EASIN_taxa) %>%
        readr::write_tsv(
          file = fs::path(path_EASIN, "new_EASIN_taxa.txt"),
          progress = FALSE)
    }
    ## Save EASIN taxa - RData ----
    ecokit::cat_time("Save EASIN taxa - RData", level = 1L)
    save(EASIN_taxa, file = path_EASIN_taxa)
  } else {
    ecokit::cat_time("Loading EASIN taxa list")
    EASIN_taxa <- ecokit::load_as(path_EASIN_taxa)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Download EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Download EASIN data")

  if (extract_data) {

    time_start_data <- lubridate::now(tzone = "CET")

    ## Prepare working in parallel ----
    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    ecokit::cat_time("Processing EASIN data", level = 1L)

    # Start downloading, allow for a maximum of `NumDownTries` trials
    ry_n <- 0

    # Start downloading ----
    repeat {
      try_n <- try_n + 1
      not_processed <- list.files(path_EASIN_interim, pattern = ".RData$") %>%
        stringr::str_remove(".RData") %>%
        setdiff(EASIN_taxa$EASINID, .)

      if (length(not_processed) == 0) {
        ecokit::cat_time("Data for all EASIN taxa were downloaded", level = 1L)
        break
      }

      ecokit::cat_time(paste0("Try number: ", ry_n), level = 1L)
      ecokit::cat_time(
        paste0(
          "There are ", length(not_processed), " EASIN taxa to be downloaded"),
        level = 2L)

      if (ry_n > n_download_attempts) {
        ecokit::cat_time(
          paste0(
            "Download failed for ", length(not_processed),
            " EASIN taxa after ", n_download_attempts, " download attempts"),
          level = 2L)
        break
      }

      data_download <- try(
        future.apply::future_lapply(
          X = not_processed, FUN = IASDT.R::EASIN_download, env_file = env_file,
          delete_chunks = delete_chunks, n_search = n_search,
          sleep_time = sleep_time,
          future.scheduling = Inf, future.seed = TRUE,
          future.packages = pkg_to_export,
          future.globals = c(
            "path_EASIN_interim", "n_search", "delete_chunks", "sleep_time",
            "env_file")),
        silent = TRUE)

      if (inherits(data_download, "try-error")) {
        next
      } else {
        break
      }

      rm(data_download, envir = environment())
      Sys.sleep(sleep_time)
    }

    ecokit::cat_diff(
      time_start_data,
      prefix = "Downloading EASIN data was finished in ", level = 1L)

    # Stopping cluster ----
    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Merging EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Merging EASIN data")

  ## Checking taxa with no data -----
  ecokit::cat_time("Checking taxa with no data", level = 1L)

  EASIN_files <- list.files(
    path_EASIN_interim, full.names = TRUE, pattern = ".RData")

  not_processed <- setdiff(
    paste0(EASIN_taxa$EASINID, ".RData"), basename(EASIN_files))

  if (length(not_processed) > 0) {
    path_not_processed <- fs::path(path_EASIN, "not_processedID.txt")
    stringr::str_remove_all(not_processed, ".RData") %>%
      cat(file = path_not_processed, sep = "\n")
    ecokit::cat_time(
      paste0(
        "There are ", length(not_processed), " not processed EASIN ID(s)."),
      level = 2L)

    ecokit::cat_time(
      paste0("EASIN IDs are saved to: ", path_not_processed),
      level = 2L)
  }

  ## Loading/merging EASIN data -----
  ecokit::cat_time("Loading or merging EASIN data", level = 1L)
  EASIN_data_original <- purrr::map_dfr(.x = EASIN_files, .f = ecokit::load_as)

  ## Save merged EASIN data - RData -----
  ecokit::cat_time("Save merged EASIN data - RData", level = 1L)
  save(
    EASIN_data_original,
    file = fs::path(path_EASIN, "EASIN_data_original.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Cleaning EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Cleaning EASIN data")
  EASIN_data <- EASIN_data_original %>%
    dplyr::mutate(Year = as.integer(Year), SpeciesName = NULL) %>%
    dplyr::rename(EASINID = SpeciesId) %>%
    # exclude observations with no spatial information or < start_year
    dplyr::filter(!is.na(WKT), Year >= start_year) %>%
    # Join with EASIN Taxa information
    dplyr::left_join(EASIN_taxa, by = "EASINID")

  ## Extract coordinates from WKT string ----
  ecokit::cat_time("Extract coordinates from WKT string", level = 1L)

  WKTs <- dplyr::distinct(EASIN_data, WKT) %>%
    dplyr::mutate(
      Points = purrr::map(
        .x = WKT,
        .f = ~ {
          # Extract POINT coordinates from WKT string
          Points <- stringr::str_extract_all(
            .x, "POINT\\s*\\(\\s*-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+\\s*\\)")[[1]]

          if (length(Points) > 0) {
            purrr::map(Points, ecokit::text_to_coordinates) %>%
              dplyr::bind_rows() %>%
              dplyr::mutate(
                dplyr::across(
                  .cols = c("Longitude", "Latitude"), .fns = ~ round(.x, 5)))
          } else {
            tibble::tibble(Longitude = NA_real_, Latitude = NA_real_)
          }
        }))

  ## Add coordinates to data and convert to sf ----
  ecokit::cat_time("Add coordinates to data and convert to sf", level = 1L)
  grid_sf <- ecokit::load_as(grid_sf) %>%
    magrittr::extract2("Grid_10_sf_s")
  land_grids <- ecokit::load_as(land_grids) %>%
    sf::st_drop_geometry() %>%
    dplyr::pull("CellCode")

  EASIN_data <- dplyr::left_join(EASIN_data, WKTs, by = "WKT") %>%
    tidyr::unnest(Points) %>%
    dplyr::filter(!is.na(Longitude) & !is.na(Latitude)) %>%
    # convert to sf object, while keeping original coordinates as columns
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    # spatial transform the data to EPSG:3035
    sf::st_transform(3035) %>%
    # Add grid ID
    sf::st_join(grid_sf) %>%
    # exclude observations not overlapping with the study area
    dplyr::filter(CellCode %in% land_grids)

  rm(land_grids, WKTs, envir = environment())

  ## Save cleaned EASIN Data - RData ----
  ecokit::cat_time("Save cleaned EASIN Data - RData", level = 1L)
  save(EASIN_data, file = fs::path(path_EASIN, "EASIN_data.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Summarizing EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Summarizing EASIN data")

  grid_r <- ecokit::load_as(grid_r, unwrap_r = TRUE)

  ## n_obs ----
  ecokit::cat_time("Number of observations per grid cell", level = 1L)

  ecokit::cat_time("All data", level = 2L)
  EASIN_n_obs <- sf::st_drop_geometry(EASIN_data) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(n_obs = n) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(grid_r, field = "n_obs") %>%
    terra::mask(grid_r) %>%
    terra::wrap()

  Path_n_obs <- fs::path(path_EASIN_summary, "EASIN_n_obs.RData")
  save(EASIN_n_obs, file = Path_n_obs)

  ### n_obs per partner ----
  ecokit::cat_time("Number of observations per partner", level = 2L)

  EASIN_n_obs_partner <- sf::st_drop_geometry(EASIN_data) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(n_obs = n) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ ecokit::replace_space(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = grid_r, field = "n_obs") %>%
          terra::mask(grid_r)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  Path_n_obs_partner <- fs::path(
    path_EASIN_summary, "EASIN_n_obs_partner.RData")
  save(EASIN_n_obs_partner, file = Path_n_obs_partner)

  ## NSpecies ----
  ecokit::cat_time("Number of species per grid cell", level = 1L)

  ecokit::cat_time("All data", level = 2L)
  EASIN_n_sp <- sf::st_drop_geometry(EASIN_data) %>%
    dplyr::distinct(CellCode, Species_name) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(n_sp = n) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(grid_r, field = "n_sp") %>%
    terra::mask(grid_r) %>%
    terra::wrap()

  Path_NSp <- fs::path(path_EASIN_summary, "EASIN_n_sp.RData")
  save(EASIN_n_sp, file = Path_NSp)

  ### n_sp per partner ----
  ecokit::cat_time("Number of species per partner", level = 2L)

  EASIN_n_sp_partner <- sf::st_drop_geometry(EASIN_data) %>%
    dplyr::distinct(CellCode, Species_name, DataPartnerName) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(n_sp = n) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ ecokit::replace_space(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = grid_r, field = "n_sp") %>%
          terra::mask(grid_r)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  path_n_sp_partner <- fs::path(
    path_EASIN_summary, "EASIN_n_sp_partner.RData")
  save(EASIN_n_sp_partner, file = path_n_sp_partner)

  if (delete_chunks) {
    fs::dir_delete(fs::path(path_EASIN_interim, "FileParts"))
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  ## Species-specific data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Species-specific data")

  time_start_data <- lubridate::now(tzone = "CET")

  EASIN_data_2 <- dplyr::group_by(EASIN_data, Species_File) %>%
    dplyr::group_split() %>%
    purrr::set_names(
      purrr::map_chr(., ~ ecokit::replace_space(.x$Species_File[1]))
    )

  seq_len(length(EASIN_data_2)) %>%
    purrr::walk(
      .f = ~ {
        # # ||||||||||||||||||||||||||||||||||
        # Data as sf object
        # # ||||||||||||||||||||||||||||||||||
        ecokit::save_as(
          object = EASIN_data_2[[.x]],
          object_name = names(EASIN_data_2)[.x],
          out_path = fs::path(
            path_EASIN_data, paste0(names(EASIN_data_2)[.x], "_data.RData"))
        )

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - sf
        # # ||||||||||||||||||||||||||||||||||
        data_grid <- sf::st_drop_geometry(EASIN_data_2[[.x]]) %>%
          dplyr::count(
            EASINID, taxon_name, Species_name, Species_name2,
            Species_File, CellCode) %>%
          dplyr::rename(n_obs = n) %>%
          dplyr::left_join(grid_sf, by = "CellCode") %>%
          sf::st_as_sf()

        ecokit::save_as(
          object = data_grid,
          object_name = paste0(names(EASIN_data_2)[.x], "_grid"),
          out_path = fs::path(
            path_EASIN_grid, paste0(names(EASIN_data_2)[.x], "_grid.RData")))

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - raster
        # # ||||||||||||||||||||||||||||||||||
        data_r <- terra::rasterize(data_grid, grid_r, field = "n_obs")
        data_r$PA <- terra::as.int(data_r > 0)
        ecokit::save_as(
          object = terra::wrap(data_r),
          object_name = paste0(names(EASIN_data_2)[.x], "_r"),
          out_path = fs::path(
            path_EASIN_r, paste0(names(EASIN_data_2)[.x], "_r.RData")))

        return(invisible(NULL))
      }, .progress = FALSE)

  rm(EASIN_data_2, envir = environment())

  ecokit::cat_diff(
    time_start_data,
    prefix = "Preparing Species-specific data was finished in ", level = 1L)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  if (plot) {
    ecokit::cat_time("Plotting")
    IASDT.R::EASIN_plot(env_file = env_file)
  }

  # # ..................................................................... ###

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing EASIN data was finished in ")

  return(invisible(NULL))
}
