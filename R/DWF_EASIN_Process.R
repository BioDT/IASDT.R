#' Process EASIN data for the `IAS-pDT`
#'
#' Extracts, processes, and visualizes data from the [European Alien Species
#' Information Network (EASIN)](https://easin.jrc.ec.europa.eu/) for the
#' Invasive Alien Species prototype Digital Twin (`IAS-pDT`). Manages taxonomy,
#' occurrence data, and plots, handling API pagination and server limits.
#' Orchestrated by `EASIN_Process()` with helpers `EASIN_Taxonomy()`,
#' `EASIN_Down()`, and `EASIN_Plot()`.
#'
#' @param ExtractTaxa Logical. If `TRUE`, extracts taxonomy using
#'   `EASIN_Taxonomy()`. Default: `TRUE`.
#' @param ExtractData Logical.If `TRUE`, downloads occurrence data with
#'   `EASIN_Down()`. Default: `TRUE`.
#' @param NDownTries Integer. Retry attempts for downloads. Default: `10`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6. The maximum number of allowed cores are 8.
#' @param SleepTime Numeric. Seconds to wait between download attempts/chunks.
#'   Default: `10`.
#' @param NSearch Integer. Records per taxonomy or data request (max 1000).
#'   Default: `1000`.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param DeleteChunks Logical. If `TRUE`, removes intermediate files. Default:
#'   `FALSE`.
#' @param StartYear Integer. Earliest year for occurrence data (excludes earlier
#'   records). Default: `1981` (aligned with CHELSA climate data).
#' @param Plot Logical. If `TRUE`, generates plots via `EASIN_Plot()`. Default:
#'   `TRUE`.
#' @param SpKey Character. EASIN taxon ID for which data is to be retrieved.
#'   This parameter cannot be `NULL`.
#' @param Timeout Integer. Download timeout in seconds. Default: `200`.
#' @param Verbose Logical. If `TRUE`, prints progress messages. Default:
#'   `FALSE`.
#' @param NSearch Integer. Number of records to attempt to retrieve per request.
#'   Default: 1000, which is the current maximum allowed by the API.
#' @param Attempts Integer. Max download attempts per chunk. Default: `10`.
#' @param SleepTime Integer. Number of seconds to pause between each data
#'   retrieval request to prevent overloading the server. Default: 5 second.
#' @param DeleteChunks Logical. Whether to delete temporary files for data
#'   chunks from the `FileParts` subdirectory. Defaults to `TRUE`.
#' @param ReturnData Logical. If `TRUE`, returns data as a dataframe; otherwise,
#'   saves to disk and returns `invisible(NULL)`. Default: `FALSE`.
#' @param Kingdom Character. Taxonomic kingdom to query. Default: `"Plantae"`.
#' @param Phylum Character. Taxonomic phylum within kingdom. Default:
#'   `"Tracheophyta"`
#' @note Uses a static RDS file with EASIN-GBIF taxonomic standardization,
#'   prepared by Marina Golivets (Feb 2024).
#' @section Functions details:
#' - **`EASIN_Process()`**: Orchestrates taxonomy extraction, data downloads,
#'   and plotting for EASIN species data.
#' - **`EASIN_Taxonomy()`**: Fetches taxonomy data in chunks via the EASIN API,
#'   filtered by kingdom and phylum. Returns a tibble.
#' - **`EASIN_Down()`**: Downloads occurrence data for a given EASIN ID,
#'   handling pagination and pauses. Returns a dataframe if `ReturnData = TRUE`,
#'   else `invisible(NULL)`.
#' - **`EASIN_Plot()`**: Creates summary plots (observations count, species
#'   count, distribution by partner) as JPEGs. Returns `invisible(NULL)`.

# # |------------------------------------------------------------------------| #
# EASIN_Process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 1

EASIN_Process <- function(
    ExtractTaxa = TRUE, ExtractData = TRUE, NDownTries = 10L, NCores = 6L,
    SleepTime = 10L, NSearch = 1000L, EnvFile = ".env", DeleteChunks = TRUE,
    StartYear = 1981L, Plot = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("ExtractTaxa", "ExtractData", "Verbose", "Plot", "DeleteChunks"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("DownTries", "NCores", "SleepTime", "NSearch", "StartYear"))

  if (NCores > 8) {
    message("Number of cores were reset from ", NCores, " to 8")
    NCores <- 8
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- TaxaInfoFile <- Name <- speciesKey <- EASINID <- taxon_name <-
    SpeciesId <- CellCode <- DataPartnerName <- Species_name <- Species_File <-
    EASIN_Ref <- Year <- WKT <- Path_EASIN <- Path_EASIN_Interim <-
    n <- Path_Grid_Ref <- Coords <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_EASIN", "DP_R_EASIN_processed", FALSE, FALSE,
    "Path_EASIN_Interim", "DP_R_EASIN_interim", FALSE, FALSE,
    "EASIN_URL", "DP_R_EASIN_url", FALSE, FALSE,
    "TaxaInfoFile", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "EASIN_Ref", "DP_R_Taxa_easin", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # General input data + Paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking input and output paths")

  Path_EASIN_DT <- IASDT.R::Path(Path_EASIN, "Sp_DT")
  Path_EASIN_Grid <- IASDT.R::Path(Path_EASIN, "Sp_Grid")
  Path_EASIN_R <- IASDT.R::Path(Path_EASIN, "Sp_R")
  Path_EASIN_Summary <- IASDT.R::Path(Path_EASIN, "Summary")

  fs::dir_create(
    c(
      Path_EASIN_DT, Path_EASIN_Grid, Path_EASIN_R, Path_EASIN_Summary,
      Path_EASIN_Interim))

  # # ||||||||||||||||||||||||||||||||||||||||||||||

  ## Grid - raster ----
  GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop("Path for the reference grid does not exist: ", GridR, call. = FALSE)
  }

  ## Grid - sf ----
  GridSf <- IASDT.R::Path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Path for the reference grid does not exist: ", GridSf, call. = FALSE)
  }

  ## Grid - sf - study area ----
  # Grid ID overlapping with study area
  LandGrids <- IASDT.R::Path(Path_Grid, "Grid_10_Land_sf.RData")
  if (!file.exists(LandGrids)) {
    stop(
      "Path for the reference grid does not exist: ", LandGrids, call. = FALSE)
  }

  ## Species list ----
  IASDT.R::CatTime("Loading species list", Level = 1)
  TaxaList <- IASDT.R::LoadAs(TaxaInfoFile)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Download EASIN taxonomy ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Extract EASIN taxonomy list")

  if (ExtractTaxa) {
    IASDT.R::CatTime("Download EASIN taxa", Level = 1)

    # Download EASIN taxa
    EASIN_Taxa_Orig <- IASDT.R::EASIN_Taxonomy(
      Kingdom = "Plantae", Phylum = "Tracheophyta", NSearch = NSearch)
    save(
      EASIN_Taxa_Orig,
      file = IASDT.R::Path(Path_EASIN, "EASIN_Taxa_Orig.RData"))

    IASDT.R::CatTime("Loading pre-standardized EASIN taxonomy", Level = 1)
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
          "Pueraria montana (Lour.) Merr. var. lobata", 2977636, "R12644")
      )

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

    IASDT.R::CatTime("Check EASIN taxa not in the list", Level = 1)
    # EASIN taxa not in the reference list
    New_EASIN_Taxa <- setdiff(EASIN_Taxa$EASIN_Name, EASIN_Ref$Name)
    if (length(New_EASIN_Taxa) > 0) {
      tibble::tibble(NewTaxa = New_EASIN_Taxa) %>%
        readr::write_tsv(
          file = IASDT.R::Path(Path_EASIN, "New_EASIN_Taxa.txt"),
          progress = FALSE)
    }
    ## Save EASIN taxa - RData ----
    IASDT.R::CatTime("Save EASIN taxa - RData", Level = 1)
    save(EASIN_Taxa, file = IASDT.R::Path(Path_EASIN, "EASIN_Taxa.RData"))
  } else {
    IASDT.R::CatTime("Loading EASIN taxa list")
    load(IASDT.R::Path(Path_EASIN, "EASIN_Taxa.RData"))
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Download EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Download EASIN data")

  if (ExtractData) {
    TimeStartData <- lubridate::now(tzone = "CET")

    ## Prepare working on parallel ----

    IASDT.R::CatTime(
      paste0("Prepare working on parallel using `", NCores, "` cores."),
      Level = 1)

    if (NCores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    # Start downloading, allow for a maximum of `NumDownTries` trials
    Try <- 0

    # Start downloading ----
    repeat {
      Try <- Try + 1
      NotProcessed <- list.files(Path_EASIN_Interim, pattern = ".RData$") %>%
        stringr::str_remove(".RData") %>%
        setdiff(EASIN_Taxa$EASINID, .)

      if (length(NotProcessed) == 0) {
        IASDT.R::CatTime("Data for all EASIN taxa were downloaded", Level = 1)
        break
      }

      IASDT.R::CatTime(paste0("Try number: ", Try), Level = 1)
      IASDT.R::CatTime(
        paste0(
          "There are ", length(NotProcessed), " EASIN taxa to be downloaded"),
        Level = 2)

      if (Try > NDownTries) {
        IASDT.R::CatTime(
          paste0(
            "Download failed for ", length(NotProcessed),
            " EASIN taxa after ", NDownTries, " download attempts"),
          Level = 2)
        break
      }

      IASDT.R::CatTime("Processing EASIN data", Level = 1)
      Down <- try(
        future.apply::future_lapply(
          X = NotProcessed, FUN = IASDT.R::EASIN_Down, EnvFile = EnvFile,
          DeleteChunks = DeleteChunks,
          NSearch = NSearch, SleepTime = SleepTime,
          future.scheduling = Inf, future.seed = TRUE,
          future.globals = c(
            "Path_EASIN_Interim", "NSearch", "DeleteChunks", "SleepTime",
            "EnvFile"),
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
      Sys.sleep(SleepTime)
    }

    IASDT.R::CatDiff(
      TimeStartData,
      Prefix = "Downloading EASIN data was finished in ", Level = 1)

    # Stopping cluster ----
    IASDT.R::CatTime("Stopping cluster", Level = 1)

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Merging EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Merging EASIN data")

  ## Checking taxa with no data -----
  IASDT.R::CatTime("Checking taxa with no data", Level = 1)

  EASIN_Files <- list.files(
    Path_EASIN_Interim, full.names = TRUE, pattern = ".RData")

  NotProcessed <- setdiff(
    paste0(EASIN_Taxa$EASINID, ".RData"), basename(EASIN_Files))

  if (length(NotProcessed) > 0) {
    Path_NotProcessed <- IASDT.R::Path(Path_EASIN, "NotProcessedID.txt")
    stringr::str_remove_all(NotProcessed, ".RData") %>%
      cat(file = Path_NotProcessed, sep = "\n")
    IASDT.R::CatTime(
      paste0(
        "There are ", length(NotProcessed), " not processed EASIN ID(s)."),
      Level = 2)

    IASDT.R::CatTime(
      paste0("EASIN IDs are saved to: ", Path_NotProcessed),
      Level = 2)
  }

  ## Loading/merging EASIN data -----
  IASDT.R::CatTime("Loading/merging EASIN data", Level = 1)
  EASIN_Data_Orig <- purrr::map_dfr(.x = EASIN_Files, .f = IASDT.R::LoadAs)

  ## Save merged EASIN data - RData -----
  IASDT.R::CatTime("Save merged EASIN data - RData", Level = 1)
  save(
    EASIN_Data_Orig, file = IASDT.R::Path(Path_EASIN, "EASIN_Data_Orig.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Cleaning EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Cleaning EASIN data")
  EASIN_Data <- EASIN_Data_Orig %>%
    dplyr::mutate(Year = as.integer(Year), SpeciesName = NULL) %>%
    dplyr::rename(EASINID = SpeciesId) %>%
    # exclude observations with no spatial information or < StartYear
    dplyr::filter(!is.na(WKT), Year >= StartYear) %>%
    # Join with EASIN Taxa information
    dplyr::left_join(EASIN_Taxa, by = "EASINID")

  ## Extract coordinates from WKT string ----
  IASDT.R::CatTime("Extract coordinates from WKT string", Level = 1)
  WKTs <- dplyr::distinct(EASIN_Data, WKT) %>%
    dplyr::mutate(Coords = purrr::map(WKT, IASDT.R::Text2Coords)) %>%
    tidyr::unnest_wider(Coords) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c("Longitude", "Latitude"), .fns = ~ round(.x, 5)))

  ## Add coordinates to data and convert to sf ----
  IASDT.R::CatTime("Add coordinates to data and convert to sf", Level = 1)
  GridSf <- IASDT.R::LoadAs(GridSf) %>%
    magrittr::extract2("Grid_10_sf_s")
  LandGrids <- IASDT.R::LoadAs(LandGrids) %>%
    sf::st_drop_geometry() %>%
    dplyr::pull("CellCode")

  EASIN_Data <- dplyr::left_join(EASIN_Data, WKTs, by = "WKT") %>%
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
  IASDT.R::CatTime("Save cleaned EASIN Data - RData", Level = 1)
  save(EASIN_Data, file = IASDT.R::Path(Path_EASIN, "EASIN_Data.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Summarizing EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Summarizing EASIN data")

  GridR <- terra::unwrap(IASDT.R::LoadAs(GridR))

  ## NObs ----
  IASDT.R::CatTime("Number of observations per grid cell", Level = 1)

  IASDT.R::CatTime("All data", Level = 2)
  EASIN_NObs <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(NObs = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(GridR, field = "NObs") %>%
    terra::mask(GridR) %>%
    terra::wrap()

  Path_NObs <- IASDT.R::Path(Path_EASIN_Summary, "EASIN_NObs.RData")
  save(EASIN_NObs, file = Path_NObs)

  ### NObs per partner ----
  IASDT.R::CatTime("Number of observations per partner", Level = 2)

  EASIN_NObs_PerPartner <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(NObs = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::ReplaceSpace(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = GridR, field = "NObs") %>%
          terra::mask(GridR)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NObs_PerPartner <- IASDT.R::Path(
    Path_EASIN_Summary, "EASIN_NObs_PerPartner.RData")
  save(EASIN_NObs_PerPartner, file = Path_NObs_PerPartner)

  ## NSpecies ----
  IASDT.R::CatTime("Number of species per grid cell", Level = 1)

  IASDT.R::CatTime("All data", Level = 2)
  EASIN_NSp <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::distinct(CellCode, Species_name) %>%
    dplyr::count(CellCode) %>%
    dplyr::rename(NSp = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(GridR, field = "NSp") %>%
    terra::mask(GridR) %>%
    terra::wrap()

  Path_NSp <- IASDT.R::Path(Path_EASIN_Summary, "EASIN_NSp.RData")
  save(EASIN_NSp, file = Path_NSp)

  ### NSp per partner ----
  IASDT.R::CatTime("Number of species per per partner", Level = 2)

  EASIN_NSp_PerPartner <- sf::st_drop_geometry(EASIN_Data) %>%
    dplyr::distinct(CellCode, Species_name, DataPartnerName) %>%
    dplyr::count(CellCode, DataPartnerName) %>%
    dplyr::rename(NSp = n) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    dplyr::group_by(DataPartnerName) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::ReplaceSpace(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~ {
        terra::rasterize(x = .x, y = GridR, field = "NSp") %>%
          terra::mask(GridR)
      }
    ) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NSp_PerPartner <- IASDT.R::Path(
    Path_EASIN_Summary, "EASIN_NSp_PerPartner.RData")
  save(EASIN_NSp_PerPartner, file = Path_NSp_PerPartner)

  if (DeleteChunks) {
    fs::dir_delete(IASDT.R::Path(Path_EASIN_Interim, "FileParts"))
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  ## Species-specific data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Species-specific data")

  TimeStartData <- lubridate::now(tzone = "CET")

  EASIN_Data2 <- dplyr::group_by(EASIN_Data, Species_File) %>%
    dplyr::group_split() %>%
    purrr::set_names(
      purrr::map_chr(., ~ IASDT.R::ReplaceSpace(.x$Species_File[1]))
    )

  seq_len(length(EASIN_Data2)) %>%
    purrr::walk(
      .f = ~ {
        # # ||||||||||||||||||||||||||||||||||
        # Data as sf object
        # # ||||||||||||||||||||||||||||||||||
        IASDT.R::SaveAs(
          InObj = EASIN_Data2[[.x]],
          OutObj = names(EASIN_Data2)[.x],
          OutPath =
            IASDT.R::Path(
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

        IASDT.R::SaveAs(
          InObj = DT_Grid, OutObj = paste0(names(EASIN_Data2)[.x], "_Grid"),
          OutPath = IASDT.R::Path(
            Path_EASIN_Grid, paste0(names(EASIN_Data2)[.x], "_Grid.RData")))

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - raster
        # # ||||||||||||||||||||||||||||||||||
        DT_R <- terra::rasterize(DT_Grid, GridR, field = "NObs")
        DT_R$PA <- terra::as.int(DT_R > 0)
        IASDT.R::SaveAs(
          InObj = terra::wrap(DT_R),
          OutObj = paste0(names(EASIN_Data2)[.x], "_R"),
          OutPath = IASDT.R::Path(
            Path_EASIN_R, paste0(names(EASIN_Data2)[.x], "_R.RData")))

        return(invisible(NULL))
      }, .progress = FALSE)

  rm(EASIN_Data2, envir = environment())

  IASDT.R::CatDiff(
    TimeStartData,
    Prefix = "Preparing Species-specific data was finished in ", Level = 1)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  if (Plot) {
    IASDT.R::CatTime("Plotting")
    IASDT.R::EASIN_Plot(EnvFile = EnvFile)
  }

  # # ..................................................................... ###

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nProcessing EASIN data was finished in ")

  return(invisible(NULL))
}
