# # |------------------------------------------------------------------------| #
# EASIN_Process ----
## |------------------------------------------------------------------------| #

#' Extracts and processes EASIN data
#'
#' This function extracts and processes data from the European Alien Species
#' Information Network ([EASIN](https://easin.jrc.ec.europa.eu/)) for vascular
#' plants. This function extracts plant species data from the EASIN database,
#' matches them with a pre-processed standardized list of taxa, and prepares
#' species-specific maps and summary maps. It also supports downloading data in
#' chunks, handling pagination, and retrying failed downloads.
#' @param ExtractTaxa Logical. If `TRUE`, the function will extract the EASIN
#'   taxonomy list using [EASIN_Taxonomy]. Default is `TRUE`.
#' @param ExtractData Logical. If `TRUE`, the function will download EASIN
#'   species occurrence data using [EASIN_Down]. Default is `TRUE`.
#' @param NDownTries Integer. Number of attempts to retry downloading data in
#'   case of failure. Default is 10.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 6.
#' @param SleepTime Numeric. Time in seconds to wait between download attempts
#'   and between chunks. Default is 10 seconds.
#' @param NSearch Integer. Number of observations or species to download during
#'   EASIN taxonomy or data extraction, respectively. Default is 1000.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param EnvFile Character. The path to the environment file containing
#'   variables required by the function. Default is ".env".
#' @param DeleteChunks Logical. If `TRUE`, the function will delete intermediate
#'   files after processing. Default is `FALSE`.
#' @param StartYear Integer. Minimum year for filtering species occurrence data.
#'   Records before this year will be excluded. Default is `1981`, which matches
#'   the year ranges of CHELSA current climate data.
#' @param Plot Logical. If `TRUE`, the function will generate summary plots of
#'   the processed data using [EASIN_Plot]. Default is `TRUE`.
#' @author Ahmed El-Gabbas
#' @name EASIN_Process
#' @return The function Returns `NULL` invisibly after completing the data
#'   extraction, processing, and optional plotting. The function saves multiple
#'   outputs to disk, including the extracted and processed EASIN data,
#'   species-specific data files, and summary statistics. The main outputs are:
#'   - `EASIN_Taxa.RData`: A dataset containing the standardized EASIN taxonomy.
#'   - `EASIN_Data.RData`: A cleaned and merged dataset of species occurrence
#'   data.
#'   - `EASIN_NObs.RData`: A rasterized dataset showing the number of
#'   observations per grid cell.
#'   - `EASIN_NObs_PerPartner.RData`: A rasterized dataset showing the number of
#'   observations per data partner.
#'   - `EASIN_NSp.RData`: A rasterized dataset showing the number of species per
#'   grid cell.
#'   - `EASIN_NSp_PerPartner.RData`: A rasterized dataset showing the number of
#'   species per data partner.
#'   -  Species-specific data files, saved as both sf and raster objects.
#' @note
#'   - The function assumes that the necessary environment variables are
#' correctly set up in the specified `.env` file. Users should ensure that all
#' required files and directories are accessible before running the function.
#'   - The function skips processing (i.e. reuse) species data or data chunks if
#' the data already exist on the raw directory. The function assumes that the
#' contents of this folder should be removed as part of the data workflow.
#' Skipping processing available data can help not to re-download already
#' available data from the EASIN server.
#'   - This function depends on the following functions: [EASIN_Taxonomy] for
#'   getting the most recent EASIN taxonomy; [EASIN_Down] for processing EASIN
#'   dataset; and [EASIN_Plot] for plotting.
#' @export

EASIN_Process <- function(
    ExtractTaxa = TRUE, ExtractData = TRUE, NDownTries = 10,
    NCores = 6, SleepTime = 10, NSearch = 1000, FromHPC = TRUE,
    EnvFile = ".env", DeleteChunks = TRUE, StartYear = 1981, Plot = TRUE) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("ExtractTaxa", "ExtractData", "FromHPC", "Verbose",
             "Plot", "DeleteChunks"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("DownTries", "NCores", "SleepTime", "NSearch", "StartYear"))

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- TaxaInfoFile <- Name <- speciesKey <- EASINID <- taxon_name <-
    SpeciesId <- CellCode <- DataPartnerName <- Species_name <- Species_File <-
    EASIN_Ref <- Year <- WKT <- Path_EASIN <- Path_EASIN_Interim <-
    n <- x <- Path_Grid_Ref <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Environment variables")
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_EASIN", "DP_R_EASIN", FALSE, FALSE,
      "Path_EASIN_Interim", "DP_R_EASIN_Interim", FALSE, FALSE,
      "TaxaInfoFile", "DP_R_TaxaInfo_RData", FALSE, TRUE,
      "EASIN_Ref", "DP_R_TaxaInfo_EASIN", FALSE, TRUE,

      # The following are needed for other called functions
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Path_EASIN_Summary", "DP_R_EASIN_Summary", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_EASIN", "DP_R_EASIN_Local", FALSE, FALSE,
      "Path_EASIN_Interim", "DP_R_EASIN_Interim_Local", FALSE, FALSE,
      "TaxaInfoFile", "DP_R_TaxaInfo_RData_Local", FALSE, TRUE,
      "EASIN_Ref", "DP_R_TaxaInfo_EASIN_Local", FALSE, TRUE,

      # The following are needed for other called functions
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Path_EASIN_Summary", "DP_R_EASIN_Summary_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # General input data + Paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking input and output paths")

  Path_EASIN_DT <- file.path(Path_EASIN, "Sp_DT")
  Path_EASIN_Grid <- file.path(Path_EASIN, "Sp_Grid")
  Path_EASIN_R <- file.path(Path_EASIN, "Sp_R")
  Path_EASIN_Summary <- file.path(Path_EASIN, "Summary")

  fs::dir_create(
    c(Path_EASIN_DT, Path_EASIN_Grid, Path_EASIN_R, Path_EASIN_Summary,
      Path_EASIN_Interim))

  # # ||||||||||||||||||||||||||||||||||||||||||||||

  ## Grid - raster ----
  GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop(
      paste0("Path for the reference grid does not exist: ", GridR),
      call. = FALSE)
  }

  ## Grid - sf ----
  GridSf <- file.path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(GridSf)) {
    stop(
      paste0("Path for the reference grid does not exist: ", GridSf),
      call. = FALSE)
  }

  ## Grid - sf - study area ----
  # Grid ID overlapping with study area
  LandGrids <- file.path(Path_Grid, "Grid_10_Land_sf.RData")
  if (!file.exists(LandGrids)) {
    stop(
      paste0("Path for the reference grid does not exist: ", LandGrids),
      call. = FALSE)
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

    EASIN_Taxa_Orig <- IASDT.R::EASIN_Taxonomy(
      # Download EASIN taxa
      BaseURL = "https://easin.jrc.ec.europa.eu/apixg/catxg",
      Kingdom = "Plantae", Phylum = "Tracheophyta", NSearch = NSearch)
    save(
      EASIN_Taxa_Orig, file = file.path(Path_EASIN, "EASIN_Taxa_Orig.RData"))


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

    IASDT.R::CatTime("Check EASIN taxa not in the list",  Level = 1)
    # EASIN taxa not in the reference list
    New_EASIN_Taxa <- setdiff(EASIN_Taxa$EASIN_Name, EASIN_Ref$Name)
    if (length(New_EASIN_Taxa) > 0) {
      tibble::tibble(NewTaxa = New_EASIN_Taxa) %>%
        readr::write_tsv(
          file = file.path(Path_EASIN, "New_EASIN_Taxa.txt"),
          progress = FALSE)
    }
    ## Save EASIN taxa to disk ----
    IASDT.R::CatTime("Save EASIN taxa to disk",  Level = 1)
    save(EASIN_Taxa, file = file.path(Path_EASIN, "EASIN_Taxa.RData"))

  } else {

    IASDT.R::CatTime("Loading EASIN taxa list")
    load(file.path(Path_EASIN, "EASIN_Taxa.RData"))

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

    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

    future::plan(future::cluster, workers = NCores, gc = TRUE)
    on.exit(future::plan(future::sequential))

    # Start downloading, allow for a maximum of `NumDownTries` trials
    Try <- 0

    # Start downloading ----
    while (TRUE) {
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
          X = NotProcessed, FUN = IASDT.R::EASIN_Down,
          Path_Raw = Path_EASIN_Interim, DeleteChunks = DeleteChunks,
          NSearch = NSearch, SleepTime = SleepTime,
          future.scheduling = Inf, future.seed = TRUE,
          future.globals = c(
            "Path_EASIN_Interim", "NSearch", "DeleteChunks", "SleepTime"),
          future.packages = c("dplyr", "jsonlite")),
        silent = TRUE)

      if (inherits(Down, "try-error")) {
        next
      }

      rm(Down)
      Sys.sleep(SleepTime)
    }

    IASDT.R::CatDiff(
      TimeStartData,
      Prefix = "Downloading EASIN data was finished in ", Level = 1)

    # Stop cluster ----
    IASDT.R::CatTime("Stop cluster", Level = 1)
    future::plan(future::sequential, gc = TRUE)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Merging EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Merging EASIN data")

  ## Loading input maps -----
  IASDT.R::CatTime("Loading input maps", Level = 1)

  EASIN_Files <- list.files(
    Path_EASIN_Interim, full.names = TRUE, pattern = ".RData")

  NotProcessed <- setdiff(
    paste0(EASIN_Taxa$EASINID, ".RData"), basename(EASIN_Files))

  if (length(NotProcessed) > 0) {
    PathNotProcessed <- file.path(Path_EASIN, "NotProcessedID.txt")
    stringr::str_remove_all(NotProcessed, ".RData") %>%
      cat(file = PathNotProcessed, sep = "\n")
    IASDT.R::CatTime(
      paste0(
        "There are ", length(NotProcessed), " not processed EASIN ID(s)."),
      Level = 2)

    IASDT.R::CatTime(
      paste0("EASIN IDs are saved to: ", PathNotProcessed), Level = 2)
  }

  ## Merging EASIN data -----
  IASDT.R::CatTime("Merging EASIN data", Level = 1)

  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 2)

  withr::local_options(future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  future::plan(future::cluster, workers = NCores, gc = TRUE)
  on.exit(future::plan(future::sequential), add = TRUE)

  EASIN_Data_Orig <- future.apply::future_lapply(
    X = EASIN_Files, FUN = IASDT.R::LoadAs,
    future.scheduling = Inf, future.seed = TRUE,
    future.globals = "EASIN_Files") %>%
    dplyr::bind_rows()

  ## Save merged EASIN data -----
  IASDT.R::CatTime("Save merged EASIN data", Level = 1)
  save(
    EASIN_Data_Orig, file = file.path(Path_EASIN, "EASIN_Data_Orig.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Cleaning EASIN data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Cleaning EASIN data")

  GridSf <- IASDT.R::LoadAs(GridSf) %>%
    magrittr::extract2("Grid_10_sf_s")
  LandGrids <- IASDT.R::LoadAs(LandGrids) %>%
    sf::st_drop_geometry() %>%
    dplyr::pull("CellCode")

  EASIN_Data <- EASIN_Data_Orig %>%
    dplyr::mutate(Year = as.integer(Year), SpeciesName = NULL) %>%
    dplyr::rename(EASINID = SpeciesId) %>%
    # exclude observations with no spatial information or < StartYear
    dplyr::filter(!is.na(WKT), Year >= StartYear) %>%
    # Join with EASIN Taxa information
    dplyr::left_join(EASIN_Taxa, by = "EASINID")

  ## Extract coordinates from WKT ----
  IASDT.R::CatTime("Extract coordinates from WKT", Level = 1)
  WKTs <- EASIN_Data %>%
    dplyr::distinct(WKT) %>%
    dplyr::pull(WKT) %>%
    future.apply::future_lapply(
      X = .,
      FUN = function(x) {
        tibble::tibble(WKT = x, IASDT.R::Text2Coords(x))
      },
      future.scheduling = Inf, future.seed = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(dplyr::across(
      .cols = c("Longitude", "Latitude"), .fns = ~round(.x, 5)))

  ## Add coordinates to data and convert to sf ----
  IASDT.R::CatTime("Add coordinates to data and convert to sf", Level = 1)
  EASIN_Data <- EASIN_Data %>%
    dplyr::left_join(WKTs, by = "WKT") %>%
    # convert to sf object, while keeping original coordinates as columns
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    # spatial transform the data to EPSG:3035
    sf::st_transform(3035) %>%
    # Add grid ID
    sf::st_join(GridSf) %>%
    # exclude observations not overlapping with the study area
    dplyr::filter(CellCode %in% LandGrids)


  future::plan(future::sequential)
  rm(LandGrids, WKTs)

  ## Save cleaned EASIN Data ----
  IASDT.R::CatTime("Save cleaned EASIN Data", Level = 1)
  save(EASIN_Data, file = file.path(Path_EASIN, "EASIN_Data.RData"))

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

  Path_NObs <- file.path(Path_EASIN_Summary, "EASIN_NObs.RData")
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
      purrr::map_chr(., ~IASDT.R::ReplaceSpace(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~{
        terra::rasterize(x = .x, y = GridR, field = "NObs") %>%
          terra::mask(GridR)
      }) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NObs_PerPartner <- file.path(
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

  Path_NSp <- file.path(Path_EASIN_Summary, "EASIN_NSp.RData")
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
      purrr::map_chr(., ~IASDT.R::ReplaceSpace(.x$DataPartnerName[1]))) %>%
    purrr::map(
      .f = ~{
        terra::rasterize(x = .x, y = GridR, field = "NSp") %>%
          terra::mask(GridR)
      }) %>%
    terra::rast() %>%
    terra::wrap()

  Path_NSp_PerPartner <- file.path(
    Path_EASIN_Summary, "EASIN_NSp_PerPartner.RData")
  save(EASIN_NSp_PerPartner, file = Path_NSp_PerPartner)

  if (DeleteChunks) {
    fs::dir_delete(file.path(Path_EASIN_Interim, "FileParts"))
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
      purrr::map_chr(., ~IASDT.R::ReplaceSpace(.x$Species_File[1])))

  seq_len(length(EASIN_Data2)) %>%
    purrr::walk(
      .f = ~{

        # # ||||||||||||||||||||||||||||||||||
        # Data as sf object
        # # ||||||||||||||||||||||||||||||||||
        IASDT.R::SaveAs(
          InObj = EASIN_Data2[[.x]],
          OutObj = names(EASIN_Data2)[.x],
          OutPath =
            file.path(
              Path_EASIN_DT, paste0(names(EASIN_Data2)[.x], "_DT.RData")))

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
          OutPath = file.path(
            Path_EASIN_Grid, paste0(names(EASIN_Data2)[.x], "_Grid.RData")))

        # # ||||||||||||||||||||||||||||||||||
        # Presence grid - raster
        # # ||||||||||||||||||||||||||||||||||
        DT_R <- terra::rasterize(DT_Grid, GridR, field = "NObs")
        DT_R$PA <- terra::as.int(DT_R > 0)
        IASDT.R::SaveAs(
          InObj = terra::wrap(DT_R),
          OutObj = paste0(names(EASIN_Data2)[.x], "_R"),
          OutPath = file.path(
            Path_EASIN_R, paste0(names(EASIN_Data2)[.x], "_R.RData")))

        return(invisible(NULL))
      }, .progress = FALSE)

  rm(EASIN_Data2)

  IASDT.R::CatDiff(
    TimeStartData,
    Prefix = "Preparing Species-specific data was finished in ", Level = 1)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  if (Plot) {
    IASDT.R::CatTime("Plotting")
    IASDT.R::EASIN_Plot(EnvFile = EnvFile, FromHPC = FromHPC)
  }

  # # ..................................................................... ###

  ## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Processing EASIN data was finished in ")

  return(invisible(NULL))
}
