## |------------------------------------------------------------------------| #
# Sampling_Efforts ----
## |------------------------------------------------------------------------| #

#' Download and Process sampling efforts data
#'
#' This script download GBIF data for all vascular plants (grouped by order) in
#' the study area. Data is converted to raster to represent the number of
#' vascular plant observation per grid
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param EnvFile Character. The path to the environment file containing
#'   variables required by the function. Default: `.env`.
#' @param Renviron Character. The path to the `.Renviron` file containing GBIF
#'   login credentials (email, user, password). Default: `.Renviron`.
#' @param RequestData Logical. If `TRUE`, the function requests data from GBIF.
#'   If `FALSE`, previously requested data is loaded from disk. Defaults to
#'   `TRUE`.
#' @param DownloadData Logical. If `TRUE`, the function downloads data and
#'   stores it on disk. If `FALSE`, it skips the download step. Defaults to
#'   `TRUE`.
#' @param Boundaries Numeric vector of length 4. Specifies geographical
#'   boundaries for the requested GBIF data in the order: Left, Right, Bottom,
#'   Top. Default: `c(-30, 50, 25, 75)`.
#' @param NCores Integer. The number of cores to use for parallel processing.
#' @param StartYear Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be requested from GBIF. Default: `1980`.
#' @note The function is expected to take substantial amount of time (> 9 hours
#'   on a Windows PC running with 6 cores). The data request is expected to take
#'   5 hours, as it includes waiting for the GBIF data to be ready.
#' @author Ahmed El-Gabbas
#' @name Sampling_Efforts
#' @export

Sampling_Efforts <- function(
    FromHPC = TRUE, EnvFile = ".env", Renviron = ".Renviron",
    RequestData = TRUE, DownloadData = TRUE, NCores = 6, StartYear = 1980,
    Boundaries = c(-30, 50, 25, 75)) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Renviron", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("FromHPC", "RequestData", "DownloadData"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "Boundaries", "StartYear"))

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Efforts <- Path_Efforts_Raw <- Path_Efforts_Interim <-
    Path_Grid <- Taxa_Stand <- EU_Bound <- orderKey <- Request <-
    DownDetails <- Size <- DownPath <- PathCSV <- order <- class <-
    TotalRecords <- CellCode <- speciesKey <- coordinateUncertaintyInMeters <-
    decimalLongitude <- decimalLatitude <- year <- UncertainKm <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime(
    "Ensure that GBIF access information is available or can be read")

  IASDT.R::Check_GBIF(Renviron = Renviron)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts", FALSE, FALSE,
      "Path_Efforts_Raw", "DP_R_Efforts_Raw", FALSE, FALSE,
      "Path_Efforts_Interim", "DP_R_Efforts_Interim", FALSE, FALSE,
      "Taxa_Stand", "DP_R_TaxaStand", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts_Local", FALSE, FALSE,
      "Path_Efforts_Raw", "DP_R_Efforts_Raw_Local", FALSE, FALSE,
      "Path_Efforts_Interim", "DP_R_Efforts_Interim_Local", FALSE, FALSE,
      "Taxa_Stand", "DP_R_TaxaStand_Local", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading input data ------
  IASDT.R::CatTime("Loading input data")

  ## Create paths -----
  Path_Efforts_Requests <- file.path(Path_Efforts, "Requests")
  Path_Efforts_Data <- file.path(Path_Efforts_Interim, "CleanedData")
  fs::dir_create(c(
    Path_Efforts, Path_Efforts_Raw, Path_Efforts_Interim, Path_Efforts_Data,
    Path_Efforts_Requests))

  ## Reference grid ----
  Grids <- Path_Grid %>%
    file.path(c("Grid_10_Land_Crop_sf.RData", "Grid_10_Land_Crop.RData"))

  if (!all(file.exists(Grids))) {
    stop(
      paste0("The following grid file(s) do not exist:\n",
             paste0(" >>> ", Grids[!file.exists(Grids)], recycle0 = "\n")),
      call. = FALSE)
  }

  ## IAS species list ----
  IAS_List <- unique(dplyr::pull(readRDS(Taxa_Stand), "speciesKey"))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Request efforts data ------
  IASDT.R::CatTime("Request efforts data")

  if (RequestData) {

    IASDT.R::CatTime("Requesting efforts data")

    Efforts_AllRequests <- IASDT.R::Efforts_Request(
      NCores = NCores, Path_Requests = Path_Efforts_Requests,
      Path_Efforts = Path_Efforts, StartYear = StartYear,
      Boundaries = Boundaries)

  } else {

    IASDT.R::CatTime("Efforts data was not requested, but loaded", Level = 1)
    Efforts_AllRequests <- IASDT.R::LoadAs(
      file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Download efforts data ------
  IASDT.R::CatTime("Download efforts data")

  if (DownloadData) {

    Efforts_AllRequests <- IASDT.R::Efforts_Download(
      NCores = NCores, Path_Raw = Path_Efforts_Raw,
      Path_Interim = Path_Efforts_Interim, Path_Efforts = Path_Efforts)

  } else {

    IASDT.R::CatTime("Efforts data was not downloaded", Level = 1)
    Efforts_AllRequests <- IASDT.R::LoadAs(
      file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Processing efforts data ------

  IASDT.R::CatTime("Processing efforts data")
  IASDT.R::Efforts_Process(
    NCores = NCores, Path_Efforts = Path_Efforts,
    Path_Interim = Path_Efforts_Interim, Path_Data = Path_Efforts_Data,
    Path_Grid = Path_Grid, IAS_List = IAS_List,
    Efforts_AllRequests = Efforts_AllRequests)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # # Plotting ----

  IASDT.R::CatTime("\nPlotting efforts data")
  IASDT.R::Efforts_Plot(Path_Efforts = Path_Efforts, EU_Bound = EU_Bound)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, CatInfo = FALSE,
    Prefix = "\nProcessing efforts data took ", ... = "\n")

  return(invisible(NULL))
}

## ............................................ ------
## ............................................ ------


## |------------------------------------------------------------------------| #
# Efforts_Request ----
## |------------------------------------------------------------------------| #

#' Request and Manage GBIF Data for Vascular Plant Orders
#'
#' This function requests GBIF data for each vascular plant order, processes it
#' in parallel, and manages the data download and storage. If data is already
#' available, it loads the data instead of making new requests.
#'
#' @param NCores Integer. The number of cores to use for parallel processing (1 to 3). Must be a positive integer, up to a maximum of 3.
#' @param Path_Requests Character. The directory path to save individual request
#'   files.
#' @param Path_Efforts Character. The directory path to save the final compiled
#'   data.
#' @param StartYear Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be requested from GBIF. Defaults to
#'   1980.
#' @param Boundaries Numeric vector of length 4. Specifies geographical
#'   boundaries for the requested GBIF data in the order: Left, Right, Bottom,
#'   Top. Defaults to c(-30, 50, 25, 75).
#' @return The function returns the GBIF data requests processed and stored in
#'   the specified directories.
#' @author Ahmed El-Gabbas
#' @name Efforts_Request
#' @export


Efforts_Request <- function(
    NCores, Path_Requests, Path_Efforts, StartYear, Boundaries) {

  # In earlier tries, requesting all vascular plants occurrences in a single
  # request returned 80 GB compressed file. The extracted "occurrences.txt" is
  # >280 GB (220M observations).
  #
  # The following makes individual request for each vascular plant order. This
  # can take up to 5 hours for the data to be ready

  .StartTimeRequest <- lubridate::now(tzone = "CET")

  if (missing(NCores) || !is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer between 1 and 3.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- DownDetails <- orderKey <- Size <- NumberDatasets <-
    TotalRecords <- Created <- Modified <- EraseAfter <- NULL

  # # ..................................................................... ###

  # Prepare working on parallel -----
  #
  # GBIF allows only 3 parallel requests. Here I wait until previous request
  # is finished.
  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- snow::makeSOCKcluster(min(NCores, 3))
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  snow::clusterEvalQ(
    cl = c1,
    expr = IASDT.R::LoadPackages(List = c("dplyr", "IASDT.R", "rgbif")))

  # # ..................................................................... ###

  # Requesting efforts Data on parallel -----
  IASDT.R::CatTime("Requesting efforts Data on parallel", Level = 1)
  IASDT.R::CatTime("This may take up to 4 hours", Level = 2)

  # Extract taxonomic info for vascular plants orders
  SelectedCols <- c(
    "class", "classKey", "order", "orderKey", "numDescendants")

  Efforts_AllRequests <- rgbif::name_backbone("Tracheophyta") %>%
    dplyr::pull("phylumKey") %>%
    rgbif::name_lookup(rank = "ORDER", higherTaxonKey = .) %>%
    # Get info on order names
    magrittr::extract2("data") %>%
    dplyr::select(tidyselect::all_of(SelectedCols)) %>%
    dplyr::mutate(
      Request = furrr::future_map(
        .x = orderKey,
        .f = ~{

          Request_ID <- paste0("Request_", .x)
          Request_Path <- file.path(Path_Requests, paste0(Request_ID, ".RData"))

          if (file.exists(Request_Path)) {
            # load previous request
            Down <- IASDT.R::LoadAs(Request_Path)
          } else {
            # Attempt the request with error handling
            tryCatch({
              # Make data request
              Down <- rgbif::occ_download(
                rgbif::pred_in("taxonKey", .x),
                # Only with coordinates & no spatial issues
                rgbif::pred("hasCoordinate", TRUE),
                rgbif::pred("hasGeospatialIssue", FALSE),
                # Only after (>=) a certain year
                rgbif::pred_gte("year", StartYear),
                # Only within specific boundaries
                rgbif::pred_within(
                  value = IASDT.R::DownBoundary(
                    Left = Boundaries[1], Right = Boundaries[2],
                    Bottom = Boundaries[3], Top = Boundaries[4])),
                format = "SIMPLE_CSV")

              IASDT.R::SaveAs(
                InObj = Down, OutObj = Request_ID, OutPath = Request_Path)

            },
            error = function(e) {
              stop(paste0(
                "Failed to request data for taxonKey ", .x, ": ",
                conditionMessage(e)),
                call. = FALSE)
            })
          }

          # Waiting for data to be ready
          rgbif::occ_download_wait(Down, quiet = TRUE)

          return(Down)
        },
        .options = furrr::furrr_options(seed = TRUE))) %>%
    dplyr::rowwise() %>%

    # Add columns for metadata
    dplyr::mutate(
      DownDetails = list(rgbif::occ_download_wait(Request, quiet = TRUE)),
      # Extract some info from metadata
      DownloadKey = DownDetails$key,
      DOI = DownDetails$doi,
      Created = DownDetails$created,
      Modified = DownDetails$modified,
      EraseAfter = DownDetails$eraseAfter,
      DownLink = DownDetails$downloadLink,
      Size = DownDetails$size,
      TotalRecords = DownDetails$totalRecords,
      NumberDatasets = DownDetails$numberDatasets,
      Status = DownDetails$status,
      # Size of data in megabytes
      Size = as.numeric(Size) / (1024 * 1024),
      # Convert some columns to numeric (double)
      dplyr::across(Size:NumberDatasets, as.numeric),
      # Convert some columns to integer
      dplyr::across(TotalRecords:NumberDatasets, as.integer),
      # Convert some columns to date type
      dplyr::across(c(Created, Modified, EraseAfter), lubridate::as_date)) %>%
    dplyr::ungroup()  %>%
    # how to cite data
    dplyr::mutate(Citation = purrr::map_chr(Request, attr, "citation"))

  # # ..................................................................... ###

  # Save efforts request data ------
  IASDT.R::CatTime("Save efforts request data", Level = 1)

  save(Efforts_AllRequests,
       file = file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeRequest, CatInfo = FALSE,
    Prefix = "Requesting efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(Efforts_AllRequests)

}

## ............................................ ------
## ............................................ ------

## |------------------------------------------------------------------------| #
# Efforts_Download ----
## |------------------------------------------------------------------------| #

#' Download and Manage GBIF Data for Vascular Plant Orders
#'
#' This function handles the downloading of GBIF data in parallel, checks the
#' validity of downloaded files, and stores the data in specified directories.
#' If data has already been downloaded, it validates the files instead of
#' downloading them again.
#' @param NCores Integer. Number of cores to use for parallel processing.
#' @param Path_Raw Character. Path where the raw downloaded data will be saved.
#' @param Path_Interim Character. Path where the interim CSV files will be
#'   saved.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved.
#' @name Efforts_Download
#' @author Ahmed El-Gabbas
#' @return A data frame (`Efforts_AllRequests`) with updated download paths and
#'   interim file paths.
#' @export

Efforts_Download <- function(NCores, Path_Raw, Path_Interim, Path_Efforts) {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  DownPath <- Request <- NULL

  # # ..................................................................... ###

  Path_Efforts_Request <- file.path(Path_Efforts, "Efforts_AllRequests.RData")
  if (!file.exists(Path_Efforts_Request)) {
    stop(
      paste0(
        "The path for the `Efforts_AllRequests` data does not exist: ",
        Path_Efforts_Request),
      call. = FALSE)
  }

  Efforts_AllRequests <- IASDT.R::LoadAs(Path_Efforts_Request)

  # # ..................................................................... ###

  ## Prepare working on parallel -----
  #
  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  snow::clusterEvalQ(
    cl = c1,
    expr = IASDT.R::LoadPackages(
      List = c("dplyr", "IASDT.R", "rgbif", "stringr")))

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  IASDT.R::CatTime("Downloading/checking efforts data", Level = 1)

  Efforts_AllRequests <- Efforts_AllRequests %>%
    dplyr::mutate(
      # Download datasets on parallel
      DownPath = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          DownFile <- file.path(Path_Raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not exist download it
          if (file.exists(DownFile)) {
            FileOkay <- system2(
              "unzip", args = c("-t", DownFile),
              stdout = TRUE, stderr = TRUE) %>%
              stringr::str_detect("No errors detected in compressed data") %>%
              any()

            if (FileOkay) {
              Success <- TRUE
            } else {
              Success <- FALSE
            }

          } else {
            Success <- FALSE
          }

          # Try downloading data for a max of 3 attempts, each with 20 mins
          # time out
          withr::local_options(list(timeout = 1200))

          Attempt <- 1
          Attempts <- 3

          while (!Success && (Attempt <= Attempts)) {
            tryCatch({
              suppressMessages(
                rgbif::occ_download_get(
                  key = .x, path = Path_Raw, overwrite = TRUE))

              ZipStatus <- system2(
                "unzip", args = c("-t", DownFile),
                stdout = TRUE, stderr = TRUE) %>%
                stringr::str_detect("No errors detected in compressed data") %>%
                any()

              # Ensure Success is only TRUE if both the zip file exists and
              # passes integrity check
              Success <- file.exists(DownFile) && ZipStatus

            },
            error = function(e) {
              if (Attempt < Attempts) {
                Attempt <- Attempt + 1
              } else {
                stop(
                  paste0(
                    "Failed to download data after ", Attempts, " attempts: ",
                    conditionMessage(e)),
                  call. = FALSE)
              }
            })
          }

          return(DownFile)

        }, .options = furrr::furrr_options(seed = TRUE)),

      PathCSV = purrr::map_chr(
        .x = DownPath,
        .f = ~{
          stringr::str_replace(basename(.x), "zip", "csv") %>%
            file.path(Path_Interim, .)
        }))

  save(Efforts_AllRequests,
       file = file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeDown, CatInfo = FALSE,
    Prefix = "Downloading efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(Efforts_AllRequests)
}


## ............................................ ------
## ............................................ ------

## |------------------------------------------------------------------------| #
# Efforts_Process ----
## |------------------------------------------------------------------------| #

#' Process and Summarize GBIF Data for vascular plants
#'
#' This function processes GBIF data for vascular plants by extracting,
#' summarizing, and saving the data. It also creates summary maps and saves them
#' in both `RData` and `TIFF` formats.
#'
#' @param NCores Numeric. The number of cores to use for parallel processing.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved.
#' @param Path_Interim Character. The directory path to save interim data.
#' @param Path_Data Character. The directory path to save detailed effort data
#'   as `RData`.
#' @param Path_Grid Character. The directory path to load the grid data.
#' @param IAS_List A list of invasive alien species keys.
#' @param Efforts_AllRequests A data frame containing the details of the GBIF
#'   download, including paths to CSV files, zip files, order, class, and total
#'   records.
#' @return The function returns and saves the GBIF data summary.
#' @author Ahmed El-Gabbas
#' @name Efforts_Process
#' @export

Efforts_Process <- function(
    NCores, Path_Efforts, Path_Interim, Path_Data, Path_Grid,
    IAS_List, Efforts_AllRequests) {

  .StartTimeProcess <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  if (missing(NCores) || !is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  if (missing(Path_Data) || !is.character(Path_Data) ||
      !dir.exists(Path_Data)) {
    stop("Path_Data must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Grid) || !is.character(Path_Grid) ||
      !dir.exists(Path_Grid)) {
    stop("Path_Grid must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Interim) || !is.character(Path_Interim) || !dir.exists(Path_Interim)) {
    stop("Path_Interim must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Efforts) || !is.character(Path_Efforts) ||
      !dir.exists(Path_Efforts)) {
    stop("Path_Efforts must be a valid directory path.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  speciesKey <- CellCode <- coordinateUncertaintyInMeters <- ObsN <-
    decimalLongitude <- decimalLatitude <- year <- UncertainKm <-
    Latitude <- Longitude <- taxonRank <- NULL

  # # ..................................................................... ###

  Path_Grid_R <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  Path_Grid_SF <- file.path(Path_Grid, "Grid_10_Land_Crop_sf.RData")

  if (!file.exists(Path_Grid_R)) {
    stop(
      paste0(
        "Reference grid was not found at the specified path: ", Path_Grid_R),
      call. = FALSE)
  }

  if (!file.exists(Path_Grid_SF)) {
    stop(
      paste0("The path for the reference grid does not exist: ", Path_Grid_SF),
      call. = FALSE)
  }

  Grid_SF <- IASDT.R::LoadAs(Path_Grid_SF)

  # # ..................................................................... ###

  # Prepare working on parallel -----
  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  invisible(
    snow::clusterEvalQ(
      cl = c1,
      IASDT.R::LoadPackages(
        List = c(
          "terra", "IASDT.R", "stringr", "fs", "sf", "readr", "dplyr"))))

  snow::clusterExport(
    cl = c1,
    list = c(
      "Path_Efforts", "Path_Interim", "Efforts_AllRequests", "Path_Grid_R",
      "Path_Data", "Grid_SF", "IAS_List"),
    envir = environment())

  # # ..................................................................... ###

  # Processing efforts data ------
  IASDT.R::CatTime("Processing efforts data", Level = 1)

  Efforts_Summary0 <- future.apply::future_lapply(
    X = seq_len(nrow(Efforts_AllRequests)),
    FUN = function(x) {

      CSV_File <- Efforts_AllRequests$PathCSV[x]
      Zip_File <- Efforts_AllRequests$DownPath[x]
      order <- Efforts_AllRequests$order[x]
      class <- Efforts_AllRequests$class[x]
      TotalRecords <- Efforts_AllRequests$TotalRecords[x]
      ClassOrder <- paste0(class, "_", order)

      Summary_Path <- file.path(
        Path_Interim, paste0("Summary_", ClassOrder, ".RData"))
      Path_DT <- file.path(Path_Data, paste0(ClassOrder, ".RData"))
      Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

      # Exclude processed Order
      if (file.exists(Summary_Path)) {
        return(Summary_Path)
      }

      # Exclude Orders without any observations
      if (TotalRecords == 0) {
        tibble::tibble(
          ObsN = 0, ObsN_Native = 0,
          NObs_R = list(NA_integer_), NObs_Native_R = list(NA_integer_),
          NSp_R = list(NA_integer_), NSp_Native_R = list(NA_integer_),
          Path_DT = NA_character_) %>%
          IASDT.R::SaveAs(OutObj = ClassOrder, Summary_Path)
        return(Summary_Path)
      }

      # # ................................... ###

      # Extract observations file
      if (!file.exists(CSV_File)) {
        paste0("unzip -qqo ", Zip_File, " -d ", Path_Interim) %>%
          system() %>%
          invisible()
      }

      # # ................................... ###

      Success <- FALSE
      Attempt <- 1
      Attempts <- 5

      while (!Success && (Attempt <= Attempts)) {
        tryCatch({

          DT <- readr::read_tsv(
            file = CSV_File,
            col_select = tidyselect::all_of(
              c("coordinateUncertaintyInMeters", "decimalLongitude",
                "decimalLatitude", "year", "speciesKey", "taxonRank")),
            progress = FALSE, show_col_types = FALSE,
            col_types = readr::cols(
              coordinateUncertaintyInMeters = readr::col_double(),
              decimalLongitude = readr::col_double(),
              decimalLatitude = readr::col_double(),
              year = readr::col_integer(),
              speciesKey = readr::col_integer(),
              taxonRank = readr::col_character(),
              .default = readr::col_double()))

          DT <- dplyr::rename(
            DT,
            UncertainKm = coordinateUncertaintyInMeters,
            Longitude = decimalLongitude, Latitude = decimalLatitude) %>%
            dplyr::mutate(
              year = as.integer(year), UncertainKm = UncertainKm / 1000) %>%
            dplyr::filter(
              !is.na(Latitude), !is.na(Longitude),
              speciesKey != "",  year > 1980,
              taxonRank %in% c("FORM", "SPECIES", "SUBSPECIES", "VARIETY"),
              UncertainKm <= 100 | is.na(UncertainKm)) %>%
            sf::st_as_sf(
              coords = c("Longitude", "Latitude"),
              crs = 4326, remove  = FALSE) %>%
            # project to 3035
            sf::st_transform(3035) %>%
            sf::st_join(Grid_SF) %>%
            dplyr::filter(magrittr::not(is.na(CellCode)))

          ObsN <- nrow(DT)

          # # ................................... ###

          if (ObsN > 0) {

            # Save order data as RData
            IASDT.R::SaveAs(InObj = DT, OutObj = ClassOrder, OutPath = Path_DT)

            # Number of observations
            NObs_R <- sf::st_drop_geometry(DT) %>%
              dplyr::count(CellCode, name = "NObs") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(y = Grid_R, field = "NObs") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NObs_", ClassOrder)) %>%
              terra::wrap()

            # # ................................... ###

            # Number of species
            NSp_R <- sf::st_drop_geometry(DT) %>%
              dplyr::distinct(CellCode, speciesKey) %>%
              dplyr::count(CellCode, name = "NSp") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(y = Grid_R, field = "NSp") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NSp_", ClassOrder)) %>%
              terra::wrap()

            # # ................................... ###

            # Data on native species
            DT_Native <- dplyr::filter(DT, !(speciesKey %in% IAS_List))

            rm(DT)

            if (nrow(DT_Native) > 0) {

              ObsN_Native <- nrow(DT_Native)

              # # ................................... ###

              # Number of observations of native species
              NObs_Native_R <- sf::st_drop_geometry(DT_Native) %>%
                dplyr::count(CellCode, name = "NObs_Native") %>%
                dplyr::left_join(Grid_SF, by = "CellCode") %>%
                sf::st_as_sf() %>%
                terra::rasterize(Grid_R, field = "NObs_Native") %>%
                terra::classify(cbind(NA, 0)) %>%
                terra::mask(Grid_R) %>%
                IASDT.R::setRastCRS() %>%
                IASDT.R::setRastVals() %>%
                stats::setNames(paste0("NObsNative_", ClassOrder)) %>%
                terra::wrap()

              # # ................................... ###

              # Number of native species
              NSp_Native_R <- sf::st_drop_geometry(DT_Native) %>%
                dplyr::distinct(CellCode, speciesKey) %>%
                dplyr::count(CellCode, name = "NSp_Native") %>%
                dplyr::left_join(Grid_SF, by = "CellCode") %>%
                sf::st_as_sf() %>%
                terra::rasterize(Grid_R, field = "NSp_Native") %>%
                terra::classify(cbind(NA, 0)) %>%
                terra::mask(Grid_R) %>%
                IASDT.R::setRastCRS() %>%
                IASDT.R::setRastVals() %>%
                stats::setNames(paste0("NSpNative_", ClassOrder)) %>%
                terra::wrap()

            } else {
              ObsN_Native <- 0
              NObs_Native_R <- NSp_Native_R <- list()
            }

          } else {
            ObsN <- ObsN_Native <- 0
            Path_DT <- NA_character_
            NObs_R <- NSp_R <- NObs_Native_R <- NSp_Native_R <- NA_real_
          }

          tibble::tibble(
            ObsN = ObsN, ObsN_Native = ObsN_Native,
            NObs_R = list(NObs_R), NObs_Native_R = list(NObs_Native_R),
            NSp_R = list(NSp_R), NSp_Native_R = list(NSp_Native_R),
            Path_DT = Path_DT) %>%
            IASDT.R::SaveAs(OutObj = ClassOrder, OutPath = Summary_Path)

          Success <- TRUE

        },
        error = function(e) {
          IASDT.R::CatTime(
            paste0("Error on attempt #", Attempt, ": ", conditionMessage(e)),
            Level = 2)

          if (Attempt < Attempts) {
            Attempt <- Attempt + 1
          } else {
            stop(
              paste0(
                "Failed to Process efforts data after ", Attempts,
                " attempts: ", conditionMessage(e)),
              call. = FALSE)
          }
        })
      }

      fs::file_delete(CSV_File)

      return(Summary_Path)
    },
    future.scheduling = Inf, future.seed = TRUE)

  Efforts_Summary <- purrr::map_dfr(Efforts_Summary0, IASDT.R::LoadAs) %>%
    dplyr::bind_cols(Efforts_AllRequests, .)

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  # # ..................................................................... ###

  # Save Efforts_Summary ----
  IASDT.R::CatTime("Save `Efforts_Summary`", Level = 1)
  save(Efforts_Summary,
       file = file.path(Path_Efforts, "Efforts_Summary.RData"))

  # # ..................................................................... ###

  # Prepare summary maps ----
  IASDT.R::CatTime("Prepare summary maps", Level = 1)
  CalcNObsNSp <- function(List, Name) {
    purrr::map(.x = unlist(List), .f = terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE) %>%
      IASDT.R::setRastCRS() %>%
      IASDT.R::setRastVals() %>%
      stats::setNames(Name)
  }


  # Exclude orders with no data
  Efforts_SummaryR <- dplyr::filter(Efforts_Summary, ObsN > 0)

  Efforts_SummaryR <- list(
    CalcNObsNSp(Efforts_SummaryR$NObs_R, "NObs"),
    CalcNObsNSp(Efforts_SummaryR$NObs_Native_R, "NObs_Native"),
    CalcNObsNSp(Efforts_SummaryR$NSp_R, "NSp"),
    CalcNObsNSp(Efforts_SummaryR$NSp_Native_R, "NSp_Native")) %>%
    terra::rast() %>%
    IASDT.R::setRastCRS() %>%
    IASDT.R::setRastVals()

  # # ..................................................................... ###

  # Save summary data as RData ----
  IASDT.R::CatTime("Save as RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Efforts_SummaryR), OutObj = "Efforts_SummaryR",
    OutPath = file.path(Path_Efforts, "Efforts_SummaryR.RData"))

  # # ..................................................................... ###

  # Save summary data as tif ----
  IASDT.R::CatTime("Save as tif", Level = 1)
  terra::writeRaster(
    Efforts_SummaryR, overwrite = TRUE,
    filename = file.path(
      Path_Efforts,
      paste0("Efforts_GBIF_", names(Efforts_SummaryR), ".tif")))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeProcess, CatInfo = FALSE,
    Prefix = "Processing Efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}


## ............................................ ------
## ............................................ ------

## |------------------------------------------------------------------------| #
# Efforts_Plot ----
## |------------------------------------------------------------------------| #

#' Plot the output of efforts maps
#'
#' This function generates and saves multiple plots of plant observation
#' efforts, both in raw and log10 scales, using provided spatial boundary and
#' summary data.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved.
#' @param EU_Bound Character. Path to file containing country boundaries.
#' @return The function saves generated plots as JPEG files in the specified
#'   directory and returns NULL invisibly.
#' @author Ahmed El-Gabbas
#' @name Efforts_Plot
#' @details This function generates and saves effort maps visualizing the number
#'   of plant observations and species, including both native and non-native
#'   species, within Europe. It produces both standard and log10-scaled plots.
#' @export

Efforts_Plot <- function(Path_Efforts, EU_Bound) {

  File_SummaryR <- file.path(Path_Efforts, "Efforts_SummaryR.RData")
  if (!file.exists(File_SummaryR)) {
    stop(paste0("Summary maps cannot be loaded: ", File_SummaryR),
         call. = FALSE)
  }

  Efforts_SummaryR <- terra::unwrap(IASDT.R::LoadAs(File_SummaryR))

  # # ..................................................................... ###

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
      plot.title = ggplot2::element_text(
        size = 14, color = "blue", face = "bold", hjust = 0.5,
        margin = ggplot2::margin(0, 0, 0, 0)),
      plot.subtitle = ggplot2::element_text(
        size = 9.5, color = "black", face = "italic", hjust = 0.5,
        vjust = 1, margin = ggplot2::margin(0, 0, 0, 0)),
      strip.text = ggplot2::element_text(size = 5.5, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "transparent", colour = "transparent"),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 8),
      legend.position	= "inside",
      legend.position.inside = c(0.925, 0.85),
      legend.title = ggplot2::element_text(
        color = "black", size = 8, face = "bold"),
      legend.spacing.x = grid::unit(0.2, "cm"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE,
      panel.spacing = grid::unit(0.05, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA))

  EurBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  MapLabel <- c(
    "Number of plant observations",
    "Number of plant observations (native species)",
    "Number of plant species",
    "Number of native species")

  # # ..................................................................... ###

  Efforts_GBIF_Plots <- purrr::map(
    .x = seq_len(4),
    .f = ~{
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(
          data = Efforts_SummaryR[[.x]], maxcell = Inf) +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(title = MapLabel[.x], fill = NULL) +
        PlottingTheme
    }) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  Efforts_GBIF_Plots_Log <- purrr::map(
    .x = seq_len(4),
    .f  = ~{
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(
          data = log10(Efforts_SummaryR[[.x]]), maxcell = Inf) +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(
          title = paste0(MapLabel[.x], " - log10 scale"), fill = NULL) +
        PlottingTheme
    }) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  plots <- list(
    list(Efforts_GBIF_Plots$NObs, Efforts_GBIF_Plots_Log$NObs),
    list(Efforts_GBIF_Plots$NObs_Native, Efforts_GBIF_Plots_Log$NObs_Native),
    list(Efforts_GBIF_Plots$NSp, Efforts_GBIF_Plots_Log$NSp),
    list(Efforts_GBIF_Plots$NSp_Native, Efforts_GBIF_Plots_Log$NSp_Native))

  filenames <- c(
    "Efforts_GBIF_NObs.jpeg", "Efforts_GBIF_NObs_Native.jpeg",
    "Efforts_GBIF_NSp.jpeg", "Efforts_GBIF_NSp_Native.jpeg")

  purrr::walk2(
    .x = plots, .y = filenames,
    .f = ~{
      patchwork::wrap_plots(.x, ncol = 2, nrow = 1) %>%
        ggplot2::ggsave(
          filename = file.path(Path_Efforts, .y),
          width = 31, height = 16, units = "cm", dpi = 600)
    })

  # # ..................................................................... ###

  return(invisible(NULL))
}
