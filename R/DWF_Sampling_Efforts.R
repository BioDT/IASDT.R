## |------------------------------------------------------------------------| #
# Sampling_Efforts ----
## |------------------------------------------------------------------------| #

#' Download and Process sampling efforts data
#'
#' This function downloads GBIF data for all vascular plants within a specified
#' geographical area (Europe), grouped by order, and converts the data to raster
#' format to represent the number of vascular plant observations and species per
#' grid cell.
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
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. See [Efforts_Split] and [Efforts_Process] for more details.
#' @note This function is expected to take a substantial amount of time (>9
#'   hours on a Windows PC with 6 cores). The data request from GBIF may take
#'   around 5 hours to be ready. The function requests GBIF data for each
#'   vascular plant order and waits for the data to be ready before processing
#'   them.
#' @return Returns `NULL` invisibly. The function generates various output
#'   files, maps, and logs, and it is designed to be used for its side effects.
#' @author Ahmed El-Gabbas
#' @name Sampling_Efforts
#' @export

Sampling_Efforts <- function(
    FromHPC = TRUE, EnvFile = ".env", Renviron = ".Renviron",
    RequestData = TRUE, DownloadData = TRUE, NCores = 6, StartYear = 1980,
    Boundaries = c(-30, 50, 25, 75), ChunkSize = 10000) {

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

  # Validate Boundaries argument
  if (length(Boundaries) != 4) {
    stop("Boundaries must be a numeric vector of length 4.")
  }

  # Validate ChunkSize
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive numeric value.")
  }

  # Validate NCores
  if (!is.numeric(NCores) || NCores <= 0 || NCores > 50) {
    stop("NCores must be a positive integer.")
  }

  # Validate StartYear
  if (!is.numeric(StartYear) || StartYear <= 1950) {
    stop("StartYear must be a positive integer after 1950")
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Efforts <- Path_Efforts_Raw <- Path_Efforts_Interim <- Path_Grid <-
    Taxa_Stand <- EU_Bound <- NULL

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
  # Create required directories
  fs::dir_create(c(
    Path_Efforts, Path_Efforts_Raw, Path_Efforts_Interim, Path_Efforts_Data,
    Path_Efforts_Requests))

  ## Reference grid ----
  Grids <- Path_Grid %>%
    file.path(c("Grid_10_Land_Crop_sf.RData", "Grid_10_Land_Crop.RData"))

  missing_grids <- Grids[!file.exists(Grids)]
  if (length(missing_grids) > 0) {
    stop(
      paste0("The following grid file(s) do not exist:\n",
             paste0(" >>> ", missing_grids, collapse = "\n")),
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
    Path_Efforts_Interim = Path_Efforts_Interim,
    Path_Efforts_Data = Path_Efforts_Data,
    Path_Grid = Path_Grid, IAS_List = IAS_List,
    Efforts_AllRequests = Efforts_AllRequests, ChunkSize = ChunkSize)

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
#' @param NCores Integer. The number of cores to use for parallel processing
#'   (between 1 and 3). Must be a positive integer.
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

  if (missing(NCores) || !is.numeric(NCores) || NCores < 1 || NCores > 3) {
    stop("`NCores` must be a positive integer between 1 and 3.", call. = FALSE)
  }

  if (!is.character(Path_Requests) || !dir.exists(Path_Requests)) {
    stop("`Path_Requests` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Efforts) || !dir.exists(Path_Efforts)) {
    stop("`Path_Efforts` must be a valid directory path.", call. = FALSE)
  }

  if (!is.numeric(StartYear) || StartYear <= 1950) {
    stop("`StartYear` must be a positive integer after 1950")
  }

  if (!is.numeric(Boundaries) || length(Boundaries) != 4) {
    stop("`Boundaries` must be a numeric vector of length 4.", call. = FALSE)
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
#' @param NCores Integer. Number of cores to use for parallel processing.  Must
#'   be a positive integer. This directory must exist or be created beforehand.
#' @param Path_Raw Character. Path where the raw downloaded data will be saved.
#'   This directory must exist or be created beforehand.
#' @param Path_Interim Character. Path where the interim CSV files will be
#'   saved. This directory must exist or be created beforehand.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved. This directory must exist or be created beforehand.
#' @name Efforts_Download
#' @author Ahmed El-Gabbas
#' @return A data frame (`Efforts_AllRequests`) with updated download paths and
#'   interim file paths.
#' @export

Efforts_Download <- function(NCores = 6, Path_Raw, Path_Interim, Path_Efforts) {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  DownPath <- Request <- NULL

  # # ..................................................................... ###

  # Validate NCores
  if (missing(NCores) || !is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  # Validate directory paths
  if (!is.character(Path_Raw) || !dir.exists(Path_Raw)) {
    stop("`Path_Raw` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Interim) || !dir.exists(Path_Interim)) {
    stop("`Path_Interim` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Efforts) || !dir.exists(Path_Efforts)) {
    stop("`Path_Efforts` must be a valid directory path.", call. = FALSE)
  }

  IASDT.R::CheckCommands("unzip")

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
            FileOkay <- tryCatch({
              system2(
                "unzip", args = c("-t", DownFile), stdout = TRUE, stderr = TRUE)
            }, error = function(e) {
              message("Error during file validation: ", conditionMessage(e))
              return(NULL)
            }) %>%
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

        }, .options = furrr::furrr_options(seed = TRUE)))

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
#' @param Path_Efforts_Interim Character. The directory path to save interim
#'   data.
#' @param Path_Efforts_Data Character. The directory path to save detailed
#'   effort data as `RData`.
#' @param Path_Grid Character. The directory path to load the grid data.
#' @param IAS_List A list of invasive alien species keys.
#' @param Efforts_AllRequests A data frame containing the details of the GBIF
#'   download, including paths to CSV files, zip files, order, class, and total
#'   records.
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. See [Efforts_Split] for more details.
#' @return The function returns and saves the GBIF data summary.
#' @author Ahmed El-Gabbas
#' @name Efforts_Process
#' @export

Efforts_Process <- function(
    NCores, Path_Efforts, Path_Efforts_Interim, Path_Efforts_Data, Path_Grid,
    IAS_List, Efforts_AllRequests, ChunkSize = 100000) {

  .StartTimeProcess <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  if (!is.numeric(NCores) || length(NCores) != 1 || NCores <= 0) {
    stop("NCores must be a single positive integer.", call. = FALSE)
  }

  if (!is.character(Path_Efforts_Data) || length(Path_Efforts_Data) != 1 ||
      !dir.exists(Path_Efforts_Data)) {
    stop(
      "Path_Efforts_Data must be a single valid directory path.",
      call. = FALSE)
  }

  if (missing(Path_Grid) || !is.character(Path_Grid) ||
      !dir.exists(Path_Grid)) {
    stop("Path_Grid must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Efforts_Interim) || !is.character(Path_Efforts_Interim) ||
      !dir.exists(Path_Efforts_Interim)) {
    stop("Path_Efforts_Interim must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Efforts) || !is.character(Path_Efforts) ||
      !dir.exists(Path_Efforts)) {
    stop("Path_Efforts must be a valid directory path.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  speciesKey <- CellCode <- ObsN <- year <- UncertainKm <- Latitude <-
    Longitude <- taxonRank <- ID <- DownPath <- Chunks <- TotalRecords <-
    ClassOrder <- Path_DT <- NULL

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
  withr::local_options(future.globals.maxSize = 8000 * 1024^2)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  invisible(
    snow::clusterEvalQ(
      cl = c1,
      IASDT.R::LoadPackages(
        List = c("terra", "IASDT.R", "stringr", "fs", "sf", "readr", "dplyr"))))

  snow::clusterExport(
    cl = c1,
    list = c(
      "Path_Efforts", "Path_Efforts_Interim", "Efforts_AllRequests",
      "Path_Grid_R", "Path_Efforts_Data", "Grid_SF", "IAS_List", "ChunkSize"),
    envir = environment())

  # # ..................................................................... ###

  # Reading data from the zipped archive -----
  IASDT.R::CatTime("Reading data from the zipped archive", Level = 1)

  Efforts_Summary <- Efforts_AllRequests %>%
    dplyr::select(
      tidyselect::all_of(
        c("class", "order", "Request", "DownLink",
          "TotalRecords", "DownPath"))) %>%
    dplyr::mutate(
      ClassOrder = purrr::map2_chr(class, order, ~paste0(.x, "_", .y)),
      Path_DT = furrr::future_pmap_chr(
        .l = list(DownPath, TotalRecords, ClassOrder),
        .f = function(DownPath, TotalRecords, ClassOrder) {

          Path_DT <- file.path(Path_Efforts_Data, paste0(ClassOrder, ".RData"))

          if (file.exists(Path_DT)) {
            return(Path_DT)
          }

          # Split data into chunks
          Chunks <- IASDT.R::Efforts_Split(
            Path_Zip = DownPath, Path_Output = Path_Efforts_Interim,
            ChunkSize = ChunkSize)

          if (TotalRecords == 0 || length(Chunks) == 0) {

            return(NA_character_)

          } else {

            DT <- readr::read_tsv(
              file = Chunks,
              col_names = c(
                "taxonRank", "Latitude", "Longitude",
                "UncertainKm", "year", "speciesKey"),
              progress = FALSE, show_col_types = FALSE,
              col_types = readr::cols(
                UncertainKm = readr::col_double(),
                Longitude = readr::col_double(),
                Latitude = readr::col_double(),
                year = readr::col_integer(),
                speciesKey = readr::col_integer(),
                taxonRank = readr::col_character(),
                .default = readr::col_double())) %>%
              dplyr::mutate(UncertainKm = UncertainKm / 1000) %>%
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

            IASDT.R::SaveAs(InObj = DT, OutObj = ClassOrder, OutPath = Path_DT)

            fs::file_delete(Chunks)
            return(Path_DT)
          }
        },
        .options = furrr::furrr_options(seed = TRUE)
      ))

  # # ..................................................................... ###

  # Number of observations and species per order ----
  IASDT.R::CatTime("Number of observations and species per order", Level = 1)

  Efforts_Summary <- Efforts_Summary %>%
    dplyr::mutate(
      NObs_NSp_R = furrr::future_map2(
        .x = Path_DT, .y = ClassOrder,
        .f = ~{

          if (is.na(.x)) {
            ObsN <- 0
          } else {
            DT <- IASDT.R::LoadAs(.x)
            ObsN <- nrow(DT)
          }

          Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

          if (ObsN == 0) {
            NObs_R <- NSp_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_R <- setNames(NObs_R, paste0("NObs_", .y)) %>%
              terra::wrap()
            NSp_R <- setNames(NSp_R, paste0("NSp_", .y)) %>%
              terra::wrap()
          } else {

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
              stats::setNames(paste0("NObs_", .y)) %>%
              terra::wrap()

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
              stats::setNames(paste0("NSp_", .y)) %>%
              terra::wrap()
          }

          invisible(gc())

          tibble::tibble(
            ObsN = ObsN, NObs_R = list(NObs_R), NSp_R = list(NSp_R)) %>%
            return()
        },
        .options = furrr::furrr_options(seed = TRUE))) %>%
    tidyr::unnest_wider("NObs_NSp_R") %>%
    tidyr::unnest(c("NObs_R", "NSp_R"))

  # # ..................................................................... ###

  # Number of observations and species per order - native species ----
  IASDT.R::CatTime(
    "Number of observations and species per order - native species", Level = 1)

  Efforts_Summary <- Efforts_Summary %>%
    dplyr::mutate(
      Native = purrr::map2(
        .x = Path_DT, .y = ClassOrder,
        .f = ~{

          if (is.na(.x)) {
            ObsN_Native <- 0
          } else {
            # Data on native species
            DT_Native <- IASDT.R::LoadAs(.x) %>%
              dplyr::filter(!(speciesKey %in% IAS_List))
            ObsN_Native <- nrow(DT_Native)
          }

          Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

          if (ObsN_Native == 0) {
            NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_Native_R <- setNames(
              NObs_Native_R, paste0("NObsNative_", .y)) %>%
              terra::wrap()
            NSp_Native_R <- setNames(
              NSp_Native_R, paste0("NSpNative_", .y)) %>%
              terra::wrap()

          } else {

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
              stats::setNames(paste0("NObsNative_", .y)) %>%
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
              stats::setNames(paste0("NSpNative_", .y)) %>%
              terra::wrap()
          }

          tibble::tibble(
            ObsN_Native = ObsN_Native, NObs_Native_R = list(NObs_Native_R),
            NSp_Native_R = list(NSp_Native_R)) %>%
            return()
        })
    ) %>%
    tidyr::unnest_wider("Native") %>%
    tidyr::unnest(c("NObs_Native_R", "NSp_Native_R"))

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

  # # ..................................................................... ###

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
#' @param Path_Efforts Character. Path to the directory where the generated
#'   plots will be saved. The directory must exist.
#' @param EU_Bound Character. Path to `RData` file containing country
#'   boundaries.
#' @return The function saves generated plots as JPEG files in the specified
#'   directory and returns NULL invisibly.
#' @author Ahmed El-Gabbas
#' @name Efforts_Plot
#' @details This function generates and saves effort maps visualizing the number
#'   of plant observations and species, including both native and non-native
#'   species, within Europe. It produces both standard and log10-scaled plots.
#' @export

Efforts_Plot <- function(Path_Efforts, EU_Bound) {

  # # ..................................................................... ###

  File_SummaryR <- file.path(Path_Efforts, "Efforts_SummaryR.RData")
  if (!file.exists(File_SummaryR)) {
    stop(paste0("Summary maps cannot be loaded: ", File_SummaryR),
         call. = FALSE)
  }

  Efforts_SummaryR <- terra::unwrap(IASDT.R::LoadAs(File_SummaryR))

  # # ..................................................................... ###

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.02, 0, 0.02, 0, "cm"),
      plot.title = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 8),
      legend.position	= "inside",
      legend.position.inside = c(0.925, 0.85),
      legend.title = ggplot2::element_text(
        color = "black", size = 8, face = "bold"),
      # legend.spacing.x = grid::unit(0.2, "cm"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA)
    )

  EurBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # # ..................................................................... ###

  Efforts_Plots <- purrr::map(
    .x = seq_len(terra::nlyr(Efforts_SummaryR)),
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
        ggplot2::labs(fill = NULL) +
        PlottingTheme
    }) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  Efforts_Plots_Log <- purrr::map(
    .x = seq_len(terra::nlyr(Efforts_SummaryR)),
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
        ggplot2::labs(fill = "log10") +
        PlottingTheme
    }) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  PlotDF <- tibble::tibble(
    Plots = list(
      list(Efforts_Plots$NObs, Efforts_Plots_Log$NObs),
      list(Efforts_Plots$NObs_Native, Efforts_Plots_Log$NObs_Native),
      list(Efforts_Plots$NSp, Efforts_Plots_Log$NSp),
      list(Efforts_Plots$NSp_Native, Efforts_Plots_Log$NSp_Native))) %>%
    dplyr::mutate(
      FileName = c(
        "Efforts_GBIF_NObs.jpeg", "Efforts_GBIF_NObs_Native.jpeg",
        "Efforts_GBIF_NSp.jpeg", "Efforts_GBIF_NSp_Native.jpeg"),
      Title =  c(
        "Number of plant observations",
        "Number of plant observations (native species)",
        "Number of plant species",
        "Number of native species"))

  purrr::walk(
    .x = seq_len(nrow(PlotDF)),
    .f = ~{
      CurrPlot <- patchwork::wrap_plots(
        PlotDF$Plots[[.x]], ncol = 2, nrow = 1) +
        patchwork::plot_annotation(
          title = PlotDF$Title[[.x]],
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 14, face = "bold", hjust = 0.5, colour = "blue",
              margin = ggplot2::margin(0.25, 0, 0.5, 0))))

      ggplot2::ggsave(
        plot = CurrPlot,
        filename = file.path(Path_Efforts, PlotDF$FileName[[.x]]),
        width = 31, height = 16.25, units = "cm", dpi = 600)
    })

  # # ..................................................................... ###

  return(invisible(NULL))
}

## ............................................ ------
## ............................................ ------

## |------------------------------------------------------------------------| #
# Efforts_Split ----
## |------------------------------------------------------------------------| #

#' Split Order data into smaller chunks
#'
#' This function extracts Order data without extracting the zipped archive. The
#' function reads CSV files inside the zipped file, selects specified columns,
#' and splits the data into smaller chunks of specified row size, saving each
#' chunk as a separate file.
#' @param Path_Zip Character. The file path to the zip archive containing the
#'   CSV file. The file must be a ZIP archive containing a single CSV file.
#' @param Path_Output Character. The directory where the split files will be
#'   saved. The directory must exist.
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. Note: Larger chunk sizes may require significant memory and
#'   processing power.
#' @return A character vector of file path(s) to the created chunk files.
#' @name Efforts_Split
#' @author Ahmed El-Gabbas
#' @export

Efforts_Split <- function(Path_Zip, Path_Output, ChunkSize = 100000) {

  ID <- Col <- NULL

  # Check if ChunkSize is valid (greater than zero)
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive number.")
  }

  # Checking required bash tools
  IASDT.R::CheckCommands(c("unzip", "cut", "sed", "split"))

  # ensure that `ChunkSize` is not formatted in scientific notation
  ChunkSize <- format(ChunkSize, scientific = FALSE)

  CSV_File <- stringr::str_replace(basename(Path_Zip), ".zip$", ".csv")
  OutPrefix <- stringr::str_replace(basename(Path_Zip), ".zip$", "_") %>%
    file.path(Path_Output, .)

  # extract column names and their numbers from the zipped file without
  # extraction read first line
  SelectedCols <- stringr::str_glue(
    "unzip -p {Path_Zip} {CSV_File} | head -n 1") %>%
    IASDT.R::System() %>%
    # Split the first row into column names. Data is tab-separated
    stringr::str_split("\t") %>%
    magrittr::extract2(1) %>%
    dplyr::tibble(Col = .) %>%
    # column number in the original data
    dplyr::mutate(ID = seq_len(dplyr::n())) %>%
    # Only keep selected columns
    dplyr::filter(
      Col %in% c(
        "taxonRank", "decimalLatitude", "decimalLongitude",
        "coordinateUncertaintyInMeters", "year", "speciesKey")) %>%
    dplyr::pull(ID) %>%
    paste0(collapse = ",")

  Command <- stringr::str_glue(
    'unzip -p {Path_Zip} {CSV_File} | cut -f{SelectedCols} -d "\t" | ',
    'sed -n "1!p" | split -l {ChunkSize} ',
    "-a 3 -d - {OutPrefix} --additional-suffix=.txt")

  Path_Chunks <- tryCatch({
    IASDT.R::System(Command, RObj = FALSE)
  }, error = function(e) {
    stop("Failed to execute system command: ", e$message)
  })

  list.files(
    Path_Output, pattern = paste0(basename(OutPrefix), ".+txt"),
    full.names = TRUE) %>%
    return()
}

## ............................................ ------
## ............................................ ------
