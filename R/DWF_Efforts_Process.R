#' Process GBIF sampling effort data for the `IAS-pDT`
#'
#' Downloads and processes GBIF sampling effort data for vascular plants in
#' Europe, supporting the Invasive Alien Species prototype Digital Twin
#' (`IAS-pDT`). Orchestrated by `Efforts_Process()`, it uses helper functions to
#' request, download, split, summarize, and visualize data at the Order level.
#' The functions prepares raster maps for the number of vascular plant
#' observations and species per grid cell.
#'
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Renviron Character. Path to `.Renviron` file with GBIF credentials
#'   (`GBIF_EMAIL`, `GBIF_USER`, `GBIF_PWD`). Default: `".Renviron"`. The
#'   credentials must be in the format:
#'    - `GBIF_EMAIL=your_email`
#'    - `GBIF_USER=your_username`
#'    - `GBIF_PWD=your_password`
#' @param Request Logical. If `TRUE` (default), requests GBIF data; otherwise,
#'   loads existing data.
#' @param Download Logical. If `TRUE` (default), downloads and saves GBIF data;
#'   otherwise, skips download. Default: `TRUE`.
#' @param Boundaries Numeric vector (length 4). GBIF data bounds (Left, Right,
#'   Bottom, Top). Default: `c(-30, 50, 25, 75)`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6, except for `Efforts_Request`, which defaults to 3 with a
#'   maximum of 3.
#' @param StartYear Integer. Earliest year for GBIF records (matches CHELSA
#'   climate data). Default: `1981`.
#' @param ChunkSize Integer. Rows per chunk file. Default: `100000`.
#' @param DeleteChunks Logical. If `TRUE` (default), deletes chunk files
#'   post-processing.
#' @param DeleteProcessed Logical. If `TRUE` (default), removes raw GBIF files
#'   after processing (>22 GB).
#' @param Path_Zip Character. Path to zip file with CSV for splitting.
#'
#' @note
#' - `Efforts_Process()` is the main entry point for processing sampling effort
#' data.
#' - Time-intensive (>9 hours on 6-core Windows PC; GBIF request ~5 hours).
#' - Detects and processes only new/missing data by order.
#'
#' @section Functions details:
#' - **`Efforts_Process()`**: Manages the workflow for requesting, downloading,
#'   processing, and plotting GBIF vascular plant data.
#' - **`Efforts_Request()`**: Requests GBIF data by order in parallel. Stores
#'   results to disk.
#' - **`Efforts_Download()`**: Downloads GBIF data, validates files, and loads
#'   existing data if available. Returns a dataframe (`Efforts_AllRequests`)
#'   with paths.
#' - **`Efforts_Split()`**: Splits zipped CSV data by order into chunks, saving
#'   each separately.
#' - **`Efforts_Summarize()`**: Processes and summarizes data into `RData` and
#'   TIFF rasters.
#' - **`Efforts_Plot()`**: Plots observation efforts (raw and log10 scales).
#' @references Data source: <https://www.gbif.org>

## |------------------------------------------------------------------------| #
# Efforts_Process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 1
#' @export

Efforts_Process <- function(
    EnvFile = ".env", Renviron = ".Renviron", Request = TRUE, Download = TRUE,
    NCores = 6L, StartYear = 1981L, Boundaries = c(-30, 50, 25, 75),
    ChunkSize = 100000L, DeleteChunks = TRUE, DeleteProcessed = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Renviron", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical", Args = c("Request", "Download"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "Boundaries", "StartYear"))

  # Validate Boundaries argument
  if (length(Boundaries) != 4) {
    stop("Boundaries must be a numeric vector of length 4.", call. = FALSE)
  }

  # Validate ChunkSize
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive numeric value.", call. = FALSE)
  }

  # Validate NCores
  if (!is.numeric(NCores) || NCores <= 0 || NCores > 50) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  # Validate StartYear
  if (!is.numeric(StartYear) || StartYear <= 1950) {
    stop("StartYear must be a positive integer after 1950", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Efforts <- Path_Raw <- Path_Interim <- Path_Grid <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime("Ensure that GBIF access information is available")
  IASDT.R::GBIF_Check(Renviron = Renviron)

  # # ..................................................................... ###

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE,
    "Path_Raw", "DP_R_Efforts_raw", FALSE, FALSE,
    "Path_Interim", "DP_R_Efforts_interim", FALSE, FALSE,
    "Taxa_Stand", "DP_R_Taxa_stand", FALSE, TRUE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  AllRequests <- IASDT.R::Path(Path_Efforts, "Efforts_AllRequests.RData")

  # # ..................................................................... ###

  # Loading input data ------
  IASDT.R::CatTime("Loading input data")

  ## Create paths -----
  Path_Efforts_Requests <- IASDT.R::Path(Path_Efforts, "Requests")
  Path_Efforts_Cleaned <- IASDT.R::Path(Path_Interim, "CleanedData")
  # Create required directories
  fs::dir_create(
    c(
      Path_Efforts, Path_Raw, Path_Interim, Path_Efforts_Cleaned,
      Path_Efforts_Requests))

  ## Reference grid ----
  Grids <- Path_Grid %>%
    IASDT.R::Path(c("Grid_10_Land_Crop_sf.RData", "Grid_10_Land_Crop.RData"))

  missing_grids <- Grids[!file.exists(Grids)]
  if (length(missing_grids) > 0) {
    stop(
      "The following grid file(s) do not exist:\n",
      paste0(" >>> ", missing_grids, collapse = "\n"), call. = FALSE)
  }

  # # ..................................................................... ###

  # Request efforts data ------
  
  if (Request) {

    IASDT.R::CatTime("Requesting efforts data")
    
    IASDT.R::Efforts_Request(
      EnvFile = EnvFile, NCores = NCores, StartYear = StartYear,
      Renviron = Renviron, Boundaries = Boundaries)

  } else {

    if (!file.exists(AllRequests)) {
      stop(
        "Efforts data was not requested and the file does not exist.",
        call. = FALSE)
    }

    IASDT.R::CatTime(
      "Efforts data was not requested, but already available on disk")

  }

  # # ..................................................................... ###

  # Download efforts data ------

  if (Download) {

    IASDT.R::CatTime("Downloading efforts data")
    IASDT.R::Efforts_Download(NCores = NCores, EnvFile = EnvFile)

  } else {

    if (!file.exists(AllRequests)) {
      stop(
        "Efforts data was not downloaded and the file does not exist.",
        call. = FALSE)
    }

    Efforts_AllRequests <- IASDT.R::LoadAs(AllRequests)

    if (!("DownPath" %in% names(Efforts_AllRequests))) {
      stop(
        "Efforts data was not downloaded and the 'DownPath' column is missing.",
        call. = FALSE)
    }

    rm(Efforts_AllRequests, envir = environment())
    invisible(gc())

    IASDT.R::CatTime("Efforts data was not downloaded")

  }

  # # ..................................................................... ###

  # Processing efforts data ------

  IASDT.R::CatTime("Processing efforts data")
  IASDT.R::Efforts_Summarize(
    EnvFile = EnvFile, NCores = NCores, ChunkSize = ChunkSize,
    DeleteChunks = DeleteChunks)

  # # ..................................................................... ###

  # # Plotting ----

  IASDT.R::CatTime("Plotting sampling efforts")
  IASDT.R::Efforts_Plot(EnvFile = EnvFile)

  # # ..................................................................... ###

  # # Cleaning up ----

  if (DeleteProcessed) {
    IASDT.R::CatTime("Cleaning up - delete downloaded GBIF data")
    fs::file_delete(list.files(Path_Raw, full.names = TRUE))
    fs::dir_delete(Path_Raw)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing efforts data took ", ... = "\n")

  return(invisible(NULL))
}
