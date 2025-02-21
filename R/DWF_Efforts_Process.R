## |------------------------------------------------------------------------| #
# Efforts_Process ----
## |------------------------------------------------------------------------| #

#' Download and Process sampling efforts data
#'
#' This function downloads GBIF data for all vascular plants within a specified
#' geographical area (Europe), grouped by order, and converts the data to raster
#' format to represent the number of vascular plant observations and species per
#' grid cell.
#' @param FromHPC Logical. Whether the processing is being done on an
#'   High-Performance Computing (HPC) environment, to adjust file paths
#'   accordingly. Default: `TRUE`.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Renviron Character. The path to the `.Renviron` file containing GBIF
#'   login credentials (email, user, password). Default: `.Renviron`.
#' @param RequestData Logical. If `TRUE`, the function requests data from GBIF.
#'   If `FALSE`, previously requested data is loaded from disk. Defaults to
#'   `TRUE`.
#' @param DownloadData Logical. If `TRUE`, the function downloads data and
#'   stores it on disk. If `FALSE`, it skips the download step. Defaults to
#'   `TRUE`.
#' @param Boundaries Numeric vector of length 4. The boundaries of the 
#'   requested GBIF data in the order: Left, Right, Bottom, Top. Defaults to 
#'   `c(-30, 50, 25, 75)`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6.
#' @param StartYear Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be requested from GBIF. Default is
#'   `1981`, which matches the year ranges of CHELSA current climate data.
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. See [Efforts_Split] and [Efforts_Summarize] for more details.
#' @param DeleteChunks Logical. Whether to remove file chunks after processing
#'   the data. Defaults to `TRUE`.
#' @param DeleteProcessed Logical. Whether to delete the raw downloaded GBIF
#'   data after processing them. This helps to free large unnecessary file space
#'   (> 22 GB). Defaults to `TRUE`.
#' @note
#' - This function is expected to take a substantial amount of time (>9
#' hours on a Windows PC with 6 cores). The data request from GBIF may take
#' around 5 hours to be ready. The function requests GBIF data for each vascular
#' plant order and waits for the data to be ready before processing them.
#' - This function should be the only function to be called to prepare sampling
#' efforts data. It calls other functions [Efforts_Request] to request data from
#' GBIF, [Efforts_Download] to download zipped archive for each vascular plant
#' order, [Efforts_Summarize] and [Efforts_Split] to process data in small
#' chunks, and [Efforts_Plot] for plotting.
#' @return Returns `NULL` invisibly. The function generates various output
#'   files, maps, and logs, and it is designed to be used for its side effects.
#' @author Ahmed El-Gabbas
#' @name Efforts_Process
#' @export

Efforts_Process <- function(
    FromHPC = TRUE, EnvFile = ".env", Renviron = ".Renviron",
    RequestData = TRUE, DownloadData = TRUE, NCores = 6L, StartYear = 1981L,
    Boundaries = c(-30, 50, 25, 75), ChunkSize = 100000L,
    DeleteChunks = TRUE, DeleteProcessed = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Renviron", "EnvFile")
  )
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("FromHPC", "RequestData", "DownloadData")
  )
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "Boundaries", "StartYear")
  )

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
    "Ensure that GBIF access information is available", Level = 1)
  IASDT.R::GBIF_Check(Renviron = Renviron)

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

  # Loading input data ------
  IASDT.R::CatTime("Loading input data")

  ## Create paths -----
  Path_Efforts_Requests <- IASDT.R::Path(Path_Efforts, "Requests")
  Path_Efforts_Data <- IASDT.R::Path(Path_Efforts_Interim, "CleanedData")
  # Create required directories
  fs::dir_create(
    c(
      Path_Efforts, Path_Efforts_Raw, Path_Efforts_Interim, Path_Efforts_Data,
      Path_Efforts_Requests))

  ## Reference grid ----
  Grids <- Path_Grid %>%
    IASDT.R::Path(c("Grid_10_Land_Crop_sf.RData", "Grid_10_Land_Crop.RData"))

  missing_grids <- Grids[!file.exists(Grids)]
  if (length(missing_grids) > 0) {
    stop(
      paste0(
        "The following grid file(s) do not exist:\n",
        paste0(" >>> ", missing_grids, collapse = "\n")
      ),
      call. = FALSE)
  }

  ## IAS species list ----
  IAS_List <- unique(dplyr::pull(readRDS(Taxa_Stand), "speciesKey"))

  # # ..................................................................... ###

  # Request efforts data ------
  IASDT.R::CatTime("Request efforts data")

  if (RequestData) {
    IASDT.R::CatTime("Requesting efforts data", Level = 1)

    Efforts_AllRequests <- IASDT.R::Efforts_Request(
      NCores = NCores, Path_Requests = Path_Efforts_Requests,
      Path_Efforts = Path_Efforts, StartYear = StartYear,
      Boundaries = Boundaries)
  } else {
    IASDT.R::CatTime("Efforts data was not requested, but loaded", Level = 1)
    Efforts_AllRequests <- IASDT.R::LoadAs(
      IASDT.R::Path(Path_Efforts, "Efforts_AllRequests.RData"))
  }

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
      IASDT.R::Path(Path_Efforts, "Efforts_AllRequests.RData"))
  }

  # # ..................................................................... ###

  # Processing efforts data ------

  IASDT.R::CatTime("Processing efforts data")
  IASDT.R::Efforts_Summarize(
    NCores = NCores, Path_Efforts = Path_Efforts,
    Path_Efforts_Interim = Path_Efforts_Interim,
    Path_Efforts_Data = Path_Efforts_Data,
    Path_Grid = Path_Grid, IAS_List = IAS_List,
    Efforts_AllRequests = Efforts_AllRequests, ChunkSize = ChunkSize,
    DeleteChunks = DeleteChunks)

  # # ..................................................................... ###

  # # Plotting ----

  IASDT.R::CatTime("Plotting sampling efforts")
  IASDT.R::Efforts_Plot(Path_Efforts = Path_Efforts, EU_Bound = EU_Bound)

  # # ..................................................................... ###

  # # Cleaning up ----

  if (DeleteProcessed) {
    IASDT.R::CatTime("Cleaning up - delete downloaded GBIF data")
    fs::file_delete(list.files(Path_Efforts_Raw, full.names = TRUE))
    fs::dir_delete(Path_Efforts_Raw)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing efforts data took ", ... = "\n")

  return(invisible(NULL))
}
