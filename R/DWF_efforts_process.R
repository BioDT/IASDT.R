#' Process GBIF sampling effort data for the `IASDT`
#'
#' Downloads and processes GBIF sampling effort data for vascular plants in
#' Europe, supporting the Invasive Alien Species Digital Twin (`IASDT`).
#' Orchestrated by `efforts_process()`, it uses helper functions to request,
#' download, split, summarise, and visualise data at the Order level. The
#' functions prepares raster maps for the number of vascular plant observations
#' and species per grid cell.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param r_environ Character. Path to `.Renviron` file with GBIF credentials
#'   (`GBIF_EMAIL`, `GBIF_USER`, `GBIF_PWD`). Default: `".Renviron"`. The
#'   credentials must be in the format:
#'    - `GBIF_EMAIL=your_email`
#'    - `GBIF_USER=your_username`
#'    - `GBIF_PWD=your_password`
#' @param request Logical. If `TRUE` (default), requests GBIF data; otherwise,
#'   loads existing data.
#' @param download Logical. If `TRUE` (default), downloads and saves GBIF data;
#'   otherwise, skips download. Default: `TRUE`.
#' @param boundaries Numeric vector (length 4). GBIF data bounds (Left, Right,
#'   Bottom, Top). Default: `c(-30, 50, 25, 75)`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6, except for `efforts_request`, which defaults to 3 with a
#'   maximum of 3.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "future::sequential", "future::multisession",
#'   "future::multicore", and "future::cluster". Defaults to
#'   `"future::multicore"` (`"future::multisession"` on Windows). See
#'   [future::plan()] and [ecokit::set_parallel()] for details.
#' @param start_year Integer. Earliest year for GBIF records (matches CHELSA
#'   climate data). Default: `1981`.
#' @param chunk_size Integer. Rows per chunk file. Default: `100000`.
#' @param delete_chunks Logical. If `TRUE` (default), deletes chunk files
#'   post-processing.
#' @param delete_processed Logical. If `TRUE` (default), removes raw GBIF files
#'   after processing (>22 GB).
#' @param path_zip Character. Path to zip file with CSV for splitting.
#'
#' @note
#' - `efforts_process()` is the main entry point for processing sampling effort
#' data.
#' - Time-intensive (>9 hours on 6-core Windows PC; GBIF request ~5 hours).
#' - Detects and processes only new/missing data by order.
#'
#' @section Functions details:
#' - **`efforts_process()`**: Manages the workflow for requesting, downloading,
#'   processing, and plotting GBIF vascular plant data.
#' - **`efforts_request()`**: Requests GBIF data by order in parallel. Stores
#'   results to disk.
#' - **`efforts_download()`**: Downloads GBIF data, validates files, and loads
#'   existing data if available. Returns a dataframe (`Efforts_AllRequests`)
#'   with paths.
#' - **`efforts_split()`**: Splits zipped CSV data by order into chunks, saving
#'   each separately.
#' - **`efforts_summarize()`**: Processes and summarises data into `RData` and
#'   TIFF rasters.
#' - **`efforts_plot()`**: Plots observation efforts (raw and log10 scales).
#' @references Data source: <https://www.gbif.org>

## |------------------------------------------------------------------------| #
# efforts_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 1
#' @export

efforts_process <- function(
    env_file = ".env", r_environ = ".Renviron", request = TRUE, download = TRUE,
    n_cores = 6L, strategy = "future::multicore", start_year = 1981L,
    boundaries = c(-30, 50, 25, 75), chunk_size = 100000L,
    delete_chunks = TRUE, delete_processed = TRUE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("r_environ", "env_file", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("request", "download"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "boundaries", "start_year"))

  # Validate boundaries argument
  if (length(boundaries) != 4) {
    ecokit::stop_ctx(
      "`boundaries` must be a numeric vector of length 4.",
      boundaries = boundaries, length_boundaries = length(boundaries),
      include_backtrace = TRUE)
  }

  # Validate chunk_size
  if (!is.numeric(chunk_size) || chunk_size <= 0) {
    ecokit::stop_ctx(
      "`chunk_size` must be a positive numeric value.",
      chunk_size = chunk_size, include_backtrace = TRUE)
  }

  # Validate n_cores
  if (!is.numeric(n_cores) || n_cores <= 0 || n_cores > 50) {
    ecokit::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (strategy == "future::sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c(
    "future::sequential", "future::multisession", "future::multicore",
    "future::cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  # Validate start_year
  if (!is.numeric(start_year) || start_year <= 1950) {
    ecokit::stop_ctx(
      "`start_year` must be a positive integer after 1950",
      start_year = start_year, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Efforts <- Path_Raw <- Path_Interim <- Path_Grid <- NULL

  # # ..................................................................... ###

  ecokit::cat_time("Ensure that GBIF access information is available")
  IASDT.R::GBIF_check(r_environ = r_environ)

  # # ..................................................................... ###

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE,
    "Path_Raw", "DP_R_Efforts_raw", FALSE, FALSE,
    "Path_Interim", "DP_R_Efforts_interim", FALSE, FALSE,
    "Taxa_Stand", "DP_R_Taxa_stand", FALSE, TRUE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  AllRequests <- fs::path(Path_Efforts, "Efforts_AllRequests.RData")

  # # ..................................................................... ###

  # Loading input data ------
  ecokit::cat_time("Loading input data")

  ## Create paths -----
  Path_Efforts_Requests <- fs::path(Path_Efforts, "Requests")
  Path_Efforts_Cleaned <- fs::path(Path_Interim, "CleanedData")
  # Create required directories
  fs::dir_create(
    c(
      Path_Efforts, Path_Raw, Path_Interim, Path_Efforts_Cleaned,
      Path_Efforts_Requests))

  ## Reference grid ----
  Grids <- Path_Grid %>%
    fs::path(c("Grid_10_Land_Crop_sf.RData", "Grid_10_Land_Crop.RData"))

  missing_grids <- Grids[!file.exists(Grids)]
  if (length(missing_grids) > 0) {
    ecokit::stop_ctx(
      paste0(
        "The following grid file(s) do not exist:\n",
        paste0(" >>> ", missing_grids, collapse = "\n")),
      missing_grids = missing_grids,
      length_missing_grids = length(missing_grids), include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # request efforts data ------

  if (request) {

    ecokit::cat_time("Requesting efforts data")

    IASDT.R::efforts_request(
      env_file = env_file, n_cores = n_cores, start_year = start_year,
      r_environ = r_environ, boundaries = boundaries)

  } else {

    if (!file.exists(AllRequests)) {
      ecokit::stop_ctx(
        "Efforts data was not requested and the file does not exist.",
        AllRequests = AllRequests, include_backtrace = TRUE)
    }

    ecokit::cat_time(
      stringr::str_glue(
        "Efforts data was not requested. Previous request information is \\
        already available on disk"))
  }

  # # ..................................................................... ###

  # download efforts data ------

  if (download) {

    ecokit::cat_time("Downloading efforts data")
    IASDT.R::efforts_download(
      n_cores = n_cores, strategy = strategy, env_file = env_file)

  } else {

    if (!file.exists(AllRequests)) {
      ecokit::stop_ctx(
        "Efforts data was not downloaded and the file does not exist.",
        AllRequests = AllRequests, include_backtrace = TRUE)
    }

    Efforts_AllRequests <- ecokit::load_as(AllRequests)

    if (!("DownPath" %in% names(Efforts_AllRequests))) {
      ecokit::stop_ctx(
        "Efforts data was not downloaded and the 'DownPath' column is missing.",
        Efforts_AllRequests = Efforts_AllRequests,
        names_Efforts_AllRequests = names(Efforts_AllRequests),
        include_backtrace = TRUE)
    }

    rm(Efforts_AllRequests, envir = environment())
    invisible(gc())

    ecokit::cat_time("Efforts data was not downloaded")

  }

  # # ..................................................................... ###

  # Processing efforts data ------

  ecokit::cat_time("Processing efforts data")
  IASDT.R::efforts_summarize(
    env_file = env_file, n_cores = n_cores, strategy = strategy,
    chunk_size = chunk_size, delete_chunks = delete_chunks)

  # # ..................................................................... ###

  # # Plotting ----

  ecokit::cat_time("Plotting sampling efforts")
  IASDT.R::efforts_plot(env_file = env_file)

  # # ..................................................................... ###

  # # Cleaning up ----

  if (delete_processed) {
    ecokit::cat_time("Cleaning up - delete downloaded GBIF data")
    fs::file_delete(list.files(Path_Raw, full.names = TRUE))
    fs::dir_delete(Path_Raw)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing efforts data took ", ... = "\n")

  return(invisible(NULL))
}
