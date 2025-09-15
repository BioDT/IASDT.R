#' Process CHELSA Climate Data for the `IASDT`
#'
#' Downloads, processes, and projects [CHELSA](https://chelsa-climate.org/)
#' climate data at the European scale for the Invasive Alien Species Digital
#' Twin (`IASDT`). Supports multiple climate scenarios, outputting data in TIFF
#' and NetCDF formats. Orchestrated by `CHELSA_process()`, with helper functions
#' `CHELSA_prepare()` and `CHELSA_project()`.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param download Logical. If `TRUE`, downloads CHELSA files. Default: `FALSE`.
#' @param download_attempts Integer. Maximum download retries. Default: `10`.
#' @param sleep Integer. Seconds to wait between download attempts. Default:
#'   `5`.
#' @param download_n_cores Integer. Number of CPU cores to use for parallel
#'   downloading of CHELSA data. Only valid if download = `TRUE`. Defaults to 4.
#' @param overwrite Logical. If `TRUE`, re-downloads existing files. Default:
#'   `FALSE`.
#' @param compression_level Integer. NetCDF compression level (1 = least, 9 =
#'   most). Default: `5`.
#' @param overwrite_processed Logical. If `TRUE`, overwrites processed files.
#'   Default: `FALSE`.
#' @param other_variables Character. Additional variables to process (e.g.,
#'   `"npp"` for Net Primary Productivity alongside 19 bioclimatic variables
#'   bio1-bio19). Use `""` for bioclimatic only. See [CHELSA_variables] for
#'   details. Default: `"npp"`.
#' @param metadata Tibble. Single-row metadata for input files, prepared by
#'   `CHELSA_prepare()`
#' @section Functions details:
#' - **`CHELSA_process()`**: Main function; optionally downloads CHELSA data,
#'   processes it to the European scale and reference grid, and saves TIFF and
#'   NetCDF outputs for 46 climate scenarios.
#' - **`CHELSA_prepare()`**: Extracts metadata from local URL text files and
#'   manages optional downloads.
#' - **`CHELSA_project()`**: Projects data to the `IASDT` reference grid with
#'   optional transformations.
#' @importFrom rlang !!
#' @note
#' - `CHELSA_prepare()` and `CHELSA_project()` are internal helpers, not for
#' direct use.
#' - Processes 19 bioclimatic variables (bio1–bio19) plus optional variables
#' (e.g., NPP) for 46 scenarios (1 current, 45 future).
#' - Time-intensive; depends on file size and compute resources.

## |------------------------------------------------------------------------| #
# CHELSA_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 1

CHELSA_process <- function(
    env_file = ".env", n_cores = 8L, strategy = "multisession",
    download = FALSE, overwrite = FALSE, download_attempts = 10L, sleep = 5L,
    other_variables = "npp", download_n_cores = 4, compression_level = 5,
    overwrite_processed = FALSE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  ecokit::cat_time("Checking input arguments")

  LogicArgs <- c("download", "overwrite", "overwrite_processed")
  ecokit::check_args(args_to_check = LogicArgs, args_type = "logical")

  NumericArgs <- c("compression_level", "sleep", "download_attempts")
  ecokit::check_args(args_to_check = NumericArgs, args_type = "numeric")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- download_n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)
  download_n_cores <- .validate_n_cores(download_n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_CHELSA_In <- Path_CHELSA_Out <- Path_Out_NC <- TimePeriod <-
    ClimateModel <- ClimateScenario <- Path_Out_tif <- Name <- FilePath <-
    File_List <- Processed_Maps <- Path_Down <-
    InputOkay <- AllOkay <- Process <- Failed <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_CHELSA_In", "DP_R_chelsa_raw", FALSE, FALSE,
    "Path_CHELSA_Out", "DP_R_chelsa_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # Ensure that the output path exists
  fs::dir_create(
    c(
      Path_CHELSA_In, Path_CHELSA_Out,
      fs::path(Path_CHELSA_Out, c("Tif", "NC", "Processed"))))

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "ecokit", "fs", "terra", "stringr", "ncdf4", "magrittr", "dplyr",
      "IASDT.R", "tibble", "purrr", "sf"),
    strategy = strategy)

  # # ..................................................................... ###

  # Prepare CHELSA metadata / download CHELSA data -----

  ecokit::cat_time("Prepare CHELSA metadata or download CHELSA data")
  TimePrepare <- lubridate::now(tzone = "CET")

  # 19 Bioclimatic variables (+ other_variables, if not empty string) * 46 CC
  # options
  CHELSA_Data <- IASDT.R::CHELSA_prepare(
    env_file = env_file, download = download, n_cores = download_n_cores,
    strategy = strategy, overwrite = overwrite,
    other_variables = other_variables, download_attempts = download_attempts,
    sleep = sleep)

  if (nrow(CHELSA_Data) == 0) {
    ecokit::stop_ctx("No CHELSA data was found", include_backtrace = TRUE)
  }

  ecokit::cat_diff(
    init_time = TimePrepare, level = 1L,
    prefix = "Prepare CHELSA metadata or download CHELSA data was finished in ")

  # # ..................................................................... ###

  # Check input CHELSA files -----

  if (isFALSE(download)) {

    # If download = TRUE, there is no need to re-check the files as these files
    # were already checked while downloading the files

    if (n_cores == 1) {
      ecokit::cat_time("Check input CHELSA files sequentially")
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::cat_time("Check input CHELSA files in parallel")
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    CHELSA_Data_Checked <- future.apply::future_lapply(
      X = seq_len(nrow(CHELSA_Data)),
      FUN = function(x) {
        tibble::tibble(Path_Down = CHELSA_Data$Path_Down[x]) %>%
          dplyr::mutate(
            InputOkay = ecokit::check_tiff(Path_Down, warning = FALSE))
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export, future.globals = "CHELSA_Data") %>%
      dplyr::bind_rows() %>%
      dplyr::filter(InputOkay == FALSE) # nolint

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("sequential", gc = TRUE)
    }

    if (nrow(CHELSA_Data_Checked) > 0) {
      readr::write_lines(
        x = dplyr::pull(CHELSA_Data_Checked, "InputOkay"),
        file = fs::path(Path_CHELSA_Out, "ProblematicTiffs.txt"))

      ecokit::stop_ctx(
        paste0(
          "Not all input tiff files are available and valid. ",
          "Check `ProblematicTiffs.txt`"),
        include_backtrace = TRUE)
    }

    # CHELSA files that will not be processed
    Diff <- setdiff(
      list.files(Path_CHELSA_In, pattern = "CHELSA.+.tif$", full.names = TRUE),
      CHELSA_Data$Path_Down)

    if (length(Diff) > 0) {
      message(
        " >> Only Bioclimatic variables and variables identified in ",
        "`other_variables`, if any, will be processed (",
        nrow(CHELSA_Data), " files)\n >> ", length(Diff),
        " files will not be processed.\n",
        " >> See `NotProcessed.txt` for the list of files")

      readr::write_lines(
        x = Diff, file = fs::path(Path_CHELSA_Out, "NotProcessed.txt"))
    }

    rm(AllOkay, Diff, envir = environment())
  }

  # # ..................................................................... ###

  # Processing CHELSA data -----

  TimeProcess <- lubridate::now(tzone = "CET")

  if (n_cores == 1) {
    ecokit::cat_time("Processing CHELSA maps sequentially")
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::cat_time("Processing CHELSA maps in parallel")
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  if (overwrite_processed) {

    ecokit::cat_time(
      "overwrite_processed = TRUE; all files will be processed", level = 1L)

    CHELSA2Process <- dplyr::select(
      .data = CHELSA_Data, Path_Down, Path_Out_NC, Path_Out_tif)

  } else {

    # Exclude previously processed files (after checking)
    ecokit::cat_time(
      "Exclude previously processed files (after checking)", level = 1L)

    check_processed <- future.apply::future_lapply(
      X = seq_len(nrow(CHELSA_Data)),
      FUN = function(x) {

        file_nc <- CHELSA_Data$Path_Out_NC[[x]]
        file_tif <- CHELSA_Data$Path_Out_tif[[x]]
        NC_exists <- file.exists(file_nc)
        tif_exists <- file.exists(file_tif)

        if (!tif_exists || !NC_exists) {
          if (NC_exists) fs::file_delete(file_nc)
          if (tif_exists) fs::file_delete(file_tif)
          return(TRUE)
        }

        NC_Okay <- suppressWarnings(
          ecokit::check_tiff(file_nc, warning = FALSE))
        Tif_Okay <- suppressWarnings(
          ecokit::check_tiff(file_tif, warning = FALSE))
        need_processing <- isFALSE(NC_Okay && Tif_Okay)

        if (need_processing) {
          fs::file_delete(c(file_nc, file_tif))
        }

        return(need_processing)
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export, future.globals = "CHELSA_Data")

    CHELSA2Process <- CHELSA_Data %>%
      dplyr::mutate(Process = unlist(check_processed)) %>%
      dplyr::filter(Process) %>%
      dplyr::select(-"Process")
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||


  # Processing CHELSA files
  ecokit::cat_time("Processing CHELSA files", level = 1L)

  if (nrow(CHELSA2Process) > 0) {

    ecokit::cat_time(
      paste0(
        nrow(CHELSA2Process), " of ", nrow(CHELSA_Data),
        " files need to be processed"),
      level = 2L)

    # process CHELSA files and return info if processing failed
    Failed2process <- future.apply::future_lapply(
      X = seq_len(nrow(CHELSA2Process)),
      FUN = function(x) {

        Path_Down <- CHELSA2Process$Path_Down[[x]]
        file_tif <- CHELSA2Process$Path_Out_tif[[x]]
        file_nc <- CHELSA2Process$Path_Out_NC[[x]]

        Tiffs_okay <- all(
          ecokit::check_tiff(file_tif, warning = FALSE),
          ecokit::check_tiff(file_nc, warning = FALSE)) %>%
          suppressWarnings()

        if (Tiffs_okay) {
          return(FALSE)
        }

        if (file.exists(file_nc)) fs::file_delete(file_nc)
        if (file.exists(file_tif)) fs::file_delete(file_tif)

        # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
        # official parameters (overriding the ones from GeoTIFF keys)
        # see: https://stackoverflow.com/questions/78007307
        terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

        Attempt <- 0
        repeat {
          Attempt <- Attempt + 1
          Try <- try(
            expr = {
              IASDT.R::CHELSA_project(
                metadata = dplyr::slice(CHELSA2Process, x), env_file = env_file,
                compression_level = compression_level) %>%
                # suppress known warning
                # https://github.com/rspatial/terra/issues/1212
                # https://github.com/rspatial/terra/issues/1832
                # https://stackoverflow.com/questions/78098166
                suppressWarnings()
              "Okay"
            },
            silent = TRUE)

          Sys.sleep(2)

          Tiffs_okay <- all(
            ecokit::check_tiff(file_tif, warning = FALSE),
            ecokit::check_tiff(file_nc, warning = FALSE)) %>%
            suppressWarnings()

          if ((!inherits(Try, "try-error") && Tiffs_okay) || Attempt >= 5) {
            break
          }

          if (inherits(Try, "try-error") &&
              isFALSE(
                ecokit::check_tiff(Path_Down, warning = FALSE))) {
            # re-download the file if it fails to be processed and downloaded
            # file is not valid
            system(
              command = CHELSA2Process$DownCommand[[x]],
              ignore.stdout = TRUE, ignore.stderr = TRUE)
          }
          if (file.exists(file_nc)) fs::file_delete(file_nc)
          if (file.exists(file_tif)) fs::file_delete(file_tif)
        }

        if (inherits(Try, "try-error")) {
          return(TRUE)
        } else if (Tiffs_okay) {
          return(FALSE)
        } else {
          return(TRUE)
        }
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c("env_file", "compression_level", "CHELSA2Process"))

    CHELSA2Process <- CHELSA2Process %>%
      dplyr::mutate(Failed = unlist(Failed2process)) %>%
      dplyr::filter(Failed)

    if (nrow(CHELSA2Process) > 0) {
      readr::write_lines(
        x = CHELSA2Process$Path_Down,
        file = fs::path(Path_CHELSA_Out, "FailedProcessing.txt"))

      ecokit::stop_ctx(
        paste0(
          "\n >> ", nrow(CHELSA2Process), " files failed to process.\n",
          " >> Check `FailedProcessing.txt` for more details"),
        CHELSA2Process = CHELSA2Process, include_backtrace = TRUE)
    }

    ecokit::cat_time("All tiff files were processed", level = 1L)

  } else {

    ecokit::cat_time("All tiff files were already processed.", level = 1L)

  }

  rm(CHELSA2Process, envir = environment())

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  ecokit::cat_diff(
    init_time = TimeProcess,
    prefix = "Processing CHELSA data was finished in ", level = 1L)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Grouping CHELSA data by time, climate model/scenario ----

  if (n_cores == 1) {
    ecokit::cat_time(
      "Grouping CHELSA data by time and climate model+scenario sequentially")
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::cat_time(
      "Grouping CHELSA data by time and climate model+scenario in parallel")
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # String to be matched to extract variable names
  SelectedVars <- c("bio", other_variables) %>%
    # Only keep non-empty strings. If `other_variables` = "", only bioclimatic
    # variables will be processed.
    stringr::str_subset(".+") %>%
    # If other_variables = "", this will return "(?i)(bio)\\d*_"
    # If other_variables = "npp", this will return "(?i)(bio|npp)\\d*_"
    # "(?i)" represents case-insensitive matching
    # \\d+_ means one or more digits followed by an underscore
    paste(collapse = "|") %>%
    paste0("(?i)(", ., ")\\d*_")

  CHELSA_Processed <- CHELSA_Data %>%
    dplyr::select(TimePeriod, ClimateModel, ClimateScenario, Path_Out_tif) %>%
    dplyr::slice(gtools::mixedorder(Path_Out_tif)) %>%
    dplyr::summarise(
      File_List = list(Path_Out_tif),
      .by = c(TimePeriod, ClimateModel, ClimateScenario)) %>%

    dplyr::mutate(
      Name = paste0(
        "R_", TimePeriod, "_", ClimateModel, "_", ClimateScenario),
      Name = stringr::str_replace(
        Name, "1981-2010_Current_Current", "Current"),
      Name = stringr::str_replace_all(Name, "-", "_"),

      FilePath = fs::path(
        Path_CHELSA_Out, "Processed", paste0(Name, ".RData")),

      Processed_Maps = furrr::future_pmap(
        .l = list(File_List, FilePath, Name),
        .f = function(File_List, FilePath, Name) {

          MapNames <- basename(File_List) %>%
            # Extract variable names
            stringr::str_extract(SelectedVars) %>%
            # remove the last underscore
            stringr::str_remove_all("_$")

          Map <- terra::rast(File_List) %>%
            stats::setNames(MapNames) %>%
            terra::subset(gtools::mixedsort(MapNames)) %>%
            ecokit::set_raster_crs(crs = "epsg:3035") %>%
            terra::toMemory() %>%
            terra::wrap()

          # save to disk
          ecokit::save_as(
            object = Map, object_name = Name, out_path = FilePath)

          return(Map)

        },
        .options = furrr::furrr_options(
          seed = TRUE, packages = pkg_to_export,
          globals = c("CHELSA_Data", "SelectedVars"))))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  save(
    CHELSA_Processed,
    file = fs::path(Path_CHELSA_Out, "CHELSA_Processed.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  CHELSA_Processed_DT <- dplyr::select(CHELSA_Processed, -Processed_Maps)

  save(
    CHELSA_Processed_DT,
    file = fs::path(Path_CHELSA_Out, "CHELSA_Processed_DT.RData"))

  readr::write_csv(
    x = CHELSA_Processed_DT,
    file = fs::path(Path_CHELSA_Out, "CHELSA_Processed_DT.csv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nProcessing CHELSA data took ")

  return(invisible(NULL))
}
