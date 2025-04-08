#' Process CHELSA Climate Data for the `IAS-pDT`
#'
#' Downloads, processes, and projects [CHELSA](https://chelsa-climate.org/)
#' climate data at the European scale for the Invasive Alien Species prototype
#' Digital Twin (`IAS-pDT`). Supports multiple climate scenarios, outputting
#' data in TIFF and NetCDF formats. Orchestrated by `CHELSA_process()`, with
#' helper functions `CHELSA_prepare()` and `CHELSA_project()`.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
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
#' - **`CHELSA_project()`**: Projects data to the IAS-pDT reference grid with
#'   optional transformations.
#'
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
    env_file = ".env", n_cores = 8L, download = FALSE, overwrite = FALSE,
    download_attempts = 10L, sleep = 5L, other_variables = "npp",
    download_n_cores = 4, compression_level = 5, overwrite_processed = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  IASDT.R::cat_time("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c("env_file", "other_variables")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = CharArgs, args_type = "character")

  LogicArgs <- c("download", "overwrite", "overwrite_processed")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = LogicArgs, args_type = "logical")

  NumericArgs <- c(
    "download_n_cores", "compression_level",
    "sleep", "download_attempts", "n_cores")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = NumericArgs, args_type = "numeric")

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs, envir = environment())

  if (n_cores < 1 || download_n_cores < 1) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_CHELSA_In <- Path_CHELSA_Out <- Path_Out_NC <- TimePeriod <-
    ClimateModel <- ClimateScenario <- Path_Out_tif <- Name <- FilePath <-
    File_List <- Processed_Maps <- Path_Down <-
    InputOkay <- AllOkay <- Process <- Failed <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  IASDT.R::cat_time("Environment variables")

  if (!file.exists(env_file)) {
    stop(
      "Path to environment variables: ", env_file, " was not found",
      call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_CHELSA_In", "DP_R_CHELSA_raw", FALSE, FALSE,
    "Path_CHELSA_Out", "DP_R_CHELSA_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Ensure that the output path exists
  fs::dir_create(
    c(
      Path_CHELSA_In, Path_CHELSA_Out,
      IASDT.R::path(Path_CHELSA_Out, c("Tif", "NC", "Processed"))))

  # # ..................................................................... ###

  # Prepare CHELSA metadata / download CHELSA data -----

  IASDT.R::cat_time("Prepare CHELSA metadata / download CHELSA data")
  TimePrepare <- lubridate::now(tzone = "CET")

  # 19 Bioclimatic variables (+ other_variables, if not empty string) * 46 CC
  # options
  CHELSA_Data <- IASDT.R::CHELSA_prepare(
    env_file = env_file, download = download, n_cores = download_n_cores,
    overwrite = overwrite, other_variables = other_variables,
    download_attempts = download_attempts, sleep = sleep)

  IASDT.R::cat_diff(
    init_time = TimePrepare, level = 1,
    prefix = "Prepare CHELSA metadata / download CHELSA data was finished in ")

  # # ..................................................................... ###

  # Check input CHELSA files -----

  if (isFALSE(download)) {

    # If download = TRUE, there is no need to re-check the files as these files
    # were already checked while downloading the files

    if (n_cores == 1) {
      IASDT.R::cat_time("Check input CHELSA files sequentially")
      future::plan("future::sequential", gc = TRUE)
    } else {
      IASDT.R::cat_time("Check input CHELSA files on parallel")
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(n_cores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      withr::defer(future::plan("future::sequential", gc = TRUE))
    }

    CHELSA_Data_Checked <- CHELSA_Data %>%
      dplyr::select(Path_Down) %>%
      dplyr::mutate(
        InputOkay = furrr::future_map_lgl(
          .x = Path_Down,
          .f = IASDT.R::check_tiff, warning = FALSE,
          .options = furrr::furrr_options(seed = TRUE, packages = "IASDT.R"))
      ) %>%
      dplyr::filter(isFALSE(InputOkay))

    if (n_cores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

    if (nrow(CHELSA_Data_Checked) > 0) {
      readr::write_lines(
        x = dplyr::pull(CHELSA_Data_Checked, "InputOkay"),
        file = IASDT.R::path(Path_CHELSA_Out, "ProblematicTiffs.txt"))

      stop(
        "Not all input tiff files are available and valid. ",
        "Check `ProblematicTiffs.txt`", call. = FALSE
      )
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
        x = Diff, file = IASDT.R::path(Path_CHELSA_Out, "NotProcessed.txt"))
    }

    rm(AllOkay, Diff, envir = environment())
  }

  # # ..................................................................... ###

  # Processing CHELSA data -----

  TimeProcess <- lubridate::now(tzone = "CET")

  if (n_cores == 1) {
    IASDT.R::cat_time("Processing CHELSA maps sequentially")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    IASDT.R::cat_time("Processing CHELSA maps on parallel")
    c1 <- snow::makeSOCKcluster(n_cores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }


  if (overwrite_processed) {

    IASDT.R::cat_time(
      "overwrite_processed = TRUE; all files will be processed", level = 1)

    CHELSA2Process <- dplyr::select(
      .data = CHELSA_Data, Path_Down, Path_Out_NC, Path_Out_tif)

  } else {

    # Exclude previously processed files (after checking)
    IASDT.R::cat_time(
      "Exclude previously processed files (after checking)", level = 1)

    CHELSA2Process <- CHELSA_Data %>%
      dplyr::select(Path_Down, Path_Out_NC, Path_Out_tif) %>%
      dplyr::mutate(
        Process = furrr::future_map2_lgl(
          .x = Path_Out_NC, .y = Path_Out_tif,
          .f = ~ {
            NC_Okay <- IASDT.R::check_tiff(.x, warning = FALSE)
            Tif_Okay <- IASDT.R::check_tiff(.y, warning = FALSE)
            return(isFALSE(NC_Okay && Tif_Okay))
          },
          .options = furrr::furrr_options(seed = TRUE, packages = "IASDT.R"))
      ) %>%
      dplyr::filter(Process) %>%
      dplyr::select(-"Process")
  }

  # Processing CHELSA files
  IASDT.R::cat_time("Processing CHELSA files", level = 1)

  if (nrow(CHELSA2Process) > 0) {
    CHELSA2Process <- CHELSA2Process %>%
      dplyr::mutate(
        Failed = furrr::future_pmap_lgl(
          .l = list(Path_Out_NC, Path_Out_tif, Path_Down),
          .f = function(Path_Out_NC, Path_Out_tif, Down = Path_Down) {
            FileMetadata <- dplyr::filter(CHELSA_Data, Path_Down == Down)

            # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
            # official parameters (overriding the ones from GeoTIFF keys)
            # see: https://stackoverflow.com/questions/78007307
            terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

            Attempt <- 0
            repeat {
              Attempt <- Attempt + 1
              Try <- try(
                IASDT.R::CHELSA_project(
                  metadata = FileMetadata, env_file = env_file,
                  compression_level = compression_level),
                silent = TRUE)

              if (!inherits(Try, "try-error") || Attempt >= 5) {
                break
              }

              if (inherits(Try, "try-error")) {
                # re-download the file if it fails to be processed
                system(
                  command = FileMetadata$DownCommand,
                  ignore.stdout = TRUE, ignore.stderr = TRUE)
              }
            }

            invisible(gc())

            Tiffs_okay <- all(
              IASDT.R::check_tiff(Path_Out_tif, warning = FALSE),
              IASDT.R::check_tiff(Path_Out_NC, warning = FALSE))

            if (inherits(Try, "try-error")) {
              return(TRUE)
            } else if (Tiffs_okay) {
              return(FALSE)
            } else {
              return(TRUE)
            }
          },
          .options = furrr::furrr_options(
            seed = TRUE,
            packages = c("dplyr", "terra", "IASDT.R", "tibble", "ncdf4"),
            globals = c("CHELSA_Data", "env_file", "compression_level"))
        )
      ) %>%
      dplyr::filter(Failed)

    if (nrow(CHELSA2Process) > 0) {
      readr::write_lines(
        x = CHELSA2Process$Path_Down,
        file = IASDT.R::path(Path_CHELSA_Out, "FailedProcessing.txt"))

      stop(
        "\n >> ", nrow(CHELSA2Process), " files failed to process.\n",
        " >> Check `FailedProcessing.txt` for more details", call. = FALSE)
    }

    IASDT.R::cat_time("All tiff files were processed", level = 1)

  } else {

    IASDT.R::cat_time("All tiff files were already processed.", level = 1)

  }

  rm(CHELSA2Process, envir = environment())

  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::cat_diff(
    init_time = TimeProcess,
    prefix = "Processing CHELSA data was finished in ", level = 1)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Group CHELSA data by time, climate model/scenario ----

  if (n_cores == 1) {
    IASDT.R::cat_time(
      "Group CHELSA data by time and climate model/scenario sequentially")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    IASDT.R::cat_time(
      "Group CHELSA data by time and climate model/scenario on parallel")
    c1 <- snow::makeSOCKcluster(n_cores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    withr::defer(future::plan("future::sequential", gc = TRUE))
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

      FilePath = IASDT.R::path(
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
            IASDT.R::set_raster_CRS() %>%
            IASDT.R::set_raster_values() %>%
            terra::wrap()

          # save to disk
          IASDT.R::save_as(
            object = Map, object_name = Name, out_path = FilePath)

          return(Map)

        },
        .options = furrr::furrr_options(
          seed = TRUE, packages = c("terra", "IASDT.R", "stringr"),
          globals = c("CHELSA_Data", "SelectedVars"))))

  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  save(
    CHELSA_Processed,
    file = IASDT.R::path(Path_CHELSA_Out, "CHELSA_Processed.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  CHELSA_Processed_DT <- dplyr::select(CHELSA_Processed, -Processed_Maps)

  save(
    CHELSA_Processed_DT,
    file = IASDT.R::path(Path_CHELSA_Out, "CHELSA_Processed_DT.RData"))

  readr::write_csv(
    x = CHELSA_Processed_DT,
    file = IASDT.R::path(Path_CHELSA_Out, "CHELSA_Processed_DT.csv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nProcessing CHELSA data took ")

  return(invisible(NULL))
}
