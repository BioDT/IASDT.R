#' Process CHELSA Climate Data for the `IASDT`
#'
#' Downloads, processes, and projects [CHELSA](https://chelsa-climate.org/)
#' climate data at the European scale for the Invasive Alien Species Digital
#' Twin (`IASDT`). Supports multiple climate scenarios, outputting data in TIFF
#' and NetCDF formats. Orchestrated by `chelsa_process()`, with helper functions
#' `chelsa_prepare()` and `chelsa_project()`.
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
#'   bio1-bio19). Use `""` for bioclimatic only. See [chelsa_variables] for
#'   details. Default: `"npp"`.
#' @param metadata Tibble. Single-row metadata for input files, prepared by
#'   `chelsa_prepare()`
#' @section Functions details:
#' - **`chelsa_process()`**: Main function; optionally downloads CHELSA data,
#'   processes it to the European scale and reference grid, and saves TIFF and
#'   NetCDF outputs for 46 climate scenarios.
#' - **`chelsa_prepare()`**: Extracts metadata from local URL text files and
#'   manages optional downloads.
#' - **`chelsa_project()`**: Projects data to the `IASDT` reference grid with
#'   optional transformations.
#' @importFrom rlang !!
#' @note
#' - `chelsa_prepare()` and `chelsa_project()` are internal helpers, not for
#' direct use.
#' - Processes 19 bioclimatic variables (bio1–bio19) plus optional variables
#' (e.g., NPP) for 46 scenarios (1 current, 45 future).
#' - Time-intensive; depends on file size and compute resources.

## |------------------------------------------------------------------------| #
# chelsa_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name chelsa_data
#' @rdname chelsa_data
#' @order 1

chelsa_process <- function(
    env_file = ".env", n_cores = 8L, strategy = "multisession",
    download = FALSE, overwrite = FALSE, download_attempts = 10L, sleep = 5L,
    other_variables = "npp", download_n_cores = 4, compression_level = 5,
    overwrite_processed = FALSE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  ecokit::cat_time("Checking input arguments")

  logical_args <- c("download", "overwrite", "overwrite_processed")
  ecokit::check_args(args_to_check = logical_args, args_type = "logical")

  numeric_args <- c("compression_level", "sleep", "download_attempts")
  ecokit::check_args(args_to_check = numeric_args, args_type = "numeric")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- download_n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)
  download_n_cores <- .validate_n_cores(download_n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_chelsa_in <- path_chelsa_out <- path_out_netcdf <- time_period <-
    climate_model <- climate_scenario <- path_out_tif <- climate_name <-
    file_path <- file_list <- processed_maps <- path_download <-
    input_okay <- process <- failed <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_chelsa_in", "DP_R_chelsa_raw", FALSE, FALSE,
    "path_chelsa_out", "DP_R_chelsa_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # Ensure that the output path exists
  fs::dir_create(
    c(
      path_chelsa_in, path_chelsa_out,
      fs::path(
        path_chelsa_out, c("chelsa_tif", "chelsa_netcdf", "projected_rdata"))))

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
  .time_prepare <- lubridate::now(tzone = "CET")

  # 19 Bioclimatic variables (+ other_variables, if not empty string) * 46 CC
  # options
  chelsa_data <- IASDT.R::chelsa_prepare(
    env_file = env_file, download = download, n_cores = download_n_cores,
    strategy = strategy, overwrite = overwrite,
    other_variables = other_variables, download_attempts = download_attempts,
    sleep = sleep)

  if (nrow(chelsa_data) == 0) {
    ecokit::stop_ctx("No CHELSA data was found", include_backtrace = TRUE)
  }

  ecokit::cat_diff(
    init_time = .time_prepare, level = 1L,
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

    chelsa_data_checked <- future.apply::future_lapply(
      X = seq_len(nrow(chelsa_data)),
      FUN = function(x) {
        tibble::tibble(path_download = chelsa_data$path_download[x]) %>%
          dplyr::mutate(
            input_okay = ecokit::check_tiff(path_download, warning = FALSE))
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export, future.globals = "chelsa_data") %>%
      dplyr::bind_rows() %>%
      dplyr::filter(input_okay == FALSE) # nolint

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("sequential", gc = TRUE)
    }

    if (nrow(chelsa_data_checked) > 0) {
      readr::write_lines(
        x = dplyr::pull(chelsa_data_checked, "input_okay"),
        file = fs::path(path_chelsa_out, "problematic_tiffs.txt"))

      ecokit::stop_ctx(
        paste0(
          "Not all input tiff files are available and valid. ",
          "Check `problematic_tiffs.txt`"),
        include_backtrace = TRUE)
    }

    # CHELSA files that will not be processed
    file_diff <- setdiff(
      list.files(path_chelsa_in, pattern = "CHELSA.+.tif$", full.names = TRUE),
      chelsa_data$path_download)

    if (length(file_diff) > 0) {
      message(
        " >> Only Bioclimatic variables and variables identified in ",
        "`other_variables`, if any, will be processed (",
        nrow(chelsa_data), " files)\n >> ", length(file_diff),
        " files will not be processed.\n",
        " >> See `NotProcessed.txt` for the list of files")

      readr::write_lines(
        x = file_diff, file = fs::path(path_chelsa_out, "NotProcessed.txt"))
    }

    rm(file_diff, envir = environment())
  }

  # # ..................................................................... ###

  # Processing CHELSA data -----

  .time_process <- lubridate::now(tzone = "CET")

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

    chelsa_to_process <- dplyr::select(
      chelsa_data, path_download, path_out_netcdf, path_out_tif)

  } else {

    # Exclude previously processed files (after checking)
    ecokit::cat_time(
      "Exclude previously processed files (after checking)", level = 1L)

    check_processed <- future.apply::future_lapply(
      X = seq_len(nrow(chelsa_data)),
      FUN = function(x) {

        file_nc <- chelsa_data$path_out_netcdf[[x]]
        file_tif <- chelsa_data$path_out_tif[[x]]
        netcdf_exists <- file.exists(file_nc)
        tif_exists <- file.exists(file_tif)

        if (!tif_exists || !netcdf_exists) {
          if (netcdf_exists) fs::file_delete(file_nc)
          if (tif_exists) fs::file_delete(file_tif)
          return(TRUE)
        }

        netcdf_okay <- suppressWarnings(
          ecokit::check_tiff(file_nc, warning = FALSE))
        tif_okay <- suppressWarnings(
          ecokit::check_tiff(file_tif, warning = FALSE))
        need_processing <- isFALSE(netcdf_okay && tif_okay)

        if (need_processing) {
          fs::file_delete(c(file_nc, file_tif))
        }

        need_processing
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export, future.globals = "chelsa_data")

    chelsa_to_process <- chelsa_data %>%
      dplyr::mutate(process = unlist(check_processed)) %>%
      dplyr::filter(process) %>%
      dplyr::select(-"process")
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||


  # Processing CHELSA files
  ecokit::cat_time("Processing CHELSA files", level = 1L)

  if (nrow(chelsa_to_process) > 0) {

    ecokit::cat_time(
      paste0(
        nrow(chelsa_to_process), " of ", nrow(chelsa_data),
        " files need to be processed"),
      level = 2L)

    # process CHELSA files and return info if processing failed
    failed_to_process <- future.apply::future_lapply(
      X = seq_len(nrow(chelsa_to_process)),
      FUN = function(x) {

        path_download <- chelsa_to_process$path_download[[x]]
        file_tif <- chelsa_to_process$path_out_tif[[x]]
        file_nc <- chelsa_to_process$path_out_netcdf[[x]]

        tiffs_okay <- all(
          ecokit::check_tiff(file_tif, warning = FALSE),
          ecokit::check_tiff(file_nc, warning = FALSE)) %>%
          suppressWarnings()

        if (tiffs_okay) {
          return(FALSE)
        }

        if (file.exists(file_nc)) fs::file_delete(file_nc)
        if (file.exists(file_tif)) fs::file_delete(file_tif)

        # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
        # official parameters (overriding the ones from GeoTIFF keys)
        # see: https://stackoverflow.com/questions/78007307
        terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

        attempt <- 0
        repeat {
          attempt <- attempt + 1
          try_n <- try(
            expr = {
              IASDT.R::chelsa_project(
                metadata = dplyr::slice(chelsa_to_process, x),
                env_file = env_file, compression_level = compression_level) %>%
                # suppress known warning
                # https://github.com/rspatial/terra/issues/1212
                # https://github.com/rspatial/terra/issues/1832
                # https://stackoverflow.com/questions/78098166
                suppressWarnings()
              "Okay"
            },
            silent = TRUE)

          Sys.sleep(2)

          tiffs_okay <- all(
            ecokit::check_tiff(file_tif, warning = FALSE),
            ecokit::check_tiff(file_nc, warning = FALSE)) %>%
            suppressWarnings()

          if ((!inherits(try_n, "try-error") && tiffs_okay) || attempt >= 5) {
            break
          }

          if (inherits(try_n, "try-error") &&
              isFALSE(
                ecokit::check_tiff(path_download, warning = FALSE))) {
            # re-download the file if it fails to be processed and downloaded
            # file is not valid
            system(
              command = chelsa_to_process$download_command[[x]],
              ignore.stdout = TRUE, ignore.stderr = TRUE)
          }
          if (file.exists(file_nc)) fs::file_delete(file_nc)
          if (file.exists(file_tif)) fs::file_delete(file_tif)
        }

        if (inherits(try_n, "try-error")) {
          TRUE
        } else if (tiffs_okay) {
          FALSE
        } else {
          TRUE
        }
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c("env_file", "compression_level", "chelsa_to_process"))

    chelsa_to_process <- chelsa_to_process %>%
      dplyr::mutate(failed = unlist(failed_to_process)) %>%
      dplyr::filter(failed)

    if (nrow(chelsa_to_process) > 0) {
      readr::write_lines(
        x = chelsa_to_process$path_download,
        file = fs::path(path_chelsa_out, "failed_processing.txt"))

      ecokit::stop_ctx(
        paste0(
          "\n >> ", nrow(chelsa_to_process), " files failed to process.\n",
          " >> Check `FailedProcessing.txt` for more details"),
        chelsa_to_process = chelsa_to_process, include_backtrace = TRUE)
    }

    ecokit::cat_time("All tiff files were processed", level = 1L)

  } else {

    ecokit::cat_time("All tiff files were already processed.", level = 1L)

  }

  rm(chelsa_to_process, envir = environment())

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  ecokit::cat_diff(
    init_time = .time_process,
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
  selected_vars <- c("bio", other_variables) %>%
    # Only keep non-empty strings. If `other_variables` = "", only bioclimatic
    # variables will be processed.
    stringr::str_subset(".+") %>%
    # If other_variables = "", this will return "(?i)(bio)\\d*_"
    # If other_variables = "npp", this will return "(?i)(bio|npp)\\d*_"
    # "(?i)" represents case-insensitive matching
    # \\d+_ means one or more digits followed by an underscore
    paste(collapse = "|") %>%
    paste0("(?i)(", ., ")\\d*_")

  chelsa_processed <- chelsa_data %>%
    dplyr::select(
      time_period, climate_model, climate_scenario, path_out_tif) %>%
    dplyr::slice(gtools::mixedorder(path_out_tif)) %>%
    dplyr::summarise(
      file_list = list(path_out_tif),
      .by = c(time_period, climate_model, climate_scenario)) %>%

    dplyr::mutate(
      climate_name = paste0(
        "R_", time_period, "_", climate_model, "_", climate_scenario),
      climate_name = stringr::str_replace(
        climate_name, "1981-2010_current_current", "current"),
      climate_name = stringr::str_replace_all(climate_name, "-", "_"),

      file_path = fs::path(
        path_chelsa_out, "projected_rdata", paste0(climate_name, ".RData")),

      processed_maps = furrr::future_pmap(
        .l = list(file_list, file_path, climate_name),
        .f = function(file_list, file_path, climate_name) {

          map_names <- basename(file_list) %>%
            # Extract variable names
            stringr::str_extract(selected_vars) %>%
            # remove the last underscore
            stringr::str_remove_all("_$")

          out_map <- terra::rast(file_list) %>%
            stats::setNames(map_names) %>%
            terra::subset(gtools::mixedsort(map_names)) %>%
            ecokit::set_raster_crs(crs = "epsg:3035") %>%
            terra::toMemory() %>%
            terra::wrap()

          # save to disk
          ecokit::save_as(
            object = out_map, object_name = climate_name, out_path = file_path)

          out_map
        },
        .options = furrr::furrr_options(
          seed = TRUE, packages = pkg_to_export,
          globals = c("chelsa_data", "selected_vars"))))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  save(
    chelsa_processed,
    file = fs::path(path_chelsa_out, "chelsa_processed.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  chelsa_processed_data <- dplyr::select(chelsa_processed, -processed_maps)

  save(
    chelsa_processed_data,
    file = fs::path(path_chelsa_out, "chelsa_processed_data.RData"))

  readr::write_csv(
    x = chelsa_processed_data,
    file = fs::path(path_chelsa_out, "chelsa_processed_data.csv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nProcessing CHELSA data took ")

  return(invisible(NULL))
}
