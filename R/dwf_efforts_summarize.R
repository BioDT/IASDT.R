## |------------------------------------------------------------------------| #
# efforts_summarize ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 4
#' @export

efforts_summarize <- function(
    env_file = ".env", n_cores = 6L, strategy = "multisession",
    chunk_size = 100000L, delete_chunks = TRUE) {

  # # ..................................................................... ###

  .start_time_process <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  speciesKey <- n_obs <- path_interim <- taxa_stand <- path_grid <-
    path_efforts <- NULL

  # # ..................................................................... ###

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_efforts", "DP_R_efforts_processed", TRUE, FALSE,
    "path_interim", "DP_R_efforts_interim", TRUE, FALSE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "taxa_stand", "DP_R_taxa_stand", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  path_efforts_cleaned <- fs::path(path_interim, "cleaned_data")

  ## NAPS species list ----
  naps_list <- readRDS(taxa_stand) %>%
    dplyr::pull("speciesKey") %>%
    unique()

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "ecokit", "stringr", "fs", "sf", "readr", "dplyr", "terra",
      "purrr", "tibble", "R.utils", "IASDT.R", "magrittr"),
    strategy = strategy)

  # # ..................................................................... ###

  path_grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  path_grid_sf <- fs::path(path_grid, "grid_10_land_sf.RData")

  if (!file.exists(path_grid_r)) {
    ecokit::stop_ctx(
      "Reference grid was not found", path_grid_r = path_grid_r,
      include_backtrace = TRUE)
  }

  if (!file.exists(path_grid_sf)) {
    ecokit::stop_ctx(
      "Reference grid file was not found", path_grid_sf = path_grid_sf,
      include_backtrace = TRUE)
  }

  grid_sf <- ecokit::load_as(path_grid_sf)

  # # ..................................................................... ###

  # Prepare working in parallel -----

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Processing data from zipped archives -----
  ecokit::cat_time("Processing data from zipped archives", level = 1L)

  if (delete_chunks) {
    ecokit::cat_time(
      "Chunk files will be deleted after finishing processing", level = 2L)
  }

  # Earlier attempts with `furrr::future_map()` failed

  path_efforts_request <- fs::path(path_efforts, "efforts_all_requests.RData")

  if (!file.exists(path_efforts_request)) {
    ecokit::stop_ctx(
      "The path for the `efforts_all_requests` data does not exist",
      path_efforts_request = path_efforts_request, include_backtrace = TRUE)
  }

  efforts_all_requests <- ecokit::load_as(path_efforts_request)

  data_paths <- future.apply::future_lapply(
    X = seq_len(nrow(efforts_all_requests)),
    FUN = function(id) {

      download_path <- efforts_all_requests$download_path[id]
      total_records <- efforts_all_requests$total_records[id]
      class <- efforts_all_requests$class[id]
      order <- efforts_all_requests$order[id]
      class_order <- paste0(class, "_", order)

      # Output path to save the data
      path_data <- fs::path(path_efforts_cleaned, paste0(class_order, ".RData"))

      # Should path_data be returned as the path of the RData file containing
      # the data or NA if there are no records in the current order or no
      # records left after processing
      returned_no_data <- (total_records == 0)

      # Check if the RData file for the current order exists and valid.
      file_okay <- ecokit::check_data(path_data, warning = FALSE)

      # Process current order data if the output file is not okay and the order
      # have observations
      if (isFALSE(file_okay) && total_records > 0) {

        if (file.exists(path_data)) {
          fs::file_delete(path_data)
        }

        # Check if previous chunk files for the current order exist and contain
        # the same total number of observations. If this is true, do not split
        # the data and use the chunk files directly; otherwise, split the data
        # first into small chunks. This helps to continue working on the same
        # data should previous function try failed.

        # List of chunks
        chunks <- list.files(
          path = path_interim, full.names = TRUE,
          pattern = stringr::str_remove(basename(download_path), ".zip"))

        # If there are chunk files on disk, count their total number of
        # observations
        if (length(chunks) > 0) {
          # Total number of lines in all chunk files
          n_lines <- sum(purrr::map_int(chunks, R.utils::countLines))

          # if there are less than the total number of records, delete the chunk
          # files and recreate them
          if (n_lines != total_records) {
            purrr::walk(chunks, file.remove)
            split_chunks <- TRUE
            rm(chunks, envir = environment())
          } else {
            split_chunks <- FALSE
          }
        } else {
          # If there is no chunk files available, split into chunks
          split_chunks <- TRUE
        }


        # Split data into chunks
        if (split_chunks) {
          chunks <- IASDT.R::efforts_split(
            path_zip = download_path, env_file = ".env",
            chunk_size = chunk_size)
        }

        # Process chunk files
        # nolint start
        accepted_ranks <- c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")
        column_names <- c(
          "taxonRank", "Latitude", "Longitude", "uncertain_km", "speciesKey")
        column_types <- readr::cols(
          uncertain_km = readr::col_double(),
          Longitude = readr::col_double(),
          Latitude = readr::col_double(),
          speciesKey = readr::col_integer(),
          taxonRank = readr::col_character(),
          .default = readr::col_double())
        # nolint end

        chunk_data <- purrr::map_dfr(
          .x = chunks,
          .f = ~ {
            readr::read_tsv(
              file = .x, col_names = column_names, progress = FALSE,
              show_col_types = FALSE, col_types = column_types) %>%
              dplyr::mutate(uncertain_km = uncertain_km / 1000) %>%
              dplyr::filter(
                !is.na(Latitude), !is.na(Longitude),
                nzchar(speciesKey), taxonRank %in% accepted_ranks,
                uncertain_km <= 100 | is.na(uncertain_km)) %>%
              sf::st_as_sf(
                coords = c("Longitude", "Latitude"),
                crs = 4326, remove = FALSE) %>%
              sf::st_transform(3035) %>%
              sf::st_join(grid_sf) %>%
              dplyr::filter(magrittr::not(is.na(CellCode)))
          }
        )

        # if there are observations after the filtering, save the data to disk
        # and return the saved path, otherwise return no path
        if (nrow(chunk_data) > 0) {
          ecokit::save_as(
            object = chunk_data, object_name = class_order,
            out_path = path_data)
        } else {
          returned_no_data <- TRUE
        }

        rm(chunk_data, envir = environment())

        # delete chunk files for the current order
        if (delete_chunks && (total_records > 0)) {
          purrr::walk(chunks, file.remove)
        }
      }

      invisible(gc())

      tibble::tibble(
        class_order = class_order,
        path_data = dplyr::if_else(returned_no_data, NA_character_, path_data),
        class = class, order = order)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c(
      "path_interim", "efforts_all_requests", "path_efforts_cleaned",
      "grid_sf", "chunk_size", "delete_chunks")) %>%
    dplyr::bind_rows()

  # # ++++++++++++++++++++++++++++++ ###

  # only selected columns from `efforts_all_requests`
  request_cols <- c(
    "class", "order", "request", "download_link",
    "total_records", "download_path")

  # join data with requests summary
  efforts_summary <- efforts_all_requests %>%
    dplyr::select(tidyselect::all_of(request_cols)) %>%
    dplyr::left_join(data_paths, by = c("class", "order"))

  rm(data_paths, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Prepare summary maps per order ----
  ecokit::cat_time("Prepare summary maps per order", level = 1L)

  summary_maps <- future.apply::future_lapply(
    X = seq_len(nrow(efforts_summary)),
    FUN = function(id) {
      path_data <- efforts_summary$path_data[id]
      class_order <- efforts_summary$class_order[id]

      if (is.na(path_data)) {
        # If there is no data for the current order
        n_obs <- n_obs_native <- 0L
      } else {
        # Load data on current order
        summary_data <- ecokit::load_as(path_data)
        # Number of observation in the cleaned data
        n_obs <- nrow(summary_data)

        # Only native species (exclude NAPS list)
        summary_data_native <- dplyr::filter(
          summary_data, !(speciesKey %in% naps_list))
        # Number of data for native species
        n_obs_native <- nrow(summary_data_native)
      }
      grid_r <- ecokit::load_as(path_grid_r, unwrap_r = TRUE)

      # # ++++++++++++++++++++++++++++++++++ ###

      # All species ----

      if (n_obs == 0) {
        # Create dummy maps for the number of species and records
        n_obs_r <- terra::classify(grid_r, cbind(1, 0)) %>%
          stats::setNames(paste0("n_obs_", class_order)) %>%
          terra::wrap()

        n_sp_r <- terra::classify(grid_r, cbind(1, 0)) %>%
          stats::setNames(paste0("n_sp_", class_order)) %>%
          terra::wrap()
      } else {
        # Number of observations
        n_obs_r <- efforts_summarize_maps(
          input_data = summary_data, n_sp = FALSE, count_name = "n_obs",
          class_order = class_order, grid_sf = grid_sf, grid_r = grid_r)

        # Number of species
        n_sp_r <- efforts_summarize_maps(
          input_data = summary_data, n_sp = TRUE, count_name = "n_sp",
          class_order = class_order, grid_sf = grid_sf, grid_r = grid_r)

        rm(summary_data, envir = environment())
      }

      # # ++++++++++++++++++++++++++++++++++ ###

      # Only native species ----

      if (n_obs_native == 0) {
        # Create dummy maps for the number of species and records

        n_obs_native_r <- terra::classify(grid_r, cbind(1, 0)) %>%
          stats::setNames(paste0("n_obs_native_", class_order)) %>%
          terra::wrap()

        n_sp_native_r <- terra::classify(grid_r, cbind(1, 0)) %>%
          stats::setNames(paste0("n_sp_native_", class_order)) %>%
          terra::wrap()
      } else {
        # Number of observations of native species
        n_obs_native_r <- efforts_summarize_maps(
          input_data = summary_data_native, n_sp = FALSE,
          count_name = "n_obs_native", class_order = class_order,
          grid_sf = grid_sf, grid_r = grid_r)

        # Number of native species
        n_sp_native_r <- efforts_summarize_maps(
          input_data = summary_data_native, n_sp = TRUE,
          count_name = "n_sp_native", class_order = class_order,
          grid_sf = grid_sf, grid_r = grid_r)

        rm(summary_data_native, envir = environment())
      }

      # # ++++++++++++++++++++++++++++++++++ ###

      tibble::tibble(
        class_order = class_order,
        n_obs = n_obs, n_obs_r = list(n_obs_r), n_sp_r = list(n_sp_r),
        n_obs_native = n_obs_native, n_obs_native_r = list(n_obs_native_r),
        n_sp_native_r = list(n_sp_native_r))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("path_grid_r", "grid_sf", "naps_list")) %>%
    dplyr::bind_rows()


  # join data with requests summary
  efforts_summary <- dplyr::left_join(
    efforts_summary, summary_maps, by = "class_order")

  rm(summary_maps, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Stopping cluster ------
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  invisible(gc())

  # # ..................................................................... ###

  # Save summary results: `efforts_summary` ----
  ecokit::cat_time("Save summary results: `efforts_summary`", level = 1L)
  save(
    efforts_summary,
    file = fs::path(path_efforts, "efforts_summary.RData"))

  # # ..................................................................... ###

  # Prepare summary maps - all sampling efforts ----
  ecokit::cat_time("Prepare summary maps - all sampling efforts", level = 1L)

  calc_n_obs_n_sp <- function(map_list, map_name) {
    purrr::map(.x = unlist(map_list), .f = terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE) %>%
      ecokit::set_raster_crs(crs = "epsg:3035") %>%
      terra::toMemory() %>%
      stats::setNames(map_name)
  }

  # # ..................................................................... ###

  # Exclude orders with no data
  efforts_summary_r <- dplyr::filter(efforts_summary, n_obs > 0)

  efforts_summary_r <- list(
    calc_n_obs_n_sp(efforts_summary_r$n_obs_r, "n_obs"),
    calc_n_obs_n_sp(efforts_summary_r$n_obs_native_r, "n_obs_native"),
    calc_n_obs_n_sp(efforts_summary_r$n_sp_r, "n_sp"),
    calc_n_obs_n_sp(efforts_summary_r$n_sp_native_r, "n_sp_native")) %>%
    terra::rast() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::toMemory()

  ## Save summary maps - `RData` ----
  ecokit::cat_time("Save summary maps", level = 1L)

  ecokit::cat_time("`RData`", level = 2L)
  ecokit::save_as(
    object = terra::wrap(efforts_summary_r), object_name = "efforts_summary_r",
    out_path = fs::path(path_efforts, "efforts_summary_r.RData"))

  ## Save summary maps - `tif` ----
  ecokit::cat_time("`tif`", level = 2L)
  terra::writeRaster(
    efforts_summary_r,
    overwrite = TRUE,
    filename = fs::path(
      path_efforts, paste0("efforts_gbif_", names(efforts_summary_r), ".tif")))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time_process,
    prefix = "Processing Efforts data took ", level = 1L)

  # # ..................................................................... ###

  return(invisible(NULL))
}

## |------------------------------------------------------------------------| #
# efforts_summarize_maps ----
## |------------------------------------------------------------------------| #

#' Summarize maps for efforts data
#'
#' This function processes spatial data (as an `sf` object), summarises it based
#' on the number of observations or distinct species, and generates a raster
#' layer.
#' @param input_data An `sf` object containing spatial data, with a column named
#'   `CellCode`.
#' @param n_sp Logical. Whether to generate distinct species counts (`TRUE`) or
#'   total observation counts (`FALSE`).
#' @param count_name Character. Name of the count field and the prefix for the
#'  final raster layer's name.
#' @param class_order Character. The class and order combination (separated by
#'   an underscore) represented in the `input_data`.
#' @param grid_sf,grid_r Reference grid in the form of simple feature and
#'   raster.
#' @return A processed `terra` raster object representing the summarised data.
#' @note This function is not intended to be used directly by the user, but
#'   only used inside the [efforts_process] and [efforts_summarize] functions.
#' @author Ahmed El-Gabbas
#' @name efforts_summarize_maps
#' @keywords internal
#' @noRd

efforts_summarize_maps <- function(
    input_data, n_sp, count_name, class_order, grid_sf, grid_r) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CellCode <- speciesKey <- NULL

  # # ..................................................................... ###

  # Validate if input_data is an sf object
  if (!inherits(input_data, "sf")) {
    ecokit::stop_ctx(
      paste0(
        "Input data must be a simple feature (sf) object. ",
        "Provided data is of type: ", paste(class(input_data), collapse = "+")),
      input_data = input_data, class_data = class(input_data),
      include_backtrace = TRUE)
  }

  # Validate if n_sp is logical
  if (!is.logical(n_sp) || length(n_sp) != 1) {
    ecokit::stop_ctx(
      paste0(
        "The parameter `n_sp` must be a single logical value (TRUE or FALSE). ",
        "Provided value is of type: ", paste(class(n_sp), collapse = "+")),
      n_sp = n_sp, class_n_sp = class(n_sp), length_n_sp = length(n_sp),
      include_backtrace = TRUE)
  }

  # Validate the count_name parameter
  if (is.null(count_name)) {
    ecokit::stop_ctx(
      "The parameter `count_name` can not be empty", count_name = count_name,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Drop geometry from input_data
  summary_data <- sf::st_drop_geometry(input_data)

  # Generate distinct species counts if n_sp is TRUE
  if (n_sp) {
    summary_data <- dplyr::distinct(summary_data, CellCode, speciesKey)
  }

  # Count observations or species, join with the grid, and rasterize
  summary_data <- summary_data %>%
    dplyr::count(CellCode, name = count_name) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(grid_r, field = count_name) %>%
    terra::classify(cbind(NA, 0)) %>%
    terra::mask(grid_r) %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::toMemory() %>%
    stats::setNames(paste0(count_name, "_", class_order)) %>%
    terra::wrap()

  summary_data
}
