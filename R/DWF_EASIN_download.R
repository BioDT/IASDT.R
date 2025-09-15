## |------------------------------------------------------------------------| #
# EASIN_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 3

EASIN_download <- function(
    species_key, timeout = 200, verbose = FALSE, env_file = ".env",
    n_search = 1000, n_attempts = 10, sleep_time = 5,
    delete_chunks = TRUE, return_data = FALSE) {

  # # ..................................................................... ###

  if (is.null(species_key)) {
    ecokit::stop_ctx(
      "species_key cannot be NULL", species_key = species_key,
      include_backtrace = TRUE)
  }

  path_EASIN <- NULL

  # # ..................................................................... ###

  # Checking arguments ----
  if (verbose) {
    ecokit::cat_time("Checking arguments")
  }

  ecokit::check_args(args_to_check = "species_key", args_type = "character")
  ecokit::check_args(
    args_to_check = c("return_data", "verbose", "delete_chunks"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c("timeout", "n_search", "sleep_time"),
    args_type = "numeric")

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "EASIN_URL", "DP_R_easin_data_url", FALSE, FALSE,
    "path_EASIN", "DP_R_easin_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999, timeout = timeout))

  # Output file for the merged datasets
  path_out <- fs::path(path_EASIN, paste0(species_key, ".RData"))

  # Ensure that the directory for temporary files exist
  fs::dir_create(fs::path(path_EASIN, "file_parts"))

  # Check if species data already available
  out_file_exists <- ecokit::check_data(path_out, warning = FALSE)

  if (out_file_exists && verbose) {
    ecokit::cat_time("Output file already exists")
  }

  # # ..................................................................... ###

  if (isFALSE(out_file_exists)) {
    # Download chunk data
    if (verbose) {
      ecokit::cat_time("Download chunk data")
    }
    chunk_n <- 0L

    repeat {
      chunk_n <- chunk_n + 1

      object_out <- paste0(
        species_key, "_", stringr::str_pad(chunk_n, width = 5, pad = "0"))

      path_part <- fs::path(
        path_EASIN, "file_parts", paste0(object_out, ".RData"))

      if (ecokit::check_data(path_part, warning = FALSE)) {
        next
      }

      # nolint start
      Skip <- (chunk_n - 1) * n_search
      # `exclude/dps/1` excludes GBIF observations
      URL <- stringr::str_glue(
        "{EASIN_URL}/{species_key}/exclude/dps/1/{Skip}/{n_search}")
      # nolint end

      download_try <- 0
      while (download_try <= n_attempts) {
        download_try <- download_try + 1

        chunk_data <- try(
          {
            chunk_data_raw <- RCurl::getURL(URL, .mapUnicode = FALSE)
            # Error <- stringr::str_detect(
            #   chunk_data_raw, pattern = "An error occurred while")
            n_obs <- stringr::str_detect(
              chunk_data_raw, pattern = "There are no results based on your")
            chunk_data_raw
          },
          silent = TRUE)

        if (n_obs) {
          break
        }

        chunk_data <- tibble::tibble(
          jsonlite::fromJSON(chunk_data, flatten = TRUE))

        if (inherits(chunk_data, "data.frame")) {
          if (verbose) {
            ecokit::cat_time(
              paste0(
                "chunk ", chunk_n, " - attempt ", download_try), level = 1L)
          }
          break
        }
      }

      if (inherits(chunk_data, "data.frame")) {
        ecokit::save_as(
          object = chunk_data, object_name = object_out, out_path = path_part)

        if (nrow(chunk_data) < n_search) {
          break
        }
      } else {
        break
      }

      # sleep at each chunk download
      Sys.sleep(sleep_time)
    }


    if (verbose) {
      ecokit::cat_time("Save taxa data")
    }

    chunk_list <- list.files(
      fs::path(path_EASIN, "file_parts"),
      paste0("^", species_key, ".+"), full.names = TRUE)

    ecokit::save_as(
      object = purrr::map_dfr(chunk_list, ecokit::load_as),
      object_name = species_key, out_path = path_out)

    if (delete_chunks) {
      if (verbose) {
        ecokit::cat_time("Delete chunks")
      }
      fs::file_delete(chunk_list)
    }
  }

  # # ..................................................................... ###

  if (return_data && file.exists(path_out)) {
    return(ecokit::load_as(path_out))
  } else {
    return(invisible(NULL))
  }
}
