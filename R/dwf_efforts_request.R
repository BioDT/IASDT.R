## |------------------------------------------------------------------------| #
# efforts_request ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 2
#' @export

efforts_request <- function(
    env_file = ".env", n_cores = 3L, strategy = "multisession",
    start_year = 1981L, r_environ = ".Renviron",
    boundaries = c(-30, 50, 25, 75)) {

  # # ..................................................................... ###

  # In earlier tries, requesting all vascular plants occurrences in a single
  # request returned 80 GB compressed file. The extracted "occurrences.txt" is
  # >280 GB (220M observations).
  #
  # The following makes individual request for each vascular plant order. This
  # can take up to 5 hours for the data to be ready

  # # ..................................................................... ###

  .start_time_request <- lubridate::now(tzone = "CET")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  if (!is.numeric(start_year) || start_year <= 1950) {
    ecokit::stop_ctx(
      "`start_year` must be a positive integer after 1950",
      start_year = start_year, include_backtrace = TRUE)
  }

  if (!is.numeric(boundaries) || length(boundaries) != 4) {
    ecokit::stop_ctx(
      "`boundaries` must be a numeric vector of length 4.",
      boundaries = boundaries, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  request <- download_details <- orderKey <- size <- n_datasets <-
    total_records <- path_efforts <- NULL
  # # ..................................................................... ###

  # Ensure that GBIF access information is available
  ecokit::check_gbif(r_environ = r_environ)

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_efforts", "DP_R_efforts_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("dplyr", "ecokit", "rgbif", "fs"), strategy = strategy)

  # # ..................................................................... ###

  # Prepare working in parallel -----

  # GBIF allows only 3 parallel requests. Here I wait until previous request
  # is finished.

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, 3), level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Requesting efforts data in parallel -----

  ecokit::cat_time(
    "Requesting efforts data in parallel (This may take up to 4 hours)",
    level = 1L)

  # Extract taxonomic info for vascular plants orders
  selected_columns <- c(
    "class", "classKey", "order", "orderKey", "numDescendants")

  efforts_all_requests <- rgbif::name_backbone("Tracheophyta") %>%
    dplyr::pull("phylumKey") %>%
    rgbif::name_lookup(rank = "ORDER", higherTaxonKey = .) %>%
    # Get info on order names
    magrittr::extract2("data") %>%
    dplyr::select(tidyselect::all_of(selected_columns)) %>%
    dplyr::mutate(
      request = furrr::future_map(
        .x = orderKey,
        .f = ~ {
          request_id <- paste0("request_", .x)
          request_path <- fs::path(
            path_efforts, "requests", paste0(request_id, ".RData"))

          if (file.exists(request_path)) {
            # load previous request
            request_down <- ecokit::load_as(request_path)
          } else {
            # Attempt the request with error handling
            tryCatch(
              {
                # Make data request
                request_down <- rgbif::occ_download(
                  rgbif::pred_in("taxonKey", .x),
                  # Only with coordinates & no spatial issues
                  rgbif::pred("hasCoordinate", TRUE),
                  rgbif::pred("hasGeospatialIssue", FALSE),
                  # Only after (>=) a certain year
                  rgbif::pred_gte("year", start_year),
                  # Only within specific boundaries
                  rgbif::pred_within(
                    value = ecokit::boundary_to_wkt(
                      left = boundaries[1], right = boundaries[2],
                      bottom = boundaries[3], top = boundaries[4])),
                  format = "SIMPLE_CSV")

                ecokit::save_as(
                  object = request_down, object_name = request_id,
                  out_path = request_path)
              },
              error = function(e) {
                ecokit::stop_ctx(
                  paste0(
                    "Failed to request data for taxonKey ", .x, ": ",
                    conditionMessage(e)),
                  include_backtrace = TRUE)
              })
          }

          # Waiting for data to be ready
          rgbif::occ_download_wait(request_down, quiet = TRUE)

          return(request_down)
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf, packages = pkg_to_export,
          globals = c("path_efforts", "boundaries", "start_year")))
    ) %>%
    dplyr::rowwise() %>%
    # Add columns for metadata
    dplyr::mutate(
      download_details = list(rgbif::occ_download_wait(request, quiet = TRUE)),
      # Extract some info from metadata
      download_key = download_details$key,
      DOI = download_details$doi,
      created_time = download_details$created,
      modified_time = download_details$modified,
      erase_after = download_details$eraseAfter,
      download_link = download_details$downloadLink,
      size = download_details$size,
      total_records = download_details$totalRecords,
      n_datasets = download_details$numberDatasets,
      status = download_details$status,
      # size of data in megabytes
      size = as.numeric(size) / (1024 * 1024),
      # Convert some columns to numeric (double)
      dplyr::across(size:n_datasets, as.numeric),
      # Convert some columns to integer
      dplyr::across(total_records:n_datasets, as.integer),
      # Convert some columns to date type
      dplyr::across(
        c("created_time", "modified_time", "erase_after"),
        lubridate::as_date)) %>%
    dplyr::ungroup() %>%
    # how to cite data
    dplyr::mutate(citation = purrr::map_chr(request, attr, "citation"))

  ecokit::cat_time("Requesting efforts data was finished", level = 2L)

  # # ..................................................................... ###

  # Save efforts request data ------
  ecokit::cat_time("Save efforts request data", level = 1L)

  save(
    efforts_all_requests,
    file = fs::path(path_efforts, "efforts_all_requests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time_request,
    prefix = "Requesting efforts data took ", level = 1L)

  # # ..................................................................... ###

  return(invisible(NULL))
}
