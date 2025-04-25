## |------------------------------------------------------------------------| #
# efforts_request ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 2
#' @export

efforts_request <- function(
    env_file = ".env", n_cores = 3L, start_year = 1981L,
    r_environ = ".Renviron", boundaries = c(-30, 50, 25, 75)) {

  # # ..................................................................... ###

  # In earlier tries, requesting all vascular plants occurrences in a single
  # request returned 80 GB compressed file. The extracted "occurrences.txt" is
  # >280 GB (220M observations).
  #
  # The following makes individual request for each vascular plant order. This
  # can take up to 5 hours for the data to be ready

  # # ..................................................................... ###

  .StartTimeRequest <- lubridate::now(tzone = "CET")

  if (missing(n_cores) || !is.numeric(n_cores) || n_cores < 1) {
    IASDT.R::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores)
  }

  if (!is.numeric(start_year) || start_year <= 1950) {
    IASDT.R::stop_ctx(
      "`start_year` must be a positive integer after 1950",
      start_year = start_year)
  }

  if (!is.numeric(boundaries) || length(boundaries) != 4) {
    IASDT.R::stop_ctx(
      "`boundaries` must be a numeric vector of length 4.",
      boundaries = boundaries)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- DownDetails <- orderKey <- Size <- NumberDatasets <-
    TotalRecords <- Path_Efforts <- NULL
  # # ..................................................................... ###

  IASDT.R::GBIF_check(r_environ = r_environ)

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Prepare working in parallel -----

  # GBIF allows only 3 parallel requests. Here I wait until previous request
  # is finished.
  IASDT.R::cat_time(
    paste0("Prepare working in parallel using ", min(n_cores, 3), " cores"),
    level = 1)

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(n_cores, 3))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Requesting efforts data in parallel -----

  IASDT.R::cat_time(
    "Requesting efforts data in parallel (This may take up to 4 hours)",
    level = 1)

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
        .f = ~ {
          Request_ID <- paste0("Request_", .x)
          Request_Path <- IASDT.R::path(
            Path_Efforts, "Requests", paste0(Request_ID, ".RData"))

          if (file.exists(Request_Path)) {
            # load previous request
            Down <- IASDT.R::load_as(Request_Path)
          } else {
            # Attempt the request with error handling
            tryCatch(
              {
                # Make data request
                Down <- rgbif::occ_download(
                  rgbif::pred_in("taxonKey", .x),
                  # Only with coordinates & no spatial issues
                  rgbif::pred("hasCoordinate", TRUE),
                  rgbif::pred("hasGeospatialIssue", FALSE),
                  # Only after (>=) a certain year
                  rgbif::pred_gte("year", start_year),
                  # Only within specific boundaries
                  rgbif::pred_within(
                    value = IASDT.R::boundary_to_WKT(
                      left = boundaries[1], right = boundaries[2],
                      bottom = boundaries[3], top = boundaries[4])),
                  format = "SIMPLE_CSV")

                IASDT.R::save_as(
                  object = Down, object_name = Request_ID,
                  out_path = Request_Path)
              },
              error = function(e) {
                IASDT.R::stop_ctx(
                  paste0(
                    "Failed to request data for taxonKey ", .x, ": ",
                    conditionMessage(e)))
              })
          }

          # Waiting for data to be ready
          rgbif::occ_download_wait(Down, quiet = TRUE)

          return(Down)
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf,
          globals = c("Path_Efforts", "boundaries", "start_year"),
          packages = c("dplyr", "IASDT.R", "rgbif"))
      )
    ) %>%
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
      dplyr::across(
        c("Created", "Modified", "EraseAfter"),
        lubridate::as_date)) %>%
    dplyr::ungroup() %>%
    # how to cite data
    dplyr::mutate(Citation = purrr::map_chr(Request, attr, "citation"))

  IASDT.R::cat_time("Requesting efforts data was finished", level = 2)

  # # ..................................................................... ###

  # Save efforts request data ------
  IASDT.R::cat_time("Save efforts request data", level = 1)

  save(
    Efforts_AllRequests,
    file = IASDT.R::path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::cat_time("Stopping cluster", level = 1)
  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTimeRequest,
    prefix = "Requesting efforts data took ", level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}
