## |------------------------------------------------------------------------| #
# Efforts_Request ----
## |------------------------------------------------------------------------| #

#' Request and Manage GBIF Data for Vascular Plant Orders
#'
#' This function requests GBIF data for each vascular plant order, processes it
#' in parallel, and manages the data download and storage. If data is already
#' available, it loads the data instead of making new requests.
#' @param NCores Integer. The number of cores to use for parallel processing
#'   (between 1 and 3). Must be a positive integer. Defaults to 3.
#' @param Path_Requests Character. The directory path to save individual request
#'   files.
#' @param Path_Efforts Character. The directory path to save the final compiled
#'   data.
#' @param StartYear Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be requested from GBIF. Default is
#'   `1981`, which matches the year ranges of CHELSA current climate data.
#' @param Boundaries Numeric vector of length 4. Specifies geographical
#'   boundaries for the requested GBIF data in the order: Left, Right, Bottom,
#'   Top. Defaults to c(-30, 50, 25, 75).
#' @return The function returns the GBIF data requests processed and stored in
#'   the specified directories.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Efforts_Process] function.
#' @author Ahmed El-Gabbas
#' @name Efforts_Request
#' @export

Efforts_Request <- function(
    NCores = 3, Path_Requests, Path_Efforts,
    StartYear = 1981, Boundaries = c(-30, 50, 25, 75)) {

  # In earlier tries, requesting all vascular plants occurrences in a single
  # request returned 80 GB compressed file. The extracted "occurrences.txt" is
  # >280 GB (220M observations).
  #
  # The following makes individual request for each vascular plant order. This
  # can take up to 5 hours for the data to be ready

  .StartTimeRequest <- lubridate::now(tzone = "CET")

  if (missing(NCores) || !is.numeric(NCores) || NCores < 1 || NCores > 3) {
    stop("`NCores` must be a positive integer between 1 and 3.", call. = FALSE)
  }

  if (!is.character(Path_Requests) || !dir.exists(Path_Requests)) {
    stop("`Path_Requests` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Efforts) || !dir.exists(Path_Efforts)) {
    stop("`Path_Efforts` must be a valid directory path.", call. = FALSE)
  }

  if (!is.numeric(StartYear) || StartYear <= 1950) {
    stop("`StartYear` must be a positive integer after 1950")
  }

  if (!is.numeric(Boundaries) || length(Boundaries) != 4) {
    stop("`Boundaries` must be a numeric vector of length 4.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- DownDetails <- orderKey <- Size <- NumberDatasets <-
    TotalRecords <- Created <- Modified <- EraseAfter <- NULL

  # # ..................................................................... ###

  # Prepare working on parallel -----
  
  # GBIF allows only 3 parallel requests. Here I wait until previous request
  # is finished.
  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 1)

  withr::local_options(future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  c1 <- snow::makeSOCKcluster(min(NCores, 3))
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  snow::clusterEvalQ(
    cl = c1,
    expr = IASDT.R::LoadPackages(List = c("dplyr", "IASDT.R", "rgbif")))

  # # ..................................................................... ###

  # Requesting efforts Data on parallel -----
  IASDT.R::CatTime("Requesting efforts Data on parallel", Level = 1)
  IASDT.R::CatTime("This may take up to 4 hours", Level = 2)

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
          Request_Path <- file.path(Path_Requests, paste0(Request_ID, ".RData"))

          if (file.exists(Request_Path)) {
            # load previous request
            Down <- IASDT.R::LoadAs(Request_Path)
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
                  rgbif::pred_gte("year", StartYear),
                  # Only within specific boundaries
                  rgbif::pred_within(
                    value = IASDT.R::DownBoundary(
                      Left = Boundaries[1], Right = Boundaries[2],
                      Bottom = Boundaries[3], Top = Boundaries[4])
                  ),
                  format = "SIMPLE_CSV"
                )

                IASDT.R::SaveAs(
                  InObj = Down, OutObj = Request_ID, OutPath = Request_Path)

              },
              error = function(e) {
                stop(
                  paste0(
                    "Failed to request data for taxonKey ", .x, ": ",
                    conditionMessage(e)
                  ),
                  call. = FALSE
                )
              }
            )
          }

          # Waiting for data to be ready
          rgbif::occ_download_wait(Down, quiet = TRUE)

          return(Down)
        },
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
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
      dplyr::across(c(Created, Modified, EraseAfter), lubridate::as_date)) %>%
    dplyr::ungroup() %>%
    # how to cite data
    dplyr::mutate(Citation = purrr::map_chr(Request, attr, "citation"))

  # # ..................................................................... ###

  # Save efforts request data ------
  IASDT.R::CatTime("Save efforts request data", Level = 1)

  save(Efforts_AllRequests,
    file = file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeRequest, 
    Prefix = "Requesting efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(Efforts_AllRequests)
}
