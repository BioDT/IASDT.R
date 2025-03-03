## |------------------------------------------------------------------------| #
# Efforts_Request ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 2
#' @export

Efforts_Request <- function(
  EnvFile = ".env", NCores = 3L, StartYear = 1981L,
  Renviron = ".Renviron", Boundaries = c(-30, 50, 25, 75)) {

  # # ..................................................................... ###

  # In earlier tries, requesting all vascular plants occurrences in a single
  # request returned 80 GB compressed file. The extracted "occurrences.txt" is
  # >280 GB (220M observations).
  #
  # The following makes individual request for each vascular plant order. This
  # can take up to 5 hours for the data to be ready

  # # ..................................................................... ###

  .StartTimeRequest <- lubridate::now(tzone = "CET")

  if (missing(NCores) || !is.numeric(NCores) || NCores < 1) {
    stop("`NCores` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(StartYear) || StartYear <= 1950) {
    stop("`StartYear` must be a positive integer after 1950", call. = FALSE)
  }

  if (!is.numeric(Boundaries) || length(Boundaries) != 4) {
    stop("`Boundaries` must be a numeric vector of length 4.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- DownDetails <- orderKey <- Size <- NumberDatasets <-
    TotalRecords <- Path_Efforts <- NULL
  # # ..................................................................... ###

  IASDT.R::GBIF_Check(Renviron = Renviron)

  # # ..................................................................... ###

  # Environment variables ----
  
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Prepare working on parallel -----

  # GBIF allows only 3 parallel requests. Here I wait until previous request
  # is finished.
  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", min(NCores, 3), "` cores."),
    Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(NCores, 3))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Requesting efforts data on parallel -----

  "Requesting efforts data on parallel (This may take up to 4 hours)" %>% 
    IASDT.R::CatTime(Level = 1)

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
          Request_Path <- IASDT.R::Path(
            Path_Efforts, "Requests", paste0(Request_ID, ".RData"))

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
                      Bottom = Boundaries[3], Top = Boundaries[4])),
                  format = "SIMPLE_CSV")

                IASDT.R::SaveAs(
                  InObj = Down, OutObj = Request_ID, OutPath = Request_Path)
              },
              error = function(e) {
                stop(
                  "Failed to request data for taxonKey ", .x, ": ",
                  conditionMessage(e), call. = FALSE)
              })
          }

          # Waiting for data to be ready
          rgbif::occ_download_wait(Down, quiet = TRUE)

          return(Down)
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf,
          globals = c("Path_Efforts", "Boundaries", "StartYear"),
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

  IASDT.R::CatTime("Requesting efforts data was finished", Level = 2)

  # # ..................................................................... ###

  # Save efforts request data ------
  IASDT.R::CatTime("Save efforts request data", Level = 1)

  save(
    Efforts_AllRequests,
    file = IASDT.R::Path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeRequest,
    Prefix = "Requesting efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}
