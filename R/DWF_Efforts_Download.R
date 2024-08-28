## |------------------------------------------------------------------------| #
# Efforts_Download ----
## |------------------------------------------------------------------------| #

#' Download and Manage GBIF Data for Vascular Plant Orders
#'
#' This function handles the downloading of GBIF data in parallel, checks the
#' validity of downloaded files, and stores the data in specified directories.
#' If data has already been downloaded, it validates the files instead of
#' downloading them again.
#' @param NCores Integer. Number of cores to use for parallel processing.  Must
#'   be a positive integer. This directory must exist or be created beforehand.
#' @param Path_Raw Character. Path where the raw downloaded data will be saved.
#'   This directory must exist or be created beforehand.
#' @param Path_Interim Character. Path where the interim CSV files will be
#'   saved. This directory must exist or be created beforehand.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved. This directory must exist or be created beforehand.
#' @name Efforts_Download
#' @author Ahmed El-Gabbas
#' @return A data frame (`Efforts_AllRequests`) with updated download paths and
#'   interim file paths.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Sampling_Efforts] function.
#' @export

Efforts_Download <- function(NCores = 6, Path_Raw, Path_Interim, Path_Efforts) {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  DownPath <- Request <- NULL

  # # ..................................................................... ###

  # Validate NCores
  if (missing(NCores) || !is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  # Validate directory paths
  if (!is.character(Path_Raw) || !dir.exists(Path_Raw)) {
    stop("`Path_Raw` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Interim) || !dir.exists(Path_Interim)) {
    stop("`Path_Interim` must be a valid directory path.", call. = FALSE)
  }

  if (!is.character(Path_Efforts) || !dir.exists(Path_Efforts)) {
    stop("`Path_Efforts` must be a valid directory path.", call. = FALSE)
  }

  IASDT.R::CheckCommands("unzip")

  # # ..................................................................... ###

  Path_Efforts_Request <- file.path(Path_Efforts, "Efforts_AllRequests.RData")
  if (!file.exists(Path_Efforts_Request)) {
    stop(
      paste0(
        "The path for the `Efforts_AllRequests` data does not exist: ",
        Path_Efforts_Request),
      call. = FALSE)
  }

  Efforts_AllRequests <- IASDT.R::LoadAs(Path_Efforts_Request)

  # # ..................................................................... ###

  ## Prepare working on parallel -----
  #
  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  snow::clusterEvalQ(
    cl = c1,
    expr = IASDT.R::LoadPackages(
      List = c("dplyr", "IASDT.R", "rgbif", "stringr")))

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  IASDT.R::CatTime("Downloading/checking efforts data", Level = 1)

  Efforts_AllRequests <- Efforts_AllRequests %>%
    dplyr::mutate(
      # Download datasets on parallel
      DownPath = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          DownFile <- file.path(Path_Raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not exist download it
          if (file.exists(DownFile)) {
            FileOkay <- tryCatch({
              system2(
                "unzip", args = c("-t", DownFile), stdout = TRUE, stderr = TRUE)
            }, error = function(e) {
              message("Error during file validation: ", conditionMessage(e))
              return(NULL)
            }) %>%
              stringr::str_detect("No errors detected in compressed data") %>%
              any()

            if (FileOkay) {
              Success <- TRUE
            } else {
              Success <- FALSE
            }

          } else {
            Success <- FALSE
          }

          # Try downloading data for a max of 3 attempts, each with 20 mins
          # time out
          withr::local_options(list(timeout = 1200))

          Attempt <- 1
          Attempts <- 3

          while (isFALSE(Success) && (Attempt <= Attempts)) {
            tryCatch({
              suppressMessages(
                rgbif::occ_download_get(
                  key = .x, path = Path_Raw, overwrite = TRUE))

              ZipStatus <- system2(
                "unzip", args = c("-t", DownFile),
                stdout = TRUE, stderr = TRUE) %>%
                stringr::str_detect("No errors detected in compressed data") %>%
                any()

              # Ensure Success is only TRUE if both the zip file exists and
              # passes integrity check
              Success <- file.exists(DownFile) && ZipStatus

            },
            error = function(e) {
              if (Attempt < Attempts) {
                Attempt <- Attempt + 1
              } else {
                stop(
                  paste0(
                    "Failed to download data after ", Attempts, " attempts: ",
                    conditionMessage(e)),
                  call. = FALSE)
              }
            })
          }

          return(DownFile)

        }, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)))

  save(Efforts_AllRequests,
       file = file.path(Path_Efforts, "Efforts_AllRequests.RData"))

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeDown, CatInfo = FALSE,
    Prefix = "Downloading efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(Efforts_AllRequests)
}
