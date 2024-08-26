## |------------------------------------------------------------------------| #
# EASIN_Down ----
## |------------------------------------------------------------------------| #

#' Extract EASIN data using EASIN ID
#'
#' This function downloads data for a specified species from the EASIN database.
#' It handles pagination automatically, downloading data in chunks until all
#' available data for the species is retrieved. The function also supports
#' pausing between requests to avoid overloading the server.
#' @param BaseURL character; the base URL for downloading EASIN data. Default
#'   is "https://easin.jrc.ec.europa.eu/apixg/geoxg".
#' @param SpKey character; the EASIN taxon ID for which data is to be retrieved.
#'   This parameter cannot be `NULL`.
#' @param Timeout integer; time in seconds before a download attempt times out.
#'   Default is 200.
#' @param Verbose logical. Indicating whether to print messages for the progress
#'   of the function
#' @param NSearch the number of records to attempt to retrieve per request.
#'   Default is 1000, which is the current maximum allowed by the API.
#' @param Attempts integer; maximum number of download attempts per data
#'   chunk. Defaults to 10.
#' @param SleepTime integer; the number of seconds to pause between each data
#'   retrieval request to prevent overloading the server. Default is 5 second.
#' @param Path_Raw character; the path where the raw data files and temporary
#'   parts will be stored. Default is `datasets/interim/EASIN`.
#' @param DeleteChunks logical, indicating whether to delete temporary files for
#'   data chunks from the `FileParts` subdirectory. Defaults to `TRUE`.
#' @param ReturnData logical; if `TRUE`, the function will return the combined
#'   data as a dataframe. If `FALSE` (default), the data will only be saved to
#'   disk and the function will invisibly return `NULL`.
#' @name EASIN_Down
#' @author Ahmed El-Gabbas
#' @return If ReturnData is `TRUE`, returns a dataframe containing all the data
#'   retrieved for the specified EASIN ID. If ReturnData is `FALSE`, returns
#'   `NULL`.
#' @export
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only called from the [EASIN_Processing] function.
#' @details This function extracts EASIN data for a given EASIN_ID, handles
#'   pagination, and ensures that data retrieval is efficient by managing
#'   retries and pauses.

EASIN_Down <- function(
    SpKey, Timeout = 200, Verbose = FALSE,
    BaseURL = "https://easin.jrc.ec.europa.eu/apixg/geoxg",
    NSearch = 1000, Attempts = 10, SleepTime = 5,
    Path_Raw = "datasets/interim/EASIN",
    DeleteChunks = TRUE, ReturnData = FALSE) {

  if (is.null(SpKey)) {
    stop("SpKey cannot be NULL", call. = FALSE)
  }

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("SpKey", "BaseURL", "Path_Raw"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("ReturnData", "Verbose", "DeleteChunks"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("Timeout", "NSearch", "Attempts", "SleepTime"))

  # # ..................................................................... ###

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999, timeout = Timeout))

  # Output file for the merged datasets
  Path_Out <- file.path(Path_Raw, paste0(SpKey, ".RData"))

  # Ensure that the directory for temporary files exist
  fs::dir_create(file.path(Path_Raw, "FileParts"))

  SpeciesOkay <- TRUE

  # Check if species data already available
  OutFileExist <- file.exists(Path_Out)

  if (OutFileExist) {
    if (IASDT.R::CheckRData(Path_Out)) {
      if (Verbose) {
        IASDT.R::CatTime("Output file already exists")
      }
    } else {
      OutFileExist <- FALSE
    }
  }

  if (!OutFileExist) {
    # See https://easin.jrc.ec.europa.eu/apixg/home/geoqueries/ for help on how
    # to formulate the URL. `exclude/dps/1/` excludes GBIF data

    # Looping over data chunks
    Skip <- 0         # Skip = 0; start at the first presence grid
    Chunk <- 0        # iteration ID

    while (TRUE) {
      Chunk <- Chunk + 1
      Obj_Out <- paste0(SpKey, "_", Chunk)
      Path_Part <- file.path(Path_Raw, "FileParts", paste0(Obj_Out, ".RData"))

      if (!file.exists(Path_Part)) {
        # download only observation in this chunk to check for error or if there
        # is no observations or
        ChunkDT_URL1 <- stringr::str_glue(
          "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/1")

        # downloading data chunks for a maximum of `Attempts` times
        DownTry <- 0
        while (TRUE) {
          DownTry <- DownTry + 1

          if (Verbose) {
            IASDT.R::CatTime(
              paste0(SpKey, ": Chunk ", Chunk, " - attempt ", DownTry))
          }

          # Download only one observation for this chunk. `.mapUnicode` = FALSE
          # to ignore encoding issues in some references
          ChunkDT1 <- try(
            RCurl::getURL(ChunkDT_URL1, .mapUnicode = FALSE), silent = TRUE)

          # Break if there is no observations and save empty tibble
          NoObs <- stringr::str_detect(
            string = ChunkDT1, pattern = "There are no results based ")
          if (NoObs) {
            SpeciesOkay <- FALSE
            IASDT.R::SaveAs(
              InObj = tibble::tibble(), OutObj = SpKey, OutPath = Path_Out)
            break
          }

          # Download the current chunk
          ChunkDT_URL <- stringr::str_glue(
            "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/{NSearch}")
          ChunkDT <- try(
            RCurl::getURL(ChunkDT_URL, .mapUnicode = FALSE), silent = TRUE)

          if (
            (inherits(ChunkDT, "try-error") ||
             stringr::str_detect(
               string = ChunkDT,
               pattern = "An error occurred while retrieving")) &&
            DownTry < Attempts) {
            next
          } else {
            break
          }

          if (DownTry >= Attempts) {
            SpeciesOkay <- FALSE
            break
          }

          # sleep at each chunk download
          Sys.sleep(SleepTime)
        }

        # If species is okay, save current chunk data as tibble
        if (SpeciesOkay) {

          ChunkDT <- jsonlite::fromJSON(ChunkDT, flatten = TRUE)

          if (nrow(ChunkDT) > 0) {
            ChunkDT <- tibble::tibble(ChunkDT)
          } else {
            ChunkDT <- tibble::tibble()
          }

          IASDT.R::SaveAs(
            InObj = ChunkDT, OutObj = Obj_Out, OutPath = Path_Part)

          if (nrow(ChunkDT) < NSearch) {
            break
          } else {
            Skip <- Skip + NSearch
          }
        } else {
          break
        }
      }
    }

    if (SpeciesOkay) {
      # Merge chunk data together in a single tibble
      Tiles <- list.files(
        path = file.path(Path_Raw, "FileParts"), full.names = TRUE,
        pattern = paste0(SpKey, "_.+.RData")) %>%
        gtools::mixedsort()

      DT <- purrr::map_dfr(Tiles, IASDT.R::LoadAs) %>%
        dplyr::distinct()

      if (Verbose) {
        IASDT.R::CatTime(
          paste0(SpKey, ": ", nrow(DT), " observations"), Level = 1)
      }

      IASDT.R::SaveAs(InObj = DT, OutObj = SpKey, OutPath = Path_Out)

      if (DeleteChunks) {
        fs::file_delete(Tiles)
      }
    }
  }

  if (ReturnData && SpeciesOkay) {
    if (OutFileExist) {
      DT <- IASDT.R::LoadAs(Path_Out)
    }
    return(DT)
  } else {
    return(invisible(NULL))
  }
}
