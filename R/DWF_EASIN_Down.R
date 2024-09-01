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
#' @return If `ReturnData` is `TRUE`, returns a dataframe containing all the data
#'   retrieved for the specified EASIN ID. If `ReturnData` is `FALSE`, returns
#'   `NULL`.
#' @export
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [EASIN_Process] function.
#' @details This function extracts EASIN data for a given EASIN_ID, handles
#'   pagination, and ensures that data retrieval is efficient by managing
#'   retries and pauses.

EASIN_Down <- function(
    SpKey, Timeout = 200, Verbose = FALSE,
    BaseURL = "https://easin.jrc.ec.europa.eu/apixg/geoxg",
    NSearch = 1000, Attempts = 10, SleepTime = 5,
    Path_Raw = "datasets/interim/EASIN",
    DeleteChunks = TRUE, ReturnData = FALSE) {

  # # ..................................................................... ###

  if (is.null(SpKey)) {
    stop("SpKey cannot be NULL", call. = FALSE)
  }

  # # ..................................................................... ###

  # Checking arguments ----
  if (Verbose) {
    IASDT.R::CatTime("Checking arguments")
  }

  AllArgs <- ls(envir = environment())
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
    Args = c("Timeout", "NSearch", "SleepTime"))

  # # ..................................................................... ###

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999, timeout = Timeout))

  # Output file for the merged datasets
  Path_Out <- file.path(Path_Raw, paste0(SpKey, ".RData"))

  # Ensure that the directory for temporary files exist
  fs::dir_create(file.path(Path_Raw, "FileParts"))

  # Check if species data already available
  OutFileExist <- file.exists(Path_Out) && IASDT.R::CheckRData(Path_Out)

  if (OutFileExist) {
    if (Verbose) {
      IASDT.R::CatTime("Output file already exists")
    }
  }

  # # ..................................................................... ###

  if (isFALSE(OutFileExist)) {

    # Download chunk data
    if (Verbose) {
      IASDT.R::CatTime("Download chunk data")
    }
    Chunk <- 0L
    PrevChunks <- list.files(file.path(Path_Raw, "FileParts"), SpKey) %>%
      stringr::str_remove_all(paste0(".RData|", SpKey, "_")) %>%
      as.integer()

    while (TRUE) {
      Chunk <- Chunk + 1

      Obj_Out <- paste0(
        SpKey, "_", stringr::str_pad(Chunk, width = 5, pad = "0"))
      Path_Part <- file.path(Path_Raw, "FileParts", paste0(Obj_Out, ".RData"))
      if (file.exists(Path_Part) && IASDT.R::CheckRData(Path_Part)) {
        next
      }
      Skip <- (Chunk - 1) * NSearch
      URL <- stringr::str_glue(
        "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/{NSearch}")

      DownTry <- 0
      while (DownTry <= Attempts) {
        DownTry <- DownTry + 1

        ChunkDT <- try({
          ChunkDT0 <- RCurl::getURL(URL, .mapUnicode = FALSE)
          Error <- stringr::str_detect(
            ChunkDT0, pattern = "An error occurred while")
          NoObs <- stringr::str_detect(
            ChunkDT0, pattern = "There are no results based on your")
          ChunkDT0
        },
        silent = TRUE)

        if (NoObs) {
          break
        }

        ChunkDT <- tibble::tibble(jsonlite::fromJSON(ChunkDT, flatten = TRUE))

        if (inherits(ChunkDT, "data.frame")) {
          if (Verbose) {
            IASDT.R::CatTime(
              paste0("Chunk ", Chunk, " - attempt ", DownTry), Level = 1)
          }
          break
        }
      }


      if (inherits(ChunkDT, "data.frame")) {
        IASDT.R::SaveAs(InObj = ChunkDT, OutObj = Obj_Out, OutPath = Path_Part)

        if (nrow(ChunkDT) < NSearch) {
          break
        }
      } else {
        break
      }

      # sleep at each chunk download
      Sys.sleep(SleepTime)
    }


    if (Verbose) {
      IASDT.R::CatTime("Save taxa data")
    }

    ChunkList <- list.files(
      file.path(Path_Raw, "FileParts"),
      paste0("^", SpKey, ".+"), full.names = TRUE)

    IASDT.R::SaveAs(
      InObj = purrr::map_dfr(ChunkList, IASDT.R::LoadAs),
      OutObj = SpKey, OutPath = Path_Out)

    if (DeleteChunks) {
      if (Verbose) {
        IASDT.R::CatTime("Delete chunks")
      }
      fs::file_delete(ChunkList)
    }
  }

  # # ..................................................................... ###

  if (ReturnData && file.exists(Path_Out)) {
    return(IASDT.R::LoadAs(Path_Out))
  } else {
    return(invisible(NULL))
  }
}
