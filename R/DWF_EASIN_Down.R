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
  OutFileExist <- file.exists(Path_Out)

  if (OutFileExist && IASDT.R::CheckRData(Path_Out)) {
    if (Verbose) {
      IASDT.R::CatTime("Output file already exists")
    }
  } else {
    OutFileExist <- FALSE
  }

  # # ..................................................................... ###

  if (isFALSE(OutFileExist)) {

    # Get the number of chunks to download

    if (Verbose) {
      IASDT.R::CatTime("Get the number of chunks to download")
    }

    Skip <- Chunk <- 0
    while (TRUE) {
      Chunk <- Chunk + 1
      if (Verbose) {
        IASDT.R::CatTime(paste0("Inspecting chunk # ", Chunk), Level = 1)
      }
      ChunkDT_URL1 <- stringr::str_glue(
        "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/1")
      NoObs <- try(
        RCurl::getURL(ChunkDT_URL1, .mapUnicode = FALSE), silent = TRUE) %>%
        stringr::str_detect(pattern = "There are no results based ")
      if (NoObs) {
        break
      }
      Skip <- Skip + NSearch
    }

    # ---------------------------------- #

    # Download chunk data
    if (Verbose) {
      IASDT.R::CatTime("Download chunk data")
    }

    DT_Chunks <- tibble::tibble(ID = seq_len(Chunk - 1)) %>%
      dplyr::mutate(
        ChunkDT = purrr::map_chr(
          .x = ID,
          .f = ~{
            Obj_Out <- paste0(
              SpKey, "_", stringr::str_pad(.x, width = 5, pad = "0"))
            Path_Part <- file.path(
              Path_Raw, "FileParts", paste0(Obj_Out, ".RData"))

            if (file.exists(Path_Part)) {
              return(Path_Part)
            }

            Skip <- (.x - 1) * NSearch
            URL <- stringr::str_glue(
              "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/{NSearch}")
            DownTry <- 0

            while (DownTry <= Attempts) {
              DownTry <- DownTry + 1
              if (Verbose) {
                IASDT.R::CatTime(
                  paste0("Chunk ", .x, " - attempt ", DownTry), Level = 1)
              }

              ChunkDT <- try(
                RCurl::getURL(URL, .mapUnicode = FALSE), silent = TRUE)

              if (
                inherits(ChunkDT, "try-error") ||
                stringr::str_detect(
                  string = ChunkDT,
                  pattern = "An error occurred while retrieving")) {
                next
              }

              ChunkDT <- jsonlite::fromJSON(ChunkDT, flatten = TRUE)

              if (nrow(ChunkDT) > 0) {
                ChunkDT <- tibble::tibble(ChunkDT)
              } else {
                ChunkDT <- tibble::tibble()
              }

              # sleep at each chunk download
              Sys.sleep(SleepTime)

              if (inherits(ChunkDT, "data.frame")) {
                break
              }
            }

            # sleep at each chunk download
            Sys.sleep(SleepTime)

            IASDT.R::SaveAs(
              InObj = ChunkDT, OutObj = Obj_Out, OutPath = Path_Part)

            return(Path_Part)
          }))


    # Merge chunk data
    if (Verbose) {
      IASDT.R::CatTime("Merge chunk data")
    }

    DT <- dplyr::pull(DT_Chunks, ChunkDT) %>%
      purrr::map_dfr(IASDT.R::LoadAs) %>%
      dplyr::distinct()

    if (Verbose) {
      IASDT.R::CatTime(
        paste0(SpKey, ": ", nrow(DT), " observations"), Level = 1)
    }

    IASDT.R::SaveAs(InObj = DT, OutObj = SpKey, OutPath = Path_Out)

    if (DeleteChunks) {
      fs::file_delete(DT_Chunks$ChunkDT)
    }
  }

  # # ..................................................................... ###

  if (ReturnData && file.exists(Path_Out)) {
    return(IASDT.R::LoadAs(Path_Out))
  } else {
    return(invisible(NULL))
  }
}
