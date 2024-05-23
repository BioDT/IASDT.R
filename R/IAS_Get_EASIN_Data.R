## |------------------------------------------------------------------------| #
# Get_EASIN_Data ----
## |------------------------------------------------------------------------| #

#' Extract EASIN data using EASIN ID
#'
#' Extract EASIN data using EASIN ID
#' @param SpKey character; EASIN taxon ID
#' @param NSearch number of observations to extract each time; default: 1000 (the current maximum per request)
#' @param SleepPart integer; time in seconds to wait between each request. Default: 1.
#' @param Path_Raw character; a valid path for the location to store the output data. A subfolder named `FileParts` will be created for temporary files for individual requests.
#' @param ReturnVal logical; should the data be returned or only saved to disk. Default: `FALSE`
#' @param SleepAfter integer; time in seconds to wait after finishing with downloading taxon data. Default: 10.
#' @name Get_EASIN_Data
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @details
#' A function to extract EASIN data for a given EASIN_ID

Get_EASIN_Data <- function(
    SpKey, NSearch = 1000, SleepPart = 1, Path_Raw = "Data/EASIN/Raw",
    ReturnVal = FALSE, SleepAfter = 5) {

  withr::local_options(list(scipen = 999, timeout = 200))

  Path_Out <- file.path(Path_Raw, paste0(SpKey, ".RData"))
  fs::dir_create(file.path(Path_Raw, "FileParts"))
  SpeciesOkay <- TRUE
  OutFileExist <- FALSE

  if (file.exists(Path_Out)) {
    OutFileExist <- TRUE
  } else {
    # Update 08.03.2024
    # The following does not work as of 08.03.2024 due to changes in the EASIN API
    # URL <- "https://easin.jrc.ec.europa.eu/apixg/geoxg/speciesid/{SpKey}/layertype/grid/skip/{Skip}/take/{NSearch}" %>%
    #    stringr::str_glue()
    # See https://easin.jrc.ec.europa.eu/apixg/home/geoqueries/ for help on how to formulate the URL
    # `exclude/dps/1/` excludes GBIF data

    # Looping over data chunks
    BaseURL <- "https://easin.jrc.ec.europa.eu/apixg/geoxg"
    Skip <- 0         # Skip = 0; start at the first presence grid
    ID <- 0           # iteration ID
    # DT <- list()      # List object to save the data

    while (TRUE) {
      ID <- ID + 1
      Obj_Out <- paste0(SpKey, "_", ID)
      Path_Part <- file.path(Path_Raw, "FileParts", paste0(Obj_Out, ".RData"))

      if (file.exists(Path_Part)) {
        Skip <- Skip + NSearch
        next
      } else {

        ChunkDT_URL <- "{BaseURL}/{SpKey}/exclude/dps/1/{Skip}/{NSearch}" %>%
          stringr::str_glue()
        DownTry <- 0

        while (TRUE) {
          DownTry <- DownTry + 1
          # .mapUnicode = FALSE to ignore encoding issues in some references
          ChunkDT <- RCurl::getURL(ChunkDT_URL, .mapUnicode = FALSE)
          ReDown <- inherits(ChunkDT, "try-error") ||
            stringr::str_detect(
              string = ChunkDT, pattern = "An error occurred while retrieving")
          if (magrittr::not(ReDown)) break()
          if (ReDown && DownTry > 10) {
            SpeciesOkay <- FALSE
            break()
          }
          Sys.sleep(1)
        }

        if (magrittr::not(SpeciesOkay)) {
          "data for {SpKey} can not be downloaded. Download failed after {DownTry} trials. Check server status." %>%
            stringr::str_glue() %>%
            IASDT.R::CatTime()
          break()
        }

        if (stringr::str_detect(
          string = ChunkDT,
          pattern = "There are no results based on your search criteria")) {
          IASDT.R::SaveAs(
            InObj = tibble::tibble(), OutObj = SpKey, OutPath = Path_Out)
          break()
        }

        ChunkDT <- tibble::tibble(jsonlite::fromJSON(ChunkDT))
        IASDT.R::SaveAs(InObj = ChunkDT, OutObj = Obj_Out, OutPath = Path_Part)

        if (nrow(ChunkDT) < NSearch) {
          break()
        } else {
          # Wait after each part
          Sys.sleep(SleepPart)
          Skip <- Skip + NSearch
        }
      }
    }
  }

  if (SpeciesOkay) {

    if (OutFileExist) {
      DT <- IASDT.R::LoadAs(Path_Out)
    } else {
      Tiles <- list.files(
        path = file.path(Path_Raw, "FileParts"), full.names = TRUE,
        pattern = paste0(SpKey, "_.+.RData")) %>%
        gtools::mixedsort()

      DT <- purrr::map_dfr(Tiles, IASDT.R::LoadAs)

      IASDT.R::SaveAs(InObj = DT, OutObj = SpKey, OutPath = Path_Out)
      fs::file_delete(Tiles)
    }

    if (ReturnVal) {
      Sys.sleep(SleepAfter)
      return(DT)
    } else {
      Sys.sleep(SleepAfter)
      return(invisible(NULL))
    }
  } else {
    Sys.sleep(SleepAfter)
    return(invisible(NULL))
  }
}
