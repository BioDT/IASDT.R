## |------------------------------------------------------------------------| #
# EASIN_Down ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 3

EASIN_Down <- function(
    SpKey, Timeout = 200, Verbose = FALSE, FromHPC = TRUE, EnvFile = ".env",
    NSearch = 1000, Attempts = 10, SleepTime = 5,
    DeleteChunks = TRUE, ReturnData = FALSE) {

  # # ..................................................................... ###

  if (is.null(SpKey)) {
    stop("SpKey cannot be NULL", call. = FALSE)
  }

  Path_EASIN <- NULL

  # # ..................................................................... ###

  # Checking arguments ----
  if (Verbose) {
    IASDT.R::CatTime("Checking arguments")
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("SpKey", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("ReturnData", "Verbose", "DeleteChunks", "FromHPC"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("Timeout", "NSearch", "SleepTime"))


  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "EASIN_URL", "DP_R_EASIN_URL", FALSE, FALSE,
      "Path_EASIN", "DP_R_EASIN_Interim", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "EASIN_URL", "DP_R_EASIN_URL", FALSE, FALSE,
      "Path_EASIN", "DP_R_EASIN_Interim_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999, timeout = Timeout))

  # Output file for the merged datasets
  Path_Out <- IASDT.R::Path(Path_EASIN, paste0(SpKey, ".RData"))

  # Ensure that the directory for temporary files exist
  fs::dir_create(IASDT.R::Path(Path_EASIN, "FileParts"))

  # Check if species data already available
  OutFileExist <- IASDT.R::CheckData(Path_Out, warning = FALSE)

  if (OutFileExist && Verbose) {
    IASDT.R::CatTime("Output file already exists")
  }

  # # ..................................................................... ###

  if (isFALSE(OutFileExist)) {
    # Download chunk data
    if (Verbose) {
      IASDT.R::CatTime("Download chunk data")
    }
    Chunk <- 0L

    repeat {
      Chunk <- Chunk + 1

      Obj_Out <- paste0(
        SpKey, "_", stringr::str_pad(Chunk, width = 5, pad = "0"))

      Path_Part <- IASDT.R::Path(
        Path_EASIN, "FileParts", paste0(Obj_Out, ".RData"))
      if (IASDT.R::CheckData(Path_Part, warning = FALSE)) {
        next
      }

      # nolint start
      Skip <- (Chunk - 1) * NSearch
      # nolint end

      URL <- stringr::str_glue(
        "{EASIN_URL}/{SpKey}/exclude/dps/1/{Skip}/{NSearch}")

      DownTry <- 0
      while (DownTry <= Attempts) {
        DownTry <- DownTry + 1

        ChunkDT <- try(
          {
            ChunkDT0 <- RCurl::getURL(URL, .mapUnicode = FALSE)
            # Error <- stringr::str_detect(
            #   ChunkDT0, pattern = "An error occurred while")
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
      IASDT.R::Path(Path_EASIN, "FileParts"),
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
