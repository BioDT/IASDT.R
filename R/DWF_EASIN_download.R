## |------------------------------------------------------------------------| #
# EASIN_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 3

EASIN_download <- function(
    species_key, timeout = 200, verbose = FALSE, env_file = ".env",
    n_search = 1000, n_attempts = 10, sleep_time = 5,
    delete_chunks = TRUE, return_data = FALSE) {

  # # ..................................................................... ###

  if (is.null(species_key)) {
    IASDT.R::stop_ctx("species_key cannot be NULL", species_key = species_key)
  }

  Path_EASIN <- NULL

  # # ..................................................................... ###

  # Checking arguments ----
  if (verbose) {
    IASDT.R::cat_time("Checking arguments")
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("species_key", "env_file"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("return_data", "verbose", "delete_chunks"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("timeout", "n_search", "sleep_time"))

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "EASIN_URL", "DP_R_EASIN_data_url", FALSE, FALSE,
    "Path_EASIN", "DP_R_EASIN_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999, timeout = timeout))

  # Output file for the merged datasets
  Path_Out <- fs::path(Path_EASIN, paste0(species_key, ".RData"))

  # Ensure that the directory for temporary files exist
  fs::dir_create(fs::path(Path_EASIN, "FileParts"))

  # Check if species data already available
  OutFileExist <- IASDT.R::check_data(Path_Out, warning = FALSE)

  if (OutFileExist && verbose) {
    IASDT.R::cat_time("Output file already exists")
  }

  # # ..................................................................... ###

  if (isFALSE(OutFileExist)) {
    # Download chunk data
    if (verbose) {
      IASDT.R::cat_time("Download chunk data")
    }
    Chunk <- 0L

    repeat {
      Chunk <- Chunk + 1

      Obj_Out <- paste0(
        species_key, "_", stringr::str_pad(Chunk, width = 5, pad = "0"))

      Path_Part <- fs::path(
        Path_EASIN, "FileParts", paste0(Obj_Out, ".RData"))

      if (IASDT.R::check_data(Path_Part, warning = FALSE)) {
        next
      }

      # nolint start
      Skip <- (Chunk - 1) * n_search
      URL <- stringr::str_glue(
        "{EASIN_URL}/{species_key}/exclude/dps/1/{Skip}/{n_search}")
      # nolint end

      DownTry <- 0
      while (DownTry <= n_attempts) {
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
          if (verbose) {
            IASDT.R::cat_time(
              paste0("Chunk ", Chunk, " - attempt ", DownTry), level = 1L)
          }
          break
        }
      }

      if (inherits(ChunkDT, "data.frame")) {
        IASDT.R::save_as(
          object = ChunkDT, object_name = Obj_Out, out_path = Path_Part)

        if (nrow(ChunkDT) < n_search) {
          break
        }
      } else {
        break
      }

      # sleep at each chunk download
      Sys.sleep(sleep_time)
    }


    if (verbose) {
      IASDT.R::cat_time("Save taxa data")
    }

    ChunkList <- list.files(
      fs::path(Path_EASIN, "FileParts"),
      paste0("^", species_key, ".+"), full.names = TRUE)

    IASDT.R::save_as(
      object = purrr::map_dfr(ChunkList, IASDT.R::load_as),
      object_name = species_key, out_path = Path_Out)

    if (delete_chunks) {
      if (verbose) {
        IASDT.R::cat_time("Delete chunks")
      }
      fs::file_delete(ChunkList)
    }
  }

  # # ..................................................................... ###

  if (return_data && file.exists(Path_Out)) {
    return(IASDT.R::load_as(Path_Out))
  } else {
    return(invisible(NULL))
  }
}
