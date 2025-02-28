## |------------------------------------------------------------------------| #
# Efforts_Download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 3
#' @export

Efforts_Download <- function(NCores = 6L, FromHPC = TRUE, EnvFile = ".env") {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- Path_Efforts <- Path_Raw <- NULL

  # # ..................................................................... ###

  # Validate NCores
  if (missing(NCores) || !is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  Commands <- "unzip"
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::CheckCommands)
  if (!all(CommandsAvail)) {
    Missing <- paste(Commands[!CommandsAvail], collapse = " + ")
    stop("The following command(s) are not available: ", Missing, call. = FALSE)
  }

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts", FALSE, FALSE,
      "Path_Raw", "DP_R_Efforts_Raw", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts_Local", FALSE, FALSE,
      "Path_Raw", "DP_R_Efforts_Raw_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  Path_Efforts_Request <- IASDT.R::Path(
    Path_Efforts, "Efforts_AllRequests.RData")

  if (!file.exists(Path_Efforts_Request)) {
    stop(
      "The path for the `Efforts_AllRequests` data does not exist: ",
      Path_Efforts_Request, call. = FALSE)
  }

  Efforts_AllRequests <- IASDT.R::LoadAs(Path_Efforts_Request)

  # # ..................................................................... ###

  ## Prepare working on parallel -----

  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  IASDT.R::CatTime("Downloading/checking efforts data", Level = 1)

  Efforts_AllRequests <- Efforts_AllRequests %>%
    dplyr::mutate(
      # Download datasets on parallel
      DownPath = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          DownFile <- IASDT.R::Path(Path_Raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not download it
          if (file.exists(DownFile)) {
            Success <- IASDT.R::CheckZip(DownFile)
            if (isFALSE(Success)) {
              fs::file_delete(DownFile)
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

              # Ensure Success is only TRUE if both the zip file exists and
              # passes integrity check
              Success <- file.exists(DownFile) && IASDT.R::CheckZip(DownFile)

            },
            error = function(e) {
              if (Attempt >= Attempts) {
                stop(
                  "Failed to download data after ", Attempts, " attempts: ",
                  conditionMessage(e), call. = FALSE)
              }
              Attempt <- Attempt + 1
            })
          }

          return(DownFile)

        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf,
          globals = "Path_Raw",
          packages = c("dplyr", "IASDT.R", "rgbif", "stringr", "fs", "withr"))))

  save(Efforts_AllRequests, file = Path_Efforts_Request)

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeDown,
    Prefix = "Downloading efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}
