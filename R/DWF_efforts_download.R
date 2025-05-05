## |------------------------------------------------------------------------| #
# efforts_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 3
#' @export

efforts_download <- function(n_cores = 6L, env_file = ".env") {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- Path_Efforts <- Path_Raw <- NULL

  # # ..................................................................... ###

  # Validate n_cores
  if (missing(n_cores) || !is.numeric(n_cores) || n_cores <= 0) {
    IASDT.R::stop_ctx("n_cores must be a positive integer.", n_cores = n_cores)
  }

  if (isFALSE(IASDT.R::check_system_command("unzip"))) {
    IASDT.R::stop_ctx("The 'unzip' command is not available")
  }

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE,
    "Path_Raw", "DP_R_Efforts_raw", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  Path_Efforts_Request <- fs::path(Path_Efforts, "Efforts_AllRequests.RData")

  if (!file.exists(Path_Efforts_Request)) {
    IASDT.R::stop_ctx(
      "The path for the `Efforts_AllRequests` data does not exist",
      Path_Efforts_Request = Path_Efforts_Request)
  }

  Efforts_AllRequests <- IASDT.R::load_as(Path_Efforts_Request)

  # # ..................................................................... ###

  ## Prepare working in parallel -----

  IASDT.R::cat_time(
    paste0("Prepare working in parallel using ", n_cores, " cores"),
    level = 1L)

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(n_cores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  IASDT.R::cat_time("Downloading & checking efforts data", level = 1L)

  Efforts_AllRequests <- Efforts_AllRequests %>%
    dplyr::mutate(
      # Download datasets in parallel
      DownPath = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          DownFile <- fs::path(Path_Raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not download it
          if (file.exists(DownFile)) {
            Success <- IASDT.R::check_zip(DownFile)
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
              Success <- file.exists(DownFile) && IASDT.R::check_zip(DownFile)

            },
            error = function(e) {
              if (Attempt >= Attempts) {
                IASDT.R::stop_ctx(
                  paste0(
                    "Failed to download data after ", Attempts, " attempts: ",
                    conditionMessage(e)))
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
  IASDT.R::cat_time("Stopping cluster", level = 1L)
  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTimeDown,
    prefix = "Downloading efforts data took ", level = 1L)

  # # ..................................................................... ###

  return(invisible(NULL))
}
