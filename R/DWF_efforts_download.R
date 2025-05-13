## |------------------------------------------------------------------------| #
# efforts_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 3
#' @export

efforts_download <- function(
    n_cores = 6L, strategy = "future::multicore", env_file = ".env") {

  .StartTimeDown <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- Path_Efforts <- NULL

  # # ..................................................................... ###

  # Validate n_cores
  if (missing(n_cores) || !is.numeric(n_cores) || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "future::sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c(
    "future::sequential", "future::multisession", "future::multicore",
    "future::cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  if (isFALSE(ecokit::check_system_command("unzip"))) {
    ecokit::stop_ctx(
      "The 'unzip' command is not available", include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Efforts", "DP_R_Efforts_processed", FALSE, FALSE,
    "Path_Raw", "DP_R_Efforts_raw", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  Path_Efforts_Request <- fs::path(Path_Efforts, "Efforts_AllRequests.RData")

  if (!file.exists(Path_Efforts_Request)) {
    ecokit::stop_ctx(
      "The path for the `Efforts_AllRequests` data does not exist",
      Path_Efforts_Request = Path_Efforts_Request, include_backtrace = TRUE)
  }

  Efforts_AllRequests <- ecokit::load_as(Path_Efforts_Request)

  # # ..................................................................... ###

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  ecokit::cat_time("Downloading & checking efforts data", level = 1L)

  if (strategy == "future::multicore") {
    pkg_to_export <- NULL
  } else {
    pkg_to_export <- c("dplyr", "ecokit", "rgbif", "stringr", "fs", "withr")
  }

  Efforts_AllRequests <- Efforts_AllRequests %>%
    dplyr::mutate(
      # Download datasets in parallel
      DownPath = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          DownFile <- fs::path(Path_Raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not download it
          if (file.exists(DownFile)) {
            Success <- ecokit::check_zip(DownFile)
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
              Success <- file.exists(DownFile) && ecokit::check_zip(DownFile)

            },
            error = function(e) {
              if (Attempt >= Attempts) {
                ecokit::stop_ctx(
                  paste0(
                    "Failed to download data after ", Attempts, " attempts: ",
                    conditionMessage(e)),
                  include_backtrace = TRUE)
              }
              Attempt <- Attempt + 1
            })
          }

          return(DownFile)

        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf,
          globals = "Path_Raw", packages = pkg_to_export)))

  save(Efforts_AllRequests, file = Path_Efforts_Request)

  # # ..................................................................... ###

  # Stopping cluster ------
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .StartTimeDown,
    prefix = "Downloading efforts data took ", level = 1L)

  # # ..................................................................... ###

  return(invisible(NULL))
}
