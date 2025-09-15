## |------------------------------------------------------------------------| #
# efforts_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 3
#' @export

efforts_download <- function(
    n_cores = 6L, strategy = "multisession", env_file = ".env") {

  .start_time_down <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Request <- path_efforts <- NULL

  # # ..................................................................... ###

  # Validate n_cores
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  if (isFALSE(ecokit::check_system_command("unzip"))) {
    ecokit::stop_ctx(
      "The 'unzip' command is not available", include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_efforts", "DP_R_efforts_processed", FALSE, FALSE,
    "path_raw", "DP_R_efforts_raw", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("dplyr", "ecokit", "rgbif", "stringr", "fs", "withr"),
    strategy = strategy)

  # # ..................................................................... ###

  path_efforts_request <- fs::path(path_efforts, "efforts_all_requests.RData")

  if (!file.exists(path_efforts_request)) {
    ecokit::stop_ctx(
      "The path for the `efforts_all_requests` data does not exist",
      path_efforts_request = path_efforts_request, include_backtrace = TRUE)
  }

  efforts_all_requests <- ecokit::load_as(path_efforts_request)

  # # ..................................................................... ###

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  # Downloading/checking efforts data ------
  ecokit::cat_time("Downloading & checking efforts data", level = 1L)

  efforts_all_requests <- efforts_all_requests %>%
    dplyr::mutate(
      # Download datasets in parallel
      download_path = furrr::future_map_chr(
        .x = Request,
        .f = ~{

          down_file <- fs::path(path_raw, paste0(as.character(.x), ".zip"))

          # Check zip file if exist, if not download it
          if (file.exists(down_file)) {
            success <- ecokit::check_zip(down_file)
            if (isFALSE(success)) {
              fs::file_delete(down_file)
            }
          } else {
            success <- FALSE
          }

          # Try downloading data for a max of 3 attempts, each with 20 mins
          # time out
          withr::local_options(list(timeout = 1200))

          attempt <- 1
          attempts <- 3

          while (isFALSE(success) && (attempt <= attempts)) {
            tryCatch({
              suppressMessages(
                rgbif::occ_download_get(
                  key = .x, path = path_raw, overwrite = TRUE))

              # Ensure success is only TRUE if both the zip file exists and
              # passes integrity check
              success <- file.exists(down_file) && ecokit::check_zip(down_file)

            },
            error = function(e) {
              if (attempt >= attempts) {
                ecokit::stop_ctx(
                  paste0(
                    "Failed to download data after ", attempts, " attempts: ",
                    conditionMessage(e)),
                  include_backtrace = TRUE)
              }
              attempt <- attempt + 1
            })
          }

          return(down_file)

        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = Inf,
          globals = "path_raw", packages = pkg_to_export)))

  save(efforts_all_requests, file = path_efforts_request)

  # # ..................................................................... ###

  # Stopping cluster ------
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time_down,
    prefix = "Downloading efforts data took ", level = 1L)

  # # ..................................................................... ###

  return(invisible(NULL))
}
