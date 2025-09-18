## |------------------------------------------------------------------------| #
# mod_fit_windows ----
## |------------------------------------------------------------------------| #

#' Fit Hmsc-HPC models on UFZ Windows Server
#'
#' This function fits Hmsc models on a UFZ Windows Server. It reads model
#' configurations from a specified path, loads environment variables, checks
#' input arguments for validity, and executes model fitting in parallel if
#' required.
#' @param path_model Character. Path to the model files. This argument can not
#'   be empty.
#' @param python_ve Character. Path to a valid Python virtual environment.
#'   Defaults to `NULL`. This argument can not be empty.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @name mod_fit_windows
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but prints messages to the
#'   console regarding the progress and completion of model fitting.
#' @export

mod_fit_windows <- function(
    path_model = NULL, python_ve = NULL, n_cores = NULL,
    strategy = "multisession") {

  # # |||||||||||||||||||||||||||||||||||
  # # Check input arguments
  # # |||||||||||||||||||||||||||||||||||

  ecokit::check_args(
    args_to_check = c("path_model", "python_ve"), args_type = "character")
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # exit the function if not running on Windows
  if (ecokit::os() != "Windows") {
    ecokit::stop_ctx(
      "This function is only for Windows OS.", OS = ecokit::os(),
      include_backtrace = TRUE)
  }

  # Check if path_model is a valid directory
  if (!fs::dir_exists(path_model)) {
    ecokit::stop_ctx(
      "Model directory does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  # Check if python_ve is a valid directory
  if (!fs::dir_exists(python_ve)) {
    ecokit::stop_ctx(
      "Python virtual environment directory does not exist",
      python_ve = python_ve, include_backtrace = TRUE)
  }

  # Check if python_ve directory contains python virtual environment
  python_exe <- fs::path(python_ve, "Scripts", "python.exe")
  if (!fs::file_exists(python_exe)) {
    ecokit::stop_ctx(
      "Python virtual environment does not exist", python_ve = python_ve,
      include_backtrace = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_mod_progress <- command_ws <- path_hmsc_ws <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # List of models to be fitted
  # # |||||||||||||||||||||||||||||||||||

  model_to_run <- fs::path(path_model, "model_info.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(path_mod_progress, command_ws) %>%
    tidyr::unnest(cols = c("path_mod_progress", "command_ws")) %>%
    dplyr::filter(!file.exists(path_mod_progress))

  if (nrow(model_to_run) > 0) {
    ecokit::cat_time(
      paste0("There are ", nrow(model_to_run), " model variants to be fitted."))

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    run_commands <- future.apply::future_lapply(
      X = seq_len(nrow(model_to_run)),
      FUN = function(x) {
        system2(
          command = path_hmsc_ws, args = model_to_run$command_ws[x],
          stdout = model_to_run$path_mod_progress[x],
          stderr = model_to_run$path_mod_progress[x])
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c("path_hmsc_ws", "model_to_run"))

    rm(run_commands, envir = environment())

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("sequential", gc = TRUE)
    }

  } else {
    ecokit::cat_time("All model variants were already fitted.")
  }

  return(invisible(NULL))
}
