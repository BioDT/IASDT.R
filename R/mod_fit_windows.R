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
#' @param python_VE Character. Path to a valid Python virtual environment.
#'   Defaults to `NULL`. This argument can not be empty.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#' @name mod_fit_windows
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but prints messages to the
#'   console regarding the progress and completion of model fitting.
#' @export

mod_fit_windows <- function(
    path_model = NULL, python_VE = NULL, n_cores = NULL) {

  # exit the function if not running on Windows
  if (ecokit::os() != "Windows") {
    ecokit::stop_ctx(
      "This function is only for Windows OS.", OS = ecokit::os(),
      include_backtrace = TRUE)
  }

  if (is.null(path_model) || is.null(python_VE) || is.null(n_cores)) {
    ecokit::stop_ctx(
      "`path_model`, `python_VE`, and `n_cores` cannot be empty",
      path_model = path_model, python_VE = python_VE, n_cores = n_cores,
      include_backtrace = TRUE)
  }

  # Check if path_model is a valid directory
  if (!fs::dir_exists(path_model)) {
    ecokit::stop_ctx(
      "Model directory does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  # Check if python_VE is a valid directory
  if (!fs::dir_exists(python_VE)) {
    ecokit::stop_ctx(
      "Python virtual environment directory does not exist",
      python_VE = python_VE, include_backtrace = TRUE)
  }

  # Check if python_VE directory contains python virtual environment
  python_exe <- fs::path(python_VE, "Scripts", "python.exe")
  if (!fs::file_exists(python_exe)) {
    ecokit::stop_ctx(
      "Python virtual environment does not exist", python_VE = python_VE,
      include_backtrace = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_ModProg <- Command_WS <- Path_Hmsc_WS <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Check input arguments
  # # |||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("Path_Hmsc_WS", "path_model"))

  ecokit::check_args(
    args_all = AllArgs, args_to_check = "n_cores", args_type = "numeric")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # List of models to be fitted
  # # |||||||||||||||||||||||||||||||||||

  Model2Run <- fs::path(path_model, "Model_Info.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(Path_ModProg, Command_WS) %>%
    tidyr::unnest(cols = c("Path_ModProg", "Command_WS")) %>%
    dplyr::filter(!file.exists(Path_ModProg))

  if (nrow(Model2Run) > 0) {
    ecokit::cat_time(
      paste0("There are ", nrow(Model2Run), " model variants to be fitted."))

    if (n_cores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = "future::multicore")
      withr::defer(future::plan("future::sequential", gc = TRUE))
    }

    RunCommands <- future.apply::future_lapply(
      X = seq_len(nrow(Model2Run)),
      FUN = function(x) {
        system2(
          command = Path_Hmsc_WS, args = Model2Run$Command_WS[x],
          stdout = Model2Run$Path_ModProg[x],
          stderr = Model2Run$Path_ModProg[x])
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c("Path_Hmsc_WS", "Model2Run"))

    rm(RunCommands, envir = environment())

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("future::sequential", gc = TRUE)
    }

  } else {
    ecokit::cat_time("All model variants were already fitted.")
  }

  return(invisible(NULL))
}
