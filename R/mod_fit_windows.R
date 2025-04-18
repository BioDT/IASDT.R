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
  if (IASDT.R::OS() != "Windows") {
    stop("This function is only for Windows OS.", call. = FALSE)
  }

  if (is.null(path_model) || is.null(python_VE) || is.null(n_cores)) {
    stop(
      "`path_model`, `python_VE`, and `n_cores` cannot be empty",
      call. = FALSE)
  }

  # Check if path_model is a valid directory
  if (!fs::dir_exists(path_model)) {
    stop("Model directory does not exist:", path_model, call. = FALSE)
  }

  # Check if python_VE is a valid directory
  if (!fs::dir_exists(python_VE)) {
    stop(
      "Python virtual environment directory does not exist:", python_VE,
      call. = FALSE)
  }

  # Check if python_VE directory contains python virtual environment
  python_exe <- IASDT.R::path(python_VE, "Scripts", "python.exe")
  if (!fs::file_exists(python_exe)) {
    stop(
      "Python virtual environment does not exist:", python_VE, call. = FALSE)
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

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("Path_Hmsc_WS", "path_model"))

  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = "n_cores", args_type = "numeric")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # List of models to be fitted
  # # |||||||||||||||||||||||||||||||||||

  Model2Run <- IASDT.R::path(path_model, "Model_Info.RData") %>%
    IASDT.R::load_as() %>%
    dplyr::select(Path_ModProg, Command_WS) %>%
    tidyr::unnest(cols = c("Path_ModProg", "Command_WS")) %>%
    dplyr::filter(!file.exists(Path_ModProg))

  if (nrow(Model2Run) > 0) {
    IASDT.R::cat_time(
      paste0("There are ", nrow(Model2Run), " model variants to be fitted."))

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
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

  } else {
    IASDT.R::cat_time("All model variants were already fitted.")
  }

  return(invisible(NULL))
}
