## |------------------------------------------------------------------------| #
# Mod_Fit_WS ----
## |------------------------------------------------------------------------| #

#' Fit Hmsc-HPC models on UFZ Windows Server
#'
#' This function fits Hmsc models on a UFZ Windows Server. It reads model
#' configurations from a specified path, loads environment variables, checks
#' input arguments for validity, and executes model fitting in parallel if
#' required.
#' @param Path_Model Character. Path to the model files.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#' @name Mod_Fit_WS
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but prints messages to the
#'   console regarding the progress and completion of model fitting.
#' @export

Mod_Fit_WS <- function(Path_Model, EnvFile = ".env", NCores = NULL) {

  if (is.null(Path_Model) || is.null(NCores)) {
    stop("Path_Model and NCores cannot be empty", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_ModProg <- Command_WS <- Path_Hmsc_WS <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (!file.exists(EnvFile)) {
    stop(
      "Path to environment variables: ", EnvFile, " was not found",
      call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Hmsc_WS", "DP_R_Hmsc_ve_win", TRUE, FALSE)

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Check input arguments
  # # |||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Path_Hmsc_WS", "Path_Model"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "NCores", Type = "numeric")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # List of models to be fitted
  # # |||||||||||||||||||||||||||||||||||

  if (!fs::dir_exists(Path_Model)) {
    stop("Model directory does not exist:", Path_Model, call. = FALSE)
  }

  Model2Run <- IASDT.R::Path(Path_Model, "Model_Info.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(Path_ModProg, Command_WS) %>%
    tidyr::unnest(cols = c("Path_ModProg", "Command_WS")) %>%
    dplyr::filter(!file.exists(Path_ModProg))

  if (nrow(Model2Run) > 0) {
    IASDT.R::CatTime(
      paste0("There are ", nrow(Model2Run), " model variants to be fitted."))

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

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

  } else {
    IASDT.R::CatTime("All model variants were already fitted.")
  }

  return(invisible(NULL))
}
