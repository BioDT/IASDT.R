## |------------------------------------------------------------------------| #
# Mod_Fit_WS ----
## |------------------------------------------------------------------------| #

#' Fit Hmsc-HPC models on UFZ Windows Server
#'
#' Fit Hmsc-HPC models on UFZ Windows Server
#'
#' @param Path_Model String. Path to the model files (without trailing slash)
#' @param Path_EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param NCores Integer. Number of parallel processes.
#' @name Mod_Fit_WS
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_Fit_WS <- function(
    Path_Model = NULL, Path_EnvFile = ".env", NCores = NULL) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_ModPorg <- Command_WS <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (file.exists(Path_EnvFile)) {
    readRenviron(Path_EnvFile)
    Path_Hmsc_WS <- Sys.getenv("DP_R_Mod_Path_VE_WS")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", Path_EnvFile, " was not found")
    stop(MSG)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # CHECK input arguments
  # # |||||||||||||||||||||||||||||||||||

  AllArgs <- ls()
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

  Model2Run <- file.path(Path_Model, "Model_Info.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(Path_ModPorg, Command_WS) %>%
    tidyr::unnest(cols = c("Path_ModPorg", "Command_WS")) %>%
    dplyr::filter(magrittr::not(file.exists(Path_ModPorg)))

  if (nrow(Model2Run) > 0) {
    IASDT.R::CatTime(
      paste0("There are ", nrow(Model2Run), " model variants to be fitted."))

    c1 <- snow::makeSOCKcluster(NCores)
    future::plan(future::cluster, workers = c1, gc = TRUE)
    snow::clusterExport(
      cl = c1, list = c("Path_Hmsc_WS", "Model2Run"), envir = environment())

    RunCommands <- future.apply::future_lapply(
      X = seq_len(nrow(Model2Run)),
      FUN = function(x) {
        system2(
          command = Path_Hmsc_WS, args = Model2Run$Command_WS[x],
          stdout = Model2Run$Path_ModPorg[x],
          stderr = Model2Run$Path_ModPorg[x])
      },
      future.scheduling = Inf, future.seed = TRUE)
  } else {
    IASDT.R::CatTime("All model variants were already fitted.")
  }

  return(invisible(NULL))
}
