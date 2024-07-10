## |------------------------------------------------------------------------| #
# PlotGelman ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks
#'
#' Plot Gelman-Rubin-Brooks. This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases. For more information, see `coda:::gelman.plot`
#'
#' @param InputCoda Path to RData file containing the coda object
#' @param Beta Logical. Run `IASDT.R::PlotGelman_Beta`?
#' @param Rho Logical. Run `IASDT.R::PlotGelman_Rho`?
#' @param Omega Logical. Run `IASDT.R::PlotGelman_Omega`?
#' @param Alpha Logical. Run `IASDT.R::PlotGelman_Alpha`
#' @param NCores Integer. Number of parallel processes.
#' @param NOmega Integer. Number of species to be sampled for the Omega parameter
#' @param PlottingAlpha Double. Plotting alpha for line transparency
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param SavePlot Logical. Save the outputs as PDF
#' @param ReturnPlots Logical. Return ggplot objects
#' @param OutPath String. Folder path to save the output figures
#' @name PlotGelman
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotGelman <- function(
    InputCoda = NULL, Beta = TRUE, Rho = TRUE, Omega = TRUE, Alpha = TRUE,
    NCores = NULL, NOmega = 1000, PlottingAlpha = 0.25, EnvFile = ".env",
    SavePlot = TRUE, ReturnPlots = FALSE, OutPath = NULL) {

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Checking arguments
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (SavePlot == FALSE && ReturnPlots == FALSE) {
    MSG <- "SavePlot & ReturnPlots can not be both FALSE"
    stop(MSG)
  }

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs, c("NCores", "NOmega", "PlottingAlpha"), "numeric")
  IASDT.R::CheckArgs(
    AllArgs, c("Beta", "Rho", "Omega", "Alpha", "SavePlot", "ReturnPlots"), "logical")

  if (SavePlot) {
    IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "OutPath", Type = "character")
    if (Beta) IASDT.R::CheckArgs(AllArgs, "Beta_File", "character")
    if (Rho) IASDT.R::CheckArgs(AllArgs, "Rho_File", "character")
    if (Omega) IASDT.R::CheckArgs(AllArgs, "Omega_File", "character")
    if (Alpha) IASDT.R::CheckArgs(AllArgs, "Alpha_File", "character")
    fs::dir_create(OutPath)
  }

  rm(AllArgs)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading coda object ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (inherits(InputCoda, "character")) {
    CodaObj <- IASDT.R::LoadAs(InputCoda)
  } else {
    CodaObj <- InputCoda
    rm(InputCoda)
    invisible(gc())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Alpha -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Alpha")
  if (Alpha) {
    PlotObj_Alpha <- CodaObj %>%
      magrittr::extract2("Alpha") %>%
      magrittr::extract2(1)
    PlotObj_Alpha <- IASDT.R::PlotGelman_Alpha(
      CodaObj = PlotObj_Alpha, NCores = NCores, PlottingAlpha = PlottingAlpha)
  } else {
    PlotObj_Alpha <- NULL
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Beta")
  if (Beta) {
    PlotObj_Beta <- IASDT.R::PlotGelman_Beta(
      CodaObj = magrittr::extract2(CodaObj, "Beta"),
      NCores = NCores, PlottingAlpha = PlottingAlpha, EnvFile = EnvFile)
    invisible(gc())
  } else {
    PlotObj_Beta <- NULL
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Omega -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Omega")
  if (Omega) {
    PlotObj_Omega <- CodaObj %>%
      magrittr::extract2("Omega") %>%
      magrittr::extract2(1)
    PlotObj_Omega <- IASDT.R::PlotGelman_Omega(
      CodaObj = PlotObj_Omega, NCores = NCores,
      PlottingAlpha = PlottingAlpha, NOmega = NOmega)
    invisible(gc())
  } else {
    PlotObj_Omega <- NULL
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Rho -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Rho")
  if (Rho) {
    PlotObj_Rho <- IASDT.R::PlotGelman_Rho(
      CodaObj = magrittr::extract2(CodaObj, "Rho"))
    invisible(gc())
  } else {
    PlotObj_Rho <- NULL
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Saving plots -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  PlotList <- list(
    Alpha = PlotObj_Alpha, Beta = PlotObj_Beta, Omega = PlotObj_Omega,
    Rho = PlotObj_Rho)

  if (SavePlot) {
    PlotList4Plot <- purrr::list_flatten(purrr::discard(PlotList, is.null))
    ggplot2::ggsave(
      filename = file.path(OutPath, "GelmanPlots.pdf"),
      plot = gridExtra::marrangeGrob(
        grobs = PlotList4Plot, nrow = 1, ncol = 1, top = NULL),
      width = 13, height = 7)
  }

  if (ReturnPlots) {
    return(PlotList)
  } else {
    return(invisible(NULL))
  }
}
