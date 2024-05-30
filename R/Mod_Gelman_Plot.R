## |------------------------------------------------------------------------| #
# Gelman_Plot ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks plots
#'
#' Plot Gelman-Rubin-Brooks plots
#'
#' @param CodaPath Path to RData file containing the coda object
#' @param Beta Logical. Run `IASDT.R::Gelman_Beta`?
#' @param Rho Logical. Run `IASDT.R::Gelman_Rho`?
#' @param Omega Logical. Run `IASDT.R::Gelman_Omega`?
#' @param Alpha Logical. Run `IASDT.R::Gelman_Alpha` (not yet available)?
#' @param NCores Integer. Number of parallel processes.
#' @param NOmega Integer. Number of species to be sampled for the Omega parameter
#' @param PlotAlpha Double. Plotting alpha for line transparency
#' @param SavePlot Logical. Save the outputs as JPEG file
#' @param OutPath String. Folder path to save the output figures
#' @param Beta_File String. File name (with extension) for saving the results of `IASDT.R::Gelman_Beta`
#' @param Rho_File String. File name (with extension) for saving the results of `IASDT.R::Gelman_Rho`
#' @param Omega_File String. File name (with extension) for saving the results of `IASDT.R::Gelman_Omega`
#' @param Alpha_File String. File name (with extension) for saving the results of `IASDT.R::Gelman_Alpha`
#' @name Gelman_Plot
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Gelman_Plot <- function(
    CodaPath = NULL, Beta = TRUE, Rho = TRUE, Omega = TRUE, Alpha = TRUE,
    NCores = NULL, NOmega = 1000, PlotAlpha = 0.25,
    SavePlot = NULL, OutPath = NULL, Beta_File = NULL, Rho_File = NULL,
    Omega_File = NULL, Alpha_File = NULL) {

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs, "CodaPath", "character")
  IASDT.R::CheckArgs(AllArgs, c("NCores", "NOmega", "PlotAlpha"), "numeric")
  IASDT.R::CheckArgs(
    AllArgs, c("Beta", "Rho", "Omega", "Alpha", "SavePlot"), "logical")

  if (SavePlot) {
    IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "OutPath", Type = "character")
    if (Beta) IASDT.R::CheckArgs(AllArgs, "Beta_File", "character")
    if (Rho) IASDT.R::CheckArgs(AllArgs, "Rho_File", "character")
    if (Omega) IASDT.R::CheckArgs(AllArgs, "Omega_File", "character")
    if (Alpha) IASDT.R::CheckArgs(AllArgs, "Alpha_File", "character")
    fs::dir_create(OutPath)
  }

  rm(AllArgs)


  IASDT.R::CatTime("Alpha")
  if (Alpha) {
    Plot_Alpha <- NULL
    # if (SavePlot) {
    #   ggplot2::ggsave(
    #     plot = Plot_Alpha, filename = file.path(OutPath, Alpha_File),
    #     width = 16, height = 9, dpi = 600)
    # }

  } else {
    Plot_Alpha <- NULL
  }

  IASDT.R::CatTime("Beta")
  if (Beta) {
    Plot_Beta <- CodaPath %>%
      IASDT.R::LoadAs() %>%
      magrittr::extract2("Beta") %>%
      IASDT.R::Gelman_Beta(NCores = NCores, PlotAlpha = PlotAlpha)

    if (SavePlot) {
      ggplot2::ggsave(
        plot = Plot_Beta, filename = file.path(OutPath, Beta_File),
        width = 16, height = 9, dpi = 600)
    }

    invisible(gc())
  } else {
    Plot_Beta <- NULL
  }

  IASDT.R::CatTime("Omega")
  if (Omega) {
    Plot_Omega <- CodaPath %>%
      IASDT.R::LoadAs() %>%
      magrittr::extract2("Omega") %>%
      magrittr::extract2(1) %>%
      IASDT.R::Gelman_Omega(NCores = NCores, PlotAlpha = PlotAlpha, NOmega = NOmega)

    if (SavePlot) {
      ggplot2::ggsave(
        plot = Plot_Omega, filename = file.path(OutPath, Omega_File),
        width = 16, height = 9, dpi = 600)
    }

    invisible(gc())
  } else {
    Plot_Omega <- NULL
  }

  IASDT.R::CatTime("Rho")
  if (Rho) {
    Plot_Rho <- CodaPath %>%
      IASDT.R::LoadAs() %>%
      magrittr::extract2("Rho") %>%
      IASDT.R::Gelman_Rho()

    if (SavePlot) {
      ggplot2::ggsave(
        plot = Plot_Rho, filename = file.path(OutPath, Rho_File),
        width = 16, height = 9, dpi = 600)
    }

    invisible(gc())
  } else {
    Plot_Rho <- NULL
  }

  return(list(
    Alpha = Plot_Alpha, Beta = Plot_Beta,
    Omega = Plot_Omega, Rho = Plot_Rho))
}
