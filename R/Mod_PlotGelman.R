## |------------------------------------------------------------------------| #
# PlotGelman ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks
#'
#' This function generates and optionally saves a series of plots showing the
#' evolution of Gelman and Rubin's shrink factor for various model parameters as
#' the number of iterations increases. It is designed to help assess the
#' convergence of Hmsc model by visualizing the shrink factor over iterations.
#' The function supports parallel processing and can handle multiple parameters
#' simultaneously.
#' @param InputCoda coda object or path to RData file containing the coda object.
#' @param Alpha Logical indicating whether to plot the Gelman-Rubin statistic
#'   for the Alpha parameter. If `TRUE` (default), the function executes the
#'   [IASDT.R::PlotGelman_Alpha] function.
#' @param Beta Logical indicating whether to plot the Gelman-Rubin statistic for
#'   the Beta parameter. If `TRUE` (default), the function executes the
#'   [IASDT.R::PlotGelman_Beta] function.
#' @param Omega Logical indicating whether to plot the Gelman-Rubin statistic
#'   for the Omega parameter. If `TRUE` (default), the function executes the
#'   [IASDT.R::PlotGelman_Omega] function.
#' @param Rho Logical indicating whether to plot the Gelman-Rubin statistic for
#'   the Rho parameter.  If `TRUE` (default), the function executes the
#'   [IASDT.R::PlotGelman_Rho] function.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param NOmega Integer specifying the number of species to be sampled for the
#'   Omega parameter.
#' @param PlottingAlpha A numeric value between 0 and 1 indicating the
#'   transparency level of the plot lines. Defaults to 0.25.
#' @param EnvFile String specifying the path to the environment variables file.
#'   The default value is ".env".
#' @param SavePlot Logical indicating whether to save the generated plots as a
#'   PDF file.
#' @param ReturnPlots String specifying the folder path where output figures
#'   should be saved. This parameter is mandatory if `SavePlot` is TRUE.
#' @name PlotGelman
#' @author Ahmed El-Gabbas
#' @return If `ReturnPlots` is `TRUE`, returns a list of ggplot objects
#'   corresponding to the generated plots. Otherwise, returns `NULL`. If
#'   `SavePlot` is `TRUE`, pdf file for the plots will be saved.
#' @export
#' @seealso
#' [IASDT.R::PlotGelman_Alpha]<br>[IASDT.R::PlotGelman_Beta]<br>
#' [IASDT.R::PlotGelman_Rho]<br>[IASDT.R::PlotGelman_Omega]

PlotGelman <- function(
    InputCoda = NULL, Alpha = TRUE, Beta = TRUE, Omega = TRUE, Rho = TRUE,
    NOmega = 1000, FromHPC = TRUE, PlottingAlpha = 0.25,
    EnvFile = ".env", SavePlot = TRUE, ReturnPlots = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments --------

  if (sum(Alpha, Beta, Omega, Rho) == 0) {
    stop("At least one of Alpha, Beta, Omega, and Rho must be `TRUE`",
         call. = FALSE)
  }

  if (is.null(InputCoda)) {
    stop("InputCoda cannot be empty", call. = FALSE)
  }

  if (isFALSE(SavePlot) && isFALSE(ReturnPlots)) {
    stop("At least one of SavePlot or ReturnPlots must be TRUE", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs, c("NOmega", "PlottingAlpha"), "numeric")
  IASDT.R::CheckArgs(
    AllArgs, c("Beta", "Rho", "Omega", "Alpha", "SavePlot", "ReturnPlots"),
    "logical")

  if (SavePlot) {
    if (Beta) IASDT.R::CheckArgs(AllArgs, "Beta_File", "character")
    if (Rho) IASDT.R::CheckArgs(AllArgs, "Rho_File", "character")
    if (Omega) IASDT.R::CheckArgs(AllArgs, "Omega_File", "character")
    if (Alpha) IASDT.R::CheckArgs(AllArgs, "Alpha_File", "character")
  }

  rm(AllArgs)

  # # ..................................................................... ###

  # Loading coda object ------

  if (inherits(InputCoda, "character")) {
    IASDT.R::CatTime("Loading coda object")
    CodaObj <- IASDT.R::LoadAs(InputCoda)
  } else {

    if (!inherits(InputCoda, "list")) {
      stop("`InputCoda` is neither character path or a list", call. = FALSE)
    }
    if (!inherits(InputCoda[[1]], "mcmc.list")) {
      stop("`InputCoda` has no mcmc.list items", call. = FALSE)
    }

    CodaObj <- InputCoda
    rm(InputCoda)
  }

  OutPath <- file.path(dirname(dirname(InputCoda)), "Model_Convergence")
  fs::dir_create(OutPath)

  # # ..................................................................... ###

  # Alpha -----

  if (Alpha) {
    IASDT.R::CatTime("Alpha")
    PlotObj_Alpha <- IASDT.R::PlotGelman_Alpha(
      CodaObj = CodaObj$Alpha[[1]], PlottingAlpha = PlottingAlpha)
  } else {
    PlotObj_Alpha <- NULL
  }

  # # ..................................................................... ###

  # Beta -----

  if (Beta) {
    IASDT.R::CatTime("Beta")
    PlotObj_Beta <- IASDT.R::PlotGelman_Beta(
      CodaObj = CodaObj$Beta, EnvFile = EnvFile,
      PlottingAlpha = PlottingAlpha, FromHPC = FromHPC)
  } else {
    PlotObj_Beta <- NULL
  }

  # # ..................................................................... ###

  # Omega -----

  if (Omega) {
    IASDT.R::CatTime("Omega")
    PlotObj_Omega <- IASDT.R::PlotGelman_Omega(
      CodaObj = CodaObj$Omega[[1]], NOmega = NOmega,
      PlottingAlpha = PlottingAlpha)
  } else {
    PlotObj_Omega <- NULL
  }

  # # ..................................................................... ###

  # Rho -----

  if (Rho && ("Rho" %in% names(CodaObj))) {
    IASDT.R::CatTime("Rho")
    PlotObj_Rho <- magrittr::extract2(CodaObj, "Rho") %>%
      IASDT.R::PlotGelman_Rho()
  } else {
    PlotObj_Rho <- NULL
  }

  # # ..................................................................... ###

  # Saving plots -----

  if (SavePlot) {
    PlotList <- list(
      Alpha = PlotObj_Alpha, Beta = PlotObj_Beta, Omega = PlotObj_Omega,
      Rho = PlotObj_Rho)

    PlotList4Plot <- purrr::list_flatten(purrr::discard(PlotList, is.null))

    if (length(PlotList4Plot) > 0) {
      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      grDevices::pdf(
        file = file.path(OutPath, "GelmanPlots.pdf"), width = 13, height = 7)
      print(gridExtra::marrangeGrob(
        grobs = PlotList4Plot, nrow = 1, ncol = 1, top = NULL))
      grDevices::dev.off()
    } else {
      warning("No plots to save")
    }

    IASDT.R::CatDiff(InitTime = .StartTime)

    if (ReturnPlots) {
      return(PlotList)
    } else {
      return(invisible(NULL))
    }
  } else {
    return(invisible(NULL))
  }
}
