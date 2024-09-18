## |------------------------------------------------------------------------| #
# PlotConvergence ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of a selected model
#'
#' This function generates and saves plots for model convergence diagnostics,
#' including rho, alpha, omega, and beta parameters. It supports parallel
#' processing for faster execution and can work with models fitted on
#' High-Performance Computing HPC environments.
#' @param Path_Coda String. Path to the coda object containing MCMC samples.
#' @param Path_FittedModel String. Path to the fitted Hmsc model object.
#' @param EnvFile String. Path to the environment file containing necessary
#'   environment variables. Default: ".env".
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param NChains Integer. Number of MCMC chains used in the model. Default: 4.
#' @param Title String. Title for the rho and alpha plots. Default: " ".
#' @param NOmega Integer. Number of species interactions to sample for omega
#'   parameter convergence diagnostics. Default: 1000.
#' @param NCores Integer. Number of cores to use for parallel processing.
#' @param NRC Numeric Vector. Specifies the number of rows and columns for
#'   arranging multiple plots on a page. Default: c(2, 3).
#' @param SavePlotData Logical. Indicates whether to save the plot data as RData
#'   files. Default: TRUE.
#' @param PagePerFile Integer. Indicates the number of pages per single pdf page
#'   for the `Omega` parameter. Default: 20.
#' @param Cols Character Vector. Specifies the colors to use for each chain in
#'   the plots. Default: c("red", "blue", "darkgreen", "darkgrey").
#' @name PlotConvergence
#' @author Ahmed El-Gabbas
#' @return The function does not return a value but generates and saves plots to
#'   disk.
#' @export

PlotConvergence <- function(
    Path_Coda = NULL, Path_FittedModel = NULL, EnvFile = ".env",
    FromHPC = TRUE, NChains = 4, Title = " ", NOmega = 1000, NCores = NULL,
    NRC = c(2, 3), SavePlotData = TRUE, PagePerFile = 20,
    Cols = c("red", "blue", "darkgreen", "darkgrey")) {


  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Coda) || is.null(Path_FittedModel) || is.null(NCores)) {
    stop(
      "Path_Coda, Path_FittedModel, and NCores cannot be empty", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpComb <- `2.5%` <- `97.5%` <- Class <- Order <- Family <- Plot <- DT <-
    IAS_ID <- Species <- Variable <- data <- PlotID <- Var <- PlotFixedY <-
    File <- Page <- Iter <- Value <- Chain <- y <- label <-
    Var_Sp <- NULL

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Check input arguments ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Path_Coda", "Path_FittedModel"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NChains", "NOmega", "NCores", "NRC"))
  rm(AllArgs)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Create path ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Create path")
  Path_Convergence <- dirname(dirname(Path_Coda)) %>%
    file.path("Model_Convergence")
  Path_Convergence_BySp <- file.path(Path_Convergence, "Beta_BySpecies")
  fs::dir_create(c(Path_Convergence, Path_Convergence_BySp))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Prepare convergence data ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare convergence data")

  IASDT.R::CatTime("Loading coda object", Level = 1)
  Coda_Obj <- IASDT.R::LoadAs(Path_Coda)

  IASDT.R::CatTime("Loading fitted model object", Level = 1)
  Model <- IASDT.R::LoadAs(Path_FittedModel)

  ModVars <- Model$covNames

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Rho ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if ("Rho" %in% names(Coda_Obj)) {
    IASDT.R::CatTime("Rho")

    FileConv_Rho <- file.path(Path_Convergence, "Convergence_Rho.RData")

    if (file.exists(FileConv_Rho)) {
      IASDT.R::CatTime("Loading plotting data", Level = 1)
      PlotObj_Rho <- IASDT.R::LoadAs(FileConv_Rho)
    } else {
      IASDT.R::CatTime("Prepare plot", Level = 1)
      PlotObj_Rho <- IASDT.R::PlotRho(
        Post = Coda_Obj, Model = Model, Title = Title, Cols = Cols)

      if (SavePlotData) {
        IASDT.R::CatTime("Save plotting data", Level = 1)
        IASDT.R::SaveAs(
          InObj = PlotObj_Rho, OutObj = "Convergence_Rho",
          OutPath = FileConv_Rho)
      }
    }

    IASDT.R::CatTime("Save plot", Level = 1)
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::pdf(
      file = file.path(Path_Convergence, "Convergence_Rho.pdf"),
      width = 18, height = 12)
    print(PlotObj_Rho)
    grDevices::dev.off()

    rm(PlotObj_Rho)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Alpha  ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Alpha")

  FileConv_Alpha <- file.path(Path_Convergence, "Convergence_Alpha.RData")

  if (file.exists(FileConv_Alpha)) {
    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Alpha <- IASDT.R::LoadAs(FileConv_Alpha)
  } else {
    IASDT.R::CatTime("Prepare plot", Level = 1)
    PlotObj_Alpha <- IASDT.R::PlotAlpha(
      Post = Coda_Obj, Model = Model, Title = Title, NRC = NRC,
      AddFooter = FALSE, AddTitle = FALSE, Cols = Cols)

    if (SavePlotData) {
      IASDT.R::CatTime("Save plotting data", Level = 1)
      IASDT.R::SaveAs(
        InObj = PlotObj_Alpha, OutObj = "Convergence_Alpha",
        OutPath = file.path(Path_Convergence, "Convergence_Alpha.RData"))
    }
  }

  IASDT.R::CatTime("Save plot", Level = 1)
  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::pdf(
    file = file.path(Path_Convergence, "Convergence_Alpha.pdf"),
    width = 18, height = 10)
  print(PlotObj_Alpha)
  grDevices::dev.off()

  Obj_Omega <- Coda_Obj$Omega[[1]]
  Obj_Beta <- Coda_Obj$Beta

  rm(Model, Coda_Obj, PlotObj_Alpha)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Omega  ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Omega")

  FileConv_Omega <- file.path(Path_Convergence, "Convergence_Omega.RData")

  if (file.exists(FileConv_Omega)) {
    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Omega <- IASDT.R::LoadAs(FileConv_Omega)
  } else {

    IASDT.R::CatTime("Coda to tibble", Level = 1)
    OmegaDF <- IASDT.R::Coda_to_tibble(
      CodaObj = Obj_Omega, Type = "omega", NOmega = NOmega,
      EnvFile = EnvFile, FromHPC = FromHPC)

    SelectedCombs <- unique(OmegaDF$SpComb)

    IASDT.R::CatTime("Prepare confidence interval data", Level = 1)
    CI <- purrr::map(.x  = Obj_Omega, .f = ~ .x[, SelectedCombs]) %>%
      coda::mcmc.list() %>%
      summary(quantiles = c(0.025, 0.975)) %>%
      magrittr::extract2("quantiles") %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "SpComb")

    IASDT.R::CatTime("Prepare working in parallel", Level = 1)
    if (NCores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
    }

    IASDT.R::CatTime("Prepare trace plots data", Level = 1)
    PlotObj_Omega <- future.apply::future_lapply(
      X = seq_len(NOmega),
      FUN = function(x) {

        CI0 <- dplyr::filter(CI, SpComb == SelectedCombs[x]) %>%
          dplyr::select(-SpComb) %>%
          round(2) %>%
          unlist()

        CI1 <- paste0(CI0, collapse = "-") %>%
          paste0("<i><b>95% credible interval:</b></i> ", .) %>%
          data.frame(x = -Inf, y = -Inf, label = .)

        PanelTitle <- data.frame(
          x = Inf, y = Inf,
          label = {
            sort(c(OmegaDF$IAS1[[x]], OmegaDF$IAS2[[x]])) %>%
              paste0("<i>", ., "</i>") %>%
              paste0(collapse = " & <br>") %>%
              paste0(., " ")
          })

        OmegaDF2 <- dplyr::filter(OmegaDF, SpComb == SelectedCombs[x]) %>%
          dplyr::pull("DT") %>%
          magrittr::extract2(1)

        ValRange <- range(OmegaDF2$Value)

        Plot <- ggplot2::ggplot(
          data = OmegaDF2,
          mapping = ggplot2::aes(
            x = Iter, y = Value, color = factor(Chain))) +
          ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
          ggplot2::geom_smooth(
            method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
          ggplot2::geom_point(alpha = 0) +
          ggplot2::geom_hline(
            yintercept = CI0, linetype = "dashed", color = "black",
            linewidth = 1) +
          ggplot2::scale_color_manual(values = Cols) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label), size = 4,
            data = CI1, inherit.aes = FALSE, hjust = 0, vjust = 0,
            lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label), size = 5,
            data = PanelTitle, inherit.aes = FALSE, colour = "blue",
            hjust = 1, vjust = 1, lineheight = 0, fill = NA,
            label.color = NA) +
          ggplot2::theme_bw() +
          ggplot2::xlab(NULL) +
          ggplot2::ylab(NULL) +
          ggplot2::theme(
            axis.text = ggplot2::element_text(size = 12),
            legend.position = "none")

        if (dplyr::between(0, ValRange[1], ValRange[2])) {
          Plot <- Plot +
            ggplot2::geom_hline(
              yintercept = 0, linetype = "solid",
              color = "green", linewidth = 0.6)
        }

        Plot <- ggExtra::ggMarginal(
          p = Plot, type = "density", margins = "y", size = 5,
          color = "steelblue4")

        return(Plot)
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c("NOmega", "CI", "SelectedCombs", "OmegaDF", "Cols"),
      future.packages = c("dplyr", "coda", "ggplot2", "ggExtra", "ggtext"))

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("sequential", gc = TRUE)
    }


    if (SavePlotData) {
      IASDT.R::CatTime("Save plot data", Level = 1)
      IASDT.R::SaveAs(
        InObj = PlotObj_Omega, OutObj = "Convergence_Omega",
        OutPath = file.path(Path_Convergence, "Convergence_Omega.RData"))
    }
    rm(Obj_Omega, OmegaDF, SelectedCombs, CI)
  }

  OmegaPlotList <- tibble::tibble(PlotID = seq_len(length(PlotObj_Omega))) %>%
    dplyr::mutate(
      File = ceiling(PlotID / (PagePerFile * NRC[2] * NRC[1])),
      Page = ceiling(PlotID / (NRC[2] * NRC[1]))) %>%
    tidyr::nest(.by = c("File", "Page"), .key = "PlotID") %>%
    dplyr::mutate(
      PlotID = purrr::pmap(
        .l = list(File, Page, PlotID),
        .f = function(File, Page, PlotID) {

          PlotTitle <- ggplot2::ggplot() +
            ggplot2::labs(
              title = paste0(
                "Convergence of the omega parameter - a sample of ",
                NOmega, " species pair"),
              subtitle = paste0(
                "   File ", File, " | Page ",
                (Page - ((File - 1) * PagePerFile)))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              plot.title = ggplot2::element_text(
                face = "bold", size = 20, hjust = 0.5),
              plot.subtitle = ggplot2::element_text(
                size = 12, colour = "grey",
                margin = ggplot2::margin(-5, 0, 0, 0)))

          PlotObj_Omega[unlist(as.vector(PlotID))] %>%
            cowplot::plot_grid(
              plotlist = ., ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
            cowplot::plot_grid(
              PlotTitle, ., ncol = 1, rel_heights = c(0.05, 1))
        }))

  rm(PlotObj_Omega)

  purrr::walk(
    .x = seq_along(unique(OmegaPlotList$File)),
    .f = ~{
      invisible({
        CurrPlotOrder <- dplyr::filter(OmegaPlotList, File == .x)
        grDevices::pdf(
          file = file.path(
            Path_Convergence, paste0("Convergence_Omega_", .x, ".pdf")),
          width = 18, height = 12)
        purrr::map(CurrPlotOrder$PlotID, grid::grid.draw, recording = FALSE)
        grDevices::dev.off()
      })
    })

  rm(OmegaPlotList, Obj_Omega)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta - 1. Prepare data ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Beta")

  FileConv_Beta <- file.path(Path_Convergence, "Convergence_Beta.RData")

  if (file.exists(FileConv_Beta)) {
    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Beta <- IASDT.R::LoadAs(FileConv_Beta)

  } else {

    IASDT.R::CatTime("Prepare trace plot data", Level = 1)

    IASDT.R::CatTime(
      paste0("Prepare working in parallel using `", NCores, "`"), Level = 2)
    if (NCores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
    }


    IASDT.R::CatTime("Prepare 95% credible interval data", Level = 2)
    CI <- summary(Obj_Beta, quantiles = c(0.025, 0.975))$quantiles %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "Var_Sp") %>%
      dplyr::rename(CI_025 = `2.5%`, CI_975 = `97.5%`)

    IASDT.R::CatTime("Coda to tibble", Level = 2)
    PlotObj_Beta <- IASDT.R::Coda_to_tibble(
      CodaObj = Obj_Beta, Type = "beta", EnvFile = EnvFile,
      FromHPC = FromHPC) %>%
      dplyr::left_join(CI, by = "Var_Sp")
    rm(CI)

    VarRanges <- dplyr::arrange(PlotObj_Beta, Variable, IAS_ID) %>%
      dplyr::select(Var = Variable, DT) %>%
      dplyr::mutate(
        Range = purrr::map(.x = DT, .f = ~ range(dplyr::pull(.x, Value)))) %>%
      dplyr::select(-DT) %>%
      tidyr::nest(data = c("Range")) %>%
      dplyr::mutate(
        Range = purrr::map(
          .x = data,
          .f = ~{
            dplyr::pull(.x, Range) %>%
              as.vector() %>%
              range() %>%
              purrr::set_names(c("Min", "Max"))
          })) %>%
      dplyr::select(-data) %>%
      tidyr::unnest_wider("Range")

    IASDT.R::CatTime("Prepare data in parallel", Level = 2)
    PlotObj_Beta <- future.apply::future_lapply(
      X = PlotObj_Beta$Var_Sp,
      FUN = function(x) {

        SubDT <- dplyr::filter(PlotObj_Beta, Var_Sp == x)
        Species <- SubDT$Species
        IAS_ID <- SubDT$IAS_ID
        DT <- SubDT$DT[[1]]

        Variable <- SubDT$Variable
        ValRange <- range(DT$Value)

        ValRangeByVar <- dplyr::filter(VarRanges, Var == Variable) %>%
          dplyr::select(-Var) %>%
          unlist()

        CI0 <- round(c(SubDT$CI_025, SubDT$CI_975), 2)
        CI1 <- paste0(CI0, collapse = "  -  ") %>%
          paste0("<i><b>95% credible interval:</b><br></i>", .) %>%
          data.frame(x = -Inf, y = -Inf, label = .)

        PanelTitle <- data.frame(
          x = Inf, y = Inf,
          label = paste0("<br><b><i>", Species, "</i></b>"))

        PanelTitle2 <- IASDT.R::GetSpeciesName(
          SpID = IAS_ID, EnvFile = EnvFile, FromHPC = FromHPC) %>%
          dplyr::select(Class, Order, Family) %>%
          unlist() %>%
          paste0(collapse = " | ") %>%
          paste0("<b>", ., "</b>") %>%
          paste0("<br>", IAS_ID) %>%
          data.frame(x = -Inf, y = Inf, label = .)

        Plot <- ggplot2::ggplot(
          data = DT,
          mapping = ggplot2::aes(x = Iter, y = Value, color = factor(Chain))) +
          ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
          ggplot2::geom_smooth(
            method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
          ggplot2::geom_point(alpha = 0) +
          ggplot2::geom_hline(
            yintercept = CI0, linetype = "dashed", color = "black",
            linewidth = 1) +
          ggplot2::scale_color_manual(values = Cols) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(limits = ValRange) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label), data = CI1,
            inherit.aes = FALSE, hjust = 0, vjust = 0, lineheight = 0,
            fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = PanelTitle, inherit.aes = FALSE, colour = "blue", hjust = 1,
            vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = PanelTitle2, inherit.aes = FALSE, hjust = 0, vjust = 1,
            lineheight = 0, fill = NA, label.color = NA) +
          ggplot2::theme_bw() +
          ggplot2::xlab(NULL) +
          ggplot2::ylab(NULL) +
          ggplot2::theme(
            legend.position = "none",
            axis.text = ggplot2::element_text(size = 12))

        suppressMessages({
          Plot2 <- Plot +
            ggplot2::scale_y_continuous(limits = ValRangeByVar)
        })

        if (dplyr::between(0, ValRangeByVar[1], ValRangeByVar[2])) {
          Plot2 <- Plot2 +
            ggplot2::geom_hline(
              yintercept = 0, linetype = "solid",
              color = "green", linewidth = 0.6)
        }

        if (dplyr::between(0, ValRange[1], ValRange[2])) {
          Plot <- Plot +
            ggplot2::geom_hline(
              yintercept = 0, linetype = "solid",
              color = "green", linewidth = 0.6)
        }

        Plot <- ggExtra::ggMarginal(
          p = Plot, type = "density", margins = "y", size = 5,
          color = "steelblue4")
        Plot2 <- ggExtra::ggMarginal(
          p = Plot2, type = "density", margins = "y", size = 5,
          color = "steelblue4")

        invisible(gc())

        return(
          tibble::tibble(
            Var_Sp = x, Plot = list(Plot), PlotFixedY = list(Plot2)))
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c(
        "PlotObj_Beta", "Cols", "VarRanges", "EnvFile", "FromHPC"),
      future.packages = c("dplyr", "coda", "ggplot2", "ggExtra", "ggtext")) %>%
      dplyr::bind_rows() %>%
      dplyr::right_join(PlotObj_Beta, by = "Var_Sp")

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("sequential", gc = TRUE)
    }


    if (SavePlotData) {
      IASDT.R::CatTime("Save trace plot data", Level = 1)
      IASDT.R::SaveAs(
        InObj = PlotObj_Beta, OutObj = "Convergence_Beta",
        OutPath = file.path(Path_Convergence, "Convergence_Beta.RData"))
    }
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta - 2. by variable ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Trace plots, grouped by variables", Level = 1)

  IASDT.R::CatTime("Preparing data", Level = 2)
  BetaTracePlots_ByVar <- dplyr::arrange(PlotObj_Beta, Variable, IAS_ID) %>%
    dplyr::select(Variable, Plot, PlotFixedY) %>%
    tidyr::nest(data = c("Plot", "PlotFixedY"))

  if (NCores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
  }

  IASDT.R::CatTime("Save plots in parallel", Level = 2)
  BetaTracePlots_ByVar0 <- future.apply::future_lapply(
    X = BetaTracePlots_ByVar$Variable,
    FUN = function(x) {

      VarName <- dplyr::case_when(
        x == "HabLog" ~ "% Habitat coverage",
        x == "RoadRailLog" ~ "Road + Rail intensity",
        x == "BiasLog" ~ "Sampling intensity",
        .default = x)

      BetaPlots <- dplyr::filter(BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(data) %>%
        magrittr::extract2(1) %>%
        dplyr::pull("Plot")

      BetaPlotsFixedY <- dplyr::filter(
        BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(data) %>%
        magrittr::extract2(1) %>%
        dplyr::pull("PlotFixedY")

      PlotTitle <- ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("Convergence of the beta parameter - ", VarName)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold", size = 24, hjust = 0.5))
      PlotTitleFixed <- ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0(
            "Convergence of the beta parameter - ", VarName,
            " (fixed y-axis range)")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold", size = 24, hjust = 0.5))

      BetaPlotList <- tibble::tibble(PlotID = seq_len(length(BetaPlots))) %>%
        dplyr::mutate(Page = ceiling(PlotID / (NRC[2] * NRC[1]))) %>%
        tidyr::nest(.by = "Page", .key = "PlotID") %>%
        dplyr::mutate(
          Plot = purrr::map(
            .x = PlotID,
            .f = ~ {
              BetaPlots[unlist(as.vector(.x))] %>%
                cowplot::plot_grid(
                  plotlist = ., ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitle, ., ncol = 1, rel_heights = c(0.04, 1))
            }),
          PlotFixedY = purrr::map(
            .x = PlotID,
            .f = ~ {
              BetaPlotsFixedY[unlist(as.vector(.x))] %>%
                cowplot::plot_grid(
                  plotlist = ., ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitleFixed, ., ncol = 1, rel_heights = c(0.04, 1))
            }),
          PlotID = NULL)

      invisible({
        grDevices::pdf(
          file = file.path(
            Path_Convergence, paste0("Convergence_Beta_", x, "_FreeY.pdf")),
          width = 18, height = 12)
        purrr::map(BetaPlotList$Plot, grid::grid.draw, recording = FALSE)
        grDevices::dev.off()

        grDevices::pdf(
          file = file.path(
            Path_Convergence, paste0("Convergence_Beta_", x, "_FixedY.pdf")),
          width = 18, height = 12)
        purrr::map(
          .x = BetaPlotList$PlotFixedY, .f = grid::grid.draw)
        grDevices::dev.off()
      })
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.globals = c(
      "BetaTracePlots_ByVar", "Path_Convergence", "NRC", "Cols"),
    future.packages = c("dplyr", "coda", "ggplot2", "ggExtra", "ggtext"))

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("sequential", gc = TRUE)
  }

  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta - 3. by species ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Trace plots, grouped by species", Level = 1)

  IASDT.R::CatTime("Preparing data", Level = 2)
  Order <- stringr::str_remove_all(ModVars, "\\(|\\)")
  BetaTracePlots_BySp <- PlotObj_Beta %>%
    dplyr::slice(order(factor(Variable, levels = Order))) %>%
    dplyr::select(Species, IAS_ID, Plot, Variable) %>%
    tidyr::nest(data = -c("Species", "IAS_ID"))

  IASDT.R::CatTime("Preparing working in parallel", Level = 2)
  if (NCores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
  }


  IASDT.R::CatTime("Save plots in parallel", Level = 2)
  BetaTracePlots_BySp0 <- future.apply::future_lapply(
    X = BetaTracePlots_BySp$Species,
    FUN = function(x) {

      # VarName <- dplyr::case_when(
      #   x == "HabLog" ~ "% Habitat coverage",
      #   x == "RoadRailLog" ~ "Road + Rail intensity",
      #   x == "BiasLog" ~ "Sampling intensity",
      #   .default = x)

      PlotTitle <- ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0(
            "Convergence of the beta parameter - <i>", x, "</i>")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            face = "bold", size = 24, hjust = 0.5))

      SpDT <- dplyr::filter(BetaTracePlots_BySp, Species == x)
      SpPlots <- dplyr::pull(SpDT, "data") %>%
        magrittr::extract2(1) %>%
        dplyr::pull("Plot") %>%
        cowplot::plot_grid(
          plotlist = ., ncol = 4, nrow = 3, align = "hv") %>%
        cowplot::plot_grid(
          PlotTitle, ., ncol = 1, rel_heights = c(0.04, 1))

      invisible({
        grDevices::pdf(
          file = file.path(
            Path_Convergence_BySp,
            paste0(
              "Convergence_Beta_", SpDT$IAS_ID, "_",
              IASDT.R::ReplaceSpace(x), ".pdf")),
          width = 18, height = 15)
        grid::grid.draw(SpPlots)
        grDevices::dev.off()
      })

      return(invisible(NULL))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.globals = c("BetaTracePlots_BySp", "Path_Convergence_BySp"),
    future.packages = c("dplyr", "coda", "ggplot2", "ggExtra", "ggtext"))

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("sequential", gc = TRUE)
  }


  rm(BetaTracePlots_BySp0)

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  return(invisible(NULL))
}
