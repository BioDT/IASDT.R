## |------------------------------------------------------------------------| #
# Convergence_Plot ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of a selected model
#'
#' The `Convergence_Plot()` function generates and saves convergence diagnostics
#' plots for the `rho`, `alpha`, `omega`, and `beta` parameters in an Hmsc
#' model. These plots help assess whether the MCMC chains have reached
#' stationarity. It supports parallel processing and can work with models fitted
#' on HPC environments.
#' @param Path_Coda Character. Path to a saved coda object containing MCMC
#'   samples.
#' @param Path_Model Character. Path to a saved Hmsc model object.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical. Whether the processing is being done on an
#'   High-Performance Computing (HPC) environment, to adjust file paths
#'   accordingly. Default: `TRUE`.
#' @param Title Character. Title for **rho** and **alpha** convergence plots.
#'   Default: " "
#' @param NOmega Integer. Number of species interactions sampled for Omega
#'   parameter diagnostics. Default: 1000L
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param NRC Numeric vector. Grid layout (rows&times;columns) for arranging
#'   alpha parameter plots. Default: `c(2, 2)`. If `NULL`, the layout is
#'   automatically determined based on the number of alpha levels.
#' @param Beta_NRC Numeric vector. The grid layout (rows&times;columns) for
#'   arranging beta parameter plots. Default: `c(3, 3)`.
#' @param SavePlotData Logical. If `TRUE` (default), saves the plot data as
#'   `.RData` files.
#' @param PagePerFile Integer. Number of plots per page in the Omega parameter
#'   output. Default: 20L.
#' @param Cols Character vector. MCMC chain colors (optional). Default: `NULL`.
#' @param Post `mcmc.list` or character. Either an MCMC object (`mcmc.list`)
#'   containing posterior samples, or a file path to a saved coda object.
#' @param Model `Hmsc` object or character. Either a fitted Hmsc model
#'   object, or a file path to a saved Hmsc model.
#' @param AddFooter Logical. If `TRUE` (default), adds a footer with page
#'   numbers to each plot.
#' @param AddTitle Logical. If `TRUE` (default), adds the main title (`Title`)
#'   to the plot.
#' @param MarginType Character. The type of marginal plot to add to
#'   the main plot. Valid options are "histogram" (default) or "density".
#'
#' @details `Convergence_Alpha()` and `Convergence_Rho()` are internal functions
#'   and should not be called directly.
#' @rdname Convergence_plots
#' @name Convergence_plots
#' @order 1
#' @author Ahmed El-Gabbas
#' @export

Convergence_Plot <- function(
    Path_Coda = NULL, Path_Model = NULL, EnvFile = ".env",
    FromHPC = TRUE, Title = " ", NOmega = 1000L, NCores = 8L,
    NRC = c(2L, 2L), Beta_NRC = c(3L, 3L), SavePlotData = TRUE,
    PagePerFile = 20L, Cols = NULL, MarginType = "histogram") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Coda) || is.null(Path_Model) || is.null(NCores)) {
    stop(
      "Path_Coda, Path_Model, and NCores cannot be empty", call. = FALSE)
  }

  if (length(MarginType) != 1) {
    stop("`MarginType` must be a single value.", call. = FALSE)
  }

  if (!MarginType %in% c("histogram", "density")) {
    stop(
      "`MarginType` must be either 'histogram' or 'density'.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpComb <- `2.5%` <- `97.5%` <- Class <- Order <- Family <- DT <- IAS_ID <-
    Species <- Variable <- data <- PlotID <- File <- Page <- Iter <- Value <-
    Chain <- y <- label <- Var_Sp <- CI_025 <- CI_975 <- Var_Min <- Var_Max <-
    Variable2 <- Plot_File <- Var_Sp2 <- Species_name <- Species_File <-
    Path_PA <- VarDesc <- Term <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  IASDT.R::CatTime("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Path_Coda", "Path_Model"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("NOmega", "NCores", "NRC"))
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # # Load species summary
  IASDT.R::CatTime("Load species summary")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  SpSummary <- IASDT.R::Path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    stop(paste0(SpSummary, " file does not exist"), call. = FALSE)
  }

  SpSummary <- readr::read_csv(SpSummary, show_col_types = FALSE) %>%
    dplyr::select(Species = Species_name, Species_File)

  # # ..................................................................... ###

  # Create path ------

  IASDT.R::CatTime("Create path")
  Path_Convergence <- dirname(dirname(Path_Coda)) %>%
    IASDT.R::Path("Model_Convergence")
  Pah_Beta_Data <- IASDT.R::Path(Path_Convergence, "Beta_Data")
  Path_Convergence_BySp <- IASDT.R::Path(Path_Convergence, "Beta_BySpecies")
  fs::dir_create(c(Path_Convergence, Path_Convergence_BySp, Pah_Beta_Data))

  # # ..................................................................... ###

  # Prepare convergence data ------

  IASDT.R::CatTime("Prepare convergence data")

  if (!file.exists(Path_Model)) {
    stop("`Path_Model` does not exist", call. = FALSE)
  }

  if (!file.exists(Path_Coda)) {
    stop("`Path_Coda` does not exist", call. = FALSE)
  }

  IASDT.R::CatTime("Loading coda object", Level = 1)
  Coda_Obj <- IASDT.R::LoadAs(Path_Coda)

  IASDT.R::CatTime("Loading fitted model object", Level = 1)
  Model <- IASDT.R::LoadAs(Path_Model)

  # Model variables
  ModVars <- Model$covNames

  # Number of chains
  NChains <- length(Model$postList)

  # Number of samples
  SampleSize <- Model$samples

  #  Plotting colours
  if (is.null(Cols)) {
    Cols <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }
  if (length(Cols) != NChains) {
    warning(
      "The length of provided colours != number of chains", .call. = FALSE)
    Cols <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }

  # # ..................................................................... ###

  # Rho ------

  if ("Rho" %in% names(Coda_Obj)) {
    IASDT.R::CatTime("Rho")

    FileConv_Rho <- IASDT.R::Path(Path_Convergence, "Convergence_Rho.RData")

    if (IASDT.R::CheckData(FileConv_Rho, warning = FALSE)) {
      IASDT.R::CatTime("Loading plotting data", Level = 1)
      PlotObj_Rho <- IASDT.R::LoadAs(FileConv_Rho)
    } else {
      IASDT.R::CatTime("Prepare plot", Level = 1)
      PlotObj_Rho <- IASDT.R::Convergence_Rho(
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
    grDevices::cairo_pdf(
      filename = IASDT.R::Path(Path_Convergence, "Convergence_Rho.pdf"),
      width = 18, height = 12, onefile = TRUE)
    plot(PlotObj_Rho)
    grDevices::dev.off()

    rm(PlotObj_Rho, envir = environment())
  }

  # # ..................................................................... ###

  # Alpha  ------

  IASDT.R::CatTime("Alpha")

  FileConv_Alpha <- IASDT.R::Path(Path_Convergence, "Convergence_Alpha.RData")

  if (file.exists(FileConv_Alpha)) {
    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Alpha <- IASDT.R::LoadAs(FileConv_Alpha)
  } else {
    IASDT.R::CatTime("Prepare plot", Level = 1)
    PlotObj_Alpha <- IASDT.R::Convergence_Alpha(
      Post = Coda_Obj, Model = Model, Title = Title, NRC = NRC,
      AddFooter = FALSE, AddTitle = FALSE, Cols = Cols)

    if (SavePlotData) {
      IASDT.R::CatTime("Save plotting data", Level = 1)
      IASDT.R::SaveAs(
        InObj = PlotObj_Alpha, OutObj = "Convergence_Alpha",
        OutPath = IASDT.R::Path(Path_Convergence, "Convergence_Alpha.RData"))
    }
  }

  IASDT.R::CatTime("Save plots", Level = 1)
  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = IASDT.R::Path(Path_Convergence, "Convergence_Alpha.pdf"),
    width = 18, height = 14, onefile = TRUE)
  print(PlotObj_Alpha)
  grDevices::dev.off()

  Obj_Omega <- Coda_Obj$Omega[[1]]
  Obj_Beta <- Coda_Obj$Beta

  rm(Model, Coda_Obj, PlotObj_Alpha, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Omega  ------

  IASDT.R::CatTime("Omega")

  FileConv_Omega <- IASDT.R::Path(Path_Convergence, "Convergence_Omega.qs2")

  if (file.exists(FileConv_Omega)) {
    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Omega <- IASDT.R::LoadAs(FileConv_Omega)
  } else {
    IASDT.R::CatTime("Coda to tibble", Level = 1)
    OmegaDF <- IASDT.R::Coda_to_tibble(
      CodaObj = Obj_Omega, Type = "omega", NOmega = NOmega,
      EnvFile = EnvFile, FromHPC = FromHPC)
    invisible(gc())
    SelectedCombs <- unique(OmegaDF$SpComb)

    IASDT.R::CatTime("Prepare confidence interval data", Level = 1)
    CI <- purrr::map(.x = Obj_Omega, .f = ~ .x[, SelectedCombs]) %>%
      coda::mcmc.list() %>%
      summary(quantiles = c(0.025, 0.975)) %>%
      magrittr::extract2("quantiles") %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "SpComb") %>%
      stats::setNames(c("SpComb", "CI_25", "CI_975"))

    OmegaDF <- dplyr::left_join(OmegaDF, CI, by = "SpComb")
    OmegaNames <- attr(Obj_Omega[[1]], "dimnames")[[2]]

    IASDT.R::CatTime("Prepare plots", Level = 1)
    PlotObj_Omega <- purrr::map_dfr(
      .x = seq_len(NOmega),
      .f = function(x) {

        CombData <- dplyr::filter(OmegaDF, SpComb == SelectedCombs[x])
        CurrPost <- purrr::map(
          .x = Obj_Omega,
          .f = ~ .x[, which(OmegaNames == CombData$SpComb)]) %>%
          coda::as.mcmc.list()

        ## Gelman convergence diagnostic
        Label_Gelman <- try(
          coda::gelman.diag(CurrPost, multivariate = FALSE),
          silent = TRUE)

        if (inherits(Label_Gelman, "try-error")) {
          Label_Gelman <- coda::gelman.diag(
            CurrPost, multivariate = FALSE, autoburnin = FALSE)
        }

        Label_Gelman <- Label_Gelman %>%
          magrittr::extract2("psrf") %>%
          magrittr::extract(1) %>%
          round(3) %>%
          paste0("<b><i>Gelman convergence diagnostic:</i></b> ", .) %>%
          data.frame(x = -Inf, y = Inf, label = .)

        ## Effective sample size
        Label_ESS <- coda::effectiveSize(CurrPost) %>%
          magrittr::divide_by(NChains) %>%
          round(1) %>%
          paste0(
            "<b><i>Mean effective sample size:</i></b> ", ., " / ", SampleSize)
        CurrCI <- dplyr::select(CombData, c("CI_25", "CI_975")) %>%
          round(2) %>%
          unlist()
        Label_CI <- CurrCI %>%
          paste0(collapse = " to ") %>%
          paste0("<b><i>95% credible interval:</i></b> ", .)
        Label_ESS_CI <- data.frame(
          x = -Inf, y = -Inf, label = paste0(Label_ESS, "<br>", Label_CI))

        Label_Panel <- c(CombData$IAS1, CombData$IAS2) %>%
          sort() %>%
          paste0("<i>", ., "</i>") %>%
          paste0(collapse = " & <br>") %>%
          paste0(., " ") %>%
          data.frame(x = Inf, y = Inf, label = .)

        Plot <- ggplot2::ggplot(
          data = CombData$DT[[1]],
          mapping = ggplot2::aes(
            x = Iter, y = Value, color = factor(Chain))) +
          ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
          ggplot2::geom_smooth(
            method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
          ggplot2::geom_point(alpha = 0) +
          ggplot2::geom_hline(
            yintercept = CurrCI, linetype = "dashed", color = "black",
            linewidth = 1) +
          # Ensure that y-axis always show 0
          ggplot2::geom_hline(
            yintercept = 0, linetype = "dashed",
            color = "transparent", linewidth = 0.6) +
          ggplot2::scale_color_manual(values = Cols) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label), size = 5,
            data = Label_Panel, inherit.aes = FALSE, colour = "blue",
            hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = Label_Gelman, inherit.aes = FALSE, size = 6, hjust = 0,
            vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = Label_ESS_CI, inherit.aes = FALSE, size = 6, hjust = 0,
            vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
          ggplot2::theme_bw() +
          ggplot2::xlab(NULL) +
          ggplot2::ylab(NULL) +
          ggplot2::theme(
            axis.text = ggplot2::element_text(size = 12),
            legend.position = "none")

        if (MarginType == "histogram") {
          Plot <- ggExtra::ggMarginal(
            p = Plot, type = MarginType, margins = "y", size = 6,
            color = "steelblue4", fill = "steelblue4", bins = 100)
        } else {
          Plot <- ggExtra::ggMarginal(
            p = Plot, type = MarginType, margins = "y", size = 6,
            color = "steelblue4")
        }
        # Making marginal background matching the plot background
        # https://stackoverflow.com/a/78196022/3652584
        Plot$layout$t[1] <- 1
        Plot$layout$r[1] <- max(Plot$layout$r)

        return(tibble::tibble(SpComb = CombData$SpComb, Plot = list(Plot)))
      }
    )

    if (SavePlotData) {
      IASDT.R::CatTime("Save plot data", Level = 1)
      IASDT.R::SaveAs(InObj = PlotObj_Omega, OutPath = FileConv_Omega)
    }
    rm(Obj_Omega, OmegaDF, SelectedCombs, CI, OmegaNames, envir = environment())
    invisible(gc())
  }

  IASDT.R::CatTime("Arrange plots", Level = 1)
  OmegaPlotList <- tibble::tibble(PlotID = seq_len(nrow(PlotObj_Omega))) %>%
    dplyr::mutate(
      File = ceiling(PlotID / (PagePerFile * NRC[2] * NRC[1])),
      Page = ceiling(PlotID / (NRC[2] * NRC[1]))) %>%
    tidyr::nest(.by = c("File", "Page"), .key = "PlotID") %>%
    dplyr::mutate(
      PlotID = purrr::map(PlotID, ~ unlist(as.vector(.x))),
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

          cowplot::plot_grid(
            plotlist = PlotObj_Omega$Plot[PlotID],
            ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
            cowplot::plot_grid(PlotTitle, ., ncol = 1, rel_heights = c(0.05, 1))
        }
      ))

  IASDT.R::CatTime("Save plots", Level = 1)
  purrr::walk(
    .x = seq_along(unique(OmegaPlotList$File)),
    .f = ~ {
      invisible({
        CurrPlotOrder <- dplyr::filter(OmegaPlotList, File == .x)
        grDevices::cairo_pdf(
          filename = IASDT.R::Path(
            Path_Convergence, paste0("Convergence_Omega_", .x, ".pdf")),
          width = 18, height = 14, onefile = TRUE)
        purrr::map(CurrPlotOrder$PlotID, grid::grid.draw, recording = FALSE)
        grDevices::dev.off()
      })
    })

  rm(OmegaPlotList, PlotObj_Omega, Obj_Omega, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Beta - 1. Prepare data ------

  IASDT.R::CatTime("Beta")

  FileConv_Beta <- IASDT.R::Path(Path_Convergence, "Convergence_Beta.RData")

  if (file.exists(FileConv_Beta)) {

    IASDT.R::CatTime("Loading plotting data", Level = 1)
    PlotObj_Beta <- IASDT.R::LoadAs(FileConv_Beta)

  } else {

    HTML1 <- "<span style='color:blue;'><b>"
    HTML2 <- "</b></span>"
    Intercept <- c(
      Variable = "Intercept",
      VarDesc = stringr::str_glue("{HTML1}Intercept{HTML2}"))

    LinearTerms <- tibble::tribble(
      ~Variable, ~VarDesc,
      "RiversLog", "% Habitat coverage",
      "RoadRailLog", "Road + Rail intensity",
      "EffortsLog", "Sampling efforts",
      "HabLog", "River length") %>%
      dplyr::mutate(
        VarDesc = stringr::str_glue(
          "{HTML1}{stringr::str_to_sentence(VarDesc)}{HTML2}"),
        VarDesc = paste0(VarDesc, " (log<sub>10</sub>(x + 0.1))")) %>%
      dplyr::bind_rows(Intercept, .)

    HTML1 <- "<span style='color:blue;'><b>"
    HTML2 <- "</b></span><span style='color:grey;'>"
    HTML3 <- "</span>"
    VarsDesc <- tibble::tribble(
      ~Variable, ~VarDesc,
      "bio1", "annual mean temperature",
      "bio2", "mean diurnal range",
      "bio3", "isothermality (bio2/bio7)",
      "bio4", "temperature seasonality",
      "bio5", "max temperature of warmest month",
      "bio6", "temperature of the coldest month",
      "bio7", "temperature annual range (bio5-bio6)",
      "bio8", "temperatures of the wettest quarter",
      "bio9", "mean temperature of driest quarter",
      "bio10", "mean temperature of warmest quarter",
      "bio11", "mean temperature of coldest quarter",
      "bio12", "annual precipitation amount",
      "bio13", "precipitation of wettest month",
      "bio14", "precipitation of driest month",
      "bio15", "precipitation seasonality",
      "bio16", "precipitation of wettest quarter",
      "bio17", "precipitation of driest quarter",
      "bio18", "precipitation of the warmest quarter",
      "bio19", "precipitation of coldest quarter",
      "npp", "net primary productivity") %>%
      tidyr::expand_grid(Term = c("L", "Q")) %>%
      dplyr::mutate(
        VarDesc = stringr::str_glue(
          "{HTML1}{stringr::str_to_sentence(Variable)}\\
        {HTML2} [{VarDesc}]{HTML3} - {Term}"),
        VarDesc = stringr::str_replace(
          VarDesc, " - L$", "&nbsp&nbsp&mdash&nbsp&nbspLinear"),
        VarDesc = stringr::str_replace(
          VarDesc, " - Q$", "&nbsp&nbsp&mdash&nbsp&nbspQuadratic"),
        Variable = paste0(Variable, "_", Term),
        Term = NULL) %>%
      dplyr::bind_rows(LinearTerms)

    IASDT.R::CatTime("Prepare trace plots", Level = 1)

    BetaNames <- attr(Obj_Beta[[1]], "dimnames")[[2]]

    IASDT.R::CatTime("Prepare 95% credible interval data", Level = 2)
    CI <- summary(Obj_Beta, quantiles = c(0.025, 0.975))$quantiles %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "Var_Sp") %>%
      dplyr::rename(CI_025 = `2.5%`, CI_975 = `97.5%`)

    IASDT.R::CatTime("Coda to tibble", Level = 2)
    Beta_DF <- IASDT.R::Coda_to_tibble(
      CodaObj = Obj_Beta, Type = "beta",
      EnvFile = EnvFile, FromHPC = FromHPC) %>%
      dplyr::left_join(CI, by = "Var_Sp")

    # Variable ranges
    IASDT.R::CatTime("Variable ranges", Level = 2)
    VarRanges <- dplyr::arrange(Beta_DF, Variable, IAS_ID) %>%
      dplyr::select(Variable, DT) %>%
      dplyr::mutate(
        Range = purrr::map(.x = DT, .f = ~ range(dplyr::pull(.x, Value)))) %>%
      dplyr::select(-DT) %>%
      tidyr::nest(data = c("Range")) %>%
      dplyr::mutate(
        Range = purrr::map(
          .x = data,
          .f = ~ {
            dplyr::pull(.x, Range) %>%
              as.vector() %>%
              range() %>%
              purrr::set_names(c("Var_Min", "Var_Max"))
          })) %>%
      dplyr::select(-data) %>%
      tidyr::unnest_wider("Range")

    # Species taxonomy
    IASDT.R::CatTime("Species taxonomy", Level = 2)
    SpeciesTaxonomy <- IASDT.R::GetSpeciesName(
      EnvFile = EnvFile, FromHPC = FromHPC) %>%
      dplyr::select(IAS_ID, Class, Order, Family)

    IASDT.R::CatTime("Prepare Beta data", Level = 2)
    Beta_DF <- Beta_DF %>%
      dplyr::left_join(VarRanges, by = "Variable") %>%
      dplyr::left_join(SpeciesTaxonomy, by = "IAS_ID") %>%
      dplyr::mutate(
        Var_Sp2 = paste0(Variable, "_", IAS_ID),
        Var_Sp_File = IASDT.R::Path(Pah_Beta_Data, paste0(Var_Sp2, ".RData")),
        DT = purrr::pmap(
          .l = list(
            Var_Sp, DT, CI_025, CI_975, Var_Min,
            Var_Max, Class, Order, Family),
          .f = function(Var_Sp, DT, CI_025, CI_975, Var_Min,
                        Var_Max, Class, Order, Family) {
            Beta_ID <- which(BetaNames == Var_Sp)
            Post <- coda::as.mcmc.list(Obj_Beta[, Beta_ID])

            Gelman <- try(
              coda::gelman.diag(Post, multivariate = FALSE),
              silent = TRUE)

            if (inherits(Gelman, "try-error")) {
              Gelman <- coda::gelman.diag(
                Post, multivariate = FALSE, autoburnin = FALSE)
            }

            ESS <- coda::effectiveSize(Post)

            list(
              DT = DT, CI_025 = CI_025, CI_975 = CI_975,
              Var_Min = Var_Min, Var_Max = Var_Max,
              Class = Class, Order = Order, Family = Family,
              Beta_ID = Beta_ID, Post = Post, Gelman = Gelman, ESS = ESS)
          }
        )) %>%
      dplyr::select(
        -tidyselect::all_of(
          c(
            "CI_025", "CI_975", "Var_Min", "Var_Max",
            "Class", "Order", "Family")))

    Beta_DF %>%
      dplyr::group_by(Var_Sp) %>%
      dplyr::group_split() %>%
      purrr::walk(
        .f = ~ IASDT.R::SaveAs(
          InObj = .x$DT[[1]], OutObj = .x$Var_Sp2, OutPath = .x$Var_Sp_File))

    Beta_DF <- dplyr::select(Beta_DF, -DT)

    rm(CI, VarRanges, SpeciesTaxonomy, Obj_Beta, envir = environment())
    invisible(gc())

    # Prepare working on parallel
    IASDT.R::CatTime("Prepare working on parallel", Level = 2)
    if (NCores > 1) {
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(min(NCores, nrow(Beta_DF)))
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    } else {
      future::plan("future::sequential", gc = TRUE)
    }

    # Prepare plots
    IASDT.R::CatTime("Prepare plots", Level = 2)
    PlotObj_Beta <- future.apply::future_lapply(
      X = seq_len(nrow(Beta_DF)),
      FUN = function(x) {
        Var_Sp <- Beta_DF$Var_Sp[x]
        Species <- Beta_DF$Species[x]
        Curr_IAS <- Beta_DF$IAS_ID[x]
        Var_Sp_File <- Beta_DF$Var_Sp_File[x]
        Plot_File <- stringr::str_replace(Var_Sp_File, ".RData$", "_Plots.qs2")

        DT_all <- IASDT.R::LoadAs(Var_Sp_File)
        DT_all$Post <- NULL
        invisible(gc())

        ## Gelman convergence diagnostic
        Label_Gelman <- round(DT_all$Gelman$psrf, 3) %>%
          paste0(collapse = " / ") %>%
          paste0("<b><i>Gelman convergence diagnostic:</i></b> ", .) %>%
          data.frame(x = Inf, y = -Inf, label = .)

        ## Effective sample size / CI
        Label_ESS <- round(DT_all$ESS / NChains) %>%
          paste0(
            "<b><i>Mean effective sample size:</i></b> ", ., " / ",
            SampleSize)
        CurrCI <- c(DT_all$CI_025, DT_all$CI_975)
        Label_CI <- paste0(round(CurrCI, 4), collapse = " to ") %>%
          paste0("<b><i>95% credible interval:</i></b> ", .)
        Label_ESS_CI <- data.frame(
          x = -Inf, y = -Inf, label = paste0(Label_ESS, "<br>", Label_CI))

        Label_Panel <- data.frame(
          x = Inf, y = Inf,
          label = paste0("<br><b><i>", Species, "</i></b>"))

        PanelTitle <- c(DT_all$Class, DT_all$Order, DT_all$Family) %>%
          paste0(collapse = " | ") %>%
          paste0("<b>", ., "</b>") %>%
          paste0("<br>", Curr_IAS) %>%
          data.frame(x = -Inf, y = Inf, label = .)

        Plot <- ggplot2::ggplot(
          data = DT_all$DT,
          mapping = ggplot2::aes(
            x = Iter, y = Value, color = factor(Chain))) +
          ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
          ggplot2::geom_smooth(
            method = "loess", formula = y ~ x,
            se = FALSE, linewidth = 0.8) +
          ggplot2::geom_point(alpha = 0) +
          ggplot2::geom_hline(
            yintercept = CurrCI, linetype = "dashed", color = "black",
            linewidth = 1) +
          # Ensure that y-axis always show 0
          ggplot2::geom_hline(
            yintercept = 0, linetype = "dashed",
            color = "transparent", linewidth = 0.6) +
          ggplot2::scale_color_manual(values = Cols) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = Label_Gelman, inherit.aes = FALSE, size = 3.5, hjust = 1,
            vjust = -0, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = Label_ESS_CI, inherit.aes = FALSE, size = 3.5, hjust = 0,
            vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = Label_Panel, inherit.aes = FALSE, colour = "blue",
            hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = PanelTitle, inherit.aes = FALSE, hjust = 0, vjust = 1,
            lineheight = 0, fill = NA, label.color = NA) +
          ggplot2::theme_bw() +
          ggplot2::xlab(NULL) +
          ggplot2::ylab(NULL) +
          ggplot2::theme(
            legend.position = "none",
            axis.text = ggplot2::element_text(size = 12))

        suppressMessages({
          Plot2 <- Plot +
            ggplot2::scale_y_continuous(
              limits = c(DT_all$Var_Min, DT_all$Var_Max))
        })

        if (MarginType == "histogram") {
          Plot_Marginal <- ggExtra::ggMarginal(
            p = Plot, type = MarginType, margins = "y", size = 6,
            color = "steelblue4", fill = "steelblue4", bins = 100)
        } else {
          Plot_Marginal <- ggExtra::ggMarginal(
            p = Plot, type = MarginType, margins = "y", size = 6,
            color = "steelblue4")
        }
        # Making marginal background matching the plot background
        # https://stackoverflow.com/a/78196022/3652584
        Plot_Marginal$layout$t[1] <- 1
        Plot_Marginal$layout$r[1] <- max(Plot_Marginal$layout$r)

        if (MarginType == "histogram") {
          Plot2_Marginal <- ggExtra::ggMarginal(
            p = Plot2, type = MarginType, margins = "y", size = 6,
            color = "steelblue4", fill = "steelblue4", bins = 100)
        } else {
          Plot2_Marginal <- ggExtra::ggMarginal(
            p = Plot2, type = MarginType, margins = "y", size = 6,
            color = "steelblue4")
        }
        # Making marginal background matching the plot background
        # https://stackoverflow.com/a/78196022/3652584
        Plot2_Marginal$layout$t[1] <- 1
        Plot2_Marginal$layout$r[1] <- max(Plot2_Marginal$layout$r)

        IASDT.R::SaveAs(
          InObj = list(
            Plot = Plot, Plot_Marginal = Plot_Marginal,
            PlotFixedY_Marginal = Plot2_Marginal),
          OutObj = paste0(Beta_DF$Var_Sp2[x], "_Plots"),
          OutPath = Plot_File)

        return(tibble::tibble(Var_Sp = Var_Sp, Plot_File = Plot_File))
      },
      future.seed = TRUE,
      future.globals = c(
        "Beta_DF", "NChains", "SampleSize", "Cols", "MarginType"),
      future.packages = c(
        "dplyr", "ggplot2", "ggtext", "magrittr", "coda", "IASDT.R")) %>%
      dplyr::bind_rows() %>%
      dplyr::left_join(Beta_DF, ., by = "Var_Sp") %>%
      dplyr::left_join(VarsDesc, by = "Variable")

    # Stopping cluster
    if (NCores > 1) {
      IASDT.R::CatTime("Stopping cluster", Level = 2)
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

    rm(Beta_DF, BetaNames, envir = environment())
    invisible(gc())

    if (SavePlotData) {
      IASDT.R::CatTime("Save trace plot data", Level = 2)
      IASDT.R::SaveAs(
        InObj = PlotObj_Beta, OutObj = "Convergence_Beta",
        OutPath = FileConv_Beta)
    }
  }

  # # ..................................................................... ###

  # Beta - 2. by variable ------

  IASDT.R::CatTime("Trace plots, grouped by variables", Level = 1)

  IASDT.R::CatTime("Preparing data", Level = 2)
  BetaTracePlots_ByVar <- dplyr::arrange(PlotObj_Beta, Variable, IAS_ID) %>%
    dplyr::select(Variable, Plot_File, VarDesc) %>%
    tidyr::nest(Plot_File = "Plot_File") %>%
    dplyr::mutate(Plot_File = purrr::map(Plot_File, unlist))

  # Prepare working on parallel
  IASDT.R::CatTime("Prepare working on parallel", Level = 2)
  if (NCores > 1) {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(BetaTracePlots_ByVar)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  } else {
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::CatTime("Save plots", Level = 2)
  BetaTracePlots_ByVar0 <- future.apply::future_lapply(
    X = BetaTracePlots_ByVar$Variable,
    FUN = function(x) {

      VarDesc <- BetaTracePlots_ByVar %>%
        dplyr::filter(Variable == x) %>%
        dplyr::pull(VarDesc)

      Plots <- dplyr::filter(BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(Plot_File) %>%
        magrittr::extract2(1) %>%
        purrr::map(
          .f = ~ {
            IASDT.R::LoadAs(.x) %>%
              magrittr::extract(c("Plot_Marginal", "PlotFixedY_Marginal"))
          })

      BetaPlots <- purrr::map(
        .x = Plots, .f = ~ magrittr::extract2(.x, "Plot_Marginal"))
      BetaPlotsFixedY <- purrr::map(
        .x = Plots, .f = ~ magrittr::extract2(.x, "PlotFixedY_Marginal"))
      rm(Plots, envir = environment())

      PlotTitle <- ggplot2::ggplot() +
        ggplot2::labs(title = VarDesc) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            size = 24, hjust = 0.5, margin = ggplot2::margin(t = 15, b = 15)))

      PlotTitleFixed <- ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0(VarDesc, " (fixed y-axis range)")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            size = 24, hjust = 0.5, margin = ggplot2::margin(t = 15, b = 15)))

      BetaPlotList <- tibble::tibble(PlotID = seq_len(length(BetaPlots))) %>%
        dplyr::mutate(Page = ceiling(PlotID / (NRC[2] * NRC[1]))) %>%
        tidyr::nest(.by = "Page", .key = "PlotID") %>%
        dplyr::mutate(
          Plot = purrr::map(
            .x = PlotID,
            .f = ~ {
              ID <- unlist(as.vector(.x))
              cowplot::plot_grid(
                plotlist = BetaPlots[ID],
                ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitle, ., ncol = 1, rel_heights = c(0.035, 1))
            }),
          PlotFixedY = purrr::map(
            .x = PlotID,
            .f = ~ {
              ID <- unlist(as.vector(.x))
              cowplot::plot_grid(
                plotlist = BetaPlotsFixedY[ID],
                ncol = NRC[2], nrow = NRC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitleFixed, .,
                  ncol = 1, rel_heights = c(0.035, 1))
            }))

      invisible({
        grDevices::cairo_pdf(
          filename = IASDT.R::Path(
            Path_Convergence, paste0("Convergence_Beta_", x, "_FreeY.pdf")),
          width = 18, height = 13, onefile = TRUE)
        purrr::walk(BetaPlotList$Plot, grid::grid.draw, recording = FALSE)
        grDevices::dev.off()

        grDevices::cairo_pdf(
          filename = IASDT.R::Path(
            Path_Convergence, paste0("Convergence_Beta_", x, "_FixedY.pdf")),
          width = 18, height = 13, onefile = TRUE)
        purrr::walk(.x = BetaPlotList$PlotFixedY, .f = grid::grid.draw)
        grDevices::dev.off()
      })

      rm(
        PlotTitle, PlotTitleFixed, BetaPlotsFixedY, BetaPlots, BetaPlotList,
        envir = environment())

      invisible(gc())
      return(NULL)
    },
    future.seed = TRUE,
    future.globals = c("BetaTracePlots_ByVar", "NRC", "Path_Convergence"),
    future.packages = c(
      "dplyr", "ggplot2", "magrittr", "purrr", "IASDT.R",
      "tibble", "tidyr", "cowplot"))

  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar, envir = environment())
  invisible(gc())

  # Stopping cluster
  if (NCores > 1) {
    IASDT.R::CatTime("Stopping cluster", Level = 2)
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Beta - 3. by species ------

  IASDT.R::CatTime("Trace plots, grouped by species", Level = 1)

  IASDT.R::CatTime("Preparing data", Level = 2)
  Order <- stringr::str_remove_all(ModVars, "\\(|\\)")
  BetaTracePlots_BySp <- PlotObj_Beta %>%
    dplyr::arrange(Species, factor(Variable, levels = Order)) %>%
    dplyr::select(Species, IAS_ID, Plot_File, Variable, VarDesc) %>%
    tidyr::nest(data = -c("Species", "IAS_ID")) %>%
    dplyr::left_join(SpSummary, by = "Species")

  # Prepare working on parallel
  IASDT.R::CatTime("Prepare working on parallel", Level = 2)
  if (NCores > 1) {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(BetaTracePlots_BySp)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  } else {
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::CatTime("Save plots", Level = 2)
  BetaTracePlots_BySp0 <- future.apply::future_lapply(
    X = BetaTracePlots_BySp$Species,
    FUN = function(x) {

      PlotTitle <- ggplot2::ggplot() +
        ggplot2::labs(title = paste0("<i>", x, "</i>")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            face = "bold", size = 24, hjust = 0.5))

      SpDT <- dplyr::filter(BetaTracePlots_BySp, Species == x)

      VarOrder <- SpDT$data[[1]]$Variable %>%
        stringr::str_subset("Intercept", negate = TRUE) %>%
        gtools::mixedsort() %>%
        c("Intercept", .)

      SpPlots <- SpDT$data[[1]] %>%
        dplyr::arrange(factor(Variable, levels = VarOrder)) %>%
        dplyr::mutate(
          Plot = purrr::map2(
            .x = Plot_File, .y = VarDesc,
            .f = ~ {
              Plot <- IASDT.R::LoadAs(.x)$Plot +
                ggplot2::ggtitle(.y) +
                ggplot2::theme(
                  plot.title = ggtext::element_markdown(
                    hjust = 0.5, size = 15, colour = "red",
                    margin = ggplot2::margin(0, 0, -2.5, 0)))

              if (MarginType == "histogram") {
                Plot <- ggExtra::ggMarginal(
                  p = Plot, type = MarginType, margins = "y", size = 6,
                  color = "steelblue4", fill = "steelblue4", bins = 100)
              } else {
                Plot <- ggExtra::ggMarginal(
                  p = Plot, type = MarginType, margins = "y", size = 6,
                  color = "steelblue4")
              }
              Plot$layout$t[1] <- 1
              Plot$layout$r[1] <- max(Plot$layout$r)
              return(Plot)
            })) %>%
        dplyr::pull("Plot")

      NumPages <- ceiling(length(SpPlots) / (Beta_NRC[1] * Beta_NRC[2]))

      SpPlots2 <- IASDT.R::SplitVector(
        Vector = seq_len(length(SpPlots)), NSplit = NumPages) %>%
        purrr::map(
          .f = ~ {
            SpPlots[.x] %>%
              cowplot::plot_grid(
                plotlist = ., ncol = Beta_NRC[2],
                nrow = Beta_NRC[1], align = "hv") %>%
              cowplot::plot_grid(
                PlotTitle, ., ncol = 1, rel_heights = c(0.03, 1))
          })

      invisible({
        grDevices::cairo_pdf(
          filename = IASDT.R::Path(
            Path_Convergence_BySp,
            paste0(
              "Convergence_Beta_", SpDT$IAS_ID, "_",
              SpDT$Species_File, ".pdf")),
          width = 23, height = 17, onefile = TRUE)
        purrr::walk(SpPlots2, grid::grid.draw)
        grDevices::dev.off()
      })

      rm(PlotTitle, envir = environment())
      return(invisible(NULL))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.globals = c(
      "MarginType", "BetaTracePlots_BySp", "Path_Convergence_BySp", "Beta_NRC"),
    future.packages = c("dplyr", "coda", "ggplot2", "ggExtra", "ggtext"))

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  rm(BetaTracePlots_BySp0, envir = environment())

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Plot model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
