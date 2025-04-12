## |------------------------------------------------------------------------| #
# convergence_plot ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of a selected model
#'
#' The `convergence_plot()` function generates and saves convergence diagnostics
#' plots for the `rho`, `alpha`, `omega`, and `beta` parameters in an Hmsc
#' model. These plots help assess whether the MCMC chains have reached
#' stationarity. It supports parallel processing and can work with models fitted
#' on HPC environments.
#' @param path_coda Character. Path to a saved coda object containing MCMC
#'   samples.
#' @param path_model Character. Path to a saved Hmsc model object.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param title Character. title for **rho** and **alpha** convergence plots.
#'   Default: " "
#' @param n_omega Integer. Number of species interactions sampled for Omega
#'   parameter diagnostics. Default: 1000L
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param n_RC Numeric vector. Grid layout (rows&times;columns) for arranging
#'   alpha parameter plots. Default: `c(2, 2)`. If `NULL`, the layout is
#'   automatically determined based on the number of alpha levels.
#' @param beta_n_RC Numeric vector. The grid layout (rows&times;columns) for
#'   arranging beta parameter plots. Default: `c(3, 3)`.
#' @param save_plotting_data Logical. If `TRUE` (default), saves the plot data
#'   as `.RData` files.
#' @param pages_per_file Integer. Number of plots per page in the Omega
#'   parameter output. Default: 20L.
#' @param chain_colors Character vector. MCMC chain colors (optional). Default:
#'   `NULL`.
#' @param posterior `mcmc.list` or character. Either an MCMC object
#'   (`mcmc.list`) containing posterior samples, or a file path to a saved coda
#'   object.
#' @param model_object `Hmsc` object or character. Either a fitted Hmsc model
#'   object, or a file path to a saved Hmsc model.
#' @param add_footer Logical. If `TRUE` (default), adds a footer with page
#'   numbers to each plot.
#' @param add_title Logical. If `TRUE` (default), adds the main title (`title`)
#'   to the plot.
#' @param margin_type Character. The type of marginal plot to add to the main
#'   plot. Valid options are "histogram" (default) or "density".
#'
#' @details `convergence_alpha()` and `convergence_rho()` are internal functions
#'   and should not be called directly.
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 1
#' @author Ahmed El-Gabbas
#' @export

convergence_plot <- function(
    path_coda = NULL, path_model = NULL, env_file = ".env", title = " ",
    n_omega = 1000L, n_cores = 8L, n_RC = c(2L, 2L), beta_n_RC = c(3L, 3L),
    save_plotting_data = TRUE, pages_per_file = 20L, chain_colors = NULL,
    margin_type = "histogram") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(path_coda) || is.null(path_model) || is.null(n_cores)) {
    stop(
      "path_coda, path_model, and n_cores cannot be empty", call. = FALSE)
  }

  if (length(margin_type) != 1) {
    stop("`margin_type` must be a single value.", call. = FALSE)
  }

  if (!margin_type %in% c("histogram", "density")) {
    stop(
      "`margin_type` must be either 'histogram' or 'density'.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpComb <- `2.5%` <- `97.5%` <- Class <- Order <- Family <- DT <- IAS_ID <-
    Species <- Variable <- data <- PlotID <- File <- Page <- Iter <- Value <-
    Chain <- y <- label <- Var_Sp <- CI_025 <- CI_975 <- Var_Min <- Var_Max <-
    Plot_File <- Var_Sp2 <- Species_name <- Species_File <-
    Path_PA <- VarDesc <- Term <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  IASDT.R::cat_time("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("path_coda", "path_model"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_omega", "n_cores", "n_RC"))
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # # Load species summary
  IASDT.R::cat_time("Load species summary")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  SpSummary <- IASDT.R::path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    stop(SpSummary, " file does not exist", call. = FALSE)
  }

  SpSummary <- readr::read_csv(SpSummary, show_col_types = FALSE) %>%
    dplyr::select(Species = Species_name, Species_File)

  # # ..................................................................... ###

  # Create path ------

  IASDT.R::cat_time("Create path")
  Path_Convergence <- dirname(dirname(path_coda)) %>%
    IASDT.R::path("Model_Convergence")
  Pah_Beta_Data <- IASDT.R::path(Path_Convergence, "Beta_Data")
  Path_Convergence_BySp <- IASDT.R::path(Path_Convergence, "Beta_BySpecies")
  fs::dir_create(c(Path_Convergence, Path_Convergence_BySp, Pah_Beta_Data))

  # # ..................................................................... ###

  # Prepare convergence data ------

  IASDT.R::cat_time("Prepare convergence data")

  if (!file.exists(path_model)) {
    stop("`path_model` does not exist", call. = FALSE)
  }

  if (!file.exists(path_coda)) {
    stop("`path_coda` does not exist", call. = FALSE)
  }

  IASDT.R::cat_time("Loading coda object", level = 1)
  Coda_Obj <- IASDT.R::load_as(path_coda)

  IASDT.R::cat_time("Loading fitted model object", level = 1)
  Model <- IASDT.R::load_as(path_model)

  # Model variables
  ModVars <- Model$covNames

  # Number of chains
  NChains <- length(Model$postList)

  # Number of samples
  SampleSize <- Model$samples

  #  Plotting colours
  if (is.null(chain_colors)) {
    chain_colors <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }
  if (length(chain_colors) != NChains) {
    warning(
      "The length of provided colours != number of chains", call. = FALSE)
    chain_colors <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }

  # # ..................................................................... ###

  # Rho ------

  if ("Rho" %in% names(Coda_Obj)) {
    IASDT.R::cat_time("Rho")

    FileConv_Rho <- IASDT.R::path(Path_Convergence, "convergence_rho.RData")

    if (IASDT.R::check_data(FileConv_Rho, warning = FALSE)) {
      IASDT.R::cat_time("Loading plotting data", level = 1)
      PlotObj_Rho <- IASDT.R::load_as(FileConv_Rho)
    } else {
      IASDT.R::cat_time("Prepare plot", level = 1)
      PlotObj_Rho <- IASDT.R::convergence_rho(
        posterior = Coda_Obj, model_object = Model, title = title,
        chain_colors = chain_colors)

      if (save_plotting_data) {
        IASDT.R::cat_time("Save plotting data", level = 1)
        IASDT.R::save_as(
          object = PlotObj_Rho, object_name = "convergence_rho",
          out_path = FileConv_Rho)
      }
    }

    IASDT.R::cat_time("Save plot", level = 1)
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = IASDT.R::path(Path_Convergence, "Convergence_Rho.pdf"),
      width = 18, height = 12, onefile = TRUE)
    plot(PlotObj_Rho)
    grDevices::dev.off()

    rm(PlotObj_Rho, envir = environment())
  }

  # # ..................................................................... ###

  # Alpha  ------

  IASDT.R::cat_time("Alpha")

  FileConv_Alpha <- IASDT.R::path(Path_Convergence, "Convergence_Alpha.RData")

  if (file.exists(FileConv_Alpha)) {
    IASDT.R::cat_time("Loading plotting data", level = 1)
    PlotObj_Alpha <- IASDT.R::load_as(FileConv_Alpha)
  } else {
    IASDT.R::cat_time("Prepare plot", level = 1)
    PlotObj_Alpha <- IASDT.R::convergence_alpha(
      posterior = Coda_Obj, model_object = Model, title = title, n_RC = n_RC,
      add_footer = FALSE, add_title = FALSE, chain_colors = chain_colors)

    if (save_plotting_data) {
      IASDT.R::cat_time("Save plotting data", level = 1)
      IASDT.R::save_as(
        object = PlotObj_Alpha, object_name = "convergence_alpha",
        out_path = IASDT.R::path(Path_Convergence, "Convergence_Alpha.RData"))
    }
  }

  IASDT.R::cat_time("Save plots", level = 1)
  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::cairo_pdf(
    filename = IASDT.R::path(Path_Convergence, "Convergence_Alpha.pdf"),
    width = 18, height = 14, onefile = TRUE)
  print(PlotObj_Alpha)
  grDevices::dev.off()

  Obj_Omega <- Coda_Obj$Omega[[1]]
  Obj_Beta <- Coda_Obj$Beta

  rm(Model, Coda_Obj, PlotObj_Alpha, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Omega  ------

  IASDT.R::cat_time("Omega")

  FileConv_Omega <- IASDT.R::path(Path_Convergence, "Convergence_Omega.qs2")

  if (file.exists(FileConv_Omega)) {
    IASDT.R::cat_time("Loading plotting data", level = 1)
    PlotObj_Omega <- IASDT.R::load_as(FileConv_Omega)
  } else {
    IASDT.R::cat_time("Coda to tibble", level = 1)
    OmegaDF <- IASDT.R::coda_to_tibble(
      coda_object = Obj_Omega, posterior_type = "omega", n_omega = n_omega,
      env_file = env_file)
    invisible(gc())
    SelectedCombs <- unique(OmegaDF$SpComb)

    IASDT.R::cat_time("Prepare confidence interval data", level = 1)
    CI <- purrr::map(.x = Obj_Omega, .f = ~ .x[, SelectedCombs]) %>%
      coda::mcmc.list() %>%
      summary(quantiles = c(0.025, 0.975)) %>%
      magrittr::extract2("quantiles") %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "SpComb") %>%
      stats::setNames(c("SpComb", "CI_25", "CI_975"))

    OmegaDF <- dplyr::left_join(OmegaDF, CI, by = "SpComb")
    OmegaNames <- attr(Obj_Omega[[1]], "dimnames")[[2]]

    IASDT.R::cat_time("Prepare plots", level = 1)
    PlotObj_Omega <- purrr::map_dfr(
      .x = seq_len(n_omega),
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
          paste(collapse = " to ") %>%
          paste0("<b><i>95% credible interval:</i></b> ", .)
        Label_ESS_CI <- data.frame(
          x = -Inf, y = -Inf, label = paste0(Label_ESS, "<br>", Label_CI))

        Label_Panel <- c(CombData$IAS1, CombData$IAS2) %>%
          sort() %>%
          paste0("<i>", ., "</i>") %>%
          paste(collapse = " & <br>") %>%
          paste0(" ") %>%
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
          ggplot2::scale_color_manual(values = chain_colors) +
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

        if (margin_type == "histogram") {
          Plot <- ggExtra::ggMarginal(
            p = Plot, type = margin_type, margins = "y", size = 6,
            color = "steelblue4", fill = "steelblue4", bins = 100)
        } else {
          Plot <- ggExtra::ggMarginal(
            p = Plot, type = margin_type, margins = "y", size = 6,
            color = "steelblue4")
        }
        # Making marginal background matching the plot background
        # https://stackoverflow.com/a/78196022/3652584
        Plot$layout$t[1] <- 1
        Plot$layout$r[1] <- max(Plot$layout$r)

        return(tibble::tibble(SpComb = CombData$SpComb, Plot = list(Plot)))
      }
    )
    IASDT.R::cat_time("Plots preparation is finished", level = 2)

    if (save_plotting_data) {
      IASDT.R::cat_time("Save plot data", level = 1)
      IASDT.R::save_as(object = PlotObj_Omega, out_path = FileConv_Omega)
    }
    rm(OmegaDF, SelectedCombs, CI, OmegaNames, envir = environment())
    invisible(gc())
  }

  IASDT.R::cat_time("Arrange plots", level = 1)
  OmegaPlotList <- tibble::tibble(PlotID = seq_len(nrow(PlotObj_Omega))) %>%
    dplyr::mutate(
      File = ceiling(PlotID / (pages_per_file * n_RC[2] * n_RC[1])),
      Page = ceiling(PlotID / (n_RC[2] * n_RC[1]))) %>%
    tidyr::nest(.by = c("File", "Page"), .key = "PlotID") %>%
    dplyr::mutate(
      PlotID = purrr::map(PlotID, ~ unlist(as.vector(.x))),
      PlotID = purrr::pmap(
        .l = list(File, Page, PlotID),
        .f = function(File, Page, PlotID) {
          PlotTitle <- ggplot2::ggplot() +
            ggplot2::labs(
              title = paste0(
                "Convergence of the omega parameter --- a sample of ",
                n_omega, " species pair"),
              subtitle = paste0(
                "   File ", File, " | Page ",
                (Page - ((File - 1) * pages_per_file)))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              plot.title = ggplot2::element_text(
                face = "bold", size = 20, hjust = 0.5),
              plot.subtitle = ggplot2::element_text(
                size = 12, colour = "grey",
                margin = ggplot2::margin(-5, 0, 0, 0)))

          cowplot::plot_grid(
            plotlist = PlotObj_Omega$Plot[PlotID],
            ncol = n_RC[2], nrow = n_RC[1], align = "hv") %>%
            cowplot::plot_grid(PlotTitle, ., ncol = 1, rel_heights = c(0.05, 1))
        }
      ))

  IASDT.R::cat_time("Save plots", level = 1)
  purrr::walk(
    .x = seq_along(unique(OmegaPlotList$File)),
    .f = ~ {
      invisible({
        CurrPlotOrder <- dplyr::filter(OmegaPlotList, File == .x)
        grDevices::cairo_pdf(
          filename = IASDT.R::path(
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

  IASDT.R::cat_time("Beta")

  FileConv_Beta <- IASDT.R::path(Path_Convergence, "Convergence_Beta.RData")

  if (file.exists(FileConv_Beta)) {

    IASDT.R::cat_time("Loading plotting data", level = 1)
    PlotObj_Beta <- IASDT.R::load_as(FileConv_Beta)

    rm(Obj_Beta, envir = environment())
    invisible(gc())

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
    HTML4 <- "\n&nbsp;&mdash;&nbsp;"
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
          VarDesc, " - L$", paste0(HTML4, "Linear")),
        VarDesc = stringr::str_replace(
          VarDesc, " - Q$", paste0(HTML4, "Quadratic")),
        Variable = paste0(Variable, "_", Term),
        Term = NULL) %>%
      dplyr::bind_rows(LinearTerms)

    IASDT.R::cat_time("Prepare trace plots", level = 1)

    BetaNames <- attr(Obj_Beta[[1]], "dimnames")[[2]]

    IASDT.R::cat_time("Prepare 95% credible interval data", level = 2)
    CI <- summary(Obj_Beta, quantiles = c(0.025, 0.975))$quantiles %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "Var_Sp") %>%
      dplyr::rename(CI_025 = `2.5%`, CI_975 = `97.5%`)

    IASDT.R::cat_time("Coda to tibble", level = 2)
    Beta_DF <- IASDT.R::coda_to_tibble(
      coda_object = Obj_Beta, posterior_type = "beta", env_file = env_file) %>%
      dplyr::left_join(CI, by = "Var_Sp")

    # Variable ranges
    IASDT.R::cat_time("Variable ranges", level = 2)
    VarRanges <- dplyr::arrange(Beta_DF, Variable, IAS_ID) %>%
      dplyr::select(Variable, DT) %>%
      dplyr::mutate(
        Range = purrr::map(.x = DT, .f = ~ range(dplyr::pull(.x, Value)))) %>%
      dplyr::select(-DT) %>%
      tidyr::nest(data = "Range") %>%
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
    IASDT.R::cat_time("Species taxonomy", level = 2)
    SpeciesTaxonomy <- IASDT.R::get_species_name(env_file = env_file) %>%
      dplyr::select(IAS_ID, Class, Order, Family)

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    IASDT.R::cat_time("Preparing data for plotting", level = 2)
    Cols2remove <- c(
      "CI_025", "CI_975", "Var_Min", "Var_Max", "Class", "Order", "Family")
    Beta_DF <- Beta_DF %>%
      dplyr::left_join(VarRanges, by = "Variable") %>%
      dplyr::left_join(SpeciesTaxonomy, by = "IAS_ID") %>%
      dplyr::mutate(
        Var_Sp2 = paste0(Variable, "_", IAS_ID),
        Var_Sp_File = IASDT.R::path(Pah_Beta_Data, paste0(Var_Sp2, ".RData")),
        Plot_File = IASDT.R::path(Pah_Beta_Data, paste0(Var_Sp2, "_Plots.qs2")),
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
          })) %>%
      dplyr::select(-tidyselect::all_of(Cols2remove))

    rm(
      CI, VarRanges, SpeciesTaxonomy, Cols2remove, Obj_Beta,
      envir = environment())
    invisible(gc())

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Prepare working on parallel
    IASDT.R::cat_time("Prepare working on parallel", level = 2)
    IASDT.R::set_parallel(n_cores = min(n_cores, nrow(Beta_DF)), level = 3)
    withr::defer(future::plan("future::sequential", gc = TRUE))

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Split data for each of variables and species combination
    IASDT.R::cat_time(
      "Split data for each of variables and species combination", level = 2)

    Beta_DF <- Beta_DF %>%
      dplyr::mutate(
        Save = furrr::future_pmap(
          .l = list(Var_Sp_File, Var_Sp2, DT),
          .f = function(Var_Sp_File, Var_Sp2, DT) {

            # try saving for a max of 5 attempts using repeat loop
            attempt <- 1
            repeat {

              if (attempt > 5) {
                stop(
                  "Maximum attempts (5) reached without success: ",
                  Var_Sp_File, call. = FALSE)
              }

              try({
                IASDT.R::save_as(
                  object = DT, object_name = Var_Sp2, out_path = Var_Sp_File)
                Sys.sleep(2)
              },
              silent = TRUE)

              if (IASDT.R::check_data(Var_Sp_File, warning = FALSE)) {
                break
              }

              # Increment attempt counter
              attempt <- attempt + 1
            }
          },
          .options = furrr::furrr_options(
            seed = TRUE, packages = c("IASDT.R", "tibble"))),
        Save = NULL,
        DT = NULL)

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Prepare plots
    IASDT.R::cat_time("Prepare plots", level = 2)

    PlotObj_Beta <- future.apply::future_lapply(
      X = seq_len(nrow(Beta_DF)),
      FUN = function(x) {

        Var_Sp <- Beta_DF$Var_Sp[x]
        Species <- Beta_DF$Species[x]
        Curr_IAS <- Beta_DF$IAS_ID[x]
        Var_Sp_File <- Beta_DF$Var_Sp_File[x]
        Plot_File <- Beta_DF$Plot_File[x]

        # check if input data exists
        if (isFALSE(IASDT.R::check_data(Var_Sp_File, warning = FALSE))) {
          stop("File ", x, ": ", Var_Sp_File, " does not exist.", call. = FALSE)
        }

        # Check if the output file already exists
        if (IASDT.R::check_data(Plot_File, warning = FALSE)) {
          return(tibble::tibble(Var_Sp = Var_Sp, Plot_File = Plot_File))
        }

        # delete file if corrupted
        if (file.exists(Plot_File)) {
          IASDT.R::system_command(
            command = paste0("rm -f ", Plot_File), R_object = FALSE,
            ignore.stdout = TRUE)
        }

        attempt <- 1

        repeat {

          if (attempt > 5) {
            stop(
              "Maximum attempts (5) reached without success: ", Var_Sp_File,
              call. = FALSE)
          }

          try({

            DT_all <- IASDT.R::load_as(Var_Sp_File)
            if (is.null(DT_all) || !is.list(DT_all)) {
              stop("Loaded data is invalid for file: ", Var_Sp_File)
            }

            DT_all$Post <- NULL
            invisible(gc())


            ## Gelman convergence diagnostic
            Label_Gelman <- round(DT_all$Gelman$psrf, 3) %>%
              paste(collapse = " / ") %>%
              paste0("<b><i>Gelman convergence diagnostic:</i></b> ", .) %>%
              data.frame(x = Inf, y = -Inf, label = .)

            ## Effective sample size / CI
            Label_ESS <- round(DT_all$ESS / NChains) %>%
              paste0(
                "<b><i>Mean effective sample size:</i></b> ",
                ., " / ", SampleSize)
            CurrCI <- c(DT_all$CI_025, DT_all$CI_975)
            Label_CI <- paste(round(CurrCI, 4), collapse = " to ") %>%
              paste0("<b><i>95% credible interval:</i></b> ", .)
            Label_ESS_CI <- data.frame(
              x = -Inf, y = -Inf, label = paste0(Label_ESS, "<br>", Label_CI))

            Label_Panel <- data.frame(
              x = Inf, y = Inf,
              label = paste0("<br><b><i>", Species, "</i></b>"))

            PanelTitle <- c(DT_all$Class, DT_all$Order, DT_all$Family) %>%
              paste(collapse = " | ") %>%
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
              ggplot2::scale_color_manual(values = chain_colors) +
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
                hjust = 1, vjust = 1, lineheight = 0, fill = NA,
                label.color = NA) +
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


            if (margin_type == "histogram") {
              Plot_Marginal <- ggExtra::ggMarginal(
                p = Plot, type = margin_type, margins = "y", size = 6,
                color = "steelblue4", fill = "steelblue4", bins = 100)
            } else {
              Plot_Marginal <- ggExtra::ggMarginal(
                p = Plot, type = margin_type, margins = "y", size = 6,
                color = "steelblue4")
            }

            # Making marginal background matching the plot background
            # https://stackoverflow.com/a/78196022/3652584
            Plot_Marginal$layout$t[1] <- 1
            Plot_Marginal$layout$r[1] <- max(Plot_Marginal$layout$r)

            suppressWarnings({
              if (margin_type == "histogram") {
                Plot2_Marginal <- ggExtra::ggMarginal(
                  p = Plot2, type = margin_type, margins = "y", size = 6,
                  color = "steelblue4", fill = "steelblue4", bins = 100)
              } else {
                Plot2_Marginal <- ggExtra::ggMarginal(
                  p = Plot2, type = margin_type, margins = "y", size = 6,
                  color = "steelblue4")
              }
            })

            # Making marginal background matching the plot background
            # https://stackoverflow.com/a/78196022/3652584
            Plot2_Marginal$layout$t[1] <- 1
            Plot2_Marginal$layout$r[1] <- max(Plot2_Marginal$layout$r)

            IASDT.R::save_as(
              object = list(
                Plot = Plot, Plot_Marginal = Plot_Marginal,
                PlotFixedY_Marginal = Plot2_Marginal),
              out_path = Plot_File)

            Sys.sleep(2)

          },
          silent = TRUE)

          if (IASDT.R::check_data(Plot_File, warning = FALSE)) {
            break
          }

          # Increment attempt counter
          attempt <- attempt + 1
        }

        # Return result if successful
        tibble::tibble(Var_Sp = Var_Sp, Plot_File = Plot_File)

      },
      future.seed = TRUE,
      future.globals = c(
        "Beta_DF", "NChains", "SampleSize", "chain_colors", "margin_type"),
      future.packages = c(
        "dplyr", "ggplot2", "ggtext", "magrittr", "stringr", "ggExtra",
        "coda", "IASDT.R", "qs2", "tibble"))

    PlotObj_Beta <- dplyr::left_join(Beta_DF, VarsDesc, by = "Variable")

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Stopping cluster
    IASDT.R::set_parallel(stop = TRUE, level = 2)

    rm(Beta_DF, BetaNames, envir = environment())
    invisible(gc())

    if (save_plotting_data) {
      IASDT.R::cat_time("Save trace plot data", level = 2)
      IASDT.R::save_as(
        object = PlotObj_Beta, object_name = "Convergence_Beta",
        out_path = FileConv_Beta)
    }
  }

  # # ..................................................................... ###

  # Beta - 2. by variable ------

  IASDT.R::cat_time("Trace plots, grouped by variables", level = 1)

  IASDT.R::cat_time("Preparing data", level = 2)
  BetaTracePlots_ByVar <- dplyr::arrange(PlotObj_Beta, Variable, IAS_ID) %>%
    dplyr::select(Variable, Plot_File, VarDesc) %>%
    tidyr::nest(Plot_File = "Plot_File") %>%
    dplyr::mutate(Plot_File = purrr::map(Plot_File, unlist))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Prepare working on parallel
  IASDT.R::set_parallel(
    n_cores = min(n_cores, nrow(BetaTracePlots_ByVar)), level = 2)
  withr::defer(future::plan("future::sequential", gc = TRUE))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  IASDT.R::cat_time("Save plots", level = 2)
  VarNames <- BetaTracePlots_ByVar$Variable

  BetaTracePlots_ByVar0 <- future.apply::future_lapply(
    X = VarNames,
    FUN = function(x) {

      VarDesc <- BetaTracePlots_ByVar %>%
        dplyr::filter(Variable == x) %>%
        dplyr::pull(VarDesc)

      Plots <- dplyr::filter(BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(Plot_File) %>%
        magrittr::extract2(1) %>%
        purrr::map(
          .f = ~ {
            IASDT.R::load_as(.x) %>%
              magrittr::extract(c("Plot_Marginal", "PlotFixedY_Marginal"))
          })

      BetaPlots <- purrr::map(
        .x = Plots, .f = magrittr::extract2, "Plot_Marginal")
      BetaPlotsFixedY <- purrr::map(
        .x = Plots, .f = magrittr::extract2, "PlotFixedY_Marginal")
      rm(Plots, envir = environment())
      invisible(gc())

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
        dplyr::mutate(Page = ceiling(PlotID / (n_RC[2] * n_RC[1]))) %>%
        tidyr::nest(.by = "Page", .key = "PlotID") %>%
        dplyr::mutate(
          Plot = purrr::map(
            .x = PlotID,
            .f = ~ {
              ID <- unlist(as.vector(.x))
              cowplot::plot_grid(
                plotlist = BetaPlots[ID],
                ncol = n_RC[2], nrow = n_RC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitle, ., ncol = 1, rel_heights = c(0.035, 1))
            }),
          PlotFixedY = purrr::map(
            .x = PlotID,
            .f = ~ {
              ID <- unlist(as.vector(.x))
              cowplot::plot_grid(
                plotlist = BetaPlotsFixedY[ID],
                ncol = n_RC[2], nrow = n_RC[1], align = "hv") %>%
                cowplot::plot_grid(
                  PlotTitleFixed, .,
                  ncol = 1, rel_heights = c(0.035, 1))
            }))

      invisible({
        grDevices::cairo_pdf(
          filename = IASDT.R::path(
            Path_Convergence, paste0("Convergence_Beta_", x, "_FreeY.pdf")),
          width = 18, height = 13, onefile = TRUE)
        purrr::walk(BetaPlotList$Plot, grid::grid.draw, recording = FALSE)
        grDevices::dev.off()

        grDevices::cairo_pdf(
          filename = IASDT.R::path(
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
    future.globals = c("BetaTracePlots_ByVar", "n_RC", "Path_Convergence"),
    future.packages = c(
      "tidyr", "dplyr", "ggplot2", "purrr", "ggtext",
      "tibble", "cowplot", "grDevices", "IASDT.R"))

  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Stopping cluster
  IASDT.R::set_parallel(stop = TRUE, level = 2)

  # # ..................................................................... ###

  # Beta - 3. by species ------

  IASDT.R::cat_time("Trace plots, grouped by species", level = 1)

  IASDT.R::cat_time("Preparing data", level = 2)
  Order <- stringr::str_remove_all(ModVars, "\\(|\\)")
  BetaTracePlots_BySp <- PlotObj_Beta %>%
    dplyr::arrange(Species, factor(Variable, levels = Order)) %>%
    dplyr::select(Species, IAS_ID, Plot_File, Variable, VarDesc) %>%
    tidyr::nest(data = -c("Species", "IAS_ID")) %>%
    dplyr::left_join(SpSummary, by = "Species")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Prepare working on parallel
  IASDT.R::set_parallel(
    n_cores = min(n_cores, nrow(BetaTracePlots_BySp)), level = 2)
  withr::defer(future::plan("future::sequential", gc = TRUE))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  IASDT.R::cat_time("Save plots", level = 2)

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
              Plot <- IASDT.R::load_as(.x)$Plot +
                ggplot2::ggtitle(.y) +
                ggplot2::theme(
                  plot.title = ggtext::element_markdown(
                    hjust = 0.5, size = 15, colour = "red",
                    margin = ggplot2::margin(0, 0, -2.5, 0)))

              if (margin_type == "histogram") {
                Plot <- ggExtra::ggMarginal(
                  p = Plot, type = margin_type, margins = "y", size = 6,
                  color = "steelblue4", fill = "steelblue4", bins = 100)
              } else {
                Plot <- ggExtra::ggMarginal(
                  p = Plot, type = margin_type, margins = "y", size = 6,
                  color = "steelblue4")
              }
              Plot$layout$t[1] <- 1
              Plot$layout$r[1] <- max(Plot$layout$r)
              return(Plot)
            })) %>%
        dplyr::pull("Plot")

      NumPages <- ceiling(length(SpPlots) / (beta_n_RC[1] * beta_n_RC[2]))

      SpPlots2 <- IASDT.R::split_vector(
        vector = seq_len(length(SpPlots)), n_splits = NumPages) %>%
        purrr::map(
          .f = ~ {
            SpPlots[.x] %>%
              cowplot::plot_grid(
                plotlist = ., ncol = beta_n_RC[2],
                nrow = beta_n_RC[1], align = "hv") %>%
              cowplot::plot_grid(
                PlotTitle, ., ncol = 1, rel_heights = c(0.03, 1))
          })

      invisible({
        grDevices::cairo_pdf(
          filename = IASDT.R::path(
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
    future.seed = TRUE,
    future.globals = c(
      "margin_type", "BetaTracePlots_BySp", "Path_Convergence_BySp",
      "beta_n_RC"),
    future.packages = c(
      "dplyr", "coda", "ggplot2", "ggExtra", "ggtext", "IASDT.R",
      "stringr", "gtools", "cowplot", "purrr", "grDevices"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Stopping cluster
  IASDT.R::set_parallel(stop = TRUE, level = 2)

  rm(BetaTracePlots_BySp0, envir = environment())

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "Plot model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
