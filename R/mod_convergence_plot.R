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
#' @param model_dir Character. A path to the model directory.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param title Character. title for **rho** and **alpha** convergence plots.
#'   Default: " "
#' @param n_omega Integer. Number of species interactions sampled for Omega
#'   parameter diagnostics. Default: 1000L
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param n_RC Numeric vector. Grid layout (rows&times;columns) for arranging
#'   alpha parameter plots. Default: `c(2, 2)`. If `NULL`, the layout is
#'   automatically determined based on the number of alpha levels.
#' @param beta_n_RC Numeric vector. The grid layout (rows&times;columns) for
#'   arranging beta parameter plots. Default: `c(3, 3)`.
#' @param pages_per_file Integer. Number of plots per page in the Omega
#'   parameter output. Default: 20L.
#' @param chain_colors Character vector. MCMC chain colours (optional). Default:
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
#' @param spatial_model Logical. Whether the model is a spatial model. If `TRUE`
#'   (default), the function will generate additional plots for the model's
#'   `Alpha` parameter.
#'
#' @details `convergence_alpha()`, `convergence_rho()`, and
#'   `convergence_beta_ranges` are internal functions and should not be called
#'   directly. The `convergence_beta_ranges` plots the convergence range of the
#'   each species beta parameters. It can be used to check if any of the chains
#'   show convergence issues; i.e., showing exceptionally high or low beta
#'   values.
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 1
#' @author Ahmed El-Gabbas
#' @export

convergence_plot <- function(
    path_coda = NULL, path_model = NULL, env_file = ".env", title = " ",
    n_omega = 1000L, n_cores = 8L, strategy = "multisession",
    n_RC = c(2L, 2L), beta_n_RC = c(3L, 3L), pages_per_file = 20L,
    chain_colors = NULL, margin_type = "histogram", spatial_model = TRUE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  if (is.null(path_coda)) {
    ecokit::stop_ctx(
      "`path_coda` cannot be empty",
      path_coda = path_coda, include_backtrace = TRUE)
  }
  if (is.null(path_model)) {
    ecokit::stop_ctx(
      "path_model cannot be empty",
      path_model = path_model, include_backtrace = TRUE)
  }

  if (length(margin_type) != 1) {
    ecokit::stop_ctx(
      "`margin_type` must be a single value.",
      margin_type = margin_type, length_margin_type = length(margin_type),
      include_backtrace = TRUE)
  }

  if (!margin_type %in% c("histogram", "density")) {
    ecokit::stop_ctx(
      "`margin_type` must be either 'histogram' or 'density'.",
      margin_type = margin_type, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpComb <- `2.5%` <- `97.5%` <- Class <- Order <- Family <- DT <- IAS_ID <-
    Species <- Variable <- data <- PlotID <- File <- Page <- Iter <- Value <-
    Chain <- y <- label <- Var_Sp <- CI_025 <- CI_975 <- Var_Min <- Var_Max <-
    Plot_File <- Var_Sp2 <- Species_name <- Species_File <- Path_PA <-
    VarDesc <- is_intercept <- ModVar <- LQ <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::cat_time("Check input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("path_coda", "path_model", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_omega", "n_cores", "n_RC"))
  rm(AllArgs, envir = environment())

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "ggplot2", "ggtext", "magrittr", "stringr", "ggExtra",
      "coda", "ecokit", "qs2", "tibble", "tidyr", "purrr", "cowplot",
      "gtools", "withr", "fs"),
    strategy = strategy)

  # # ..................................................................... ###

  # # Load species summary
  ecokit::cat_time("Load species summary")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  SpSummary <- fs::path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    ecokit::stop_ctx(
      "SpSummary file does not exist", SpSummary = SpSummary,
      include_backtrace = TRUE)
  }

  SpSummary <- readr::read_csv(
    file = SpSummary, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(Species = Species_name, Species_File)

  # # ..................................................................... ###

  # Create path ------

  ecokit::cat_time("Create path")
  Path_Convergence <- dirname(dirname(path_coda)) %>%
    fs::path("Model_Convergence")
  Path_Beta_Data <- fs::path(Path_Convergence, "Beta_Data")
  Path_Convergence_BySp <- fs::path(Path_Convergence, "Beta_BySpecies")
  fs::dir_create(c(Path_Convergence, Path_Convergence_BySp, Path_Beta_Data))

  # # ..................................................................... ###

  # Prepare convergence data ------

  ecokit::cat_time("Prepare convergence data")

  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "`path_model` does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  if (!file.exists(path_coda)) {
    ecokit::stop_ctx(
      "`path_coda` does not exist", path_coda = path_coda,
      include_backtrace = TRUE)
  }

  ecokit::cat_time("Loading coda object", level = 1L)
  Coda_Obj <- ecokit::load_as(path_coda)
  names_coda <- names(Coda_Obj)

  ecokit::cat_time("Loading fitted model object", level = 1L)
  Model <- ecokit::load_as(path_model)

  # Model variables
  ModVars <- Model$covNames

  # Number of chains
  NChains <- length(Model$postList)

  # Number of samples
  SampleSize <- Model$samples

  #  Plotting colours
  define_chain_colors <- FALSE

  if (is.null(chain_colors)) {
    define_chain_colors <- TRUE
  } else if (length(chain_colors) != NChains) {
    define_chain_colors <- TRUE
    warning(
      "The length of provided colours != number of chains", call. = FALSE)
  }

  if (define_chain_colors) {
    # minimum value of n colours in RColorBrewer::brewer.pal is 3.
    # black and grey will be used anyway
    if (NChains >= 4) {
      chain_colors <- c(
        "black", "grey60",
        RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
    } else if (NChains == 3) {
      chain_colors <- c("black", "grey60", "red")
    } else if (NChains == 2) {
      chain_colors <- c("black", "grey60")
    }
  }

  # # ..................................................................... ###

  # Rho ------

  if ("Rho" %in% names_coda) {

    ecokit::cat_time("Rho")
    FileConv_Rho <- fs::path(Path_Convergence, "convergence_rho.RData")

    if (ecokit::check_data(FileConv_Rho, warning = FALSE)) {
      ecokit::cat_time("Loading plotting data", level = 1L)
      PlotObj_Rho <- ecokit::load_as(FileConv_Rho)
    } else {
      ecokit::cat_time("Prepare plot", level = 1L)
      PlotObj_Rho <- IASDT.R::convergence_rho(
        posterior = Coda_Obj, model_object = Model, title = title,
        chain_colors = chain_colors)

      ecokit::cat_time("Save plotting data", level = 1L)
      ecokit::save_as(
        object = PlotObj_Rho, object_name = "convergence_rho",
        out_path = FileConv_Rho)
    }

    ecokit::cat_time("Save plot", level = 1L)
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(Path_Convergence, "Convergence_Rho.pdf"),
      width = 18, height = 12, onefile = TRUE)
    plot(PlotObj_Rho)
    grDevices::dev.off()

    rm(PlotObj_Rho, envir = environment())
  }

  # # ..................................................................... ###

  # Alpha  ------

  if (("Alpha" %in% names_coda) && spatial_model) {

    ecokit::cat_time("Alpha")
    FileConv_Alpha <- fs::path(Path_Convergence, "Convergence_Alpha.RData")

    # Ensure that all latent factors of the model are plotted
    n_lf <- ncol(Coda_Obj$Alpha[[1]][[1]])
    n_RC_alpha <- n_RC
    if ((n_RC_alpha[1] * n_RC_alpha[2]) < n_lf) {
      if (n_RC_alpha[1] == 1) {
        n_RC_alpha[1] <- 2
      } else if (n_RC_alpha[1] == 2) {
        n_RC_alpha[1] <- 3
      }
    }

    if (ecokit::check_data(FileConv_Alpha, warning = FALSE)) {
      ecokit::cat_time("Loading plotting data", level = 1L)
      PlotObj_Alpha <- ecokit::load_as(FileConv_Alpha)
    } else {
      ecokit::cat_time("Prepare plot", level = 1L)
      PlotObj_Alpha <- IASDT.R::convergence_alpha(
        posterior = Coda_Obj, model_object = Model, title = title,
        n_RC = n_RC_alpha, add_footer = FALSE, add_title = FALSE,
        chain_colors = chain_colors)

      ecokit::cat_time("Save plotting data", level = 1L)
      ecokit::save_as(
        object = PlotObj_Alpha, object_name = "convergence_alpha",
        out_path = FileConv_Alpha)
    }

    ecokit::cat_time("Save plots", level = 1L)
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(Path_Convergence, "Convergence_Alpha.pdf"),
      width = 18, height = 14, onefile = TRUE)
    print(PlotObj_Alpha)
    grDevices::dev.off()

    Obj_Omega <- Coda_Obj$Omega[[1]]
    Obj_Beta <- Coda_Obj$Beta

    rm(Model, Coda_Obj, PlotObj_Alpha, envir = environment())
    invisible(gc())
  }

  # # ..................................................................... ###

  # Omega  ------

  if ("Omega" %in% names_coda) {

    ecokit::cat_time("Omega")

    FileConv_Omega <- fs::path(Path_Convergence, "Convergence_Omega.qs2")

    if (ecokit::check_data(FileConv_Omega, warning = FALSE)) {

      ecokit::cat_time("Loading plotting data", level = 1L)
      PlotObj_Omega <- ecokit::load_as(FileConv_Omega)

    } else {

      ecokit::cat_time("Coda to tibble", level = 1L)
      OmegaDF <- IASDT.R::coda_to_tibble(
        coda_object = Obj_Omega, posterior_type = "omega", n_omega = n_omega,
        env_file = env_file)
      invisible(gc())
      SelectedCombs <- unique(OmegaDF$SpComb)

      ecokit::cat_time("Prepare confidence interval data", level = 1L)
      CI <- purrr::map(.x = Obj_Omega, .f = ~ .x[, SelectedCombs]) %>%
        coda::mcmc.list() %>%
        summary(quantiles = c(0.025, 0.975)) %>%
        magrittr::extract2("quantiles") %>%
        as.data.frame() %>%
        tibble::as_tibble(rownames = "SpComb") %>%
        stats::setNames(c("SpComb", "CI_25", "CI_975"))

      OmegaDF <- dplyr::left_join(OmegaDF, CI, by = "SpComb")
      OmegaNames <- attr(Obj_Omega[[1]], "dimnames")[[2]]

      # Prepare omega plots
      ecokit::cat_time("Prepare omega plots", level = 1L)

      PlotObj_Omega <- purrr::map_dfr(
        .x = seq_len(n_omega),
        .f = function(x) {

          temp_file <- fs::file_temp(
            pattern = paste0("plot_", x, "_"), ext = "pdf")
          grDevices::cairo_pdf(temp_file)
          on.exit({
            grDevices::dev.off()
            try(fs::file_delete(temp_file), silent = TRUE)
          },
          add = TRUE)

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
              "<b><i>Mean effective sample size:</i></b> ",
              ., " / ", SampleSize)
          CurrCI <- dplyr::select(CombData, c("CI_25", "CI_975")) %>%
            unlist() %>%
            round(2)
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
              hjust = 1, vjust = 1, lineheight = 0, fill = NA,
              label.color = NA) +
            ggtext::geom_richtext(
              mapping = ggplot2::aes(x = x, y = y, label = label),
              data = Label_Gelman, inherit.aes = FALSE, size = 6, hjust = 0,
              vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
            ggtext::geom_richtext(
              mapping = ggplot2::aes(x = x, y = y, label = label),
              data = Label_ESS_CI, inherit.aes = FALSE, size = 6, hjust = 0,
              vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
            ggplot2::labs(x = NULL, y = NULL) +
            ggplot2::theme_bw() +
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

      ecokit::cat_time("Save plot data", level = 1L)
      ecokit::save_as(object = PlotObj_Omega, out_path = FileConv_Omega)

      rm(OmegaDF, SelectedCombs, CI, OmegaNames, envir = environment())
      invisible(gc())
    }


    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    ecokit::cat_time("Arrange plots", level = 1L)
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
                  "Convergence of the omega parameter ---  a sample of ",
                  n_omega, " species pairs"),
                subtitle = paste0(
                  "   File ", File, " | Page ",
                  (Page - ((File - 1) * pages_per_file)))) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                text = ggplot2::element_text(family = "sans"),
                plot.title = ggtext::element_markdown(
                  face = "bold", size = 20, hjust = 0.5),
                plot.subtitle = ggplot2::element_text(
                  size = 12, colour = "grey",
                  margin = ggplot2::margin(-5, 0, 0, 0)))

            cowplot::plot_grid(
              plotlist = PlotObj_Omega$Plot[PlotID],
              ncol = n_RC[2], nrow = n_RC[1], align = "hv") %>%
              cowplot::plot_grid(
                PlotTitle, ., ncol = 1, rel_heights = c(0.05, 1))

          }
        ))

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    ecokit::cat_time("Save plots", level = 1L)
    ecokit::cat_time(
      paste0(
        "Saving omega convergence plots for ", n_omega,
        " species associations to ", length(unique(OmegaPlotList$File)),
        " pdf files"),
      level = 2L, cat_timestamp = FALSE)

    purrr::walk(
      .x = seq_along(unique(OmegaPlotList$File)),
      .f = ~ {
        invisible({
          CurrPlotOrder <- dplyr::filter(OmegaPlotList, File == .x)
          grDevices::cairo_pdf(
            filename = fs::path(
              Path_Convergence, paste0("Convergence_Omega_", .x, ".pdf")),
            width = 18, height = 14, onefile = TRUE)
          purrr::map(CurrPlotOrder$PlotID, grid::grid.draw, recording = FALSE)
          grDevices::dev.off()
        })
      })

    rm(OmegaPlotList, PlotObj_Omega, Obj_Omega, envir = environment())
    invisible(gc())

  }

  # # ..................................................................... ###

  # Beta - 1. Prepare data ------

  ecokit::cat_time("Beta")


  arrange_vars <- function(df) {
    ecokit::arrange_alphanum(df, Variable) %>%
      dplyr::mutate(is_intercept = (Variable == "Intercept")) %>%
      dplyr::arrange(dplyr::desc(is_intercept)) %>%
      dplyr::select(-is_intercept)
  }

  FileConv_Beta <- fs::path(Path_Convergence, "Convergence_Beta.qs2")

  if (ecokit::check_data(FileConv_Beta, warning = FALSE, n_threads = 5)) {

    ecokit::cat_time("Loading plotting data", level = 1L)
    PlotObj_Beta <- ecokit::load_as(FileConv_Beta, n_threads = 5)

    rm(Obj_Beta, envir = environment())
    invisible(gc())

  } else {

    VarsDesc <- tibble::tribble(
      ~Variable, ~VarDesc,
      "(Intercept)", "Intercept",
      "RiversLog", "River length (log<sub>10</sub>)",
      "RoadRailLog", "Road + Rail intensity (log<sub>10</sub>)",
      "EffortsLog", "Sampling efforts (log<sub>10</sub>)",
      "HabLog", "% Habitat coverage",
      "bio1", "bio1 (annual mean temperature)",
      "bio2", "bio2 (mean diurnal range)",
      "bio3", "bio3 (isothermality)",
      "bio4", "bio4 (temperature seasonality)",
      "bio5", "bio5 (max temperature of warmest month)",
      "bio6", "bio6 (temperature of the coldest month)",
      "bio7", "bio7 (temperature annual range)",
      "bio8", "bio8 (temperatures of the wettest quarter)",
      "bio9", "bio9 (mean temperature of driest quarter)",
      "bio10", "bio10 (mean temperature of warmest quarter)",
      "bio11", "bio11 (mean temperature of coldest quarter)",
      "bio12", "bio12 (annual precipitation amount)",
      "bio13", "bio13 (precipitation of wettest month)",
      "bio14", "bio14 (precipitation of driest month)",
      "bio15", "bio15 (precipitation seasonality)",
      "bio16", "bio16 (precipitation of wettest quarter)",
      "bio17", "bio17 (precipitation of driest quarter)",
      "bio18", "bio18 (precipitation of the warmest quarter)",
      "bio19", "bio19 (precipitation of coldest quarter)",
      "npp", "npp (net primary productivity)",
      "soil", "soil bulk density",
      "wetness", "topographic wetness index")

    VarsDesc <- tibble::tibble(ModVar = ModVars) %>%
      dplyr::mutate(
        Variable = stringr::str_remove_all(
          ModVar, "stats::poly\\(|, degree = [0-9], raw = TRUE\\)[0-9]"),
        LQ = dplyr::case_when(
          stringr::str_detect(ModVar, "raw = TRUE\\)1") ~ "L",
          stringr::str_detect(ModVar, "raw = TRUE\\)2") ~ "Q",
          .default = "L_only")) %>%
      dplyr::left_join(VarsDesc, by = "Variable") %>%
      dplyr::mutate(
        Variable = purrr::map2_chr(
          .x = Variable, .y = LQ,
          .f = ~ {
            dplyr::case_when(
              .y == "L" ~ paste0(.x, "_L"),
              .y == "Q" ~ paste0(.x, "_Q"),
              .default = .x)
          }),
        VarDesc = purrr::map2_chr(
          .x = VarDesc, .y = LQ,
          .f = ~{
            desc <- stringr::str_to_sentence(.x)
            desc <- dplyr::case_when(
              (.y == "L_only" && .x != "Intercept") ~
                paste0(desc, "\n&nbsp;&mdash;&nbsp;Linear"),
              .y == "L" ~ paste0(desc, "\n&nbsp;&mdash;&nbsp;Linear"),
              .y == "Q" ~ paste0(desc, "\n&nbsp;&mdash;&nbsp;Quadratic"),
              .default = .x)

            stringr::str_glue(
              "<span style='color:blue;'><b>{desc}</b></span>")
          }),
        Variable = stringr::str_replace(
          Variable, "\\(Intercept\\)", "Intercept"),
        LQ = NULL) %>%
      arrange_vars()

    ecokit::cat_time("Prepare trace plots", level = 1L)

    BetaNames <- attr(Obj_Beta[[1]], "dimnames")[[2]]

    ecokit::cat_time("Prepare 95% credible interval data", level = 2L)
    CI <- summary(Obj_Beta, quantiles = c(0.025, 0.975))$quantiles %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "Var_Sp") %>%
      dplyr::rename(CI_025 = `2.5%`, CI_975 = `97.5%`)

    ecokit::cat_time("Coda to tibble", level = 2L)
    Beta_DF <- IASDT.R::coda_to_tibble(
      coda_object = Obj_Beta, posterior_type = "beta", env_file = env_file) %>%
      dplyr::left_join(CI, by = "Var_Sp")

    # Variable ranges
    ecokit::cat_time("Variable ranges", level = 2L)
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
      ecokit::arrange_alphanum(Variable) %>%
      tidyr::unnest_wider("Range") %>%
      arrange_vars()

    # Species taxonomy
    ecokit::cat_time("Species taxonomy", level = 2L)
    SpeciesTaxonomy <- IASDT.R::get_species_name(env_file = env_file) %>%
      dplyr::select(IAS_ID, Class, Order, Family)

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    ecokit::cat_time("Preparing plotting data", level = 2L)
    Cols2remove <- c(
      "CI_025", "CI_975", "Var_Min", "Var_Max", "Class", "Order", "Family")

    Beta_DF <- Beta_DF %>%
      dplyr::left_join(VarRanges, by = "Variable") %>%
      dplyr::left_join(SpeciesTaxonomy, by = "IAS_ID") %>%
      dplyr::mutate(
        Var_Sp2 = paste0(Variable, "_", IAS_ID),
        Var_Sp_File = fs::path(Path_Beta_Data, paste0(Var_Sp2, ".RData")),
        Plot_File = fs::path(Path_Beta_Data, paste0(Var_Sp2, "_Plots.qs2")),
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

    # Prepare working in parallel
    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, nrow(Beta_DF)), level = 2L,
        strategy = strategy, cat_timestamp = FALSE, future_max_size = 1500L)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Split data for each of variables and species combination
    ecokit::cat_time(
      "Split data for each of variables and species combination", level = 2L)

    Beta_DF2 <- future.apply::future_lapply(
      X = seq_len(nrow(Beta_DF)),
      FUN = function(x) {

        Var_Sp_File <- Beta_DF$Var_Sp_File[[x]]

        if (ecokit::check_data(Var_Sp_File, warning = FALSE)) {
          return(NULL)
        }

        # try saving for a max of 5 attempts using repeat loop
        attempt <- 1
        repeat {

          if (attempt > 5) {
            ecokit::stop_ctx(
              "Maximum attempts (5) reached without success: ",
              Var_Sp_File = Var_Sp_File, include_backtrace = TRUE)
          }

          try({
            ecokit::save_as(
              object = Beta_DF$DT[[x]], object_name = Beta_DF$Var_Sp2[[x]],
              out_path = Var_Sp_File)
            Sys.sleep(2)
          },
          silent = TRUE)

          if (ecokit::check_data(Var_Sp_File, warning = FALSE)) {
            break
          }

          # Increment attempt counter
          attempt <- attempt + 1
        }
      },
      future.seed = TRUE, future.globals = "Beta_DF",
      future.packages = pkg_to_export)

    rm(Beta_DF2, envir = environment())

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Prepare beta plots
    ecokit::cat_time("Prepare beta plots", level = 2L)

    PlotObj_Beta <- future.apply::future_lapply(
      X = seq_len(nrow(Beta_DF)),
      FUN = function(x) {

        Var_Sp <- Beta_DF$Var_Sp[x]
        Species <- Beta_DF$Species[x]
        Curr_IAS <- Beta_DF$IAS_ID[x]
        Var_Sp_File <- Beta_DF$Var_Sp_File[x]
        Plot_File <- Beta_DF$Plot_File[x]

        # check if input data exists
        if (isFALSE(ecokit::check_data(Var_Sp_File, warning = FALSE))) {
          ecokit::stop_ctx(
            "File does not exist.", x = x, Var_Sp_File = Var_Sp_File,
            include_backtrace = TRUE)
        }

        # Check if the output file already exists
        if (ecokit::check_data(Plot_File, warning = FALSE)) {
          return(tibble::tibble(Var_Sp = Var_Sp, Plot_File = Plot_File))
        }

        # delete file if corrupted
        if (file.exists(Plot_File)) {
          fs::file_delete(Plot_File)
        }

        temp_file <- fs::file_temp(ext = "pdf")
        grDevices::cairo_pdf(temp_file)

        on.exit({
          grDevices::dev.off()
          try(fs::file_delete(temp_file), silent = TRUE)
        },
        add = TRUE)

        attempt <- 1

        repeat {

          if (attempt > 5) {
            ecokit::stop_ctx(
              "Maximum attempts (5) reached without success: ",
              Var_Sp_File = Var_Sp_File, include_backtrace = TRUE)
          }

          try({

            DT_all <- ecokit::load_as(Var_Sp_File)
            if (is.null(DT_all) || !is.list(DT_all)) {
              ecokit::stop_ctx(
                "Loaded data is invalid", Var_Sp_File = Var_Sp_File,
                DT_all = DT_all, class_DT_all = DT_all,
                include_backtrace = TRUE)
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
              ggplot2::labs(x = NULL, y = NULL) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                text = ggplot2::element_text(family = "sans"),
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

            ecokit::save_as(
              object = list(
                Plot = Plot, Plot_Marginal = Plot_Marginal,
                PlotFixedY_Marginal = Plot2_Marginal),
              out_path = Plot_File)

            Sys.sleep(2)

          },
          silent = TRUE)

          if (ecokit::check_data(Plot_File, warning = FALSE)) {
            break
          }

          # Increment attempt counter
          attempt <- attempt + 1
        }

        # Return result if successful
        tibble::tibble(Var_Sp = Var_Sp, Plot_File = Plot_File)

      },
      future.seed = TRUE, future.packages = pkg_to_export,
      future.globals = c(
        "Beta_DF", "NChains", "SampleSize", "chain_colors", "margin_type")) %>%
      ecokit::quiet_device()

    PlotObj_Beta <- dplyr::left_join(Beta_DF, VarsDesc, by = "Variable")

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    # Stopping cluster
    if (n_cores > 1) {
      ecokit::set_parallel(
        stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
      future::plan("sequential", gc = TRUE)
    }

    rm(Beta_DF, BetaNames, envir = environment())
    invisible(gc())

    ecokit::cat_time("Save trace plot data", level = 2L)
    ecokit::save_as(
      object = PlotObj_Beta, object_name = "Convergence_Beta",
      out_path = FileConv_Beta, n_threads = 5)

  }

  # # ..................................................................... ###

  # plot minimum and maximum value of each beta parameter
  ecokit::cat_time(
    "plot minimum and maximum value of each beta parameter", level = 1L)

  IASDT.R::convergence_beta_ranges(model_dir = dirname(path_model))

  # # ..................................................................... ###

  # Beta - 2. by variable ------

  ecokit::cat_time("Trace plots, grouped by variables", level = 1L)

  ecokit::cat_time("Preparing data", level = 2L)
  BetaTracePlots_ByVar <- dplyr::arrange(PlotObj_Beta, Variable, IAS_ID) %>%
    dplyr::select(Variable, Plot_File, VarDesc) %>%
    tidyr::nest(Plot_File = "Plot_File") %>%
    dplyr::mutate(Plot_File = purrr::map(Plot_File, unlist)) %>%
    arrange_vars()

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(BetaTracePlots_ByVar)), level = 2L,
      future_max_size = 1500, strategy = strategy, cat_timestamp = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  ecokit::cat_time("Save plots, grouped by variables", level = 2L)
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
            ecokit::load_as(.x) %>%
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
          text = ggplot2::element_text(family = "sans"),
          plot.title = ggtext::element_markdown(
            size = 24, hjust = 0.5, margin = ggplot2::margin(t = 15, b = 15)))

      if (!stringr::str_detect(VarDesc, "\n&nbsp;&mdash;&nbsp;")) {
        VarDesc <- paste0(VarDesc, "  ---  ")
      }

      PlotTitleFixed <- ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0(VarDesc, " (fixed y-axis range)")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(family = "sans"),
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

      grDevices::cairo_pdf(
        filename = fs::path(
          Path_Convergence, paste0("Convergence_Beta_", x, "_FreeY.pdf")),
        width = 18, height = 13, onefile = TRUE)
      purrr::walk(BetaPlotList$Plot, grid::grid.draw, recording = FALSE)
      grDevices::dev.off()

      grDevices::cairo_pdf(
        filename = fs::path(
          Path_Convergence, paste0("Convergence_Beta_", x, "_FixedY.pdf")),
        width = 18, height = 13, onefile = TRUE)
      purrr::walk(.x = BetaPlotList$PlotFixedY, .f = grid::grid.draw)
      grDevices::dev.off()

      rm(
        PlotTitle, PlotTitleFixed, BetaPlotsFixedY, BetaPlots, BetaPlotList,
        envir = environment())

      invisible(gc())
      return(NULL)
    },
    future.seed = TRUE, future.packages = pkg_to_export,
    future.globals = c("BetaTracePlots_ByVar", "n_RC", "Path_Convergence")) %>%
    ecokit::quiet_device()

  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Stopping cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Beta - 3. by species ------

  ecokit::cat_time("Trace plots, grouped by species", level = 1L)

  ecokit::cat_time("Preparing data", level = 2L)
  Order <- stringr::str_remove_all(ModVars, "\\(|\\)")
  BetaTracePlots_BySp <- PlotObj_Beta %>%
    dplyr::arrange(Species, factor(Variable, levels = Order)) %>%
    dplyr::select(Species, IAS_ID, Plot_File, Variable, VarDesc) %>%
    tidyr::nest(data = -c("Species", "IAS_ID")) %>%
    dplyr::left_join(SpSummary, by = "Species")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(BetaTracePlots_BySp)), level = 2L,
      future_max_size = 1500, strategy = strategy, cat_timestamp = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  ecokit::cat_time("Save plots, grouped by species", level = 2L)

  BetaTracePlots_BySp0 <- future.apply::future_lapply(
    X = BetaTracePlots_BySp$Species,
    FUN = function(x) {

      PlotTitle <- ggplot2::ggplot() +
        ggplot2::labs(title = paste0("<i>", x, "</i>")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(family = "sans"),
          plot.title = ggtext::element_markdown(
            face = "bold", size = 24, hjust = 0.5))

      SpDT <- dplyr::filter(BetaTracePlots_BySp, Species == x)

      VarOrder <- SpDT$data[[1]]$Variable %>%
        stringr::str_subset("Intercept", negate = TRUE) %>%
        gtools::mixedsort() %>%
        c("Intercept", .)

      temp_file <- fs::file_temp(ext = "pdf")
      grDevices::cairo_pdf(temp_file)

      SpPlots <- SpDT$data[[1]] %>%
        dplyr::arrange(factor(Variable, levels = VarOrder)) %>%
        dplyr::mutate(
          Plot = purrr::map2(
            .x = Plot_File, .y = VarDesc,
            .f = ~ {
              Plot <- ecokit::load_as(.x)$Plot +
                ggplot2::ggtitle(.y) +
                ggplot2::theme(
                  text = ggplot2::element_text(family = "sans"),
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

      grDevices::dev.off()
      try(fs::file_delete(temp_file), silent = TRUE)

      NumPages <- ceiling(length(SpPlots) / (beta_n_RC[1] * beta_n_RC[2]))

      SpPlots2 <- ecokit::split_vector(
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

      grDevices::cairo_pdf(
        filename = fs::path(
          Path_Convergence_BySp,
          paste0(
            "Convergence_Beta_", SpDT$IAS_ID, "_",
            SpDT$Species_File, ".pdf")),
        width = 23, height = 17, onefile = TRUE)
      purrr::walk(SpPlots2, grid::grid.draw)
      grDevices::dev.off()

      rm(PlotTitle, envir = environment())
      return(invisible(NULL))
    },
    future.seed = TRUE, future.packages = pkg_to_export,
    future.globals = c(
      "margin_type", "BetaTracePlots_BySp", "Path_Convergence_BySp",
      "beta_n_RC")) %>%
    ecokit::quiet_device()

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  # Stopping cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
    future::plan("sequential", gc = TRUE)
  }

  rm(BetaTracePlots_BySp0, envir = environment())

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plot model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# plot_Beta_ranges
## |------------------------------------------------------------------------| #

#' @rdname convergence_plots
#' @name convergence_plots
#' @order 4
#' @author Ahmed El-Gabbas
#' @export

convergence_beta_ranges <- function(model_dir) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Variable <- chain <- Var_Sp_File <- Range <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  if (!is.character(model_dir) || length(model_dir) != 1 ||
      !nzchar(model_dir)) {
    ecokit::stop_ctx(
      "The specified model_dir is not a character of length 1 and nchar > 1.",
      model_dir = model_dir, include_backtrace = TRUE)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "The specified model_dir does not exist.", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  # Construct path to beta parameter data
  beta_data_path <- fs::path(
    dirname(model_dir), "Model_Convergence", "Convergence_Beta.qs2")

  # Check if the beta data file exists
  if (!file.exists(beta_data_path)) {
    ecokit::stop_ctx(
      "The beta parameter data file does not exist.",
      beta_data_path = beta_data_path, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Load beta parameter information -----

  Beta_ranges <- ecokit::load_as(beta_data_path, n_threads = 5)

  # number of chains
  n_chains <- length(Beta_ranges$DT[[1]]$Post) # nolint: object_length_linter

  fct_levels <- unique(
    c("Intercept", gtools::mixedsort(unique(Beta_ranges$Variable))))

  Beta_ranges <- Beta_ranges %>%
    dplyr::mutate(
      Variable = forcats::fct(Variable, levels = fct_levels),
      # Calculate the min and max of the beta values for each species
      Range = purrr::map(
        .x = Var_Sp_File,
        .f = ~ {
          Post <- ecokit::load_as(.x)$Post
          purrr::map(
            .x = seq_len(n_chains),
            .f = ~ {
              tibble::tibble(
                chain = .x,
                min = min(Post[[.x]]),
                max = max(Post[[.x]]))
            }) %>%
            dplyr::bind_rows()
        })) %>%
    tidyr::unnest(Range) %>%
    dplyr::select(Variable, chain, min, max)

  plot_subtitle <- stringr::str_glue(
    "Values of beta parameters (<span style='color:red'>minimum</span> and \\
    <span style='color:blue'>maximum</span> per species) across chains")

  # Construct path for saving the plot
  plot_path <- fs::path(
    dirname(model_dir), "Model_Postprocessing", "Beta_min_max_Habitat.jpeg")
  fs::dir_create(dirname(plot_path))

  Beta_plot <- Beta_ranges %>%
    ggplot2::ggplot(ggplot2::aes(x = (chain - 0.125), y = min)) +
    ggplot2::geom_point(pch = 20, colour = "red", size = 0.8, alpha = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(x = (chain + 0.125), y = max), show.legend = FALSE,
      pch = 20, colour = "blue", size = 0.8, alpha = 0.5) +
    ggplot2::scale_x_continuous(breaks = ecokit::integer_breaks(n_chains)) +
    ggplot2::facet_wrap(~Variable, scales = "free_y", ncol = 5, nrow = 4) +
    ggplot2::labs(
      title = "Convergence of beta parameters", subtitle = plot_subtitle,
      x = "Chain",
      y = stringr::str_glue(
        "<span style='color:red'>Minimum</span> and \\
      <span style='color:blue'>maximum</span> beta values")) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      plot.subtitle = ggtext::element_markdown(),
      axis.title = ggtext::element_markdown(face = "bold"))

  ragg::agg_jpeg(
    filename = plot_path, width = 30, height = 20, res = 600,
    quality = 100, units = "cm")
  print(Beta_plot)
  grDevices::dev.off()

  invisible(NULL)
}
