## |------------------------------------------------------------------------| #
# convergence_plot_all ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of multiple modelling alternatives
#'
#' This function generates and saves a series of diagnostic plots to assess the
#' convergence of Hmsc models across multiple modelling alternatives. It checks
#' model convergence using trace plots and Gelman-Rubin diagnostics for key
#' model parameters.
#'
#' @param model_dir Character. Path to the root directory of the fitted model.
#'   The convergence outputs will be saved to the `Model_Convergence_All`
#'   subdirectory.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param spatial_model Logical. Whether the model is a spatial model. If `TRUE`
#'   (default), the function will generate additional plots for the model's
#'   `Alpha` parameter.
#' @name convergence_plot_all
#' @inheritParams convergence_plots
#' @author Ahmed El-Gabbas
#' @return The function returns `invisible(NULL)` and does not return any value,
#'   but saves a series of diagnostic plots in the specified path.
#' @export

convergence_plot_all <- function(
    model_dir = NULL, n_omega = 1000L, margin_type = "histogram",
    spatial_model = TRUE, n_cores = NULL, strategy = "multisession",
    n_rc_alpha = c(2L, 3L)) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` must not be NULL",
      model_dir = model_dir, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  GPP_Thin <- M_Name_Fit <- rL <- M_thin <- M_samples <- Omega_Gelman <-
    Omega_ESS <- Beta_Gelman <- Beta_ESS <- ESS2 <- Path_Trace_Alpha <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::cat_time("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "margin_type", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_omega")
  rm(AllArgs, envir = environment())

  if (length(margin_type) != 1) {
    ecokit::stop_ctx(
      "`margin_type` must be a single value.", margin_type = margin_type,
      include_backtrace = TRUE)
  }

  if (!margin_type %in% c("histogram", "density")) {
    ecokit::stop_ctx(
      "`margin_type` must be either 'histogram' or 'density'.",
      margin_type = margin_type, include_backtrace = TRUE)
  }

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "sf", "Hmsc", "coda", "magrittr", "ggplot2", "ecokit",
      "purrr", "stringr", "fs", "IASDT.R", "withr", "ggtext"),
    strategy = strategy)

  # # ..................................................................... ###

  # Prepare/load convergence data ------

  ecokit::cat_time("Prepare or load convergence data")

  path_conv_all <- fs::path(model_dir, "Model_Convergence_All")
  fs::dir_create(path_conv_all)

  Model_Info <- fs::path(model_dir, "Model_Info.RData")
  if (!file.exists(Model_Info)) {
    ecokit::stop_ctx(
      "Model info file does not exist", Model_Info = Model_Info,
      include_backtrace = TRUE)
  }
  Model_Info <- ecokit::load_as(Model_Info)

  # Extract number of chains
  NChains <- length(Model_Info$Chain[[1]])

  # # ..................................................................... ###

  PrepConvergence <- function(ID) {

    temp_file <- fs::file_temp(ext = "pdf")
    grDevices::cairo_pdf(temp_file)
    on.exit({
      grDevices::dev.off()
      try(fs::file_delete(temp_file), silent = TRUE)
    }, add = TRUE)

    path_coda <- Model_Info$Path_Coda[[ID]]
    M_Name_Fit <- Model_Info$M_Name_Fit[[ID]]
    Tree <- Model_Info$Tree[[ID]]

    # Prepare traceplot ----

    if (isFALSE(fs::file_exists(path_coda))) {

      Path_Trace_Rho <- Path_Trace_Alpha <- Beta_Gelman <-
        Beta_ESS <- Omega_Gelman <- Omega_ESS  <- NULL

    } else {

      Obj_Rho <- paste0(M_Name_Fit, "_trace_rho")
      Path_Trace_Rho <- fs::path(path_conv_all, paste0(Obj_Rho, ".RData"))

      Obj_Alpha <- paste0(M_Name_Fit, "_trace_alpha")
      Path_Trace_Alpha <- fs::path(path_conv_all, paste0(Obj_Alpha, ".RData"))

      Obj_Beta <- paste0(M_Name_Fit, "_beta")
      Path_Beta <- fs::path(path_conv_all, paste0(Obj_Beta, ".RData"))

      Obj_Omega <- paste0(M_Name_Fit, "_omega")
      Path_Omega <- fs::path(path_conv_all, paste0(Obj_Omega, ".RData"))

      Coda_Obj <- ecokit::load_as(path_coda)

      # names of coda
      names_coda <- names(Coda_Obj)
      # Number of chains
      n_chains <- length(Coda_Obj$Beta)

      # Number of samples
      n_samples <- nrow(Coda_Obj$Beta[[1]])

      # Rho -----
      if (!file.exists(Path_Trace_Rho) && ("Rho" %in% names_coda)) {
        if (Tree == "Tree") {
          RhoTitle <- stringr::str_remove_all(
            string = basename(path_coda), pattern = "_Tree|_Coda|.RData$|.qs2")

          PlotObj_Rho <- IASDT.R::convergence_rho(
            posterior = Coda_Obj, title = RhoTitle, margin_type = margin_type,
            n_chains = n_chains, n_samples = n_samples)

          ecokit::save_as(
            object = PlotObj_Rho, object_name = Obj_Rho,
            out_path = Path_Trace_Rho)

          rm(PlotObj_Rho, envir = environment())

        } else {
          Path_Trace_Rho <- NA_character_
        }
      } else {
        Path_Trace_Rho <- NA_character_
      }

      if ("Rho" %in% names_coda) Coda_Obj$Rho <- NULL

      # Alpha -----
      if (!file.exists(Path_Trace_Alpha) && ("Alpha" %in% names_coda) &&
          spatial_model) {
        PlotObj_Alpha <- IASDT.R::convergence_alpha(
          posterior = Coda_Obj,
          title = stringr::str_remove_all(
            basename(path_coda), "_Tree|_Coda|.RData$|.qs2"),
          n_rc_alpha = n_rc_alpha, margin_type = margin_type,
          n_chains = n_chains, n_samples = n_samples)

        ecokit::save_as(
          object = PlotObj_Alpha, object_name = Obj_Alpha,
          out_path = Path_Trace_Alpha)

        rm(PlotObj_Alpha, envir = environment())
      } else {
        Path_Trace_Alpha <- NA_character_
      }
      if ("Alpha" %in% names_coda) Coda_Obj$Alpha <- NULL

      # Beta -----
      if (file.exists(Path_Beta)) {

        Beta_conv <- ecokit::load_as(Path_Beta)
        Beta_Gelman <- Beta_conv$Beta_Gelman
        Beta_ESS <- Beta_conv$Beta_ESS
        rm(Beta_conv, envir = environment())

      } else {

        Beta <- magrittr::extract2(Coda_Obj, "Beta")

        # BETA - effectiveSize
        Beta_ESS <- coda::effectiveSize(Beta)

        # BETA - gelman.diag
        Beta_Gelman <- coda::gelman.diag(Beta, multivariate = FALSE) %>%
          magrittr::extract2("psrf") %>%
          as.data.frame() %>%
          dplyr::pull(1) %>%
          magrittr::set_names(NULL)

        Beta_conv <- list(Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS)
        save(Beta_conv, file = Path_Beta)
        rm(Beta_conv, Beta, envir = environment())
      }
      if ("Beta" %in% names_coda) Coda_Obj$Beta <- NULL

      # Omega -----

      if ("Omega" %in% names_coda) {

        if (file.exists(Path_Omega)) {
          Omega_conv <- ecokit::load_as(Path_Omega)
          Omega_ESS <- Omega_conv$Omega_ESS
          Omega_Gelman <- Omega_conv$Omega_Gelman
          rm(Omega_conv, envir = environment())

        } else {
          Omega <- magrittr::extract2(Coda_Obj, "Omega") %>%
            magrittr::extract2(1)

          # OMEGA - effectiveSize
          Omega_ESS <- coda::effectiveSize(Omega)

          # OMEGA - gelman.diag
          sel <- sample.int(
            n = dim(Omega[[1]])[2], size = min(n_omega, dim(Omega[[1]])[2]))
          Omega_Gelman <- purrr::map(.x = Omega, .f = ~ .x[, sel]) %>%
            coda::gelman.diag(multivariate = FALSE) %>%
            magrittr::extract2("psrf") %>%
            as.data.frame() %>%
            dplyr::pull(1) %>%
            magrittr::set_names(NULL)

          Omega_conv <- list(Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS)
          save(Omega_conv, file = Path_Omega)
          rm(sel, Omega_conv, Omega, envir = environment())
        }

      } else {
        Omega_Gelman <- Omega_ESS <- data.frame()
      }
      rm(Coda_Obj, envir = environment())

    }

    invisible(gc())

    return(
      list(
        Path_Trace_Alpha = Path_Trace_Alpha,
        Path_Trace_Rho = Path_Trace_Rho,
        Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS,
        Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS))
  }

  # # ..................................................................... ###

  # Processing convergence data -----

  path_data <- fs::path(path_conv_all, "convergence_data.RData")

  if (ecokit::check_data(path_data, warning = FALSE)) {

    ecokit::cat_time("Loading convergence data", level = 1L)
    convergence_data <- ecokit::load_as(path_data)

  } else {

    ecokit::cat_time("Processing convergence data", level = 1L)
    n_cores <- min(n_cores, nrow(Model_Info))
    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 2L,
        future_max_size = 800L, strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    ecokit::cat_time(
      "Starting processing convergence data", level = 2L, cat_timestamp = FALSE)

    convergence_data <- Model_Info %>%
      dplyr::mutate(
        Plots = future.apply::future_lapply(
          X = seq_len(nrow(Model_Info)), FUN = PrepConvergence,
          future.scheduling = Inf, future.seed = TRUE,
          future.packages = pkg_to_export,
          future.globals = c(
            "Model_Info", "path_conv_all", "n_omega", "PrepConvergence",
            "spatial_model", "margin_type"))) %>%
      ecokit::quiet_device()

    convergence_data <- convergence_data %>%
      dplyr::select(tidyselect::all_of(c("M_Name_Fit", "Plots"))) %>%
      tidyr::unnest_wider("Plots") %>%
      # Arrange data alphanumerically by model name
      dplyr::arrange(gtools::mixedorder(M_Name_Fit))

    ecokit::save_as(
      object = convergence_data, object_name = "convergence_data",
      out_path = path_data)

    if (n_cores > 1) {
      ecokit::set_parallel(
        stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
      future::plan("sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  # Plotting theme -----

  Theme <-  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    axis.title = ggtext::element_markdown(
      size = 20, colour = "darkgrey", face = "bold"),
    axis.text = ggplot2::element_text(size = 16),
    title = ggtext::element_markdown(size = 20, face = "bold", color = "blue"),
    axis.text.y = ggtext::element_markdown(
      hjust = 0, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5)),
    panel.spacing = ggplot2::unit(0.75, "lines"))

  Label <- ggplot2::as_labeller(c(
    `2000` = "2000 samples",
    `1000` = "1000 samples",
    `3000` = "3000 samples",
    `4000` = "4000 samples",
    `5000` = "5000 samples",
    Tree = "Phylogenetic (taxonomic) tree",
    NoTree = "No phylogenetic (taxonomic) tree"))

  # # ..................................................................... ###

  # Alpha - trace plots ------
  if (!all(is.na(convergence_data$Path_Trace_Alpha)) && spatial_model) {
    Convergence_DT_alpha <- dplyr::filter(
      convergence_data, !is.na(Path_Trace_Alpha))

    ecokit::cat_time("Alpha - trace plots", level = 2L)
    grDevices::cairo_pdf(
      filename = fs::path(path_conv_all, "TracePlots_Alpha.pdf"),
      width = 18, height = 12, onefile = TRUE)
    purrr::walk(
      .x = Convergence_DT_alpha$Path_Trace_Alpha,
      .f = purrr::safely(~{
        gridExtra::grid.arrange(ecokit::load_as(.x)[[1]])
      }))
    grDevices::dev.off()

    rm(Convergence_DT_alpha, envir = environment())
  }

  # # ..................................................................... ###

  # Rho - trace plots ------
  ecokit::cat_time("Rho - trace plots", level = 2L)
  Path_Trace_Rho <- NULL

  Plot <- convergence_data %>%
    dplyr::filter(stringr::str_detect(M_Name_Fit, "_Tree_")) %>%
    dplyr::mutate(
      Rho = purrr::map_if(
        .x = Path_Trace_Rho,
        .p = ~is.na(.x),
        .f = ~grid::grid.rect(gp = grid::gpar(col = "white")),
        .else = ~ecokit::load_as(.x)))

  if (nrow(convergence_data) > 1) {
    layout_matrix <- matrix(seq_len(2 * 2), nrow = 2, byrow = TRUE)
    Plot <- Plot$Rho %>%
      gridExtra::marrangeGrob(
        bottom = bquote(paste0("page ", g, " of ", npages)),
        top = grid::textGrob(
          label = "Convergence of the rho parameter",
          gp = grid::gpar(fontface = "bold", fontsize = 20)),
        nrow = 2, ncol = 2, layout_matrix = layout_matrix)
  } else {
    Plot <- Plot$Rho[[1]]
  }

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::cairo_pdf(
    filename = fs::path(path_conv_all, "trace_plots_rho_phylogenetic.pdf"),
    width = 18, height = 15, onefile = TRUE)
  invisible(plot(Plot))
  grDevices::dev.off()

  # # ..................................................................... ###

  if (!all(is.na(convergence_data$Omega_Gelman))) {

    # Omega - Gelman convergence ------
    ecokit::cat_time("Omega - Gelman convergence", level = 2L)

    Plot_Path <- fs::path(path_conv_all, paste0("convergence_omega_gelman.pdf"))

    Plot_Title <- paste0(
      "Gelman convergence diagnostic --- Omega (", n_omega, " samples)")

    Plot <- convergence_data %>%
      dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rL", "Tree", "M_thin", "M_samples", "Omega_Gelman"))) %>%
      # Add rL as NA if not present
      dplyr::mutate(
        rL = dplyr::coalesce(rL, NA_real_),
        GPP_Thin = dplyr::if_else(
          is.na(rL),
          paste0("non-spatial | Th", M_thin),
          paste0("GPP", rL, " | Th", M_thin)),
        GPP_Thin = factor(
          GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
      tidyr::unnest("Omega_Gelman")

    Plot <- Plot %>%
      ggplot2::ggplot(
        ggplot2::aes(GPP_Thin, Omega_Gelman), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::scale_y_log10() +
      ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
      ggplot2::labs(title = Plot_Title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(
        "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
      ggplot2::coord_flip(expand = FALSE) +
      Theme

    # Suppress the following message: Scale for y is already present. Adding
    # another scale for y, which will replace the existing scale.
    Plot2 <- Plot +
      ggplot2::ylab(
        paste0(
          "Gelman and Rubin's convergence diagnostic ",
          "<sub>(only values between 0.9 and 1.1)</sub>")) +
      ggplot2::ylim(c(0.9, 1.1)) +
      Theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
    print(Plot)
    print(Plot2)
    grDevices::dev.off()
    rm(Plot, Plot2, envir = environment())

    # # .................................................................. ###

    # Omega - Effective sample size -----
    ecokit::cat_time("Omega - Effective sample size", level = 2L)

    Plot_Path <- fs::path(path_conv_all, paste0("convergence_omega_ess.pdf"))

    Plot_Title <- paste0(
      "Effective sample size --- Omega (", n_omega, " samples)")

    Plot <- convergence_data %>%
      dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rL", "Tree", "M_thin", "M_samples", "Omega_ESS"))) %>%
      # Add rL as NA if not present
      dplyr::mutate(
        rL = dplyr::coalesce(rL, NA_real_),
        GPP_Thin = dplyr::if_else(
          is.na(rL),
          paste0("non-spatial | Th", M_thin),
          paste0("GPP", rL, " | Th", M_thin)),
        GPP_Thin = factor(
          GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
      tidyr::unnest("Omega_ESS") %>%
      ggplot2::ggplot(
        ggplot2::aes(GPP_Thin, Omega_ESS), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
      ggplot2::labs(title = Plot_Title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(paste0("Effective sample size (", NChains, " chains)")) +
      ggplot2::coord_flip(expand = FALSE) +
      Theme

    Plot2 <- convergence_data %>%
      dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rL", "Tree", "M_thin", "M_samples", "Omega_ESS"))) %>%
      # Add rL as NA if not present
      dplyr::mutate(
        rL = dplyr::coalesce(rL, NA_real_),
        GPP_Thin = dplyr::if_else(
          is.na(rL),
          paste0("non-spatial | Th", M_thin),
          paste0("GPP", rL, " | Th", M_thin)),
        GPP_Thin = factor(
          GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
      tidyr::unnest("Omega_ESS") %>%
      dplyr::mutate(ESS2 = (100 * Omega_ESS / (M_samples * NChains)))
    Plot2 <- Plot2 %>%
      ggplot2::ggplot(ggplot2::aes(GPP_Thin, ESS2), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
      ggplot2::labs(title = Plot_Title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Mean effective sample size (%)") +
      ggplot2::coord_flip(expand = FALSE) +
      Theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
    print(Plot)
    print(Plot2)
    grDevices::dev.off()
    rm(Plot, Plot2, envir = environment())
  }

  # # ..................................................................... ###

  # Beta - Gelman convergence ------
  ecokit::cat_time("Beta - Gelman convergence", level = 2L)

  Plot_Title <- paste0("Gelman convergence diagnostic --- Beta")

  Plot_Path <- fs::path(path_conv_all, paste0("convergence_beta_gelman.pdf"))

  Plot <- convergence_data %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rL", "Tree", "M_thin", "M_samples", "Beta_Gelman"))) %>%
    # Add rL as NA if not present
    dplyr::mutate(
      rL = dplyr::coalesce(rL, NA_real_),
      GPP_Thin = dplyr::if_else(
        is.na(rL),
        paste0("non-spatial | Th", M_thin),
        paste0("GPP", rL, " | Th", M_thin)),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_Gelman")
  Plot <- Plot %>%
    ggplot2::ggplot(
      ggplot2::aes(GPP_Thin, Beta_Gelman), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(
      "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- Plot +
    ggplot2::ylab(
      paste0(
        "Gelman and Rubin's convergence diagnostic ",
        "<sub>(only values between 0.9 and 1.1)</sub>")) +
    ggplot2::ylim(c(0.9, 1.1)) +
    Theme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
  print(Plot)
  print(Plot2)
  grDevices::dev.off()
  rm(Plot, Plot2, envir = environment())

  # # ..................................................................... ###

  # Beta - Effective sample size -----
  ecokit::cat_time("Beta - Effective sample size", level = 2L)

  Plot_Path <- fs::path(path_conv_all, paste0("convergence_beta_ess.pdf"))
  Plot_Title <- "Effective sample size --- Beta"

  Plot <- convergence_data %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rL", "Tree", "M_thin", "M_samples", "Beta_ESS"))) %>%
    # Add rL as NA if not present
    dplyr::mutate(
      rL = dplyr::coalesce(rL, NA_real_),
      GPP_Thin = dplyr::if_else(
        is.na(rL),
        paste0("non-spatial | Th", M_thin),
        paste0("GPP", rL, " | Th", M_thin)),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_ESS")
  Plot <- Plot %>%
    ggplot2::ggplot(
      ggplot2::aes(GPP_Thin, Beta_ESS), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size (", NChains, " chains)")) +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- convergence_data %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rL", "Tree", "M_thin", "M_samples", "Beta_ESS"))) %>%
    # Add rL as NA if not present
    dplyr::mutate(
      rL = dplyr::coalesce(rL, NA_real_),
      GPP_Thin = dplyr::if_else(
        is.na(rL),
        paste0("non-spatial | Th", M_thin),
        paste0("GPP", rL, " | Th", M_thin)),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_ESS") %>%
    dplyr::mutate(ESS2 = (100 * Beta_ESS / (M_samples * NChains)))
  Plot2 <- Plot2 %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, ESS2), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  # # ..................................................................... ###

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
  print(Plot)
  print(Plot2)
  grDevices::dev.off()
  rm(Plot, Plot2, envir = environment())

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plotting model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
