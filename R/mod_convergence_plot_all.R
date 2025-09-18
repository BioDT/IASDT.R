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
#'   The convergence outputs will be saved to the `model_convergence_all`
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

  # Check input arguments ------

  ecokit::cat_time("Check input arguments")
  ecokit::check_args(
    args_to_check = c("model_dir", "margin_type"), args_type = "character")
  ecokit::check_args(
    args_to_check = c("n_omega", "n_rc_alpha"),
    args_type = "numeric", arg_length = c(1L, 2L))
  ecokit::check_args(args_to_check = "spatial_model", args_type = "logical")

  if (!margin_type %in% c("histogram", "density")) {
    ecokit::stop_ctx(
      "`margin_type` must be either 'histogram' or 'density'.",
      margin_type = margin_type, include_backtrace = TRUE)
  }

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  gpp_thin <- m_name_fit <- rl <- m_thin <- m_samples <- omega_gelman <-
    omega_ess <- beta_gelman <- beta_ess <- ess_2 <- path_trace_alpha <- NULL

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

  path_conv_all <- fs::path(model_dir, "model_convergence_all")
  fs::dir_create(path_conv_all)

  model_info <- fs::path(model_dir, "model_info.RData")
  if (!file.exists(model_info)) {
    ecokit::stop_ctx(
      "Model info file does not exist", model_info = model_info,
      include_backtrace = TRUE)
  }
  model_info <- ecokit::load_as(model_info)

  # Extract number of chains
  n_chains <- length(model_info$chain[[1]])

  # # ..................................................................... ###

  prepare_convergence <- function(id) {

    temp_file <- fs::file_temp(ext = "pdf")
    grDevices::cairo_pdf(temp_file)
    on.exit({
      grDevices::dev.off()
      try(fs::file_delete(temp_file), silent = TRUE)
    }, add = TRUE)

    path_coda <- model_info$path_coda[[id]]
    m_name_fit <- model_info$m_name_fit[[id]]
    tree <- model_info$tree[[id]]

    # Prepare traceplot ----

    if (isFALSE(fs::file_exists(path_coda))) {

      path_trace_rho <- path_trace_alpha <- beta_gelman <-
        beta_ess <- omega_gelman <- omega_ess  <- NULL

    } else {

      obj_rho <- paste0(m_name_fit, "_trace_rho")
      path_trace_rho <- fs::path(path_conv_all, paste0(obj_rho, ".RData"))

      obj_alpha <- paste0(m_name_fit, "_trace_alpha")
      path_trace_alpha <- fs::path(path_conv_all, paste0(obj_alpha, ".RData"))

      obj_beta <- paste0(m_name_fit, "_beta")
      path_beta <- fs::path(path_conv_all, paste0(obj_beta, ".RData"))

      obj_omega <- paste0(m_name_fit, "_omega")
      path_omega <- fs::path(path_conv_all, paste0(obj_omega, ".RData"))

      coda_obj <- ecokit::load_as(path_coda)

      # names of coda
      names_coda <- names(coda_obj)
      # Number of chains
      n_chains <- length(coda_obj$Beta)

      # Number of samples
      n_samples <- nrow(coda_obj$Beta[[1]])

      # Rho -----
      if (!file.exists(path_trace_rho) && ("Rho" %in% names_coda)) {
        if (tree == "tree") {
          rho_title <- stringr::str_remove_all(
            string = basename(path_coda), pattern = "_tree|_coda|.RData$|.qs2")

          plot_obj_rho <- IASDT.R::convergence_rho(
            posterior = coda_obj, title = rho_title, margin_type = margin_type,
            n_chains = n_chains, n_samples = n_samples)

          ecokit::save_as(
            object = plot_obj_rho, object_name = obj_rho,
            out_path = path_trace_rho)

          rm(plot_obj_rho, envir = environment())

        } else {
          path_trace_rho <- NA_character_
        }
      } else {
        path_trace_rho <- NA_character_
      }

      if ("Rho" %in% names_coda) coda_obj$Rho <- NULL

      # Alpha -----
      if (!file.exists(path_trace_alpha) && ("Alpha" %in% names_coda) &&
          spatial_model) {
        plot_obj_alpha <- IASDT.R::convergence_alpha(
          posterior = coda_obj,
          title = stringr::str_remove_all(
            basename(path_coda), "_tree|_coda|.RData$|.qs2"),
          n_rc_alpha = n_rc_alpha, margin_type = margin_type,
          n_chains = n_chains, n_samples = n_samples)

        ecokit::save_as(
          object = plot_obj_alpha, object_name = obj_alpha,
          out_path = path_trace_alpha)

        rm(plot_obj_alpha, envir = environment())
      } else {
        path_trace_alpha <- NA_character_
      }
      if ("Alpha" %in% names_coda) coda_obj$Alpha <- NULL

      # Beta -----
      if (file.exists(path_beta)) {

        beta_conv <- ecokit::load_as(path_beta)
        beta_gelman <- beta_conv$beta_gelman
        beta_ess <- beta_conv$beta_ess
        rm(beta_conv, envir = environment())

      } else {

        beta <- magrittr::extract2(coda_obj, "Beta")

        # BETA - effectiveSize
        beta_ess <- coda::effectiveSize(beta)

        # BETA - gelman.diag
        beta_gelman <- coda::gelman.diag(beta, multivariate = FALSE) %>%
          magrittr::extract2("psrf") %>%
          as.data.frame() %>%
          dplyr::pull(1) %>%
          magrittr::set_names(NULL)

        beta_conv <- list(beta_gelman = beta_gelman, beta_ess = beta_ess)
        save(beta_conv, file = path_beta)
        rm(beta_conv, beta, envir = environment())
      }
      if ("Beta" %in% names_coda) coda_obj$Beta <- NULL

      # Omega -----

      if ("Omega" %in% names_coda) {

        if (file.exists(path_omega)) {
          omega_conv <- ecokit::load_as(path_omega)
          omega_ess <- omega_conv$omega_ess
          omega_gelman <- omega_conv$omega_gelman
          rm(omega_conv, envir = environment())

        } else {
          omega <- magrittr::extract2(coda_obj, "Omega") %>%
            magrittr::extract2(1)

          # OMEGA - effectiveSize
          omega_ess <- coda::effectiveSize(omega)

          # OMEGA - gelman.diag
          sel <- sample.int(
            n = dim(omega[[1]])[2], size = min(n_omega, dim(omega[[1]])[2]))
          omega_gelman <- purrr::map(.x = omega, .f = ~ .x[, sel]) %>%
            coda::gelman.diag(multivariate = FALSE) %>%
            magrittr::extract2("psrf") %>%
            as.data.frame() %>%
            dplyr::pull(1) %>%
            magrittr::set_names(NULL)

          omega_conv <- list(omega_gelman = omega_gelman, omega_ess = omega_ess)
          save(omega_conv, file = path_omega)
          rm(sel, omega_conv, omega, envir = environment())
        }

      } else {
        omega_gelman <- omega_ess <- data.frame()
      }
      rm(coda_obj, envir = environment())

    }

    invisible(gc())

    list(
      path_trace_alpha = path_trace_alpha,
      path_trace_rho = path_trace_rho,
      beta_gelman = beta_gelman, beta_ess = beta_ess,
      omega_gelman = omega_gelman, omega_ess = omega_ess)
  }

  # # ..................................................................... ###

  # Processing convergence data -----

  path_data <- fs::path(path_conv_all, "convergence_data.RData")

  if (ecokit::check_data(path_data, warning = FALSE)) {

    ecokit::cat_time("Loading convergence data", level = 1L)
    convergence_data <- ecokit::load_as(path_data)

  } else {

    ecokit::cat_time("Processing convergence data", level = 1L)
    n_cores <- min(n_cores, nrow(model_info))
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

    convergence_data <- model_info %>%
      dplyr::mutate(
        plots = future.apply::future_lapply(
          X = seq_len(nrow(model_info)), FUN = prepare_convergence,
          future.scheduling = Inf, future.seed = TRUE,
          future.packages = pkg_to_export,
          future.globals = c(
            "model_info", "path_conv_all", "n_omega", "prepare_convergence",
            "spatial_model", "margin_type"))) %>%
      ecokit::quiet_device()

    convergence_data <- convergence_data %>%
      dplyr::select(tidyselect::all_of(c("m_name_fit", "plots"))) %>%
      tidyr::unnest_wider("plots") %>%
      # Arrange data alphanumerically by model name
      dplyr::arrange(gtools::mixedorder(m_name_fit))

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

  plot_theme <-  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    axis.title = ggtext::element_markdown(
      size = 20, colour = "darkgrey", face = "bold"),
    axis.text = ggplot2::element_text(size = 16),
    title = ggtext::element_markdown(size = 20, face = "bold", color = "blue"),
    axis.text.y = ggtext::element_markdown(
      hjust = 0, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5)),
    panel.spacing = ggplot2::unit(0.75, "lines"))

  label <- ggplot2::as_labeller(c(
    `2000` = "2000 samples",
    `1000` = "1000 samples",
    `3000` = "3000 samples",
    `4000` = "4000 samples",
    `5000` = "5000 samples",
    tree = "Phylogenetic (taxonomic) tree",
    no_tree = "No phylogenetic (taxonomic) tree"))

  # # ..................................................................... ###

  # Alpha - trace plots ------
  if (!all(is.na(convergence_data$path_trace_alpha)) && spatial_model) {
    convergence_data_alpha <- dplyr::filter(
      convergence_data, !is.na(path_trace_alpha))

    ecokit::cat_time("Alpha - trace plots", level = 2L)
    grDevices::cairo_pdf(
      filename = fs::path(path_conv_all, "trace_plots_alpha.pdf"),
      width = 18, height = 12, onefile = TRUE)
    purrr::walk(
      .x = convergence_data_alpha$path_trace_alpha,
      .f = purrr::safely(~{
        gridExtra::grid.arrange(ecokit::load_as(.x)[[1]])
      }))
    grDevices::dev.off()

    rm(convergence_data_alpha, envir = environment())
  }

  # # ..................................................................... ###

  # Rho - trace plots ------
  ecokit::cat_time("Rho - trace plots", level = 2L)
  path_trace_rho <- NULL

  plot <- convergence_data %>%
    dplyr::filter(stringr::str_detect(m_name_fit, "_tree_")) %>%
    dplyr::mutate(
      Rho = purrr::map_if(
        .x = path_trace_rho,
        .p = ~is.na(.x),
        .f = ~grid::grid.rect(gp = grid::gpar(col = "white")),
        .else = ~ecokit::load_as(.x)))

  if (nrow(convergence_data) > 1) {
    layout_matrix <- matrix(seq_len(2 * 2), nrow = 2, byrow = TRUE)
    plot <- plot$Rho %>%
      gridExtra::marrangeGrob(
        bottom = bquote(paste0("page ", g, " of ", npages)),
        top = grid::textGrob(
          label = "Convergence of the rho parameter",
          gp = grid::gpar(fontface = "bold", fontsize = 20)),
        nrow = 2, ncol = 2, layout_matrix = layout_matrix)
  } else {
    plot <- plot$Rho[[1]]
  }

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::cairo_pdf(
    filename = fs::path(path_conv_all, "trace_plots_rho_phylogenetic.pdf"),
    width = 18, height = 15, onefile = TRUE)
  invisible(plot(plot))
  grDevices::dev.off()

  # # ..................................................................... ###

  if (!all(is.na(convergence_data$omega_gelman))) {

    # Omega - Gelman convergence ------
    ecokit::cat_time("Omega - Gelman convergence", level = 2L)

    plot_path <- fs::path(path_conv_all, paste0("convergence_omega_gelman.pdf"))

    plot_title <- paste0(
      "Gelman convergence diagnostic --- Omega (", n_omega, " samples)")

    plot <- convergence_data %>%
      dplyr::left_join(model_info, by = "m_name_fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rl", "tree", "m_thin", "m_samples", "omega_gelman"))) %>%
      # Add rl as NA if not present
      dplyr::mutate(
        rl = dplyr::coalesce(rl, NA_real_),
        gpp_thin = dplyr::if_else(
          is.na(rl),
          paste0("nonspatial | th", m_thin),
          paste0("gpp", rl, " | th", m_thin)),
        gpp_thin = factor(
          gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
      tidyr::unnest("omega_gelman")

    plot <- plot %>%
      ggplot2::ggplot(
        ggplot2::aes(gpp_thin, omega_gelman), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::scale_y_log10() +
      ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
      ggplot2::labs(title = plot_title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(
        "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
      ggplot2::coord_flip(expand = FALSE) +
      plot_theme

    # Suppress the following message: Scale for y is already present. Adding
    # another scale for y, which will replace the existing scale.
    plot_2 <- plot +
      ggplot2::ylab(
        paste0(
          "Gelman and Rubin's convergence diagnostic ",
          "<sub>(only values between 0.9 and 1.1)</sub>")) +
      ggplot2::coord_cartesian(ylim = c(0.9, 1.1)) +
      plot_theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = plot_path, width = 18, height = 12, onefile = TRUE)
    print(plot)
    print(plot_2)
    grDevices::dev.off()
    rm(plot, plot_2, envir = environment())

    # # .................................................................. ###

    # Omega - Effective sample size -----
    ecokit::cat_time("Omega - Effective sample size", level = 2L)

    plot_path <- fs::path(path_conv_all, paste0("convergence_omega_ess.pdf"))

    plot_title <- paste0(
      "Effective sample size --- Omega (", n_omega, " samples)")

    plot <- convergence_data %>%
      dplyr::left_join(model_info, by = "m_name_fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rl", "tree", "m_thin", "m_samples", "omega_ess"))) %>%
      # Add rl as NA if not present
      dplyr::mutate(
        rl = dplyr::coalesce(rl, NA_real_),
        gpp_thin = dplyr::if_else(
          is.na(rl),
          paste0("nonspatial | th", m_thin),
          paste0("gpp", rl, " | th", m_thin)),
        gpp_thin = factor(
          gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
      tidyr::unnest("omega_ess") %>%
      ggplot2::ggplot(
        ggplot2::aes(gpp_thin, omega_ess), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
      ggplot2::labs(title = plot_title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(paste0("Effective sample size (", n_chains, " chains)")) +
      ggplot2::coord_flip(expand = FALSE) +
      plot_theme

    plot_2 <- convergence_data %>%
      dplyr::left_join(model_info, by = "m_name_fit") %>%
      dplyr::select(
        tidyselect::any_of(
          c("rl", "tree", "m_thin", "m_samples", "omega_ess"))) %>%
      # Add rl as NA if not present
      dplyr::mutate(
        rl = dplyr::coalesce(rl, NA_real_),
        gpp_thin = dplyr::if_else(
          is.na(rl),
          paste0("nonspatial | th", m_thin),
          paste0("gpp", rl, " | th", m_thin)),
        gpp_thin = factor(
          gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
      tidyr::unnest("omega_ess") %>%
      dplyr::mutate(ess_2 = (100 * omega_ess / (m_samples * n_chains)))
    plot_2 <- plot_2 %>%
      ggplot2::ggplot(ggplot2::aes(gpp_thin, ess_2), environment = emptyenv()) +
      ggplot2::geom_violin() +
      ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
      ggplot2::labs(title = plot_title) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Mean effective sample size (%)") +
      ggplot2::coord_flip(expand = FALSE) +
      plot_theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = plot_path, width = 18, height = 12, onefile = TRUE)
    print(plot)
    print(plot_2)
    grDevices::dev.off()
    rm(plot, plot_2, envir = environment())
  }

  # # ..................................................................... ###

  # Beta - Gelman convergence ------
  ecokit::cat_time("Beta - Gelman convergence", level = 2L)

  plot_title <- paste0("Gelman convergence diagnostic --- Beta")

  plot_path <- fs::path(path_conv_all, paste0("convergence_beta_gelman.pdf"))

  plot <- convergence_data %>%
    dplyr::left_join(model_info, by = "m_name_fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rl", "tree", "m_thin", "m_samples", "beta_gelman"))) %>%
    # Add rl as NA if not present
    dplyr::mutate(
      rl = dplyr::coalesce(rl, NA_real_),
      gpp_thin = dplyr::if_else(
        is.na(rl),
        paste0("nonspatial | th", m_thin),
        paste0("gpp", rl, " | th", m_thin)),
      gpp_thin = factor(
        gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
    tidyr::unnest("beta_gelman")
  plot <- plot %>%
    ggplot2::ggplot(
      ggplot2::aes(gpp_thin, beta_gelman), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
    ggplot2::labs(title = plot_title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(
      "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
    ggplot2::coord_flip(expand = FALSE) +
    plot_theme

  plot_2 <- plot +
    ggplot2::ylab(
      paste0(
        "Gelman and Rubin's convergence diagnostic ",
        "<sub>(only values between 0.9 and 1.1)</sub>")) +
    ggplot2::coord_cartesian(ylim = c(0.9, 1.1)) +
    plot_theme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = plot_path, width = 18, height = 12, onefile = TRUE)
  print(plot)
  print(plot_2)
  grDevices::dev.off()
  rm(plot, plot_2, envir = environment())

  # # ..................................................................... ###

  # Beta - Effective sample size -----
  ecokit::cat_time("Beta - Effective sample size", level = 2L)

  plot_path <- fs::path(path_conv_all, paste0("convergence_beta_ess.pdf"))
  plot_title <- "Effective sample size --- Beta"

  plot <- convergence_data %>%
    dplyr::left_join(model_info, by = "m_name_fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rl", "tree", "m_thin", "m_samples", "beta_ess"))) %>%
    # Add rl as NA if not present
    dplyr::mutate(
      rl = dplyr::coalesce(rl, NA_real_),
      gpp_thin = dplyr::if_else(
        is.na(rl),
        paste0("nonspatial | th", m_thin),
        paste0("gpp", rl, " | th", m_thin)),
      gpp_thin = factor(
        gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
    tidyr::unnest("beta_ess")
  plot <- plot %>%
    ggplot2::ggplot(
      ggplot2::aes(gpp_thin, beta_ess), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
    ggplot2::labs(title = plot_title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size (", n_chains, " chains)")) +
    ggplot2::coord_flip(expand = FALSE) +
    plot_theme

  plot_2 <- convergence_data %>%
    dplyr::left_join(model_info, by = "m_name_fit") %>%
    dplyr::select(
      tidyselect::any_of(
        c("rl", "tree", "m_thin", "m_samples", "beta_ess"))) %>%
    # Add rl as NA if not present
    dplyr::mutate(
      rl = dplyr::coalesce(rl, NA_real_),
      gpp_thin = dplyr::if_else(
        is.na(rl),
        paste0("nonspatial | th", m_thin),
        paste0("gpp", rl, " | th", m_thin)),
      gpp_thin = factor(
        gpp_thin, levels = rev(gtools::mixedsort(unique(gpp_thin))))) %>%
    tidyr::unnest("beta_ess") %>%
    dplyr::mutate(ess_2 = (100 * beta_ess / (m_samples * n_chains)))
  plot_2 <- plot_2 %>%
    ggplot2::ggplot(ggplot2::aes(gpp_thin, ess_2), environment = emptyenv()) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(tree ~ m_samples, labeller = label) +
    ggplot2::labs(title = plot_title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    plot_theme

  # # ..................................................................... ###

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = plot_path, width = 18, height = 12, onefile = TRUE)
  print(plot)
  print(plot_2)
  grDevices::dev.off()
  rm(plot, plot_2, envir = environment())

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plotting model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
