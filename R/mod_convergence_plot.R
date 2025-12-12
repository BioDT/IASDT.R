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
#' @param n_rc List of 3 numeric vectors representing the number of rows and
#'   columns for grid layout of the convergence plots of alpha, omega, and beta
#'   parameters. .
#' @param n_rc_alpha Numeric vector of length 2. Number of rows and columns for
#'   the convergence plots of the `alpha` parameter. Default: `c(2L, 3L)`.
#' @param pages_per_file Integer. Number of plots per page in the Omega
#'   parameter output. Default: 20L.
#' @param chain_colors Character vector. MCMC chain colours (optional). Default:
#'   `NULL`.
#' @param posterior `mcmc.list` or character. Either an MCMC object
#'   (`mcmc.list`) containing posterior samples, or a file path to a saved coda
#'   object.
#' @param add_footer Logical. If `TRUE` (default), adds a footer with page
#'   numbers to each plot.
#' @param add_title Logical. If `TRUE` (default), adds the main title (`title`)
#'   to the plot.
#' @param margin_type Character. The type of marginal plot to add to the main
#'   plot. Valid options are "histogram" (default) or "density".
#' @param spatial_model Logical. Whether the model is a spatial model. If `TRUE`
#'   (default), the function will generate additional plots for the model's
#'   `Alpha` parameter.
#' @param future_max_size	Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details.
#' @param beta_data Data frame. Beta parameter summary data frame.
#' @param n_chains Integer. Number of MCMC chains.
#' @param n_samples Integer. Number of MCMC samples.
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
    path_coda = NULL, env_file = ".env", title = " ", n_omega = 1000L,
    n_cores = 8L, strategy = "multisession", future_max_size = 2000L,
    n_rc = list(alpha = c(2L, 3L), omega = c(2L, 2L), beta = c(3L, 3L)),
    pages_per_file = 20L, chain_colors = NULL, margin_type = "histogram",
    spatial_model = TRUE) {

  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::cat_time("Check input arguments")
  ecokit::check_args(args_to_check = "spatial_model", args_type = "logical")
  ecokit::check_args(
    args_to_check = c("path_coda", "title", "margin_type"),
    args_type = "character")
  ecokit::check_args(
    args_to_check = c("n_omega", "future_max_size", "pages_per_file"),
    args_type = "numeric")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  if (!margin_type %in% c("histogram", "density")) {
    ecokit::stop_ctx(
      "`margin_type` must be either 'histogram' or 'density'.",
      margin_type = margin_type, include_backtrace = TRUE)
  }

  # Validate n_rc
  if (!is.list(n_rc) || length(n_rc) != 3 ||
      !all(sapply(n_rc, is.numeric)) || !all(sapply(n_rc, length) == 2)) {
    ecokit::stop_ctx(
      "`n_rc` must be a list of three numeric vectors of length 2.",
      n_rc = n_rc, include_backtrace = TRUE)
  }
  if (any(sapply(n_rc, function(x) any(x <= 0)))) {
    ecokit::stop_ctx(
      "`n_rc` values must be positive integers.",
      n_rc = n_rc, include_backtrace = TRUE)
  }

  # Validate n_omega
  if (n_omega <= 0) {
    ecokit::stop_ctx(
      "`n_omega` must be a positive integer.",
      n_omega = n_omega, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  species_combs <- `2.5%` <- `97.5%` <- class <- order <- family <- data <-
    ias_id <- variable <- data <- plot_id <- file <- page <- iter <- value <-
    chain <- y <- label <- var_sp_2 <- species_name <- species_file <-
    path_pa <- is_intercept <- var_sp_file <- description <-
    linear_quadratic <- NULL

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "ggplot2", "ggtext", "magrittr", "stringr", "ggExtra",
      "coda", "ecokit", "qs2", "tibble", "tidyr", "purrr", "cowplot",
      "gtools", "withr", "fs", "grid", "ggpubr"),
    strategy = strategy)

  # # ..................................................................... ###

  # # Load species summary
  ecokit::cat_time("Load species summary")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_pa", "DP_R_pa", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  sp_summary <- fs::path(path_pa, "sp_pa_summary_df.RData")
  if (!file.exists(sp_summary)) {
    ecokit::stop_ctx(
      "sp_summary file does not exist", sp_summary = sp_summary,
      include_backtrace = TRUE)
  }

  sp_summary <- readr::read_csv(
    file = sp_summary, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(species = species_name, species_file)

  # # ..................................................................... ###

  # Create path ------

  ecokit::cat_time("Create path")
  dir_convergence <- fs::path(dirname(dirname(path_coda)), "model_convergence")
  dir_beta_data <- fs::path(dir_convergence, "beta_data")
  dir_convergence_by_sp <- fs::path(dir_convergence, "beta_by_species")
  fs::dir_create(c(dir_convergence, dir_convergence_by_sp, dir_beta_data))

  # # ..................................................................... ###

  # Prepare convergence data ------

  ecokit::cat_time("Prepare convergence data")

  if (!file.exists(path_coda)) {
    ecokit::stop_ctx(
      "`path_coda` does not exist", path_coda = path_coda,
      include_backtrace = TRUE)
  }

  ecokit::cat_time("Loading coda object", level = 1L)
  coda_obj <- ecokit::load_as(path_coda)
  names_coda <- names(coda_obj)

  # Number of chains
  n_chains <- length(coda_obj$Beta)

  # Number of samples
  n_samples <- nrow(coda_obj$Beta[[1]])

  #  Plotting colours
  define_chain_colors <- FALSE

  if (is.null(chain_colors)) {
    define_chain_colors <- TRUE
  } else if (length(chain_colors) != n_chains) {
    define_chain_colors <- TRUE
    warning(
      "The length of provided colours != number of chains", call. = FALSE)
  }

  if (define_chain_colors) {
    # minimum value of n colours in RColorBrewer::brewer.pal is 3.
    # black and grey will be used anyway
    if (n_chains >= 4) {
      chain_colors <- c(
        "black", "grey60",
        RColorBrewer::brewer.pal(n = n_chains - 2, name = "Set1"))
    } else if (n_chains == 3) {
      chain_colors <- c("black", "grey60", "red")
    } else if (n_chains == 2) {
      chain_colors <- c("black", "grey60")
    }
  }

  # # ..................................................................... ###

  # Model variables ------

  arrange_vars <- function(df) {
    ecokit::arrange_alphanum(df, variable) %>%
      dplyr::mutate(is_intercept = (variable == "Intercept")) %>%
      dplyr::arrange(dplyr::desc(is_intercept)) %>%
      dplyr::select(-is_intercept)
  }

  vars_desc <- tibble::tribble(
    ~variable, ~description,
    "Intercept", "Intercept",
    "rivers_log", "River length (log<sub>10</sub>)",
    "road_rail_log", "Road + Rail intensity (log<sub>10</sub>)",
    "efforts_log", "Sampling efforts (log<sub>10</sub>)",
    "habitat_log", "% Habitat coverage",
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

  model_terms <- attr(coda_obj$Beta[[1]], "dimnames") %>%
    unlist() %>%
    stringr::str_remove_all("^B\\[|, sp_\\d{4}\\]$") %>%
    unique() %>%
    sort()
  model_vars <- model_terms %>%
    stringr::str_remove_all("stats::poly|, degree = 2, raw = TRUE\\)[0-9]") %>%
    stringr::str_remove_all("\\(|\\)") %>%
    stringr::str_trim()
  model_vars <- dplyr::bind_cols(variable = model_vars, term = model_terms) %>%
    dplyr::mutate(
      linear_quadratic = dplyr::case_when(
        stringr::str_detect(term, "raw = TRUE\\)1") ~ "L",
        stringr::str_detect(term, "raw = TRUE\\)2") ~ "Q",
        .default = "L_only")) %>%
    dplyr::left_join(vars_desc, by = "variable") %>%
    dplyr::mutate(
      variable = purrr::map2_chr(
        .x = variable, .y = linear_quadratic,
        .f = ~ {
          dplyr::case_when(
            .y == "L" ~ paste0(.x, "_L"),
            .y == "Q" ~ paste0(.x, "_Q"),
            .default = .x)
        }),
      description = purrr::map2_chr(
        .x = description, .y = linear_quadratic,
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
      variable = stringr::str_replace(
        variable, "\\(Intercept\\)", "Intercept"),
      linear_quadratic = NULL) %>%
    arrange_vars()

  # # ..................................................................... ###

  # Rho ------

  if ("Rho" %in% names_coda) {

    ecokit::cat_time("Rho")
    plot_obj_rho <- IASDT.R::convergence_rho(
      posterior = coda_obj, title = title, chain_colors = chain_colors,
      n_chains = n_chains, n_samples = n_samples)

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(dir_convergence, "convergence_rho.pdf"),
      width = 18, height = 12, onefile = TRUE)
    plot(plot_obj_rho)
    grDevices::dev.off()

    coda_obj$Rho <- NULL
    rm(plot_obj_rho, envir = environment())
    invisible(gc())
  }

  # # ..................................................................... ###

  # Alpha  ------

  if (("Alpha" %in% names_coda) && spatial_model) {

    ecokit::cat_time("Alpha")

    # Ensure that all latent factors of the model are plotted
    n_lf <- ncol(coda_obj$Alpha[[1]][[1]])
    n_rc_alpha <- n_rc$alpha
    if ((n_rc_alpha[1] * n_rc_alpha[2]) < n_lf) {
      if (n_rc_alpha[1] == 1) {
        n_rc_alpha[1] <- 2
      } else if (n_rc_alpha[1] == 2) {
        n_rc_alpha[1] <- 3
      }
    }

    plot_obj_alpha <- IASDT.R::convergence_alpha(
      posterior = coda_obj, title = title, n_rc_alpha = n_rc_alpha,
      add_footer = FALSE, add_title = FALSE, chain_colors = chain_colors,
      n_chains = n_chains, n_samples = n_samples)

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(dir_convergence, "convergence_alpha.pdf"),
      width = 18, height = 14, onefile = TRUE)
    print(plot_obj_alpha)
    grDevices::dev.off()

    rm(plot_obj_alpha, n_rc_alpha, envir = environment())
  }

  if ("Alpha" %in% names_coda) coda_obj$Alpha <- NULL

  obj_omega <- coda_obj$Omega[[1]]
  obj_beta <- coda_obj$Beta
  rm(coda_obj, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Omega  ------

  if ("Omega" %in% names_coda) {

    ecokit::cat_time("Omega")

    ecokit::cat_time("Coda to tibble", level = 1L)
    omega_df <- IASDT.R::coda_to_tibble(
      coda_object = obj_omega, posterior_type = "omega", n_omega = n_omega,
      env_file = env_file)
    invisible(gc())
    selected_combinations <- unique(omega_df$species_combs)

    ecokit::cat_time("Prepare confidence interval data", level = 1L)
    ci <- purrr::map(.x = obj_omega, .f = ~ .x[, selected_combinations]) %>%
      coda::mcmc.list() %>%
      summary(quantiles = c(0.025, 0.975)) %>%
      magrittr::extract2("quantiles") %>%
      as.data.frame() %>%
      tibble::as_tibble(rownames = "species_combs") %>%
      stats::setNames(c("species_combs", "ci_25", "ci_975"))

    omega_df <- dplyr::left_join(omega_df, ci, by = "species_combs")
    omega_names <- attr(obj_omega[[1]], "dimnames")[[2]]

    # Prepare omega plots
    ecokit::cat_time("Prepare omega plots", level = 1L)

    plot_obj_omega <- purrr::map_dfr(
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

        comb_data <- dplyr::filter(
          omega_df, species_combs == selected_combinations[x])
        curr_post <- purrr::map(
          .x = obj_omega,
          .f = ~ .x[, which(omega_names == comb_data$species_combs)]) %>%
          coda::as.mcmc.list()

        ## gelman convergence diagnostic
        label_gelman <- try(
          coda::gelman.diag(curr_post, multivariate = FALSE), silent = TRUE)

        if (inherits(label_gelman, "try-error")) {
          label_gelman <- coda::gelman.diag(
            curr_post, multivariate = FALSE, autoburnin = FALSE)
        }

        label_gelman <- round(label_gelman$psrf[[1]], 3) %>%
          paste0("<b><i>Gelman convergence diagnostic:</i></b> ", .) %>%
          data.frame(x = -Inf, y = Inf, label = .)

        ## Effective sample size
        label_ess <- round((coda::effectiveSize(curr_post) / n_chains), 1) %>%
          paste0(
            "<b><i>Mean effective sample size:</i></b> ", ., " / ", n_samples)
        curr_ci <- dplyr::select(comb_data, c("ci_25", "ci_975")) %>%
          unlist() %>%
          round(2)
        label_ci <- paste0(
          "<b><i>95% credible interval:</i></b> ",
          paste(curr_ci, collapse = " to "))
        label_ess_ci <- data.frame(
          x = -Inf, y = -Inf, label = paste0(label_ess, "<br>", label_ci))

        label_panel <- sort(c(comb_data$naps_1, comb_data$naps_2)) %>%
          paste0("<i>", ., "</i>") %>%
          paste(collapse = " & <br>") %>%
          paste0(" ") %>%
          data.frame(x = Inf, y = Inf, label = .)

        plot <- ggplot2::ggplot(
          data = comb_data$data[[1]], environment = emptyenv(),
          mapping = ggplot2::aes(
            x = iter, y = value, color = factor(chain))) +
          ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
          ggplot2::geom_smooth(
            method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
          ggplot2::geom_point(alpha = 0) +
          ggplot2::geom_hline(
            yintercept = curr_ci, linetype = "dashed", color = "black",
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
            data = label_panel, inherit.aes = FALSE, colour = "blue",
            hjust = 1, vjust = 1, lineheight = 0, fill = NA,
            label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = label_gelman, inherit.aes = FALSE, size = 6, hjust = 0,
            vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
          ggtext::geom_richtext(
            mapping = ggplot2::aes(x = x, y = y, label = label),
            data = label_ess_ci, inherit.aes = FALSE, size = 6, hjust = 0,
            vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            axis.text = ggplot2::element_text(size = 12),
            legend.position = "none")

        if (margin_type == "histogram") {
          plot <- ggExtra::ggMarginal(
            p = plot, type = margin_type, margins = "y", size = 6,
            color = "steelblue4", fill = "steelblue4", bins = 100)
        } else {
          plot <- ggExtra::ggMarginal(
            p = plot, type = margin_type, margins = "y", size = 6,
            color = "steelblue4")
        }
        # Making marginal background matching the plot background
        # https://stackoverflow.com/a/78196022/3652584
        plot$layout$t[1] <- 1
        plot$layout$r[1] <- max(plot$layout$r)

        tibble::tibble(
          species_combs = comb_data$species_combs, plot = list(plot))
      }
    )

    rm(
      omega_df, obj_omega, selected_combinations, ci, omega_names,
      envir = environment())
    invisible(gc())

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    ecokit::cat_time("Arrange plots", level = 1L)
    n_rc_omega <- n_rc$omega

    omega_plot_list <- tibble::tibble(plot_id = seq_along(plot_obj_omega)) %>%
      dplyr::mutate(
        file = ceiling(
          plot_id / (pages_per_file * n_rc_omega[2] * n_rc_omega[1])),
        page = ceiling(plot_id / (n_rc_omega[2] * n_rc_omega[1]))) %>%
      tidyr::nest(.by = c("file", "page"), .key = "plot_id") %>%
      dplyr::mutate(
        plot_id = purrr::map(plot_id, ~ unlist(as.vector(.x))),
        plot_id = purrr::pmap(
          .l = list(file, page, plot_id),
          .f = function(file, page, plot_id) {

            plot_title <- ggplot2::ggplot(environment = emptyenv()) +
              ggplot2::labs(
                title = paste0(
                  "Convergence of the omega parameter ---  a sample of ",
                  n_omega, " species pairs"),
                subtitle = paste0(
                  "   file ", file, " | page ",
                  (page - ((file - 1) * pages_per_file)))) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                text = ggplot2::element_text(family = "sans"),
                plot.title = ggtext::element_markdown(
                  face = "bold", size = 20, hjust = 0.5),
                plot.subtitle = ggplot2::element_text(
                  size = 12, colour = "grey",
                  margin = ggplot2::margin(-5, 0, 0, 0)))

            ecokit::quietly({
              plot <- cowplot::plot_grid(
                plotlist = plot_obj_omega$plot[plot_id],
                ncol = n_rc_omega[2], nrow = n_rc_omega[1], align = "hv")
            },
            "Removed [0-9]+ rows containing non-finite outside the scale range")

            plot %>%
              cowplot::plot_grid(
                plot_title, ., ncol = 1, rel_heights = c(0.05, 1))
          }
        )
      )

    # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

    ecokit::cat_time("Save plots", level = 1L)
    ecokit::cat_time(
      paste0(
        "Saving omega convergence plots for ", n_omega,
        " species associations to ", length(unique(omega_plot_list$file)),
        " pdf files"),
      level = 2L, cat_timestamp = FALSE)

    purrr::walk(
      .x = seq_along(unique(omega_plot_list$file)),
      .f = ~ {
        invisible({
          curr_plot_order <- dplyr::filter(omega_plot_list, file == .x)
          file_omega <- fs::path(
            dir_convergence, paste0("convergence_omega_", .x, ".pdf"))
          grDevices::cairo_pdf(
            filename = file_omega, width = 18,
            height = 13.5, onefile = TRUE)
          purrr::map(
            curr_plot_order$plot_id, grid::grid.draw, recording = FALSE)
          grDevices::dev.off()
        })
      })

    rm(omega_plot_list, plot_obj_omega, n_rc_omega, envir = environment())
    invisible(gc())
  }

  # ..................................................................... ###

  # Beta - 1. Prepare data ------

  ecokit::cat_time("Beta")

  ecokit::cat_time("Prepare trace plots", level = 1L)
  beta_names <- attr(obj_beta[[1]], "dimnames")[[2]]

  ecokit::cat_time("Prepare 95% credible interval data", level = 2L)
  ci <- summary(obj_beta, quantiles = c(0.025, 0.975))$quantiles %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "var_sp") %>%
    dplyr::rename(ci_025 = `2.5%`, ci_975 = `97.5%`)

  ecokit::cat_time("Coda to tibble", level = 2L)
  beta_data <- IASDT.R::coda_to_tibble(
    coda_object = obj_beta, posterior_type = "beta", env_file = env_file) %>%
    dplyr::left_join(ci, by = "var_sp")

  # variable ranges
  ecokit::cat_time("variable ranges", level = 2L)
  vars_ranges <- dplyr::arrange(beta_data, variable, ias_id) %>%
    dplyr::select(variable, data) %>%
    dplyr::mutate(
      Range = purrr::map(.x = data, .f = ~ range(dplyr::pull(.x, value)))) %>%
    dplyr::select(-data) %>%
    tidyr::nest(data = "Range") %>%
    dplyr::mutate(
      Range = purrr::map(
        .x = data,
        .f = ~ {
          dplyr::pull(.x, Range) %>%
            as.vector() %>%
            range() %>%
            purrr::set_names(c("var_min", "var_max"))
        })) %>%
    dplyr::select(-data) %>%
    ecokit::arrange_alphanum(variable) %>%
    tidyr::unnest_wider("Range") %>%
    arrange_vars()

  # Species taxonomy
  ecokit::cat_time("Species taxonomy", level = 2L)
  species_taxonomy <- IASDT.R::get_species_name(env_file = env_file) %>%
    dplyr::select(ias_id, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  ecokit::cat_time("Preparing plotting data", level = 2L)

  # untar the tar file containing data (if exists) to the beta_data directory
  beta_data_tar <- fs::path(dir_beta_data, "beta_data.tar")
  if (file.exists(beta_data_tar)) {
    utils::untar(beta_data_tar, exdir = dir_beta_data)
  }

  beta_data <- beta_data %>%
    dplyr::left_join(vars_ranges, by = "variable") %>%
    dplyr::left_join(species_taxonomy, by = "ias_id") %>%
    dplyr::left_join(model_vars, by = "variable") %>%
    dplyr::mutate(
      var_sp_2 = paste0(variable, "_", ias_id),
      var_sp_file = fs::path(dir_beta_data, paste0(var_sp_2, ".qs2")))

  purrr::walk(
    .x = seq_len(nrow(beta_data)),
    .f = ~ {

      option_data <- beta_data[.x, ]

      if (ecokit::check_data(option_data$var_sp_file, warning = FALSE)) {
        return(NULL)
      }

      beta_id <- which(beta_names == option_data$var_sp)
      post_mcmc <- coda::as.mcmc.list(obj_beta[, beta_id])

      gelman <- try(
        coda::gelman.diag(post_mcmc, multivariate = FALSE), silent = TRUE)
      if (inherits(gelman, "try-error")) {
        gelman <- coda::gelman.diag(
          post_mcmc, multivariate = FALSE, autoburnin = FALSE)
      }

      ess <- as.vector(coda::effectiveSize(post_mcmc))

      output <- option_data %>%
        dplyr::mutate(
          gelman = list(gelman), ess = ess,
          beta_id = beta_id, post_mcmc = list(post_mcmc)) %>%
        dplyr::rename(plotting_data = data)
      ecokit::save_as(object = output, out_path = option_data$var_sp_file)
      return(NULL)
    })

  columns_to_remove <- c("data", "ci_025", "ci_975", "var_min", "var_max")
  beta_data <- dplyr::select(beta_data, -tidyselect::all_of(columns_to_remove))

  rm(
    ci, vars_ranges, species_taxonomy, columns_to_remove, obj_beta, beta_names,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # minimum and maximum value of each beta parameter -----
  ecokit::cat_time(
    "plot minimum and maximum value of each beta parameter", level = 1L)

  IASDT.R::convergence_beta_ranges(
    model_dir = fs::path(dirname(dirname(path_coda)), "model_fitted"),
    beta_data = beta_data, n_chains = n_chains)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

  plot_beta <- function(
    data_option, fixed_y = TRUE, add_taxonomy = TRUE,
    margin_type = "histogram", n_chains = NULL, n_samples = NULL,
    chain_colors = NULL, plot_title = NULL) {

    x <- NULL
    var_sp_file <- data_option$var_sp_file

    # check if input data exists
    if (isFALSE(ecokit::check_data(var_sp_file, warning = FALSE))) {
      ecokit::stop_ctx(
        "file does not exist.", var_sp_file = var_sp_file,
        include_backtrace = TRUE)
    }

    temp_file <- fs::file_temp(ext = "pdf")
    grDevices::cairo_pdf(temp_file)
    on.exit({
      grDevices::dev.off()
      try(fs::file_delete(temp_file), silent = TRUE)
    },
    add = TRUE)

    data_all <- ecokit::load_as(var_sp_file)
    if (is.null(data_all) || !inherits(data_all, "data.frame")) {
      ecokit::stop_ctx(
        "Loaded data is invalid", var_sp_file = var_sp_file,
        data_all = data_all, class_data_all = class(data_all),
        include_backtrace = TRUE)
    }

    ## gelman convergence diagnostic
    label_gelman <- round(data_all$gelman[[1]]$psrf, 3) %>%
      paste(collapse = " / ") %>%
      paste0("<b><i>gelman convergence diagnostic:</i></b> ", .) %>%
      data.frame(x = Inf, y = -Inf, label = .)

    ## Effective sample size / ci
    label_ess <- paste0(
      "<b><i>Mean effective sample size:</i></b> ",
      round(data_all$ess / n_chains), " / ", n_samples)

    curr_ci <- c(data_all$ci_025, data_all$ci_975)
    label_ci <- paste(round(curr_ci, 4), collapse = " to ") %>%
      paste0("<b><i>95% credible interval:</i></b> ", .)
    label_ess_ci <- data.frame(
      x = -Inf, y = -Inf, label = paste0(label_ess, "<br>", label_ci))

    plot <- ggplot2::ggplot(
      data = data_all$plotting_data[[1]], environment = emptyenv(),
      mapping = ggplot2::aes(x = iter, y = value, color = factor(chain))) +
      ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
      ggplot2::geom_smooth(
        method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
      ggplot2::geom_point(alpha = 0) +
      ggplot2::geom_hline(
        yintercept = curr_ci, linetype = "dashed", color = "black",
        linewidth = 1) +
      # Ensure that y-axis always show 0
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "transparent", linewidth = 0.6) +
      ggplot2::scale_color_manual(values = chain_colors) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggtext::geom_richtext(
        mapping = ggplot2::aes(x = x, y = y, label = label),
        data = label_gelman, inherit.aes = FALSE, size = 3.5, hjust = 1,
        vjust = -1.5, lineheight = 0, fill = NA, label.color = NA) +
      ggtext::geom_richtext(
        mapping = ggplot2::aes(x = x, y = y, label = label),
        data = label_ess_ci, inherit.aes = FALSE, size = 3.5, hjust = 0,
        vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "sans"),
        legend.position = "none",
        axis.text = ggplot2::element_text(size = 12))

    if (add_taxonomy) {
      label_panel <- data.frame(
        x = Inf, y = Inf,
        label = paste0("<br><b><i>", data_all$species, "</i></b>"))
      panel_title <- c(data_all$class, data_all$order, data_all$family) %>%
        paste(collapse = " | ") %>%
        paste0("<b>", ., "</b>") %>%
        paste0("<br>", data_all$ias_id) %>%
        data.frame(x = -Inf, y = Inf, label = .)

      plot <- plot +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = label_panel, inherit.aes = FALSE, colour = "blue", hjust = 1,
          vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = panel_title, inherit.aes = FALSE, hjust = 0, vjust = 1,
          lineheight = 0, fill = NA, label.color = NA)
    }

    if (fixed_y) {
      vars_limits <- c(data_all$var_min, data_all$var_max)
      ecokit::quietly({
        plot <- plot + ggplot2::scale_y_continuous(limits = vars_limits)
      },
      "Scale for \\w+ is already present.")
    }

    if (!is.null(plot_title)) {
      plot <- plot +
        ggplot2::ggtitle(plot_title) +
        ggplot2::theme(
          plot.title.position = "panel",
          plot.margin = ggplot2::margin(1, 0, 1, 0, unit = "pt"),
          plot.title = ggtext::element_markdown(
            hjust = 0.5, size = 11, margin = ggplot2::margin(-0.5, 0, -2.5, 0)))
    }

    if (margin_type == "histogram") {
      plot <- ggExtra::ggMarginal(
        p = plot, type = margin_type, margins = "y", size = 6,
        color = "steelblue4", fill = "steelblue4", bins = 100)
    } else {
      plot <- ggExtra::ggMarginal(
        p = plot, type = margin_type, margins = "y", size = 6,
        color = "steelblue4")
    }

    return(plot)
  }

  # *******************************************************************

  # Prepare beta plots
  ecokit::cat_time("Prepare beta plots", level = 2L)

  # Beta - 2. by variable ------

  ecokit::cat_time("Trace plots, grouped by variables", level = 1L)

  plots_by_variable <- beta_data %>%
    dplyr::mutate(variable = factor(variable, levels = model_vars$variable)) %>%
    dplyr::arrange(ias_id) %>%
    dplyr::select(variable, var_sp_file, description) %>%
    tidyr::nest(var_data = -c("variable", "description")) %>%
    dplyr::mutate(variable = as.character(variable))

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 2L, strategy = strategy,
      future_max_size = future_max_size)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  plots_by_variable2 <- future.apply::future_lapply(
    X =  seq_len(nrow(plots_by_variable)),
    FUN = function(x) {

      plot_id <- NULL
      n_rc_beta <- n_rc$beta

      plot_title_orig <- plots_by_variable$description[x]
      data_list <- unname(unlist(plots_by_variable$var_data[[x]]))
      title_theme <- ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(family = "sans"),
          plot.title = ggtext::element_markdown(
            size = 24, hjust = 0.5, margin = ggplot2::margin(t = 15, b = 15)))

      plots <- purrr::map(
        .x = data_list,
        .f = function(y) {
          plot_beta(
            data_option = ecokit::load_as(y), fixed_y = FALSE,
            n_chains = n_chains, n_samples = n_samples,
            chain_colors = chain_colors) %>%
            ggpubr::as_ggplot()
        })

      plot_title <- ggplot2::ggplot(environment = emptyenv()) +
        ggplot2::labs(title = plot_title_orig) +
        title_theme

      plot_list <- tibble::tibble(plot_id = seq_along(data_list)) %>%
        dplyr::mutate(
          page = ceiling(plot_id / (n_rc_beta[2] * n_rc_beta[1]))) %>%
        tidyr::nest(.by = "page", .key = "plot_id") %>%
        dplyr::mutate(
          plot = purrr::map(
            .x = plot_id,
            .f = function(z) {
              cowplot::plot_grid(
                plotlist = plots[unlist(as.vector(z))],
                ncol = n_rc_beta[2], nrow = n_rc_beta[1], align = "hv") %>%
                cowplot::plot_grid(
                  plot_title, ., ncol = 1, rel_heights = c(0.035, 1))
            }))

      file_free_y <- fs::path(
        dir_convergence,
        paste0(
          "convergence_beta_", plots_by_variable$variable[x], "_free_y.pdf"))
      grDevices::cairo_pdf(
        filename = file_free_y, width = 18, height = 13, onefile = TRUE)
      purrr::walk(plot_list$plot, grid::grid.draw, recording = FALSE)
      grDevices::dev.off()

      rm(plot_title, plots, plot_list, envir = environment())
      invisible(gc())

      # **********

      if (!stringr::str_detect(plot_title_orig, "\n&nbsp;&mdash;&nbsp;")) {
        plot_title_orig <- paste0(plot_title_orig, "  ---  ")
      }

      plot_title <- ggplot2::ggplot(environment = emptyenv()) +
        ggplot2::labs(
          title = paste0(plot_title_orig, " (fixed y-axis range)")) +
        title_theme

      plots <- purrr::map(
        .x = data_list,
        .f = function(y) {
          plot_beta(
            data_option = ecokit::load_as(y), fixed_y = TRUE,
            n_chains = n_chains, n_samples = n_samples,
            chain_colors = chain_colors) %>%
            ggpubr::as_ggplot()
        })

      plot_list <- tibble::tibble(plot_id = seq_along(data_list)) %>%
        dplyr::mutate(
          page = ceiling(plot_id / (n_rc_beta[2] * n_rc_beta[1]))) %>%
        tidyr::nest(.by = "page", .key = "plot_id") %>%
        dplyr::mutate(
          plot = purrr::map(
            .x = plot_id,
            .f = function(z) {
              cowplot::plot_grid(
                plotlist = plots[unlist(as.vector(z))],
                ncol = n_rc_beta[2], nrow = n_rc_beta[1], align = "hv") %>%
                cowplot::plot_grid(
                  plot_title, ., ncol = 1, rel_heights = c(0.035, 1))
            }))

      file_fixed_y <- fs::path(
        dir_convergence,
        paste0(
          "convergence_beta_", plots_by_variable$variable[x], "_fixed_y.pdf"))

      grDevices::cairo_pdf(
        filename = file_fixed_y, width = 18, height = 13, onefile = TRUE)
      purrr::walk(plot_list$plot, grid::grid.draw, recording = FALSE)
      grDevices::dev.off()

      invisible(NULL)
    },
    future.seed = TRUE, future.packages = pkg_to_export,
    future.globals = c(
      "beta_trace_plots_by_var", "n_rc", "dir_convergence", "n_chains",
      "chain_colors", "n_samples", "plot_beta", "plots_by_variable"))

  rm(plots_by_variable2, envir = environment())
  invisible(gc())

  # Stopping cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Beta - 3. by species ------

  ecokit::cat_time("Trace plots, grouped by species", level = 1L)
  cols_to_keep <- c(
    "species", "variable", "class", "order", "family", "ias_id",
    "var_sp_file", "description")

  plots_by_species <- beta_data %>%
    dplyr::arrange(variable, ias_id) %>%
    dplyr::select(tidyselect::all_of(cols_to_keep)) %>%
    tidyr::nest(var_data = -c("class", "order", "family", "species", "ias_id"))

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 2L, strategy = strategy,
      future_max_size = future_max_size)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  plots_by_species2 <- future.apply::future_lapply(
    X =  seq_len(nrow(plots_by_species)),
    FUN = function(x) {

      n_rc_beta <- n_rc$beta

      n_per_page <- n_rc_beta[1] * n_rc_beta[2]
      title_text <- paste0("<i>", plots_by_species$species[[x]], "</i>")
      subtitle_text <- paste0(
        "**Class:** ", plots_by_species$class[[x]], " / **Order:** ",
        plots_by_species$order[[x]], " / **Family:** ",
        plots_by_species$family[[x]])

      variable_order <- model_vars$variable %>%
        stringr::str_subset("Intercept", negate = TRUE) %>%
        gtools::mixedsort() %>%
        c("Intercept", .)

      species_data <- plots_by_species$var_data[[x]] %>%
        dplyr::mutate(variable = factor(variable, levels = variable_order)) %>%
        dplyr::arrange(variable)
      n_vars <- nrow(species_data)

      plots <- purrr::map(
        .x = seq_len(nrow(species_data)),
        .f = function(y) {
          plot0 <- plot_beta(
            data_option = ecokit::load_as(species_data$var_sp_file[y]),
            fixed_y = FALSE, n_chains = n_chains, n_samples = n_samples,
            chain_colors = chain_colors, add_taxonomy = FALSE,
            plot_title = species_data$description[y])

          # Making marginal background matching the plot background
          # https://stackoverflow.com/a/78196022/3652584
          plot0$layout$t[1] <- 1
          plot0$layout$r[1] <- max(plot0$layout$r)
          plot0
        })

      plot_title <- ggplot2::ggplot() +
        ggtext::geom_richtext(
          ggplot2::aes(x = 0.5, y = 1, label = title_text),
          size = 8, hjust = 0.5, vjust = 1, fill = NA, label.color = NA,
          fontface = "bold") +
        ggtext::geom_richtext(
          ggplot2::aes(x = 0, y = 0.95, label = subtitle_text),
          size = 5, hjust = 0, vjust = 1, fill = NA, label.color = NA,
          colour = "darkblue") +
        ggtext::geom_richtext(
          ggplot2::aes(x = 0.9, y = 0.95, label = plots_by_species$ias_id[[x]]),
          size = 5, hjust = 0, vjust = 1, fill = NA, label.color = NA,
          colour = "red", fontface = "bold") +
        ggplot2::theme_void() +
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

      plot_list <- split(
        seq_len(n_vars), ceiling(seq_along(plots) / n_per_page)) %>%
        purrr::map(
          .f = ~ {
            cowplot::plot_grid(
              plotlist = plots[.x], ncol = n_rc_beta[2],
              nrow = n_rc_beta[1], align = "hv") %>%
              cowplot::plot_grid(
                plot_title, ., ncol = 1, rel_heights = c(0.05, 1))
          })
      rm(plot_title, plots, envir = environment())

      file_free_y <- fs::path(
        dir_convergence_by_sp,
        paste0("convergence_beta_", plots_by_species$species[x], ".pdf"))
      grDevices::cairo_pdf(
        filename = file_free_y, width = 18, height = 13, onefile = TRUE)
      purrr::walk(plot_list, grid::grid.draw, recording = FALSE)
      grDevices::dev.off()

      invisible(gc())
      invisible(NULL)
    },
    future.seed = TRUE, future.packages = pkg_to_export,
    future.globals = c(
      "beta_trace_plots_by_var", "n_rc", "dir_convergence", "n_chains",
      "chain_colors", "n_samples", "plot_beta", "plots_by_species"))

  rm(plots_by_species2, envir = environment())
  invisible(gc())

  # Stopping cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # tar all qs2 files in the `dir_beta_data` directory
  ecokit::cat_time("Tar all .qs2 files", level = 1L)

  qs2_files <- list.files(path = dir_beta_data, pattern = "\\.qs2$")
  qs2_files_2 <- paste(qs2_files, collapse = " ")

  # Command to create the tar file
  tar_command <- stringr::str_glue(
    "cd {fs::path_abs(dir_beta_data)}; tar -cf {basename(beta_data_tar)} \\
    -b 2048 {qs2_files_2}")

  # Create tar file
  system(tar_command)

  # Change the permission of the tar file
  Sys.chmod(beta_data_tar, "755", use_umask = FALSE)

  # Delete all .qs2 files
  try(fs::file_delete(fs::path(dir_beta_data, qs2_files)), silent = TRUE)
  rm(qs2_files_2, envir = environment())

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plot model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# plot_beta_ranges
## |------------------------------------------------------------------------| #

#' @rdname convergence_plots
#' @name convergence_plots
#' @order 4
#' @author Ahmed El-Gabbas
#' @export

convergence_beta_ranges <- function(
    model_dir = NULL, beta_data = NULL, n_chains = NULL) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  variable <- var_sp_file <- max_val <- min_val <- chain <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  if (!is.character(model_dir) || length(model_dir) != 1 ||
      !nzchar(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `model_dir`: must be a character vector of length 1 ",
        "with at least one character."),
      model_dir = model_dir, include_backtrace = TRUE)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "The specified model_dir does not exist.", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!is.data.frame(beta_data)) {
    ecokit::stop_ctx(
      "`beta_data` is not a data frame object",
      class_beta_object = class(beta_data), include_backtrace = TRUE)
  }
  if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains < 1 ||
      n_chains %% 1 != 0) {
    ecokit::stop_ctx(
      "n_chains must be a single positive integer.",
      n_chains = n_chains, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Load beta parameter information -----

  # calculate minimum and maximum value for combination of species and cv
  keep_vars <- c("variable", "chain", "species", "min_val", "max_val")
  beta_ranges <- beta_data %>%
    dplyr::mutate(
      min_max = purrr::map(
        .x = var_sp_file,
        .f = ~ {
          ecokit::load_as(.x) %>%
            dplyr::pull("plotting_data") %>%
            magrittr::extract2(1) %>%
            dplyr::reframe(
              min_val = min(value), max_val = max(value), .by = chain) %>%
            dplyr::mutate(chain = as.integer(chain))
        })) %>%
    tidyr::unnest("min_max") %>%
    dplyr::select(tidyselect::all_of(keep_vars)) %>%
    dplyr::summarise(
      min_val = min(min_val), max_val = max(max_val),
      .by = c("variable", "chain", "species"))

  fct_levels <- unique(
    c("Intercept", gtools::mixedsort(unique(beta_ranges$variable))))

  beta_ranges <- beta_ranges %>%
    dplyr::mutate(variable = forcats::fct(variable, levels = fct_levels))

  plot_subtitle <- stringr::str_glue(
    "Values of beta parameters (<span style='color:red'>minimum</span> and \\
    <span style='color:blue'>maximum</span> per species) across chains")

  # Construct path for saving the plot
  plot_path <- fs::path(
    dirname(model_dir), "model_postprocessing", "beta_min_max.jpeg")
  fs::dir_create(dirname(plot_path))

  beta_plot <- beta_ranges %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = (chain - 0.125), y = min_val), environment = emptyenv()) +
    ggplot2::geom_point(pch = 20, colour = "red", size = 0.8, alpha = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(x = (chain + 0.125), y = max_val), show.legend = FALSE,
      pch = 20, colour = "blue", size = 0.8, alpha = 0.5) +
    ggplot2::scale_x_continuous(breaks = ecokit::integer_breaks(n_chains)) +
    ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 5, nrow = 4) +
    ggplot2::labs(
      title = "Convergence of beta parameters", subtitle = plot_subtitle,
      x = "chain",
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
  print(beta_plot)
  grDevices::dev.off()

  invisible(NULL)
}
