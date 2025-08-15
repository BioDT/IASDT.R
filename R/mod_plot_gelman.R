## |------------------------------------------------------------------------| #
# plot_gelman ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks
#'
#' The `PlotGelman_*()` functions generate plots visualising the evolution of
#' the Gelman-Rubin-Brooks shrink factor for different model parameters as the
#' number of iterations increases. These plots help assess whether MCMC chains
#' have converged to a common distribution. Each plot includes: median (solid
#' line) and 97.5<sup>th</sup> percentile (dashed line) of the shrink factor and
#' a dashed horizontal line at 1.1, representing the common convergence
#' threshold. The primary function for users is `plot_gelman()`, which
#' internally calls:
#' - `plot_gelman_alpha()`: Plots shrink factor for the **Alpha** parameter
#' - `plot_gelman_beta()`: Plots shrink factor for the **Beta** parameters
#' - `plot_gelman_omega()`: Plots shrink factor for the **Omega** parameter
#' - `plot_gelman_rho()`: Plots shrink factor for the **Rho** parameter
#' @param path_coda Character or `mcmc.list`. Path to a file containing the coda
#'   object, or an `mcmc.list` object representing MCMC samples.
#' @param alpha,beta,omega,rho Logical. If `TRUE`, plots the Gelman-Rubin
#'   statistic for the respective model parameters (alpha, beta, omega, or rho).
#'   Default: `TRUE` for all parameters.
#' @param n_omega Integer. Number of species sampled for the omega parameter.
#'   Default: 1000L.
#' @param plotting_alpha Numeric. Transparency level (alpha) for plot lines (0 =
#'   fully transparent, 1 = fully opaque). Default: 0.25.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param return_plots Character. Path to the folder where the output plots will
#'   be saved.
#' @param coda_object `mcmc.list`. An MCMC sample object containing posterior
#'   distributions from an Hmsc model.
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 1
#' @export
#' @author Ahmed El-Gabbas

plot_gelman <- function(
    path_coda = NULL, alpha = TRUE, beta = TRUE, omega = TRUE, rho = TRUE,
    n_omega = 1000L, plotting_alpha = 0.25, env_file = ".env",
    return_plots = FALSE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments --------

  if (sum(alpha, beta, omega, rho) == 0) {
    ecokit::stop_ctx(
      "At least one of `alpha`, `beta`, `omega`, and `rho` must be `TRUE`",
      alpha = alpha, beta = beta, omega = omega, rho = rho,
      include_backtrace = TRUE)
  }

  if (is.null(path_coda)) {
    ecokit::stop_ctx(
      "path_coda cannot be empty", path_coda = path_coda,
      include_backtrace = TRUE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_omega", "plotting_alpha"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("beta", "rho", "omega", "alpha", "return_plots"))

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Loading coda object ------

  if (inherits(path_coda, "character")) {

    ecokit::cat_time("Loading coda object")
    coda_object <- ecokit::load_as(path_coda)

  } else {

    if (!inherits(path_coda, "list")) {
      ecokit::stop_ctx(
        "`path_coda` is neither character path or a list",
        path_coda = path_coda, class_path_coda = class(path_coda),
        include_backtrace = TRUE)
    }
    if (!inherits(path_coda[[1]], "mcmc.list")) {
      ecokit::stop_ctx(
        "`path_coda` has no mcmc.list items",
        path_coda = path_coda, class_path_coda = class(path_coda),
        include_backtrace = TRUE)
    }

    coda_object <- path_coda
    rm(path_coda, envir = environment())
  }
  names_coda <- names(coda_object)

  # # ..................................................................... ###

  out_path <- fs::path(dirname(dirname(path_coda)), "Model_Convergence")
  fs::dir_create(out_path)

  # # ..................................................................... ###

  # alpha -----

  if (alpha && ("Alpha" %in% names_coda)) {
    ecokit::cat_time("alpha")
    PlotObj_Alpha <- plot_gelman_alpha(
      coda_object = coda_object$Alpha[[1]], plotting_alpha = plotting_alpha)
  } else {
    PlotObj_Alpha <- NULL
  }

  # # ..................................................................... ###

  # beta -----

  if (beta) {
    ecokit::cat_time("beta")
    PlotObj_Beta <- IASDT.R::plot_gelman_beta(
      coda_object = coda_object$Beta, env_file = env_file,
      plotting_alpha = plotting_alpha)
  } else {
    PlotObj_Beta <- NULL
  }

  # # ..................................................................... ###

  # omega -----

  if (omega && ("Omega" %in% names_coda)) {
    ecokit::cat_time("omega")
    PlotObj_Omega <- IASDT.R::plot_gelman_omega(
      coda_object = coda_object$Omega[[1]], n_omega = n_omega,
      plotting_alpha = plotting_alpha)
  } else {
    PlotObj_Omega <- NULL
  }

  # # ..................................................................... ###

  # rho -----

  if (rho && ("Rho" %in% names_coda)) {
    ecokit::cat_time("rho")
    PlotObj_Rho <- IASDT.R::plot_gelman_rho(coda_object$Rho)
  } else {
    PlotObj_Rho <- NULL
  }

  # # ..................................................................... ###

  # Saving plots as PDF -----

  PlotList <- list(
    Alpha = PlotObj_Alpha, Beta = PlotObj_Beta,
    Omega = PlotObj_Omega, Rho = PlotObj_Rho)

  PlotList4Plot <- purrr::list_flatten(purrr::discard(PlotList, is.null))

  if (length(PlotList4Plot) > 0) {
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(out_path, "GelmanPlots.pdf"),
      width = 13, height = 7, onefile = TRUE)
    purrr::walk(PlotList4Plot, grid::grid.draw)
    grDevices::dev.off()
  } else {
    warning("No plots to save", call. = FALSE)
  }

  # # ..................................................................... ###

  # Saving plots as qs2 -----
  ecokit::save_as(
    object = PlotList, object_name = "GelmanPlots", n_threads = 5L,
    out_path = fs::path(out_path, "GelmanPlots.qs2"))

  # # ..................................................................... ###

  ecokit::cat_diff(init_time = .start_time)

  if (return_plots) {
    return(PlotList)
  } else {
    return(invisible(NULL))
  }
}

# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# plot_gelman_alpha ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 2
#' @author Ahmed El-Gabbas

plot_gelman_alpha <- function(coda_object, plotting_alpha = 0.25) {

  # # ..................................................................... ###

  if (!inherits(coda_object, "mcmc.list")) {
    ecokit::stop_ctx(
      "`coda_object` has to be of class mcmc.list",
      coda_object = coda_object, class_coda_object = class(coda_object),
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var_LV <- Type <- Iter <- Plot <- data <- Median <- NULL

  # # ..................................................................... ###

  Gelman_Alpha_DT <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Alpha_Preplot1 <- lapply(coda_object, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Alpha_Preplot2 <- try(
          gelman_preplot(
            x = Alpha_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Alpha_Preplot2, "try-error")) {
          Alpha_Preplot2 <- gelman_preplot(
            x = Alpha_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Alpha_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::filter(!is.nan(Median)) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Var_LV = x)
      })

  Gelman_Alpha_Plot <- Gelman_Alpha_DT %>%
    dplyr::mutate(
      group = paste0(Var_LV, "_", Type),
      Var_LV = purrr::map_chr(
        .x = Var_LV, .f = stringr::str_remove_all,
        pattern = "Alpha1\\[|\\]")) %>%
    tidyr::nest(.by = "Var_LV") %>%
    dplyr::mutate(
      Plot = purrr::map2(
        .x = Var_LV, .y = data,
        .f = ~{
          ggplot2::ggplot(data = .y) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(
                x = Iter, y = ShrinkFactor, group = group, color = Type),
              alpha = plotting_alpha) +
            ggplot2::scale_color_manual(
              values = c(Median = "red", Q97_5 = "black")) +
            ggplot2::geom_hline(
              yintercept = 1.1, linetype = "dashed", col = "darkgrey",
              linewidth = 0.8) +
            ggplot2::facet_grid(~ Type, labeller = ggplot2::label_parsed) +
            ggplot2::scale_x_continuous(
              limits = range(Gelman_Alpha_DT$Iter), expand = c(0, 0)) +
            ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
            ggplot2::labs(
              title = "Gelman-Rubin-Brooks plot --- alpha",
              subtitle = .x,
              caption = paste0(
                "This plot shows the evolution of Gelman and Rubin's shrink ",
                "factor as the number of iterations increases.")) +
            ggplot2::xlab(NULL) +
            ggplot2::ylab("Shrink factor") +
            ggplot2::theme_bw() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(5, 20, 5, 5),
              legend.position = "none",
              legend.background = ggplot2::element_blank(),
              legend.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(size = 12, face = "bold"),
              axis.title = ggplot2::element_text(
                size = 12, colour = "darkgrey", face = "bold"),
              axis.text = ggplot2::element_text(size = 12),
              plot.title = ggtext::element_markdown(
                size = 18, face = "bold", color = "blue"),
              plot.subtitle = ggplot2::element_text(
                size = 12, face = "italic", color = "darkgrey"),
              panel.spacing = ggplot2::unit(0.85, "lines"),
              plot.caption = ggplot2::element_text(
                color = "darkgrey", face = "italic", size = 10))
        })
    ) %>%
    dplyr::pull(Plot)

  return(Gelman_Alpha_Plot)
}


# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# plot_gelman_beta ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 3
#' @author Ahmed El-Gabbas

plot_gelman_beta <- function(
    coda_object, env_file = ".env", plotting_alpha = 0.25) {

  # # ..................................................................... ###

  if (is.null(coda_object)) {
    ecokit::stop_ctx(
      "`coda_object` cannot be empty", coda_object = coda_object,
      include_backtrace = TRUE)
  }

  if (!inherits(coda_object, "mcmc.list")) {
    ecokit::stop_ctx(
      "`coda_object` has to be of class mcmc.list",
      coda_object = coda_object, class_coda_object = class(coda_object),
      include_backtrace = TRUE)
  }

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Var_Sp <- ShrinkFactor <- group <- NULL

  # # ..................................................................... ###

  Beta_Coda <- IASDT.R::coda_to_tibble(
    coda_object = coda_object, posterior_type = "beta", env_file = env_file)

  NVars <- length(unique(Beta_Coda$Variable))
  NSp <- length(unique(Beta_Coda$Species))
  SubTitle <- paste0(NVars, " covariates - ", NSp, " species")
  rm(Beta_Coda, envir = environment())

  Gelman_Beta_Vals <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Beta_Preplot1 <- lapply(coda_object, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Beta_Preplot2 <- try(
          gelman_preplot(
            x = Beta_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Beta_Preplot2, "try-error")) {
          Beta_Preplot2 <- gelman_preplot(
            x = Beta_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Beta_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Var_Sp = x)
      }) %>%
    dplyr::mutate(group = paste0(Var_Sp, "_", Type))

  # # ..................................................................... ###

  Gelman_Beta_Plot <- ggplot2::ggplot(data = Gelman_Beta_Vals) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = plotting_alpha) +
    ggplot2::scale_color_manual(values = c(Median = "red", Q97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(c(Median = "Median", Q97_5 = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot --- beta", subtitle = SubTitle,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's shrink factor as ",
        "the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme(
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  # # ..................................................................... ###

  return(Gelman_Beta_Plot)
}

# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# plot_gelman_omega ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 4
#' @author Ahmed El-Gabbas

plot_gelman_omega <- function(
    coda_object, n_omega = 1000L, plotting_alpha = 0.25) {

  # # ..................................................................... ###

  if (is.null(coda_object)) {
    ecokit::stop_ctx(
      "`coda_object` cannot be empty", coda_object = coda_object,
      include_backtrace = TRUE)
  }

  if (!is.numeric(n_omega) || n_omega <= 0 || n_omega != as.integer(n_omega)) {
    ecokit::stop_ctx(
      "`n_omega` must be a positive integer.", n_omega = n_omega,
      include_backtrace = TRUE)
  }
  n_omega <- as.integer(n_omega)

  if (!inherits(coda_object, "mcmc.list")) {
    ecokit::stop_ctx(
      "`coda_object` has to be of class mcmc.list",
      coda_object = coda_object, class_coda_object = class(coda_object),
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Sp_comb <- ShrinkFactor <- group <- NULL

  # # ..................................................................... ###

  Gelman_OmegaDT <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sample(min(n_omega, length(.))) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Omega_Preplot1 <- lapply(coda_object, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Omega_Preplot2 <- try(
          gelman_preplot(
            x = Omega_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Omega_Preplot2, "try-error")) {
          Omega_Preplot2 <- gelman_preplot(
            x = Omega_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Omega_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Sp_comb = x)
      }) %>%
    dplyr::mutate(group = paste0(Sp_comb, "_", Type))

  # # ..................................................................... ###

  Gelman_Omega_Plot <- ggplot2::ggplot(data = Gelman_OmegaDT) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = plotting_alpha) +
    ggplot2::scale_color_manual(values = c(Median = "red", Q97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(c(Median = "Median", Q97_5 = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      title = paste0(
        "Gelman-Rubin-Brooks plot --- omega --- ",
        n_omega, " species combination samples"),
      subtitle = NULL,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's shrink factor as ",
        "the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 20, 5, 5),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  return(Gelman_Omega_Plot)
}


# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# plot_gelman_rho ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 5
#' @author Ahmed El-Gabbas

plot_gelman_rho <- function(coda_object) {

  if (is.null(coda_object)) {
    ecokit::stop_ctx(
      "`coda_object` cannot be empty", coda_object = coda_object,
      include_backtrace = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- ShrinkFactor <- NULL

  Gelman_Rho_Plot <- try(
    gelman_preplot(
      x = coda_object, bin.width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = TRUE),
    silent = TRUE)

  if (inherits(Gelman_Rho_Plot, "try-error")) {
    Gelman_Rho_Plot <- gelman_preplot(
      x = coda_object, bin.width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = FALSE)
  }

  Gelman_Rho_Plot <- Gelman_Rho_Plot %>%
    magrittr::extract2("shrink") %>%
    tibble::as_tibble(rownames = "Iter") %>%
    purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
    dplyr::mutate(Iter = as.integer(Iter)) %>%
    tidyr::pivot_longer(
      cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
    dplyr::arrange(Type, Iter) %>%
    dplyr::mutate(Type = factor(Type)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = Iter, y = ShrinkFactor, color = Type)) +
    ggplot2::scale_color_manual(values = c(Median = "red", Q97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(c(Median = "Median", Q97_5 = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot --- rho",
      subtitle = NULL,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's ",
        "shrink factor as the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 20, 5, 5),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  return(Gelman_Rho_Plot)
}
