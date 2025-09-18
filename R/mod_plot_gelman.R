## |------------------------------------------------------------------------| #
# plot_gelman ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks
#'
#' The `plot_gelman_*()` functions generate plots visualising the evolution of
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
#' @param path_coda Character. Path to a file containing the coda object,
#'   representing MCMC samples.
#' @param alpha,beta,omega,rho Logical. If `TRUE`, plots the Gelman-Rubin
#'   statistic for the respective model parameters (alpha, beta, omega, or rho).
#'   Default: `TRUE` for all parameters.
#' @param n_omega Integer. Number of species sampled for the omega parameter.
#'   Default: 1000L.
#' @param plotting_alpha Numeric. Transparency level (alpha) for plot lines (0 =
#'   fully transparent, 1 = fully opaque). Default: 0.25.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param coda_object `mcmc.list`. An MCMC sample object containing posterior
#'   distributions from an Hmsc model.
#' @rdname plot_gelman
#' @name plot_gelman
#' @order 1
#' @export
#' @author Ahmed El-Gabbas

plot_gelman <- function(
    path_coda = NULL, alpha = TRUE, beta = TRUE, omega = TRUE, rho = TRUE,
    n_omega = 1000L, plotting_alpha = 0.25, env_file = ".env") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments --------
  ecokit::check_args(args_to_check = "path_coda", args_type = "character")
  ecokit::check_args(
    args_to_check = c("n_omega", "plotting_alpha"), args_type = "numeric")
  ecokit::check_args(
    args_to_check = c("beta", "rho", "omega", "alpha"), args_type = "logical")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (sum(alpha, beta, omega, rho) == 0) {
    ecokit::stop_ctx(
      "At least one of `alpha`, `beta`, `omega`, and `rho` must be `TRUE`",
      alpha = alpha, beta = beta, omega = omega, rho = rho,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Loading coda object ------

  ecokit::cat_time("Loading coda object")
  coda_object <- ecokit::load_as(path_coda)
  names_coda <- names(coda_object)

  # # ..................................................................... ###

  out_path <- fs::path(dirname(dirname(path_coda)), "model_convergence")
  fs::dir_create(out_path)

  # # ..................................................................... ###

  # alpha -----

  if (alpha && ("Alpha" %in% names_coda)) {
    ecokit::cat_time("alpha")
    plot_obj_alpha <- plot_gelman_alpha(
      coda_object = coda_object$Alpha[[1]], plotting_alpha = plotting_alpha)
  } else {
    plot_obj_alpha <- NULL
  }

  # # ..................................................................... ###

  # beta -----

  if (beta) {
    ecokit::cat_time("beta")
    plot_obj_beta <- IASDT.R::plot_gelman_beta(
      coda_object = coda_object$Beta, env_file = env_file,
      plotting_alpha = plotting_alpha)
  } else {
    plot_obj_beta <- NULL
  }

  # # ..................................................................... ###

  # omega -----

  if (omega && ("Omega" %in% names_coda)) {
    ecokit::cat_time("omega")
    plot_obj_omega <- IASDT.R::plot_gelman_omega(
      coda_object = coda_object$Omega[[1]], n_omega = n_omega,
      plotting_alpha = plotting_alpha)
  } else {
    plot_obj_omega <- NULL
  }

  # # ..................................................................... ###

  # rho -----

  if (rho && ("Rho" %in% names_coda)) {
    ecokit::cat_time("rho")
    plot_obj_rho <- IASDT.R::plot_gelman_rho(coda_object$Rho)
  } else {
    plot_obj_rho <- NULL
  }

  # # ..................................................................... ###

  # Saving plots as PDF -----

  plot_list <- list(
    Alpha = plot_obj_alpha, Beta = plot_obj_beta,
    Omega = plot_obj_omega, Rho = plot_obj_rho)

  plot_list_for_plot <- purrr::list_flatten(purrr::discard(plot_list, is.null))

  if (length(plot_list_for_plot) > 0) {
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = fs::path(out_path, "gelman_plots.pdf"),
      width = 13, height = 7, onefile = TRUE)
    purrr::walk(plot_list_for_plot, grid::grid.draw)
    grDevices::dev.off()
  } else {
    warning("No plots to save", call. = FALSE)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(init_time = .start_time)

  return(invisible(NULL))
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

  ecokit::check_args(args_to_check = "plotting_alpha", args_type = "numeric")

  if (!inherits(coda_object, "mcmc.list")) {
    ecokit::stop_ctx(
      "`coda_object` has to be of class mcmc.list",
      coda_object = coda_object, class_coda_object = class(coda_object),
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  var_lv <- type <- iter <- plot <- data <- median <- NULL

  # # ..................................................................... ###

  gelman_alpha_data <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        alpha_pre_plot_1 <- lapply(coda_object, function(y) {
          y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        alpha_pre_plot_2 <- try(
          gelman_preplot(
            x = alpha_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(alpha_pre_plot_2, "try-error")) {
          alpha_pre_plot_2 <- gelman_preplot(
            x = alpha_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        alpha_pre_plot_2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "iter") %>%
          purrr::set_names(c("iter", "median", "q_97_5")) %>%
          dplyr::filter(!is.nan(median)) %>%
          dplyr::mutate(iter = as.integer(iter)) %>%
          tidyr::pivot_longer(
            cols = -iter, names_to = "type", values_to = "shrink_factor") %>%
          dplyr::arrange(type, iter) %>%
          dplyr::mutate(type = factor(type), var_lv = x)
      })

  gelman_alpha_plot <- gelman_alpha_data %>%
    dplyr::mutate(
      group = paste0(var_lv, "_", type),
      var_lv = purrr::map_chr(
        .x = var_lv, .f = stringr::str_remove_all,
        pattern = "Alpha1\\[|\\]")) %>%
    tidyr::nest(.by = "var_lv") %>%
    dplyr::mutate(
      plot = purrr::map2(
        .x = var_lv, .y = data,
        .f = ~{
          plot <- ggplot2::ggplot(data = .y, environment = emptyenv()) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(
                x = iter, y = shrink_factor, group = group, color = type),
              alpha = plotting_alpha) +
            ggplot2::scale_color_manual(
              values = c(median = "red", q_97_5 = "black")) +
            ggplot2::geom_hline(
              yintercept = 1.1, linetype = "dashed", col = "darkgrey",
              linewidth = 0.8) +
            ggplot2::facet_grid(~ type, labeller = ggplot2::label_parsed) +
            ggplot2::scale_x_continuous(
              limits = range(gelman_alpha_data$iter), expand = c(0, 0)) +
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
          plot

        })
    ) %>%
    dplyr::pull(plot)

  return(gelman_alpha_plot)
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

  ecokit::check_args(args_to_check = "plotting_alpha", args_type = "numeric")

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
  iter <- type <- var_sp <- shrink_factor <- group <- NULL

  # # ..................................................................... ###

  beta_coda <- IASDT.R::coda_to_tibble(
    coda_object = coda_object, posterior_type = "beta", env_file = env_file)

  n_vars <- length(unique(beta_coda$variable))
  n_sp <- length(unique(beta_coda$species))
  plot_subtitle <- paste0(n_vars, " covariates - ", n_sp, " species")
  rm(beta_coda, envir = environment())

  gelman_beta_vals <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        beta_pre_plot_1 <- lapply(coda_object, function(y) {
          y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        beta_pre_plot_2 <- try(
          gelman_preplot(
            x = beta_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(beta_pre_plot_2, "try-error")) {
          beta_pre_plot_2 <- gelman_preplot(
            x = beta_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        beta_pre_plot_2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "iter") %>%
          purrr::set_names(c("iter", "median", "q_97_5")) %>%
          dplyr::mutate(iter = as.integer(iter)) %>%
          tidyr::pivot_longer(
            cols = -iter, names_to = "type", values_to = "shrink_factor") %>%
          dplyr::arrange(type, iter) %>%
          dplyr::mutate(type = factor(type), var_sp = x)
      }) %>%
    dplyr::mutate(group = paste0(var_sp, "_", type))

  # # ..................................................................... ###

  gelman_beta_plot <- ggplot2::ggplot(
    data = gelman_beta_vals, environment = emptyenv()) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = iter, y = shrink_factor, group = group, color = type),
      alpha = plotting_alpha) +
    ggplot2::scale_color_manual(values = c(median = "red", q_97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~type,
      labeller = ggplot2::as_labeller(c(median = "median", q_97_5 = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot --- beta", subtitle = plot_subtitle,
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

  return(gelman_beta_plot)
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

  ecokit::check_args(
    args_to_check = c("plotting_alpha", "n_omega"), args_type = "numeric")

  if (n_omega <= 0 || n_omega != as.integer(n_omega)) {
    ecokit::stop_ctx(
      "`n_omega` must be a positive integer.", n_omega = n_omega,
      include_backtrace = TRUE)
  }
  n_omega <- as.integer(n_omega)

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

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  iter <- type <- sp_comb <- shrink_factor <- group <- NULL

  # # ..................................................................... ###

  gelman_omega_data <- magrittr::extract2(coda_object, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sample(min(n_omega, length(.))) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        omega_pre_plot_1 <- lapply(coda_object, function(y) {
          y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        omega_pre_plot_2 <- try(
          gelman_preplot(
            x = omega_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(omega_pre_plot_2, "try-error")) {
          omega_pre_plot_2 <- gelman_preplot(
            x = omega_pre_plot_1,
            bin_width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        omega_pre_plot_2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "iter") %>%
          purrr::set_names(c("iter", "median", "q_97_5")) %>%
          dplyr::mutate(iter = as.integer(iter)) %>%
          tidyr::pivot_longer(
            cols = -iter, names_to = "type", values_to = "shrink_factor") %>%
          dplyr::arrange(type, iter) %>%
          dplyr::mutate(type = factor(type), sp_comb = x)
      }) %>%
    dplyr::mutate(group = paste0(sp_comb, "_", type))

  # # ..................................................................... ###

  gelman_omega_plot <- ggplot2::ggplot(
    data = gelman_omega_data, environment = emptyenv()) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = iter, y = shrink_factor, group = group, color = type),
      alpha = plotting_alpha) +
    ggplot2::scale_color_manual(values = c(median = "red", q_97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~type,
      labeller = ggplot2::as_labeller(c(median = "median", q_97_5 = "97.5%"))) +
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

  return(gelman_omega_plot)
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
  iter <- type <- shrink_factor <- NULL

  gelman_rho_plot <- try(
    gelman_preplot(
      x = coda_object, bin_width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = TRUE),
    silent = TRUE)

  if (inherits(gelman_rho_plot, "try-error")) {
    gelman_rho_plot <- gelman_preplot(
      x = coda_object, bin_width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = FALSE)
  }

  gelman_rho_plot <- gelman_rho_plot %>%
    magrittr::extract2("shrink") %>%
    tibble::as_tibble(rownames = "iter") %>%
    purrr::set_names(c("iter", "median", "q_97_5")) %>%
    dplyr::mutate(iter = as.integer(iter)) %>%
    tidyr::pivot_longer(
      cols = -iter, names_to = "type", values_to = "shrink_factor") %>%
    dplyr::arrange(type, iter) %>%
    dplyr::mutate(type = factor(type))
  gelman_rho_plot <- ggplot2::ggplot(
    gelman_rho_plot, environment = emptyenv()) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = iter, y = shrink_factor, color = type)) +
    ggplot2::scale_color_manual(values = c(median = "red", q_97_5 = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~type,
      labeller = ggplot2::as_labeller(c(median = "median", q_97_5 = "97.5%"))) +
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

  gelman_rho_plot
}
