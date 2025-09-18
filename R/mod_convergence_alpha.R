## |------------------------------------------------------------------------| #
# convergence_alpha ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 2
#' @author Ahmed El-Gabbas

convergence_alpha <- function(
    posterior = NULL, title = NULL, n_rc_alpha = c(2L, 3L), add_footer = TRUE,
    add_title = TRUE, chain_colors = NULL, margin_type = "histogram",
    n_chains = NULL, n_samples = NULL) {

  # Checking arguments
  ecokit::check_args(
    args_to_check = c("title", "margin_type"), args_type = "character")
  ecokit::check_args(
    args_to_check = c("n_rc_alpha", "n_samples", "n_chains"),
    args_type = "numeric", arg_length = c(2L, 1L, 1L))
  ecokit::check_args(
    args_to_check = c("add_footer", "add_title"), args_type = "logical")

  temp_file <- fs::file_temp(ext = "pdf")
  grDevices::cairo_pdf(temp_file)
  on.exit({
    grDevices::dev.off()
    try(fs::file_delete(temp_file), silent = TRUE)
  },
  add = TRUE)

  if (is.null(posterior)) {
    ecokit::stop_ctx(
      "`posterior` cannot be empty", posterior = posterior,
      include_backtrace = TRUE)
  }

  if (!margin_type %in% c("histogram", "density")) {
    ecokit::stop_ctx(
      "`margin_type` must be either 'histogram' or 'density'.",
      margin_type = margin_type, include_backtrace = TRUE)
  }

  if (n_chains < 1 || n_samples < 1) {
    ecokit::stop_ctx(
      "`n_chains` and `n_samples` should be positive numeric values",
      n_chains = n_chains, n_samples = n_samples)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  factor <- value <- point_est <- upper_ci <- NULL

  # # ..................................................................... ###

  # Load coda object
  if (inherits(posterior, "character")) {
    posterior <- ecokit::load_as(posterior)
  }

  if (!("Alpha" %in% names(posterior))) {
    ecokit::stop_ctx(
      "`posterior` object does not contain 'Alpha'",
      names_posterior = names(posterior), include_backtrace = TRUE)
  }
  posterior <- posterior$Alpha[[1]]

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

  # Number of latent factors
  n_lf <- ncol(posterior[[1]])

  # # ..................................................................... ###

  ## Gelman convergence diagnostic

  gelman <- coda::gelman.diag(x = posterior, multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    as.data.frame() %>%
    stats::setNames(c("point_est", "upper_ci")) %>%
    round(2) %>%
    dplyr::mutate(gelman = paste0(point_est, " / ", upper_ci)) %>%
    dplyr::pull(gelman)

  ## Effective sample size
  ess <- coda::effectiveSize(x = posterior) %>%  # nolint: object_name_linter
    magrittr::divide_by(n_chains) %>%
    round(1)

  ## quantiles
  CI <- summary(posterior, quantiles = c(0.025, 0.975)) %>% # nolint: object_name_linter
    magrittr::extract2("quantiles") %>%
    matrix(ncol = 2L) %>%
    magrittr::divide_by(1000L) %>%
    round(3)

  alpha_data <- IASDT.R::coda_to_tibble(      # nolint: object_name_linter
    coda_object = posterior, posterior_type = "alpha") %>%
    dplyr::mutate(
      factor2 = purrr::map_int(
        .x = factor,
        .f = ~ {
          as.character(.x) %>%
            stringr::str_remove("factor") %>%
            as.integer()
        }),
      value = value / 1000)

  if (is.null(n_rc_alpha)) {
    n_rc_alpha <- dplyr::case_when(
      n_lf == 1 ~ c(1, 1), n_lf == 2 ~ c(1, 2),
      n_lf == 3 ~ c(1, 3), n_lf == 4 ~ c(2, 2),
      .default = c(2, 3))
  }

  # # ..................................................................... ###

  plots <- purrr::map(
    .x = seq_len(n_lf),
    .f = ~ {

      temp_file <- fs::file_temp(ext = "pdf")
      grDevices::cairo_pdf(temp_file)
      on.exit({
        grDevices::dev.off()
        try(fs::file_delete(temp_file), silent = TRUE)
      },
      add = TRUE)

      ess_0 <- paste0(
        "<b><i>Mean effective sample size:</i></b> ", ess[.x],
        " / ", n_samples, " samples")

      ci_0 <- paste(round(CI[.x, ], 2), collapse = " to ") %>%
        paste0("<b><i>95% credible interval:</i></b> ", ., " km")

      ess_ci <- data.frame(
        x = -Inf, y = -Inf, label = paste0(ess_0, "<br>", ci_0))

      gelman_0 <- paste0(
        "<b><i>Gelman convergence diagnostic:</i></b> ", gelman[.x])

      title2 <- data.frame(x = Inf, y = Inf, label = gelman_0)
      title3 <- data.frame(
        x = -Inf, y = Inf, label = paste0("<b>Factor", .x, "</b>"))

      plot_data <- dplyr::filter(alpha_data, factor2 == .x)

      # Pre-calculate smoothed lines
      summary_data <- plot_data %>%
        dplyr::group_by(chain) %>%
        dplyr::mutate(Smoothed = pmax(predict(loess(value ~ iter)), 0))

      plot <- ggplot2::ggplot(
        data = plot_data, environment = emptyenv(),
        mapping = ggplot2::aes(
          x = iter, y = value, color = chain, group = chain)) +
        ggplot2::geom_line(linewidth = 0.125, alpha = 0.6) +
        ggplot2::geom_line(
          data = summary_data, linewidth = 0.8, alpha = 0.6,
          mapping = ggplot2::aes(x = iter, y = Smoothed, color = chain)) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI[.x, ], linetype = "dashed", color = "black",
          linewidth = 1) +
        ggplot2::scale_color_manual(values = chain_colors) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
          limits = c(0, max(alpha_data$value) * 1.05),
          expand = c(0, 0), oob = scales::squish) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), data = title2,
          inherit.aes = FALSE, size = 4, hjust = 1, vjust = 1, lineheight = 0,
          fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), data = title3,
          inherit.aes = FALSE, size = 5, hjust = -0.1, vjust = 1,
          color = "blue", lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = ess_ci, inherit.aes = FALSE, size = 4,
          hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
        ggplot2::labs(x = NULL, y = "Distance (km)") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          text = ggplot2::element_text(family = "sans"),
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 12))

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

      return(plot)
    })

  layout_matrix <- matrix(
    seq_len(n_rc_alpha[1] * n_rc_alpha[2]), nrow = n_rc_alpha[1], byrow = TRUE)

  if (add_title && add_footer) {
    plots <- gridExtra::marrangeGrob(
      plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = n_rc_alpha[1], ncol = n_rc_alpha[2], layout_matrix = layout_matrix)
  }

  if (add_title && isFALSE(add_footer)) {
    plots <- gridExtra::marrangeGrob(
      plots, bottom = NULL,
      top = grid::textGrob(
        label = title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = n_rc_alpha[1], ncol = n_rc_alpha[2], layout_matrix = layout_matrix)
  }

  if (isFALSE(add_title) && add_footer) {
    plots <- gridExtra::marrangeGrob(
      plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = NULL, nrow = n_rc_alpha[1], ncol = n_rc_alpha[2],
      layout_matrix = layout_matrix)
  }

  if (isFALSE(add_title) && isFALSE(add_footer)) {
    plots <- cowplot::plot_grid(
      plotlist = plots, ncol = n_rc_alpha[2], nrow = n_rc_alpha[1],
      align = "hv", byrow = TRUE)
  }

  return(plots)
}
