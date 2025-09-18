## |------------------------------------------------------------------------| #
# convergence_rho ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 3
#' @author Ahmed El-Gabbas

convergence_rho <- function(
    posterior = NULL, title = NULL, chain_colors = NULL,
    margin_type = "histogram", n_chains = NULL, n_samples = NULL) {

  temp_file <- fs::file_temp(ext = "pdf")
  grDevices::cairo_pdf(temp_file)
  on.exit({
    grDevices::dev.off()
    try(fs::file_delete(temp_file), silent = TRUE)
  },
  add = TRUE)

  if (is.null(posterior) || is.null(title)) {
    ecokit::stop_ctx(
      "`posterior` and `title` cannot be empty",
      posterior = posterior, title = title, include_backtrace = TRUE)
  }

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

  if (!is.numeric(n_chains) || !is.numeric(n_samples) ||
      length(n_chains) != 1 || length(n_samples) != 1 ||
      n_chains < 1 || n_samples < 1) {
    ecokit::stop_ctx(
      "n_chains and n_samples should be positive numeric vectors of length 1",
      n_chains = n_chains, n_samples = n_samples)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  chain <- id <- value <- x <- y <- label <- NULL

  # # ..................................................................... ###

  # Load coda object
  if (inherits(posterior, "character")) {
    posterior <- ecokit::load_as(posterior)
  }
  posterior <- posterior$Rho

  # # ..................................................................... ###

  ## Effective sample size
  ess <- coda::effectiveSize(posterior) %>%
    magrittr::divide_by(n_chains) %>%
    round(1) %>%
    paste0("<b><i>Mean effective sample size:</i></b> ", ., " / ", n_samples)

  ci <- summary(posterior, quantiles = c(0.025, 0.975))$quantiles
  ci_2 <- paste0(
    "<b><i>95% credible interval:</i></b> ", paste(ci, collapse = " to "))

  rho_data <- purrr::map(
    .x = posterior, .f = tibble::as_tibble, rownames = "id") %>%
    dplyr::bind_rows(.id = "chain") %>%
    dplyr::rename(value = "var1") %>%
    dplyr::mutate(chain = factor(chain), id = as.integer(id))

  ## Gelman convergence diagnostic
  gelman <- try(
    coda::gelman.diag(x = posterior, multivariate = FALSE),
    silent = TRUE)
  if (inherits(gelman, "try-error")) {
    gelman <- try(
      coda::gelman.diag(
        x = posterior, multivariate = FALSE, autoburnin = FALSE),
      silent = TRUE)
  }
  gelman <- gelman %>%
    magrittr::extract2("psrf") %>%
    round(2) %>%
    paste(collapse = " / ") %>%
    paste0(
      '<span style="color:blue"><b><i>',
      "Gelman convergence diagnostic:</i></b></span> ", .)
  title2 <- data.frame(x = Inf, y = Inf, label = gelman)

  ess_ci <- data.frame(
    x = -Inf, y = -Inf, label = paste0(ess, "<br>", ci_2))

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

  plot <- ggplot2::ggplot(
    data = rho_data, environment = emptyenv(),
    mapping = ggplot2::aes(x = id, y = value, color = factor(chain))) +
    ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0) +
    ggplot2::geom_hline(
      yintercept = ci, linetype = "dashed", color = "black", linewidth = 1) +
    # Ensure that y-axis always show 0
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(values = chain_colors) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0), limits = c(NA, NA), oob = scales::rescale_none) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = title2, inherit.aes = FALSE, size = 6,
      hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = ess_ci, inherit.aes = FALSE, size = 6,
      hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      legend.position = "none", axis.text = ggplot2::element_text(size = 14))

  if (margin_type == "histogram") {
    plot_1 <- ggExtra::ggMarginal(
      p = plot, type = margin_type, margins = "y", size = 6,
      color = "steelblue4", fill = "steelblue4", bins = 100)
  } else {
    plot_1 <- ggExtra::ggMarginal(
      p = plot, type = margin_type, margins = "y", size = 6,
      color = "steelblue4")
  }

  # Making marginal background matching the plot background
  # https://stackoverflow.com/a/78196022/3652584
  plot_1$layout$t[1] <- 1
  plot_1$layout$r[1] <- max(plot_1$layout$r)

  return(plot_1)
}
