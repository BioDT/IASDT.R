## |------------------------------------------------------------------------| #
# convergence_alpha ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 2
#' @author Ahmed El-Gabbas

convergence_alpha <- function(
    posterior = NULL, model_object = NULL, title = NULL, n_RC = NULL,
    add_footer = TRUE, add_title = TRUE, chain_colors = NULL,
    margin_type = "histogram") {

  temp_file <- fs::file_temp(ext = "pdf")
  grDevices::cairo_pdf(temp_file)
  on.exit({
    grDevices::dev.off()
    try(fs::file_delete(temp_file), silent = TRUE)
  },
  add = TRUE)

  if (is.null(posterior) || is.null(model_object)) {
    ecokit::stop_ctx(
      "`posterior` and `model_object` cannot be empty",
      model_object = model_object, posterior = posterior,
      include_backtrace = TRUE)
  }

  if (length(margin_type) != 1) {
    ecokit::stop_ctx(
      "`margin_type` must be a single string.",
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
  Factor <- Value <- PointEst <- UpperCI <- NULL

  # # ..................................................................... ###

  # Checking arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character", args_to_check = "title")

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


  # Load model object
  if (inherits(model_object, "character")) {
    model_object <- ecokit::load_as(model_object)
  }

  SampleSize <- model_object$samples      # nolint: object_name_linter
  NChains <- length(model_object$postList)

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

  rm(model_object, envir = environment())

  # Number of latent factors
  NLV <- ncol(posterior[[1]])

  # # ..................................................................... ###

  ## Gelman convergence diagnostic

  Gelman <- coda::gelman.diag(x = posterior, multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    as.data.frame() %>%
    stats::setNames(c("PointEst", "UpperCI")) %>%
    round(2) %>%
    dplyr::mutate(Gelman = paste0(PointEst, " / ", UpperCI)) %>%
    dplyr::pull(Gelman)

  ## Effective sample size
  ESS <- coda::effectiveSize(x = posterior) %>%  # nolint: object_name_linter
    magrittr::divide_by(NChains) %>%
    round(1)

  ## quantiles
  CI <- summary(posterior, quantiles = c(0.025, 0.975)) %>% # nolint: object_name_linter
    magrittr::extract2("quantiles") %>%
    matrix(ncol = 2) %>%
    magrittr::divide_by(1000) %>%
    round(3)

  AlphaDF <- IASDT.R::coda_to_tibble(      # nolint: object_name_linter
    coda_object = posterior, posterior_type = "alpha") %>%
    dplyr::mutate(
      Factor2 = purrr::map_int(
        .x = Factor,
        .f = ~ {
          as.character(.x) %>%
            stringr::str_remove("factor") %>%
            as.integer()
        }),
      Value = Value / 1000)

  if (is.null(n_RC)) {
    n_RC <- dplyr::case_when(
      NLV == 1 ~ c(1, 1), NLV == 2 ~ c(1, 2),
      NLV == 3 ~ c(1, 3), NLV == 4 ~ c(2, 2),
      .default = c(2, 3))
  }

  # # ..................................................................... ###

  Plots <- purrr::map(
    .x = seq_len(NLV),
    .f = ~ {
      ESS0 <- paste0(
        "<b><i>Mean effective sample size:</i></b> ", ESS[.x],
        " / ", SampleSize, " samples")

      CI0 <- paste(round(CI[.x, ], 2), collapse = " to ") %>%
        paste0("<b><i>95% credible interval:</i></b> ", ., " km")

      ESS_CI <- data.frame(
        x = -Inf, y = -Inf, label = paste0(ESS0, "<br>", CI0))

      Gelman0 <- paste0(
        "<b><i>Gelman convergence diagnostic:</i></b> ", Gelman[.x])

      title2 <- data.frame(x = Inf, y = Inf, label = Gelman0)
      title3 <- data.frame(
        x = -Inf, y = Inf, label = paste0("<b>Factor", .x, "</b>"))

      PlotDT <- dplyr::filter(AlphaDF, Factor2 == .x)

      # Pre-calculate smoothed lines
      summary_data <- PlotDT %>%
        group_by(Chain) %>%
        mutate(Smoothed = pmax(predict(loess(Value ~ Iter)), 0))

      Plot <- ggplot2::ggplot(
        data = PlotDT,
        mapping = ggplot2::aes(
          x = Iter, y = Value, color = Chain, group = Chain)) +
        ggplot2::geom_line(linewidth = 0.125, alpha = 0.6) +
        ggplot2::geom_line(
          data = summary_data, linewidth = 0.8, alpha = 0.6,
          mapping = ggplot2::aes(x = Iter, y = Smoothed, color = Chain)) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI[.x, ], linetype = "dashed", color = "black",
          linewidth = 1) +
        ggplot2::scale_color_manual(values = chain_colors) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
          limits = c(0, max(AlphaDF$Value) * 1.05),
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
          data = ESS_CI, inherit.aes = FALSE, size = 4,
          hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
        ggplot2::labs(x = NULL, y = "Distance (km)") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          text = ggplot2::element_text(family = "sans"),
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 12))

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

      return(Plot)
    })

  layout_matrix <- matrix(
    seq_len(n_RC[1] * n_RC[2]), nrow = n_RC[1], byrow = TRUE)

  if (add_title && add_footer) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = n_RC[1], ncol = n_RC[2], layout_matrix = layout_matrix)
  }

  if (add_title && isFALSE(add_footer)) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = NULL,
      top = grid::textGrob(
        label = title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = n_RC[1], ncol = n_RC[2], layout_matrix = layout_matrix)
  }

  if (isFALSE(add_title) && add_footer) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = NULL, nrow = n_RC[1], ncol = n_RC[2], layout_matrix = layout_matrix)
  }

  if (isFALSE(add_title) && isFALSE(add_footer)) {
    Plots <- cowplot::plot_grid(
      plotlist = Plots, ncol = n_RC[2], nrow = n_RC[1],
      align = "hv", byrow = TRUE)
  }

  return(Plots)
}
