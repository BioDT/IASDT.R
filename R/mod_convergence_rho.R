## |------------------------------------------------------------------------| #
# convergence_rho ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname convergence_plots
#' @name convergence_plots
#' @order 3
#' @author Ahmed El-Gabbas

convergence_rho <- function(
    posterior = NULL, model_object = NULL, title = NULL,
    chain_colors = NULL, margin_type = "histogram") {

  temp_file <- fs::file_temp(ext = "pdf")
  grDevices::cairo_pdf(temp_file)
  on.exit({
    grDevices::dev.off()
    try(fs::file_delete(temp_file), silent = TRUE)
  },
  add = TRUE)

  if (is.null(posterior) || is.null(model_object) || is.null(title)) {
    ecokit::stop_ctx(
      "`posterior`, `model_object`, and `title` cannot be empty",
      posterior = posterior, model_object = model_object, title = title,
      include_backtrace = TRUE)
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

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chain <- ID <- Value <- x <- y <- label <- NULL

  # # ..................................................................... ###

  # Load coda object
  if (inherits(posterior, "character")) {
    posterior <- ecokit::load_as(posterior)
  }
  posterior <- posterior$Rho

  # # ..................................................................... ###

  # Load model object
  if (inherits(model_object, "character")) {
    model_object <- ecokit::load_as(model_object)
  }

  # # ..................................................................... ###

  SampleSize <- model_object$samples
  NChains <- length(model_object$postList)
  rm(model_object, envir = environment())

  ## Effective sample size
  ESS <- coda::effectiveSize(posterior) %>%
    magrittr::divide_by(NChains) %>%
    round(1) %>%
    paste0("<b><i>Mean effective sample size:</i></b> ", ., " / ", SampleSize)

  CI <- summary(posterior, quantiles = c(0.025, 0.975))$quantiles
  CI2 <- paste0(
    "<b><i>95% credible interval:</i></b> ", paste(CI, collapse = " to "))

  RhoDF <- purrr::map(
    .x = posterior, .f = tibble::as_tibble, rownames = "ID") %>%
    dplyr::bind_rows(.id = "Chain") %>%
    dplyr::rename(Value = "var1") %>%
    dplyr::mutate(Chain = factor(Chain), ID = as.integer(ID))

  ## Gelman convergence diagnostic
  Gelman <- try(
    coda::gelman.diag(x = posterior, multivariate = FALSE),
    silent = TRUE)
  if (inherits(Gelman, "try-error")) {
    Gelman <- try(
      coda::gelman.diag(
        x = posterior, multivariate = FALSE, autoburnin = FALSE),
      silent = TRUE)
  }
  Gelman <- Gelman %>%
    magrittr::extract2("psrf") %>%
    round(2) %>%
    paste(collapse = " / ") %>%
    paste0(
      '<span style="color:blue"><b><i>',
      "Gelman convergence diagnostic:</i></b></span> ", .)
  title2 <- data.frame(x = Inf, y = Inf, label = Gelman)

  ESS_CI <- data.frame(
    x = -Inf, y = -Inf, label = paste0(ESS, "<br>", CI2))

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

  Plot <- ggplot2::ggplot(
    data = RhoDF,
    mapping = ggplot2::aes(x = ID, y = Value, color = factor(Chain))) +
    ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0) +
    ggplot2::geom_hline(
      yintercept = CI, linetype = "dashed", color = "black", linewidth = 1) +
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
      data = ESS_CI, inherit.aes = FALSE, size = 6,
      hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      legend.position = "none", axis.text = ggplot2::element_text(size = 14))

  if (margin_type == "histogram") {
    Plot1 <- ggExtra::ggMarginal(
      p = Plot, type = margin_type, margins = "y", size = 6,
      color = "steelblue4", fill = "steelblue4", bins = 100)
  } else {
    Plot1 <- ggExtra::ggMarginal(
      p = Plot, type = margin_type, margins = "y", size = 6,
      color = "steelblue4")
  }

  # Making marginal background matching the plot background
  # https://stackoverflow.com/a/78196022/3652584
  Plot1$layout$t[1] <- 1
  Plot1$layout$r[1] <- max(Plot1$layout$r)

  return(Plot1)
}
