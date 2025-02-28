## |------------------------------------------------------------------------| #
# Convergence_Rho ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Convergence_plots
#' @name Convergence_plots
#' @order 3
#' @author Ahmed El-Gabbas

Convergence_Rho <- function(
    Post = NULL, Model = NULL, Title = NULL,
    Cols = NULL, MarginType = "histogram") {

  # # ..................................................................... ###

  if (is.null(Post) || is.null(Model) || is.null(Title)) {
    stop("Post, Model, and Title cannot be empty", call. = FALSE)
  }

  if (length(MarginType) != 1) {
    stop("`MarginType` must be a single value.", call. = FALSE)
  }

  if (!MarginType %in% c("histogram", "density")) {
    stop(
      "`MarginType` must be either 'histogram' or 'density'.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chain <- ID <- Value <- x <- y <- label <- NULL

  # # ..................................................................... ###

  # Load coda object
  if (inherits(Post, "character")) {
    Post <- IASDT.R::LoadAs(Post)
  }
  Post <- Post$Rho

  # # ..................................................................... ###

  # Load model object
  if (inherits(Model, "character")) {
    Model <- IASDT.R::LoadAs(Model)
  }

  # # ..................................................................... ###

  SampleSize <- Model$samples
  NChains <- length(Model$postList)
  rm(Model, envir = environment())

  ## Effective sample size
  ESS <- coda::effectiveSize(Post) %>%
    magrittr::divide_by(NChains) %>%
    round(1) %>%
    paste0("<b><i>Mean effective sample size:</i></b> ", ., " / ", SampleSize)

  CI <- summary(Post, quantiles = c(0.025, 0.975))$quantiles
  CI2 <- paste0(
    "<b><i>95% credible interval:</i></b> ", paste(CI, collapse = " to "))

  RhoDF <- purrr::map(.x = Post, .f = tibble::as_tibble, rownames = "ID") %>%
    dplyr::bind_rows(.id = "Chain") %>%
    dplyr::rename(Value = "var1") %>%
    dplyr::mutate(Chain = factor(Chain), ID = as.integer(ID))

  ## Gelman convergence diagnostic
  Gelman <- try(
    coda::gelman.diag(x = Post, multivariate = FALSE),
    silent = TRUE)
  if (inherits(Gelman, "try-error")) {
    Gelman <- try(
      coda::gelman.diag(x = Post, multivariate = FALSE, autoburnin = FALSE),
      silent = TRUE)
  }
  Gelman <- Gelman %>%
    magrittr::extract2("psrf") %>%
    round(2) %>%
    paste(collapse = " / ") %>%
    paste0(
      '<span style="color:blue"><b><i>',
      "Gelman convergence diagnostic:</i></b></span> ", .)
  Title2 <- data.frame(x = Inf, y = Inf, label = Gelman)

  ESS_CI <- data.frame(
    x = -Inf, y = -Inf, label = paste0(ESS, "<br>", CI2))

  #  Plotting colours
  if (is.null(Cols)) {
    Cols <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }
  if (length(Cols) != NChains) {
    warning(
      "The length of provided colours != number of chains", call. = FALSE)
    Cols <- c(
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
    ggplot2::scale_color_manual(values = Cols) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0), limits = c(NA, NA), oob = scales::rescale_none) +
    ggplot2::theme_bw() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = Title2, inherit.aes = FALSE, size = 6,
      hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = ESS_CI, inherit.aes = FALSE, size = 6,
      hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none", axis.text = ggplot2::element_text(size = 14))

  if (MarginType == "histogram") {
    Plot1 <- ggExtra::ggMarginal(
      p = Plot, type = MarginType, margins = "y", size = 6,
      color = "steelblue4", fill = "steelblue4", bins = 100)
  } else {
    Plot1 <- ggExtra::ggMarginal(
      p = Plot, type = MarginType, margins = "y", size = 6,
      color = "steelblue4")
  }

  # Making marginal background matching the plot background
  # https://stackoverflow.com/a/78196022/3652584
  Plot1$layout$t[1] <- 1
  Plot1$layout$r[1] <- max(Plot1$layout$r)

  return(Plot1)
}
