## |------------------------------------------------------------------------| #
# PlotRho ----
## |------------------------------------------------------------------------| #

#' Plot convergence traceplots for the rho parameter
#'
#' This function generates and plots convergence traceplots for the rho
#' parameter of an Hmsc model. It visualizes the trace and density of the rho
#' parameter across different chains, providing insights into the convergence
#' and distribution of the parameter estimates.
#' @param Post A `coda` object containing MCMC samples of the rho parameter or a
#'   character string specifying the path to such an object.
#' @param Model A fitted Hmsc model object or a character string specifying the
#'   path to such an object.
#' @param Title A character string specifying the title of the plot.
#' @param Cols A character vector specifying the colors to be used for the lines
#'   representing each chain in the plot. Defaults to c("red", "blue",
#'   "darkgreen", "darkgrey").
#' @name PlotRho
#' @author Ahmed El-Gabbas
#' @return A ggplot object representing the traceplot of the rho parameter,
#'   including annotations for the Gelman-Rubin diagnostic, effective sample
#'   size, and credible intervals.
#' @export

PlotRho <- function(
    Post, Model, Title, Cols = c("red", "blue", "darkgreen", "darkgrey")) {

  if (is.null(Post) || is.null(Model) || is.null(Title)) {
    stop("Post, Model, and Title cannot be empty", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chain <- ID <- Value <- x <- y <- label <- NULL

  # Load coda object
  if (inherits(Post, "character")) {
    Post <- IASDT.R::LoadAs(Post)
  }
  Post <- Post$Rho

  # Load model object
  if (inherits(Model, "character")) {
    Model <- IASDT.R::LoadAs(Model)
  }

  SampleSize <- Model$samples
  NChains <- length(Model$postList)
  rm(Model)

  ## Gelman convergence diagnostic
  Gelman <- coda::gelman.diag(x = Post, multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    magrittr::extract(1) %>%
    round(3) %>%
    paste0("<i>Gelman convergence diagnostic:</i> ", .)

  ## Effective sample size
  ESS <- coda::effectiveSize(Post) %>%
    magrittr::divide_by(NChains) %>%
    round(1) %>%
    paste0("<i>Mean effective sample size:</i> ", ., " / ", SampleSize)

  ## quantiles
  CI <- summary(Post, quantiles = c(0.25, 0.75))$quantiles
  CI2 <- paste0("<i>50% credible interval:</i> ", CI, collapse = " - ")

  RhoDF <- purrr::map(.x = Post, .f = tibble::as_tibble, rownames = "ID") %>%
    dplyr::bind_rows(.id = "Chain") %>%
    dplyr::rename(Value = "var1") %>%
    dplyr::mutate(Chain = factor(Chain), ID = as.integer(ID))

  Title2 <- data.frame(
    x = Inf, y = Inf,
    label = paste0('<b><span style="color:blue">',
                   Title, "</span></b><br>", Gelman))
  ESS_CI <- data.frame(
    x = -Inf, y = -Inf, label = paste0(ESS, "<br>", CI2))

  Plot <- ggplot2::ggplot(
    data = RhoDF,
    mapping = ggplot2::aes(x = ID, y = Value, color = factor(Chain))) +
    ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0) +
    ggplot2::geom_hline(
      yintercept = CI, linetype = "dashed", color = "black", linewidth = 1) +
    ggplot2::scale_color_manual(values = Cols) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = Title2, inherit.aes = FALSE, size = 6,
      hjust = 1, vjust = 1, lineheight  = 0, fill = NA, label.color = NA) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = ESS_CI, inherit.aes = FALSE, size = 6,
      hjust = 0, vjust = 0, lineheight  = 0, fill = NA, label.color = NA) +
    ggplot2::theme(
      legend.position = "none", axis.text = ggplot2::element_text(size = 14))

  Plot1 <- Plot %>%
    ggExtra::ggMarginal(
      type = "density", margins = "y", size = 5, color = "steelblue4")

  return(Plot1)
}
