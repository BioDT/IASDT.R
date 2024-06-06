## |------------------------------------------------------------------------| #
# Plot_Rho ----
## |------------------------------------------------------------------------| #

#' Plot convergence traceplots for the rho parameter
#'
#' Plot convergence traceplots for the rho parameter
#'
#' @param Post Coda object or path to it
#' @param Model Fitted model object or path to it
#' @param Title String. Plotting title
#' @param Cols Colours for lines for each chain
#' @name Plot_Rho
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Plot_Rho <- function(
    Post = NULL, Model = NULL, Title = NULL,
    Cols = c("red", "blue", "darkgreen", "darkgrey")) {

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
  invisible(gc())

  ## Gelman convergence diagnostic
  Gelman <- Post %>%
    coda::gelman.diag(multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    magrittr::extract(1) %>%
    round(3) %>%
    paste0("<i>Gelman convergence diagnostic:</i> ", .)

  ## Effective sample size
  ESS <- Post %>%
    coda::effectiveSize() %>%
    magrittr::divide_by(NChains) %>%
    round(1) %>%
    paste0("<i>Mean effective sample size:</i> ", ., " / ", SampleSize)

  ## quantiles
  CI <- summary(Post, quantiles = c(0.25, 0.75))$quantiles
  CI2 <- CI %>%
    paste0(collapse = " - ") %>%
    paste0("<i>50% credible interval:</i> ", .)

  RhoDF <- Post %>%
    purrr::map(tibble::as_tibble, rownames = "ID") %>%
    dplyr::bind_rows(.id = "Chain") %>%
    dplyr::rename(Value = "var1") %>%
    dplyr::mutate(Chain = factor(Chain), ID = as.integer(ID))

  Title2 <- data.frame(
    x = Inf, y = Inf,
    label = paste0('<b><span style="color:blue">',
                   Title, "</span></b><br>", Gelman))
  ESS_CI <- data.frame(
    x = -Inf, y = -Inf, label = paste0(ESS, "<br>", CI2))

  Plot <- RhoDF %>%
    ggplot2::ggplot(
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
      data = Title2, inherit.aes = FALSE,
      hjust = 1, vjust = 1, lineheight  = 0, fill = NA, label.color = NA) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(x = x, y = y, label = label),
      data = ESS_CI, inherit.aes = FALSE,
      hjust = 0, vjust = 0, lineheight  = 0, fill = NA, label.color = NA) +
    ggplot2::theme(legend.position = "none")

  Plot1 <- Plot %>%
    ggExtra::ggMarginal(
      type = "density", margins = "y", size = 5, color = "steelblue4")
  return(Plot1)
}
