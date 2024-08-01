## |------------------------------------------------------------------------| #
# PlotGelman_Rho ----
## |------------------------------------------------------------------------| #

#' Creates a Gelman-Rubin-Brooks plot for the `rho` parameter.
#'
#' This function generates a Gelman-Rubin-Brooks plot for the `rho` parameter of
#' a coda object. The plot includes lines for the median and the 97.5th
#' percentile of the shrink factor, with a dashed line at 1.1 indicating the
#' threshold for convergence. This function is not planned to be used in
#' isolation but rather within [IASDT.R::PlotGelman].
#' @param CodaObj An object of class `mcmc.list`, representing the Hmsc samples.
#' @name PlotGelman_Rho
#' @author Ahmed El-Gabbas
#' @return A ggplot object representing the Gelman-Rubin-Brooks plot for the
#'   `rho` parameter.
#' @export

PlotGelman_Rho <- function(CodaObj) {

  if (is.null(CodaObj)) {
    stop("CodaObj cannot be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- ShrinkFactor <- NULL

  Gelman_Rho_Plot <- CodaObj %>%
    gelman.preplot(
      bin.width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = TRUE) %>%
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
    ggplot2::scale_color_manual(
      values = c("Median" = "red", "Q97_5" = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(
        c(`Median` = "Median", `Q97_5` = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot - Rho",
      subtitle = NULL,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's ",
        "shrink factor as the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme(
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  return(Gelman_Rho_Plot)
}
