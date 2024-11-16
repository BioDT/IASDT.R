## |------------------------------------------------------------------------| #
# PlotGelman_Alpha ----
## |------------------------------------------------------------------------| #

#' Creates a Gelman-Rubin-Brooks plot for `alpha` parameters from an Hmsc
#' object.
#'
#' This function generates a Gelman-Rubin-Brooks plot for `alpha` parameters. It
#' uses parallel processing to speed up the computation. The plot includes lines
#' for the median and the 97.5th percentile of the shrink factor, with a dashed
#' line at 1.1 indicating the threshold for convergence. This function is not
#' planned to be used in isolation but rather within [IASDT.R::PlotGelman].
#' @param CodaObj An object of class `mcmc.list`, representing the Hmsc samples.
#' @param PlottingAlpha A double specifying the alpha (transparency) level for
#'   the plot lines. Default is 0.25.
#' @name PlotGelman_Alpha
#' @author Ahmed El-Gabbas
#' @return A ggplot object(s) representing the Gelman-Rubin-Brooks plot for the
#'   `beta` parameters.
#' @export

PlotGelman_Alpha <- function(CodaObj, PlottingAlpha = 0.25) {

  # # ..................................................................... ###

  if (!inherits(CodaObj, "mcmc.list")) {
    stop("CodaObj has to be of class mcmc.list", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var_LV <- Type <- Iter <- Plot <- data <- Median <- NULL

  # # ..................................................................... ###

  Gelman_Alpha_DT <- magrittr::extract2(CodaObj, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Alpha_Preplot1 <- lapply(CodaObj, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Alpha_Preplot2 <- try(
          gelman.preplot(
            x = Alpha_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Alpha_Preplot2, "try-error")) {
          Alpha_Preplot2 <- gelman.preplot(
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

  Xlim <- range(Gelman_Alpha_DT$Iter)

  Gelman_Alpha_Plot <- Gelman_Alpha_DT %>%
    dplyr::mutate(
      group = paste0(Var_LV, "_", Type),
      Var_LV = purrr::map_chr(
        .x = Var_LV, .f = ~stringr::str_remove_all(.x, "Alpha1\\[|\\]"))) %>%
    tidyr::nest(.by = "Var_LV") %>%
    dplyr::mutate(
      Plot = purrr::map2(
        .x = Var_LV, .y = data,
        .f = ~{
          .y %>%
            ggplot2::ggplot() +
            ggplot2::geom_line(
              mapping = ggplot2::aes(
                x = Iter, y = ShrinkFactor, group = group, color = Type),
              alpha = PlottingAlpha) +
            ggplot2::scale_color_manual(
              values = c("Median" = "red", "Q97_5" = "black")) +
            ggplot2::geom_hline(
              yintercept = 1.1, linetype = "dashed", col = "darkgrey",
              linewidth = 0.8) +
            ggplot2::facet_grid(~ Type, labeller = ggplot2::label_parsed) +
            ggplot2::scale_x_continuous(limits = Xlim, expand = c(0, 0)) +
            ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
            ggplot2::labs(
              title = "Gelman-Rubin-Brooks plot - Alpha",
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
              plot.title = ggplot2::element_text(
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
