## |------------------------------------------------------------------------| #
# PlotGelman_Alpha ----
## |------------------------------------------------------------------------| #

#' Gelman-Rubin-Brooks plot `alpha` parameters
#'
#' Gelman-Rubin-Brooks plot `alpha` parameters
#'
#' @param CodaObj an mcmc object
#' @param NCores Integer. Number of parallel processes.
#' @param PlottingAlpha Double. Plotting alpha for line transparency
#' @name PlotGelman_Alpha
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export


PlotGelman_Alpha <- function(
    CodaObj = NULL, NCores = NULL, PlottingAlpha = 0.25) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Plot <- Var_LV <- Type <- Iter <- dplyr <- coda <- tibble <-
    magrittr <- data <- NULL

  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  on.exit(snow::stopCluster(c1), add = TRUE)

  AlphaNames <- CodaObj %>%
    magrittr::extract2(1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort()

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, tibble, magrittr)))
  snow::clusterExport(cl = c1, list = "CodaObj", envir = environment())

  Gelman_Alpha_Plot <- snow::parLapply(
    cl = c1, x = AlphaNames,
    fun = function(x) {
      lapply(CodaObj, function(Y) {
        Y[, x, drop = TRUE]
      }) %>%
        coda::mcmc.list() %>%
        coda:::gelman.preplot(
          bin.width = 10, max.bins = 50, confidence = 0.95,
          transform = FALSE, autoburnin = TRUE) %>%
        magrittr::extract2("shrink") %>%
        tibble::as_tibble(rownames = "Iter") %>%
        purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
        dplyr::mutate(Iter = as.integer(Iter)) %>%
        tidyr::pivot_longer(
          cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
        dplyr::arrange(Type, Iter) %>%
        dplyr::mutate(Type = factor(Type), Var_LV = x)
    }) %>%
    dplyr::bind_rows() %>%
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
            ggplot2::coord_cartesian(expand = FALSE) +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = "Gelman-Rubin-Brooks plot - Alpha",
              subtitle = .x,
              caption = "This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases.") +
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
