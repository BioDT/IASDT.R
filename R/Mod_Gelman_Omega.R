## |------------------------------------------------------------------------| #
# Gelman_Omega ----
## |------------------------------------------------------------------------| #

#' Gelman-Rubin-Brooks plot for all omega parameters
#'
#' Gelman-Rubin-Brooks plot for all omega parameters
#'
#' @param CodaObj an mcmc object
#' @param NCores Integer. Number of parallel processes.
#' @param NOmega Integer. Number of species to be sampled for the Omega parameter
#' @param PlotAlpha Double. Plotting alpha for line transparency
#' @name Gelman_Omega
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Gelman_Omega <- function(
    CodaObj = NULL, NCores = NULL, NOmega = 1000, PlotAlpha = 0.25) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Sp_comb <- dplyr <- coda <- tibble <- magrittr <-
    ShrinkFactor <- group <- NULL

  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  OmegaNames <- CodaObj %>%
    magrittr::extract2(1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sample(NOmega) %>%
    sort()

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, tibble, magrittr)))
  snow::clusterExport(
    cl = c1, list = c("CodaObj"), envir = environment())

  Gelman_Omega <- snow::parLapply(
    cl = c1, x = OmegaNames,
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
        dplyr::mutate(
          # Colour = dplyr::if_else(Type == "Median", "black", "red"),
          # LType = dplyr::if_else(Type == "Median", "dashed", "solid"),
          Type = factor(Type), Sp_comb = x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(group = paste0(Sp_comb, "_", Type))

  Gelman_Omega_Plot <- Gelman_Omega %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = PlotAlpha) +
    ggplot2::scale_color_manual(
      values = c("Median" = "red", "Q97_5" = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey", linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(
        c(`Median` = "Median", `Q97_5` = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Gelman-Rubin-Brooks plot - Omega - ",
                     NOmega, " species combination samples"),
      subtitle = NULL,
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

  snow::stopCluster(c1)
  return(Gelman_Omega_Plot)
}
