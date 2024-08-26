## |------------------------------------------------------------------------| #
# PlotGelman_Omega ----
## |------------------------------------------------------------------------| #

#' Creates a Gelman-Rubin-Brooks plot for `omega` parameters.
#'
#' This function generates a Gelman-Rubin-Brooks plot specifically for `omega`
#' parameters using the provided Hmsc model object. It is designed to help
#' assess the convergence of MCMC simulations by plotting the shrink factor over
#' iterations for a subset of species' omega parameters. This function is not
#' planned to be used in isolation but rather within [IASDT.R::PlotGelman].
#' @param CodaObj n object of class `mcmc.list` representing the MCMC chains.
#' @param NCores An integer specifying the number of cores to use for parallel
#'   processing.
#' @param NOmega An optional integer indicating the number of species' omega
#'   parameters to sample and plot. Defaults to 1000.
#' @param PlottingAlpha A numeric value between 0 and 1 indicating the
#'   transparency level of the plot lines. Defaults to 0.25.
#' @name PlotGelman_Omega
#' @author Ahmed El-Gabbas
#' @return A ggplot object representing the Gelman-Rubin-Brooks plot for the
#'   sampled omega parameters.
#' @export

PlotGelman_Omega <- function(
    CodaObj, NCores, NOmega = 1000, PlottingAlpha = 0.25) {

  if (is.null(CodaObj) || is.null(NCores)) {
    stop("CodaObj and NCores cannot be empty", call. = FALSE)
  }

  if (!is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }
 if (!is.numeric(NOmega) || NOmega <= 0) {
    stop("NOmega must be a positive integer.", call. = FALSE)
  }

  if (!inherits(CodaObj, "mcmc.list")) {
    stop("CodaObj has to be of class mcmc.list", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Sp_comb <- ShrinkFactor <- group <- NULL

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  OmegaNames <- magrittr::extract2(CodaObj, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sample(min(NOmega, length(.))) %>%
    sort()

  invisible(snow::clusterEvalQ(
    cl = c1,
    IASDT.R::LoadPackages(List = c("dplyr", "coda", "tibble", "magrittr"))))
  snow::clusterExport(cl = c1, list = "CodaObj", envir = environment())

  Gelman_OmegaDT <- snow::parLapply(
    cl = c1, x = OmegaNames,
    fun = function(x) {
      lapply(CodaObj, function(Y) {
        Y[, x, drop = TRUE]
      }) %>%
        coda::mcmc.list() %>%
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
        dplyr::mutate(Type = factor(Type), Sp_comb = x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(group = paste0(Sp_comb, "_", Type))

  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)

  Gelman_Omega_Plot <- Gelman_OmegaDT %>%
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
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's shrink factor as ",
        "the number of iterations increases.")) +
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

  return(Gelman_Omega_Plot)
}
