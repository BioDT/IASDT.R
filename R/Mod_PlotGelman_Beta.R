## |------------------------------------------------------------------------| #
# PlotGelman_Beta ----
## |------------------------------------------------------------------------| #

#' Gelman-Rubin-Brooks plot for `beta` parameters
#'
#' This function generates a Gelman-Rubin-Brooks plot for `beta` parameters from a Hmsc object. It is used to assess the convergence of MCMC simulations by plotting the evolution of the shrink factor over iterations. The plot helps in identifying whether the chains have converged to a common distribution. The plot includes lines for the median and the 97.5th percentile of the shrink factor, with a dashed line at 1.1 indicating the threshold for convergence.
#'
#' @param CodaObj An object of class `mcmc.list`, representing the Hmsc samples.
#' @param NCores Integer indicating the number of parallel processes to be used for computation.
#' @param EnvFile String specifying the path to read the environment variables. Default value is `.env`.
#' @param PlottingAlpha Double indicating the alpha (transparency) level for the plotting lines. Default value is 0.25.
#' @name PlotGelman_Beta
#' @author Ahmed El-Gabbas
#' @return A ggplot object representing the Gelman-Rubin-Brooks plot for the `beta` parameters.
#' @export

PlotGelman_Beta <- function(
    CodaObj, NCores, EnvFile = ".env", PlottingAlpha = 0.25) {


  if (is.null(CodaObj) || is.null(NCores)) {
    stop("CodaObj and NCores cannot be empty")
  }

  if (!is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.")
  }

  if (magrittr::not(inherits(CodaObj, "mcmc.list"))) {
    stop("CodaObj has to be of class mcmc.list")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Var_Sp <- dplyr <- coda <- tibble <- magrittr <-
    ShrinkFactor <- group <- NULL

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(snow::stopCluster(c1), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  Beta_Coda <- IASDT.R::Coda_to_tibble(
    CodaObj = CodaObj, Type = "beta", EnvFile = EnvFile)
  NVars <- length(unique(Beta_Coda$Variable))
  NSp <- length(unique(Beta_Coda$Species))
  SubTitle <- paste0(NVars, " covariates - ", NSp, " species")
  rm(Beta_Coda)

  BetaNames <- CodaObj %>%
    magrittr::extract2(1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort()

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, tibble, magrittr)))
  snow::clusterExport(cl = c1, list = "CodaObj", envir = environment())

  Gelman_Beta_Vals <- snow::parLapply(
    cl = c1, x = BetaNames,
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
          Type = factor(Type), Var_Sp = x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(group = paste0(Var_Sp, "_", Type))

  Gelman_Beta_Plot <- Gelman_Beta_Vals %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = PlottingAlpha) +
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
      title = "Gelman-Rubin-Brooks plot - Beta", subtitle = SubTitle,
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

  return(Gelman_Beta_Plot)
}
