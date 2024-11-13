## |------------------------------------------------------------------------| #
# PlotAlpha ----
## |------------------------------------------------------------------------| #

#' Plot convergence traceplots for the alpha parameter
#'
#' This function generates and plots convergence traceplots for the alpha
#' parameter from an Hmsc model, using various diagnostics and visualizations to
#' assess convergence. It supports customization of plot aesthetics and layout.
#' @param Post A `coda` object or a string path to a saved `coda` object
#'   containing the posterior samples for the alpha parameter.
#' @param Model A fitted model object or a string path to a saved model object
#'   containing the samples and metadata necessary for plotting.
#' @param Title A string specifying the main title of the plot.
#' @param NRC  numeric vector of length 2 specifying the number of rows and
#'   columns per page of plots. If `NULL`, it is automatically determined based
#'   on the number of levels of the alpha parameter.
#' @param AddFooter A logical value indicating whether to add a footer with page
#'   numbers to each plot page. Defaults to `TRUE`.
#' @param AddTitle A logical value indicating whether to add the main title
#'   (specified by `Title`) to the plot. Defaults to `TRUE`.
#' @param Cols Character vector for chain colours (optional). Default: `NULL`.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param MarginType A string specifying the type of marginal plot to add to the
#'   main plot. Must be either "histogram" or "density". Default: "histogram".
#' @name PlotAlpha
#' @author Ahmed El-Gabbas
#' @return A grid of `ggplot2` objects representing the traceplots for the alpha
#'   parameter, arranged according to the specified number of rows and columns
#'   per page. If `AddTitle` or `AddFooter` is `TRUE`, additional text elements
#'   are included as specified.
#' @export

PlotAlpha <- function(
    Post = NULL, Model = NULL, Title = NULL, NRC = NULL, AddFooter = TRUE,
    AddTitle = TRUE, Cols = NULL, FromHPC = TRUE, MarginType = "histogram") {

  # # ..................................................................... ###

  if (is.null(Post) || is.null(Model)) {
    stop("Post and Model cannot be empty", call. = FALSE)
  }

  if (length(MarginType) != 1) {
    stop("`MarginType` must be a single string.")
  }

  if (!MarginType %in% c("histogram", "density")) {
    stop("`MarginType` must be either 'histogram' or 'density'.")
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Factor <- Value <- PointEst <- UpperCI <- NULL

  # # ..................................................................... ###

  # Checking arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = c("Title"))

  # # ..................................................................... ###

  # Load coda object
  if (inherits(Post, "character")) {
    Post <- IASDT.R::LoadAs(Post)
  }

  if ("Alpha" %in% names(Post)) {
    Post <- Post$Alpha[[1]]
  } else {
    stop("Post object does not contain 'Alpha'", call. = FALSE)
  }

  # Load model object
  if (inherits(Model, "character")) {
    Model <- IASDT.R::LoadAs(Model)
  }

  SampleSize <- Model$samples
  NChains <- length(Model$postList)

  #  Plotting colours
  if (is.null(Cols)) {
    Cols <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }
  if (length(Cols) != NChains) {
    warning(
      "The length of provided colours != number of chains", .call. = FALSE)
    Cols <- c(
      "black", "grey60",
      RColorBrewer::brewer.pal(n = NChains - 2, name = "Set1"))
  }

  rm(Model)

  # Number of latent factors
  NLV <- ncol(Post[[1]])

  # # ..................................................................... ###

  ## Gelman convergence diagnostic

  Gelman <- coda::gelman.diag(x = Post, multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    as.data.frame() %>%
    stats::setNames(c("PointEst", "UpperCI")) %>%
    round(2) %>%
    dplyr::mutate(Gelman = paste0(PointEst, " / ", UpperCI)) %>%
    dplyr::pull(Gelman)

  ## Effective sample size
  ESS <- coda::effectiveSize(x = Post) %>%
    magrittr::divide_by(NChains) %>%
    round(1)

  ## quantiles
  CI <- summary(Post, quantiles = c(0.025, 0.975)) %>%
    magrittr::extract2("quantiles") %>%
    matrix(ncol = 2) %>%
    magrittr::divide_by(1000) %>%
    round(3)

  AlphaDF <- IASDT.R::Coda_to_tibble(
    CodaObj = Post, Type = "alpha", FromHPC = FromHPC) %>%
    dplyr::mutate(
      Factor2 = purrr::map_int(
        .x = Factor,
        .f = ~ {
          as.character(.x) %>%
            stringr::str_remove("factor") %>%
            as.integer()
        })) %>%
    dplyr::mutate(Value = Value / 1000)

  if (is.null(NRC)) {
    NRC <- dplyr::case_when(
      NLV == 1 ~ c(1, 1), NLV == 2 ~ c(1, 2),
      NLV == 3 ~ c(1, 3), NLV == 4 ~ c(2, 2),
      .default = c(2, 3))
  }

  # # ..................................................................... ###

  Plots <- purrr::map(
    .x = seq_len(NLV),
    .f = ~ {
      ESS0 <- paste0(
        "<b><i>Mean effective sample size:</i></b> ", ESS[.x],
        " / ", SampleSize, " samples")

      CI0 <- paste0(round(CI[.x, ], 2), collapse = " to ") %>%
        paste0("<b><i>95% credible interval:</i></b> ", ., " km")

      ESS_CI <- data.frame(
        x = -Inf, y = -Inf, label = paste0(ESS0, "<br>", CI0))

      Gelman0 <- paste0(
        "<b><i>Gelman convergence diagnostic:</i></b> ", Gelman[.x])

      Title2 <- data.frame(x = Inf, y = Inf, label = Gelman0)
      Title3 <- data.frame(
        x = -Inf, y = Inf, label = paste0("<b>Factor", .x, "</b>"))

      PlotDT <- dplyr::filter(AlphaDF, Factor2 == .x)

      # Pre-calculate smoothed lines
      summary_data <- PlotDT %>%
        group_by(Chain) %>%
        mutate(Smoothed = pmax(predict(loess(Value ~ Iter)), 0))

      Plot <- ggplot2::ggplot(
        data = PlotDT,
        mapping = ggplot2::aes(
          x = Iter, y = Value, color = Chain, group = Chain)) +
        ggplot2::geom_line(linewidth = 0.125, alpha = 0.6) +
        ggplot2::geom_line(
          data = summary_data, linewidth = 0.8, alpha = 0.6,
          mapping = ggplot2::aes(x = Iter, y = Smoothed, color = Chain)) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI[.x, ], linetype = "dashed", color = "black",
          linewidth = 1) +
        ggplot2::scale_color_manual(values = Cols) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
          limits = c(0, max(AlphaDF$Value) * 1.05),
          expand = c(0, 0), oob = scales::squish) +
        ggplot2::theme_bw() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Distance (km)") +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), data = Title2,
          inherit.aes = FALSE, size = 4, hjust = 1, vjust = 1, lineheight = 0,
          fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), data = Title3,
          inherit.aes = FALSE, size = 5, hjust = -0.1, vjust = 1,
          color = "blue", lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = ESS_CI, inherit.aes = FALSE, size = 4,
          hjust = 0, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 12))


      if (MarginType == "histogram") {
        Plot <- ggExtra::ggMarginal(
          p = Plot, type = MarginType, margins = "y", size = 6,
          color = "steelblue4", fill = "steelblue4", bins = 100)
      } else {
        Plot <- ggExtra::ggMarginal(
          p = Plot, type = MarginType, margins = "y", size = 6,
          color = "steelblue4")
      }

      # Making marginal background matching the plot background
      # https://stackoverflow.com/a/78196022/3652584
      Plot$layout$t[1] <- 1
      Plot$layout$r[1] <- max(Plot$layout$r)

      return(Plot)
    })


  if (AddTitle && AddFooter) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = Title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = NRC[1], ncol = NRC[2])
  }

  if (AddTitle && isFALSE(AddFooter)) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = NULL,
      top = grid::textGrob(
        label = Title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = NRC[1], ncol = NRC[2])
  }

  if (isFALSE(AddTitle) && AddFooter) {
    Plots <- gridExtra::marrangeGrob(
      Plots, bottom = bquote(paste0("page ", g, " of ", npages)),
      top = NULL, nrow = NRC[1], ncol = NRC[2])
  }

  if (isFALSE(AddTitle) && isFALSE(AddFooter)) {
    Plots <- cowplot::plot_grid(
      plotlist = Plots, ncol = NRC[2], nrow = NRC[1], align = "hv")
  }

  return(Plots)
}
