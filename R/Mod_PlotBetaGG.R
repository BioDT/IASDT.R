## |------------------------------------------------------------------------| #
# PlotBetaGG ----
## |------------------------------------------------------------------------| #

#' Heatmaps of parameter estimates or posterior support values of species'
#' environmental responses (Beta parameters)
#'
#' This function generates heatmaps of parameter estimates or posterior support
#' values for species' environmental responses, represented by Beta parameters.
#' It is designed to visualize how species (Y) respond to various covariates (X)
#' using `ggplot2` for plotting. The function is an adaptation of
#' [Hmsc::plotBeta], focusing on `ggplot2`-based visualizations.
#' @param Path_Model String. The path to the `.RData` file containing the
#'   selected model.
#' @param supportLevel Numeric. The threshold for posterior support used in
#'   plotting. Values above this threshold (and below 1 - threshold) are
#'   considered significant and will be plotted. The default value is 0.95,
#'   indicating 95% posterior support. For more information, see
#'   [Hmsc::plotBeta]
#' @param PlotWidth,PlotHeight Numeric. The width and height of the plot in
#'   centimeters. Default is `18` cm x `20` cm.
#' @return The function does not return a value but saves heatmap plots as JPEG
#'   files in a directory related to the model's path.
#' @name PlotBetaGG
#' @export

PlotBetaGG <- function(
    Path_Model = NULL, supportLevel = 0.95, PlotWidth = 18, PlotHeight = 20) {

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be empty", call. = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Out path
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Path_Out <- dirname(dirname(Path_Model)) %>%
    file.path("Model_Postprocessing", "Parameters_Summary")
  fs::dir_create(Path_Out)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading model object
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading model object")
  Model <- IASDT.R::LoadAs(Path_Model)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plot phylogenetic tree -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("phylogenetic tree plot")

  Tree <- Model$phyloTree

  if (length(Tree$edge.length) == 2 * nrow(Tree$edge)) {
    Tree$edge.length <- rep(1, length(Tree$edge.length) / 2)
  }

  PhyloPlot <- ggtree::ggtree(
    tr = Tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 2) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.2, 0, 0, 0), "lines"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plotting theme ----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plotting theme")

  Theme <- ggplot2::theme(
    legend.title = ggtext::element_markdown(),
    legend.spacing = ggplot2::unit(0, "cm"),
    legend.key.size = ggplot2::unit(0.5, "cm"),
    legend.key.width = ggplot2::unit(0.4, "cm"),
    legend.box.margin = ggplot2::margin(0, -7, 0, 1),
    legend.box.spacing = ggplot2::unit(0, "pt"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(-0.25, -2.25, 1.25, -1), "lines"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # getPostEstimate -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("getPostEstimate")
  # Calculates mean, support and other posterior quantities for a specified
  # model parameter
  post <- Hmsc::getPostEstimate(hM = Model, parName = "Beta")

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Sign ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("1. sign")
  Plot_SignD <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0) %>%
    magrittr::multiply_by(sign(post$mean))

  rownames(Plot_SignD) <- dplyr::case_when(
    Model$covNames == "(Intercept)" ~ "\n\nIntercept\n",
    Model$covNames == "BiasLog" ~ "\n\nSampling\nintensity",
    Model$covNames == "RoadRailLog" ~ "\n\nRoad&Rail\nintensity",
    .default = paste0("\n\n", Model$covNames, "\n"))

  PosSign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  NegSign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  LegendTitle <- paste0('<span style="font-size: 12pt"><b>Beta</b></span>",
                        "<br><span style="font-size: 9pt">(sign)</span>')

  Plot_Sign <- (
    Plot_SignD %>%
      sign(x = .) %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      replace(., . == "1", PosSign) %>%
      replace(., . == "-1", NegSign) %>%
      replace(., . == "0", NA_character_) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.75, width = 12, font.size = 3, hjust = 0.5) +
      ggplot2::scale_fill_manual(
        values = c("red", "blue"), na.value = "transparent",
        breaks = c(PosSign, NegSign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 6))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Sign)),
    rel_widths = c(0.91, 0.09))
  
  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Sign.jpeg"),
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100, 
    res = 600)
  Plot
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Mean -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("2. mean")

  Plot_MeanD <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0) %>%
    magrittr::multiply_by(post$mean)

  rownames(Plot_MeanD) <- dplyr::case_when(
    Model$covNames == "(Intercept)" ~ "\n\nIntercept\n",
    Model$covNames == "BiasLog" ~ "\n\nSampling\nintensity",
    Model$covNames == "RoadRailLog" ~ "\n\nRoad&Rail\nintensity",
    .default = paste0("\n\n", Model$covNames, "\n"))

  LegendTitle <- paste0('<span style="font-size: 12pt"><b>Beta</span><br>",
                        "<span style="font-size: 9pt">(mean)</span>')

  Plot_Mean <- (
    Plot_MeanD %>%
      t() %>%
      as.data.frame() %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.75, width = 12, font.size = 3, hjust = 0.5) +
      ggplot2::scale_fill_gradientn(
        na.value = "transparent", colours = colorRamps::matlab.like(200)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 6))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(0.9, 0.1))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Mean1.jpeg"),
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100, 
    res = 600)
  Plot
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Mean - without intercept -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("3. Mean - without intercept")

  Plot_MeanD <- Plot_MeanD[-1, ]

  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 9pt">(mean)</span><br><br>',
    '<span style="font-size: 7pt">[excl.<br/>Intercept]</span>')

  Plot_Mean <- (
    Plot_MeanD %>%
      t() %>%
      as.data.frame() %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.5, width = 12, font.size = 3, hjust = 0.5) +
      ggplot2::scale_fill_gradientn(
        na.value = "transparent", colours = colorRamps::matlab.like(200)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 6))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(0.9, 0.1))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Mean2.jpeg"),
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100, 
    res = 600)
  Plot
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Support ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("4. support")

  Plot_SupportD <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0) %>%
    magrittr::multiply_by((2 * post$support - 1))

  rownames(Plot_SupportD) <- dplyr::case_when(
    Model$covNames == "(Intercept)" ~ "\n\nIntercept\n",
    Model$covNames == "BiasLog" ~ "\n\nSampling\nintensity",
    Model$covNames == "RoadRailLog" ~ "\n\nRoad&Rail\nintensity",
    .default = paste0("\n\n", Model$covNames, "\n"))

  LegendTitle <- paste0('<span style="font-size: 12pt"><b>Beta</span><br>',
                        '<span style="font-size: 7pt">(support)</span>')

  Plot_Support <- (
    Plot_SupportD %>%
      t() %>%
      # abs() %>%
      as.data.frame() %>%
      replace(., . == 0, NA_real_) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.75, width = 12, font.size = 3, hjust = 0.5) +
      ggplot2::scale_fill_gradientn(
        na.value = "transparent", colours = colorRamps::matlab.like(200)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 6))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Support + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Support)),
    rel_widths = c(0.9, 0.1))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Support.jpeg"),
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100, 
    res = 600)
  Plot
  grDevices::dev.off()

  return(invisible(NULL))
}
