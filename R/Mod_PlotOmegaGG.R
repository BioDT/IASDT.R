## |------------------------------------------------------------------------| #
# PlotOmegaGG ----
## |------------------------------------------------------------------------| #

#' Creates heatmaps of parameter estimates or posterior support values for species' residual association (`Omega` parameters).
#'
#' This function generates heatmaps to visualize the parameter estimates or posterior support values of species' residual associations (`Omega` parameters). It is designed to work with model output files and produces two types of visualizations: one indicating the sign (positive or negative) of the associations and another showing the mean values of these associations.
#'
#' @param Path_Model String. Specifies the path to the `.RData` file containing the selected model.
#' @param supportLevel Numeric. The threshold for posterior support values used to determine which associations are strong enough to be plotted. Only associations with posterior support exceeding this threshold (or falling below 1 - threshold for negative associations) will be visualized. Defaults to 0.95.
#' @param PlotWidth,PlotHeight Integer. Specifies the width and height of the generated plot in centimeters. Defaults to `22` x `20`.
#' @return Generates two JPEG files containing the heatmaps of Omega parameter: signs and mean values. These files are saved in a directory named 'Model_Postprocessing' within the parent directory of the provided model file path. The function itself returns `NULL` invisibly.
#' @name PlotOmegaGG
#' @export

PlotOmegaGG <- function(
    Path_Model, supportLevel = 0.95, PlotWidth = 22, PlotHeight = 20) {

  # devtools::install_github("YuLab-SMU/ggtree")

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be empty")
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Out path
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Path_Out <- Path_Model %>%
    dirname() %>%
    dirname() %>%
    file.path("Model_Postprocessing")
  fs::dir_create(Path_Out)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading model object
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading model object")
  Model <- IASDT.R::LoadAs(Path_Model)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plot phylogenetic tree -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Phylogenetic tree plot")

  Tree <- Model$phyloTree
  if (length(Tree$edge.length) == 2 * nrow(Tree$edge)) {
    Tree$edge.length <- rep(1, length(Tree$edge.length) / 2)
  }

  PhyloPlot <- ggtree::ggtree(
    tr = Tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 2) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.2, 0, 0, 0), "lines"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plotting theme ----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plotting theme")

  Theme <- ggplot2::theme(
    # legend.title = ggplot2::element_text(size = 7, face = "bold"),
    legend.title = ggtext::element_markdown(),
    legend.spacing = ggplot2::unit(0, "cm"),
    legend.key.size = ggplot2::unit(0.5, "cm"),
    legend.key.width = ggplot2::unit(0.4, "cm"),
    legend.box.margin = ggplot2::margin(0, -8, 0, 0),
    legend.box.spacing = ggplot2::unit(0, "pt"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(-2.1, -2.75, 1.25, -1.65), "lines"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # computeAssociations -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("ComputeAssociations")
  post <- Hmsc::computeAssociations(Model)[[1]]

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Sign ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("1. sign")

  Support <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0)
  PostMean <- post$mean
  PostMean[!Support] <- NA_real_

  PosSign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  NegSign <- '<span style="font-size: 8pt"><b>  \U2212  </b></span>'
  LegendTitle <- '<span style="font-size: 12pt"><b>Beta</b></span><br><span style="font-size: 9pt">(sign)</span>'

  Plot_Sign <- (
    sign(PostMean) %>%
      as.data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      # replace diagonals with NA
      replace(., col(.) == row(.), NA_character_) %>%
      replace(., . == "1", PosSign) %>%
      replace(., . == "-1", NegSign) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 2.75, width = 12, font.size = 2,
        colnames_angle = 90, hjust = 1) +
      ggplot2::scale_fill_manual(
        values = c("red", "blue"), na.value = "transparent",
        breaks = c(PosSign, NegSign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 8))) %>%
    # suppress the message: Scale for fill is already present. Adding another scale for fill, which will replace the existing scale.
    suppressMessages()

  cowplot::plot_grid(
    (Plot_Sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Sign)),
    rel_widths = c(0.91, 0.09)) %>%
    ggplot2::ggsave(
      filename = file.path(Path_Out, "Parameter_Omega_Sign.jpeg"),
      width = PlotWidth, height = PlotHeight, units = "cm", dpi = 600)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Mean -----
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("2. mean")

  LegendTitle <- '<span style="font-size: 12pt"><b>Beta</span><br><span style="font-size: 9pt">(mean)</span>'

  Plot_Mean <- (
    PostMean %>%
      as.data.frame() %>%
      # replace diagonals with NA
      replace(., col(.) == row(.), NA_real_) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 2.75, width = 12, font.size = 2,
        colnames_angle = 90, hjust = 1) +
      ggplot2::scale_fill_gradientn(
        na.value = "transparent", colours = colorRamps::matlab.like(200),
        labels = scales::number_format(accuracy = 0.1)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 8))) %>%
    # suppress the message: Scale for fill is already present. Adding another scale for fill, which will replace the existing scale.
    suppressMessages()

  cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(0.91, 0.09)) %>%
    ggplot2::ggsave(
      filename = file.path(Path_Out, "Parameter_Omega_Mean.jpeg"),
      width = PlotWidth, height = PlotHeight, units = "cm", dpi = 600)

  return(invisible(NULL))
}
