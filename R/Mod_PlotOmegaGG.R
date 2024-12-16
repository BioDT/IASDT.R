## |------------------------------------------------------------------------| #
# PlotOmegaGG ----
## |------------------------------------------------------------------------| #

#' Creates heatmaps of parameter estimates or posterior support values for
#' species' residual association (`Omega` parameters).
#'
#' This function generates heatmaps to visualize the parameter estimates or
#' posterior support values of species' residual associations (`Omega`
#' parameters). It is designed to work with model output files and produces two
#' types of visualizations: one indicating the sign (positive or negative) of
#' the associations and another showing the mean values of these associations.
#' @param Path_Model String. Path to the fitted Hmsc model object.
#' @param supportLevel Numeric. The threshold for posterior support values used
#'   to determine which associations are strong enough to be plotted. Only
#'   associations with posterior support exceeding this threshold (or falling
#'   below 1 - threshold for negative associations) will be visualized. Defaults
#'   to 0.95.
#' @param PlotWidth,PlotHeight Integer. Specifies the width and height of the
#'   generated plot in centimeters. Defaults to `26` x `22.5`.
#' @return Generates two JPEG files containing the heatmaps of Omega parameter:
#'   signs and mean values. These files are saved in a directory named
#'   'Model_Postprocessing' within the parent directory of the provided model
#'   file path. The function itself returns `NULL` invisibly.
#' @name PlotOmegaGG
#' @export

PlotOmegaGG <- function(
    Path_Model, supportLevel = 0.95, PlotWidth = 26, PlotHeight = 22.5) {

  # # ..................................................................... ###

  # Set null device for `cairo`. This is to properly render the plots using
  # ggtext - https://github.com/wilkelab/cowplot/issues/73
  cowplot::set_null_device("cairo")

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be empty", call. = FALSE)
  }

  # # ..................................................................... ###

  # Out path -----

  Path_Out <- dirname(dirname(Path_Model)) %>%
    file.path("Model_Postprocessing", "Parameters_Summary")
  fs::dir_create(Path_Out)

  # # ..................................................................... ###

  # Loading model object ------
  IASDT.R::CatTime("Loading model object")

  Model <- IASDT.R::LoadAs(Path_Model)

  # # ..................................................................... ###

  # Plot phylogenetic tree -----
  IASDT.R::CatTime("Phylogenetic tree plot")

  Tree <- Model$phyloTree
  # Remove the 'Sp_' prefix from tip labels
  Tree$tip.label <- stringr::str_remove(Tree$tip.label, "^Sp_")
  if (length(Tree$edge.length) == 2 * nrow(Tree$edge)) {
    Tree$edge.length <- rep(1, length(Tree$edge.length) / 2)
  }

  PhyloPlot <- ggtree::ggtree(
    tr = Tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 1) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.2, 0, 0, 0), "lines"))

  # # ..................................................................... ###

  # Plotting theme ----
  IASDT.R::CatTime("Plotting theme")

  Theme <- ggplot2::theme(
    legend.title = ggtext::element_markdown(),
    legend.spacing = ggplot2::unit(0, "cm"),
    legend.key.size = ggplot2::unit(0.65, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.box.margin = ggplot2::margin(0, -30, 0, -15),
    legend.box.spacing = ggplot2::unit(0, "pt"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, -3, 0.25, -2), "lines"))

  # # ..................................................................... ###

  # ComputeAssociations -----
  IASDT.R::CatTime("Compute associations")

  post <- Hmsc::computeAssociations(Model)[[1]]

  # # ..................................................................... ###

  # Sign ------
  IASDT.R::CatTime("1. sign")

  Support <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0)
  # remove prefix "Sp_" from co-occurrence labels
  dimnames(Support)[[2]] <- dimnames(Support)[[1]] <- stringr::str_remove(
    dimnames(Support)[[1]], "^Sp_")

  PostMean <- post$mean
  PostMean[!Support] <- NA_real_
  # remove prefix "Sp_" from co-occurrence labels
  dimnames(PostMean)[[2]] <- dimnames(PostMean)[[1]] <- stringr::str_remove(
    dimnames(PostMean)[[1]], "^Sp_")

  PosSign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  NegSign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Omega</b></span><br>',
    '<span style="font-size: 9pt">(sign)</span>')

  Plot_Sign <- (
    sign(PostMean) %>%
      as.data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      # replace diagonals with NA
      replace(., col(.) == row(.), NA_character_) %>%
      replace(., . == "1", PosSign) %>%
      replace(., . == "-1", NegSign) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.75, width = 12, font.size = 0.75,
        colnames_offset_y = -1, colnames_angle = 90, hjust = 0.5) +
      ggplot2::scale_fill_manual(
        values = c("red", "blue"), na.value = "transparent",
        breaks = c(PosSign, NegSign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 8))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()


  Plot <- cowplot::plot_grid(
    (Plot_Sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Sign)),
    rel_widths = c(1, 0.09))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Omega_Sign.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  IASDT.R::CatTime("2. mean")

  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Omega</span><br>',
    '<span style="font-size: 9pt">(mean)</span>')

  Plot_Mean <- (
    PostMean %>%
      as.data.frame() %>%
      # replace diagonals with NA
      replace(., col(.) == row(.), NA_real_) %>%
      ggtree::gheatmap(
        PhyloPlot, ., offset = 0.75, width = 12, font.size = 0.75,
        colnames_offset_y = -1, colnames_angle = 90, hjust = 1) +
      ggplot2::scale_fill_gradientn(
        na.value = "transparent", colours = colorRamps::matlab.like(200),
        labels = scales::number_format(accuracy = 0.1)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 8))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(1, 0.09))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Omega_Mean.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime, Prefix = "Plotting took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
