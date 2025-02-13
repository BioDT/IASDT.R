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
#' @param Path_Model String. Path to the fitted Hmsc model object.
#' @param supportLevel Numeric. The threshold for posterior support used in
#'   plotting. Values above this threshold (and below 1 - threshold) are
#'   considered significant and will be plotted. The default value is 0.95,
#'   indicating 95% posterior support. For more information, see
#'   [Hmsc::plotBeta]
#' @param PlotWidth,PlotHeight Numeric. The width and height of the plot in
#'   centimeters. Default is `25` cm x `35` cm.
#' @return The function does not return a value but saves heatmap plots as JPEG
#'   files in a directory related to the model's path.
#' @author Ahmed El-Gabbas
#' @name PlotBetaGG
#' @export

PlotBetaGG <- function(
    Path_Model = NULL, supportLevel = 0.95, PlotWidth = 25, PlotHeight = 35) {

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

  # Loading model object -----
  IASDT.R::CatTime("Loading model object")

  if (file.exists(Path_Model)) {
    Model <- IASDT.R::LoadAs(Path_Model)
  } else {
    stop("Model file not found", call. = FALSE)
  }

  # # ..................................................................... ###

  # Plot phylogenetic tree -----
  IASDT.R::CatTime("phylogenetic tree plot")

  Tree <- Model$phyloTree
  # remove prefix "Sp_" from tip labels
  Tree$tip.label <- stringr::str_remove(Tree$tip.label, "^Sp_")
  if (length(Tree$edge.length) == 2 * nrow(Tree$edge)) {
    Tree$edge.length <- rep(1, length(Tree$edge.length) / 2)
  }

  PhyloPlot <- ggtree::ggtree(
    tr = Tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 1) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.3, 0, 0, 0), "lines"))

  # # ..................................................................... ###

  # Plotting theme ----
  IASDT.R::CatTime("Plotting theme")

  Theme <- ggplot2::theme(
    legend.title = ggtext::element_markdown(),
    legend.spacing = ggplot2::unit(0, "cm"),
    legend.key.size = ggplot2::unit(0.65, "cm"),
    legend.key.width = ggplot2::unit(0.65, "cm"),
    legend.box.margin = ggplot2::margin(0, -20, 0, -30),
    legend.box.spacing = ggplot2::unit(0, "pt"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, -1, 1.75, -2), "lines"))

  # # ..................................................................... ###

  # getPostEstimate -----
  IASDT.R::CatTime("getPostEstimate")

  # Calculates mean, support and other posterior quantities for a specified
  # model parameter
  post <- Hmsc::getPostEstimate(hM = Model, parName = "Beta")

  # # ..................................................................... ###

  CovNames <- Model$covNames %>%
    stringr::str_remove("stats::poly\\(") %>%
    stringr::str_replace_all(", degree = 2, raw = TRUE\\)", "_") %>%
    stringr::str_replace_all("_1", "\n(L)") %>%
    stringr::str_replace_all("_2", "\n(Q)")
  ColNames <- dplyr::case_when(
    CovNames == "(Intercept)" ~ "\n\nIntercept\n",
    CovNames == "EffortsLog" ~ "\n\nSampling\nefforts",
    CovNames == "RoadRailLog" ~ "\n\n\nRoad &\nRail\nintensity",
    CovNames == "RiversLog" ~ "\n\nRiver\nlength",
    CovNames == "HabLog" ~ "\n\nHabitat\ncoverage",
    .default = paste0("\n\n", CovNames))
  SupportMatrix <- (post$support > supportLevel) %>%
    magrittr::add(post$support < (1 - supportLevel)) %>%
    magrittr::is_greater_than(0)

  # Legend colours
  Palette_Positive <- grDevices::colorRampPalette(
    c("#FFE4B2", "#FFC85C", "#FF8C00", "#E25822", "#800000"))(100)
  Palette_Negative <- grDevices::colorRampPalette(
    c("#E6FFFF", "#91D8F7", "#4682B4", "#083D77", "#001F3F"))(100) %>%
    rev()

  rm(Model, envir = environment())

  # # ..................................................................... ###

  # Support ------
  IASDT.R::CatTime("1. support")

  Plot_SupportD <- (SupportMatrix * (2 * post$support - 1)) %>%
    t() %>%
    as.data.frame() %>%
    replace(., . == 0, NA_real_)
  # remove prefix "Sp_" from species labels
  dimnames(Plot_SupportD)[[1]] <- stringr::str_remove(
    dimnames(Plot_SupportD)[[1]], "^Sp_")
  colnames(Plot_SupportD) <- ColNames

  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 7pt">(support)</span>')

  Plot_SupportD_Positive <- Plot_SupportD_Negative <- Plot_SupportD
  Plot_SupportD_Positive[Plot_SupportD_Positive < 0] <- NA_real_
  Plot_SupportD_Negative[Plot_SupportD_Negative >= 0] <- NA_real_

  suppressMessages(
    {
      Plot_Support1 <- ggtree::gheatmap(
        PhyloPlot, Plot_SupportD_Negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = LegendTitle) +
        ggnewscale::new_scale_fill()

      Plot_Support <- ggtree::gheatmap(
        Plot_Support1, Plot_SupportD_Positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggtree::scale_x_ggtree() +
        ggplot2::coord_cartesian(clip = "off")  +
        Theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  Plot <- cowplot::plot_grid(
    (Plot_Support + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Support)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Support.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Sign ------

  IASDT.R::CatTime("1. sign")

  PosSign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  NegSign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</b></span>',
    '<br><span style="font-size: 9pt">(sign)</span>')

  Plot_SignD <- (SupportMatrix * sign(post$mean)) %>%
    sign(x = .) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate_all(as.character) %>%
    replace(., . == "1", PosSign) %>%
    replace(., . == "-1", NegSign) %>%
    replace(., . == "0", NA_character_)

  # remove prefix "Sp_" from species labels
  dimnames(Plot_SignD)[[1]] <- stringr::str_remove(
    dimnames(Plot_SignD)[[1]], "^Sp_")
  colnames(Plot_SignD) <- ColNames

  Plot_Sign <- (
    ggtree::gheatmap(
      PhyloPlot, Plot_SignD, offset = -0.85, width = 12,
      font.size = 2.5, hjust = 0.5) +
      ggplot2::scale_fill_manual(
        values = c("#E25822", "#083D77"), na.value = "transparent",
        breaks = c(PosSign, NegSign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = LegendTitle) +
      Theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 6))) %>%
    suppressMessages()

  Plot <- cowplot::plot_grid(
    (Plot_Sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Sign)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Sign.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  IASDT.R::CatTime("2. mean")

  Plot_MeanD <- (SupportMatrix * post$mean) %>%
    t() %>%
    as.data.frame() %>%
    replace(., . == 0, NA_real_)
  # remove prefix "Sp_" from species labels
  dimnames(Plot_MeanD)[[1]] <- stringr::str_remove(
    dimnames(Plot_MeanD)[[1]], "^Sp_")
  colnames(Plot_MeanD) <- ColNames

  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 9pt">(mean)</span>')

  Plot_MeanD_Positive <- Plot_MeanD_Negative <- Plot_MeanD
  Plot_MeanD_Positive[Plot_MeanD_Positive < 0] <- NA_real_
  Plot_MeanD_Negative[Plot_MeanD_Negative >= 0] <- NA_real_

  suppressMessages(
    {
      Plot_Mean1 <- ggtree::gheatmap(
        PhyloPlot, Plot_MeanD_Negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggtree::scale_x_ggtree() +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = LegendTitle) +
        ggnewscale::new_scale_fill()

      Plot_Mean <- ggtree::gheatmap(
        Plot_Mean1, Plot_MeanD_Positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggplot2::coord_cartesian(clip = "off")  +
        Theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Mean1.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean - without intercept -----
  IASDT.R::CatTime("3. Mean - without intercept")

  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 9pt">(mean)</span><br>',
    '<span style="font-size: 7pt">[excl. Intercept]</span>')

  Plot_MeanD_Positive <- Plot_MeanD_Positive[, -1]
  Plot_MeanD_Negative <- Plot_MeanD_Negative[, -1]

  suppressMessages(
    {
      Plot_Mean1 <- ggtree::gheatmap(
        PhyloPlot, Plot_MeanD_Negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggtree::scale_x_ggtree() +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = LegendTitle) +
        ggnewscale::new_scale_fill()

      Plot_Mean <- ggtree::gheatmap(
        Plot_Mean1, Plot_MeanD_Positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggplot2::coord_cartesian(clip = "off")  +
        Theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = file.path(Path_Out, "Parameter_Beta_Mean2.jpeg"), res = 600,
    width = PlotWidth, height = PlotHeight, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime, Prefix = "Plotting took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
