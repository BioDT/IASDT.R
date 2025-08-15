## |------------------------------------------------------------------------| #
# mod_heatmap_beta ----
## |------------------------------------------------------------------------| #

#' Heatmaps for the `beta` and `omega` parameters of the Hmsc model
#'
#' The `mod_heatmap_beta()` and `mod_heatmap_omega()` functions generate
#' heatmaps using `ggplot2` to visualise parameter estimates or posterior
#' support values for species' environmental responses (`beta` parameters, which
#' describes how species (*Y*) respond to various covariates (*X*); see
#' [Hmsc::plotBeta]) and residual associations (`omega` parameter),
#' respectively.
#' @param path_model Character. Path to the fitted Hmsc model object.
#' @param support_level Numeric. The posterior support threshold for determining
#'   which values are considered significant in the heatmap. Defaults to 0.95,
#'   indicating 95% posterior support. Values above this threshold (or below 1 -
#'   threshold for negative associations) are considered significant and will be
#'   plotted (see [Hmsc::plotBeta]).
#' @param width,height Integer. The width and height of the generated
#'   heatmaps in centimetres. Defaults to 26&times;22.5 for `omega`; 25&times;35
#'   for `beta`.
#' @return Both functions do not return a value but saves heatmap plots as JPEG
#'   files in the `Model_Postprocessing/Parameters_Summary` subdirectory.
#' @details The functions exports three types of visualisations (see
#'   [Hmsc::plotBeta]):
#' - `Mean`: posterior mean estimate,
#' - `Support`: statistical support level, measured by the posterior
#'   probability for a positive or negative response,
#' - `Sign`: indicates whether the response is positive, negative, or neither
#'   of these based on the chosen `support_level`.
#'
#' For the `omega` parameter, the `mod_heatmap_omega()` function generates two
#' JPEG files: signs and mean values. While for the `beta` parameter, the
#' `mod_heatmap_beta()` function generates four JPEG files : support, signs,
#' mean values (including and excluding the intercept).
#' @export
#' @name parameter_heatmap
#' @rdname parameter_heatmap
#' @order 1
#' @author Ahmed El-Gabbas. The `mod_heatmap_beta()` function is adapted from
#'   [Hmsc::plotBeta]

mod_heatmap_beta <- function(
    path_model = NULL, support_level = 0.95, width = 25, height = 35) {

  # # ..................................................................... ###

  # Set null device for `cairo`. This is to properly render the plots using
  # ggtext - https://github.com/wilkelab/cowplot/issues/73
  cowplot::set_null_device("cairo")

  .start_time <- lubridate::now(tzone = "CET")

  if (is.null(path_model)) {
    ecokit::stop_ctx(
      "`path_model` cannot be empty", path_model = path_model,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # support_level has to be between 0 and 1
  if (support_level < 0 || support_level > 1) {
    ecokit::stop_ctx(
      "`support_level` has to be a numeric value between 0 and 1",
      support_level = support_level, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Out path -----

  Path_Out <- dirname(dirname(path_model)) %>%
    fs::path("Model_Postprocessing", "Parameters_Summary")
  fs::dir_create(Path_Out)

  # # ..................................................................... ###

  # Loading model object -----
  ecokit::cat_time("Loading model object")

  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model file not found", path_model = path_model, include_backtrace = TRUE)
  }
  Model <- ecokit::load_as(path_model)

  # # ..................................................................... ###

  # Plot phylogenetic tree -----
  ecokit::cat_time("phylogenetic tree plot")

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
  ecokit::cat_time("Plotting theme")

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

  # Computing posterior quantities -----
  ecokit::cat_time("Computing posterior quantities")

  # Calculates mean, support and other posterior quantities for a specified
  # model parameter
  post <- Hmsc::getPostEstimate(hM = Model, parName = "Beta")

  # # ..................................................................... ###

  CovNames <- Model$covNames %>%
    stringr::str_remove("stats::poly\\(") %>%
    stringr::str_replace_all(", degree = 2, raw = TRUE\\)", "_") %>%
    stringr::str_replace_all("_1", "\n(L)") %>%
    stringr::str_replace_all("_2", "\n(Q)")
  var_order <- c(1, 1 + gtools::mixedorder(CovNames[-1]))
  CovNames <- CovNames[var_order]

  ColNames <- dplyr::case_when(
    CovNames == "(Intercept)" ~ "Intercept",
    CovNames == "EffortsLog" ~ "Sampling\nefforts",
    CovNames == "RoadRailLog" ~ "Road &\nRail\nintensity",
    CovNames == "RiversLog" ~ "River\nlength",
    CovNames == "HabLog" ~ "Habitat\ncoverage",
    CovNames == "soil" ~ "Soil\nbulk\ndensity",
    CovNames == "wetness" ~ "Topographic\nwetness\nindex",
    .default = CovNames) %>%
    paste0("\n\n\n", .)
  # Pad each string at the end so each has max_newlines
  n_newlines <- stringr::str_count(ColNames, "\n")
  ColNames <- paste0(
    ColNames, stringr::str_dup("\n", max(n_newlines) - n_newlines))

  post_support <- post$support[var_order, ]
  post_mean <- post$mean[var_order, ]

  SupportMatrix <- (post_support > support_level) %>%
    magrittr::add(post_support < (1 - support_level)) %>%
    magrittr::is_greater_than(0)

  # Legend colours
  Palette_Positive <- grDevices::colorRampPalette(
    c("#FFE4B2", "#FFC85C", "#FF8C00", "#E25822", "#800000"))(100)
  Palette_Negative <- grDevices::colorRampPalette(
    c("#E6FFFF", "#91D8F7", "#4682B4", "#083D77", "#001F3F"))(100) %>%
    rev()

  rm(Model, post, envir = environment())

  # # ..................................................................... ###

  # Support ------
  ecokit::cat_time("1. support", level = 1L)

  Plot_SupportD <- (SupportMatrix * (2 * post_support - 1)) %>%
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
    filename = fs::path(Path_Out, "Parameter_Beta_Support.jpeg"),
    res = 600, width = width, height = height, units = "cm",
    quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Sign ------

  ecokit::cat_time("2. sign", level = 1L)

  PosSign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  NegSign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  LegendTitle <- paste0(
    '<span style="font-size: 12pt"><b>Beta</b></span>',
    '<br><span style="font-size: 9pt">(sign)</span>')

  Plot_SignD <- (SupportMatrix * sign(post_mean)) %>%
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
    filename = fs::path(Path_Out, "Parameter_Beta_Sign.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  ecokit::cat_time("3. mean", level = 1L)

  Plot_MeanD <- (SupportMatrix * post_mean) %>%
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
    filename = fs::path(Path_Out, "Parameter_Beta_Mean1.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean - without intercept -----
  ecokit::cat_time("4. Mean - without intercept", level = 1L)

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
    filename = fs::path(Path_Out, "Parameter_Beta_Mean2.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  ecokit::cat_diff(init_time = .start_time, prefix = "Plotting took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}

## |------------------------------------------------------------------------| #
# mod_heatmap_omega ----
## |------------------------------------------------------------------------| #

#' @export
#' @name parameter_heatmap
#' @rdname parameter_heatmap
#' @order 2

mod_heatmap_omega <- function(
    path_model = NULL, support_level = 0.95, width = 26, height = 22.5) {

  # # ..................................................................... ###

  # Set null device for `cairo`. This is to properly render the plots using
  # ggtext - https://github.com/wilkelab/cowplot/issues/73
  cowplot::set_null_device("cairo")

  .start_time <- lubridate::now(tzone = "CET")

  if (is.null(path_model)) {
    ecokit::stop_ctx(
      "`path_model` cannot be empty", path_model = path_model,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # support_level has to be between 0 and 1
  if (support_level < 0 || support_level > 1) {
    ecokit::stop_ctx(
      "`support_level` has to be a numeric value between 0 and 1",
      support_level = support_level, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Out path -----

  Path_Out <- dirname(dirname(path_model)) %>%
    fs::path("Model_Postprocessing", "Parameters_Summary")
  fs::dir_create(Path_Out)

  # # ..................................................................... ###

  # Loading model object ------
  ecokit::cat_time("Loading model object")

  Model <- ecokit::load_as(path_model)

  if (length(Model$rLNames) == 0) {
    ecokit::cat_time("There is no Omega parameters in this model")
    return(invisible(NULL))
  }

  # # ..................................................................... ###

  # Plot phylogenetic tree -----
  ecokit::cat_time("Phylogenetic tree plot")

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
  ecokit::cat_time("Plotting theme")

  Theme <- ggplot2::theme(
    legend.title = ggtext::element_markdown(),
    legend.spacing = ggplot2::unit(0, "cm"),
    legend.key.size = ggplot2::unit(0.65, "cm"),
    legend.key.width = ggplot2::unit(0.65, "cm"),
    legend.box.margin = ggplot2::margin(0, -30, 0, -15),
    legend.box.spacing = ggplot2::unit(0, "pt"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, -3, 0.25, -2), "lines"))

  # # ..................................................................... ###

  # ComputeAssociations -----
  ecokit::cat_time("Compute associations")

  post <- Hmsc::computeAssociations(Model)[[1]]
  post_support <- post$support
  post_mean <- post$mean
  rm(Model, post, envir = environment())

  # # ..................................................................... ###

  # Sign ------
  ecokit::cat_time("1. sign", level = 1L)

  Support <- (post_support > support_level) %>%
    magrittr::add(post_support < (1 - support_level)) %>%
    magrittr::is_greater_than(0)
  # remove prefix "Sp_" from co-occurrence labels
  dimnames(Support)[[2]] <- dimnames(Support)[[1]] <- stringr::str_remove(
    dimnames(Support)[[1]], "^Sp_")

  PostMean <- post_mean
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
      # Replace diagonal elements (self-associations) with NA for clarity in
      # the heatmap
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
    filename = fs::path(Path_Out, "Parameter_Omega_Sign.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  ecokit::cat_time("2. mean", level = 1L)

  # Legend colours
  Palette_Positive <- grDevices::colorRampPalette(
    c("#FFE4B2", "#FFC85C", "#FF8C00", "#E25822", "#800000"))(100)
  Palette_Negative <- grDevices::colorRampPalette(
    c("#E6FFFF", "#91D8F7", "#4682B4", "#083D77", "#001F3F"))(100) %>%
    rev()

  LegendTitle <- paste0(
    '<span style="font-size: 11pt"><b>Omega</span><br>',
    '<span style="font-size: 8pt">(mean)</span>')

  PostMeanD <- as.data.frame(PostMean) %>%
    # replace diagonals with NA
    replace(., col(.) == row(.), NA_real_)

  PostMeanD_Positive <- PostMeanD_Negative <- PostMeanD
  PostMeanD_Positive[PostMeanD_Positive < 0] <- NA_real_
  PostMeanD_Negative[PostMeanD_Negative >= 0] <- NA_real_

  suppressMessages(
    {
      Plot_Mean1 <- ggtree::gheatmap(
        PhyloPlot, PostMeanD_Negative, offset = 0.75, width = 12,
        font.size = 0.75, colnames_offset_y = -1, colnames_angle = 90,
        hjust = 1) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Negative,
          guide = ggplot2::guide_colorbar(order = 0),
          labels = scales::number_format(accuracy = 0.1)) +
        ggplot2::labs(fill = LegendTitle) +
        ggnewscale::new_scale_fill()

      Plot_Mean <- ggtree::gheatmap(
        Plot_Mean1, PostMeanD_Positive, offset = 0.75, width = 12,
        font.size = 0.75, colnames_offset_y = -1, colnames_angle = 90,
        hjust = 1) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = Palette_Positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL,
          labels = scales::number_format(accuracy = 0.1)) +
        ggtree::scale_x_ggtree() +
        ggplot2::coord_cartesian(clip = "off")  +
        Theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  Plot <- cowplot::plot_grid(
    (Plot_Mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(Plot_Mean)),
    rel_widths = c(1, 0.09))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_Out, "Parameter_Omega_Mean.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  ecokit::cat_diff(init_time = .start_time, prefix = "Plotting took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
