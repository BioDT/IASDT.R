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
#' @param path_model Character. Path to the fitted `Hmsc` model object.
#' @param support_level Numeric. The posterior support threshold for determining
#'   which values are considered significant in the heatmap. Defaults to 0.95,
#'   indicating 95% posterior support. Values above this threshold (or below 1 -
#'   threshold for negative associations) are considered significant and will be
#'   plotted (see [Hmsc::plotBeta]).
#' @param width,height Integer. The width and height of the generated
#'   heatmaps in centimetres. Defaults to 26&times;22.5 for `omega`; 25&times;35
#'   for `beta`.
#' @return Both functions do not return a value but saves heatmap plots as JPEG
#'   files in the `model_postprocessing/parameters_summary` subdirectory.
#' @details The functions exports three types of visualisations (see
#'   [Hmsc::plotBeta]):
#' - `mean`: posterior mean estimate,
#' - `support`: statistical support level, measured by the posterior
#'   probability for a positive or negative response,
#' - `sign`: indicates whether the response is positive, negative, or neither
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

  path_out <- dirname(dirname(path_model)) %>%
    fs::path("model_postprocessing", "parameters_summary")
  fs::dir_create(path_out)

  # # ..................................................................... ###

  # Loading model object -----
  ecokit::cat_time("Loading model object")

  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model file not found", path_model = path_model, include_backtrace = TRUE)
  }
  model_obj <- ecokit::load_as(path_model)

  # # ..................................................................... ###

  # plot phylogenetic tree -----
  ecokit::cat_time("phylogenetic tree plot")

  tree <- model_obj$phyloTree
  # remove prefix "sp_" from tip labels
  tree$tip.label <- stringr::str_remove(tree$tip.label, "^sp_")
  if (length(tree$edge.length) == 2 * nrow(tree$edge)) {
    tree$edge.length <- rep(1, length(tree$edge.length) / 2)
  }

  phylo_plot <- ggtree::ggtree(
    tr = tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 1) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.3, 0, 0, 0), "lines"))

  # # ..................................................................... ###

  # Plotting theme ----
  ecokit::cat_time("Plotting theme")

  plot_theme <- ggplot2::theme(
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
  post <- Hmsc::getPostEstimate(hM = model_obj, parName = "Beta")

  # # ..................................................................... ###

  cov_names <- model_obj$covNames %>%
    stringr::str_remove("stats::poly\\(") %>%
    stringr::str_replace_all(", degree = 2, raw = TRUE\\)", "_") %>%
    stringr::str_replace_all("_1", "\n(L)") %>%
    stringr::str_replace_all("_2", "\n(Q)")
  var_order <- c(1, 1 + gtools::mixedorder(cov_names[-1]))
  cov_names <- cov_names[var_order]

  col_names <- dplyr::case_when(
    cov_names == "(Intercept)" ~ "Intercept",
    cov_names == "efforts_log" ~ "Sampling\nefforts",
    cov_names == "road_rail_log" ~ "Road &\nRail\nintensity",
    cov_names == "rivers_log" ~ "River\nlength",
    cov_names == "habitat_log" ~ "Habitat\ncoverage",
    cov_names == "soil" ~ "Soil\nbulk\ndensity",
    cov_names == "wetness" ~ "Topographic\nwetness\nindex",
    .default = cov_names) %>%
    paste0("\n\n\n", .)
  # Pad each string at the end so each has max_newlines
  n_newlines <- stringr::str_count(col_names, "\n")
  col_names <- paste0(
    col_names, stringr::str_dup("\n", max(n_newlines) - n_newlines))

  post_support <- post$support[var_order, ]
  post_mean <- post$mean[var_order, ]

  support_matrix <- (post_support > support_level) %>%
    magrittr::add(post_support < (1 - support_level)) %>%
    magrittr::is_greater_than(0)

  # Legend colours
  palette_positive <- grDevices::colorRampPalette(
    c("#FFE4B2", "#FFC85C", "#FF8C00", "#E25822", "#800000"))(100)
  palette_negative <- grDevices::colorRampPalette(
    c("#E6FFFF", "#91D8F7", "#4682B4", "#083D77", "#001F3F"))(100) %>%
    rev()

  rm(model_obj, post, envir = environment())

  # # ..................................................................... ###

  # support ------
  ecokit::cat_time("1. support", level = 1L)

  plot_support_d <- (support_matrix * (2 * post_support - 1)) %>%
    t() %>%
    as.data.frame() %>%
    replace(., . == 0, NA_real_)
  # remove prefix "sp_" from species labels
  dimnames(plot_support_d)[[1]] <- stringr::str_remove(
    dimnames(plot_support_d)[[1]], "^sp_")
  colnames(plot_support_d) <- col_names

  legend_title <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 7pt">(support)</span>')

  plot_support_d_positive <- plot_support_d_negative <- plot_support_d
  plot_support_d_positive[plot_support_d_positive < 0] <- NA_real_
  plot_support_d_negative[plot_support_d_negative >= 0] <- NA_real_

  suppressMessages(
    {
      plot_support1 <- ggtree::gheatmap(
        phylo_plot, plot_support_d_negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = legend_title) +
        ggnewscale::new_scale_fill()

      plot_support <- ggtree::gheatmap(
        plot_support1, plot_support_d_positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggtree::scale_x_ggtree() +
        ggplot2::coord_cartesian(clip = "off")  +
        plot_theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  plot <- cowplot::plot_grid(
    (plot_support + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_support)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_beta_support.jpeg"),
    res = 600, width = width, height = height, units = "cm",
    quality = 100)
  print(plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Sign ------

  ecokit::cat_time("2. sign", level = 1L)

  pos_sign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  neg_sign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  legend_title <- paste0(
    '<span style="font-size: 12pt"><b>Beta</b></span>',
    '<br><span style="font-size: 9pt">(sign)</span>')

  plot_sign_data <- (support_matrix * sign(post_mean)) %>%
    sign(x = .) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate_all(as.character) %>%
    replace(., . == "1", pos_sign) %>%
    replace(., . == "-1", neg_sign) %>%
    replace(., . == "0", NA_character_)

  # remove prefix "sp_" from species labels
  dimnames(plot_sign_data)[[1]] <- stringr::str_remove(
    dimnames(plot_sign_data)[[1]], "^sp_")
  colnames(plot_sign_data) <- col_names

  plot_sign <- (
    ggtree::gheatmap(
      phylo_plot, plot_sign_data, offset = -0.85, width = 12,
      font.size = 2.5, hjust = 0.5) +
      ggplot2::scale_fill_manual(
        values = c("#E25822", "#083D77"), na.value = "transparent",
        breaks = c(pos_sign, neg_sign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = legend_title) +
      plot_theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 6))) %>%
    suppressMessages()

  plot <- cowplot::plot_grid(
    (plot_sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_sign)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_beta_sign.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  ecokit::cat_time("3. mean", level = 1L)

  plot_mean_data <- (support_matrix * post_mean) %>%
    t() %>%
    as.data.frame() %>%
    replace(., . == 0, NA_real_)
  # remove prefix "sp_" from species labels
  dimnames(plot_mean_data)[[1]] <- stringr::str_remove(
    dimnames(plot_mean_data)[[1]], "^sp_")
  colnames(plot_mean_data) <- col_names

  legend_title <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 9pt">(mean)</span>')

  plot_mean_data_positive <- plot_mean_data_negative <- plot_mean_data
  plot_mean_data_positive[plot_mean_data_positive < 0] <- NA_real_
  plot_mean_data_negative[plot_mean_data_negative >= 0] <- NA_real_

  suppressMessages(
    {
      plot_mean_1 <- ggtree::gheatmap(
        phylo_plot, plot_mean_data_negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggtree::scale_x_ggtree() +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = legend_title) +
        ggnewscale::new_scale_fill()

      plot_mean <- ggtree::gheatmap(
        plot_mean_1, plot_mean_data_positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggplot2::coord_cartesian(clip = "off")  +
        plot_theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  plot <- cowplot::plot_grid(
    (plot_mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_mean)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_beta_mean_1.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean - without intercept -----
  ecokit::cat_time("4. Mean - without intercept", level = 1L)

  legend_title <- paste0(
    '<span style="font-size: 12pt"><b>Beta</span><br>',
    '<span style="font-size: 9pt">(mean)</span><br>',
    '<span style="font-size: 7pt">[excl. Intercept]</span>')

  plot_mean_data_positive <- plot_mean_data_positive[, -1]
  plot_mean_data_negative <- plot_mean_data_negative[, -1]

  suppressMessages(
    {
      plot_mean_1 <- ggtree::gheatmap(
        phylo_plot, plot_mean_data_negative, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggtree::scale_x_ggtree() +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_negative,
          guide = ggplot2::guide_colorbar(order = 0)) +
        ggplot2::labs(fill = legend_title) +
        ggnewscale::new_scale_fill()

      plot_mean <- ggtree::gheatmap(
        plot_mean_1, plot_mean_data_positive, offset = -0.85, width = 12,
        font.size = 2.5, hjust = 0.5) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL) +
        ggplot2::coord_cartesian(clip = "off")  +
        plot_theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  plot <- cowplot::plot_grid(
    (plot_mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_mean)),
    rel_widths = c(0.94, 0.06))

  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_beta_mean_2.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(plot)
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

  path_out <- dirname(dirname(path_model)) %>%
    fs::path("model_postprocessing", "parameters_summary")
  fs::dir_create(path_out)

  # # ..................................................................... ###

  # Loading model object ------
  ecokit::cat_time("Loading model object")

  model_obj <- ecokit::load_as(path_model)

  if (length(model_obj$rLNames) == 0) {
    ecokit::cat_time("There is no Omega parameters in this model")
    return(invisible(NULL))
  }

  # # ..................................................................... ###

  # plot phylogenetic tree -----
  ecokit::cat_time("Phylogenetic tree plot")

  tree <- model_obj$phyloTree
  # Remove the 'sp_' prefix from tip labels
  tree$tip.label <- stringr::str_remove(tree$tip.label, "^sp_")
  if (length(tree$edge.length) == 2 * nrow(tree$edge)) {
    tree$edge.length <- rep(1, length(tree$edge.length) / 2)
  }

  phylo_plot <- ggtree::ggtree(
    tr = tree, branch.length = "none", ladderize = FALSE, linewidth = 0.25) +
    ggtree::geom_tiplab(size = 1) +
    ggtree::theme_tree() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.2, 0, 0, 0), "lines"))

  # # ..................................................................... ###

  # Plotting theme ----
  ecokit::cat_time("Plotting theme")

  plot_theme <- ggplot2::theme(
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

  post <- Hmsc::computeAssociations(model_obj)[[1]]
  post_support <- post$support
  post_mean <- post$mean
  rm(model_obj, post, envir = environment())

  # # ..................................................................... ###

  # Sign ------
  ecokit::cat_time("1. sign", level = 1L)

  support <- (post_support > support_level) %>%
    magrittr::add(post_support < (1 - support_level)) %>%
    magrittr::is_greater_than(0)
  # remove prefix "sp_" from co-occurrence labels
  dimnames(support)[[2]] <- dimnames(support)[[1]] <- stringr::str_remove(
    dimnames(support)[[1]], "^sp_")

  post_mean[!support] <- NA_real_
  # remove prefix "sp_" from co-occurrence labels
  dimnames(post_mean)[[2]] <- dimnames(post_mean)[[1]] <- stringr::str_remove(
    dimnames(post_mean)[[1]], "^sp_")

  pos_sign <- '<span style="font-size: 8pt"><b>  +  </b></span>'
  neg_sign <- '<span style="font-size: 8pt"><b>  \u2212  </b></span>'
  legend_title <- paste0(
    '<span style="font-size: 12pt"><b>Omega</b></span><br>',
    '<span style="font-size: 9pt">(sign)</span>')

  plot_sign <- (
    sign(post_mean) %>%
      as.data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      # Replace diagonal elements (self-associations) with NA for clarity in
      # the heatmap
      replace(., col(.) == row(.), NA_character_) %>%
      replace(., . == "1", pos_sign) %>%
      replace(., . == "-1", neg_sign) %>%
      ggtree::gheatmap(
        phylo_plot, ., offset = 0.75, width = 12, font.size = 0.75,
        colnames_offset_y = -1, colnames_angle = 90, hjust = 0.5) +
      ggplot2::scale_fill_manual(
        values = c("red", "blue"), na.value = "transparent",
        breaks = c(pos_sign, neg_sign)) +
      ggtree::scale_x_ggtree() +
      ggplot2::coord_cartesian(clip = "off")  +
      ggplot2::labs(fill = legend_title) +
      plot_theme +
      ggplot2::theme(legend.text = ggtext::element_markdown(size = 8))) %>%
    # suppress the message: Scale for fill is already present. Adding another
    # scale for fill, which will replace the existing scale.
    suppressMessages()

  plot <- cowplot::plot_grid(
    (plot_sign + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_sign)),
    rel_widths = c(1, 0.09))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_omega_sign.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Mean -----
  ecokit::cat_time("2. mean", level = 1L)

  # Legend colours
  palette_positive <- grDevices::colorRampPalette(
    c("#FFE4B2", "#FFC85C", "#FF8C00", "#E25822", "#800000"))(100)
  palette_negative <- grDevices::colorRampPalette(
    c("#E6FFFF", "#91D8F7", "#4682B4", "#083D77", "#001F3F"))(100) %>%
    rev()

  legend_title <- paste0(
    '<span style="font-size: 11pt"><b>Omega</span><br>',
    '<span style="font-size: 8pt">(mean)</span>')

  post_mean_d <- as.data.frame(post_mean) %>%
    # replace diagonals with NA
    replace(., col(.) == row(.), NA_real_)

  post_mean_d_positive <- post_mean_d_negative <- post_mean_d
  post_mean_d_positive[post_mean_d_positive < 0] <- NA_real_
  post_mean_d_negative[post_mean_d_negative >= 0] <- NA_real_

  suppressMessages(
    {
      plot_mean_1 <- ggtree::gheatmap(
        phylo_plot, post_mean_d_negative, offset = 0.75, width = 12,
        font.size = 0.75, colnames_offset_y = -1, colnames_angle = 90,
        hjust = 1) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_negative,
          guide = ggplot2::guide_colorbar(order = 0),
          labels = scales::number_format(accuracy = 0.1)) +
        ggplot2::labs(fill = legend_title) +
        ggnewscale::new_scale_fill()

      plot_mean <- ggtree::gheatmap(
        plot_mean_1, post_mean_d_positive, offset = 0.75, width = 12,
        font.size = 0.75, colnames_offset_y = -1, colnames_angle = 90,
        hjust = 1) +
        ggplot2::scale_fill_gradientn(
          na.value = "transparent", colours = palette_positive,
          guide = ggplot2::guide_colorbar(order = 1), name = NULL,
          labels = scales::number_format(accuracy = 0.1)) +
        ggtree::scale_x_ggtree() +
        ggplot2::coord_cartesian(clip = "off")  +
        plot_theme +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
    })

  plot <- cowplot::plot_grid(
    (plot_mean + ggplot2::theme(legend.position = "none")),
    ggpubr::as_ggplot(ggpubr::get_legend(plot_mean)),
    rel_widths = c(1, 0.09))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_out, "parameter_omega_mean.jpeg"), res = 600,
    width = width, height = height, units = "cm", quality = 100)
  print(plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  ecokit::cat_diff(init_time = .start_time, prefix = "Plotting took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
