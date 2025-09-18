## |------------------------------------------------------------------------| #
# clc_plot ------
## |------------------------------------------------------------------------| #

#' Plot Corine Land Cover Maps
#'
#' This function plots Corine Land Cover maps with percentage coverage and saves
#' the plots as JPEG files.
#'
#' @param clc_name Character. Name of the Corine Land Cover map to plot. This
#'   has to be one of `perc_cover_synhab`, `perc_cover_clc_l1`,
#'   `perc_cover_clc_l2`, `perc_cover_clc_l3` or `perc_cover_eunis2019`.
#' @param clc_map A tibble created within the [IASDT.R::clc_process] containing
#'   CLC summary maps.
#' @param eu_map `sf` object. Map of EU boundaries to overlay.
#' @param crosswalk Data frame. Contains the crosswalk between CLC codes and
#'   their labels.
#' @param path_jpeg Character. Directory path where JPEG files will be saved.
#' @name clc_plot
#' @noRd
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @return NULL. Plots are saved as JPEG files.
#' @note This function is marked as internal and not intended for direct use by
#'   end users.

clc_plot <- function(clc_name, clc_map, eu_map, crosswalk, path_jpeg) {
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  level <- map_crop <- cw_id <- label <- name <- NULL

  # # ..................................................................... ###

  if (is.null(clc_name) || is.null(eu_map) || is.null(crosswalk) ||
      is.null(path_jpeg)) {
    ecokit::stop_ctx(
      "`clc_name`, `eu_map`, `crosswalk`, and `path_jpeg` can not be empty",
      clc_name = clc_name, eu_map = eu_map, crosswalk = crosswalk,
      path_jpeg = path_jpeg, include_backtrace = TRUE)
  }

  clc_name2 <- stringr::str_remove(clc_name, "perc_cover_")
  path_jpeg_2 <- fs::path(path_jpeg, clc_name2)
  path_jpeg_2_free <- fs::path(path_jpeg, paste0(clc_name2, "_free_legend"))

  clc_map_r <- dplyr::filter(clc_map, name == clc_name) %>%
    dplyr::pull(map_crop) %>%
    magrittr::extract2(1) %>%
    terra::unwrap()

  labels <- stringr::str_remove_all(clc_name, "perc_cover_|_crop")
  class_order <- stringr::str_remove_all(names(clc_map_r), paste0(labels, "_"))
  labels <- crosswalk %>%
    dplyr::select(
      tidyselect::matches(
        match = paste0("^", labels, "$|^", labels, "_Label"))) %>%
    dplyr::distinct() %>%
    stats::setNames(c("level", "label")) %>%
    dplyr::arrange(factor(level, levels = class_order)) %>%
    dplyr::mutate(cw_id = seq_len(dplyr::n())) %>%
    dplyr::select(cw_id, level, label)

  prefix <- stringr::str_remove_all(clc_name, "perc_cover_|_crop") %>%
    stringr::str_replace_all("clc_", "CLC  ")

  file_prefix <- stringr::str_remove_all(clc_name, "perc_cover_|_crop") %>%
    stringr::str_replace_all("clc_L", "CLC")

  # determine which layers will be plotted in each figure (4 columns * 2 rows)
  split_vector <- seq_len(terra::nlyr(clc_map_r)) %>%
    split(., ceiling(seq_along(.) / 8))

  # Plotting boundaries
  x_lim <- c(2600000, 6550000)
  y_lim <- c(1450000, 5420000)

  last_update <- paste0("Last update: ", format(Sys.Date(), "%d %B %Y"))

  out_path <- paste0(
    "summary_perc_cover_", file_prefix, "_",
    seq_len(length(split_vector)), ".jpeg") %>%
    fs::path(path_jpeg_2, .)

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.05, 0, 0, 0, "cm"),
      plot.title = ggplot2::element_text(
        size = 7, color = "blue", hjust = 0,
        margin = ggplot2::margin(2, 0, 2, 0)),
      strip.text = ggplot2::element_text(size = 6, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 4),
      axis.text.y = ggplot2::element_text(size = 4, hjust = 0.5, angle = 90),
      axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.04, "cm"),
      panel.spacing = grid::unit(0.3, "lines"),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      panel.border = ggplot2::element_blank(),
      legend.position = "none")

  plot_theme_2 <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.25, 0, 0, 0.05, "cm"),
      plot.title = ggplot2::element_text(
        size = 12, color = "blue", face = "bold", hjust = 0,
        margin = ggplot2::margin(0, 0, 0, 0)),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 0.5, angle = 90),
      axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.04, "cm"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.text = ggplot2::element_text(size = 8),
      legend.box.background = ggplot2::element_rect(colour = "transparent"),
      legend.background = ggplot2::element_rect(
        colour = "transparent", fill = "transparent"),
      plot.tag.position = c(0.99, 0.992),
      plot.tag = ggplot2::element_text(colour = "grey", size = 10, hjust = 1))

  maps <- purrr::map_dfr(
    .x = seq_len(length(split_vector)),
    .f = ~ {

      plots <- purrr::map_dfr(
        .x = split_vector[[.x]],
        .f = function(YY) {

          curr_map <- clc_map_r[[YY]]
          curr_map_nozero <- terra::classify(curr_map, cbind(0, NA))

          map_title <- labels$label[[YY]] %>%
            # split long title text into multiple lines when necessary
            stringi::stri_wrap(55) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(map_title, "\n", negate = TRUE)) {
            map_title <- paste0(map_title, "\n")
          }

          # create ggplot object for each layer
          curr_map_plot <- ggplot2::ggplot() +
            tidyterra::geom_spatraster(data = curr_map) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.2)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(title = map_title, fill = NULL) +
            plot_theme

          curr_plot_nozero <- ggplot2::ggplot() +
            tidyterra::geom_spatraster(data = curr_map_nozero) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.2)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(title = map_title, fill = NULL) +
            plot_theme

          level_text <- stringr::str_replace_all(
            string = labels$level[YY], pattern = "\\.", replacement = "_") %>%
            stringr::str_remove("_$")

          tile_path <- paste0(
            "perc_cover_", file_prefix, "_", level_text, "_",
            labels$label[YY], ".jpeg") %>%
            stringr::str_replace_all("/", "_") %>%
            stringr::str_replace_all(" ", "_") %>%
            stringr::str_to_lower() %>%
            fs::path(path_jpeg_2, .)
          tile_path_nozero <- stringr::str_replace(
            tile_path, ".jpeg$", "_nozero.jpeg")

          title_lab <- paste0(
            file_prefix, " \u2014 ", labels$level[YY],
            ".", labels$label[[YY]]) %>%
            stringr::str_replace("CLC", "CLC \u2014 level ") %>%
            stringr::str_replace(" - ", " \u2014 ") %>%
            stringr::str_replace_all("\\.\\.", ".") %>%
            stringr::str_replace_all("\\.([^0-9])", " - \\1")
            stringi::stri_wrap(75) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(title_lab, "\n", negate = TRUE)) {
            title_lab <- paste0(title_lab, "\n")
          }

          plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(),
              color = "grey60", linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map) +
            ggplot2::geom_sf(
              eu_map, fill = "transparent", mapping = ggplot2::aes(),
              color = "grey75", linewidth = 0.25, inherit.aes = TRUE) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", "viridis::plasma", limits = c(0, 100)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(
              title = title_lab, fill = "%", tag = last_update,
              x = NULL, y = NULL) +
            plot_theme_2
          ragg::agg_jpeg(
            filename = tile_path, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot)
          grDevices::dev.off()

          plot_nozero <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(),
              color = "grey60", linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map_nozero) +
            ggplot2::geom_sf(
              eu_map, fill = "transparent", mapping = ggplot2::aes(),
              color = "grey75", linewidth = 0.25, inherit.aes = TRUE) +
            ggtext::geom_richtext(
              mapping = ggplot2::aes(x = x, y = y, label = label),
              data = data.frame(
                x = -Inf, y = Inf, stringsAsFactors = FALSE,
                label = "Only grid cells with > 0% habitat coverage"),
              inherit.aes = FALSE, size = 5, hjust = 0, show.legend = FALSE,
              vjust = 0.9, lineheight = 0, fill = NA, label.color = NA,
              colour = "darkgrey") +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", "viridis::plasma", limits = c(0, 100)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(
              title = title_lab, fill = "%", tag = last_update,
              x = NULL, y = NULL) +
            plot_theme_2
          ragg::agg_jpeg(
            filename = tile_path_nozero, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot_nozero)
          grDevices::dev.off()

          # |||||||||||||||||||||||||||||||||||||||||||||||

          tile_path_free <- paste0(
            "perc_cover_", file_prefix, "_", level_text, "_",
            labels$label[YY], ".jpeg") %>%
            stringr::str_replace_all("/", "_") %>%
            stringr::str_replace_all(" ", "_") %>%
            stringr::str_to_lower() %>%
            stringr::str_sub(1, 35) %>%
            fs::path(path_jpeg_2_free, .)
          tile_path_free_nozero <- stringr::str_replace(
            tile_path_free, ".jpeg$", "_nozero.jpeg")

          plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map) +
            ggplot2::geom_sf(
              eu_map, fill = "transparent", mapping = ggplot2::aes(),
              color = "grey75", linewidth = 0.25, inherit.aes = TRUE) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma") +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(
              title = title_lab, fill = "%", tag = last_update,
              x = NULL, y = NULL) +
            plot_theme_2
          ragg::agg_jpeg(
            filename = tile_path_free, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot)
          grDevices::dev.off()

          plot_nozero <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              eu_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map_nozero) +
            ggplot2::geom_sf(
              eu_map, fill = "transparent", mapping = ggplot2::aes(),
              color = "grey75", linewidth = 0.25, inherit.aes = TRUE) +
            ggtext::geom_richtext(
              mapping = ggplot2::aes(x = x, y = y, label = label),
              data = data.frame(
                x = -Inf, y = Inf, stringsAsFactors = FALSE,
                label = "Only grid cells with > 0% habitat coverage"),
              inherit.aes = FALSE, size = 5, hjust = 0, show.legend = FALSE,
              vjust = 0.9, lineheight = 0, fill = NA, label.color = NA,
              colour = "darkgrey") +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma") +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(
              title = title_lab, fill = "%", tag = last_update,
              x = NULL, y = NULL) +
            plot_theme_2
          ragg::agg_jpeg(
            filename = tile_path_free_nozero, width = 25, height = 23,
            res = 600, quality = 100, units = "cm")
          print(plot_nozero)
          grDevices::dev.off()

          tibble::tibble(
            page = .x, plot_ID = YY,
            plot = list(curr_map_plot), plot_nozero = list(curr_plot_nozero))
        }) %>%
        dplyr::bind_rows()

      return(plots)

    })

  rm(
    maps, path_jpeg_2_free, plot_theme, plot_theme_2, prefix,
    out_path, x_lim, y_lim, last_update, env = environment())

  # Multiple panels per file
  common_legend <- cowplot::get_legend(
    (ggplot2::ggplot() +
      tidyterra::geom_spatraster(
        data = terra::rast(clc_map_r[[1]]), maxcell = terra::ncell(clc_map_r)) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", palette = "viridis::plasma",
        limits = c(0, 100)) +
      ggplot2::theme(
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.text = ggplot2::element_text(size = 6),
        legend.background = ggplot2::element_rect(fill = "transparent")) +
      ggplot2::labs(fill = NULL))) %>%
    suppressWarnings()

  # arrange map tiles together into figures (4 columns * 2 rows)
  purrr::walk(
    .x = seq_len(length(split_vector)),
    .f = ~ {

      # main title of the figure - {("\u00D7")} prints the multiplication symbol
      main_title <- stringr::str_glue(
        "Percent coverage of {prefix} per 10\u00D710 km grid cell") %>%
        as.character()
      main_title <- cowplot::ggdraw() +
        cowplot::draw_label(main_title, fontface = "bold", hjust = 0.5) +
        cowplot::draw_label(
          last_update, fontface = "italic", color = "grey",
          x = 0.935, size = 3) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

      plot <- dplyr::filter(maps, page == .x) %>%
        dplyr::pull(plot)
      plot <- cowplot::plot_grid(plotlist = plot, ncol = 4, nrow = 2) %>%
        cowplot::plot_grid(common_legend, rel_widths = c(4, 0.2)) %>%
        cowplot::plot_grid(main_title, ., ncol = 1, rel_heights = c(0.05, 1))
      path_jpeg_file <- out_path[.x]
      ragg::agg_jpeg(
        filename = path_jpeg_file, width = 28, height = 15, res = 600,
        quality = 100, units = "cm")
      print(plot)
      grDevices::dev.off()

      plot_nozero <- dplyr::filter(maps, page == .x) %>%
        dplyr::pull(plot_nozero)
      plot_nozero <- cowplot::plot_grid(
        plotlist = plot_nozero, ncol = 4, nrow = 2) %>%
        cowplot::plot_grid(common_legend, rel_widths = c(4, 0.2)) %>%
        cowplot::plot_grid(main_title, ., ncol = 1, rel_heights = c(0.05, 1))
      path_jpeg_file <- stringr::str_replace(
        out_path[.x], "summary_jpeg", "summary_jpeg_nozero")
      ragg::agg_jpeg(
        filename = path_jpeg_file, width = 28, height = 15, res = 600,
        quality = 100, units = "cm")
      print(plot_nozero)
      grDevices::dev.off()
    })

  rm(common_legend, env = environment())

  invisible(NULL)
}
