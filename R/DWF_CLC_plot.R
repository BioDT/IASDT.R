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
#' @param EU_map `sf` object. Map of EU boundaries to overlay.
#' @param crosswalk Data frame. Contains the crosswalk between CLC codes and
#'   their labels.
#' @param path_jpeg Character. Directory path where JPEG files will be saved.
#' @param path_jpeg_free Character. Directory path where additional JPEG files
#'   with free legend will be saved (not bounded between 0 and 1).
#' @name clc_plot
#' @noRd
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @return NULL. Plots are saved as JPEG files.
#' @note This function is marked as internal and not intended for direct use by
#'   end users.

clc_plot <- function(
    clc_name, clc_map, EU_map, crosswalk, path_jpeg, path_jpeg_free) {
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  level <- map_crop <- ID <- label <- Name <- NULL

  # # ..................................................................... ###

  if (is.null(clc_name) || is.null(EU_map) || is.null(crosswalk) ||
      is.null(path_jpeg) || is.null(path_jpeg_free)) {
    ecokit::stop_ctx(
      paste0(
        "`clc_name`, `EU_map`, `crosswalk`, `path_jpeg`, and ",
        "`path_jpeg_free` can not be empty"),
      clc_name = clc_name, EU_map = EU_map, crosswalk = crosswalk,
      path_jpeg = path_jpeg, path_jpeg_free = path_jpeg_free,
      include_backtrace = TRUE)
  }

  clc_map_r <- dplyr::filter(clc_map, Name == clc_name) %>%
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
    dplyr::mutate(ID = seq_len(dplyr::n())) %>%
    dplyr::select(ID, level, label)

  prefix <- stringr::str_remove_all(clc_name, "perc_cover_|_crop") %>%
    stringr::str_replace_all("clc_", "CLC  ")

  file_prefix <- stringr::str_remove_all(clc_name, "perc_cover_|_crop") %>%
    stringr::str_replace_all("clc_L", "CLC")

  ecokit::cat_time(prefix, level = 1L)

  # determine which layers will be plotted in each figure (4 columns * 2 rows)
  split_vector <- seq_len(terra::nlyr(clc_map_r)) %>%
    split(., ceiling(seq_along(.) / 8))

  # nolint start

  # Plotting boundaries
  x_lim <- c(2600000, 6550000)
  y_lim <- c(1450000, 5420000)

  last_update <- paste0("Last update: ", format(Sys.Date(), "%d %B %Y"))

  out_path <- paste0(
    "summary_perc_cover_", file_prefix, "_",
    seq_len(length(split_vector)), ".jpeg") %>%
    fs::path(path_jpeg, .)

  maps <- purrr::map_dfr(
    .x = seq_len(length(split_vector)),
    .f = ~ {

      plots <- purrr::map_dfr(
        .x = split_vector[[.x]],
        .f = function(YY) {

          curr_map <- clc_map_r[[YY]]
          curr_map_no_zero <- terra::classify(curr_map, cbind(0, NA))

          ecokit::cat_time(paste0(labels$label[[YY]]), level = 2L)
          map_title <- labels$label[[YY]] %>%
            # split long title text into multiple lines when necessary
            stringi::stri_wrap(55) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(map_title, "\n", negate = TRUE)) {
            map_title <- paste0(map_title, "\n")
          }

          plot_theme <-  ggplot2::theme_bw() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(0.05, 0, 0, 0, "cm"),
              plot.title = ggplot2::element_text(
                size = 7, color = "blue", hjust = 0,
                margin = ggplot2::margin(2, 0, 2, 0)),
              strip.text = ggplot2::element_text(size = 6, face = "bold"),
              axis.text.x = ggplot2::element_text(size = 4),
              axis.text.y = ggplot2::element_text(
                size = 4, hjust = 0.5, angle = 90),
              axis.ticks = ggplot2::element_line(
                colour = "blue", linewidth = 0.25),
              axis.ticks.length = grid::unit(0.04, "cm"),
              panel.spacing = grid::unit(0.3, "lines"),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
              panel.grid.major = ggplot2::element_line(linewidth = 0.25),
              panel.border = ggplot2::element_blank(),
              legend.position = "none")

          # create ggplot object for each layer
          curr_map_plot <- ggplot2::ggplot() +
            tidyterra::geom_spatraster(data = curr_map) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.2)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
            ggplot2::labs(title = map_title, fill = NULL) +
            plot_theme

          curr_plot_no_zero <- ggplot2::ggplot() +
            tidyterra::geom_spatraster(data = curr_map_no_zero) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(), color = "grey60",
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
            fs::path(path_jpeg, .)

          tile_path_no_zero <- stringr::str_replace(
            tile_path, ".jpeg", "_no_zero.jpeg")

          theme_2 <- ggplot2::theme_minimal() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(0.25, 0, 0, 0.05, "cm"),
              plot.title = ggplot2::element_text(
                size = 12, color = "blue", face = "bold", hjust = 0,
                margin = ggplot2::margin(0, 0, 0, 0)),
              axis.text.x = ggplot2::element_text(size = 8),
              axis.text.y = ggplot2::element_text(
                size = 8, hjust = 0.5, angle = 90),
              axis.ticks = ggplot2::element_line(
                colour = "blue", linewidth = 0.25),
              axis.ticks.length = grid::unit(0.04, "cm"),
              legend.box.margin = ggplot2::margin(0, 0, 0, 0),
              legend.key.size = grid::unit(0.8, "cm"),
              legend.key.width = grid::unit(0.6, "cm"),
              legend.text = ggplot2::element_text(size = 8),
              legend.box.background = ggplot2::element_rect(
                colour = "transparent"),
              legend.background = ggplot2::element_rect(
                colour = "transparent", fill = "transparent"),
              plot.tag.position = c(0.99, 0.992),
              plot.tag = ggplot2::element_text(
                colour = "grey", size = 7, hjust = 1))

          title_lab <- paste0(
            file_prefix, " \u2014 ", labels$level[YY],
            ".", labels$label[[YY]]) %>%
            stringr::str_replace("CLC", "CLC \u2014 level ") %>%
            stringr::str_replace(" - ", " \u2014 ") %>%
            stringr::str_replace("\\.\\.", ".") %>%
            stringi::stri_wrap(75) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(title_lab, "\n", negate = TRUE)) {
            title_lab <- paste0(title_lab, "\n")
          }

          plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(),
              color = "grey60", linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map) +
            ggplot2::geom_sf(
              EU_map, fill = "transparent", mapping = ggplot2::aes(),
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
            theme_2
          ragg::agg_jpeg(
            filename = tile_path, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot)
          grDevices::dev.off()

          plot_no_zero <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(),
              color = "grey60", linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map_no_zero) +
            ggplot2::geom_sf(
              EU_map, fill = "transparent", mapping = ggplot2::aes(),
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
            theme_2
          ragg::agg_jpeg(
            filename = tile_path_no_zero, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot_no_zero)
          grDevices::dev.off()

          # |||||||||||||||||||||||||||||||||||||||||||||||

          tile_path_free <- paste0(
            "perc_cover_", file_prefix, "_", labels$label[YY], ".jpeg") %>%
            stringr::str_replace_all("/", "_") %>%
            fs::path(path_jpeg_free, .)
          tile_path_free_no_zero <- stringr::str_replace(
            tile_path_free, ".jpeg", "_no_zero.jpeg")

          plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map) +
            ggplot2::geom_sf(
              EU_map, fill = "transparent", mapping = ggplot2::aes(),
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
            theme_2
          ragg::agg_jpeg(
            filename = tile_path_free, width = 25, height = 23, res = 600,
            quality = 100, units = "cm")
          print(plot)
          grDevices::dev.off()

          plot_no_zero <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = curr_map) +
            ggplot2::geom_sf(
              EU_map, fill = "transparent", mapping = ggplot2::aes(),
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
            theme_2
          ragg::agg_jpeg(
            filename = tile_path_free_no_zero, width = 25, height = 23,
            res = 600, quality = 100, units = "cm")
          print(plot_no_zero)
          grDevices::dev.off()

          tibble::tibble(
            page = .x, plot_ID = YY,
            plot = list(curr_map_plot), plot_no_zero = list(curr_plot_no_zero))
        }) %>%
        dplyr::bind_rows()

      return(plots)

    })

  ecokit::cat_time(paste0(prefix, " - Multiple panels per file "), level = 1L)

  common_legend <- cowplot::get_legend(
    (ggplot2::ggplot() +
       tidyterra::geom_spatraster(
         data = terra::rast(clc_map_r[[1]]),
         maxcell = terra::ncell(clc_map_r)) +
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
  # nolint end


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
      path_jpeg <- out_path[.x]
      ragg::agg_jpeg(
        filename = path_jpeg, width = 28, height = 15, res = 600,
        quality = 100, units = "cm")
      print(plot)
      grDevices::dev.off()

      plot_no_zero <- dplyr::filter(maps, page == .x) %>%
        dplyr::pull(plot_no_zero)
      plot_no_zero <- cowplot::plot_grid(
        plotlist = plot_no_zero, ncol = 4, nrow = 2) %>%
        cowplot::plot_grid(common_legend, rel_widths = c(4, 0.2)) %>%
        cowplot::plot_grid(main_title, ., ncol = 1, rel_heights = c(0.05, 1))
      path_jpeg <- stringr::str_replace(
        out_path[.x], "summary_jpeg", "summary_jpeg_no_zero")
      ragg::agg_jpeg(
        filename = path_jpeg, width = 28, height = 15, res = 600,
        quality = 100, units = "cm")
      print(plot_no_zero)
      grDevices::dev.off()
    })

  return(invisible(NULL))
}
