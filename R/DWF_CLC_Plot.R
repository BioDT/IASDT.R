## |------------------------------------------------------------------------| #
# CLC_Plot ------
## |------------------------------------------------------------------------| #

#' Plot Corine Land Cover Maps
#'
#' This function plots Corine Land Cover maps with percentage coverage and saves
#' the plots as JPEG files.
#'
#' @param CLC_Name Character. Name of the Corine Land Cover map to plot. This
#'   has to be one of `PercCov_SynHab`, `PercCov_CLC_L1`, `PercCov_CLC_L2`,
#'   `PercCov_CLC_L3` or `PercCov_EUNIS_2019`.
#' @param CLC_Map A tibble created within the [IASDT.R::CLC_Process] containing
#'   CLC summary maps.
#' @param EU_Map `sf` object. Map of EU boundaries to overlay.
#' @param CrossWalk Data frame. Contains the crosswalk between CLC codes and
#'   their labels.
#' @param Path_JPEG Character. Directory path where JPEG files will be saved.
#' @param Path_JPEG_Free Character. Directory path where additional JPEG files
#'   with free legend will be saved (not bounded between 0 and 1).
#' @name CLC_Plot
#' @noRd
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @return NULL. Plots are saved as JPEG files.
#' @note This function is marked as internal and not intended for direct use by
#'   end users.

CLC_Plot <- function(
    CLC_Name, CLC_Map, EU_Map, CrossWalk, Path_JPEG, Path_JPEG_Free) {
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Level <- Map_Crop <- ID <- Label <- Name <- NULL

  # # ..................................................................... ###

  if (is.null(CLC_Name) || is.null(EU_Map) || is.null(CrossWalk) ||
    is.null(Path_JPEG) || is.null(Path_JPEG_Free)) {
    stop(
      paste0(
        "`CLC_Name`, `EU_Map`, `CrossWalk`, `Path_JPEG`, and ",
        "`Path_JPEG_Free` can not be empty"),
      call. = FALSE)
  }

  CLC_MapR <- dplyr::filter(CLC_Map, Name == CLC_Name) %>%
    dplyr::pull(Map_Crop) %>%
    magrittr::extract2(1) %>%
    terra::unwrap()

  Labels <- stringr::str_remove_all(CLC_Name, "PercCov_|_Crop")
  ClassOrder <- stringr::str_remove_all(names(CLC_MapR), paste0(Labels, "_"))
  Labels <- CrossWalk %>%
    dplyr::select(
      tidyselect::matches(
        match = paste0("^", Labels, "$|^", Labels, "_Label"))) %>%
    dplyr::distinct() %>%
    stats::setNames(c("Level", "Label")) %>%
    dplyr::arrange(factor(Level, levels = ClassOrder)) %>%
    dplyr::mutate(ID = seq_len(dplyr::n())) %>%
    dplyr::select(ID, Level, Label)

  Prefix <- stringr::str_remove_all(CLC_Name, "PercCov_|_Crop") %>%
    stringr::str_replace_all("CLC_", "CLC  ")

  FilePrefix <- stringr::str_remove_all(CLC_Name, "PercCov_|_Crop") %>%
    stringr::str_replace_all("CLC_L", "CLC")

  IASDT.R::CatTime(Prefix, Level = 1)

  # determine which layers will be plotted in each figure (4 columns * 2 rows)
  split_vector <- seq_len(terra::nlyr(CLC_MapR)) %>%
    split(., ceiling(seq_along(.) / 8))

  # Plotting boundaries
  Xlim <- c(2600000, 6550000)
  Ylim <- c(1450000, 5420000)

  LastUpdate <- paste0("Last update: ", format(Sys.Date(), "%d %B %Y"))

  OutPath <- paste0(
    "Summary_PercCover_", FilePrefix, "_",
    seq_len(length(split_vector)), ".jpeg") %>%
    file.path(Path_JPEG, .)

  MAPS <- purrr::map(
    .x = seq_len(length(split_vector)),
    .f = ~ {
      purrr::map(
        .x = split_vector[[.x]],
        .f = function(YY) {
          CurrMap <- CLC_MapR[[YY]]

          IASDT.R::CatTime(paste0(Labels$Label[[YY]]), Level = 2)
          MapTitle <- Labels$Label[[YY]] %>%
            # split long title text into multiple lines when necessary
            stringi::stri_wrap(55) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(MapTitle, "\n", negate = TRUE)) {
            MapTitle <- paste0(MapTitle, "\n")
          }

          # create ggplot object for each layer
          CurrMapPlot <- ggplot2::ggplot() +
            tidyterra::geom_spatraster(data = CurrMap) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::geom_sf(
              EU_Map,
              mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.2)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
            ggplot2::theme_bw() +
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
              legend.position = "none") +
            ggplot2::labs(title = MapTitle, fill = NULL)

          LevelTxt <- stringr::str_replace_all(
            string = Labels$Level[YY], pattern = "\\.", replacement = "_") %>%
            stringr::str_remove("_$")

          TilePath <- paste0(
            "PercCover_", FilePrefix, "_", LevelTxt, "_",
            Labels$Label[YY], ".jpeg") %>%
            stringr::str_replace_all("/", "_") %>%
            file.path(Path_JPEG, .)

          Theme2 <- ggplot2::theme_minimal() +
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

          TitleLab <- paste0(
            FilePrefix, " \U2014 ", Labels$Level[YY],
            ".", Labels$Label[[YY]]) %>%
            stringr::str_replace("CLC", "CLC \U2014 Level ") %>%
            stringr::str_replace(" - ", " \U2014 ") %>%
            stringr::str_replace("\\.\\.", ".") %>%
            stringi::stri_wrap(75) %>%
            stringr::str_c(collapse = "\n")

          if (stringr::str_detect(TitleLab, "\n", negate = TRUE)) {
            TitleLab <- paste0(TitleLab, "\n")
          }

          (ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_Map, mapping = ggplot2::aes(),
              color = "grey60", linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = CurrMap) +
            ggplot2::geom_sf(
              EU_Map, fill = "transparent",
              mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", "viridis::plasma",
              limits = c(0, 100)) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)),
              limits = Xlim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)),
              limits = Ylim) +
            ggplot2::labs(
              title = TitleLab,
              fill = NULL,
              tag = LastUpdate) +
            Theme2) %>%
            ggplot2::ggsave(
              filename = TilePath, width = 25, height = 23, units = "cm",
              dpi = 600, create.dir = TRUE)

          TilePathFree <- paste0(
            "PercCover_", FilePrefix, "_", Labels$Label[YY], ".jpeg") %>%
            stringr::str_replace_all("/", "_") %>%
            file.path(Path_JPEG_Free, .)

          (ggplot2::ggplot() +
            ggplot2::geom_sf(
              EU_Map, mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE,
              fill = scales::alpha("grey80", 0.4)) +
            tidyterra::geom_spatraster(data = CurrMap) +
            ggplot2::geom_sf(
              EU_Map, fill = "transparent",
              mapping = ggplot2::aes(), color = "grey60",
              linewidth = 0.25, inherit.aes = TRUE) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", palette = "viridis::plasma") +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
            ggplot2::labs(
              title = TitleLab,
              fill = NULL,
              tag = LastUpdate) +
            Theme2) %>%
            ggplot2::ggsave(
              filename = TilePathFree, width = 25, height = 23,
              units = "cm", dpi = 600, create.dir = TRUE)

          return(CurrMapPlot)
        }
      )
    }
  )

  IASDT.R::CatTime(paste0(Prefix, " - Multiple panels per file "), Level = 1)

  CommonLegend <- cowplot::get_legend(
    (ggplot2::ggplot() +
      tidyterra::geom_spatraster(
        data = terra::rast(CLC_MapR[[1]]), maxcell = terra::ncell(CLC_MapR)) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", palette = "viridis::plasma",
        limits = c(0, 100)) +
      ggplot2::theme(
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.text = ggplot2::element_text(size = 6),
        legend.background = ggplot2::element_rect(fill = "transparent")) +
      ggplot2::labs(fill = NULL))
  ) %>%
    suppressWarnings()

  # arrange map tiles together into figures (4 columns * 2 rows)
  purrr::walk(
    .x = seq_along(MAPS),
    .f = ~ {
      # main title of the figure - {("\U00D7")} prints the multiplication symbol
      MainTitle <- stringr::str_glue(
        "Percent coverage of {Prefix} per 10\u00D710 km grid cell") %>%
        as.character()

      MainTitle <- cowplot::ggdraw() +
        cowplot::draw_label(MainTitle, fontface = "bold", hjust = 0.5) +
        cowplot::draw_label(
          LastUpdate, fontface = "italic", color = "grey",
          x = 0.935, size = 3) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

      cowplot::plot_grid(plotlist = MAPS[[.x]], ncol = 4, nrow = 2) %>%
        cowplot::plot_grid(CommonLegend, rel_widths = c(4, .2)) %>%
        cowplot::plot_grid(
          MainTitle, ., ncol = 1, rel_heights = c(0.05, 1)) %>%
        ggplot2::ggsave(
          filename = OutPath[.x], width = 28, height = 15, units = "cm",
          dpi = 600, create.dir = TRUE)
    }
  )

  return(invisible(NULL))
}
