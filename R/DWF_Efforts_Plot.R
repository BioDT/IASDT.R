## |------------------------------------------------------------------------| #
# Efforts_Plot ----
## |------------------------------------------------------------------------| #

#' Plot the output of efforts maps
#'
#' This function generates and saves multiple plots of plant observation
#' efforts, both in raw and log10 scales, using provided spatial boundary and
#' summary data.
#' @param Path_Efforts Character. Path to the directory where the generated
#'   plots will be saved. The directory must exist.
#' @param EU_Bound Character. Path to `RData` file containing country
#'   boundaries.
#' @return The function saves generated plots as JPEG files in the specified
#'   directory and returns NULL invisibly.
#' @author Ahmed El-Gabbas
#' @name Efforts_Plot
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Efforts_Process] function.
#' @details This function generates and saves effort maps visualizing the number
#'   of plant observations and species, including both native and non-native
#'   species, within Europe. It produces both standard and log10-scaled plots.
#' @export

Efforts_Plot <- function(Path_Efforts, EU_Bound) {


  # # ..................................................................... ###

  File_SummaryR <- file.path(Path_Efforts, "Efforts_SummaryR.RData")
  if (!file.exists(File_SummaryR)) {
    stop(
      paste0("Summary maps cannot be loaded: ", File_SummaryR),
      call. = FALSE
    )
  }

  Efforts_SummaryR <- terra::unwrap(IASDT.R::LoadAs(File_SummaryR))

  # # ..................................................................... ###

  # nolint start
  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.02, 0, 0.02, 0, "cm"),
      plot.title = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 8),
      legend.position = "inside",
      legend.position.inside = c(0.925, 0.85),
      legend.title = ggplot2::element_text(
        color = "black", size = 8, face = "bold"
      ),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA))

  EurBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")
  # nolint end

  # # ..................................................................... ###

  Efforts_Plots <- purrr::map(
    .x = seq_len(terra::nlyr(Efforts_SummaryR)),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98"
        ) +
        tidyterra::geom_spatraster(
          data = Efforts_SummaryR[[.x]], maxcell = Inf
        ) +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE
        ) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma"
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)
        ) +
        ggplot2::labs(fill = NULL) +
        PlottingTheme
    }
  ) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  Efforts_Plots_Log <- purrr::map(
    .x = seq_len(terra::nlyr(Efforts_SummaryR)),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98"
        ) +
        tidyterra::geom_spatraster(
          data = log10(Efforts_SummaryR[[.x]]), maxcell = Inf
        ) +
        ggplot2::geom_sf(
          data = EurBound, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE
        ) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma"
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)
        ) +
        ggplot2::labs(fill = "log10") +
        PlottingTheme
    }
  ) %>%
    stats::setNames(names(Efforts_SummaryR))

  # # ..................................................................... ###

  PlotDF <- tibble::tibble(
    Plots = list(
      list(Efforts_Plots$NObs, Efforts_Plots_Log$NObs),
      list(Efforts_Plots$NObs_Native, Efforts_Plots_Log$NObs_Native),
      list(Efforts_Plots$NSp, Efforts_Plots_Log$NSp),
      list(Efforts_Plots$NSp_Native, Efforts_Plots_Log$NSp_Native)
    )
  ) %>%
    dplyr::mutate(
      FileName = c(
        "Efforts_GBIF_NObs.jpeg", "Efforts_GBIF_NObs_Native.jpeg",
        "Efforts_GBIF_NSp.jpeg", "Efforts_GBIF_NSp_Native.jpeg"
      ),
      Title = c(
        "Number of plant observations",
        "Number of plant observations (native species)",
        "Number of plant species",
        "Number of native species"
      )
    )

  purrr::walk(
    .x = seq_len(nrow(PlotDF)),
    .f = ~ {
      CurrPlot <- patchwork::wrap_plots(
        PlotDF$Plots[[.x]], ncol = 2, nrow = 1) +
        patchwork::plot_annotation(
          title = PlotDF$Title[[.x]],
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 14, face = "bold", hjust = 0.5, colour = "blue",
              margin = ggplot2::margin(0.25, 0, 0.5, 0))))

      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      grDevices::jpeg(
        filename = file.path(Path_Efforts, PlotDF$FileName[[.x]]),
        width = 31, height = 16.25, units = "cm", quality = 100, res = 600)
      print(CurrPlot)
      grDevices::dev.off()

    }
  )

  # # ..................................................................... ###

  return(invisible(NULL))
}
