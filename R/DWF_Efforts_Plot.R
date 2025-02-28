## |------------------------------------------------------------------------| #
# Efforts_Plot ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 6
#' @export

Efforts_Plot <- function(FromHPC = TRUE, EnvFile = ".env") {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Efforts <- EU_Bound <- NULL

  # # ..................................................................... ###

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts", FALSE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts_Local", FALSE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  File_SummaryR <- IASDT.R::Path(Path_Efforts, "Efforts_SummaryR.RData")
  if (!file.exists(File_SummaryR)) {
    stop("Summary maps cannot be loaded: ", File_SummaryR, call. = FALSE)
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
        filename = IASDT.R::Path(Path_Efforts, PlotDF$FileName[[.x]]),
        width = 31, height = 16.25, units = "cm", quality = 100, res = 600)
      print(CurrPlot)
      grDevices::dev.off()

    }
  )

  # # ..................................................................... ###

  return(invisible(NULL))
}
