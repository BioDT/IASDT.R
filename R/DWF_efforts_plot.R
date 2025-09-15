## |------------------------------------------------------------------------| #
# efforts_plot ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 6
#' @export

efforts_plot <- function(env_file = ".env") {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_efforts <- EU_boundaries <- NULL

  # # ..................................................................... ###

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_efforts", "DP_R_efforts_processed", FALSE, FALSE,
    "EU_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  file_summary_r <- fs::path(path_efforts, "efforts_summary_r.RData")
  if (!file.exists(file_summary_r)) {
    ecokit::stop_ctx(
      "Summary maps cannot be loaded: ", file_summary_r = file_summary_r,
      include_backtrace = TRUE)
  }
  efforts_summary <- ecokit::load_as(file_summary_r, unwrap_r = TRUE)

  # # ..................................................................... ###

  # nolint start
  plot_theme <- ggplot2::theme_bw() +
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
        color = "black", size = 8, face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA))

  EU_boundaries <- ecokit::load_as(EU_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")
  # nolint end

  # # ..................................................................... ###

  efforts_plots <- purrr::map(
    .x = seq_len(terra::nlyr(efforts_summary)),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EU_boundaries, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(
          data = efforts_summary[[.x]], maxcell = Inf) +
        ggplot2::geom_sf(
          data = EU_boundaries, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(fill = NULL) +
        plot_theme
    }) %>%
    stats::setNames(names(efforts_summary))

  # # ..................................................................... ###

  efforts_plots_Log <- purrr::map(
    .x = seq_len(terra::nlyr(efforts_summary)),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = EU_boundaries, mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(
          data = log10(efforts_summary[[.x]]), maxcell = Inf) +
        ggplot2::geom_sf(
          data = EU_boundaries, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6550000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(fill = "log10") +
        plot_theme
    }) %>%
    stats::setNames(names(efforts_summary))

  # # ..................................................................... ###

  plot_data <- tibble::tibble(
    plots = list(
      list(efforts_plots$n_obs, efforts_plots_Log$n_obs),
      list(efforts_plots$n_obs_native, efforts_plots_Log$n_obs_native),
      list(efforts_plots$n_sp, efforts_plots_Log$n_sp),
      list(efforts_plots$n_sp_native, efforts_plots_Log$n_sp_native))) %>%
    dplyr::mutate(
      file_name = c(
        "Efforts_GBIF_n_obs.jpeg", "Efforts_GBIF_n_obs_native.jpeg",
        "Efforts_GBIF_n_sp.jpeg", "Efforts_GBIF_n_sp_native.jpeg"
      ),
      Title = c(
        "Number of plant observations",
        "Number of plant observations (native species)",
        "Number of plant species",
        "Number of native species"))

  purrr::walk(
    .x = seq_len(nrow(plot_data)),
    .f = ~ {
      CurrPlot <- patchwork::wrap_plots(
        plot_data$plots[[.x]], ncol = 2, nrow = 1) +
        patchwork::plot_annotation(
          title = plot_data$Title[[.x]],
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 14, face = "bold", hjust = 0.5, colour = "blue",
              margin = ggplot2::margin(0.25, 0, 0.5, 0))))

      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      ragg::agg_jpeg(
        filename = fs::path(path_efforts, plot_data$file_name[[.x]]),
        width = 31, height = 16.25, res = 600, quality = 100, units = "cm")
      print(CurrPlot)
      grDevices::dev.off()

    }
  )

  # # ..................................................................... ###

  return(invisible(NULL))
}
