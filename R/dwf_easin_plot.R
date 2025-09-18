## |------------------------------------------------------------------------| #
# easin_plot ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name easin_data
#' @rdname easin_data
#' @order 4

easin_plot <- function(env_file = ".env") {

  # # ..................................................................... ###

  .plot_start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_easin_summary <- eu_boundaries <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Environment variables", level = 1L)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE,
    "path_easin", "DP_R_easin_processed", TRUE, FALSE,
    "path_easin_Interim", "DP_R_easin_interim", TRUE, FALSE,
    "path_easin_summary", "DP_R_easin_summary", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # |||||||||||||||||||||||||||||||||||
  # # Input maps ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Loading input maps", level = 1L)

  ## Country boundaries ----
  ecokit::cat_time("Country boundaries", level = 2L)
  eu_borders <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  ## EASIN summary maps -----
  ecokit::cat_time("EASIN summary maps", level = 2L)
  path_n_species <- fs::path(path_easin_summary, "easin_n_sp.RData")
  path_n_sp_partner <- fs::path(path_easin_summary, "easin_n_sp_partner.RData")

  path_n_obs <- fs::path(path_easin_summary, "easin_n_obs.RData")
  path_n_obs_partner <- fs::path(
    path_easin_summary, "easin_n_obs_partner.RData")

  path_summary_maps <- c(
    path_n_species, path_n_sp_partner, path_n_obs, path_n_obs_partner)
  summary_maps_missing <- !file.exists(path_summary_maps)

  if (any(summary_maps_missing)) {
    ecokit::stop_ctx(
      "Missing summary input files",
      missing_files = path_summary_maps[which(summary_maps_missing)],
      include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  ## NObs + NSp ----

  ecokit::cat_time("Number of species & observations", level = 1L)

  plot_easin_all <- function(
    path_map, plot_title, eu_borders, add_tag = FALSE, plot_legend = FALSE) {

    last_update <- paste0(
      "<b>Last update:</b></br><i>", format(Sys.Date(), "%d %B %Y"), "</i>")

    plot_map <- log10(ecokit::load_as(path_map, unwrap_r = TRUE))
    n_cells <- terra::ncell(plot_map)

    plot_theme <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 8, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0)),
        legend.key.size = grid::unit(0.3, "cm"),
        legend.key.width = grid::unit(0.3, "cm"),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 4),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.title = ggtext::element_markdown(
          color = "blue", size = 6, face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.92, 0.85),
        axis.text.x = ggplot2::element_text(size = 4),
        axis.text.y = ggplot2::element_text(size = 4, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        panel.spacing = grid::unit(0.3, "lines"),
        panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25),
        panel.border = ggplot2::element_blank(),
        plot.tag.position = c(0.85, 0.9999),
        plot.tag = ggtext::element_markdown(colour = "grey", size = 5))

    plot_obj <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = plot_map, maxcell = n_cells) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma") +
      ggplot2::geom_sf(
        eu_borders, inherit.aes = TRUE,
        mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.04, fill = scales::alpha("grey80", 0.2)) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)),
        limits = c(2600000, 6700000)) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)),
        limits = c(1450000, 5420000)) +
      plot_theme

    if (add_tag) {
      plot_obj <- plot_obj +
        ggplot2::labs(
          title = plot_title, fill = "log<sub>10</sub>",
          tag = last_update)
    } else {
      plot_obj <- plot_obj +
        ggplot2::labs(title = plot_title, fill = "log<sub>10</sub>")
    }
    return(plot_obj)
  }

  ### Number of observations ----
  ecokit::cat_time("Number of observations", level = 2L)
  plot_n_obs <- plot_easin_all(
    path_map = path_n_obs, plot_title = "Number of observations",
    eu_borders = eu_borders, add_tag = FALSE, plot_legend = FALSE)

  ### Number of species ----
  ecokit::cat_time("Number of species", level = 2L)
  plot_n_species <- plot_easin_all(
    path_map = path_n_species, plot_title = "Number of species",
    eu_borders = eu_borders, add_tag = TRUE, plot_legend = TRUE)

  ### Combine maps ----
  ecokit::cat_time("Merge maps side by side and save as JPEG", level = 2L)

  plot_obj <- ggpubr::ggarrange(
    plot_n_obs,
    (ggplot2::ggplot() + ggplot2::theme_void()),
    plot_n_species,
    widths = c(1, 0, 1), nrow = 1) +
    patchwork::plot_annotation(
      title = "EASIN data",
      theme = ggplot2::theme(
        plot.margin = ggplot2::margin(0.1, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 9, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0))))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_easin_summary, "easin_data.jpeg"),
    width = 20, height = 10.3, res = 600, quality = 100, units = "cm")
  print(plot_obj)
  grDevices::dev.off()

  rm(plot_n_species, plot_n_obs, plot_obj, envir = environment())

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Number of species/observations per partner ----

  ecokit::cat_time("Number of species & observations per partner", level = 1L)

  plot_easin_partner <- function(path_map, file_prefix, plot_title) {

    last_update <- paste0(
      "<b>Last update:</b> <i>", format(Sys.Date(), "%d %B %Y"), "</i>")

    plot_theme2 <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.125, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 14, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(2, 0, 4, 0)),
        panel.spacing = grid::unit(0.15, "lines"),
        panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25),
        panel.border = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 10, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 4),
        axis.text.y = ggplot2::element_text(size = 4, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        legend.text = ggplot2::element_text(size = 6),
        legend.title = ggtext::element_markdown(
          color = "blue", size = 9, face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.97, 0.90),
        legend.key.size = grid::unit(0.35, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.box.spacing = grid::unit(0, "pt"),
        plot.tag.position = c(0.92, 0.975),
        plot.tag = ggtext::element_markdown(colour = "grey", size = 9)
      )

    plot_map <- log10(ecokit::load_as(path_map, unwrap_r = TRUE))
    n_cells <- terra::ncell(plot_map)

    legend_limit <- c(
      min(terra::global(plot_map, min, na.rm = TRUE), na.rm = TRUE),
      max(terra::global(plot_map, max, na.rm = TRUE), na.rm = TRUE))

    for (i in seq_len(2)) {
      start_layer <- (i - 1) * 8 + 1
      end_layer <- min(i * 8, terra::nlyr(plot_map))

      plot_obj <- ggplot2::ggplot() +
        tidyterra::geom_spatraster(
          data = plot_map[[start_layer:end_layer]], maxcell = n_cells) +
        ggplot2::facet_wrap(~lyr, nrow = 2, ncol = 4) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma", limits = legend_limit) +
        ggplot2::geom_sf(
          eu_borders, inherit.aes = TRUE,
          mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.04, fill = scales::alpha("grey80", 0.2)) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6700000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(
          title = paste0(plot_title, "[p", i, "]"),
          fill = "log<sub>10</sub>",
          tag = last_update) +
        plot_theme2

      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      ragg::agg_jpeg(
        filename = fs::path(
          path_easin_summary, paste0(file_prefix, "_p", i, ".jpeg")),
        width = 30, height = 16.5, res = 600, quality = 100, units = "cm")
      print(plot_obj)
      grDevices::dev.off()

    }

    invisible(NULL)
  }

  ### Number of observations per partner ----
  ecokit::cat_time("Number of observations per partner", level = 2L)
  plot_easin_partner(
    path_map = path_n_obs_partner,
    file_prefix = "easin_n_obs_per_partner",
    plot_title = "EASIN data - Number of observations per data partner ")


  ### Number of species per partner ----
  ecokit::cat_time("Number of species per partner", level = 2L)
  plot_easin_partner(
    path_map = path_n_sp_partner,
    file_prefix = "easin_n_sp_per_partner",
    plot_title = "EASIN data - Number of species per data partner ")

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .plot_start_time,
    prefix = "Plotting EASIN data was finished in ", level = 1L)

  invisible(NULL)
}
