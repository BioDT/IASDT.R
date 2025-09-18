# # |------------------------------------------------------------------------| #
# gbif_species_data ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name gbif_data
#' @rdname gbif_data
#' @order 4

gbif_species_data <- function(
    species = NULL, env_file = ".env", verbose = TRUE, plot_tag = NULL) {

  # # ..................................................................... ###

  # Checking arguments ----
  if (verbose) {
    ecokit::cat_time("Checking arguments")
  }

  ecokit::check_args(
    args_to_check = c("species", "plot_tag"), args_type = "character")
  ecokit::check_args(args_to_check = "verbose", args_type = "logical")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_gbif <- CellCode <- path_grid <- eu_boundaries <- publisher_type <-
    x <- y <- NULL

  # # ..................................................................... ###

  # Environment variables
  if (verbose) {
    ecokit::cat_time("Environment variables")
  }

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_gbif", "DP_R_gbif_processed", FALSE, FALSE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # # grid_10_land_sf
  grid_sf <- fs::path(path_grid, "grid_10_land_sf.RData")
  if (!file.exists(grid_sf)) {
    ecokit::stop_ctx(
      "Reference grid file (sf) not found", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }
  grid_sf <- ecokit::load_as(grid_sf)

  # grid_10_land_crop
  grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Reference grid file not found", grid_r = grid_r,
      include_backtrace = TRUE)
  }
  grid_r <- ecokit::load_as(grid_r, unwrap_r = TRUE)

  eu_borders <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # # ..................................................................... ###

  path_species_data <- fs::path(path_gbif, "species_data")
  path_jpeg <- fs::path(path_gbif, "species_jpeg_maps")
  path_grids <- fs::path(path_gbif, "species_grids")
  path_raster <- fs::path(path_gbif, "species_raster")
  fs::dir_create(c(path_species_data, path_jpeg, path_grids, path_raster))

  if (verbose) {
    ecokit::cat_time(species, level = 1L)
  }

  sp_name <- ecokit::replace_space(species) %>%
    # replace non-ascii multiplication symbol with x
    stringr::str_replace_all("\u00D7", "x") %>%
    stringr::str_replace_all("-", "")
  sp_data_file <- fs::path(path_species_data, paste0(sp_name, ".RData"))

  if (file.exists(sp_data_file)) {

    sp_data <- ecokit::load_as(sp_data_file)

    # ****************************************************
    # species data --- grid
    # ****************************************************

    sp_grid <- sf::st_drop_geometry(sp_data) %>%
      # showing the number of observations per grid cell
      dplyr::count(CellCode, publisher_type, name = "count", sort = TRUE) %>%
      # add species name as column
      dplyr::mutate(species_name = species) %>%
      # add grid polygon to the data
      dplyr::left_join(grid_sf, by = "CellCode") %>%
      # convert to sf object
      sf::st_as_sf()

    obj_name_grid <- paste0(sp_name, "_grid")
    ecokit::save_as(
      object = sp_grid, object_name = obj_name_grid,
      out_path = fs::path(path_grids, paste0(obj_name_grid, ".RData")))

    # ****************************************************
    # species data --- raster
    # ****************************************************

    sp_r_all <- terra::rasterize(
      x = sp_grid, y = grid_r, field = "count", fun = "sum", na.rm = TRUE) %>%
      terra::mask(grid_r) %>%
      stats::setNames(sp_name) %>%
      ecokit::set_raster_crs(crs = "epsg:3035")
    sp_r_cz <- dplyr::filter(sp_grid, publisher_type == "citizen_science")
    if (nrow(sp_r_cz) > 0) {
      sp_r_cz <- sp_r_cz %>%
        terra::rasterize(y = grid_r, field = "count", fun = "sum") %>%
        terra::mask(grid_r) %>%
        stats::setNames(paste0(sp_name, "_cz")) %>%
        ecokit::set_raster_crs(crs = "epsg:3035")
    } else {
      sp_r_cz <- terra::classify(grid_r, cbind(1, 0)) %>%
        stats::setNames(paste0(sp_name, "_cz"))
    }

    sp_r_others <- dplyr::filter(sp_grid, publisher_type == "others")
    if (nrow(sp_r_others) > 0) {
      sp_r_others <- sp_r_others %>%
        terra::rasterize(y = grid_r, field = "count", fun = "sum") %>%
        terra::mask(grid_r) %>%
        stats::setNames(paste0(sp_name, "_others")) %>%
        ecokit::set_raster_crs(crs = "epsg:3035")
    } else {
      sp_r_others <- terra::classify(grid_r, cbind(1, 0)) %>%
        stats::setNames(paste0(sp_name, "_others"))
    }

    sp_r <- c(sp_r_all, sp_r_cz, sp_r_others)

    obj_name_raster <- paste0(sp_name, "_raster")
    sp_data_r_file <- fs::path(path_raster, paste0(obj_name_raster, ".RData"))
    ecokit::save_as(
      object = terra::wrap(sp_r), object_name = obj_name_raster,
      out_path = sp_data_r_file)

    # ****************************************************
    ## Plotting species data
    # ****************************************************

    plot_theme <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.25, 0, 0, 0.1, "cm"),
        plot.title = ggplot2::element_text(
          size = 14, color = "blue", face = "bold.italic", hjust = 0.5,
          margin = ggplot2::margin(2, 0, 4, 0)),
        strip.text = ggplot2::element_text(size = 6, face = "bold"),
        legend.key.size = grid::unit(0.7, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 8),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.935, 0.9),
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        panel.spacing = grid::unit(0.3, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.1, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        plot.tag.position = c(0.99, 0.9875),
        plot.tag = ggplot2::element_text(colour = "grey", size = 9, hjust = 1),
        panel.ontop = TRUE, panel.background = ggplot2::element_rect(fill = NA))

    # Plotting limits
    x_limit <- c(2600000, 6700000)
    y_limit <- c(1450000, 5420000)

    # Number of observations
    sp_plot_file <- fs::path(path_jpeg, paste0(sp_name, ".jpeg"))
    n_grids <- stringr::str_glue(
      "{scales::label_comma()(nrow(sp_data))} ",
      "observations\n{scales::label_comma()(nrow(sp_grid))} grids")
    sp_plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        eu_borders,
        mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(data = sp_r_all, maxcell = Inf) +
      ggplot2::geom_sf(
        eu_borders,
        mapping = ggplot2::aes(), color = "grey40",
        linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = ecokit::integer_breaks()) +
      ggplot2::scale_x_continuous(
        limits = x_limit, expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(
        limits = y_limit, expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(
        title = stringr::str_replace_all(species, "_", " "),
        fill = "# observations", tag = plot_tag) +
      plot_theme

    sp_plot <- cowplot::ggdraw(sp_plot) +
      cowplot::draw_label(
        n_grids, x = 0.0775, y = 0.95, size = 8, color = "red") +
      cowplot::draw_label(
        "GBIF",
        x = 0.05, y = 0.03, size = 12, color = "red", fontface = "bold")

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = sp_plot_file, width = 25, height = 25, res = 600,
      quality = 100, units = "cm")
    print(sp_plot)
    grDevices::dev.off()

    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # citizen science vs others

    sp_plot_file_cz <- fs::path(path_jpeg, paste0(sp_name, "_cz.jpeg"))
    sp_grid_cz <- sp_grid %>%
      dplyr::filter(publisher_type == "citizen_science") %>%
      dplyr::mutate(source = "Citizen science only")
    sp_grid_others <- sp_grid %>%
      dplyr::filter(publisher_type == "others") %>%
      dplyr::mutate(
        source = "Any publisher",
        line_col = dplyr::if_else(
          CellCode %in% sp_grid_cz$CellCode, "red", "transparent"),
        line_width = dplyr::if_else(CellCode %in% sp_grid_cz$CellCode, 0.1, 0))

    # Create a dummy dataframe for legend purposes
    legend_df <- c("Citizen science only", "Any publisher")
    legend_df <- data.frame(
      publisher_type = factor(legend_df, levels = legend_df),
      x = c(Inf, Inf), y = c(Inf, Inf))

    sp_plot_cz <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        eu_borders, mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.1, fill = "grey95", show.legend = FALSE) +
      ggplot2::geom_sf(
        sp_grid_cz, mapping = ggplot2::aes(), fill = "red",
        color = "transparent", linewidth = 0.2, inherit.aes = FALSE,
        show.legend = FALSE) +
      ggplot2::geom_sf(
        sp_grid_others, mapping = ggplot2::aes(),
        color = sp_grid_others$line_col, linewidth = sp_grid_others$line_width,
        fill = "blue", inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::geom_sf(
        eu_borders, mapping = ggplot2::aes(), color = "grey40",
        linewidth = 0.075, fill = "transparent", inherit.aes = FALSE,
        show.legend = FALSE) +
      ggplot2::geom_tile(
        data = legend_df, ggplot2::aes(x = x, y = y, fill = publisher_type),
        width = 1, height = 1, show.legend = TRUE) +
      ggplot2::scale_fill_manual(
        name = NULL,
        values = c("Citizen science only" = "red", "Any publisher" = "blue"),
        labels = c("Citizen science only", "Any publisher")) +
      ggplot2::scale_x_continuous(
        limits = x_limit, expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(
        limits = y_limit, expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(
        title = stringr::str_replace_all(species, "_", " "), tag = plot_tag) +
      plot_theme +
      ggplot2::theme(
        legend.position.inside = c(0.90, 0.95),
        legend.text = ggplot2::element_text(size = 10))

    ragg::agg_jpeg(
      filename = sp_plot_file_cz, width = 25, height = 25, res = 600,
      quality = 100, units = "cm")
    print(sp_plot_cz)
    grDevices::dev.off()

  }
  return(invisible(NULL))
}
