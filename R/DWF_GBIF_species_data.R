# # |------------------------------------------------------------------------| #
# GBIF_species_data ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 4

GBIF_species_data <- function(
    species = NULL, env_file = ".env", verbose = TRUE, plot_tag = NULL) {

  # # ..................................................................... ###

  # Checking arguments ----
  if (verbose) {
    ecokit::cat_time("Checking arguments")
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("env_file", "species", "plot_tag"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "verbose")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_GBIF <- CellCode <- Path_Grid <- EU_Bound <- NULL

  # # ..................................................................... ###

  # Environment variables
  if (verbose) {
    ecokit::cat_time("Environment variables")
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # # Grid_10_Land_Crop_sf
  GridSf <- fs::path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    ecokit::stop_ctx(
      "Reference grid file (sf) not found", GridSf = GridSf,
      include_backtrace = TRUE)
  }
  GridSf <- ecokit::load_as(GridSf)

  # Grid_10_Land_Crop
  GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    ecokit::stop_ctx(
      "Reference grid file not found", GridR = GridR, include_backtrace = TRUE)
  }
  GridR <- ecokit::load_as(GridR, unwrap_r = TRUE)

  EuroBound <- ecokit::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # # ..................................................................... ###

  Path_SpData <- fs::path(Path_GBIF, "Sp_Data")
  path_JPEG <- fs::path(Path_GBIF, "Sp_JPEG_Maps")
  Path_Grids <- fs::path(Path_GBIF, "Sp_Grids")
  Path_Raster <- fs::path(Path_GBIF, "Sp_Raster")
  fs::dir_create(c(Path_SpData, path_JPEG, Path_Grids, Path_Raster))

  if (verbose) {
    ecokit::cat_time(species, level = 1L)
  }

  SpName <- ecokit::replace_space(species) %>%
    # replace non-ascii multiplication symbol with x
    stringr::str_replace_all("\u00D7", "x") %>%
    stringr::str_replace_all("-", "")
  SpPath <- fs::path(Path_SpData, paste0(SpName, ".RData"))

  if (file.exists(SpPath)) {

    SpData <- ecokit::load_as(SpPath)

    # ****************************************************
    # species data --- grid
    # ****************************************************

    Obj_Name_grid <- paste0(SpName, "_Grid")
    SpGrid <- sf::st_drop_geometry(SpData) %>%
      # showing the number of observations per grid cell
      dplyr::count(CellCode, name = "Count", sort = TRUE) %>%
      # add species name as column
      dplyr::mutate(Species_name = species) %>%
      # add grid polygon to the data
      dplyr::left_join(GridSf, by = "CellCode") %>%
      # convert to sf object
      sf::st_as_sf()

    FilePath_grid <- fs::path(Path_Grids, paste0(Obj_Name_grid, ".RData"))

    ecokit::save_as(
      object = SpGrid, object_name = Obj_Name_grid, out_path = FilePath_grid)

    # ****************************************************
    # species data --- raster
    # ****************************************************

    Obj_Name_Raster <- paste0(SpName, "_Raster")
    Sp_R <- terra::rasterize(
      x = SpGrid, y = GridR, field = "Count", fun = "max") %>%
      terra::mask(GridR) %>%
      stats::setNames(SpName) %>%
      ecokit::set_raster_crs(crs = "epsg:3035")

    FilePath_R <- fs::path(Path_Raster, paste0(Obj_Name_Raster, ".RData"))
    ecokit::save_as(
      object = terra::wrap(Sp_R), object_name = Obj_Name_Raster,
      out_path = FilePath_R)

    # ****************************************************
    ## Plotting species data
    # ****************************************************

    PlottingTheme <- ggplot2::theme_bw() +
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

    FilePath_Plot <- fs::path(path_JPEG, paste0(SpName, ".jpeg"))

    NDataGrids <- stringr::str_glue(
      "{scales::label_comma()(nrow(SpData))} ",
      "observations\n{scales::label_comma()(nrow(SpGrid))} grids")

    SpPlot <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        EuroBound,
        mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(data = Sp_R, maxcell = Inf) +
      ggplot2::geom_sf(
        EuroBound,
        mapping = ggplot2::aes(), color = "grey40",
        linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = ecokit::integer_breaks()) +
      ggplot2::scale_x_continuous(
        limits = c(2600000, 6700000),
        expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(
        limits = c(1450000, 5420000),
        expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(
        title = stringr::str_replace_all(species, "_", " "),
        fill = "# observations", tag = plot_tag) +
      PlottingTheme

    SpPlot <- cowplot::ggdraw(SpPlot) +
      cowplot::draw_label(
        NDataGrids, x = 0.0775, y = 0.95, size = 8, color = "red") +
      cowplot::draw_label(
        "GBIF",
        x = 0.05, y = 0.03, size = 12, color = "red", fontface = "bold")

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = FilePath_Plot, width = 25, height = 25, res = 600,
      quality = 100, units = "cm")
    print(SpPlot)
    grDevices::dev.off()

  }
  return(invisible(NULL))
}
