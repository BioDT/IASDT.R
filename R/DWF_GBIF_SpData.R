# # |------------------------------------------------------------------------| #
# GBIF_SpData ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 5

GBIF_SpData <- function(
    Species = NULL, EnvFile = ".env", Verbose = TRUE, PlotTag = NULL) {

  # # ..................................................................... ###

  # Checking arguments ----
  if (Verbose) {
    IASDT.R::CatTime("Checking arguments")
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("EnvFile", "Species", "PlotTag"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "Verbose")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_GBIF <- CellCode <- Path_Grid <- EU_Bound <- NULL

  # # ..................................................................... ###

  # Environment variables
  if (Verbose) {
    IASDT.R::CatTime("Environment variables")
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # # Grid_10_Land_Crop_sf
  GridSf <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Reference grid file (sf) not found at: ", GridSf, call. = FALSE)
  }
  GridSf <- IASDT.R::LoadAs(GridSf)

  # Grid_10_Land_Crop
  GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop("Reference grid file not found at: ", GridR, call. = FALSE)
  }
  GridR <- terra::unwrap(IASDT.R::LoadAs(GridR))

  EuroBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # # ..................................................................... ###

  Path_SpData <- IASDT.R::Path(Path_GBIF, "Sp_Data")
  Path_JPEG <- IASDT.R::Path(Path_GBIF, "Sp_JPEG_Maps")
  Path_Grids <- IASDT.R::Path(Path_GBIF, "Sp_Grids")
  Path_Raster <- IASDT.R::Path(Path_GBIF, "Sp_Raster")
  fs::dir_create(c(Path_SpData, Path_JPEG, Path_Grids, Path_Raster))

  if (Verbose) {
    IASDT.R::CatTime(Species, Level = 1)
  }

  SpName <- IASDT.R::ReplaceSpace(Species) %>%
    # replace non-ascii multiplication symbol with x
    stringr::str_replace_all("\u00D7", "x") %>%
    stringr::str_replace_all("-", "")
  SpPath <- IASDT.R::Path(Path_SpData, paste0(SpName, ".RData"))

  if (file.exists(SpPath)) {

    SpData <- IASDT.R::LoadAs(SpPath)

    # ****************************************************
    # species data --- grid
    # ****************************************************

    Obj_Name_grid <- paste0(SpName, "_Grid")
    SpGrid <- sf::st_drop_geometry(SpData) %>%
      # showing the number of observations per grid cell
      dplyr::count(CellCode, name = "Count", sort = TRUE) %>%
      # add species name as column
      dplyr::mutate(Species_name = Species) %>%
      # add grid polygon to the data
      dplyr::left_join(GridSf, by = "CellCode") %>%
      # convert to sf object
      sf::st_as_sf()

    FilePath_grid <- IASDT.R::Path(Path_Grids, paste0(Obj_Name_grid, ".RData"))

    IASDT.R::SaveAs(
      InObj = SpGrid, OutObj = Obj_Name_grid, OutPath = FilePath_grid)

    # ****************************************************
    # species data --- raster
    # ****************************************************

    Obj_Name_Raster <- paste0(SpName, "_Raster")
    Sp_R <- terra::rasterize(
      x = SpGrid, y = GridR, field = "Count", fun = "max") %>%
      terra::mask(GridR) %>%
      stats::setNames(SpName) %>%
      IASDT.R::setRastCRS()

    FilePath_R <- IASDT.R::Path(Path_Raster, paste0(Obj_Name_Raster, ".RData"))
    IASDT.R::SaveAs(
      InObj = terra::wrap(Sp_R), OutObj = Obj_Name_Raster, OutPath = FilePath_R
    )

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

    FilePath_Plot <- IASDT.R::Path(Path_JPEG, paste0(SpName, ".jpeg"))

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
        breaks = IASDT.R::integer_breaks()) +
      ggplot2::scale_x_continuous(
        limits = c(2600000, 6700000),
        expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(
        limits = c(1450000, 5420000),
        expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(
        title = stringr::str_replace_all(Species, "_", " "),
        fill = "# observations", tag = PlotTag) +
      PlottingTheme

    SpPlot <- cowplot::ggdraw(SpPlot) +
      cowplot::draw_label(
        NDataGrids, x = 0.0775, y = 0.95, size = 8, color = "red") +
      cowplot::draw_label(
        "GBIF",
        x = 0.05, y = 0.03, size = 12, color = "red", fontface = "bold")

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::jpeg(
      filename = FilePath_Plot,
      width = 25, height = 25, units = "cm", quality = 100, res = 600)
    print(SpPlot)
    grDevices::dev.off()

  }
  return(invisible(NULL))
}
