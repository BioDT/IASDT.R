# # |------------------------------------------------------------------------| #
# road_intensity ----
## |------------------------------------------------------------------------| #

#' Calculate road intensity per grid cell
#'
#' This function downloads, processes, and analyses [GRIP global roads
#' data](https://www.globio.info/download-grip-dataset) ([Meijer et al.
#' 2018](https://iopscience.iop.org/article/10.1088/1748-9326/aabd42/meta)). The
#' function calculates the total road lengths and the distance to the nearest
#' road per grid cell (for any road type and per road type).
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @return `NULL`. The function outputs processed files to the specified
#'   directories.
#' @note
#' - The function downloads the most recent version of Global Roads Inventory
#'   Project (`GRIP`) data from the URL specified in the environment variable
#'   `DP_R_road_url`. Original data format is a zipped file containing global
#'   road data in the form of `fgdb` (`EPSG:3246`).
#' - On LUMI HPC, loading the `libarchive` module is necessary to use the
#'   `archive` R package: `module load libarchive/3.6.2-cpeGNU-23.09`
#' - The distance to roads is calculated by determining the distance from each
#'   grid
#'   cell to the nearest grid cell that overlaps with a road (not to the nearest
#'   road line). Note that this is different from calculating the actual
#'   distance to the nearest road line, which is computationally intensive and
#'   not performed in this function.
#' @name road_intensity
#' @export
#' @author Ahmed El-Gabbas

road_intensity <- function(env_file = ".env") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_road <- path_road_raw <- path_road_interim <- road_type <-
    EU_boundaries <- VarName <- road_URL <- path_grid <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_road", "DP_R_road_processed", FALSE, FALSE,
    "path_road_raw", "DP_R_road_raw", FALSE, FALSE,
    "path_road_interim", "DP_R_road_interim", FALSE, FALSE,
    "road_URL", "DP_R_road_url", FALSE, FALSE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "EU_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  fs::dir_create(c(path_road, path_road_raw, path_road_interim))

  grid_ref <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_ref)) {
    ecokit::stop_ctx(
      "The reference grid file does not exist", grid_ref = grid_ref,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # download road data ------
  ecokit::cat_time("Download road data")

  withr::local_options(timeout = 1200)
  path_download_file <- fs::path(path_road_raw, basename(road_URL))

  # Check if zip file is a valid file
  if (file.exists(path_download_file)) {
    success <- ecokit::check_zip(path_download_file)
    if (isFALSE(success)) {
      fs::file_delete(path_download_file)
    }
  } else {
    success <- FALSE
  }

  # Try downloading data for a max of 5 attempts
  attempt <- 1
  attempts <- 5

  while (isFALSE(success) && attempt <= attempts) {
    download <- try(
      expr = {
        utils::download.file(
          url = road_URL, destfile = path_download_file,
          mode = "wb", quiet = TRUE) %>%
          suppressWarnings()

        success <- ecokit::check_zip(path_download_file)
        success
      }, silent = TRUE)

    if (inherits(download, "try-error")) {
      success <- FALSE
    }

    attempt <- attempt + 1
  }

  if (isFALSE(success)) {
    ecokit::stop_ctx(
      paste0("Failed to download road data after ", attempts, " attempts"),
      road_URL = road_URL, include_backtrace = TRUE)
  }

  rm(download, envir = environment())

  # # .................................... ###

  ecokit::cat_time("Extracting files")
  archive::archive_extract(
    archive = path_download_file, dir = path_road_interim) %>%
    suppressMessages()
  rm(path_download_file, envir = environment())

  # # ..................................................................... ###

  # Processing GRIP road data ------
  ecokit::cat_time("Processing GRIP road data")

  ## Load, crop, and project GRIP data -----
  ecokit::cat_time("Load, crop, and project GRIP data", level = 1L)

  road_gdb_files <- list.files(
    path = path_road_interim, pattern = ".gdb$", full.names = TRUE)

  if (length(road_gdb_files) == 0) {
    ecokit::stop_ctx(
      "No `.gdb` files found in the directory after extraction: ",
      path_road_interim = path_road_interim, road_gdb_files = road_gdb_files,
      include_backtrace = TRUE)
  }

  road_sf <- road_gdb_files[1] %>%
    # convert to `SpatVector`, returning `SpatVectorProxy`
    terra::vect(proxy = TRUE) %>%
    # Query the data only in Europe
    terra::query(extent = terra::ext(-30, 50, 25, 75)) %>%
    # project to EPSG:3035
    terra::project("EPSG:3035") %>%
    # convert to sf object
    sf::st_as_sf()

  invisible(gc())

  # # ..................................... ###

  ## Save - RData ----
  ecokit::cat_time("Save projected data - RData", level = 1L)
  save(road_sf, file = fs::path(path_road, "road_sf.RData"))

  # # ..................................... ###

  ## One file per road type ----
  ecokit::cat_time("Save RData file per road type", level = 1L)

  tibble::tribble(
    ~road_type, ~VarName,
    1, "Highways",
    2, "Primary",
    3, "Secondary",
    4, "Tertiary",
    5, "Local") %>%
    dplyr::mutate(
      A = purrr::walk2(
        .x = road_type, .y = VarName,
        .f = ~ {
          ecokit::cat_time(paste0(.x, " - ", .y), level = 2L)
          dplyr::filter(road_sf, GP_RTP %in% .x) %>%
            dplyr::select(-GP_RTP) %>%
            terra::vect() %>%
            terra::wrap() %>%
            ecokit::save_as(
              object_name = paste0("road_sf_", .x, "_", .y),
              out_path = fs::path(
                path_road, paste0("road_sf_", .x, "_", .y, ".RData")))
          invisible(gc())
          return(invisible(NULL))
        })
    ) %>%
    invisible()

  rm(road_sf, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Road length ------
  ecokit::cat_time("Road length")

  ecokit::cat_time("Calculate Road length per road type", level = 1L)

  grid_ref <- ecokit::load_as(grid_ref, unwrap_r = TRUE)

  extract_road_summary <- function(road_type, VarName, Function = "length") {
    summary_map <- list.files(
      path = path_road, full.names = TRUE,
      pattern = paste0("^road_sf_", road_type, "_.+RData")) %>%
      ecokit::load_as(unwrap_r = TRUE) %>%
      terra::rasterizeGeom(y = grid_ref, fun = Function, unit = "km") %>%
      terra::mask(mask = grid_ref) %>%
      stats::setNames(paste0(road_type, "_", VarName))
    return(summary_map)
  }

  ecokit::cat_time("1 - Highways", level = 2L)
  grip_1 <- extract_road_summary(road_type = 1, VarName = "Highways")

  ecokit::cat_time("2 - Primary", level = 2L)
  grip_2 <- extract_road_summary(road_type = 2, VarName = "Primary")

  ecokit::cat_time("3 - Secondary", level = 2L)
  grip_3 <- extract_road_summary(road_type = 3, VarName = "Secondary")

  ecokit::cat_time("4 - Tertiary", level = 2L)
  grip_4 <- extract_road_summary(road_type = 4, VarName = "Tertiary")

  ecokit::cat_time("5 - Local", level = 2L)
  grip_5 <- extract_road_summary(road_type = 5, VarName = "Local")

  ecokit::cat_time("All roads", level = 2L)
  road_length <- (grip_1 + grip_2 + grip_3 + grip_4 + grip_5) %>%
    stats::setNames("All") %>%
    c(grip_1, grip_2, grip_3, grip_4, grip_5, .) %>%
    # Ensure that values are read from memory
    terra::toMemory()

  ecokit::cat_time("Save road length - tif", level = 1L)
  terra::writeRaster(
    x = road_length, overwrite = TRUE,
    filename = fs::path(
      path_road, paste0("road_length_", names(road_length), ".tif")))

  ecokit::cat_time("Save road length - RData", level = 1L)
  ecokit::save_as(
    object = terra::wrap(road_length),
    object_name = "road_length",
    out_path = fs::path(path_road, "road_length.RData"))

  # # ..................................................................... ###

  # Distance to roads ------
  ecokit::cat_time("Distance to roads")
  # This calculates the distance from each grid cell to the nearest grid cell
  # overlapping with a road This can be different than calculating the actual
  # distance to nearest road line; which is expected to take too much time to
  # calculate

  ecokit::cat_time("Calculate distance to roads", level = 1L)

  # suppress progress bar
  terra::terraOptions(progress = 0)

  road_distance <- purrr::map(
    .x = as.list(road_length),
    .f = ~ {
      ecokit::cat_time(names(.x), level = 2L)
      Road_Points <- terra::as.points(terra::classify(.x, cbind(0, NA)))
      terra::distance(x = .x, y = Road_Points, unit = "km") %>%
        terra::mask(grid_ref) %>%
        stats::setNames(paste0("road_distance_", names(.x))) %>%
        # Ensure that values are read from memory
        terra::toMemory()
    }) %>%
    terra::rast()

  ecokit::cat_time("Save distance to road - tif", level = 1L)
  terra::writeRaster(
    x = road_distance, overwrite = TRUE,
    filename = fs::path(path_road, paste0(names(road_distance), ".tif")))

  ecokit::cat_time("Save distance to road - RData", level = 1L)
  ecokit::save_as(
    object = terra::wrap(road_distance), object_name = "road_distance",
    out_path = fs::path(path_road, "road_distance.RData"))

  # # ..................................................................... ###

  # Plotting ------
  ecokit::cat_time("Plotting")

  EU_boundaries <- ecokit::load_as(EU_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # nolint start
  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
      plot.title = ggplot2::element_text(
        size = 12, color = "blue", face = "bold", hjust = 0.5,
        margin = ggplot2::margin(0, 0, 0, 0)),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.key.width = grid::unit(0.4, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 8),
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.85),
      legend.title = ggplot2::element_text(
        color = "black", size = 8, face = "bold"),
      legend.spacing.x = grid::unit(0.2, "cm"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE,
      panel.spacing = grid::unit(0.05, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA))
  # nolint end

  # # ..................................... ###

  ## Road length -----
  ecokit::cat_time("Road length", level = 1L)

  plots_length <- purrr::map(
    .x = terra::as.list(road_length),
    .f = ~ {
      Road <- log10(terra::classify(.x, cbind(0, NA)))
      Title <- names(.x) %>%
        stringr::str_remove("road_distance_") %>%
        stringr::str_replace("_", " - ") %>%
        paste(" roads")

      ggplot2::ggplot() +
        ggplot2::geom_sf(
          EU_boundaries,
          mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(data = Road, maxcell = Inf) +
        ggplot2::geom_sf(
          EU_boundaries,
          mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6700000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(title = Title, fill = "log10") +
        plot_theme
    }
  )

  plots_length <- patchwork::wrap_plots(plots_length, ncol = 3, nrow = 2) +
    patchwork::plot_annotation(
      title = "Road length per grid cell",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0.5, 0))))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_road, "road_length.jpeg"),
    width = 30, height = 21, res = 600, quality = 100, units = "cm")
  print(plots_length)
  grDevices::dev.off()

  rm(plots_length, envir = environment())

  # # ..................................... ###

  ## Distance to road ------
  ecokit::cat_time("Distance to road", level = 1L)
  plots_distance <- purrr::map(
    .x = terra::as.list(road_distance),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          EU_boundaries,
          mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(data = .x, maxcell = Inf) +
        ggplot2::geom_sf(
          EU_boundaries,
          mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6700000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(title = paste0(names(.x), " roads"), fill = "km") +
        plot_theme
    }
  )

  plots_distance <- patchwork::wrap_plots(plots_distance, ncol = 3, nrow = 2) +
    patchwork::plot_annotation(
      title = "Distance to nearest road",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 1, 0))))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_road, "road_distance.jpeg"),
    width = 30, height = 21, res = 600, quality = 100, units = "cm")
  print(plots_distance)
  grDevices::dev.off()

  rm(plots_distance, envir = environment())

  # ..................................................................... ###

  # Cleanup ------
  ecokit::cat_time("Cleanup")

  # Delete extracted GRIP files
  list.files(path_road_interim, full.names = TRUE, pattern = "^GRIP") %>%
    fs::file_delete()

  fs::dir_delete(c(path_road_interim, path_road_raw))

  # ..................................................................... ###

  # Function Summary ----

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing road data was finished in ", ... = "\n")

  return(invisible(NULL))
}
