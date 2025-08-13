# # |------------------------------------------------------------------------| #
# river_length ----
## |------------------------------------------------------------------------| #

#' Calculate the length of rivers in each Strahler order per grid cell
#'
#' This function processes EU-Hydro River Network Database to calculate the
#' length of rivers in each Strahler number. The Strahler number is used as an
#' index for river network classification, with higher numbers representing
#' larger, more significant river segments. The function reads and processes
#' zip-compressed geographic data (GPKG files), extracts relevant information
#' about river segments, computes the length of rivers for each Strahler order
#' per grid cell, and outputs the results both as raster files and RData
#' objects. The calculated length represents the total length of rivers in each
#' Strahler number or larger (e.g., for STRAHLER_5, the length of rivers with
#' Strahler values of 5 or higher).
#'
#' @name river_length
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param cleanup Logical indicating whether to clean up temporary files from
#'   the Interim directory after finishing calculations. Default: `FALSE`.
#' @return `NULL`. The function outputs processed files to the specified
#'   directories.
#' @details The data provides at pan-European level a photo-interpreted river
#'   network, consistent of surface interpretation of water bodies (lakes and
#'   wide rivers), and a drainage model (also called Drainage Network), derived
#'   from EU-DEM, with catchments and drainage lines and nodes.
#'
#'   - **Data source**: EU-Hydro River Network Database v013 |
#'   - **Temporal extent**: 2006-2012; **Format**: Vector (GPKG); **Size**: 4 GB
#' @export
#' @author Ahmed El-Gabbas
#' @references
#' - DOI: <https://doi.org/10.2909/393359a7-7ebd-4a52-80ac-1a18d5f3db9c>
#' - Download link: <https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database>

river_length <- function(env_file = ".env", cleanup = FALSE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # |||||||||||||||||||||||||||||||||||
  # Checking arguments ----
  ecokit::cat_time("Checking arguments")
  # # |||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character", args_to_check = "env_file")
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "cleanup")
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Rivers_Raw <- Path_Rivers_Interim <- Path_Grid <- path <- File <-
    size_original <- size <- size_on_disk <- Rivers <- STRAHLER <-
    Path_Rivers <- Path_Rivers_Zip <- NULL

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Environment variables ----
  ecokit::cat_time("Environment variables")
  # # |||||||||||||||||||||||||||||||||||

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Rivers", "DP_R_Rivers_processed", FALSE, FALSE,
    "Path_Rivers_Raw", "DP_R_Rivers_raw", FALSE, FALSE,
    "Path_Rivers_Interim", "DP_R_Rivers_interim", FALSE, FALSE,
    "Path_Rivers_Zip", "DP_R_Rivers_zip", FALSE, TRUE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  fs::dir_create(c(Path_Rivers, Path_Rivers_Raw, Path_Rivers_Interim))

  RefGrid <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(RefGrid)) {
    ecokit::stop_ctx(
      "The reference grid file does not exist", RefGrid = RefGrid,
      include_backtrace = TRUE)
  }
  RefGrid <- ecokit::load_as(RefGrid, unwrap_r = TRUE)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Unzipping river network data ----
  ecokit::cat_time("Unzipping river network data")
  # # |||||||||||||||||||||||||||||||||||

  # Check if ZIP file exists and not empty
  if (!file.exists(Path_Rivers_Zip) || fs::file_size(Path_Rivers_Zip) == 0) {
    ecokit::stop_ctx(
      "The river network ZIP file is missing or empty",
      Path_Rivers_Zip = Path_Rivers_Zip, include_backtrace = TRUE)
  }

  # Check the integrity of the ZIP file
  if (!ecokit::check_zip(Path_Rivers_Zip)) {
    ecokit::stop_ctx(
      "The river network ZIP file is corrupted",
      Path_Rivers_Zip = Path_Rivers_Zip, include_backtrace = TRUE)
  }


  # List files inside ZIP not in the interim directory, and extract them
  Rivers2extract <- archive::archive(Path_Rivers_Zip) %>%
    dplyr::select(path, size_original = size) %>%
    dplyr::mutate(
      # Use working relative path
      path = fs::path(Path_Rivers_Interim, path),
      # Check if file exists in the interim directory
      exists = file.exists(path),
      # File size on disk
      size_on_disk = dplyr::if_else(exists, fs::file_size(path), 0),
      # Check for file size difference
      diff = as.numeric(size_on_disk) - size_original) %>%
    # Filter files not in the interim directory
    dplyr::filter(diff != 0) %>%
    # Select the path of the files to extract
    dplyr::pull(path)

  # Extract only the missing or outdated files
  if (length(Rivers2extract) > 0) {
    ecokit::cat_time("Extracting files", level = 1L, cat_timestamp = FALSE)
    archive::archive_extract(
      Path_Rivers_Zip, dir = Path_Rivers_Interim,
      files = basename(Rivers2extract)) %>%
      suppressMessages()
  } else {
    message("  >>>  All river files are already extracted and up-to-date.")
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Processing river network files (gpkg) ----
  ecokit::cat_time("Processing river network files")
  # # |||||||||||||||||||||||||||||||||||

  River_Lengths <- list.files(
    Path_Rivers_Interim, pattern = ".+_GPKG.zip", full.names = TRUE) %>%
    purrr::map(
      .f = ~ {

        # Print file name without extension
        stringr::str_remove_all(
          basename(.x), "^euhydro_|_v.+_GPKG.zip$") %>%
          stringr::str_to_sentence() %>%
          ecokit::cat_time(level = 1L)

        # File path for processed sf data
        File_Sf <- stringr::str_replace_all(.x, "_GPKG.zip$", "_sf.RData")

        if (ecokit::check_data(File_Sf, warning = FALSE)) {

          # Load data from sf file if already exists and valid
          ecokit::cat_time(
            "GPKG file was already processed as sf",
            level = 2L, cat_timestamp = FALSE)
          River_sf <- ecokit::load_as(File_Sf)

        } else {

          # Extracting the GPKG file, if not already extracted

          # List of files inside the ZIP
          gpkg_Files <- archive::archive(.x) %>%
            # Filter files matching the pattern "euhydro_.+gpkg"
            dplyr::filter(
              stringr::str_detect(basename(path), "^euhydro_.+gpkg$")) %>%
            dplyr::select(path, size_original = size) %>%
            dplyr::mutate(
              # Use working relative path
              path2 = fs::path(Path_Rivers_Interim, basename(path)),
              # Check if file exists in the interim directory
              exists = file.exists(path2),
              # File size on disk
              size_on_disk = dplyr::if_else(exists, fs::file_size(path2), 0),
              # Check for file size difference
              diff = as.numeric(size_on_disk) - size_original)

          # Extract only the missing or outdated files
          gpkg_2_Extract <- dplyr::filter(gpkg_Files, diff != 0)

          if (nrow(gpkg_2_Extract) > 0) {

            ecokit::cat_time(
              "Extracting gpkg file", level = 2L, cat_timestamp = FALSE)

            # `archive::archive_extract` will always keep the file structure
            # inside of the ZIP file, which is not what needed here.
            zip::unzip(
              zipfile = .x, files = gpkg_2_Extract$path,
              exdir = Path_Rivers_Interim, junkpaths = TRUE)

          } else {

            ecokit::cat_time(
              "gpkg file was already extracted",
              level = 2L, cat_timestamp = FALSE)

          }

          # # .................................. #

          # Read the gpkg file as sf object
          ecokit::cat_time(
            "Read the gpkg files as sf", level = 2L, cat_timestamp = FALSE)

          River_sf <- sf::st_read(
            # Read the gpkg files
            dsn = gpkg_Files$path2, quiet = TRUE,
            # Select only `STRAHLER` (STRAHLER number) and `Shape` (line
            # geometry) columns from the `River_Net_l` layer
            query = "SELECT STRAHLER, Shape FROM River_Net_l") %>%
            suppressWarnings() %>%
            # Add the file name as a column
            dplyr::mutate(File = basename(gpkg_Files$path2), .before = 1) %>%
            # Merge lines with the same STRAHLER value
            tidyr::nest(River = -c(STRAHLER, File)) %>%
            # Filter out rows with NA STRAHLER values
            dplyr::filter(!is.na(STRAHLER)) %>%
            # Arrange by STRAHLER value in descending order
            dplyr::arrange(-STRAHLER) %>%
            dplyr::mutate(
              # Accumulate the River column to include all line segments of the
              # respective STRAHLER value or larger
              River = purrr::accumulate(River, dplyr::bind_rows),
              STRAHLER = paste0("STRAHLER_", STRAHLER))

          # Save the processed data as an sf object
          ecokit::cat_time(
            "Saving processed data as sf", level = 2L, cat_timestamp = FALSE)
          ecokit::save_as(
            object = River_sf, object_name = "River_sf", out_path = File_Sf)
        }

        # # .................................. #
        # # .................................. #

        # Rasterize the river network data
        ecokit::cat_time("Rasterizing", level = 2L, cat_timestamp = FALSE)

        River_sf %>%
          dplyr::mutate(
            # Calculate length of rivers per grid cell
            River = purrr::map2(
              .x = River, .y = STRAHLER,
              .f = ~ {
                terra::rasterizeGeom(
                  # Convert to spatVector
                  x = terra::vect(.x), y = RefGrid,
                  fun = "length", unit = "km") %>%
                  stats::setNames(.y)
              })) %>%
          # Arrange by STRAHLER value in ascending order
          dplyr::arrange(STRAHLER) %>%
          # Convert length column to multiple columns
          tidyr::pivot_wider(names_from = STRAHLER, values_from = River)

      }) %>%
    # Combine all the processed rivers into a single tibble
    dplyr::bind_rows()

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Merging rasterized gpkg data ----
  ecokit::cat_time("Merging rasterized gpkg data")
  # # |||||||||||||||||||||||||||||||||||

  River_Lengths <- River_Lengths %>%
    # Calculate the total length of rivers of all files
    dplyr::summarise(
      dplyr::across(
        .cols = -File,
        .fns = ~{
          purrr::discard(.x, is.null) %>%
            terra::rast() %>%
            sum() %>%
            list()
        })) %>%
    # Convert columns into two columns for STRAHLER value and River length
    tidyr::pivot_longer(
      cols = tidyselect::everything(), names_to = "STRAHLER",
      values_to = "Rivers") %>%
    # Change the layer name to the STRAHLER value
    dplyr::mutate(
      Rivers = purrr::map2(.x = Rivers, .y = STRAHLER, stats::setNames)) %>%
    # Extract `Rivers` column
    dplyr::pull(Rivers) %>%
    # Combine all the rivers into a single spatRaster
    terra::rast() %>%
    # Mask with reference grid
    terra::mask(RefGrid)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Saving ------
  ecokit::cat_time("Saving")
  # # |||||||||||||||||||||||||||||||||||

  # Save as RData
  ecokit::cat_time("Save as RData", level = 1L, cat_timestamp = FALSE)
  ecokit::save_as(
    object = terra::wrap(River_Lengths), object_name = "river_length",
    out_path = fs::path(Path_Rivers, "River_Lengths.RData"))

  # Save as tif files
  ecokit::cat_time("Save as tif files", level = 1L, cat_timestamp = FALSE)
  terra::writeRaster(
    x = River_Lengths,
    filename = fs::path(
      Path_Rivers, paste0("Rivers_", names(River_Lengths), ".tif")),
    overwrite = TRUE)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ------
  ecokit::cat_time("Plotting")
  # # |||||||||||||||||||||||||||||||||||

  RiverPlots <- purrr::map(
    .x = as.list(River_Lengths),
    .f = ~ {
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = .x, maxcell = Inf) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", palette = "viridis::plasma", name = NULL) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6500000), oob = scales::oob_keep) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(
          title = stringr::str_replace(names(.x), "_", " \u2265 ")) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0, 0.05, 0, 0.05, "cm"),
          plot.title = ggplot2::element_text(
            size = 14, color = "grey60", face = "bold", hjust = 0.5,
            margin = ggplot2::margin(t = 3, r = 0, b = 3, l = 0)),
          legend.key.size = grid::unit(0.75, "cm"),
          legend.key.width = grid::unit(0.5, "cm"),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.775),
          legend.background = ggplot2::element_rect(fill = "transparent"),
          legend.text = ggplot2::element_text(size = 10),
          legend.box.spacing = grid::unit(0, "pt"),
          axis.text.x = ggplot2::element_text(size = 7),
          axis.text.y = ggplot2::element_text(
            size = 7, hjust = 0.5, angle = 90),
          axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
          axis.ticks.length = grid::unit(0.04, "cm"),
          axis.title = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(
            linewidth = 0.075, colour = "grey40", linetype = 2),
          panel.border = ggplot2::element_blank(),
          panel.ontop = TRUE,
          panel.background = ggplot2::element_rect(fill = NA))
    })

  RiverPlots2 <- patchwork::wrap_plots(RiverPlots, ncol = 3, byrow = TRUE) +
    patchwork::plot_annotation(
      title = "Total River Length by Strahler Order",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(t = 2.5, r = 0, b = 5, l = 0))))

  ragg::agg_jpeg(
    filename = fs::path(Path_Rivers, "river_length.jpeg"),
    width = 30, height = 33, res = 600, quality = 100, units = "cm")
  print(RiverPlots2)
  grDevices::dev.off()

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Clean up temporary files ------
  # # |||||||||||||||||||||||||||||||||||

  if (cleanup) {
    ecokit::cat_time("Cleaning up interim files")

    try(
      expr = {
        file_paths <- list.files(
          path = ecokit::normalize_path(Path_Rivers_Interim),
          pattern = "_sf.RData$|.gpkg$|_GPKG.zip$", full.names = TRUE)
        fs::file_delete(file_paths)
      },
      silent = TRUE)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Function Summary ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing river data was finished in ", ... = "\n")

  # # ..................................................................... ###

  return(invisible(NULL))
}
