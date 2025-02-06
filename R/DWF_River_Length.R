# # |------------------------------------------------------------------------| #
# River_Length ----
## |------------------------------------------------------------------------| #

#' Calculate the length of rivers in each Strahler order
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
#' STRAHLER values of 5 or higher).
#'
#' @name River_Length
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param EnvFile Character. The path to the environment file containing
#'   variables required by the function. Default is ".env".
#' @param Filename Character. The name of the ZIP file containing the river
#'   network data. Default is "EU_hydro_gpkg_eu.zip".
#' @param Cleanup Logical indicating whether to clean up temporary files from
#'   the Interim directory after finishing calculations. Default: `FALSE`.
#' @return `NULL`. The function outputs processed files to the specified
#'   directories.
#' @details The data provides at pan-European level a photo-interpreted river
#'   network, consistent of surface interpretation of water bodies (lakes and
#'   wide rivers), and a drainage model (also called Drainage Network), derived
#'   from EU-DEM, with catchments and drainage lines and nodes.
#'
#'   - **Data source**: EU-Hydro River Network Database v013 | **[DOI](https://doi.org/10.2909/393359a7-7ebd-4a52-80ac-1a18d5f3db9c)** | **[Download link](https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database)**
#'   - **Temporal extent**: 2006-2012; **Format**: Vector (GPKG); **Size**: 4 GB
#' @export
#' @author Ahmed El-Gabbas

River_Length <- function(
    FromHPC = TRUE, EnvFile = ".env", Filename = "EU_hydro_gpkg_eu.zip",
    Cleanup = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # |||||||||||||||||||||||||||||||||||
  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")
  # # |||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("EnvFile", "Filename"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical", Args = c("FromHPC", "Cleanup"))

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Rivers_Raw <- Path_Rivers_Interim <- Path_Grid <- path <- File <-
    size_original <- size <- size_on_disk <- Rivers <- STRAHLER <-
    Path_Rivers <- NULL

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Environment variables ----
  IASDT.R::CatTime("Environment variables")
  # # |||||||||||||||||||||||||||||||||||

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Rivers", "DP_R_Rivers", FALSE, FALSE,
      "Path_Rivers_Raw", "DP_R_Rivers_Raw", FALSE, FALSE,
      "Path_Rivers_Interim", "DP_R_Rivers_Interim", FALSE, FALSE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Rivers", "DP_R_Rivers", FALSE, FALSE,
      "Path_Rivers_Raw", "DP_R_Rivers_Raw", FALSE, FALSE,
      "Path_Rivers_Interim", "DP_R_Rivers_Interim", FALSE, FALSE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  fs::dir_create(c(Path_Rivers, Path_Rivers_Raw, Path_Rivers_Interim))

  RefGrid <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(RefGrid)) {
    stop(
      paste0("The reference grid file does not exist: ", RefGrid),
      call. = FALSE)
  }
  RefGrid <- terra::unwrap(IASDT.R::LoadAs(RefGrid))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Unzipping river network data ----
  IASDT.R::CatTime("Unzipping river network data")
  # # |||||||||||||||||||||||||||||||||||

  # Path to the ZIP file
  Path_Rivers_Zip <- file.path(Path_Rivers_Raw, Filename)

  # Check if ZIP file exists and not empty
  if (!file.exists(Path_Rivers_Zip) || fs::file_size(Path_Rivers_Zip) == 0) {
    stop(
      paste0(
        "The river network ZIP file is missing or empty: ", Path_Rivers_Zip),
      call. = FALSE)
  }

  # List files inside ZIP not in the interim directory, and extract them
  Rivers2extract <- archive::archive(Path_Rivers_Zip) %>%
    dplyr::select(path, size_original = size) %>%
    dplyr::mutate(
      # Use working relative path
      path = file.path(Path_Rivers_Interim, path),
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
    IASDT.R::CatTime("Extracting files", Level = 1, Time = FALSE)
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
  IASDT.R::CatTime("Processing river network files")
  # # |||||||||||||||||||||||||||||||||||

  River_Lengths <- list.files(
    Path_Rivers_Interim, pattern = ".+_GPKG.zip", full.names = TRUE) %>%
    purrr::map(
      .f = ~ {

        # Print file name without extension
        stringr::str_remove_all(
          basename(.x), "^euhydro_|_v.+_GPKG.zip$") %>%
          stringr::str_to_sentence() %>%
          IASDT.R::CatTime(Level = 1)

        # File path for processed sf data
        File_Sf <- stringr::str_replace_all(.x, "_GPKG.zip$", "_sf.RData")

        if (IASDT.R::CheckData(File_Sf, warning = FALSE)) {

          # Load data from sf file if already exists and valid
          IASDT.R::CatTime(
            "GPKG file was already processed as sf", Level = 2, Time = FALSE)
          River_sf <- IASDT.R::LoadAs(File_Sf)

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
              path2 = file.path(Path_Rivers_Interim, basename(path)),
              # Check if file exists in the interim directory
              exists = file.exists(path2),
              # File size on disk
              size_on_disk = dplyr::if_else(exists, fs::file_size(path2), 0),
              # Check for file size difference
              diff = as.numeric(size_on_disk) - size_original)

          # Extract only the missing or outdated files
          gpkg_2_Extract <- dplyr::filter(gpkg_Files, diff != 0)

          if (nrow(gpkg_2_Extract) > 0) {

            IASDT.R::CatTime("Extracting gpkg file", Level = 2, Time = FALSE)

            # `archive::archive_extract` will always keep the file structure
            # inside of the ZIP file, which is not what needed here.
            zip::unzip(
              zipfile = .x, files = gpkg_2_Extract$path,
              exdir = Path_Rivers_Interim, junkpaths = TRUE)

          } else {

            IASDT.R::CatTime(
              "gpkg file was already extracted", Level = 2, Time = FALSE)

          }

          # # .................................. #

          # Read the gpkg file as sf object
          IASDT.R::CatTime("Read the gpkg files as sf", Level = 2, Time = FALSE)

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
          IASDT.R::CatTime(
            "Saving processed data as sf", Level = 2, Time = FALSE)
          IASDT.R::SaveAs(
            InObj = River_sf, OutObj = "River_sf", OutPath = File_Sf)
        }

        # # .................................. #
        # # .................................. #

        # Rasterize the river network data
        IASDT.R::CatTime("Rasterizing", Level = 2, Time = FALSE)

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
  IASDT.R::CatTime("Merging rasterized gpkg data")
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
  IASDT.R::CatTime("Saving")
  # # |||||||||||||||||||||||||||||||||||

  # Save as RData
  IASDT.R::CatTime("Save as RData", Level = 1, Time = FALSE)
  IASDT.R::SaveAs(
    InObj = terra::wrap(River_Lengths), OutObj = "River_Length",
    OutPath = file.path(Path_Rivers, "River_Lengths.RData"))

  # Save as tif files
  IASDT.R::CatTime("Save as tif files", Level = 1, Time = FALSE)
  terra::writeRaster(
    x = River_Lengths,
    filename = file.path(
      Path_Rivers, paste0("Rivers_", names(River_Lengths), ".tif")),
    overwrite = TRUE)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Plotting ------
  IASDT.R::CatTime("Plotting")
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
    filename = file.path(Path_Rivers, "River_Length.jpeg"),
    width = 30, height = 33, res = 600, quality = 100, units = "cm")
  print(RiverPlots2)
  grDevices::dev.off()

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Clean up temporary files ------
  # # |||||||||||||||||||||||||||||||||||

  if (Cleanup) {
    IASDT.R::CatTime("Cleaning up interim files")

    try(
      expr = {
        file_paths <- list.files(
          path = normalizePath(Path_Rivers_Interim, winslash = "/"),
          pattern = "_sf.RData$|.gpkg$|_GPKG.zip$", full.names = TRUE)
        fs::file_delete(file_paths)
      },
      silent = TRUE)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Function Summary ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing river data was finished in ", ... = "\n")

  # # ..................................................................... ###

  return(invisible(NULL))
}
