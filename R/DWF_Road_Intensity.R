# # |------------------------------------------------------------------------| #
# Road_Intensity ----
## |------------------------------------------------------------------------| #

#' Calculate road intensity per grid cell
#'
#' This function downloads, processes, and analyzes [GRIP global roads
#' data](https://www.globio.info/download-grip-dataset) ([Meijer et al.
#' 2018](https://iopscience.iop.org/article/10.1088/1748-9326/aabd42/meta)). The
#' function calculates the total road lengths and the distance to the nearest
#' road per grid cell (for any road type and per road type).
#'
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @return `NULL`. The function outputs processed files to the specified
#'   directories.
#' @note
#' - The function downloads the most recent version of Global Roads Inventory
#'   Project (`GRIP`) data from the URL specified in the environment variable
#'   `DP_R_Roads_url`. Original data format is a zipped file containing global
#'   road data in the form of `fgdb` (`EPSG:3246`).
#' - On LUMI HPC, loading the `libarchive` module is necessary to use the
#'   `archive` R package: `module load libarchive/3.6.2-cpeGNU-23.09`
#' - The distance to roads is calculated by determining the distance from each
#'   grid
#'   cell to the nearest grid cell that overlaps with a road (not to the nearest
#'   road line). Note that this is different from calculating the actual
#'   distance to the nearest road line, which is computationally intensive and
#'   not performed in this function.
#' @name Road_Intensity
#' @export
#' @author Ahmed El-Gabbas

Road_Intensity <- function(EnvFile = ".env") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "Download")

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Roads <- Path_Roads_Raw <- Path_Roads_Interim <- RoadType <-
    EU_Bound <- VarName <- Road_URL <- Path_Grid <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Roads", "DP_R_Roads_processed", FALSE, FALSE,
    "Path_Roads_Raw", "DP_R_Roads_raw", FALSE, FALSE,
    "Path_Roads_Interim", "DP_R_Roads_interim", FALSE, FALSE,
    "Road_URL", "DP_R_Roads_url", FALSE, FALSE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  fs::dir_create(c(Path_Roads, Path_Roads_Raw, Path_Roads_Interim))

  RefGrid <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(RefGrid)) {
    stop("The reference grid file does not exist: ", RefGrid, call. = FALSE)
  }

  # # ..................................................................... ###

  # Download road data ------
  IASDT.R::CatTime("Download road data")

  withr::local_options(timeout = 1200)
  Path_DownFile <- IASDT.R::Path(Path_Roads_Raw, basename(Road_URL))

  # Check if zip file is a valid file
  if (file.exists(Path_DownFile)) {
    Success <- IASDT.R::CheckZip(Path_DownFile)
    if (isFALSE(Success)) {
      fs::file_delete(Path_DownFile)
    }
  } else {
    Success <- FALSE
  }

  # Try downloading data for a max of 5 attempts
  Attempt <- 1
  Attempts <- 5

  while (isFALSE(Success) && Attempt <= Attempts) {
    Down <- try(
      expr = {
        utils::download.file(
          url = Road_URL, destfile = Path_DownFile,
          mode = "wb", quiet = TRUE) %>%
          suppressWarnings()

        Success <- IASDT.R::CheckZip(Path_DownFile)
        Success
      }, silent = TRUE)

    if (inherits(Down, "try-error")) {
      Success <- FALSE
    }

    Attempt <- Attempt + 1
  }

  if (isFALSE(Success)) {
    stop(
      "Failed to download road data after ", Attempts, " attempts:\n",
      Road_URL, call. = FALSE)
  }

  rm(Down, envir = environment())

  # # .................................... ###

  IASDT.R::CatTime("Extracting files")
  archive::archive_extract(
    archive = Path_DownFile, dir = Path_Roads_Interim) %>%
    suppressMessages()
  rm(Path_DownFile, envir = environment())

  # # ..................................................................... ###

  # Processing GRIP road data ------
  IASDT.R::CatTime("Processing GRIP road data")

  ## Load, crop, and project GRIP data -----
  IASDT.R::CatTime("Load, crop, and project GRIP data", Level = 1)

  Road_GDB_Files <- list.files(
    path = Path_Roads_Interim, pattern = ".gdb$", full.names = TRUE)

  if (length(Road_GDB_Files) == 0) {
    stop(
      "No `.gdb` files found in the directory after extraction: ",
      Path_Roads_Interim, call. = FALSE)
  }

  Road_sf <- Road_GDB_Files[1] %>%
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
  IASDT.R::CatTime("Save projected data - RData", Level = 1)
  save(Road_sf, file = IASDT.R::Path(Path_Roads, "Road_sf.RData"))

  # # ..................................... ###

  ## One file per road type ----
  IASDT.R::CatTime("Save RData file per road type", Level = 1)

  tibble::tribble(
    ~RoadType, ~VarName,
    1, "Highways",
    2, "Primary",
    3, "Secondary",
    4, "Tertiary",
    5, "Local") %>%
    dplyr::mutate(
      A = purrr::walk2(
        .x = RoadType, .y = VarName,
        .f = ~ {
          IASDT.R::CatTime(paste0(.x, " - ", .y), Level = 2)
          dplyr::filter(Road_sf, GP_RTP %in% .x) %>%
            dplyr::select(-GP_RTP) %>%
            terra::vect() %>%
            terra::wrap() %>%
            IASDT.R::SaveAs(
              OutObj = paste0("Road_sf_", .x, "_", .y),
              OutPath = IASDT.R::Path(
                Path_Roads, paste0("Road_sf_", .x, "_", .y, ".RData")))
          invisible(gc())
          return(invisible(NULL))
        }
      )
    ) %>%
    invisible()

  rm(Road_sf, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Road length ------
  IASDT.R::CatTime("Road length")


  IASDT.R::CatTime("Calculate Road length per road type", Level = 1)

  RefGrid <- terra::unwrap(IASDT.R::LoadAs(RefGrid))

  ExtractRoadSummary <- function(RoadType, Function = "length", VarName) {
    SummMap <- list.files(
      path = Path_Roads, full.names = TRUE,
      pattern = paste0("^Road_sf_", RoadType, "_.+RData")) %>%
      IASDT.R::LoadAs() %>%
      terra::unwrap() %>%
      terra::rasterizeGeom(y = RefGrid, fun = Function, unit = "km") %>%
      terra::mask(mask = RefGrid) %>%
      stats::setNames(paste0(RoadType, "_", VarName))
    return(SummMap)
  }

  IASDT.R::CatTime("1 - Highways", Level = 2)
  GRIP_1 <- ExtractRoadSummary(RoadType = 1, VarName = "Highways")

  IASDT.R::CatTime("2 - Primary", Level = 2)
  GRIP_2 <- ExtractRoadSummary(RoadType = 2, VarName = "Primary")

  IASDT.R::CatTime("3 - Secondary", Level = 2)
  GRIP_3 <- ExtractRoadSummary(RoadType = 3, VarName = "Secondary")

  IASDT.R::CatTime("4 - Tertiary", Level = 2)
  GRIP_4 <- ExtractRoadSummary(RoadType = 4, VarName = "Tertiary")

  IASDT.R::CatTime("5 - Local", Level = 2)
  GRIP_5 <- ExtractRoadSummary(RoadType = 5, VarName = "Local")

  IASDT.R::CatTime("All roads", Level = 2)
  Road_Length <- (GRIP_1 + GRIP_2 + GRIP_3 + GRIP_4 + GRIP_5) %>%
    stats::setNames("All") %>%
    c(GRIP_1, GRIP_2, GRIP_3, GRIP_4, GRIP_5, .) %>%
    # Ensure that values are read from memory
    IASDT.R::setRastVals()

  IASDT.R::CatTime("Save road length - tif", Level = 1)
  terra::writeRaster(
    x = Road_Length, overwrite = TRUE,
    filename = IASDT.R::Path(
      Path_Roads, paste0("Road_Length_", names(Road_Length), ".tif")))

  IASDT.R::CatTime("Save road length - RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Road_Length),
    OutObj = "Road_Length",
    OutPath = IASDT.R::Path(Path_Roads, "Road_Length.RData"))

  # # ..................................................................... ###

  # Distance to roads ------
  IASDT.R::CatTime("Distance to roads")
  # This calculates the distance from each grid cell to the nearest grid cell
  # overlapping with a road This can be different than calculating the actual
  # distance to nearest road line; which is expected to take too much time to
  # calculate

  IASDT.R::CatTime("Calculate distance to roads", Level = 1)

  # suppress progress bar
  terra::terraOptions(progress = 0)

  Road_Distance <- purrr::map(
    .x = as.list(Road_Length),
    .f = ~ {
      IASDT.R::CatTime(names(.x), Level = 2)
      Road_Points <- terra::as.points(terra::classify(.x, cbind(0, NA)))
      terra::distance(x = .x, y = Road_Points, unit = "km") %>%
        terra::mask(RefGrid) %>%
        stats::setNames(paste0("Road_Distance_", names(.x))) %>%
        # Ensure that values are read from memory
        IASDT.R::setRastVals()
    }
  ) %>%
    terra::rast()

  IASDT.R::CatTime("Save distance to road - tif", Level = 1)
  terra::writeRaster(
    x = Road_Distance, overwrite = TRUE,
    filename = IASDT.R::Path(Path_Roads, paste0(names(Road_Distance), ".tif")))

  IASDT.R::CatTime("Save distance to road - RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Road_Distance), OutObj = "Road_Distance",
    OutPath = IASDT.R::Path(Path_Roads, "Road_Distance.RData"))

  # # ..................................................................... ###

  # Plotting ------
  IASDT.R::CatTime("Plotting")

  EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  # nolint start
  PlottingTheme <- ggplot2::theme_bw() +
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
  IASDT.R::CatTime("Road length", Level = 1)

  Plots_Length <- purrr::map(
    .x = terra::as.list(Road_Length),
    .f = ~ {
      Road <- log10(terra::classify(.x, cbind(0, NA)))
      Title <- names(.x) %>%
        stringr::str_remove("Road_Distance_") %>%
        stringr::str_replace("_", " - ") %>%
        paste(" roads")

      ggplot2::ggplot() +
        ggplot2::geom_sf(
          EU_Bound,
          mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(data = Road, maxcell = Inf) +
        ggplot2::geom_sf(
          EU_Bound,
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
        PlottingTheme
    }
  )

  Plots_Length <- patchwork::wrap_plots(Plots_Length, ncol = 3, nrow = 2) +
    patchwork::plot_annotation(
      title = "Road length per grid cell",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0.5, 0))))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = IASDT.R::Path(Path_Roads, "Road_Length.jpeg"),
    width = 30, height = 21, units = "cm", quality = 100, res = 600)
  print(Plots_Length)
  grDevices::dev.off()

  rm(Plots_Length, envir = environment())

  # # ..................................... ###

  ## Distance to roads ------
  IASDT.R::CatTime("Distance to roads", Level = 1)
  Plots_Distance <- purrr::map(
    .x = terra::as.list(Road_Distance),
    .f = ~ {
      ggplot2::ggplot() +
        ggplot2::geom_sf(
          EU_Bound,
          mapping = ggplot2::aes(), color = "grey75",
          linewidth = 0.075, fill = "grey98") +
        tidyterra::geom_spatraster(data = .x, maxcell = Inf) +
        ggplot2::geom_sf(
          EU_Bound,
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
        PlottingTheme
    }
  )

  Plots_Distance <- patchwork::wrap_plots(Plots_Distance, ncol = 3, nrow = 2) +
    patchwork::plot_annotation(
      title = "Distance to nearest road",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 1, 0))))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = IASDT.R::Path(Path_Roads, "Road_Distance.jpeg"),
    width = 30, height = 21, units = "cm", quality = 100, res = 600)
  print(Plots_Distance)
  grDevices::dev.off()

  rm(Plots_Distance, envir = environment())

  # ..................................................................... ###

  # Cleanup ------
  IASDT.R::CatTime("Cleanup")

  # Delete extracted GRIP files
  list.files(Path_Roads_Interim, full.names = TRUE, pattern = "^GRIP") %>%
    fs::file_delete()

  fs::dir_delete(c(Path_Roads_Interim, Path_Roads_Raw))

  # ..................................................................... ###

  # Function Summary ----

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing road data was finished in ", ... = "\n")

  return(invisible(NULL))
}
