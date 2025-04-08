#' Process GBIF occurrence data for the `IAS-pDT`
#'
#' Extracts, processes, and visualizes occurrence data from the [Global
#' Biodiversity Information Facility (GBIF)](https://www.gbif.org) for the
#' Invasive Alien Species prototype Digital Twin (`IAS-pDT`). Orchestrated by
#' `GBIF_process()`, it requests, downloads, cleans, chunks, and maps species
#' data using helper functions.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param r_environ Character. Path to `.Renviron` file with GBIF credentials
#'   (`GBIF_EMAIL`, `GBIF_USER`, `GBIF_PWD`). Default: `".Renviron"`. The
#'   credentials must be in the format:
#'    - `GBIF_EMAIL=your_email`
#'    - `GBIF_USER=your_username`
#'    - `GBIF_PWD=your_password`
#' @param request Logical. If `TRUE` (default), requests GBIF data; otherwise,
#'   loads from disk.
#' @param download Logical. If `TRUE` (default), downloads and saves GBIF data.
#' @param split_chunks Logical. If `TRUE` (default), splits data into chunks for
#'   easier processing.
#' @param chunk_size Integer. Records per data chunk. Default: `50000`.
#' @param boundaries Numeric vector (length 4). GBIF data bounds (Left, Right,
#'   Bottom, Top). Default: `c(-30, 50, 25, 75)`.
#' @param start_year Integer. Earliest year for GBIF records (matches CHELSA
#'   climate data). Default: `1981`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6.
#' @param delete_chunks Logical. If `TRUE` (default), deletes chunk files.
#' @param chunk_file Character. Path of chunk file for processing.
#'
#' @param max_uncertainty Numeric. Maximum spatial uncertainty in kilometers.
#'   Default: `10`.
#' @param start_year Integer. Earliest collection year to be included. Default
#'   is 1981.
#' @param save_RData Logical. If `TRUE` (default), saves chunk data as `.RData`.
#' @param return_data If `TRUE`, returns chunk data; otherwise,
#'   `invisible(NULL)`. Default: `FALSE`.
#' @param overwrite Logical. If `TRUE`, reprocesses existing `.RData` chunks.
#'   Default: `FALSE`. This helps to continue working on previously processed
#'   chunks if the previous try failed, e.g. due to memory issue.
#' @param species Character. Species name for processing.
#' @param verbose Logical. If `TRUE` (default), prints progress messages.
#' @param plot_tag Character. Tag for plot titles.
#'
#' @section Functions details:
#'
#' - **`GBIF_process()`**: Orchestrates GBIF data requests, downloads,
#'   processing, and mapping. Saves `RData`, Excel, and JPEG summary files.
#' - **`GBIF_check()`**: Verifies GBIF credentials in environment or
#'   `.Renviron`. Returns `TRUE` if valid, else `FALSE`.
#' - **`GBIF_download()`**: Requests and downloads GBIF data (if `download =
#'   TRUE`),  using the specified criteria (taxa, coordinates, time period, and
#'   boundaries), splits into small chunks (if `split_chunks = TRUE`), and saves
#'   metadata. Returns `invisible(NULL)`.
#' - **`GBIF_read_chunk()`**: Filters chunk data (spatial/temporal, e.g.,
#'   spatial uncertainty, collection year, coordinate precision, and taxonomic
#'   rank), select relevant columns, and saves as `.RData` (if `save_RData =
#'   TRUE`) or returns it (if `return_data = TRUE`). Skips if `.RData` exists
#'   and `overwrite = FALSE`.
#' - **`GBIF_species_data()`**: Converts species-specific data to `sf` and
#'   raster formats, generating distribution maps.
#' @note Relies on a static RDS file listing IAS species, GBIF keys, and
#'   metadata, standardized by Marina Golivets (Feb 2024).

# # |------------------------------------------------------------------------| #
# GBIF_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 1

GBIF_process <- function(
    env_file = ".env", r_environ = ".Renviron", n_cores = 6L, request = TRUE,
    download = TRUE, split_chunks = TRUE, overwrite = FALSE,
    delete_chunks = TRUE, chunk_size = 50000L, boundaries = c(-30, 50, 25, 75),
    start_year = 1981L) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----

  IASDT.R::cat_time("Checking arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("env_file", "Renviron"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "request", "download", "split_chunks", "overwrite", "delete_chunks"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("start_year", "boundaries", "chunk_size", "n_cores"))

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
    future.seed = TRUE)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_GBIF <- Path_GBIF_Interim <- Species_name <- institutionCode <-
    TaxaInfo <- species <- CellCode <- iNaturalist <- IAS_ID <- taxon_name <-
    Species_name2 <- Species_File <- Others <- Path_Grid <- SummMap <-
    SummMapGG <- Title <- LegendLabel <- n <- EU_Bound <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "Path_GBIF_Interim", "DP_R_GBIF_interim", FALSE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "CountryCodes", "DP_R_Countrycodes", FALSE, TRUE,
    "TaxaInfo", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())
  # # ..................................................................... ###

  # request / download GBIF data and split data into chunks -----

  IASDT.R::GBIF_download(
    env_file = env_file, r_environ = r_environ, request = request,
    download = download, split_chunks = split_chunks, chunk_size = chunk_size,
    boundaries = boundaries, start_year = start_year)

  # # ..................................................................... ###

  GBIF_Metadata <- IASDT.R::path(Path_GBIF, "GBIF_Metadata.RData")
  if (!file.exists(GBIF_Metadata)) {
    stop("GBIF metadata file does not exist: ", GBIF_Metadata, call. = FALSE)
  }
  GBIF_Metadata <- IASDT.R::load_as(GBIF_Metadata)

  TaxaList <- IASDT.R::load_as(TaxaInfo)
  Path_SpData <- IASDT.R::path(Path_GBIF, "Sp_Data")
  fs::dir_create(Path_SpData)

  # Grid_10_Land_Crop_sf
  GridSf <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Reference grid (sf) file not found at: ", GridSf, call. = FALSE)
  }
  GridSf <- IASDT.R::load_as(GridSf)

  # Grid_10_Land_Crop
  GridR <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop("Reference grid file not found at: ", GridR, call. = FALSE)
  }
  GridR <- terra::unwrap(IASDT.R::load_as(GridR))

  EuroBound <- IASDT.R::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  invisible(gc())

  # # ..................................................................... ###

  # Processing data chunks -----

  IASDT.R::cat_time("Processing data chunks")
  .StartTimeChunks <- lubridate::now(tzone = "CET")

  ChunkList <- list.files(
    path = Path_GBIF_Interim, pattern = "Chunk_.+.txt", full.names = TRUE)
  ChunkListRData <- stringr::str_replace_all(ChunkList, ".txt$", ".RData")

  IASDT.R::cat_time(
    paste0("Prepare working on parallel using ", n_cores, " cores."),
    level = 1)

  if (n_cores == 1) {
    terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(n_cores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)

    # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
    # official parameters (overriding the ones from GeoTIFF keys)
    # see: https://stackoverflow.com/questions/78007307
    snow::clusterEvalQ(
      cl = c1, expr = terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")) %>%
      invisible()

    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  IASDT.R::cat_time(
    "Processing chunks on parallel, save each as RData files",
    level = 1)

  GBIF_Data <- future.apply::future_lapply(
    X = ChunkList,
    FUN = function(x) {
      IASDT.R::GBIF_read_chunk(
        chunk_file = x, env_file = env_file, save_RData = TRUE,
        return_data = FALSE, overwrite = overwrite)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = c(
      "IASDT.R", "purrr", "tibble", "terra", "tidyr", "dplyr", "sf"),
    future.globals = c("env_file", "overwrite"))

  IASDT.R::cat_diff(
    init_time = .StartTimeChunks, prefix = "Finished in ", level = 2)

  IASDT.R::cat_time("Reading processed chunks into a single dataset", level = 1)
  if (all(file.exists(ChunkListRData))) {
    GBIF_Data <- future.apply::future_lapply(
      X = ChunkListRData, FUN = IASDT.R::load_as,
      future.scheduling = Inf, future.seed = TRUE) %>%
      dplyr::bind_rows() %>%
      # merge with taxa standardization results
      dplyr::left_join(TaxaList, by = "speciesKey") %>%
      # arrange by species name
      dplyr::arrange(Species_name)
  } else {
    ChunkList[which(!file.exists(ChunkListRData))] %>%
      paste0(" >> ", ., collapse = "\n") %>%
      paste0("The following chunks were not processed\n", .) %>%
      stop(call. = FALSE)
  }

  IASDT.R::cat_diff(
    init_time = .StartTimeChunks,
    prefix = "Processing data chunks was finished in ", level = 2)

  IASDT.R::cat_time(
    paste0(
      "A total of ", format(nrow(GBIF_Data), big.mark = ","), " observations"),
    level = 2)

  IASDT.R::cat_time("Stopping cluster", level = 1)
  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::cat_time("Saving `GBIF_Data` to disk", level = 1)
  save(GBIF_Data, file = IASDT.R::path(Path_GBIF, "GBIF_Data.RData"))

  rm(TaxaList, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # iNaturalist unique grids ----

  IASDT.R::cat_time("iNaturalist unique grids")

  IASDT.R::cat_time("Prepare input data for summarizing", level = 1)
  iNaturalist_Others <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(species, CellCode, institutionCode) %>%
    dplyr::mutate(
      institutionCode = tidyr::replace_na(institutionCode, ""),
      institutionCode = institutionCode == "iNaturalist",
      institutionCode = dplyr::if_else(
        institutionCode, "iNaturalist", "Others"))

  IASDT.R::cat_time(
    "Number of unique iNaturalist grid cells per species", level = 1)
  iNaturalist_Unique <- iNaturalist_Others %>%
    dplyr::group_by(species, CellCode) %>%
    tidyr::pivot_wider(
      names_from = institutionCode, values_from = institutionCode,
      values_fn = unique) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(Others), iNaturalist == "iNaturalist") %>%
    dplyr::count(species) %>%
    dplyr::rename(iNaturalist_Unique = n)

  IASDT.R::cat_time(
    "Number of grid cells for iNaturalist and other data sources",
    level = 1)

  iNaturalist_Count <- iNaturalist_Others %>%
    dplyr::count(species, institutionCode) %>%
    tidyr::pivot_wider(names_from = institutionCode, values_from = n) %>%
    dplyr::select(species, iNaturalist, Others) %>%
    dplyr::left_join(iNaturalist_Unique, by = "species") %>%
    dplyr::select(species, Others, tidyselect::everything())

  IASDT.R::cat_time("Save iNaturalist summary data", level = 1)
  save(
    iNaturalist_Count,
    file = IASDT.R::path(Path_GBIF, "iNaturalist_Count.RData"))

  rm(
    iNaturalist_Others, iNaturalist_Unique, iNaturalist_Count,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Covert GBIF data to presence-only grids ----

  IASDT.R::cat_time("Covert GBIF data to presence-only grids")
  GBIF_Grid <- sf::st_drop_geometry(GBIF_Data) %>%
    # only unique combinations between species and grid ID
    dplyr::distinct(
      IAS_ID, taxon_name, Species_name, Species_name2,
      Species_File, CellCode) %>%
    # add polygons for grids
    dplyr::left_join(GridSf, by = "CellCode") %>%
    # ensure that the data is sf object
    sf::st_as_sf()

  IASDT.R::cat_time("Save presence-only grids", level = 1)
  save(GBIF_Grid, file = IASDT.R::path(Path_GBIF, "GBIF_Grid.RData"))
  invisible(gc())

  # # ..................................................................... ###

  # Number of observations or grid cells per species ----

  IASDT.R::cat_time("# observations or grid cells per species")

  IASDT.R::cat_time("# observations per species", level = 1)
  SpNObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumOcc = n)

  IASDT.R::cat_time("# grids per species", level = 1)
  SpNGrids <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumGrids = n)

  IASDT.R::cat_time(
    "Merge number of observations and grids per species", level = 1)
  GBIF_NObsNGrid <- dplyr::full_join(
    SpNObs, SpNGrids,
    by = c("IAS_ID", "taxon_name", "Species_name"))

  IASDT.R::cat_time("Save as RData", level = 1)
  save(GBIF_NObsNGrid, file = IASDT.R::path(Path_GBIF, "GBIF_NObsNGrid.RData"))

  IASDT.R::cat_time("Save as xlsx", level = 1)
  writexl::write_xlsx(
    x = GBIF_NObsNGrid, path = IASDT.R::path(Path_GBIF, "GBIF_NObsNGrid.xlsx"))

  rm(GBIF_NObsNGrid, SpNObs, SpNGrids, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Summary maps ----
  IASDT.R::cat_time("Summary maps")

  ## Number of observations per grid  ----

  IASDT.R::cat_time("Number of observations per grid cell", level = 1)
  GBIF_NObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumObs") %>%
    IASDT.R::set_raster_CRS() %>%
    terra::wrap()

  IASDT.R::cat_time(
    "Number of observations per grid cell (log scale)", level = 1)
  GBIF_NObs_log <- terra::unwrap(GBIF_NObs) %>%
    log10() %>%
    IASDT.R::set_raster_CRS() %>%
    terra::wrap()

  IASDT.R::cat_time("Save as RData", level = 2)
  save(
    GBIF_NObs, GBIF_NObs_log,
    file = IASDT.R::path(Path_GBIF, "GBIF_NObs.RData"))

  invisible(gc())

  # # ................................... ###

  ## Number of species per grid cell -----

  IASDT.R::cat_time("Number of IAS per grid cell", level = 1)
  GBIF_NSp <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumSpecies") %>%
    IASDT.R::set_raster_CRS() %>%
    terra::wrap()

  IASDT.R::cat_time("Number of IAS per grid cell (log)", level = 1)
  GBIF_NSp_Log <- terra::unwrap(GBIF_NSp) %>%
    log10() %>%
    IASDT.R::set_raster_CRS() %>%
    terra::wrap()

  IASDT.R::cat_time("Save as RData", level = 2)
  save(
    GBIF_NSp, GBIF_NSp_Log,
    file = IASDT.R::path(Path_GBIF, "GBIF_NSp.RData"))

  rm(GridSf, GridR, envir = environment())
  invisible(gc())

  # # ................................... ###

  ## Plot_GBIF_Summary -----
  IASDT.R::cat_time("Plotting summary maps", level = 1)

  GBIF_date <- GBIF_Metadata$StatusDetailed$modified %>%
    lubridate::date() %>%
    format("%d %B %Y")

  GBIF_DOI <- GBIF_Metadata$StatusDetailed$doi
  plot_tag <- format(Sys.Date(), "%d %B %Y")
  plot_tag <- stringr::str_glue(
    "DOI: {GBIF_DOI} ({GBIF_date})\nLast update: {plot_tag}")

  rm(GBIF_date, GBIF_DOI, GBIF_Grid, GBIF_Metadata, envir = environment())

  Plot_GBIF_Summary <- function(
    RstrMap, Title, LegendLabel = NULL, EU_map = EuroBound) {

    # Plotting limits
    Xlim <- c(2600000, 6700000)
    Ylim <- c(1450000, 5420000)

    PlottingTheme <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0.1, 0, 0.1, "cm"),
        plot.title = ggplot2::element_text(
          size = 12, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0)),
        strip.text = ggplot2::element_text(size = 6, face = "bold"),
        legend.key.size = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.92, 0.75),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 8),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.title = ggplot2::element_text(
          color = "blue", size = 7, face = "bold", hjust = 0.5),
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        panel.spacing = grid::unit(0.3, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.1, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        panel.ontop = TRUE, panel.background = ggplot2::element_rect(fill = NA))

    OutMap <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        EU_map, mapping = ggplot2::aes(), color = "grey30", linewidth = 0.1,
        fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(data = terra::trim(RstrMap), maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = IASDT.R::integer_breaks()) +
      ggplot2::geom_sf(
        EU_map, mapping = ggplot2::aes(), color = "grey40", linewidth = 0.075,
        fill = "transparent", inherit.aes = TRUE) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
      PlottingTheme +
      ggplot2::labs(title = Title, fill = LegendLabel)

    return(OutMap)
  }

  # a common main title of the figure
  MainTitle <- cowplot::ggdraw() +
    cowplot::draw_label(
      label = "Summary of IAS data (GBIF)", fontface = "bold", hjust = 0.5) +
    cowplot::draw_label(
      label = plot_tag, fontface = "italic", color = "grey", x = 0.98,
      y = 0.35, size = 8, hjust = 1) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0.4, 0))

  # Plotting summary maps
  Plot <- tibble::tibble(
    # tibble for plotting information
    SummMap = list(
      terra::unwrap(GBIF_NObs) / 1000, terra::unwrap(GBIF_NObs_log),
      terra::unwrap(GBIF_NSp) / 100, terra::unwrap(GBIF_NSp_Log)),
    Title = c(
      "Number of observations", "Number of observations (log10)",
      "Number of IAS", "Number of IAS (log10)"),
    LegendLabel = c("\u00D7 1000", "log10", "\u00D7 100", "Log10")) %>%
    # add the ggplot object as column to the data
    dplyr::mutate(
      SummMapGG = purrr::pmap(
        .l = list(SummMap, Title, LegendLabel),
        .f = function(SummMap, Title, LegendLabel) {
          Plot_GBIF_Summary(SummMap, Title, LegendLabel)
        }
      )
    ) %>%
    dplyr::pull(SummMapGG) %>%
    # Plot the four panels together on a single figure
    cowplot::plot_grid(plotlist = ., ncol = 2, nrow = 2) %>%
    # add the common title
    cowplot::plot_grid(MainTitle, ., ncol = 1, rel_heights = c(0.035, 1))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = IASDT.R::path(Path_GBIF, "GBIF_Summary.jpeg"),
    width = 25, height = 25.8, units = "cm", quality = 100, res = 600)
  print(Plot)
  grDevices::dev.off()

  rm(EuroBound, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Species-specific data ----
  IASDT.R::cat_time("Species-specific data")

  SpList <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name) %>%
    sort()

  # Species data --- sf ----
  IASDT.R::cat_time("Split species data - sf", level = 1)

  # On LUMI, extracting species data in parallel using `furrr::future_walk` took
  # much longer time than working sequentially `purrr::walk` due to the
  # existence of the very large `GBIF_Data` object. Although it is possible to
  # integrate following chunk in the `GBIF_species_data` function, it would be
  # much faster to implement the following sequentially, then process species
  # maps in parallel later on, after the deletion of the `GBIF_Data` object
  purrr::walk(
    .x = SpList,
    .f = ~ {
      # Filter data on this species
      SpName <- IASDT.R::replace_space(.x) %>%
        # replace non-ascii multiplication symbol with x
        stringr::str_replace_all("\u00D7", "x") %>%
        stringr::str_replace_all("-", "")

      SpData <- dplyr::filter(GBIF_Data, Species_name == .x)
      OutFileSF <- IASDT.R::path(Path_SpData, paste0(SpName, ".RData"))

      # Save if there is data
      if (nrow(SpData) > 0) {
        IASDT.R::save_as(
          object = SpData, object_name = SpName, out_path = OutFileSF)
      }
      return(invisible(NULL))
    },
    .progress = FALSE)

  rm(
    GBIF_Data, GBIF_NObs, GBIF_NObs_log, GBIF_NSp, GBIF_NSp_Log,
    envir = environment())
  invisible(gc())


  # Grid / raster / plotting ----
  IASDT.R::cat_time("Split species data - grid/raster/plot", level = 1)

  IASDT.R::cat_time(
    paste0("Prepare working on parallel using ", n_cores, " cores."),
    level = 1)

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(n_cores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  IASDT.R::cat_time("Splitting species data on parallel", level = 2)
  furrr::future_walk(
    .x = SpList, .f = IASDT.R::GBIF_species_data, env_file = env_file,
    verbose = FALSE, plot_tag = plot_tag,
    .options = furrr::furrr_options(seed = TRUE, packages = "dplyr")
  )

  IASDT.R::cat_time("Stopping cluster", level = 2)
  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Clean up chunk files ----
  if (delete_chunks) {
    IASDT.R::cat_time("Clean up - remove temporary chunk files")
    list.files(Path_GBIF_Interim, full.names = TRUE) %>%
      fs::file_delete()
    fs::dir_delete(Path_GBIF_Interim)
  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nProcessing GBIF data was finished in ")

  return(invisible(NULL))
}
