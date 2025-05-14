#' Process GBIF occurrence data for the `IASDT`
#'
#' Extracts, processes, and visualises occurrence data from the [Global
#' Biodiversity Information Facility (GBIF)](https://www.gbif.org) for the
#' Invasive Alien Species Digital Twin (`IASDT`). Orchestrated by
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
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "future::sequential", "future::multisession",
#'   "future::multicore", and "future::cluster". Defaults to
#'   `"future::multicore"` (`"future::multisession"` on Windows). See
#'   [future::plan()] and [ecokit::set_parallel()] for details.
#' @param delete_chunks Logical. If `TRUE` (default), deletes chunk files.
#' @param chunk_file Character. Path of chunk file for processing.
#'
#' @param max_uncertainty Numeric. Maximum spatial uncertainty in kilometres.
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
    env_file = ".env", r_environ = ".Renviron", n_cores = 6L,
    strategy = "future::multicore", request = TRUE, download = TRUE,
    split_chunks = TRUE, overwrite = FALSE, delete_chunks = TRUE,
    chunk_size = 50000L, boundaries = c(-30, 50, 25, 75), start_year = 1981L) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----

  ecokit::cat_time("Checking arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("env_file", "Renviron", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "request", "download", "split_chunks", "overwrite", "delete_chunks"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("start_year", "boundaries", "chunk_size", "n_cores"))

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "future::sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c(
    "future::sequential", "future::multisession", "future::multicore",
    "future::cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

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
  ecokit::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "Path_GBIF_Interim", "DP_R_GBIF_interim", FALSE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "CountryCodes", "DP_R_Countrycodes", FALSE, TRUE,
    "TaxaInfo", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())
  # # ..................................................................... ###

  # request / download GBIF data and split data into chunks -----

  IASDT.R::GBIF_download(
    env_file = env_file, r_environ = r_environ, request = request,
    download = download, split_chunks = split_chunks, chunk_size = chunk_size,
    boundaries = boundaries, start_year = start_year)

  # # ..................................................................... ###

  GBIF_Metadata <- fs::path(Path_GBIF, "GBIF_Metadata.RData")
  if (!file.exists(GBIF_Metadata)) {
    ecokit::stop_ctx(
      "GBIF metadata file does not exist", GBIF_Metadata = GBIF_Metadata,
      include_backtrace = TRUE)
  }
  GBIF_Metadata <- ecokit::load_as(GBIF_Metadata)

  TaxaList <- ecokit::load_as(TaxaInfo)
  Path_SpData <- fs::path(Path_GBIF, "Sp_Data")
  fs::dir_create(Path_SpData)

  # Grid_10_Land_Crop_sf
  GridSf <- fs::path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    ecokit::stop_ctx(
      "Reference grid (sf) file not found", GridSf = GridSf,
      include_backtrace = TRUE)
  }
  GridSf <- ecokit::load_as(GridSf)

  # Grid_10_Land_Crop
  GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    ecokit::stop_ctx(
      "Reference grid file not found", GridR = GridR, include_backtrace = TRUE)
  }
  GridR <- terra::unwrap(ecokit::load_as(GridR))

  EuroBound <- ecokit::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  invisible(gc())

  # # ..................................................................... ###

  # Processing data chunks -----

  ecokit::cat_time("\nProcessing data chunks")
  .StartTimeChunks <- lubridate::now(tzone = "CET")

  ChunkList <- list.files(
    path = Path_GBIF_Interim, pattern = "Chunk_.+.txt", full.names = TRUE)
  ChunkListRData <- stringr::str_replace_all(ChunkList, ".txt$", ".RData")

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  ecokit::cat_time(
    "Processing chunks in parallel, save each as RData files", level = 1L)

  if (strategy == "future::multicore") {
    pkg_to_export <- NULL
  } else {
    pkg_to_export <- c(
      "IASDT.R", "purrr", "tibble", "terra", "tidyr", "dplyr", "readr",
      "stringr", "bit64", "tidyselect", "fs", "sf", "ecokit")
  }

  GBIF_Data <- future.apply::future_lapply(
    X = ChunkList,
    FUN = function(x) {
      IASDT.R::GBIF_read_chunk(
        chunk_file = x, env_file = env_file, save_RData = TRUE,
        return_data = FALSE, overwrite = overwrite)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "overwrite"))

  ecokit::cat_diff(
    init_time = .StartTimeChunks, prefix = "Finished in ", level = 2L)

  ecokit::cat_time("Reading processed chunks into a single dataset", level = 1L)

  if (all(file.exists(ChunkListRData))) {
    GBIF_Data <- future.apply::future_lapply(
      X = ChunkListRData, FUN = ecokit::load_as,
      future.scheduling = Inf, future.seed = TRUE) %>%
      dplyr::bind_rows() %>%
      # merge with taxa standardization results
      dplyr::left_join(TaxaList, by = "speciesKey") %>%
      # arrange by species name
      dplyr::arrange(Species_name)
  } else {
    Msg <- ChunkList[which(!file.exists(ChunkListRData))] %>%
      paste0(" >> ", ., collapse = "\n") %>%
      paste0("The following chunks were not processed\n", .)
    ecokit::stop_ctx(Msg, include_backtrace = TRUE)
  }

  ecokit::cat_diff(
    init_time = .StartTimeChunks,
    prefix = "Processing data chunks was finished in ", level = 2L)

  ecokit::cat_time(
    paste0(
      "A total of ", format(nrow(GBIF_Data), big.mark = ","), " observations"),
    level = 2L, cat_timestamp = FALSE)


  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("future::sequential", gc = TRUE)
  }

  ecokit::cat_time("Saving `GBIF_Data` to disk", level = 1L)
  ecokit::save_as(
    object = GBIF_Data, out_path = fs::path(Path_GBIF, "GBIF_Data.qs2"))

  rm(TaxaList, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # iNaturalist unique grids ----

  ecokit::cat_time("iNaturalist unique grids")

  ecokit::cat_time("Prepare input data for summarising", level = 1L)
  iNaturalist_Others <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(species, CellCode, institutionCode) %>%
    dplyr::mutate(
      institutionCode = tidyr::replace_na(institutionCode, ""),
      institutionCode = institutionCode == "iNaturalist",
      institutionCode = dplyr::if_else(
        institutionCode, "iNaturalist", "Others"))

  ecokit::cat_time(
    "Number of unique iNaturalist grid cells per species", level = 1L)
  iNaturalist_Unique <- iNaturalist_Others %>%
    dplyr::group_by(species, CellCode) %>%
    tidyr::pivot_wider(
      names_from = institutionCode, values_from = institutionCode,
      values_fn = unique) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(Others), iNaturalist == "iNaturalist") %>%
    dplyr::count(species) %>%
    dplyr::rename(iNaturalist_Unique = n)

  ecokit::cat_time(
    "Number of grid cells for iNaturalist and other data sources",
    level = 1L)

  iNaturalist_Count <- iNaturalist_Others %>%
    dplyr::count(species, institutionCode) %>%
    tidyr::pivot_wider(names_from = institutionCode, values_from = n) %>%
    dplyr::select(species, iNaturalist, Others) %>%
    dplyr::left_join(iNaturalist_Unique, by = "species") %>%
    dplyr::select(species, Others, tidyselect::everything())

  ecokit::cat_time("Save iNaturalist summary data", level = 1L)
  save(
    iNaturalist_Count,
    file = fs::path(Path_GBIF, "iNaturalist_Count.RData"))

  rm(
    iNaturalist_Others, iNaturalist_Unique, iNaturalist_Count,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Covert GBIF data to presence-only grids ----

  ecokit::cat_time("Covert GBIF data to presence-only grids")
  GBIF_Grid <- sf::st_drop_geometry(GBIF_Data) %>%
    # only unique combinations between species and grid ID
    dplyr::distinct(
      IAS_ID, taxon_name, Species_name, Species_name2,
      Species_File, CellCode) %>%
    # add polygons for grids
    dplyr::left_join(GridSf, by = "CellCode") %>%
    # ensure that the data is sf object
    sf::st_as_sf()

  ecokit::cat_time("Save presence-only grids", level = 1L)
  save(GBIF_Grid, file = fs::path(Path_GBIF, "GBIF_Grid.RData"))
  invisible(gc())

  # # ..................................................................... ###

  # Number of observations or grid cells per species ----

  ecokit::cat_time("# observations or grid cells per species")

  ecokit::cat_time("# observations per species", level = 1L)
  SpNObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumOcc = n)

  ecokit::cat_time("# grids per species", level = 1L)
  SpNGrids <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumGrids = n)

  ecokit::cat_time(
    "Merge number of observations and grids per species", level = 1L)
  GBIF_NObsNGrid <- dplyr::full_join(
    SpNObs, SpNGrids,
    by = c("IAS_ID", "taxon_name", "Species_name"))

  ecokit::cat_time("Save as RData", level = 1L)
  save(GBIF_NObsNGrid, file = fs::path(Path_GBIF, "GBIF_NObsNGrid.RData"))

  ecokit::cat_time("Save as xlsx", level = 1L)
  writexl::write_xlsx(
    x = GBIF_NObsNGrid, path = fs::path(Path_GBIF, "GBIF_NObsNGrid.xlsx"))

  rm(GBIF_NObsNGrid, SpNObs, SpNGrids, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Summary maps ----
  ecokit::cat_time("Summary maps")

  ## Number of observations per grid  ----

  ecokit::cat_time("Number of observations per grid cell", level = 1L)
  GBIF_NObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumObs") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time(
    "Number of observations per grid cell (log scale)", level = 1L)
  GBIF_NObs_log <- terra::unwrap(GBIF_NObs) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Save as RData", level = 2L)
  save(GBIF_NObs, GBIF_NObs_log, file = fs::path(Path_GBIF, "GBIF_NObs.RData"))

  invisible(gc())

  # # ................................... ###

  ## Number of species per grid cell -----

  ecokit::cat_time("Number of IAS per grid cell", level = 1L)
  GBIF_NSp <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumSpecies") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Number of IAS per grid cell (log)", level = 1L)
  GBIF_NSp_Log <- terra::unwrap(GBIF_NSp) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Save as RData", level = 2L)
  save(GBIF_NSp, GBIF_NSp_Log, file = fs::path(Path_GBIF, "GBIF_NSp.RData"))

  rm(GridSf, GridR, envir = environment())
  invisible(gc())

  # # ................................... ###

  ## Plot_GBIF_Summary -----
  ecokit::cat_time("Plotting summary maps", level = 1L)

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
        breaks = ecokit::integer_breaks()) +
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
  ragg::agg_jpeg(
    filename = fs::path(Path_GBIF, "GBIF_Summary.jpeg"),
    width = 25, height = 25.8, res = 600, quality = 100, units = "cm")
  print(Plot)
  grDevices::dev.off()

  rm(EuroBound, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Species-specific data ----
  ecokit::cat_time("Species-specific data")

  SpList <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name) %>%
    sort()

  # Species data --- sf ----
  ecokit::cat_time("Split species data - sf", level = 1L)

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
      SpName <- ecokit::replace_space(.x) %>%
        # replace non-ascii multiplication symbol with x
        stringr::str_replace_all("\u00D7", "x") %>%
        stringr::str_replace_all("-", "")

      SpData <- dplyr::filter(GBIF_Data, Species_name == .x)
      OutFileSF <- fs::path(Path_SpData, paste0(SpName, ".RData"))

      # Save if there is data
      if (nrow(SpData) > 0) {
        ecokit::save_as(
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
  ecokit::cat_time("Split species data - grid + raster + plot", level = 1L)

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 2L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  ecokit::cat_time("Splitting species data in parallel", level = 2L)

  if (strategy == "future::multicore") {
    pkg_to_export <- NULL
  } else {
    pkg_to_export <- c(
      "dplyr", "ecokit", "cowplot", "ggplot2", "terra", "tidyterra", "fs",
      "tibble", "magrittr", "stringr", "sf", "ragg", "cowplot", "grid")
  }

  sp_data_tmp <- future.apply::future_lapply(
    X = SpList,
    FUN = IASDT.R::GBIF_species_data, env_file = env_file,
    verbose = FALSE, plot_tag = plot_tag,
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export, future.globals = c("env_file", "plot_tag"))
  rm(sp_data_tmp, envir = environment())

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Clean up chunk files ----
  if (delete_chunks) {
    ecokit::cat_time("Clean up - remove temporary chunk files")
    list.files(Path_GBIF_Interim, full.names = TRUE) %>%
      fs::file_delete()
    fs::dir_delete(Path_GBIF_Interim)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nProcessing GBIF data was finished in ")

  return(invisible(NULL))
}
