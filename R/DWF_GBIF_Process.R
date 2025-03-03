#' Process GBIF occurrence data for the `IAS-pDT`
#'
#' Extracts, processes, and visualizes occurrence data from the [Global
#' Biodiversity Information Facility (GBIF)](https://www.gbif.org) for the
#' Invasive Alien Species prototype Digital Twin (`IAS-pDT`). Orchestrated by
#' `GBIF_Process()`, it requests, downloads, cleans, chunks, and maps species
#' data using helper functions.
#'
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Renviron Character. The path to the `.Renviron` file containing GBIF
#'   login credentials. Defaults to `.Renviron` in the current working
#'   directory.
#' @param Renviron Character. Path to `.Renviron` file with GBIF credentials
#'   (`GBIF_EMAIL`, `GBIF_USER`, `GBIF_PWD`). Default: `".Renviron"`. The
#'   credentials must be in the format:
#'    - `GBIF_EMAIL=your_email`
#'    - `GBIF_USER=your_username`
#'    - `GBIF_PWD=your_password`
#' @param Request Logical. If `TRUE` (default), requests GBIF data; otherwise,
#'   loads from disk.
#' @param Download Logical. If `TRUE` (default), downloads and saves GBIF data.
#' @param SplitChunks Logical. If `TRUE` (default), splits data into chunks for
#'   easier processing.
#' @param ChunkSize Integer. Records per data chunk. Default: `50000`.
#' @param Boundaries Numeric vector (length 4). GBIF data bounds (Left, Right,
#'   Bottom, Top). Default: `c(-30, 50, 25, 75)`.
#' @param StartYear Integer. Earliest year for GBIF records (matches CHELSA
#'   climate data). Default: `1981`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6.
#' @param DeleteChunks Logical. If `TRUE` (default), deletes chunk files.
#' @param ChunkFile Character. Path of chunk file for processing.
#'
#' @param MaxUncert Numeric. Maximum spatial uncertainty in kilometers. Default:
#'   `10`.
#' @param StartYear Integer. Earliest collection year to be included. Default is
#'   1981.
#' @param SaveRData Logical. If `TRUE` (default), saves chunk data as `.RData`.
#' @param ReturnData If `TRUE`, returns chunk data; otherwise,
#'   `invisible(NULL)`. Default: `FALSE`.
#' @param Overwrite Logical. If `TRUE`, reprocesses existing `.RData` chunks.
#'   Default: `FALSE`. This helps to continue working on previously processed
#'   chunks if the previous try failed, e.g. due to memory issue.
#' @param Species Character. Species name for processing.
#' @param Verbose Logical. If `TRUE` (default), prints progress messages.
#' @param PlotTag Character. Tag for plot titles.
#'
#' @section Functions details:
#'
#' - **`GBIF_Process()`**: Orchestrates GBIF data requests, downloads,
#'   processing, and mapping. Saves `RData`, Excel, and JPEG summary files.
#' - **`GBIF_Check()`**: Verifies GBIF credentials in environment or
#'   `.Renviron`. Returns `TRUE` if valid, else `FALSE`.
#' - **`GBIF_Download()`**: Requests and downloads GBIF data (if `Download =
#'   TRUE`),  using the specified criteria (taxa, coordinates, time period, and
#'   boundaries), splits into small chunks (if `SplitChunks = TRUE`), and saves
#'   metadata. Returns `invisible(NULL)`.
#' - **`GBIF_ReadChunk()`**: Filters chunk data (spatial/temporal, e.g.,
#'   spatial uncertainty, collection year, coordinate precision, and taxonomic
#'   rank), select relevant columns, and saves as `.RData` (if `SaveRData =
#'   TRUE`) or returns it (if `ReturnData = TRUE`). Skips if `.RData` exists and
#'   `Overwrite = FALSE`.
#' - **`GBIF_SpData()`**: Converts species-specific data to `sf` and raster
#'   formats, generating distribution maps.
#' @note Relies on a static RDS file listing IAS species, GBIF keys, and
#'   metadata, standardized by Marina Golivets (Feb 2024).

# # |------------------------------------------------------------------------| #
# GBIF_Process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 1

GBIF_Process <- function(
    EnvFile = ".env", Renviron = ".Renviron", NCores = 6L, Request = TRUE,
    Download = TRUE, SplitChunks = TRUE, Overwrite = FALSE,
    DeleteChunks = TRUE, ChunkSize = 50000L, Boundaries = c(-30, 50, 25, 75),
    StartYear = 1981L) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----

  IASDT.R::CatTime("Checking arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("EnvFile", "Renviron1"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("Request", "Download", "SplitChunks", "Overwrite", "DeleteChunks"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("StartYear", "Boundaries", "ChunkSize", "NCores"))

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
  IASDT.R::CatTime("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "Path_GBIF_Interim", "DP_R_GBIF_interim", FALSE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "CountryCodes", "DP_R_Countrycodes", FALSE, TRUE,
    "TaxaInfo", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())
  # # ..................................................................... ###

  # Request / download GBIF data and split data into chunks -----

  IASDT.R::GBIF_Download(
    EnvFile = EnvFile, Renviron = Renviron, Request = Request,
    Download = Download, SplitChunks = SplitChunks, ChunkSize = ChunkSize,
    Boundaries = Boundaries, StartYear = StartYear)

  # # ..................................................................... ###

  GBIF_Metadata <- IASDT.R::Path(Path_GBIF, "GBIF_Metadata.RData")
  if (!file.exists(GBIF_Metadata)) {
    stop("GBIF metadata file does not exist: ", GBIF_Metadata, call. = FALSE)
  }
  GBIF_Metadata <- IASDT.R::LoadAs(GBIF_Metadata)

  TaxaList <- IASDT.R::LoadAs(TaxaInfo)
  Path_SpData <- IASDT.R::Path(Path_GBIF, "Sp_Data")
  fs::dir_create(Path_SpData)

  # Grid_10_Land_Crop_sf
  GridSf <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Reference grid (sf) file not found at: ", GridSf, call. = FALSE)
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

  invisible(gc())

  # # ..................................................................... ###

  # Processing data chunks -----

  IASDT.R::CatTime("Processing data chunks")
  .StartTimeChunks <- lubridate::now(tzone = "CET")

  ChunkList <- list.files(
    path = Path_GBIF_Interim, pattern = "Chunk_.+.txt", full.names = TRUE)
  ChunkListRData <- stringr::str_replace_all(ChunkList, ".txt$", ".RData")

  IASDT.R::CatTime(
    paste0("Prepare working on parallel using ", NCores, " cores."),
    Level = 1)

  if (NCores == 1) {
    terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(NCores)
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

  IASDT.R::CatTime(
    "Processing chunks on parallel, save each as RData files",
    Level = 1)

  GBIF_Data <- future.apply::future_lapply(
    X = ChunkList,
    FUN = function(x) {
      IASDT.R::GBIF_ReadChunk(
        ChunkFile = x, EnvFile = EnvFile, SaveRData = TRUE, ReturnData = FALSE,
        Overwrite = Overwrite)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = c(
      "IASDT.R", "purrr", "tibble", "terra", "tidyr", "dplyr", "sf"),
    future.globals = c("EnvFile", "Overwrite"))

  IASDT.R::CatDiff(
    InitTime = .StartTimeChunks, Prefix = "Finished in ", Level = 2)

  IASDT.R::CatTime("Reading processed chunks into a single dataset", Level = 1)
  if (all(file.exists(ChunkListRData))) {
    GBIF_Data <- future.apply::future_lapply(
      X = ChunkListRData, FUN = IASDT.R::LoadAs,
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

  IASDT.R::CatDiff(
    InitTime = .StartTimeChunks,
    Prefix = "Processing data chunks was finished in ", Level = 2)

  IASDT.R::CatTime(
    paste0(
      "A total of ", format(nrow(GBIF_Data), big.mark = ","), " observations"),
    Level = 2)

  IASDT.R::CatTime("Stopping cluster", Level = 1)
  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::CatTime("Saving `GBIF_Data` to disk", Level = 1)
  save(GBIF_Data, file = IASDT.R::Path(Path_GBIF, "GBIF_Data.RData"))

  rm(TaxaList, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # iNaturalist unique grids ----

  IASDT.R::CatTime("iNaturalist unique grids")

  IASDT.R::CatTime("Prepare input data for summarizing", Level = 1)
  iNaturalist_Others <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(species, CellCode, institutionCode) %>%
    dplyr::mutate(
      institutionCode = tidyr::replace_na(institutionCode, ""),
      institutionCode = institutionCode == "iNaturalist",
      institutionCode = dplyr::if_else(
        institutionCode, "iNaturalist", "Others"))

  IASDT.R::CatTime(
    "Number of unique iNaturalist grid cells per species", Level = 1)
  iNaturalist_Unique <- iNaturalist_Others %>%
    dplyr::group_by(species, CellCode) %>%
    tidyr::pivot_wider(
      names_from = institutionCode, values_from = institutionCode,
      values_fn = unique) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(Others), iNaturalist == "iNaturalist") %>%
    dplyr::count(species) %>%
    dplyr::rename(iNaturalist_Unique = n)

  IASDT.R::CatTime(
    "Number of grid cells for iNaturalist and other data sources",
    Level = 1)

  iNaturalist_Count <- iNaturalist_Others %>%
    dplyr::count(species, institutionCode) %>%
    tidyr::pivot_wider(names_from = institutionCode, values_from = n) %>%
    dplyr::select(species, iNaturalist, Others) %>%
    dplyr::left_join(iNaturalist_Unique, by = "species") %>%
    dplyr::select(species, Others, tidyselect::everything())

  IASDT.R::CatTime("Save iNaturalist summary data", Level = 1)
  save(
    iNaturalist_Count,
    file = IASDT.R::Path(Path_GBIF, "iNaturalist_Count.RData"))

  rm(
    iNaturalist_Others, iNaturalist_Unique, iNaturalist_Count,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Covert GBIF data to presence-only grids ----

  IASDT.R::CatTime("Covert GBIF data to presence-only grids")
  GBIF_Grid <- sf::st_drop_geometry(GBIF_Data) %>%
    # only unique combinations between species and grid ID
    dplyr::distinct(
      IAS_ID, taxon_name, Species_name, Species_name2,
      Species_File, CellCode) %>%
    # add polygons for grids
    dplyr::left_join(GridSf, by = "CellCode") %>%
    # ensure that the data is sf object
    sf::st_as_sf()

  IASDT.R::CatTime("Save presence-only grids", Level = 1)
  save(GBIF_Grid, file = IASDT.R::Path(Path_GBIF, "GBIF_Grid.RData"))
  invisible(gc())

  # # ..................................................................... ###

  # Number of observations or grid cells per species ----

  IASDT.R::CatTime("# observations or grid cells per species")

  IASDT.R::CatTime("# observations per species", Level = 1)
  SpNObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumOcc = n)

  IASDT.R::CatTime("# grids per species", Level = 1)
  SpNGrids <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(IAS_ID, taxon_name, Species_name) %>%
    dplyr::rename(GBIF_NumGrids = n)

  IASDT.R::CatTime(
    "Merge number of observations and grids per species",
    Level = 1)
  GBIF_NObsNGrid <- dplyr::full_join(
    SpNObs, SpNGrids,
    by = c("IAS_ID", "taxon_name", "Species_name"))

  IASDT.R::CatTime("Save as RData", Level = 1)
  save(GBIF_NObsNGrid, file = IASDT.R::Path(Path_GBIF, "GBIF_NObsNGrid.RData"))

  IASDT.R::CatTime("Save as xlsx", Level = 1)
  writexl::write_xlsx(
    x = GBIF_NObsNGrid, path = IASDT.R::Path(Path_GBIF, "GBIF_NObsNGrid.xlsx"))

  rm(GBIF_NObsNGrid, SpNObs, SpNGrids, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Summary maps ----
  IASDT.R::CatTime("Summary maps")

  ## Number of observations per grid  ----

  IASDT.R::CatTime("Number of observations per grid cell", Level = 1)
  GBIF_NObs <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumObs") %>%
    IASDT.R::setRastCRS() %>%
    terra::wrap()

  IASDT.R::CatTime(
    "Number of observations per grid cell (log scale)", Level = 1)
  GBIF_NObs_log <- terra::unwrap(GBIF_NObs) %>%
    log10() %>%
    IASDT.R::setRastCRS() %>%
    terra::wrap()

  IASDT.R::CatTime("Save as RData", Level = 2)
  save(
    GBIF_NObs, GBIF_NObs_log,
    file = IASDT.R::Path(Path_GBIF, "GBIF_NObs.RData"))

  invisible(gc())

  # # ................................... ###

  ## Number of species per grid cell -----

  IASDT.R::CatTime("Number of IAS per grid cell", Level = 1)
  GBIF_NSp <- sf::st_drop_geometry(GBIF_Grid) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(GridSf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = GridR, field = "n") %>%
    terra::mask(GridR) %>%
    stats::setNames("NumSpecies") %>%
    IASDT.R::setRastCRS() %>%
    terra::wrap()

  IASDT.R::CatTime("Number of IAS per grid cell (log)", Level = 1)
  GBIF_NSp_Log <- terra::unwrap(GBIF_NSp) %>%
    log10() %>%
    IASDT.R::setRastCRS() %>%
    terra::wrap()

  IASDT.R::CatTime("Save as RData", Level = 2)
  save(
    GBIF_NSp, GBIF_NSp_Log,
    file = IASDT.R::Path(Path_GBIF, "GBIF_NSp.RData"))

  rm(GridSf, GridR, envir = environment())
  invisible(gc())

  # # ................................... ###

  ## Plot_GBIF_Summary -----
  IASDT.R::CatTime("Plotting summary maps", Level = 1)

  GBIF_date <- GBIF_Metadata$StatusDetailed$modified %>%
    lubridate::date() %>%
    format("%d %B %Y")

  GBIF_DOI <- GBIF_Metadata$StatusDetailed$doi
  PlotTag <- format(Sys.Date(), "%d %B %Y")
  PlotTag <- stringr::str_glue(
    "DOI: {GBIF_DOI} ({GBIF_date})\nLast update: {PlotTag}")

  rm(GBIF_date, GBIF_DOI, GBIF_Grid, GBIF_Metadata, envir = environment())

  Plot_GBIF_Summary <- function(
    RstrMap, Title, LegendLabel = NULL, EU_Map = EuroBound) {

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
        EU_Map,
        mapping = ggplot2::aes(), color = "grey30", linewidth = 0.1,
        fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(data = terra::trim(RstrMap), maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = IASDT.R::integer_breaks()) +
      ggplot2::geom_sf(
        EU_Map,
        mapping = ggplot2::aes(), color = "grey40", linewidth = 0.075,
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
      label = PlotTag, fontface = "italic", color = "grey", x = 0.98,
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
    filename = IASDT.R::Path(Path_GBIF, "GBIF_Summary.jpeg"),
    width = 25, height = 25.8, units = "cm", quality = 100, res = 600)
  print(Plot)
  grDevices::dev.off()

  rm(EuroBound, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Species-specific data ----
  IASDT.R::CatTime("Species-specific data")

  SpList <- sf::st_drop_geometry(GBIF_Data) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name) %>%
    sort()

  # Species data --- sf ----
  IASDT.R::CatTime("Split species data - sf", Level = 1)

  # On LUMI, extracting species data in parallel using `furrr::future_walk` took
  # much longer time than working sequentially `purrr::walk` due to the
  # existence of the very large `GBIF_Data` object. Although it is possible to
  # integrate following chunk in the `GBIF_SpData` function, it would be much
  # faster to implement the following sequentially, then process species maps in
  # parallel later on, after the deletion of the `GBIF_Data` object
  purrr::walk(
    .x = SpList,
    .f = ~ {
      # Filter data on this species
      SpName <- IASDT.R::ReplaceSpace(.x) %>%
        # replace non-ascii multiplication symbol with x
        stringr::str_replace_all("\u00D7", "x") %>%
        stringr::str_replace_all("-", "")

      SpData <- dplyr::filter(GBIF_Data, Species_name == .x)
      OutFileSF <- IASDT.R::Path(Path_SpData, paste0(SpName, ".RData"))

      # Save if there is data
      if (nrow(SpData) > 0) {
        IASDT.R::SaveAs(InObj = SpData, OutObj = SpName, OutPath = OutFileSF)
      }
      return(invisible(NULL))
    },
    .progress = FALSE)

  rm(
    GBIF_Data, GBIF_NObs, GBIF_NObs_log, GBIF_NSp, GBIF_NSp_Log,
    envir = environment())
  invisible(gc())


  # Grid / raster / plotting ----
  IASDT.R::CatTime("Split species data - grid/raster/plot", Level = 1)

  IASDT.R::CatTime(
    paste0("Prepare working on parallel using ", NCores, " cores."),
    Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  IASDT.R::CatTime("Splitting species data on parallel", Level = 2)
  furrr::future_walk(
    .x = SpList, .f = IASDT.R::GBIF_SpData, EnvFile = EnvFile,
    Verbose = FALSE, PlotTag = PlotTag,
    .options = furrr::furrr_options(seed = TRUE, packages = "dplyr")
  )

  IASDT.R::CatTime("Stopping cluster", Level = 2)
  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Clean up chunk files ----
  if (DeleteChunks) {
    IASDT.R::CatTime("Clean up - remove temporary chunk files")
    list.files(Path_GBIF_Interim, full.names = TRUE) %>%
      fs::file_delete()
    fs::dir_delete(Path_GBIF_Interim)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nProcessing GBIF data was finished in ")

  return(invisible(NULL))
}
