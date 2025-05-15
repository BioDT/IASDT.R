#' Process and map Invasive Alien Species (IAS) data for the `IASDT`
#'
#' Processes and visualises Invasive Alien Species (IAS) distribution data from
#' GBIF, EASIN, and eLTER for the Invasive Alien Species Digital Twin (`IASDT`).
#' Merges pre-processed data, creates presence-absence rasters, summarises
#' distributions, and generates maps using helper functions.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 6.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "future::sequential", "future::multisession",
#'   "future::multicore", and "future::cluster". Defaults to
#'   `"future::multicore"` (`"future::multisession"` on Windows). See
#'   [future::plan()] and [ecokit::set_parallel()] for details.
#' @param species Character. Species name for distribution mapping.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default:
#'   `FALSE`.
#'
#' @section Functions details:
#' - **`IAS_process()`**: Merges pre-processed GBIF ([GBIF_process]), EASIN
#'   ([EASIN_process]), and eLTER ([eLTER_process]) data (run these first).
#'   Outputs `SpatRaster` distribution rasters, summary tables, and JPEG maps
#'   using `IAS_distribution()` and `IAS_plot()`.
#' - **`IAS_distribution()`**: Generates presence-absence maps (`.RData`,
#'   `.tif`)  for a species, including all grid cells in the study area and a
#'   set excluding cultivated/casual-only countries. Returns a tibble with
#'   presence counts (total, by source) and summary statistics for
#'   biogeographical regions
#' - **`IAS_plot()`**: Creates JPEG distribution maps from GBIF, EASIN, and
#'   eLTER data using `ggplot2`.

# # |------------------------------------------------------------------------| #
# IAS_process ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name IAS_data
#' @rdname IAS_data
#' @order 1

IAS_process <- function(
    env_file = ".env", n_cores = 6L, strategy = "future::multicore") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("env_file", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_cores")

  rm(AllArgs, envir = environment())

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
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

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- Path_TaxaInfo_RData <- taxon_name <- Path_TaxaStand <-
    Path_HabAff <- Species <- PA_Map <- Species_name <- NCells_All <-
    synhab_name <- IAS_ID <- Count <- Hab <- EU_Bound <- Threshold <- NSp <-
    taxon_name <- PA_Masked_Map <- NCells_Naturalized <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_GBIF", "DP_R_GBIF_processed", TRUE, FALSE,
    "Path_EASIN", "DP_R_EASIN_processed", TRUE, FALSE,
    "Path_eLTER", "DP_R_eLTER_processed", FALSE, TRUE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "Path_PA", "DP_R_PA", FALSE, FALSE,
    "Path_HabAff", "DP_R_HabAff", FALSE, TRUE,
    "Path_TaxaStand", "DP_R_Taxa_stand", FALSE, TRUE,
    "Path_TaxaInfo_RData", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "lubridate", "IASDT.R", "purrr", "stringr", "readr", "fs",
      "sf", "terra", "readxl", "tidyr", "tidyselect", "ggplot2", "ggtext",
      "grid", "tidyterra", "cowplot", "scales", "tibble", "magrittr", "ragg",
      "gtools", "ecokit"),
    strategy = strategy)

  # # ..................................................................... ###

  # Reading input data and create directories ------

  ecokit::cat_time("Reading input data and create directories")

  Path_PA_JPEG <- fs::path(Path_PA, "JPEG_Maps")
  fs::dir_create(Path_PA_JPEG)

  # last update info
  LastUpdate <- stringr::str_glue(
    'Last update: {format(Sys.Date(), "%d %B %Y")}')

  ## Standardized taxonomy -----
  ecokit::cat_time("Standardized taxonomy", level = 1L)
  # list of original taxonomy (including dummy ID column `taxon_id_`)
  TaxaList_Original <- readRDS(Path_TaxaStand) %>%
    dplyr::select(tidyselect::all_of(c("taxon_id_", "taxon_name")))

  # # .................................... ###

  ## TaxaInfo ----
  ecokit::cat_time("TaxaInfo", level = 1L)
  TaxaList <- ecokit::load_as(Path_TaxaInfo_RData)
  TaxaList_Distinct <- dplyr::distinct(TaxaList, taxon_name, IAS_ID)

  # # .................................... ###

  ## Species habitat affinity -----
  ecokit::cat_time("Species habitat affinity", level = 1L)
  SynHab_List <- tibble::tribble(
    ~SynHab_code, ~SynHab_name,
    "1", "Forests",
    "2", "Open_forests",
    "3", "Scrub",
    "4", "Grasslands",
    "4a", "Natural_grasslands",
    "4b", "Human_maintained_grasslands",
    "5", "Sandy",
    "6", "Rocky",
    "7", "Dryland",
    "8", "Saline",
    "9", "Riparian",
    "10", "Wetland",
    "11", "Aquatic",
    "12", "Man_made",
    "12a", "Ruderal_habitats",
    "12b", "Agricultural_habitats")

  Sp_SynHab <- readRDS(Path_HabAff) %>%
    dplyr::rename(SynHab_name = synhab_name) %>%
    dplyr::left_join(SynHab_List, by = "SynHab_name") %>%
    dplyr::mutate(Count = TRUE) %>%
    tidyr::pivot_wider(
      names_from = c("SynHab_code", "SynHab_name"),
      names_prefix = "Hab_", values_from = Count) %>%
    # add taxon_name column
    dplyr::left_join(TaxaList_Original, by = "taxon_id_") %>%
    # add IAS_ID column
    dplyr::left_join(TaxaList_Distinct, by = "taxon_name") %>%
    dplyr::select(-tidyselect::all_of(c("taxon_id_", "taxon_name"))) %>%
    dplyr::select(IAS_ID, gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::distinct() %>%
    dplyr::arrange(IAS_ID)

  rm(TaxaList_Original, TaxaList_Distinct, envir = environment())

  # # ..................................................................... ###

  # Species-specific data ------

  ecokit::cat_time("Species-specific data")
  .StartTimeDist <- lubridate::now(tzone = "CET")

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  # # .................................... ###

  ## Species-specific data in parallel ----
  ecokit::cat_time("Species-specific data in parallel", level = 1L)

  Sp_PA_Data <- future.apply::future_lapply(
    X = sort(unique(TaxaList$Species_name)),
    FUN = function(x) {

      # file name
      sp_file <- dplyr::filter(TaxaList, Species_name == x)$Species_File
      file_Summary <- fs::path(
        Path_PA, "SpSummary", paste0(sp_file, "_Summary.RData"))
      file_PA <- fs::path(Path_PA, "RData", paste0(sp_file, "_PA.RData"))
      file_tif_All <- fs::path(Path_PA, "tif", paste0(sp_file, "_All.tif"))
      file_tif_Masked <- fs::path(
        Path_PA, "tif", paste0(sp_file, "_Masked.tif"))
      file_jpeg <- fs::path(Path_PA, "JPEG_Maps", paste0(sp_file, ".jpeg"))
      files_all <- c(
        file_Summary, file_PA, file_tif_All, file_tif_Masked, file_jpeg)

      # Maximum attempts
      max_attempts <- 5
      attempt <- 1

      repeat {

        # Run the distribution function
        Species_Data <- IASDT.R::IAS_distribution(
          species = x, env_file = env_file, verbose = FALSE)

        # allow some time for the files to be created
        Sys.sleep(2)

        # Check for the existence and validity of all files
        all_okay <- all(
          ecokit::check_data(file_PA, warning = FALSE),
          ecokit::check_data(file_Summary, warning = FALSE),
          ecokit::check_tiff(file_tif_All, warning = FALSE),
          ecokit::check_tiff(file_tif_Masked, warning = FALSE),
          fs::file_exists(file_jpeg))

        if (all_okay) {
          break
        }

        # delete files for unsuccessful try
        purrr::walk(
          .x = files_all,
          .f = ~ if (fs::file_exists(.x)) {
            fs::file_delete(.x)
          })

        # Increment the attempt counter
        attempt <- attempt + 1

        # Exit if max attempts reached
        if (attempt > max_attempts) {
          warning("Maximum attempts reached for species: ", x, call. = FALSE)
          break
        }
      }

      return(Species_Data)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "Path_PA", "TaxaList")) %>%
    dplyr::bind_rows()

  # # .................................... ###

  ## Stopping cluster ----
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("future::sequential", gc = TRUE)
  }

  ecokit::cat_diff(
    init_time = .StartTimeDist,
    prefix = "Processing Species-specific data took ",
    msg_n_lines = 1, level = 2L)

  # # ..................................................................... ###

  # Merge Species-specific summary info -----
  ecokit::cat_time("Merge Species-specific summary info")

  Sp_PA_Data <- Sp_PA_Data %>%
    # Split number of grid cell per country / biogeographical region as separate
    # column
    tidyr::unnest_wider(c("BioRegs_DT", "BioRegsMask_DT", "SpCountry")) %>%
    dplyr::mutate(
      # replace NA for biogeographical regions and country analysis with 0
      dplyr::across(
        tidyselect::matches("^CNT_|^BioReg_"),
        ~ dplyr::if_else(is.na(.x), 0L, .x))) %>%
    # change column order
    dplyr::select(
      Species:NCells_Naturalized, PA_Map, PA_Masked_Map,
      tidyselect::matches("^GBIF"),
      tidyselect::matches("^EASIN"),
      tidyselect::matches("^eLTER"),
      tidyselect::matches("^iNat"),
      tidyselect::matches("^BioReg_"),
      tidyselect::matches("^BioRegSumm_"),
      tidyselect::matches("^BioReg_Masked_"),
      tidyselect::matches("^BioRegMaskSumm_"),
      tidyselect::matches("^CNT_"),
      tidyselect::everything()) %>%
    dplyr::rename(Species_name = Species) %>%
    dplyr::full_join(
      y = dplyr::distinct(
        dplyr::select(TaxaList, -tidyselect::all_of("speciesKey"))),
      by = "Species_name") %>%
    dplyr::select(
      tidyselect::all_of(
        c(
          "IAS_ID", "taxon_name", "Species_name",
          "Species_name2", "Species_File")),
      tidyselect::everything()
    ) %>%
    dplyr::left_join(Sp_SynHab, by = "IAS_ID")

  rm(TaxaList, envir = environment())

  ## Save summary results -----
  ecokit::cat_time("Save summary results", level = 1L)

  ecokit::cat_time("`RData`", level = 2L)
  save(Sp_PA_Data, file = fs::path(Path_PA, "Sp_PA_Data.RData"))

  ecokit::cat_time("Summary data without maps", level = 1L)
  ecokit::cat_time("csv format", level = 2L)
  Sp_PA_Summary_DF <- Sp_PA_Data %>%
    dplyr::select(
      -tidyselect::all_of(
        c("GBIF_R", "EASIN_R", "eLTER_R", "PA_Map", "PA_Masked_Map")))

  readr::write_excel_csv(
    Sp_PA_Summary_DF, file = fs::path(Path_PA, "Sp_PA_Summary_DF.csv"),
    progress = FALSE)

  ecokit::cat_time("RData format", level = 2L)
  save(
    Sp_PA_Summary_DF,
    file = fs::path(Path_PA, "Sp_PA_Summary_DF.RData"))

  # # ..................................................................... ###

  # Summary of merged data ------
  ecokit::cat_time("Summary of merged data")

  ## Number of IAS per grid cell -----
  ecokit::cat_time("# IAS per grid cell", level = 1L)

  IAS_NumSp <- dplyr::filter(Sp_PA_Data, NCells_All > 0) %>%
    dplyr::pull(PA_Map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("IAS_NumSp") %>%
    # Ensure that values are read from memory
    terra::toMemory()

  # save as RData
  ecokit::save_as(
    object = terra::wrap(IAS_NumSp), object_name = "IAS_NumSp",
    out_path = fs::path(Path_PA, "IAS_NumSp.RData"))

  # save as tif
  raster::writeRaster(
    x = IAS_NumSp, overwrite = TRUE,
    filename = fs::path(Path_PA, "IAS_NumSp.tif"))

  # # .................................... ###

  ## Number of IAS per grid cell - masked data -----

  IAS_NumSp_Masked <- dplyr::filter(Sp_PA_Data, NCells_Naturalized > 0) %>%
    dplyr::pull(PA_Masked_Map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("IAS_NumSp_Masked") %>%
    # Ensure that values are read from memory
    terra::toMemory()

  # save as RData
  ecokit::save_as(
    object = terra::wrap(IAS_NumSp_Masked), object_name = "IAS_NumSp_Masked",
    out_path = fs::path(Path_PA, "IAS_NumSp_Masked.RData"))

  # save as tif
  raster::writeRaster(
    x = IAS_NumSp_Masked, overwrite = TRUE,
    filename = fs::path(Path_PA, "IAS_NumSp_Masked.tif"))

  # # .................................... ###

  ## Plotting summary of IAS data -----
  ecokit::cat_time("Plotting summary of IAS data", level = 1L)

  EUBound <- ecokit::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")
  MapLimX <- c(2600000, 6550000)
  MapLimY <- c(1450000, 5410000)

  # # +++++++++++++++++++++++++++++++++ ###

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0.1, "cm"),
      plot.title = ggplot2::element_text(
        size = 10, color = "blue", face = "bold", hjust = 0,
        margin = ggplot2::margin(2, 0, 2, 0)),
      plot.subtitle = ggplot2::element_text(
        size = 8, color = "darkgrey", face = "italic", hjust = 0,
        margin = ggplot2::margin(1, 0, 1, 0)),
      strip.text = ggplot2::element_text(size = 6, face = "bold"),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.5, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 6),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.title = ggplot2::element_text(
        color = "blue", size = 6, face = "bold", hjust = 0),
      legend.position = "inside",
      legend.position.inside = c(0.92, 0.825),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.text.y = ggplot2::element_text(size = 7, hjust = 0.5, angle = 90),
      axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.04, "cm"),
      panel.spacing = grid::unit(0.3, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.1, colour = "grey40", linetype = 2),
      panel.border = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE, panel.background = ggplot2::element_rect(fill = NA))

  # # +++++++++++++++++++++++++++++++++ ###

  Plot_Nsp <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = terra::classify(IAS_NumSp, cbind(0, NA)), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid",
      subtitle = "All data (including cultivated or casual observations)",
      fill = "# IAS") +
    PlottingTheme

  Plot_Nsp_log <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = log10(terra::classify(IAS_NumSp, cbind(0, NA))), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid (log10 scale)",
      subtitle = "All data (including cultivated or casual observations)",
      fill = "# IAS\n(log10)") +
    PlottingTheme

  Plot <- cowplot::plot_grid(Plot_Nsp, Plot_Nsp_log, ncol = 2, nrow = 1) +
    cowplot::draw_label(
      label = LastUpdate, color = "grey", x = 0.99, y = 0.98,
      size = 7, hjust = 1)

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_PA, "IAS_NumSpecies.jpeg"),
    width = 30, height = 15.5, res = 600, quality = 100, units = "cm")
  print(Plot)
  grDevices::dev.off()

  # # +++++++++++++++++++++++++++++++++ ###

  Plot_Nsp_Masked <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = terra::classify(IAS_NumSp_Masked, cbind(0, NA)), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid",
      subtitle = "Excluding cultivated or casual observations",
      fill = "# IAS") +
    PlottingTheme

  Plot_Nsp_Masked_log <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = log10(terra::classify(IAS_NumSp_Masked, cbind(0, NA))),
      maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid (log10 scale)",
      subtitle = "Excluding cultivated or casual observations",
      fill = "# IAS\n(log10)") +
    PlottingTheme

  Plot <- cowplot::plot_grid(
    Plot_Nsp_Masked, Plot_Nsp_Masked_log, ncol = 2, nrow = 1) +
    cowplot::draw_label(
      label = LastUpdate, color = "grey", x = 0.99, y = 0.98,
      size = 7, hjust = 1)

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_PA, "IAS_NumSpecies_Masked.jpeg"),
    width = 30, height = 15.5, res = 600, quality = 100, units = "cm")
  print(Plot)
  grDevices::dev.off()

  # # +++++++++++++++++++++++++++++++++ ###

  rm(
    IAS_NumSp, Plot_Nsp, Plot_Nsp_log, PlottingTheme, EUBound, MapLimX,
    MapLimY, Plot_Nsp_Masked_log, Plot_Nsp_Masked,
    envir = environment())

  # # .................................... ###

  ## # Number of IAS to be used in the model - per habitat type ------

  Threshold_Theme <- ggplot2::theme(
    plot.margin = ggplot2::margin(0.05, 0.15, 0.05, 0, "cm"),
    plot.title = ggplot2::element_text(
      size = 11, color = "black", hjust = 0.05, face = "italic",
      margin = ggplot2::margin(0.05, 0, 0.05, 0, "cm")),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(
      size = 6, face = "bold", margin = ggplot2::margin(0, 0, 0, 0, "cm")),
    legend.box.spacing = grid::unit(-0.05, "cm"),
    legend.key.spacing.x = grid::unit(0.75, "cm"),
    legend.key = ggplot2::element_rect(fill = scales::alpha("white", 0)),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.1, "cm"),
    axis.title.y.left = ggplot2::element_text(size = 9, vjust = -1),
    axis.title.x = ggplot2::element_text(size = 9, vjust = 1),
    axis.text.y.left = ggplot2::element_text(size = 7),
    axis.text.x.bottom = ggplot2::element_text(size = 7),
    axis.text.y.right = ggplot2::element_blank(),
    axis.text.x.top = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      linewidth = 0.1, colour = "grey40", linetype = 2),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank())

  # # +++++++++++++++++++++++++++++++++ ###

  NSp_DT <- dplyr::select(Sp_PA_Data, Species_name, NCells_All)
  NSp_DT <- purrr::map_dfr(
    .x = sort(unique(NSp_DT$NCells_All)),
    .f = ~ tibble::tibble(
      Threshold = .x, NSp = sum(NSp_DT$NCells_All >= .x, na.rm = TRUE))) %>%
    dplyr::filter(Threshold < 450)

  NSp_DT_Masked <- dplyr::select(Sp_PA_Data, Species_name, NCells_Naturalized)
  NSp_DT_Masked <- purrr::map_dfr(
    .x = sort(unique(NSp_DT_Masked$NCells_Naturalized)),
    .f = ~ tibble::tibble(
      Threshold = .x,
      NSp = sum(NSp_DT_Masked$NCells_Naturalized >= .x, na.rm = TRUE))) %>%
    dplyr::filter(Threshold < 450)

  # # +++++++++++++++++++++++++++++++++ ###

  NSp_Hab_DT <- purrr::map_dfr(
    .x = c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats"),
    .f = function(Hab) {

      # nolint start
      CurrDT <- dplyr::filter(Sp_PA_Data, !!as.symbol(Hab))
      # nolint end

      purrr::map_dfr(
        sort(unique(Sp_PA_Data$NCells_All)),
        ~ tibble::tibble(
          Hab = Hab, Threshold = .x,
          NSp = sum(CurrDT$NCells_All >= .x, na.rm = TRUE)))
    }
  ) %>%
    dplyr::mutate(
      Hab = stringr::str_remove(Hab, "^Hab_"),
      Hab = stringr::str_replace(Hab, "_", " - "),
      Hab = stringr::str_replace_all(Hab, "_", " "),
      Hab = forcats::fct_inorder(Hab)) %>%
    dplyr::slice(gtools::mixedorder(Hab)) %>%
    dplyr::filter(Threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = Threshold, y = NSp, group = Hab, colour = Hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = NSp_DT, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = Threshold, y = NSp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "Threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(
      title = "All data (including cultivated or casual observations)") +
    Threshold_Theme +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1, label.position = "bottom",
        override.aes = list(linewidth = 1.5)))

  NSp_Hab_Masked_DT <- purrr::map_dfr(
    .x = c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats"),
    .f = function(Hab) {

      # nolint start
      CurrDT <- dplyr::filter(Sp_PA_Data, !!as.symbol(Hab))
      # nolint end

      purrr::map_dfr(
        sort(unique(Sp_PA_Data$NCells_Naturalized)),
        ~ tibble::tibble(
          Hab = Hab, Threshold = .x,
          NSp = sum(CurrDT$NCells_Naturalized >= .x, na.rm = TRUE)))
    }
  ) %>%
    dplyr::mutate(
      Hab = stringr::str_remove(Hab, "^Hab_"),
      Hab = stringr::str_replace(Hab, "_", " - "),
      Hab = stringr::str_replace_all(Hab, "_", " "),
      Hab = forcats::fct_inorder(Hab)) %>%
    dplyr::slice(gtools::mixedorder(Hab)) %>%
    dplyr::filter(Threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = Threshold, y = NSp, group = Hab, colour = Hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = NSp_DT_Masked, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = Threshold, y = NSp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "Threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(title = "Excluding cultivated or casual observations") +
    Threshold_Theme +
    ggplot2::theme(legend.position = "none")

  Legend <- ggpubr::as_ggplot(ggpubr::get_legend(NSp_Hab_DT))
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "Number of IAS to be used in the models based on the arbitrary ",
        "selection of # of presence grid cells per species"),
      fontface = "bold", colour = "blue")

  Plot <- cowplot::plot_grid(
    title,
    cowplot::plot_grid(
      (NSp_Hab_DT + ggplot2::theme(legend.position = "none")),
      NSp_Hab_Masked_DT,
      ncol = 2, nrow = 1),
    Legend,
    ncol = 1, rel_heights = c(0.05, 1, 0.05))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_PA, "IAS_NSp_threshold_Hab.jpeg"),
    width = 30, height = 17, res = 600, quality = 100, units = "cm")
  print(Plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Function Summary ----
  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing species data was finished in ", ... = "\n")

  return(invisible(NULL))
}
