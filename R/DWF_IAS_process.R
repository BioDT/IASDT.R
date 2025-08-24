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
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
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
#'   set excluding cultivated/casual-only countries. Returns a file path to a
#'   tibble with presence counts (total, by source) and summary statistics for
#'   biogeographical regions
#' - **`IAS_plot()`**: Creates JPEG distribution maps from GBIF, EASIN, and
#'   eLTER data using `ggplot2`.

# # |------------------------------------------------------------------------| #
# IAS_process ----
# # |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name IAS_data
#' @rdname IAS_data
#' @order 1

IAS_process <- function(
    env_file = ".env", n_cores = 6L, strategy = "multisession") {

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- Path_TaxaInfo_RData <- taxon_name <- Path_TaxaStand <- img_valid <-
    Path_HabAff <- Species <- PA_Map <- Species_name <- NCells_All <-
    synhab_name <- IAS_ID <- Count <- Hab <- EU_Bound <- Threshold <- NSp <-
    taxon_name <- PA_Masked_Map <- NCells_Naturalized <- Species_File <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

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
  invisible(gc())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "lubridate", "IASDT.R", "purrr", "stringr", "readr", "fs",
      "sf", "terra", "readxl", "tidyr", "tidyselect", "ggplot2", "ggtext",
      "grid", "cowplot", "scales", "tibble", "magrittr", "ragg", "qs2",
      "gtools", "ecokit"),
    strategy = strategy)

  # # ..................................................................... ###

  # Reading input data and create directories ------

  ecokit::cat_time("Reading input data and create directories")

  Path_PA_Summary <- fs::path(Path_PA, "PA_summary")
  Path_PA_tif <- fs::path(Path_PA, "PA_tif")
  Path_PA_RData <- fs::path(Path_PA, "PA_RData")
  Path_PA_JPEG <- fs::path(Path_PA, "Distribution_JPEG")

  # list of paths; create if not exist
  Paths_All <- list(
    summary = Path_PA_Summary, PA_tif = Path_PA_tif,
    PA_RData = Path_PA_RData, jpeg = Path_PA_JPEG)
  purrr::walk(Paths_All, fs::dir_create)

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
  invisible(gc())

  Sp_taxa <- TaxaList %>%
    dplyr::arrange(IAS_ID) %>%
    dplyr::distinct(Species_name, Species_File)

  # # ..................................................................... ###

  # Species-specific data ------

  ecokit::cat_time("Species-specific data")
  .StartTimeDist <- lubridate::now(tzone = "CET")

  # # .................................... ###

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # .................................... ###

  ## Species-specific data in parallel ----
  ecokit::cat_time("Species-specific data in parallel", level = 1L)

  ecokit::cat_time(
    paste0("There are ", nrow(Sp_taxa), " species to process"), level = 2L)

  Sp_PA_Data <- future.apply::future_lapply(
    X = seq_len(nrow(Sp_taxa)),
    FUN = function(x) {

      # species file name
      sp_file <- Sp_taxa$Species_File[x]
      sp_Name <- Sp_taxa$Species_name[x]

      if (length(sp_file) != 1 || is.na(sp_file) || !nzchar(sp_file) ||
          length(sp_Name) != 1 || is.na(sp_Name) || !nzchar(sp_Name)) {
        ecokit::stop_ctx(
          "Species file name not unique or empty", sp_Name = sp_Name,
          sp_file = sp_file)
      }

      # output paths
      file_summary <- fs::path(Paths_All$summary, paste0(sp_file, ".qs2"))
      file_PA <- fs::path(Paths_All$PA_RData, paste0(sp_file, "_PA.RData"))
      file_tif_All <- fs::path(Paths_All$PA_tif, paste0(sp_file, "_All.tif"))
      file_tif_Masked <- fs::path(
        Paths_All$PA_tif, paste0(sp_file, "_Masked.tif"))

      all_okay <- all(
        ecokit::check_data(file_summary, warning = FALSE),
        ecokit::check_data(file_PA, warning = FALSE),
        ecokit::check_tiff(file_tif_All, warning = FALSE),
        ecokit::check_tiff(file_tif_Masked, warning = FALSE))

      if (all_okay) {
        return(file_summary)
      }

      # Maximum attempts
      max_attempts <- 8
      attempt <- 1

      repeat {

        # Run the distribution function
        Species_Data <- tryCatch(
          expr = {
            IASDT.R::IAS_distribution(
              species = sp_Name, env_file = env_file, verbose = FALSE)
          },
          error = function(e) {
            NULL
          },
          warning = function(w) {
            NULL
          })

        if (is.null(Species_Data) || !is.data.frame(Species_Data)) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            break
          }
          next
        }

        if (ncol(Species_Data) != 2 && nrow(Species_Data) != 1) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            break
          }
          next
        }

        Sys.sleep(0.5)

        # species without data
        if (is.na(Species_Data$PA_summary)) {
          break
        }

        # Check for the existence and validity of all files
        all_okay <- all(
          ecokit::check_data(file_summary, warning = FALSE),
          ecokit::check_data(file_PA, warning = FALSE),
          ecokit::check_tiff(file_tif_All, warning = FALSE),
          ecokit::check_tiff(file_tif_Masked, warning = FALSE))

        if (all_okay) {
          break
        }

        attempt <- attempt + 1
        if (attempt > max_attempts) {
          break
        }
        invisible(gc())
      }

      return(Species_Data$PA_summary)

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "Paths_All", "Sp_taxa"))

  invisible(gc())

  # # .................................... ###

  ## Read processed species-specific data -----
  ecokit::cat_time("Read processed species-specific data", level = 1L)

  Sp_PA_Data <- unlist(Sp_PA_Data) %>%
    # remove NA object from a list; for species with no observations
    purrr::discard(is.na) %>%
    future.apply::future_lapply(
      FUN = ecokit::load_as, future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export) %>%
    dplyr::bind_rows()

  # # .................................... ###

  ## Stopping cluster ----

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  # # .................................... ###

  ecokit::cat_time("Identify species with no data", level = 1L)

  # Species with no data
  sp_no_data <- setdiff(TaxaList$IAS_ID, as.integer(Sp_PA_Data$species_ID))
  sp_no_data <- dplyr::filter(TaxaList, IAS_ID %in% sp_no_data)

  if (nrow(sp_no_data) > 0) {
    readr::write_tsv(
      x = sp_no_data, file = fs::path(Path_PA, "species_no_data.csv"),
      col_names = TRUE)
    ecokit::cat_time(
      paste0(
        length(unique(Sp_PA_Data$species_ID)),  " species were processed; ",
        nrow(sp_no_data), " species have no data. ",
        "See `species_no_data.csv` for details"),
      level = 2L)
  } else {
    ecokit::cat_time(
      paste0("All ", length(TaxaList$Species_name), " species were processed"),
      level = 2L)
  }

  # # .................................... ###

  ecokit::cat_diff(
    init_time = .StartTimeDist,
    prefix = "Processing Species-specific data took ",
    msg_n_lines = 1, level = 2L)

  rm(sp_no_data, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Merge species summary info -----
  ecokit::cat_time("Merge species summary info")

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
      tidyselect::everything()) %>%
    dplyr::left_join(Sp_SynHab, by = "IAS_ID")

  rm(TaxaList, envir = environment())

  ## Save summary results -----
  ecokit::cat_time("Save summary results", level = 1L)

  ecokit::cat_time("`qs2`", level = 2L)
  save(Sp_PA_Data, file = fs::path(Path_PA, "Sp_PA_Data.qs2"))

  ecokit::cat_time("Summary data without maps", level = 1L)
  ecokit::cat_time("csv format", level = 2L)
  Sp_PA_Summary_DF <- Sp_PA_Data %>%
    dplyr::select(tidyselect::where(~ !is.list(.x)))

  readr::write_excel_csv(
    x = Sp_PA_Summary_DF, file = fs::path(Path_PA, "Sp_PA_Summary_DF.csv"),
    progress = FALSE)

  ecokit::cat_time("RData format", level = 2L)
  save(Sp_PA_Summary_DF, file = fs::path(Path_PA, "Sp_PA_Summary_DF.RData"))

  # only keep selected columns for summary plots
  Sp_PA_Data <- Sp_PA_Data %>%
    dplyr::select(
      Species_name, NCells_All, NCells_Naturalized,
      PA_Map, PA_Masked_Map, tidyselect::starts_with("Hab_"))

  rm(Sp_PA_Summary_DF, envir = environment())
  invisible(gc())

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

  Sp_PA_Data <- dplyr::select(Sp_PA_Data, -PA_Map)
  invisible(gc())

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

  Sp_PA_Data <- dplyr::select(Sp_PA_Data, -PA_Masked_Map)
  invisible(gc())

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

  rm(IAS_NumSp, Plot, envir = environment())
  invisible(gc())

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
    Plot_Nsp, Plot_Nsp_log, PlottingTheme, EUBound, MapLimX,
    MapLimY, Plot_Nsp_Masked_log, Plot_Nsp_Masked, IAS_NumSp_Masked, Plot,
    envir = environment())
  invisible(gc())

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
    }) %>%
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
    }) %>%
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

  rm(
    Sp_PA_Data, NSp_Hab_DT, NSp_Hab_Masked_DT, Legend, title, Plot,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting species distribution -----
  ecokit::cat_time("Plotting species distribution")
  .StartTimePlot <- lubridate::now(tzone = "CET")

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # .................................... ###

  ## Plotting species distribution in parallel ----
  ecokit::cat_time("Plotting species distribution in parallel", level = 1L)

  Sp_plots <- future.apply::future_lapply(
    X = seq_len(nrow(Sp_taxa)),
    FUN = function(x) {

      # Maximum attempts
      max_attempts <- 4
      attempt <- 1

      repeat {
        sp_plot <- tryCatch(
          IASDT.R::IAS_plot(
            species = Sp_taxa$Species_name[x], env_file = env_file))

        if (inherits(sp_plot, "try-error") || !is.character(sp_plot) ||
            !fs::file_exists(sp_plot)) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            break
          }
          next
        }

        img_valid <- ecokit::check_image(sp_plot, warning = FALSE)

        if (img_valid) {
          break
        }

        attempt <- attempt + 1
        if (attempt > max_attempts) {
          break
        }
      }

      tibble::tibble(
        species = Sp_taxa$Species_name[x],
        img_file = sp_plot, img_valid = img_valid)

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "Paths_All", "Sp_taxa"))


  # validate that all plots were created
  Sp_plots <- dplyr::bind_rows(Sp_plots)
  invalid_plots <- dplyr::filter(Sp_plots, !img_valid)
  valid_plots <- dplyr::filter(Sp_plots, img_valid)

  if (nrow(invalid_plots) > 0) {
    readr::write_tsv(
      x = invalid_plots, file = fs::path(Path_PA, "species_invalid_plots.csv"),
      col_names = TRUE)

    ecokit::cat_time(
      paste0(
        "Distribution maps for ", nrow(valid_plots),
        " species were processed; ", nrow(invalid_plots),
        " species have no data or failed to be plotted.",
        "See `species_invalid_plots.csv` for details"),
      level = 2L)
  } else {
    ecokit::cat_time(
      paste0(
        "Distribution maps for all ", nrow(valid_plots),
        " species were plotted"),
      level = 2L)
  }

  rm(Sp_plots, valid_plots, invalid_plots, envir = environment())

  # # .................................... ###

  ## Stopping cluster ----

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }
  # # .................................... ###

  ecokit::cat_diff(
    init_time = .StartTimePlot,
    prefix = "Plotting Species distribution took ",
    msg_n_lines = 1, level = 2L)

  # # ..................................................................... ###

  # Function Summary ----
  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing species data was finished in ", ... = "\n")

  return(invisible(NULL))
}
