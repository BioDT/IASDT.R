#' Process and map Naturalized Alien Plant Species (NAPS) data for the `IASDT`
#'
#' Processes and visualises Naturalized Alien Plant Species (NAPS) distribution
#' data from GBIF, EASIN, and eLTER for the Invasive Alien Species Digital Twin
#' (`IASDT`).
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
#' @param dist_citizen Numeric. Distance in km for spatial filtering of citizen
#'   science data in GBIF. Default is 100 km.
#'
#' @section Functions details:
#' - **`naps_process()`**: Merges pre-processed GBIF ([gbif_process]), EASIN
#'   ([easin_process]), and eLTER ([elter_process]) data (run these first).
#'   Outputs `SpatRaster` distribution rasters, summary tables, and JPEG maps
#'   using `naps_distribution()` and `naps_plot()`.
#' - **`naps_distribution()`**: Generates presence-absence maps (`.RData`,
#'   `.tif`)  for a species, including all grid cells in the study area and a
#'   set excluding cultivated/casual-only countries. Returns a file path to a
#'   tibble with presence counts (total, by source) and summary statistics for
#'   biogeographical regions
#' - **`naps_plot()`**: Creates JPEG distribution maps from GBIF, EASIN, and
#'   eLTER data using `ggplot2`.

# # |------------------------------------------------------------------------| #
# naps_process ----
# # |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name naps_data
#' @rdname naps_data
#' @order 1

naps_process <- function(
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
  path_pa <- path_taxa_info_r <- taxon_name <- img_valid <- pa_cz_excluded <-
    path_hab_affinity <- species <- pa_map <- species_name <- n_cells_all <-
    synhab_name <- ias_id <- count <- hab <- eu_boundaries <- threshold <-
    n_sp <- pa_masked_map <- n_cells_naturalized <- species_file <-
    synhab_code <- pa_cz_excluded_masked <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_gbif", "DP_R_gbif_processed", TRUE, FALSE,
    "path_easin", "DP_R_easin_processed", TRUE, FALSE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE,
    "path_pa", "DP_R_pa", FALSE, FALSE,
    "path_hab_affinity", "DP_R_hab_affinity", FALSE, TRUE,
    "path_taxa_info_r", "DP_R_taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)

  rm(env_vars_to_read, envir = environment())
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

  path_pa_summary <- fs::path(path_pa, "pa_summary")
  path_pa_tif <- fs::path(path_pa, "pa_tif")
  path_pa_r <- fs::path(path_pa, "pa_raster")
  path_pa_jpeg <- fs::path(path_pa, "distribution_jpeg")

  # list of paths; create if not exist
  paths_all <- list(
    summary = path_pa_summary, pa_tif = path_pa_tif,
    pa_rdata = path_pa_r, jpeg = path_pa_jpeg)
  purrr::walk(paths_all, fs::dir_create)

  # last update info
  last_update <- stringr::str_glue(
    'Last update: {format(Sys.Date(), "%d %B %Y")}')

  # # .................................... ###

  ## Taxa info ----
  ecokit::cat_time("Taxa info", level = 1L)
  taxa_list <- ecokit::load_as(path_taxa_info_r)
  taxa_list_distinct <- dplyr::distinct(
    taxa_list, taxon_name, ias_id, species_name)

  # # .................................... ###

  ## Species habitat affinity -----
  ecokit::cat_time("Species habitat affinity", level = 1L)
  synhab_list <- tibble::tribble(
    ~synhab_code, ~synhab_name,
    "1", "forests",
    "2", "open_forests",
    "3", "scrub",
    "4", "grasslands",
    "4a", "natural_grasslands",
    "4b", "human_maintained_grasslands",
    "5", "sandy",
    "6", "rocky",
    "7", "dryland",
    "8", "saline",
    "9", "riparian",
    "10", "wetland",
    "11", "aquatic",
    "12", "man_made",
    "12a", "ruderal_habitats",
    "12b", "agricultural_habitats") %>%
    dplyr::mutate(
      synhab_name_full = paste0("hab_", synhab_code, "_", synhab_name))

  species_synhab <- readRDS(path_hab_affinity) %>%
    dplyr::mutate(synhab_name = stringr::str_to_lower(synhab_name)) %>%
    dplyr::left_join(synhab_list, by = "synhab_name") %>%
    dplyr::mutate(count = TRUE, synhab_name_full = NULL) %>%
    ecokit::arrange_alphanum(synhab_code) %>%
    tidyr::pivot_wider(
      names_from = c("synhab_code", "synhab_name"),
      names_prefix = "hab_", values_from = count) %>%
    # add ias_id column
    dplyr::full_join(taxa_list_distinct, by = "taxon_name") %>%
    dplyr::select(
      ias_id, species_name, taxon_name,
      gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::starts_with("hab"), ~tidyr::replace_na(.x, FALSE))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(ias_id)

  species_synhab_new <- dplyr::filter(species_synhab, is.na(ias_id))
  if (nrow(species_synhab_new) > 0) {
    species_synhab_new_f <- fs::path(path_pa, "species_synhab_new.txt")
    ecokit::cat_time(
      paste0(
        nrow(species_synhab_new),
        " species do not exist in taxa standardization data"),
      level = 2L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0("see `", species_synhab_new_f, "` for more details"),
      level = 2L, cat_timestamp = FALSE)
    cat(
      sort(species_synhab_new$taxon_name),
      file = species_synhab_new_f, sep = "\n")
  }

  species_synhab_missing <- dplyr::filter(
    species_synhab, dplyr::if_all(dplyr::starts_with("hab"), isFALSE))
  if (nrow(species_synhab_missing) > 0) {
    species_synhab_missing_f <- fs::path(path_pa, "species_synhab_missing.txt")
    ecokit::cat_time(
      paste0(
        nrow(species_synhab_missing),
        " species do not exist in taxa standardization data"),
      level = 2L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0("see `", species_synhab_missing_f, "`` for more details"),
      level = 2L, cat_timestamp = FALSE)
    cat(
      sort(species_synhab_missing$taxon_name),
      file = species_synhab_missing_f, sep = "\n")
  }

  rm(
    taxa_list_distinct, species_synhab_new, species_synhab_missing,
    envir = environment())
  invisible(gc())

  species_taxa <- dplyr::arrange(taxa_list, ias_id) %>%
    dplyr::distinct(species_name, species_file)

  # # ..................................................................... ###

  # Species-specific data ------

  ecokit::cat_time("Species-specific data")
  .start_time_distribution <- lubridate::now(tzone = "CET")

  # # .................................... ###

  ## Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 1000L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # .................................... ###

  ## Species-specific data in parallel ----
  ecokit::cat_time("Species-specific data in parallel", level = 1L)

  ecokit::cat_time(
    paste0("There are ", nrow(species_taxa), " species to process"), level = 2L)

  species_pa_data <- future.apply::future_lapply(
    X = seq_len(nrow(species_taxa)),
    FUN = function(x) {

      # species file name
      sp_file <- species_taxa$species_file[x]
      sp_name <- species_taxa$species_name[x]

      if (length(sp_file) != 1 || is.na(sp_file) || !nzchar(sp_file) ||
          length(sp_name) != 1 || is.na(sp_name) || !nzchar(sp_name)) {
        ecokit::stop_ctx(
          "Species file name not unique or empty", sp_name = sp_name,
          sp_file = sp_file)
      }

      # output paths
      file_summary <- fs::path(
        paths_all$summary, paste0(sp_file, "_summary.RData"))
      file_pa <- fs::path(paths_all$pa_rdata, paste0(sp_file, "_pa.RData"))
      file_tif_all <- fs::path(paths_all$pa_tif, paste0(sp_file, "_all.tif"))
      file_tif_masked <- fs::path(
        paths_all$pa_tif, paste0(sp_file, "_masked.tif"))

      all_okay <- all(
        ecokit::check_data(file_summary, warning = FALSE),
        ecokit::check_data(file_pa, warning = FALSE),
        ecokit::check_tiff(file_tif_all, warning = FALSE),
        ecokit::check_tiff(file_tif_masked, warning = FALSE))

      if (all_okay) {
        return(file_summary)
      }

      # Maximum attempts
      max_attempts <- 8
      attempt <- 1

      repeat {

        # Run the distribution function
        species_data <- tryCatch(
          expr = {
            IASDT.R::naps_distribution(
              species = sp_name, env_file = env_file, verbose = FALSE)
          },
          error = function(e) {
            NULL
          },
          warning = function(w) {
            NULL
          })

        if (is.null(species_data) || !is.data.frame(species_data)) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            break
          }
          next
        }

        if (ncol(species_data) != 2 && nrow(species_data) != 1) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            break
          }
          next
        }

        Sys.sleep(0.5)

        # species without data
        if (is.na(species_data$pa_summary)) {
          break
        }

        # Check for the existence and validity of all files
        all_okay <- all(
          ecokit::check_data(file_summary, warning = FALSE),
          ecokit::check_data(file_pa, warning = FALSE),
          ecokit::check_tiff(file_tif_all, warning = FALSE),
          ecokit::check_tiff(file_tif_masked, warning = FALSE))

        if (all_okay) {
          break
        }

        attempt <- attempt + 1
        if (attempt > max_attempts) {
          break
        }
        invisible(gc())
      }

      species_data$pa_summary

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c(
      "env_file", "paths_all", "species_taxa", "update_citizen"))

  invisible(gc())

  # # .................................... ###

  ## Read processed species-specific data -----
  ecokit::cat_time("Read processed species-specific data", level = 1L)

  species_pa_data <- unlist(species_pa_data) %>%
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
  invisible(gc())

  # # .................................... ###

  ecokit::cat_time("Identify species with no data", level = 1L)

  # Species with no data
  sp_no_data <- setdiff(
    taxa_list$ias_id, as.integer(species_pa_data$ias_id))
  sp_no_data <- dplyr::filter(taxa_list, ias_id %in% sp_no_data)

  if (nrow(sp_no_data) > 0) {
    readr::write_tsv(
      x = sp_no_data, file = fs::path(path_pa, "species_no_data.csv"),
      col_names = TRUE)
    ecokit::cat_time(
      paste0(
        length(unique(species_pa_data$ias_id)),  " species were processed"),
      level = 2L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        nrow(sp_no_data), " species have no data. ",
        "See `species_no_data.csv` for details"),
      level = 2L, cat_timestamp = FALSE)
  } else {
    ecokit::cat_time(
      paste0("All ", length(taxa_list$species_name), " species were processed"),
      level = 2L)
  }

  # # .................................... ###

  ecokit::cat_diff(
    init_time = .start_time_distribution,
    prefix = "Processing species-specific data took ",
    msg_n_lines = 1, level = 2L)

  rm(sp_no_data, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Merge species summary info -----
  ecokit::cat_time("Merge species summary info")

  unnest_columns <- c(
    "bioreg_data", "bioreg_data_mask", "species_country",
    "species_country_masked")

  species_pa_data <- species_pa_data %>%
    # Split # grid cell per country / biogeographical region as separate column
    tidyr::unnest_wider(tidyselect::all_of(unnest_columns)) %>%
    dplyr::mutate(
      # replace NA for biogeographical regions and country analysis with 0
      dplyr::across(
        tidyselect::matches("^cnt_|^bioreg_"),
        ~ dplyr::if_else(is.na(.x), 0L, .x))) %>%
    dplyr::rename(species_name = species) %>%
    dplyr::left_join(
      species_synhab, by = c("ias_id", "species_name", "taxon_name"))

  rm(taxa_list, envir = environment())

  ## Save summary results -----
  ecokit::cat_time("Save summary results", level = 1L)

  ecokit::cat_time("`qs2`", level = 2L)
  save(species_pa_data, file = fs::path(path_pa, "species_pa_data.qs2"))

  ecokit::cat_time("Summary data without maps", level = 1L)
  ecokit::cat_time("csv format", level = 2L)
  sp_pa_summary_df <- species_pa_data %>%
    dplyr::select(tidyselect::where(~ !is.list(.x)))

  readr::write_excel_csv(
    x = sp_pa_summary_df, file = fs::path(path_pa, "sp_pa_summary_df.csv"),
    progress = FALSE)

  ecokit::cat_time("RData format", level = 2L)
  save(sp_pa_summary_df, file = fs::path(path_pa, "sp_pa_summary_df.RData"))

  # only keep selected columns for summary plots
  selected_cols <- c(
    "ias_id", "species_name", "n_cells_all", "n_cells_naturalized",
    "n_cells_cz_excl", "n_cells_cz_excl_masked", "pa_map", "pa_masked_map",
    "mask_keep", "pa_cz_excluded", "pa_cz_excluded_masked")
  species_pa_data <- species_pa_data %>%
    dplyr::select(
      tidyselect::all_of(selected_cols), tidyselect::starts_with("hab_"))

  rm(sp_pa_summary_df, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Summary of merged data ------
  ecokit::cat_time("Summary of merged data")

  ## Number of NAPS per grid cell -----
  ecokit::cat_time("# NAPS per grid cell", level = 1L)

  naps_n_sp <- dplyr::filter(species_pa_data, n_cells_all > 0) %>%
    dplyr::pull(pa_map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("naps_n_sp") %>%
    terra::toMemory()

  species_pa_data <- dplyr::select(species_pa_data, -pa_map)
  invisible(gc())

  # save as RData
  ecokit::save_as(
    object = terra::wrap(naps_n_sp), object_name = "naps_n_sp",
    out_path = fs::path(path_pa, "naps_n_sp.RData"))

  # save as tif
  raster::writeRaster(
    x = naps_n_sp, overwrite = TRUE,
    filename = fs::path(path_pa, "naps_n_sp.tif"))

  # # .................................... ###

  ## Number of NAPS per grid cell - masked data -----

  naps_n_sp_masked <- species_pa_data %>%
    dplyr::filter(n_cells_naturalized > 0) %>%
    dplyr::pull(pa_masked_map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("naps_n_sp_masked") %>%
    terra::toMemory()

  species_pa_data <- dplyr::select(species_pa_data, -pa_masked_map)
  invisible(gc())

  # save as RData
  ecokit::save_as(
    object = terra::wrap(naps_n_sp_masked), object_name = "naps_n_sp_masked",
    out_path = fs::path(path_pa, "naps_n_sp_masked.RData"))

  # save as tif
  raster::writeRaster(
    x = naps_n_sp_masked, overwrite = TRUE,
    filename = fs::path(path_pa, "naps_n_sp_masked.tif"))

  # # .................................... ###

  ## Number of excluded citizen science grid cells -----

  naps_n_sp_cz <- species_pa_data %>%
    dplyr::pull(pa_cz_excluded) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("naps_n_sp_cz") %>%
    terra::toMemory()
  ecokit::save_as(
    object = terra::wrap(naps_n_sp_cz), object_name = "naps_n_sp_cz",
    out_path = fs::path(path_pa, "naps_n_sp_cz.RData"))
  raster::writeRaster(
    x = naps_n_sp_cz, overwrite = TRUE,
    filename = fs::path(path_pa, "naps_n_sp_cz.tif"))

  naps_n_sp_cz_masked <- species_pa_data %>%
    dplyr::pull(pa_cz_excluded_masked) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("naps_n_sp_cz_masked") %>%
    terra::toMemory()
  ecokit::save_as(
    object = terra::wrap(naps_n_sp_cz_masked),
    object_name = "naps_n_sp_cz_masked",
    out_path = fs::path(path_pa, "naps_n_sp_cz_masked.RData"))
  raster::writeRaster(
    x = naps_n_sp_cz_masked, overwrite = TRUE,
    filename = fs::path(path_pa, "naps_n_sp_cz_masked.tif"))

  species_pa_data <- dplyr::select(
    species_pa_data, -pa_cz_excluded, -pa_cz_excluded_masked)
  invisible(gc())

  # # .................................... ###

  ## Plotting summary of NAPS data -----
  ecokit::cat_time("Plotting summary of NAPS data", level = 1L)

  eu_bound <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")
  limit_x <- c(2600000, 6550000)
  limit_y <- c(1450000, 5410000)

  # # +++++++++++++++++++++++++++++++++ ###

  plot_theme <- ggplot2::theme_bw() +
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

  # Function to create all NAPS plots using the `create_naps_plot` function

  create_naps_plot <- function(
    rast, filename, p_title, p_subtitle,
    p_title_log, p_subtitle_log) {

    plot_n_sp <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        eu_bound, mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(
        data = terra::classify(rast, cbind(0, NA)), maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = ecokit::integer_breaks()) +
      ggplot2::geom_sf(
        eu_bound, mapping = ggplot2::aes(), color = "grey40",
        linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = limit_x) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = limit_y) +
      ggplot2::labs(
        title = p_title, subtitle = p_subtitle, fill = "# NAPS") +
      plot_theme

    plot_n_sp_log <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        eu_bound, mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(
        data = log10(terra::classify(rast, cbind(0, NA))), maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = ecokit::integer_breaks()) +
      ggplot2::geom_sf(
        eu_bound, mapping = ggplot2::aes(), color = "grey40",
        linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = limit_x) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = limit_y) +
      ggplot2::labs(
        title = p_title_log, subtitle = p_subtitle_log,
        fill = "# NAPS\n(log10)") +
      plot_theme

    plot <- cowplot::plot_grid(
      plot_n_sp, plot_n_sp_log, ncol = 2, nrow = 1) +
      cowplot::draw_label(
        label = last_update, color = "grey", x = 0.99, y = 0.98,
        size = 7, hjust = 1)

    ragg::agg_jpeg(
      filename = fs::path(path_pa, filename),
      width = 30, height = 15.5, res = 600, quality = 100, units = "cm")
    print(plot)
    grDevices::dev.off()

    rm(rast, plot, envir = environment())
    invisible(gc())
  }

  create_naps_plot(
    rast = naps_n_sp,
    filename = "naps_n_sp.jpeg",
    p_title = "Number of NAPS per 10\u00D710 km grid",
    p_subtitle = "All data (including cultivated or casual observations)",
    p_title_log = "Number of NAPS per 10\u00D710 km grid (log10 scale)",
    p_subtitle_log = "All data (including cultivated or casual observations)")

  # Masked plots (excluding cultivated or casual observations)
  create_naps_plot(
    rast = naps_n_sp_masked,
    filename = "naps_n_sp_masked.jpeg",
    p_title = "Number of NAPS per 10\u00D710 km grid",
    p_subtitle = "Excluding cultivated or casual observations",
    p_title_log = "Number of NAPS per 10\u00D710 km grid (log10 scale)",
    p_subtitle_log = "Excluding cultivated or casual observations")

  # Excluded NAPS (citizen science data) plots
  create_naps_plot(
    rast = naps_n_sp_cz,
    filename = "naps_n_sp_cz.jpeg",
    p_title =
      "Number of excluded NAPS (citizen science data) per 10\u00D710 km grid",
    p_subtitle = "All data (including cultivated or casual observations)",
    p_title_log = paste0(
      "Number of excluded NAPS (citizen science data) per 10\u00D710 km",
      "grid (log10 scale)"),
    p_subtitle_log = "All data (including cultivated or casual observations)")

  # Excluded NAPS masked (citizen science data, excluding cultivated/casual)
  # plots
  create_naps_plot(
    rast = naps_n_sp_cz_masked,
    filename = "naps_n_sp_cz_masked.jpeg",
    p_title =
      "Number of excluded NAPS (citizen science data) per 10\u00D710 km grid",
    p_subtitle = "Excluding cultivated or casual observations",
    p_title_log = paste0(
      "Number of excluded NAPS (citizen science data) per 10\u00D710 km",
      "grid (log10 scale)"),
    p_subtitle_log = "Excluding cultivated or casual observations")

  # # +++++++++++++++++++++++++++++++++ ###

  rm(
    eu_bound, limit_x, limit_y, naps_n_sp,  naps_n_sp_masked, naps_n_sp_cz,
    naps_n_sp_cz_masked, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## # Number of NAPS to be used in the model - per habitat type ------

  threshold_theme <- ggplot2::theme(
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

  n_sp_data <- dplyr::select(species_pa_data, species_name, n_cells_all)
  n_sp_data <- purrr::map_dfr(
    .x = sort(unique(n_sp_data$n_cells_all)),
    .f = ~ tibble::tibble(
      threshold = .x,
      n_sp = sum(n_sp_data$n_cells_all >= .x, na.rm = TRUE))) %>%
    dplyr::filter(threshold < 450)

  n_sp_data_masked <- dplyr::select(
    species_pa_data, species_name, n_cells_naturalized)
  n_sp_data_masked <- purrr::map_dfr(
    .x = sort(unique(n_sp_data_masked$n_cells_naturalized)),
    .f = ~ tibble::tibble(
      threshold = .x,
      n_sp = sum(n_sp_data_masked$n_cells_naturalized >= .x, na.rm = TRUE))) %>%
    dplyr::filter(threshold < 450)

  # # +++++++++++++++++++++++++++++++++ ###

  selected_hab <- c(
    "hab_1_forests", "hab_2_open_forests", "hab_3_scrub",
    "hab_4a_natural_grasslands", "hab_4b_human_maintained_grasslands",
    "hab_10_wetland", "hab_12a_ruderal_habitats",
    "hab_12b_agricultural_habitats")

  n_sp_hab_data <- purrr::map_dfr(
    .x = selected_hab,
    .f = function(hab) {

      # nolint start
      curr_data <- dplyr::filter(species_pa_data, !!as.symbol(hab))
      # nolint end

      purrr::map_dfr(
        sort(unique(species_pa_data$n_cells_all)),
        ~ tibble::tibble(
          hab = hab, threshold = .x,
          n_sp = sum(curr_data$n_cells_all >= .x, na.rm = TRUE)))}) %>%
    dplyr::mutate(
      hab = stringr::str_remove(hab, "^hab_"),
      hab = stringr::str_replace(hab, "_", " - "),
      hab = stringr::str_replace_all(hab, "_", " "),
      hab = forcats::fct_inorder(hab)) %>%
    dplyr::slice(gtools::mixedorder(hab)) %>%
    dplyr::filter(threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = threshold, y = n_sp, group = hab, colour = hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = n_sp_data, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = threshold, y = n_sp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(
      title = "All data (including cultivated or casual observations)") +
    threshold_theme +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1, label.position = "bottom",
        override.aes = list(linewidth = 1.5)))

  n_sp_hab_masked_data <- purrr::map_dfr(
    .x = selected_hab,
    .f = function(hab) {

      # nolint start
      curr_data <- dplyr::filter(species_pa_data, !!as.symbol(hab))
      # nolint end

      purrr::map_dfr(
        sort(unique(species_pa_data$n_cells_naturalized)),
        ~ tibble::tibble(
          hab = hab, threshold = .x,
          n_sp = sum(curr_data$n_cells_naturalized >= .x, na.rm = TRUE)))
    }) %>%
    dplyr::mutate(
      hab = stringr::str_remove(hab, "^hab_"),
      hab = stringr::str_replace(hab, "_", " - "),
      hab = stringr::str_replace_all(hab, "_", " "),
      hab = forcats::fct_inorder(hab)) %>%
    dplyr::slice(gtools::mixedorder(hab)) %>%
    dplyr::filter(threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = threshold, y = n_sp, group = hab, colour = hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = n_sp_data_masked, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = threshold, y = n_sp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(title = "Excluding cultivated or casual observations") +
    threshold_theme +
    ggplot2::theme(legend.position = "none")

  plot_legend <- ggpubr::as_ggplot(ggpubr::get_legend(n_sp_hab_data))
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "Number of NAPS to be used in the models based on the arbitrary ",
        "selection of # of presence grid cells per species"),
      fontface = "bold", colour = "blue")

  plot <- cowplot::plot_grid(
    plot_title,
    cowplot::plot_grid(
      (n_sp_hab_data + ggplot2::theme(legend.position = "none")),
      n_sp_hab_masked_data,
      ncol = 2, nrow = 1),
    plot_legend,
    ncol = 1, rel_heights = c(0.05, 1, 0.05))

  ragg::agg_jpeg(
    filename = fs::path(path_pa, "naps_n_sp_threshold_hab.jpeg.jpeg"),
    width = 30, height = 17, res = 600, quality = 100, units = "cm")
  print(plot)
  grDevices::dev.off()

  rm(
    species_pa_data, n_sp_hab_data, n_sp_hab_masked_data,
    plot_legend, plot_title, plot, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting species distribution -----
  ecokit::cat_time("Plotting species distribution")
  .start_time_plot <- lubridate::now(tzone = "CET")

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

  species_plots <- future.apply::future_lapply(
    X = seq_len(nrow(species_taxa)),
    FUN = function(x) {

      # Maximum attempts
      max_attempts <- 4
      attempt <- 1

      repeat {
        sp_plot <- tryCatch(
          IASDT.R::naps_plot(
            species = species_taxa$species_name[x], env_file = env_file))

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
        species = species_taxa$species_name[x],
        img_file = sp_plot, img_valid = img_valid)

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "species_taxa"))


  # validate that all plots were created
  species_plots <- dplyr::bind_rows(species_plots)
  invalid_plots <- dplyr::filter(species_plots, !img_valid)
  valid_plots <- dplyr::filter(species_plots, img_valid)

  if (nrow(invalid_plots) > 0) {
    readr::write_tsv(
      x = invalid_plots, file = fs::path(path_pa, "species_invalid_plots.csv"),
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

  rm(species_plots, valid_plots, invalid_plots, envir = environment())

  # # .................................... ###

  ## Stopping cluster ----

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }
  # # .................................... ###

  ecokit::cat_diff(
    init_time = .start_time_plot,
    prefix = "Plotting species distribution took ",
    msg_n_lines = 1, level = 2L)

  # # ..................................................................... ###

  # Function Summary ----
  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing species data was finished in ", ... = "\n")

  return(invisible(NULL))
}
