#' Process GBIF occurrence data for the `IASDT`
#'
#' Extracts, processes, and visualises occurrence data from the [Global
#' Biodiversity Information Facility (GBIF)](https://www.gbif.org) for the
#' Invasive Alien Species Digital Twin (`IASDT`). Orchestrated by
#' `gbif_process()`, it requests, downloads, cleans, chunks, and maps species
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
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param delete_chunks Logical. If `TRUE` (default), deletes chunk files.
#' @param chunk_file Character. Path of chunk file for processing.
#'
#' @param max_uncertainty Numeric. Maximum spatial uncertainty in kilometres.
#'   Default: `10`.
#' @param start_year Integer. Earliest collection year to be included. Default
#'   is 1981.
#' @param save_rdata Logical. If `TRUE` (default), saves chunk data as `.RData`.
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
#' - **`gbif_process()`**: Orchestrates GBIF data requests, downloads,
#'   processing, and mapping. Saves `RData`, Excel, and JPEG summary files.
#' - **`gbif_download()`**: Requests and downloads GBIF data (if `download =
#'   TRUE`), using the specified criteria (taxa, coordinates, time period, and
#'   boundaries), splits into small chunks (if `split_chunks = TRUE`), and saves
#'   metadata. Returns `invisible(NULL)`.
#' - **`gbif_read_chunk()`**: Filters chunk data (spatial/temporal, e.g.,
#'   spatial uncertainty, collection year, coordinate precision, and taxonomic
#'   rank), select relevant columns, and saves as `.RData` (if `save_rdata =
#'   TRUE`) or returns it (if `return_data = TRUE`). Skips if `.RData` exists
#'   and `overwrite = FALSE`.
#' - **`gbif_species_data()`**: Converts species-specific data to `sf` and
#'   raster formats, generating distribution maps.
#' @note Relies on a static RDS file listing IAS species, GBIF keys, and
#'   metadata, standardized by Marina Golivets (Feb 2024).

# # |------------------------------------------------------------------------| #
# gbif_process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name gbif_data
#' @rdname gbif_data
#' @order 1

gbif_process <- function(
    env_file = ".env", r_environ = ".Renviron", n_cores = 6L,
    strategy = "multisession", request = TRUE, download = TRUE,
    split_chunks = TRUE, overwrite = FALSE, delete_chunks = TRUE,
    chunk_size = 50000L, boundaries = c(-30, 50, 25, 75), start_year = 1981L) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----

  ecokit::cat_time("Checking arguments")
  ecokit::check_args(args_to_check = "r_environ", args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "request", "download", "split_chunks", "overwrite", "delete_chunks"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c("boundaries", "start_year", "chunk_size"),
    args_type = "numeric", arg_length = c(4L, 1L, 1L))

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
    future.seed = TRUE)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_gbif <- path_gbif_interim <- species_name <- n_obs <- publisher_type <-
    taxa_info <- species <- CellCode <- ias_id <- taxon_name <- others <-
    species_name2 <- species_file <- path_grid <- summ_map <- summ_map_gg <-
    plot_title <- legend_label <- n <- eu_borders <- publisher <-
    Longitude_3035 <- Latitude_3035 <- sp_data <- partner_map <-
    partner_id <- citizen_science <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_gbif", "DP_R_gbif_processed", FALSE, FALSE,
    "path_gbif_interim", "DP_R_gbif_interim", FALSE, FALSE,
    "eu_borders", "DP_R_country_boundaries", FALSE, TRUE,
    "taxa_info", "DP_R_taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "IASDT.R", "purrr", "tibble", "terra", "tidyr", "dplyr", "readr",
      "stringr", "bit64", "tidyselect", "fs", "sf", "ecokit", "paletteer",
      "ggplot2", "tidyterra", "magrittr", "ragg", "cowplot", "grid", "scales",
      "viridis"),
    strategy = strategy)

  # # ..................................................................... ###

  # request / download GBIF data and split data into chunks -----

  IASDT.R::gbif_download(
    env_file = env_file, r_environ = r_environ, request = request,
    download = download, split_chunks = split_chunks, chunk_size = chunk_size,
    boundaries = boundaries, start_year = start_year)

  # # ..................................................................... ###

  gbif_metadata <- fs::path(path_gbif, "gbif_metadata.RData")
  if (!file.exists(gbif_metadata)) {
    ecokit::stop_ctx(
      "GBIF metadata file does not exist", gbif_metadata = gbif_metadata,
      include_backtrace = TRUE)
  }
  gbif_metadata <- ecokit::load_as(gbif_metadata)

  taxa_list <- ecokit::load_as(taxa_info)
  path_sp_data <- fs::path(path_gbif, "species_data")
  path_publishers <- fs::path(path_gbif, "n_species_per_publisher")
  fs::dir_create(c(path_sp_data, path_publishers))

  # grid_10_land_sf
  grid_sf <- fs::path(path_grid, "grid_10_land_sf.RData")
  if (!file.exists(grid_sf)) {
    ecokit::stop_ctx(
      "Reference grid (sf) file not found", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }
  grid_sf <- ecokit::load_as(grid_sf)

  # grid_10_land_crop
  grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Reference grid file not found", grid_r = grid_r,
      include_backtrace = TRUE)
  }
  grid_r <- ecokit::load_as(grid_r, unwrap_r = TRUE)

  eu_borders <- ecokit::load_as(eu_borders) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  invisible(gc())

  # # ..................................................................... ###

  # Processing data chunks -----

  ecokit::cat_time("\nProcessing data chunks")
  .start_time_chunks <- lubridate::now(tzone = "CET")

  chunk_list <- list.files(
    path = path_gbif_interim, pattern = "Chunk_.+.txt", full.names = TRUE)
  chunk_list_rdata <- stringr::str_replace_all(chunk_list, ".txt$", ".RData")

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time(
    "Processing chunks in parallel, save each as RData files", level = 1L)

  gbif_data <- future.apply::future_lapply(
    X = chunk_list,
    FUN = function(x) {
      IASDT.R::gbif_read_chunk(
        chunk_file = x, env_file = env_file, save_rdata = TRUE,
        return_data = FALSE, overwrite = overwrite)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("env_file", "overwrite"))

  ecokit::cat_diff(
    init_time = .start_time_chunks, prefix = "Finished in ", level = 2L)

  ecokit::cat_time("Reading processed chunks into a single dataset", level = 1L)

  if (all(file.exists(chunk_list_rdata))) {
    gbif_data <- future.apply::future_lapply(
      X = chunk_list_rdata, FUN = ecokit::load_as,
      future.scheduling = Inf, future.seed = TRUE) %>%
      dplyr::bind_rows() %>%
      # merge with taxa standardization results
      dplyr::left_join(taxa_list, by = "speciesKey") %>%
      # arrange by species name
      dplyr::arrange(species_name)
  } else {
    stop_message <- chunk_list[which(!file.exists(chunk_list_rdata))] %>%
      paste0(" >> ", ., collapse = "\n") %>%
      paste0("The following chunks were not processed\n", .)
    ecokit::stop_ctx(stop_message, include_backtrace = TRUE)
  }

  ecokit::cat_diff(
    init_time = .start_time_chunks,
    prefix = "Processing data chunks was finished in ", level = 2L)

  ecokit::cat_time(
    paste0(
      "A total of ", format(nrow(gbif_data), big.mark = ","), " observations"),
    level = 2L, cat_timestamp = FALSE)


  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  ecokit::cat_time("Saving `gbif_data` to disk", level = 1L)
  ecokit::save_as(
    object = gbif_data, out_path = fs::path(path_gbif, "gbif_data.qs2"),
    n_threads = 5L)

  rm(taxa_list, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Citizen science unique grids ----

  ecokit::cat_time("Citizen science unique grids")

  ecokit::cat_time("Prepare input data for summarising", level = 1L)
  citizen_others <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::distinct(species, CellCode, publisher_type)

  ecokit::cat_time(
    "Number of unique citizen science grid cells per species", level = 1L)
  citizen_unique <- citizen_others %>%
    dplyr::group_by(species, CellCode) %>%
    tidyr::pivot_wider(
      names_from = publisher_type, values_from = publisher_type,
      values_fn = unique) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(others) & citizen_science == "citizen_science") %>%
    dplyr::count(species) %>%
    dplyr::rename(citizen_science_unique = n)

  ecokit::cat_time(
    "Number of grid cells for citizen science and other data sources",
    level = 1L)
  citizen_count <- dplyr::count(citizen_others, species, publisher_type) %>%
    tidyr::pivot_wider(names_from = publisher_type, values_from = n) %>%
    dplyr::select(species, citizen_science, others) %>%
    dplyr::left_join(citizen_unique, by = "species") %>%
    dplyr::select(species, others, tidyselect::everything())

  ecokit::cat_time("Save citizen science summary data", level = 1L)
  save(citizen_count, file = fs::path(path_gbif, "citizen_count.RData"))

  rm(citizen_others, citizen_unique, citizen_count, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Covert GBIF data to presence-only grids ----

  ecokit::cat_time("Covert GBIF data to presence-only grids")
  gbif_grid <- sf::st_drop_geometry(gbif_data) %>%
    # only unique combinations between species and grid ID
    dplyr::distinct(
      ias_id, taxon_name, species_name, species_name2,
      species_file, CellCode, publisher_type) %>%
    # add polygons for grids
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    # ensure that the data is sf object
    sf::st_as_sf()

  ecokit::cat_time("Save presence-only grids", level = 1L)
  save(gbif_grid, file = fs::path(path_gbif, "gbif_grid.RData"))
  invisible(gc())

  # # ..................................................................... ###

  # Number of observations or grid cells per species ----

  ecokit::cat_time("# observations or grid cells per species")

  ecokit::cat_time("# observations per species", level = 1L)
  species_n_obs <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::count(ias_id, taxon_name, species_name) %>%
    dplyr::rename(gbif_n_obs = n)

  ecokit::cat_time("# grids per species", level = 1L)
  species_n_grids <- dplyr::select(gbif_grid, -publisher_type) %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(ias_id, taxon_name, species_name, CellCode) %>%
    dplyr::count(ias_id, taxon_name, species_name) %>%
    dplyr::rename(gbif_n_grids = n)

  ecokit::cat_time(
    "Merge number of observations and grids per species", level = 1L)
  gbif_n_obs_n_grid <- dplyr::full_join(
    species_n_obs, species_n_grids,
    by = c("ias_id", "taxon_name", "species_name"))

  ecokit::cat_time("Save as RData", level = 1L)
  save(gbif_n_obs_n_grid, file = fs::path(path_gbif, "gbif_n_obs_n_grid.RData"))

  ecokit::cat_time("Save as xlsx", level = 1L)
  writexl::write_xlsx(
    x = gbif_n_obs_n_grid, path = fs::path(path_gbif, "gbif_n_obs_n_grid.xlsx"))

  rm(gbif_n_obs_n_grid, species_n_obs, species_n_grids, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Summary maps ----
  ecokit::cat_time("Summary maps")

  ## Number of observations per grid  ----

  ecokit::cat_time("Number of observations per grid cell", level = 1L)
  gbif_n_obs <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = grid_r, field = "n") %>%
    terra::mask(grid_r) %>%
    stats::setNames("NumObs") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()
  gbif_n_obs_log <- terra::unwrap(gbif_n_obs) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time(
    "Number of citizen science observations per grid cell", level = 1L)
  gbif_cz_n_obs <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::filter(publisher_type == "citizen_science") %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = grid_r, field = "n") %>%
    terra::mask(grid_r) %>%
    stats::setNames("NumObs") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()
  gbif_cz_n_obs_log <- terra::unwrap(gbif_cz_n_obs) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Save as RData", level = 2L)
  save(
    gbif_n_obs, gbif_n_obs_log, gbif_cz_n_obs, gbif_cz_n_obs_log,
    file = fs::path(path_gbif, "gbif_n_obs.RData"))

  invisible(gc())

  # # ................................... ###

  ## Number of species per grid cell -----

  ecokit::cat_time("Number of NAPS per grid cell", level = 1L)
  gbif_n_species <- sf::st_drop_geometry(gbif_grid) %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = grid_r, field = "n") %>%
    terra::mask(grid_r) %>%
    stats::setNames("n_sp") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()
  gbif_n_species_log <- terra::unwrap(gbif_n_species) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Number of CZ NAPS per grid cell", level = 1L)
  gbif_cz_n_species <- sf::st_drop_geometry(gbif_grid) %>%
    dplyr::filter(publisher_type == "citizen_science") %>%
    dplyr::count(CellCode) %>%
    dplyr::left_join(grid_sf, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::vect() %>%
    terra::rasterize(y = grid_r, field = "n") %>%
    terra::mask(grid_r) %>%
    stats::setNames("n_sp") %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()
  gbif_cz_n_species_log <- terra::unwrap(gbif_cz_n_species) %>%
    log10() %>%
    ecokit::set_raster_crs(crs = "epsg:3035") %>%
    terra::wrap()

  ecokit::cat_time("Save as RData", level = 2L)
  save(
    gbif_n_species, gbif_n_species_log,
    gbif_cz_n_species, gbif_cz_n_species_log,
    file = fs::path(path_gbif, "gbif_n_species.RData"))

  rm(grid_sf, envir = environment())
  invisible(gc())

  # # ................................... ###

  ## plot_gbif_summary -----

  ecokit::cat_time("Plotting summary maps", level = 1L)

  gbif_date <- gbif_metadata$gbif_status$modified %>%
    lubridate::date() %>%
    format("%d %B %Y")

  gbif_doi <- gbif_metadata$gbif_status$doi
  plot_tag <- format(Sys.Date(), "%d %B %Y")
  plot_tag <- stringr::str_glue(
    "DOI: {gbif_doi} ({gbif_date})\nLast update: {plot_tag}")

  rm(gbif_date, gbif_doi, gbif_grid, gbif_metadata, envir = environment())

  # Plotting limits
  x_limit <- c(2600000, 6700000)
  y_limit <- c(1450000, 5420000)

  plotting_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0.1, 0, 0.1, "cm"),
      plot.title = ggplot2::element_text(
        size = 12, color = "blue", face = "bold",
        hjust = 0, margin = ggplot2::margin(0, 0, 0, 0)),
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


  plot_gbif_summary <- function(
    r_map, plot_title, legend_label = NULL, eu_map = eu_borders,
    plot_theme = plotting_theme, xlim = x_limit, ylim = y_limit) {

    out_map <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        eu_map, mapping = ggplot2::aes(), color = "grey30", linewidth = 0.1,
        fill = "grey95", inherit.aes = TRUE) +
      tidyterra::geom_spatraster(data = terra::trim(r_map), maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma",
        breaks = ecokit::integer_breaks()) +
      ggplot2::geom_sf(
        eu_map, mapping = ggplot2::aes(), color = "grey40", linewidth = 0.075,
        fill = "transparent", inherit.aes = TRUE) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = xlim) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = ylim) +
      plot_theme +
      ggplot2::labs(title = plot_title, fill = legend_label)

    return(out_map)
  }

  # a common main title of the figure
  main_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      label = "Summary of naturalised alien plant species (NAPS) data (GBIF)",
      fontface = "bold", hjust = 0.5) +
    cowplot::draw_label(
      label = plot_tag, fontface = "italic", color = "grey", x = 0.98,
      y = 0.35, size = 8, hjust = 1) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0.4, 0))

  # Plotting summary maps
  gbif_summary_plot <- tibble::tibble(
    # tibble for plotting information
    summ_map = list(
      terra::unwrap(gbif_n_obs) / 1000, terra::unwrap(gbif_n_obs_log),
      terra::unwrap(gbif_n_species) / 100, terra::unwrap(gbif_n_species_log)),
    plot_title = c(
      "Number of observations", "Number of observations (log10)",
      "Number of NAPS", "Number of NAPS (log10)"),
    legend_label = c("\u00D7 1000", "log10", "\u00D7 100", "Log10")) %>%
    # add the ggplot object as column to the data
    dplyr::mutate(
      summ_map_gg = purrr::pmap(
        .l = list(summ_map, plot_title, legend_label),
        .f = function(summ_map, plot_title, legend_label) {
          plot_gbif_summary(
            r_map = summ_map, plot_title = plot_title,
            legend_label = legend_label)
        })) %>%
    dplyr::pull(summ_map_gg) %>%
    # Plot the four panels together on a single figure
    cowplot::plot_grid(plotlist = ., ncol = 2, nrow = 2) %>%
    # add the common title
    cowplot::plot_grid(main_title, ., ncol = 1, rel_heights = c(0.035, 1))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_gbif, "gbif_summary.jpeg"),
    width = 25, height = 25.8, res = 600, quality = 100, units = "cm")
  print(gbif_summary_plot)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Species-specific data ----
  ecokit::cat_time("Species-specific data")

  species_list <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::distinct(species_name) %>%
    dplyr::pull(species_name) %>%
    sort()

  # Species data --- sf ----
  ecokit::cat_time("Split species data - sf", level = 1L)

  # On LUMI, extracting species data in parallel using `furrr::future_walk` took
  # much longer time than working sequentially `purrr::walk` due to the
  # existence of the very large `gbif_data` object. Although it is possible to
  # integrate following chunk in the `gbif_species_data` function, it would be
  # much faster to implement the following sequentially, then process species
  # maps in parallel later on, after the deletion of the `gbif_data` object
  purrr::walk(
    .x = species_list,
    .f = ~ {
      # Filter data on this species
      sp_data <- dplyr::filter(gbif_data, species_name == .x)

      # Save if there is data
      if (nrow(sp_data) > 0) {
        sp_name <- ecokit::replace_space(.x) %>%
          # replace non-ascii multiplication symbol with x
          stringr::str_replace_all("\u00D7", "x") %>%
          stringr::str_replace_all("-", "")
        outfile_sf <- fs::path(path_sp_data, paste0(sp_name, ".RData"))
        ecokit::save_as(
          object = sp_data, object_name = sp_name, out_path = outfile_sf)
      }
      return(invisible(NULL))
    })

  # # ..................................................................... ###

  # Number of species per publisher --------
  ecokit::cat_time("Number of species per publisher")

  n_species_publisher <- sf::st_drop_geometry(gbif_data) %>%
    dplyr::distinct(
      ias_id, publisher, CellCode, Longitude_3035, Latitude_3035) %>%
    tidyr::nest(sp_data = -publisher) %>%
    dplyr::mutate(n_obs = purrr::map_int(sp_data, nrow)) %>%
    dplyr::arrange(dplyr::desc(n_obs)) %>%
    dplyr::mutate(
      partner_id = dplyr::row_number(),
      partner_map = purrr::map(
        .x = sp_data,
        .f = ~{
          dplyr::summarise(
            .x, n_species = length(unique(ias_id)),
            .by = c(CellCode, Longitude_3035, Latitude_3035)) %>%
            sf::st_as_sf(
              coords = c("Longitude_3035", "Latitude_3035"), crs = 3035) %>%
            dplyr::select(-CellCode) %>%
            terra::rasterize(grid_r, field = "n_species")
        })) %>%
    dplyr::select(-sp_data, -n_obs)
  invisible(gc())

  n_species_publisher <- n_species_publisher %>%
    dplyr::mutate(
      partner_map_gg = purrr::pmap(
        .l = list(publisher, partner_map, partner_id),
        .f = function(publisher, partner_map, partner_id) {

          partner_id2 <- stringr::str_pad(partner_id, width = 3, pad = "0")

          publisher_plot <- ggplot2::ggplot() +
            ggplot2::geom_sf(
              eu_borders, mapping = ggplot2::aes(), color = "grey30",
              linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
            tidyterra::geom_spatraster(
              data = terra::trim(partner_map), maxcell = Inf) +
            paletteer::scale_fill_paletteer_c(
              na.value = "transparent", "viridis::plasma",
              breaks = ecokit::integer_breaks()) +
            ggplot2::geom_sf(
              eu_borders, mapping = ggplot2::aes(), color = "grey40",
              linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
            ggplot2::scale_x_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = x_limit) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0, 0)), limits = y_limit) +
            plotting_theme +
            ggplot2::labs(title = publisher, fill = "n_species")

          ragg::agg_jpeg(
            filename = fs::path(
              path_publishers,
              paste("n_sp_per_publisher_", partner_id2, ".jpeg")),
            width = 25, height = 25, res = 600, quality = 100, units = "cm")
          print(publisher_plot)
          grDevices::dev.off()

          invisible(return(NULL))
        }))

  n_species_publisher <- n_species_publisher %>%
    dplyr::mutate(
      partner_map_gg = NULL, partner_id = NULL,
      partner_map = purrr::map(partner_map, terra::wrap))
  ecokit::save_as(
    object = n_species_publisher, object_name = "n_species_publisher",
    out_path = fs::path(path_gbif, "n_species_publisher.RData"))

  rm(
    gbif_data, gbif_n_obs, gbif_n_obs_log, gbif_n_species, gbif_n_species_log,
    n_species_publisher, grid_r, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Grid / raster / plotting ----
  ecokit::cat_time("Split species data - grid + raster + plot", level = 1L)

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 2L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time("Splitting species data in parallel", level = 2L)

  sp_data_tmp <- ecokit::quietly(
    expr = future.apply::future_lapply(
      X = species_list,
      FUN = IASDT.R::gbif_species_data, env_file = env_file,
      verbose = FALSE, plot_tag = plot_tag,
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c("env_file", "plot_tag")),
    "Ignoring empty aesthetics")
  rm(sp_data_tmp, envir = environment())

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Clean up chunk files ----
  if (delete_chunks) {
    ecokit::cat_time("Clean up - remove temporary chunk files")
    list.files(path_gbif_interim, full.names = TRUE) %>%
      fs::file_delete()
    fs::dir_delete(path_gbif_interim)
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nProcessing GBIF data was finished in ")

  return(invisible(NULL))
}
