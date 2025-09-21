# # |------------------------------------------------------------------------| #
# naps_distribution ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name naps_data
#' @rdname naps_data
#' @order 2

naps_distribution <- function(
    species = NULL, env_file = ".env", verbose = FALSE, dist_citizen = 100L) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")
  ecokit::info_chunk(species, verbose = verbose)

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::check_args(args_to_check = "verbose", args_type = "logical")

  if (is.null(species) || is.na(species) || species == "") {
    ecokit::stop_ctx(
      "`species` cannot be empty", species = species, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid <- path_taxa_info <- path_gbif <- path_elter <- path_easin <-
    path_bioreg <- path_grid_ref <- path_pa <- biogreg <- pa_masked <-
    country <- species_name2 <- n <- gbif <- easin <- elter <- species_name <-
    pa <- `status-decision` <- path_taxa_country <- country <- taxon_name <-
    status_decision <- path_taxa_info_rdata <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables", verbose = verbose)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_grid_ref", "DP_R_grid_raw", TRUE, FALSE,
    "path_gbif", "DP_R_gbif_processed", TRUE, FALSE,
    "path_easin", "DP_R_easin_processed", TRUE, FALSE,
    "path_elter", "DP_R_elter_processed", FALSE, TRUE,
    "path_pa", "DP_R_pa", FALSE, FALSE,
    "path_taxa_country", "DP_R_taxa_country", FALSE, TRUE,
    "path_taxa_info_rdata", "DP_R_taxa_info_rdata", FALSE, TRUE,
    "path_taxa_info", "DP_R_taxa_info", FALSE, TRUE,
    "path_bioreg", "DP_R_bioreg_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Preparing input data ----
  ecokit::cat_time("Preparing input data", verbose = verbose)

  # # ................................ ###

  ## Current species info
  ecokit::cat_time("Current species info", level = 1L, verbose = verbose)

  species2 <- ecokit::replace_space(species)
  ias_id <- readr::read_tsv(
    file = path_taxa_info, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(species_name2 == species2)
  species_file <- unique(ias_id$species_file)
  ias_id <- unique(ias_id$ias_id)
  taxon_name2 <- ecokit::load_as(path_taxa_info_rdata) %>%
    dplyr::filter(species_name == species) %>%
    dplyr::pull("taxon_name") %>%
    unique()
  gbif_keys <- ecokit::load_as(path_taxa_info_rdata) %>%
    dplyr::filter(species_name == species) %>%
    dplyr::pull("speciesKey")

  # # ................................ ###

  ## Check directories ----
  ecokit::cat_time("Check directories", level = 1L, verbose = verbose)

  ecokit::cat_time(
    "Check and create directories", level = 2L, verbose = verbose)
  path_pa_summary <- fs::path(path_pa, "pa_summary")
  path_pa_tif <- fs::path(path_pa, "pa_tif")
  path_pa_r <- fs::path(path_pa, "pa_raster")
  path_pa_jpeg <- fs::path(path_pa, "distribution_jpeg")

  c(path_pa_summary, path_pa_tif, path_pa_r, path_pa_jpeg) %>%
    purrr::walk(fs::dir_create)

  path_jpeg <- fs::path(path_pa_jpeg, paste0(species_file, ".jpeg"))
  path_summary <- fs::path(
    path_pa_summary, paste0(species_file, "_summary.RData"))

  # # ................................ ###

  ecokit::cat_time(
    "Check path for easin and GBIF data", level = 2L, verbose = verbose)
  path_gbif_data <- fs::path(path_gbif, "species_data")
  path_easin <- fs::path(path_easin, "species_data")

  if (!dir.exists(path_gbif_data)) {
    ecokit::stop_ctx(
      "Required path for GBIF data do not exist",
      path_gbif_data = path_gbif_data, include_backtrace = TRUE)
  }

  if (!dir.exists(path_easin)) {
    ecokit::stop_ctx(
      "Required path for easin data do not exist", path_easin = path_easin,
      include_backtrace = TRUE)
  }

  # # ................................ ###

  ecokit::cat_time("eLTER data", level = 1L, verbose = verbose)
  elter_data <- ecokit::load_as(path_elter) %>%
    dplyr::filter(species_name == species)

  # # ................................ ###

  # Loading reference grids ----
  ecokit::cat_time("Loading reference grids", level = 1L, verbose = verbose)

  grid_paths <- fs::path(
    path_grid,
    c(
      "grid_10_land_sf_country.RData", "grid_10_land_crop.RData",
      "grid_10_land_sf.RData"))

  if (!all(file.exists(grid_paths))) {
    ecokit::stop_ctx(
      paste0(
        "grid files do not exist: \n  >>> ",
        paste(grid_paths[!file.exists(grid_paths)], collapse = "\n  >>> ")),
      grid_paths = grid_paths,
      missing_grid_paths = grid_paths[!file.exists(grid_paths)],
      include_backtrace = TRUE)
  }

  ### sf - 100 km ----
  ecokit::cat_time("sf - 100 km", level = 2L, verbose = verbose)
  grid_100_file <- fs::path(path_grid_ref, "Grid_100_sf.RData")
  if (!file.exists(grid_100_file)) {
    ecokit::stop_ctx(
      "grid file does not exist: Grid_100_sf.RData",
      grid_100_file = grid_100_file, include_backtrace = TRUE)
  }
  grid_100_sf <- ecokit::load_as(grid_100_file) %>%
    magrittr::extract2("Grid_100_sf_s")

  ### sf - 10 km ----
  ecokit::cat_time("sf - 10 km - CNT", level = 2L, verbose = verbose)
  grid_10_country <- fs::path(
    path_grid, "grid_10_land_sf_country.RData") %>%
    ecokit::load_as()

  ### raster - 10 km ----
  ecokit::cat_time("raster - 10 km", level = 2L, verbose = verbose)
  reference_grid <- fs::path(path_grid, "grid_10_land_crop.RData") %>%
    ecokit::load_as(unwrap_r = TRUE)

  ### raster - 100 km ----
  ecokit::cat_time("raster - 100 km", level = 2L, verbose = verbose)
  grid_100_land <- fs::path(path_grid, "grid_10_land_sf.RData") %>%
    ecokit::load_as() %>%
    sf::st_geometry() %>%
    sf::st_centroid() %>%
    sf::st_filter(x = grid_100_sf, join = sf::st_within)

  grid_100_empty <- dplyr::slice(grid_100_sf, 0)

  rm(grid_100_sf, envir = environment())
  invisible(gc())

  # # ................................ ###

  ## `pa_100km` function - species distribution at 100 km ----
  ecokit::cat_time(
    "`pa_100km` function - species distribution at 100 km", level = 1L,
    verbose = verbose)

  pa_100km <- function(species_sf, grid100 = grid_100_land) {

    sf::st_filter(x = grid100, y = species_sf, join = sf::st_within) %>%
      dplyr::select("geometry") %>%
      tibble::tibble() %>%
      dplyr::distinct() %>%
      sf::st_as_sf()
  }

  # # ................................ ###

  ## Exclude countries with only cultivated or casual observations ----
  ecokit::cat_time(
    "Exclude countries with only cultivated or casual observations",
    level = 1L, verbose = verbose)

  valid_countries <- c(
    "Albania", "Andorra", "Austria", "Belgium", "Bosnia and Herzegovina",
    "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia",
    "Faroes", "Finland", "France", "Germany", "Greece", "Guernsey",
    "Hungary", "Iceland", "Ireland", "Isle of Man", "Italy", "Jersey",
    "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Malta",
    "Monaco", "Montenegro", "Netherlands", "North Macedonia", "Norway",
    "Poland", "Portugal", "Romania", "Serbia", "Slovakia", "Slovenia",
    "Spain", "Sweden", "Switzerland", "United Kingdom", "Vatican City") %>%
    sort()

  # List of countries to be excluded for the current species
  countries_to_exclude <- readxl::read_xlsx(
    path = path_taxa_country, sheet = 1) %>%
    dplyr::rename(status_decision = `status-decision`) %>%
    # filter data on the current species and only keep countries with
    # `naturalized` status
    dplyr::filter(
      taxon_name == taxon_name2, status_decision != "naturalized") %>%
    # rename countries to match country names of the countries boundaries
    dplyr::mutate(
      country = dplyr::case_when(
        country == "United_Kingdom" ~ "United Kingdom",
        country == "North_Macedonia" ~ "North Macedonia",
        country == "Isle_of_Man" ~ "Isle of Man",
        country == "Bosnia_and_Herzegovina" ~ "Bosnia and Herzegovina",
        .default = country)) %>%
    # Exclude Turkey (if exist), as it is currently not in the current study
    # area
    dplyr::filter(country != "Turkey") %>%
    dplyr::pull("country")

  if (length(countries_to_exclude) > 0 &&
      !all(countries_to_exclude %in% valid_countries)) {
    is_valid_country <- countries_to_exclude %in% valid_countries
    invalid_countries <- countries_to_exclude[!is_valid_country]
    ecokit::stop_ctx(
      "Some of the to be excluded countries are invalid",
      invalid_countries = invalid_countries,
      n_invalid_countries = length(invalid_countries),
      include_backtrace = TRUE)
  }

  # Mask grid to exclude countries - `TRUE` for grid cells to be considered as
  # presence if present in any of the data source; `FALSE` for grid cells need
  # to be masked as 0 in species distribution maps (1 becomes 0)
  if (length(countries_to_exclude) > 0) {

    ecokit::cat_time(
      paste0(
        "There are ", length(countries_to_exclude), " countries to exclude:"),
      level = 2L, cat_timestamp = FALSE, verbose = verbose)

    ecokit::cat_time(
      paste(sort(countries_to_exclude), collapse = " + "),
      level = 3L, cat_timestamp = FALSE, verbose = verbose)

    mask_keep <- grid_10_country %>%
      dplyr::mutate(keep = !(country %in% countries_to_exclude)) %>%
      dplyr::select("keep") %>%
      terra::rasterize(y = reference_grid, field = "keep") %>%
      terra::as.bool() %>%
      terra::toMemory() %>%
      stats::setNames("mask_keep")
  } else {
    ecokit::cat_time(
      "No countries are excluded", level = 2L,
      cat_timestamp = FALSE, verbose = verbose)
    mask_keep <- terra::as.bool(reference_grid) %>%
      terra::toMemory() %>%
      stats::setNames("mask_keep")
  }

  rm(grid_10_country, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Preparing species data from the 3 sources -----
  ecokit::cat_time(
    "Preparing species data from the 3 sources", verbose = verbose)

  ## 1. GBIF -----
  ecokit::cat_time("1. GBIF", level = 1L, verbose = verbose)

  gbif_keys <- paste(gbif_keys, collapse = "_")

  # Path for the current species GBIF data
  path_gbif_data <- fs::path(path_gbif_data, paste0(species_file, ".RData"))
  path_gbif_r <- fs::path(
    path_gbif, "species_raster", paste0(species_file, "_raster.RData"))

  if (all(file.exists(path_gbif_data, path_gbif_r))) {
    # GBIF data as sf object
    gbif_data <- ecokit::load_as(path_gbif_data)

    # GBIF data as raster map
    gbif_r <- ecokit::load_as(path_gbif_r, unwrap_r = TRUE) %>%
      # convert raster map into binary (1/0)
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames(c("gbif_all", "gbif_cz", "gbif_others"))

    # presence 100 km grid
    gbif_grid100 <- pa_100km(gbif_data)

    rm(gbif_data, envir = environment())
    invisible(gc())

    # Map to export out of this function
    gbif_r_out <- stats::setNames(
      gbif_r, paste0(species_file, "_", c("all", "cz", "others"))) %>%
      terra::wrap() %>%
      list()

    gbif_path <- path_gbif_data
  } else {
    gbif_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("gbif")

    gbif_r_out <- stats::setNames(gbif_r, species_file) %>%
      terra::wrap() %>%
      list()

    gbif_grid100 <- grid_100_empty
    gbif_path <- NA_character_
  }

  # # .................................... ###

  ## 2. easin -----
  ecokit::cat_time("2. easin", level = 1L, verbose = verbose)

  path_easin_data <- fs::path(path_easin, paste0(species_file, "_data.RData"))

  if (file.exists(path_easin_data)) {
    easin_data <- ecokit::load_as(path_easin_data)

    # raster map
    easin_r <- dplyr::select(easin_data, "species_name") %>%
      terra::rasterize(reference_grid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames("easin")

    easin_grid100 <- pa_100km(easin_data)

    rm(easin_data, envir = environment())
    invisible(gc())

    easin_r_out <- stats::setNames(easin_r, species_file) %>%
      terra::wrap() %>%
      list()

    easin_path <- path_easin_data
  } else {
    easin_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("easin")

    easin_r_out <- stats::setNames(easin_r, species_file) %>%
      terra::wrap() %>%
      list()

    easin_grid100 <- grid_100_empty
    easin_path <- NA_character_
  }

  # # .................................... ###

  ## 3. eLTER -----
  ecokit::cat_time("3. eLTER", level = 1L, verbose = verbose)

  if (nrow(elter_data) > 0) {
    elter_r <- dplyr::select(elter_data, "species_name") %>%
      terra::rasterize(reference_grid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames("elter")

    elter_grid100 <- pa_100km(elter_data)

    elter_r_out <- stats::setNames(elter_r, species_file) %>%
      terra::wrap() %>%
      list()
  } else {
    elter_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("elter")

    elter_r_out <- stats::setNames(elter_r, species_file) %>%
      terra::wrap() %>%
      list()

    elter_grid100 <- grid_100_empty
  }

  rm(elter_data, grid_100_land, grid_100_empty, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ## 4. Merging data from the 3 data sources -----
  ecokit::cat_time(
    "Merging data from the 3 data sources", level = 1L, verbose = verbose)

  species_pa <- sum(gbif_r$gbif_others, easin_r, elter_r) %>%
    ecokit::raster_to_pres_abs() %>%
    stats::setNames("pa") %>%
    terra::mask(reference_grid) %>%
    update_citizen(
      ecokit::raster_to_pres_abs(gbif_r$gbif_cz), dist_km = dist_citizen)

  species_pa$pa_masked <- (species_pa$pa * mask_keep)

  species_pa <- c(gbif_r, easin_r, elter_r, species_pa, mask_keep) %>%
    # Ensure that values are read from memory
    terra::toMemory()

  # citizen science grid cells that are excluded from the data
  species_pa$pa_cz_excluded <- all(species_pa$gbif_all == 1, species_pa$pa == 0)

  species_pa$pa_cz_excluded_masked <-
    (species_pa$pa_cz_excluded * species_pa$mask_keep)

  maps_order <- c(
    "gbif_all", "gbif_cz", "gbif_others", "easin", "elter", "pa",
    "pa_cz_excluded", "mask_keep", "pa_masked", "pa_cz_excluded_masked")
  species_pa <- terra::subset(species_pa, maps_order)

  rm(mask_keep, envir = environment())
  invisible(gc())

  # number of cells with values
  n_cells_all <- terra::global(species_pa$pa == 1, sum, na.rm = TRUE) %>%
    as.integer()
  n_cells_naturalized <- terra::global(
    x = species_pa$pa_masked == 1, sum, na.rm = TRUE) %>%
    as.integer()
  n_cells_cz_excl <- terra::global(
    x = species_pa$pa_cz_excluded == 1, sum, na.rm = TRUE) %>%
    as.integer()
  n_cells_cz_excl_masked <- terra::global(
    x = species_pa$pa_cz_excluded_masked == 1, sum, na.rm = TRUE) %>%
    as.integer()

  # # ..................................................................... ###

  # Processing species data ----
  ecokit::cat_time("Processing species data", verbose = verbose)

  # If there is no presence grid cell for the current species, return NA early
  if (n_cells_all == 0 && n_cells_cz_excl == 0) {
    return(tibble::tibble(species = species, pa_summary = NA_character_))
  }

  # # .................................... ###

  ## Save maps ------

  ecokit::cat_time("Save species data", level = 1L, verbose = verbose)

  ### RData -----
  ecokit::cat_time("`.RData`", level = 2L, verbose = verbose)
  ecokit::save_as(
    object = terra::wrap(species_pa),
    object_name = paste0(species_file, "_pa"),
    out_path = fs::path(path_pa_r, paste0(species_file, "_pa.RData")))

  ### tif - all presence grid cells -----
  ecokit::cat_time(
    "`.tif` - all presence grid cells", level = 2L, verbose = verbose)
  terra::writeRaster(
    x = stats::setNames(species_pa$pa, species_file), overwrite = TRUE,
    filename = fs::path(path_pa_tif, paste0(species_file, "_all.tif")))

  ### tif - Excluding cultivated or casual observations -----
  ecokit::cat_time(
    "`.tif` - Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)
  terra::writeRaster(
    x = stats::setNames(species_pa$pa_masked, paste0(species_file, "_masked")),
    overwrite = TRUE,
    filename = fs::path(path_pa_tif, paste0(species_file, "_masked.tif")))

  # # .................................... ###

  ## Biogeographical regions ----

  # number of biogeographical regions and grid cells per region
  ecokit::cat_time(
    "Analysis per biogeographical regions", level = 1L, verbose = verbose)

  ### All -----
  ecokit::cat_time("all data", level = 2L, verbose = verbose)

  bioreg_r <- fs::path(path_bioreg, "bioreg_r.RData")
  if (isFALSE(ecokit::check_data(bioreg_r, warning = FALSE))) {
    ecokit::stop_ctx(
      "Required file for biogeographical regions does not exist",
      bioreg_r = bioreg_r, include_backtrace = TRUE)
  }
  bioreg_r <- ecokit::load_as(bioreg_r, unwrap_r = TRUE)

  # name of Biogeographical regions
  bioreg_names <- c(
    "alpine", "arctic", "atlantic", "blackSea", "boreal",
    "continental", "mediterranean", "pannonian", "steppic") %>%
    paste0("bioreg_", .)

  species_bioregs <- terra::classify(species_pa$pa, cbind(0, NA)) %>%
    terra::mask(mask = ., x = bioreg_r) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(biogreg = "short_name") %>%
    dplyr::count(biogreg) %>%
    tidyr::pivot_wider(names_from = biogreg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("bioreg_", .x)) %>%
    ecokit::add_missing_columns(0L, bioreg_names) %>%
    # Sort columns
    dplyr::select(tidyselect::all_of(bioreg_names)) %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::all_of(bioreg_names), ~tidyr::replace_na(.x, FALSE)))

  bioreg_n <- sum(species_bioregs > 0)

  ### Masked -----
  ecokit::cat_time(
    "Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)

  bioreg_names2 <- stringr::str_replace(
    string = bioreg_names, pattern = "bioreg_", replacement = "bioreg_masked_")

  species_bioregs_masked <- terra::classify(
    species_pa$pa_masked, cbind(0, NA)) %>%
    terra::mask(mask = ., x = bioreg_r) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(biogreg = "short_name") %>%
    dplyr::count(biogreg) %>%
    tidyr::pivot_wider(names_from = biogreg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("bioreg_masked_", .x)) %>%
    ecokit::add_missing_columns(0L, bioreg_names2) %>%
    dplyr::select(tidyselect::all_of(bioreg_names2)) %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::all_of(bioreg_names2), ~tidyr::replace_na(.x, FALSE)))

  rm(bioreg_r, envir = environment())
  invisible(gc())

  if (nrow(species_bioregs_masked) > 0) {
    bioreg_n_mask <- sum(species_bioregs_masked > 0)
  } else {
    bioreg_n_mask <- 0
  }

  # # .................................... ###

  ## Number of presence grid cells per data source -----
  ecokit::cat_time(
    "# presence grid cells per data source", level = 1L, verbose = verbose)

  # presence grid cells per data type
  summary_maps <- (species_pa$pa == 1 & gbif_r$gbif_all == 1) %>%
    ecokit::raster_to_pres_abs() %>%
    terra::mask(reference_grid) %>%
    stats::setNames("gbif") %>%
    c(easin_r, elter_r)
  summary_maps_masked <- summary_maps * species_pa$mask_keep
  invisible(gc())

  r_values <- as.data.frame(summary_maps, na.rm = TRUE) %>%
    tibble::tibble() %>%
    # Keep only grid cells present for any data source
    dplyr::filter_all(dplyr::any_vars(. > 0))

  r_values_masked <- as.data.frame(summary_maps_masked, na.rm = TRUE) %>%
    tibble::tibble() %>%
    dplyr::filter_all(dplyr::any_vars(. > 0))

  n_gbif <- as.integer(terra::global(summary_maps$gbif, sum, na.rm = TRUE))
  n_gbif_unique <- dplyr::filter(
    r_values, gbif == 1, easin == 0, elter == 0) %>%
    dplyr::pull("gbif") %>%
    sum()

  n_gbif_masked <- terra::global(
    summary_maps_masked$gbif, sum, na.rm = TRUE) %>%
    as.integer()
  n_gbif_unique_masked <- dplyr::filter(
    r_values_masked, gbif == 1, easin == 0, elter == 0) %>%
    dplyr::pull("gbif") %>%
    sum()

  n_easin <- as.integer(terra::global(summary_maps$easin, sum, na.rm = TRUE))
  n_easin_unique <- dplyr::filter(
    r_values, gbif == 0, easin == 1, elter == 0) %>%
    dplyr::pull("easin") %>%
    sum()

  n_easin_masked <- terra::global(
    summary_maps_masked$easin, sum, na.rm = TRUE) %>%
    as.integer()
  n_easin_unique_masked <- dplyr::filter(
    r_values_masked, gbif == 0, easin == 1, elter == 0) %>%
    dplyr::pull("easin") %>%
    sum()

  n_elter <- as.integer(terra::global(summary_maps$elter, sum, na.rm = TRUE))
  n_elter_unique <- dplyr::filter(
    r_values, gbif == 0, easin == 0, elter == 1) %>%
    dplyr::pull("elter") %>%
    sum()

  n_elter_masked <- terra::global(
    summary_maps_masked$elter, sum, na.rm = TRUE) %>%
    as.integer()
  n_elter_unique_masked <- dplyr::filter(
    r_values_masked, gbif == 0, easin == 0, elter == 1) %>%
    dplyr::pull("elter") %>%
    sum()

  rm(
    r_values, r_values_masked, gbif_r, easin_r, elter_r,
    summary_maps, summary_maps_masked, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## Number of presence grid cells per country -----
  ecokit::cat_time(
    "# presence grid cells per country", level = 1L, verbose = verbose)

  country_list <- paste0("cnt_", ecokit::replace_space(valid_countries))

  grid_country <- fs::path(path_grid, "grid_10_land_sf_country.RData") %>%
    ecokit::load_as() %>%
    dplyr::select("country")

  species_country <- terra::as.points(species_pa$pa) %>%
    sf::st_as_sf() %>%
    dplyr::filter(pa == 1) %>%
    dplyr::select(-"pa") %>%
    sf::st_join(grid_country) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble() %>%
    dplyr::mutate(country = ecokit::replace_space(country)) %>%
    dplyr::count(country) %>%
    tidyr::pivot_wider(names_from = country, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("cnt_", .x)) %>%
    ecokit::add_missing_columns(0L, country_list) %>%
    dplyr::select(tidyselect::all_of(country_list))

  species_country_masked <- terra::as.points(species_pa$pa_masked) %>%
    sf::st_as_sf() %>%
    dplyr::filter(pa_masked == 1) %>%
    dplyr::select(-"pa_masked") %>%
    sf::st_join(grid_country) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble() %>%
    dplyr::mutate(country = ecokit::replace_space(country)) %>%
    dplyr::count(country) %>%
    tidyr::pivot_wider(names_from = country, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("cnt_", .x)) %>%
    ecokit::add_missing_columns(0L, country_list) %>%
    dplyr::select(tidyselect::all_of(country_list)) %>%
    dplyr::rename_with(~paste0(.x, "_masked"))

  rm(grid_country, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Prepare and export species summary info -------
  ecokit::cat_time("Prepare and export species summary info", verbose = verbose)

  binary_r_out <- stats::setNames(species_pa$pa, species_file) %>%
    terra::wrap() %>%
    list()

  binary_r_masked_out <- stats::setNames(species_pa$pa_masked, species_file) %>%
    terra::wrap() %>%
    list()

  binary_r_mask_keep <- stats::setNames(species_pa$mask_keep, species_file) %>%
    terra::wrap() %>%
    list()

  integer_columns <- c(
    "n_cells_all", "n_cells_naturalized",
    "gbif", "gbif_unique", "gbif_masked", "gbif_masked_unique",
    "easin", "easin_unique", "easin_masked", "easin_masked_unique",
    "elter", "elter_unique", "elter_masked", "elter_masked_unique",
    "bioreg_n", "bioreg_n_mask")

  results <- tibble::tibble(

    # ||||||||||||||||||||||||||||||||||
    # Species information
    # ||||||||||||||||||||||||||||||||||
    ias_id = ias_id,
    species = species,
    taxon_name = taxon_name2,
    species_file = species_file,
    gbif_keys = gbif_keys,

    # ||||||||||||||||||||||||||||||||||
    # Number of presence grid cells
    # ||||||||||||||||||||||||||||||||||
    # merged PA data (excl. isolated citizen science data)
    n_cells_all = n_cells_all,
    # merged PA data (excl. countries with cultivated and casual observations)
    n_cells_naturalized = n_cells_naturalized,
    # number of excluded citizen science grid cells (full Europe)
    n_cells_cz_excl = n_cells_cz_excl,
    # number of excluded citizen science grid cells (only countries with
    # naturalized status)
    n_cells_cz_excl_masked = n_cells_cz_excl_masked,

    # ||||||||||||||||||||||||||||||||||
    # GBIF
    # ||||||||||||||||||||||||||||||||||
    # number of grid cells in the final PA map from GBIF
    gbif = n_gbif,
    # number of unique grid cells from GBIF in the final PA map
    gbif_unique = n_gbif_unique,
    # number of grid cells in the final PA map from GBIF (masked)
    gbif_masked = n_gbif_masked,
    # number of unique grid cells from GBIF in the final PA map (masked)
    gbif_masked_unique = n_gbif_unique_masked,
    # path for GBIF input data
    gbif_path = gbif_path,
    # GBIF data at 100 km grid; for plotting
    gbif_grid100 = list(gbif_grid100),
    # SpatRasters for species' GBIF data
    gbif_r = gbif_r_out,

    # ||||||||||||||||||||||||||||||||||
    # EASIN
    # ||||||||||||||||||||||||||||||||||
    # number of grid cells in the final PA map from EASIN
    easin = n_easin,
    # number of unique grid cells from EASIN in the final PA map
    easin_unique = n_easin_unique,
    # number of grid cells in final PA map from EASIN (masked)
    easin_masked = n_easin_masked,
    # number of unique grid cells from EASIN in final PA map (masked)
    easin_masked_unique = n_easin_unique_masked,
    # path for EASIN input data
    easin_path = easin_path,
    # EASIN data at 100 km grid; for plotting
    easin_grid100 = list(easin_grid100),
    # SpatRasters for species' EASIN data
    easin_r = easin_r_out,

    # ||||||||||||||||||||||||||||||||||
    # eLTER
    # ||||||||||||||||||||||||||||||||||
    # number of grid cells in the final PA map from EASIN
    elter = n_elter,
    # number of unique grid cells from EASIN in the final PA map
    elter_unique = n_elter_unique,
    # number of grid cells in final PA map from EASIN (masked)
    elter_masked = n_elter_masked,
    # number of unique grid cells from EASIN in final PA map (masked)
    elter_masked_unique = n_elter_unique_masked,
    # EASIN data at 100 km grid; for plotting
    elter_grid100 = list(elter_grid100),
    # SpatRasters for species' eLTER data
    elter_r = elter_r_out,

    # ||||||||||||||||||||||||||||||||||
    # SpatRasters for distribution maps
    # ||||||||||||||||||||||||||||||||||
    # Final PA map for the species
    pa_map = binary_r_out,
    # Final PA map for the species (masked)
    pa_masked_map = binary_r_masked_out,
    # mask for included/excluded countries
    mask_keep = binary_r_mask_keep,
    # Excluded citizen science data
    pa_cz_excluded = list(terra::wrap(species_pa$pa_cz_excluded)),
    pa_cz_excluded_masked = list(terra::wrap(species_pa$pa_cz_excluded_masked)),

    # ||||||||||||||||||||||||||||||||||
    # Country summary
    # ||||||||||||||||||||||||||||||||||
    # summary for the number of grid cell per country
    species_country = list(species_country),
    # summary for the number of grid cell per country (masked)
    species_country_masked = list(species_country_masked),
    # list of excluded countries
    countries_to_exclude = list(countries_to_exclude),

    # ||||||||||||||||||||||||||||||||||
    # Biogeographical regions
    # ||||||||||||||||||||||||||||||||||
    # number of presence grid cells in each of biogeographical region
    bioreg_data = list(species_bioregs),
    # unique number of biogeographical regions in presence grid cells
    bioreg_n = bioreg_n,
    # number of presence grid cells in each of biogeographical region (masked)
    bioreg_data_mask = list(species_bioregs_masked),
    # unique number of biogeographical regions in presence grid cells (masked)
    bioreg_n_mask = bioreg_n_mask,

    # Path for species distribution map
    path_jpeg = path_jpeg) %>%
    dplyr::mutate_at(integer_columns, ~ as.integer(.))

  rm(species_pa, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, verbose = verbose,
    prefix = "\nProcessing species data was finished in ", ... = "\n")

  # save species data
  ecokit::save_as(
    object = results, object_name = paste0(species_file, "_summary"),
    out_path = path_summary)

  return(tibble::tibble(species = species, pa_summary = path_summary))

}

# # |------------------------------------------------------------------------| #
# update_citizen ----
## |------------------------------------------------------------------------| #

#' Update species distribution raster with citizen science data within a
#' distance threshold
#'
#' Adds presence grid cells from a citizen science raster (`r_cz`) to an
#' existing species distribution raster (`r_sp`) only if the citizen science
#' cell is within `dist_km` kilometers of any existing cell in `r_sp`, but not
#' within 500 meters of any existing cell in `r_sp`. This avoids including
#' isolated citizen science presences that are too far from known presences.
#'
#' @param r_sp SpatRaster. Binary raster (0, 1, NA) representing merged/primary
#'   species distribution.
#' @param r_cz SpatRaster. Binary raster (0, 1, NA) representing citizen science
#'   presences.
#' @param dist_km numeric. Maximum distance (in kilometers) to allow new citizen
#'   science presences to be added to `r_sp`.
#' @return SpatRaster. Updated raster with citizen science presences added under
#'   the specified constraints.
#' @keywords internal
#' @noRd

update_citizen <- function(r_sp, r_cz, dist_km = 100) {

  if (!inherits(r_sp, "SpatRaster") || !inherits(r_cz, "SpatRaster")) {
    ecokit::stop_ctx(
      "r_sp and r_cz must be SpatRaster objects", include_backtrace = TRUE,
      class_r_sp = class(r_sp), class_r_cz = class(r_cz))
  }
  if (terra::nlyr(r_sp) != 1 || terra::nlyr(r_cz) != 1) {
    ecokit::stop_ctx(
      "r_sp and r_cz must be single layer SpatRaster objects",
      include_backtrace = TRUE,
      nlyr_r_sp = terra::nlyr(r_sp), nlyr_r_cz = terra::nlyr(r_cz))
  }
  if (!all(terra::ext(r_sp) == terra::ext(r_cz))) {
    ecokit::stop_ctx(
      "r_sp and r_cz must have the same extent",
      include_backtrace = TRUE,
      ext_r_sp = as.character(terra::ext(r_sp)),
      ext_r_cz = as.character(terra::ext(r_cz)))
  }
  if (!all(terra::res(r_sp) == terra::res(r_cz))) {
    ecokit::stop_ctx(
      "r_sp and r_cz must have the same resolution",
      include_backtrace = TRUE)
  }
  if (!all(unique(terra::values(r_sp)) %in% c(NaN, NA, 0, 1))) {
    ecokit::stop_ctx(
      "r_sp must be a binary raster with values of 0, 1, or NA",
      include_backtrace = TRUE)
  }
  if (!all(unique(terra::values(r_cz)) %in% c(NaN, NA, 0, 1))) {
    ecokit::stop_ctx(
      "r_cz must be a binary raster with values of 0, 1, or NA",
      include_backtrace = TRUE)
  }

  # Convert r_sp and r_cz to sf points
  r_sp_sf <- terra::classify(r_sp == 1, cbind(0, NA)) %>%
    stats::setNames("r_sp") %>%
    as.data.frame(xy = TRUE) %>%
    dplyr::filter(r_sp == 1)
  if (nrow(r_sp_sf) == 0) {
    return(r_sp)
  }
  r_sp_sf <- sf::st_as_sf(r_sp_sf, coords = c("x", "y"), crs = 3035L)
  r_cz_sf <- terra::classify(r_cz == 1, cbind(0, NA)) %>%
    stats::setNames("r_cz") %>%
    as.data.frame(xy = TRUE) %>%
    dplyr::filter(r_cz == 1)
  if (nrow(r_cz_sf) == 0) {
    return(r_sp)
  }

  r_cz_sf <- sf::st_as_sf(r_cz_sf, coords = c("x", "y"), crs = 3035L)

  # Identify points in r_cz_sf within dist_km of any point in r_sp_sf, but not
  # within 500 m of any point in r_sp_sf
  points_to_add <- sf::st_filter(
    x = r_cz_sf, y = r_sp_sf,
    .predicate = sf::st_is_within_distance, dist = dist_km * 1000L)

  if (nrow(points_to_add) == 0) {
    return(r_sp)
  }

  # exclude overlapped points with r_sp_sf
  points_to_add <- dplyr::filter(
    points_to_add,
    rowSums(
      sf::st_is_within_distance(
        points_to_add, r_sp_sf, dist = 500, sparse = FALSE)) == 0)

  if (nrow(points_to_add) == 0) {
    return(r_sp)
  }

  # Rasterize the points to add back to r_sp
  terra::rasterize(points_to_add, r_sp, field = 1, update = TRUE) %>%
    stats::setNames(names(r_sp))
}
