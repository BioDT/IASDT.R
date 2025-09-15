# # |------------------------------------------------------------------------| #
# IAS_distribution ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name IAS_data
#' @rdname IAS_data
#' @order 2

IAS_distribution <- function(
    species = NULL, env_file = ".env", verbose = FALSE) {

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
  path_grid <- path_taxa_info <- path_GBIF <- path_eLTER <- path_EASIN <-
    path_bioreg <- path_grid_ref <- path_PA <- biogreg <-
    country <- species_name2 <- species_name <- n <- GBIF <- EASIN <- eLTER <-
    PA <- `status-decision` <- path_taxa_country <- country <-
    gbif_key <- status_decision <- path_taxa_info_rdata <- NULL

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
    "path_GBIF", "DP_R_gbif_processed", TRUE, FALSE,
    "path_EASIN", "DP_R_easin_processed", TRUE, FALSE,
    "path_eLTER", "DP_R_elter_processed", FALSE, TRUE,
    "path_PA", "DP_R_PA", FALSE, FALSE,
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

  Species2 <- ecokit::replace_space(species)
  ias_id <- readr::read_tsv(
    file = path_taxa_info, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(species_name2 == Species2)
  species_file <- unique(ias_id$species_file)
  ias_id <- dplyr::pull(ias_id, "ias_id") %>%
    unique() %>%
    stringr::str_pad(width = 4, pad = "0")
  GBIF_keys <- ecokit::load_as(path_taxa_info_rdata) %>%
    dplyr::filter(species_name == species) %>%
    dplyr::pull("speciesKey")

  # # ................................ ###

  ## Check directories ----
  ecokit::cat_time("Check directories", level = 1L, verbose = verbose)

  ecokit::cat_time(
    "Check and create directories", level = 2L, verbose = verbose)
  path_PA_summary <- fs::path(path_PA, "PA_summary")
  path_PA_tif <- fs::path(path_PA, "PA_tif")
  path_PA_r <- fs::path(path_PA, "PA_r")
  path_PA_JPEG <- fs::path(path_PA, "distribution_jpeg")

  c(path_PA_summary, path_PA_tif, path_PA_r, path_PA_JPEG) %>%
    purrr::walk(fs::dir_create)

  out_JPEG <- fs::path(path_PA_JPEG, paste0(species_file, ".jpeg"))
  out_qs2 <- fs::path(path_PA_summary, paste0(species_file, ".qs2"))

  # # ................................ ###

  ecokit::cat_time(
    "Check path for EASIN and GBIF data", level = 2L, verbose = verbose)
  path_GBIF_data <- fs::path(path_GBIF, "species_data")
  path_EASIN <- fs::path(path_EASIN, "species_data")

  if (!dir.exists(path_GBIF_data)) {
    ecokit::stop_ctx(
      "Required path for GBIF data do not exist",
      path_GBIF_data = path_GBIF_data, include_backtrace = TRUE)
  }

  if (!dir.exists(path_EASIN)) {
    ecokit::stop_ctx(
      "Required path for EASIN data do not exist", path_EASIN = path_EASIN,
      include_backtrace = TRUE)
  }

  # # ................................ ###

  ecokit::cat_time("eLTER data", level = 1L, verbose = verbose)
  eLTER_data <- ecokit::load_as(path_eLTER) %>%
    dplyr::filter(species_name == species)

  # # ................................ ###

  # Loading reference grids ----
  ecokit::cat_time("Loading reference grids", level = 1L, verbose = verbose)
  grid_100_file <- fs::path(path_grid_ref, "Grid_100_sf.RData")
  if (!file.exists(grid_100_file)) {
    ecokit::stop_ctx(
      "grid file does not exist: Grid_100_sf.RData",
      grid_100_file = grid_100_file, include_backtrace = TRUE)
  }

  grid_paths <- fs::path(
    path_grid,
    c(
      "grid_10_land_crop_sf_country.RData", "grid_10_land_crop.RData",
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
  grid_100_sf <- ecokit::load_as(grid_100_file) %>%
    magrittr::extract2("Grid_100_sf_s")

  ### sf - 10 km ----
  ecokit::cat_time("sf - 10 km - CNT", level = 2L, verbose = verbose)
  grid_10_country <- fs::path(
    path_grid, "grid_10_land_crop_sf_country.RData") %>%
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

  ## `PA_100km` function - species distribution at 100 km ----
  ecokit::cat_time(
    "`PA_100km` function - species distribution at 100 km", level = 1L,
    verbose = verbose)

  PA_100km <- function(species_sf, grid100 = grid_100_land) {
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

  # List of countries to be excluded for the current species
  countries_to_exclude <- readxl::read_xlsx(
    path = path_taxa_country, sheet = 1) %>%
    dplyr::rename(status_decision = `status-decision`) %>%
    # filter data on the current species and only keep countries with `cult+cas`
    # or `delete` status
    dplyr::filter(
      gbif_key %in% GBIF_keys,
      status_decision %in% c("cult+cas", "delete")) %>%
    # rename countries to match country names of the countries boundaries
    dplyr::mutate(
      country = dplyr::case_when(
        country == "United_Kingdom" ~ "United Kingdom",
        country == "North_Macedonia" ~ "North Macedonia",
        country == "Isle_of_Man" ~ "Isle of Man",
        country == "Bosnia_and_Herzegovina" ~ "Bosnia and Herzegovina",
        .default = country)) %>%
    # Exclude Turkey, as it is currently not in the current study area
    dplyr::filter(country != "Turkey") %>%
    dplyr::pull("country")

  ecokit::cat_time(
    paste0(
      "There are ", length(countries_to_exclude), " countries to exclude:"),
    level = 2L, cat_timestamp = FALSE, verbose = verbose)

  # Mask grid to exclude countries - `TRUE` for grid cells to be considered as
  # presence if present in any of the data source; `FALSE` for grid cells need
  # to be masked as 0 in species distribution maps (1 becomes 0)
  if (length(countries_to_exclude) > 0) {

    ecokit::cat_time(
      paste(sort(countries_to_exclude), collapse = " + "),
      level = 3L, cat_timestamp = FALSE, verbose = verbose)

    mask_keep <- grid_10_country %>%
      dplyr::mutate(Keep = !(country %in% countries_to_exclude)) %>%
      dplyr::select("Keep") %>%
      terra::rasterize(y = reference_grid, field = "Keep") %>%
      terra::as.bool() %>%
      terra::toMemory() %>%
      stats::setNames("mask_keep")
  } else {
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

  GBIF_keys <- paste(GBIF_keys, collapse = "_")

  # Path for the current species GBIF data
  path_GBIF_data <- fs::path(path_GBIF_data, paste0(species_file, ".RData"))
  path_GBIF_r <- fs::path(
    path_GBIF, "species_raster", paste0(species_file, "_raster.RData"))

  if (all(file.exists(path_GBIF_data, path_GBIF_r))) {
    # GBIF data as sf object
    GBIF_data <- ecokit::load_as(path_GBIF_data)

    # GBIF data as raster map
    GBIF_r <- ecokit::load_as(path_GBIF_r) %>%
      terra::unwrap() %>%
      # convert raster map into binary (1/0)
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames("GBIF")

    # presence 100 km grid
    GBIF_grid100 <- PA_100km(GBIF_data)

    rm(GBIF_data, envir = environment())
    invisible(gc())

    # Map to export out of this function
    GBIF_r_out <- stats::setNames(GBIF_r, species_file) %>%
      terra::wrap() %>%
      list()

    GBIF_path <- path_GBIF_data
  } else {
    GBIF_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("GBIF")

    GBIF_r_out <- stats::setNames(GBIF_r, species_file) %>%
      terra::wrap() %>%
      list()

    GBIF_grid100 <- grid_100_empty
    GBIF_path <- NA_character_
  }

  # # .................................... ###

  ## 2. EASIN -----
  ecokit::cat_time("2. EASIN", level = 1L, verbose = verbose)

  path_EASIN_data <- fs::path(path_EASIN, paste0(species_file, "_data.RData"))

  if (file.exists(path_EASIN_data)) {
    EASIN_data <- ecokit::load_as(path_EASIN_data)

    # raster map
    EASIN_r <- dplyr::select(EASIN_data, "species_name") %>%
      terra::rasterize(reference_grid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames("EASIN")

    EASIN_grid100 <- PA_100km(EASIN_data)

    rm(EASIN_data, envir = environment())
    invisible(gc())

    EASIN_r_out <- stats::setNames(EASIN_r, species_file) %>%
      terra::wrap() %>%
      list()

    EASIN_path <- path_EASIN_data
  } else {
    EASIN_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("EASIN")

    EASIN_r_out <- stats::setNames(EASIN_r, species_file) %>%
      terra::wrap() %>%
      list()

    EASIN_grid100 <- grid_100_empty
    EASIN_path <- NA_character_
  }

  # # .................................... ###

  ## 3. eLTER -----
  ecokit::cat_time("3. eLTER", level = 1L, verbose = verbose)

  if (nrow(eLTER_data) > 0) {
    eLTER_r <- dplyr::select(eLTER_data, "species_name") %>%
      terra::rasterize(reference_grid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(reference_grid) %>%
      stats::setNames("eLTER")

    eLTER_grid100 <- PA_100km(eLTER_data)

    eLTER_r_out <- stats::setNames(eLTER_r, species_file) %>%
      terra::wrap() %>%
      list()
  } else {
    eLTER_r <- terra::classify(reference_grid, cbind(1, 0)) %>%
      stats::setNames("eLTER")

    eLTER_r_out <- stats::setNames(eLTER_r, species_file) %>%
      terra::wrap() %>%
      list()

    eLTER_grid100 <- grid_100_empty
  }

  rm(eLTER_data, grid_100_land, grid_100_empty, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ## 4. Merging data from the 3 data sources -----
  ecokit::cat_time(
    "Merging data from the 3 data sources", level = 1L, verbose = verbose)

  species_PA <- sum(GBIF_r, EASIN_r, eLTER_r) %>%
    ecokit::raster_to_pres_abs() %>%
    stats::setNames("PA") %>%
    terra::mask(reference_grid)
  species_PA$PA_masked <- (species_PA$PA * mask_keep)
  species_PA <- c(GBIF_r, EASIN_r, eLTER_r, species_PA, mask_keep) %>%
    # Ensure that values are read from memory
    terra::toMemory()

  rm(reference_grid, mask_keep, envir = environment())
  invisible(gc())

  # number of cells with values
  PA_n_cells_all <- terra::global(species_PA$PA == 1, sum, na.rm = TRUE) %>%
    as.integer()
  PA_n_cells_naturalized <- terra::global(
    x = species_PA$PA_masked == 1, sum, na.rm = TRUE) %>%
    as.integer()

  # # ..................................................................... ###

  # Processing species data ----
  ecokit::cat_time("Processing species data", verbose = verbose)

  # If there is no presence grid cell for the current species, return NA early
  if (PA_n_cells_all == 0) {
    return(tibble::tibble(species = species, PA_summary = NA_character_))
  }

  # # .................................... ###

  ## Save maps ------

  ecokit::cat_time("Save species data", level = 1L, verbose = verbose)

  ### RData -----
  ecokit::cat_time("`.RData`", level = 2L, verbose = verbose)
  path_rData <- fs::path(path_PA_r, paste0(species_file, "_PA.RData"))
  ecokit::save_as(
    object = terra::wrap(species_PA), object_name = paste0(species_file, "_PA"),
    out_path = path_rData)

  ### tif - all presence grid cells -----
  ecokit::cat_time(
    "`.tif` - all presence grid cells", level = 2L, verbose = verbose)
  terra::writeRaster(
    x = species_PA$PA, overwrite = TRUE,
    filename = fs::path(path_PA_tif, paste0(species_file, "_all.tif")))

  ### tif - Excluding cultivated or casual observations -----
  ecokit::cat_time(
    "`.tif` - Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)
  terra::writeRaster(
    x = species_PA$PA_masked, overwrite = TRUE,
    filename = fs::path(path_PA_tif, paste0(species_file, "_masked.tif")))

  ### NetCDF -----
  # NOT IMPLEMENTED YET

  # # .................................... ###

  ## Biogeographical regions ----
  ecokit::cat_time(
    "Analysis per biogeographical regions", level = 1L, verbose = verbose)

  ### All -----
  ecokit::cat_time("all data", level = 2L, verbose = verbose)

  # Number of grid per each biogeographical region
  #
  # number of biogeographical regions per species and minimum / maximum / mean
  # number of grid cells per biogeographical regions

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

  species_bioregs <- terra::classify(species_PA$PA, cbind(0, NA)) %>%
    terra::mask(mask = ., x = bioreg_r) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(biogreg = "short_name") %>%
    dplyr::count(biogreg) %>%
    tidyr::pivot_wider(names_from = biogreg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("bioreg_", .x)) %>%
    ecokit::add_missing_columns(0L, bioreg_names) %>%
    # Sort columns
    dplyr::select(tidyselect::all_of(bioreg_names))

  bioreg_summ_min <- min(species_bioregs, na.rm = TRUE)
  bioreg_summ_max <- max(species_bioregs, na.rm = TRUE)
  bioreg_summ_mean <- unlist(species_bioregs) %>%
    mean(na.rm = TRUE) %>%
    round(1)
  bioreg_summ_n <- sum(species_bioregs > 0)

  ### Masked -----
  ecokit::cat_time(
    "Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)

  bioreg_names2 <- stringr::str_replace(
    string = bioreg_names, pattern = "bioreg_", replacement = "bioreg_masked_")

  sp_bioreg_masked <- terra::classify(species_PA$PA_masked, cbind(0, NA)) %>%
    terra::mask(mask = ., x = bioreg_r) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(biogreg = "short_name") %>%
    dplyr::count(biogreg) %>%
    tidyr::pivot_wider(names_from = biogreg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("bioreg_masked_", .x)) %>%
    ecokit::add_missing_columns(0L, bioreg_names2) %>%
    dplyr::select(tidyselect::all_of(bioreg_names2))

  rm(bioreg_r, envir = environment())
  invisible(gc())

  if (nrow(sp_bioreg_masked) > 0) {
    bioreg_mask_summ_min <- min(sp_bioreg_masked, na.rm = TRUE)
    bioreg_mask_summ_max <- max(sp_bioreg_masked, na.rm = TRUE)
    bioreg_mask_summ_mean <- unlist(sp_bioreg_masked) %>%
      mean(na.rm = TRUE) %>%
      round(1)
    bioreg_mask_summ_n <- sum(sp_bioreg_masked > 0)
  } else {
    bioreg_mask_summ_min <- 0
    bioreg_mask_summ_max <- 0
    bioreg_mask_summ_mean <- 0
    bioreg_mask_summ_n <- 0
  }

  # # .................................... ###

  ## Number of presence grid cells per data provider -----
  ecokit::cat_time(
    "# presence grid cells per data provider", level = 1L, verbose = verbose)

  # presence grid cells per data type
  r_values <- c(GBIF_r, EASIN_r, eLTER_r) %>%
    as.data.frame(na.rm = TRUE) %>%
    tibble::tibble() %>%
    # Keep only grid cells present for any data provider
    dplyr::filter_all(dplyr::any_vars(. > 0))

  r_values_masked <- c(GBIF_r, EASIN_r, eLTER_r) %>%
    magrittr::multiply_by(species_PA$PA_masked) %>%
    as.data.frame(na.rm = TRUE) %>%
    tibble::tibble() %>%
    dplyr::filter_all(dplyr::any_vars(. > 0))

  n_GBIF <- as.integer(terra::global(GBIF_r, sum, na.rm = TRUE))
  n_GBIF_unique <- dplyr::filter(
    r_values, GBIF == 1, EASIN == 0, eLTER == 0) %>%
    dplyr::pull("GBIF") %>%
    sum()

  n_GBIF_masked <- terra::global(
    (species_PA$GBIF * species_PA$mask_keep), sum, na.rm = TRUE) %>%
    as.integer()
  n_GBIF_unique_masked <- dplyr::filter(
    r_values_masked, GBIF == 1, EASIN == 0, eLTER == 0) %>%
    dplyr::pull("GBIF") %>%
    sum()

  n_EASIN <- as.integer(terra::global(EASIN_r, sum, na.rm = TRUE))
  n_EASIN_unique <- dplyr::filter(
    r_values, GBIF == 0, EASIN == 1, eLTER == 0) %>%
    dplyr::pull("EASIN") %>%
    sum()

  n_EASIN_masked <- terra::global(
    (species_PA$EASIN * species_PA$mask_keep), sum, na.rm = TRUE) %>%
    as.integer()
  n_EASIN_unique_masked <- dplyr::filter(
    r_values_masked, GBIF == 0, EASIN == 1, eLTER == 0) %>%
    dplyr::pull("EASIN") %>%
    sum()

  n_eLTER <- as.integer(terra::global(eLTER_r, sum, na.rm = TRUE))
  n_eLTER_unique <- dplyr::filter(
    r_values, GBIF == 0, EASIN == 0, eLTER == 1) %>%
    dplyr::pull("eLTER") %>%
    sum()

  n_eLTER_masked <- terra::global(
    (species_PA$eLTER * species_PA$mask_keep), sum, na.rm = TRUE) %>%
    as.integer()
  n_eLTER_unique_masked <- dplyr::filter(
    r_values_masked, GBIF == 0, EASIN == 0, eLTER == 1) %>%
    dplyr::pull("eLTER") %>%
    sum()

  rm(r_values, r_values_masked, GBIF_r, EASIN_r, eLTER_r, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## Number of presence grid cells per country -----
  ecokit::cat_time(
    "# presence grid cells per country", level = 1L, verbose = verbose)

  country_list <- c(
    "Albania", "Austria", "Belgium", "Bosnia and Herzegovina",
    "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia",
    "Finland", "France", "Germany", "Greece", "Hungary", "Iceland",
    "Ireland", "Isle of Man", "Italy", "Latvia", "Liechtenstein",
    "Lithuania", "Luxembourg", "Monaco", "Montenegro", "Netherlands",
    "North Macedonia", "Norway", "Poland", "Portugal", "Romania",
    "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland",
    "United Kingdom", "Andorra", "Faroes", "Guernsey", "Jersey",
    "Malta", "Vatican City") %>%
    sort() %>%
    ecokit::replace_space() %>%
    paste0("CNT_", .)

  grid_country <- fs::path(path_grid, "grid_10_land_crop_sf_country.RData") %>%
    ecokit::load_as() %>%
    dplyr::select("country")

  species_country <- terra::as.points(species_PA$PA) %>%
    sf::st_as_sf() %>%
    dplyr::filter(PA == 1) %>%
    dplyr::select(-"PA") %>%
    sf::st_join(grid_country) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble() %>%
    dplyr::mutate(country = ecokit::replace_space(country)) %>%
    dplyr::count(country) %>%
    tidyr::pivot_wider(names_from = country, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("CNT_", .x)) %>%
    ecokit::add_missing_columns(0L, country_list) %>%
    dplyr::select(tidyselect::all_of(country_list))

  rm(grid_country, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## Number of unique iNaturalist grid cells -----
  ecokit::cat_time(
    "# unique iNaturalist grid cells", level = 1L, verbose = verbose)

  iNatur_data <- fs::path(path_GBIF, "iNaturalist_count.RData") %>%
    ecokit::load_as() %>%
    dplyr::filter(species == !!species)

  if (nrow(iNatur_data) > 0) {
    iNatur_unique <- iNatur_data$iNaturalist_unique
    iNatur_percent <- round(100 * (iNatur_data$iNaturalist_unique / n_GBIF), 1)
  } else {
    iNatur_unique <- iNatur_percent <- 0
  }

  rm(iNatur_data, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Prepare and export species summary info -------
  ecokit::cat_time("Prepare and export species summary info", verbose = verbose)

  binary_r_out <- stats::setNames(species_PA$PA, species_file) %>%
    terra::wrap() %>%
    list()

  binary_r_masked_out <- stats::setNames(species_PA$PA_masked, species_file) %>%
    terra::wrap() %>%
    list()

  binary_r_mask_keep <- stats::setNames(species_PA$mask_keep, species_file) %>%
    terra::wrap() %>%
    list()

  rm(species_PA, envir = environment())
  invisible(gc())

  integer_columns <- c(
    "n_cells_all", "n_cells_naturalized",
    "GBIF", "GBIF_unique", "GBIF_masked", "GBIF_masked_unique",
    "EASIN", "EASIN_unique", "EASIN_masked", "EASIN_masked_unique",
    "eLTER", "eLTER_unique", "eLTER_masked", "eLTER_masked_unique",
    "bioreg_summ_min", "bioreg_summ_max", "bioreg_summ_n",
    "bioreg_mask_summ_min", "bioreg_mask_summ_max",
    "bioreg_mask_summ_n", "iNaturalist_unique")

  Results <- tibble::tibble(
    species = species,
    GBIF_keys = GBIF_keys,
    species_id = ias_id,
    n_cells_all = PA_n_cells_all,
    n_cells_naturalized = PA_n_cells_naturalized,

    GBIF = n_GBIF,
    GBIF_unique = n_GBIF_unique,
    GBIF_masked = n_GBIF_masked,
    GBIF_masked_unique = n_GBIF_unique_masked,
    GBIF_path = GBIF_path,
    GBIF_grid100 = list(GBIF_grid100),
    GBIF_r = GBIF_r_out,

    EASIN = n_EASIN,
    EASIN_unique = n_EASIN_unique,
    EASIN_masked = n_EASIN_masked,
    EASIN_masked_unique = n_EASIN_unique_masked,
    EASIN_path = EASIN_path,
    EASIN_grid100 = list(EASIN_grid100),
    EASIN_r = EASIN_r_out,

    eLTER = n_eLTER,
    eLTER_unique = n_eLTER_unique,
    eLTER_masked = n_eLTER_masked,
    eLTER_masked_unique = n_eLTER_unique_masked,
    eLTER_grid100 = list(eLTER_grid100),
    eLTER_r = eLTER_r_out,

    PA_map = binary_r_out,
    PA_masked_map = binary_r_masked_out,
    mask_keep = binary_r_mask_keep,

    species_country = list(species_country),
    bioreg_data = list(species_bioregs),
    bioreg_summ_min = bioreg_summ_min,
    bioreg_summ_max = bioreg_summ_max,
    bioreg_summ_mean = bioreg_summ_mean,
    bioreg_summ_n = bioreg_summ_n,
    bioreg_mask_summ_data = list(species_bioregs_masked),
    bioreg_mask_summ_min = bioreg_mask_summ_min,
    bioreg_mask_summ_max = bioreg_mask_summ_max,
    bioreg_mask_summ_mean = bioreg_mask_summ_mean,
    bioreg_mask_summ_n = bioreg_mask_summ_n,
    iNaturalist_unique = iNatur_unique,
    iNatur_percent = iNatur_percent,
    countries_to_exclude = list(countries_to_exclude),
    path_JPEG = out_JPEG) %>%
    dplyr::mutate_at(integer_columns, ~ as.integer(.))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing species data was finished in ", ... = "\n",
    verbose = verbose)

  # save species data
  ecokit::save_as(object = Results, out_path = out_qs2)

  return(tibble::tibble(species = species, PA_summary = out_qs2))

}

# # |------------------------------------------------------------------------| #
# update_citizen ----
## |------------------------------------------------------------------------| #

#' Update species distribution raster with citizen science data within a
#' distance threshold
#'
#' Adds presence grid cells from a citizen science raster (`r1`) to an existing
#' species distribution raster (`r2`) only if the citizen science cell is within
#' `dist_km` kilometers of any existing cell in `r2`, but not within 500 meters
#' of any existing cell in `r2`. This avoids including isolated citizen science
#' presences that are too far from known presences.
#'
#' @param r1 SpatRaster. Binary raster (0, 1, NA) representing citizen science
#'   presences.
#' @param r2 SpatRaster. Binary raster (0, 1, NA) representing merged/primary
#'   species distribution.
#' @param dist_km numeric. Maximum distance (in kilometers) to allow new citizen
#'   science presences to be added to `r2`.
#' @return SpatRaster. Updated raster with citizen science presences added under
#'   the specified constraints.
#' @keywords internal
#' @noRd

update_citizen <- function(r1, r2, dist_km = 100) {

  if (!inherits(r1, "SpatRaster") || !inherits(r2, "SpatRaster")) {
    ecokit::stop_ctx(
      "r1 and r2 must be SpatRaster objects", include_backtrace = TRUE,
      class_r1 = class(r1), class_r2 = class(r2))
  }
  if (terra::nlyr(r1) != 1 || terra::nlyr(r2) != 1) {
    ecokit::stop_ctx(
      "r1 and r2 must be single layer SpatRaster objects",
      include_backtrace = TRUE,
      nlyr_r1 = terra::nlyr(r1), nlyr_r2 = terra::nlyr(r2))
  }
  if (!all(terra::ext(r1) == terra::ext(r2))) {
    ecokit::stop_ctx(
      "r1 and r2 must have the same extent",
      include_backtrace = TRUE, ext_r1 = as.character(terra::ext(r1)),
      ext_r2 = as.character(terra::ext(r2)))
  }
  if (!all(terra::res(r1) == terra::res(r2))) {
    ecokit::stop_ctx(
      "r1 and r2 must have the same resolution",
      include_backtrace = TRUE)
  }
  if (!all(unique(terra::values(r1)) %in% c(NaN, 0, 1))) {
    ecokit::stop_ctx(
      "r1 must be a binary raster with values of 0, 1, or NA",
      include_backtrace = TRUE)
  }
  if (!all(unique(terra::values(r2)) %in% c(NaN, 0, 1))) {
    ecokit::stop_ctx(
      "r2 must be a binary raster with values of 0, 1, or NA",
      include_backtrace = TRUE)
  }

  # Convert r1 and r2 to sf points
  r1_sf <- terra::classify(r1 == 1, cbind(0, NA)) %>%
    stats::setNames("r1") %>%
    as.data.frame(xy = TRUE) %>%
    dplyr::filter(r1 == 1) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(r1))

  r2_sf <- terra::classify(r2 == 1, cbind(0, NA)) %>%
    stats::setNames("r2") %>%
    as.data.frame(xy = TRUE) %>%
    dplyr::filter(r2 == 1) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(r1))

  # Identify points in r1_sf within dist_km of any point in r2_sf
  # but not within 500 m of any point in r2_sf
  points_to_add <- sf::st_filter(
    x = r1_sf, y = r2_sf,
    .predicate = sf::st_is_within_distance, dist = dist_km * 1000)
  if (nrow(points_to_add) == 0) {
    return(r2)
  }

  # exclude overlapped points with r2_sf
  points_to_add <- dplyr::filter(
    points_to_add,
    rowSums(
      sf::st_is_within_distance(
        points_to_add, r2_sf, dist = 500, sparse = FALSE)) == 0)
  if (nrow(points_to_add) == 0) {
    return(r2)
  }

  # Rasterize the points to add back to r2
  terra::rasterize(points_to_add, r2, field = 1, update = TRUE)
}
