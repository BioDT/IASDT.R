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

  if (is.null(species) || is.na(species) || species == "") {
    ecokit::stop_ctx(
      "`species` cannot be empty", species = species, include_backtrace = TRUE)
  }

  ecokit::info_chunk(species, verbose = verbose)

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::cat_time("Checking arguments", verbose = verbose)

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)
  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("species", "env_file"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "verbose")
  rm(AllArgs, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- Path_TaxaInfo <- Path_GBIF <- Path_eLTER <- Path_EASIN <-
    Path_BioReg <- Path_Grid_Ref <- Path_PA <- Path_TaxaInfo_RData <- BiogReg <-
    Country <- Species_name2 <- Species_name <- n <- GBIF <- EASIN <- eLTER <-
    PA <- `status-decision` <- Path_TaxaCNT <- country <-
    gbif_key <- status_decision <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables", verbose = verbose)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", TRUE, FALSE,
    "Path_EASIN", "DP_R_EASIN_processed", TRUE, FALSE,
    "Path_eLTER", "DP_R_eLTER_processed", FALSE, TRUE,
    "Path_PA", "DP_R_PA", FALSE, FALSE,
    "Path_TaxaCNT", "DP_R_Taxa_country", FALSE, TRUE,
    "Path_TaxaInfo_RData", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "Path_TaxaInfo", "DP_R_Taxa_info", FALSE, TRUE,
    "Path_BioReg", "DP_R_BioReg_processed", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Preparing input data ----
  ecokit::cat_time("Preparing input data", verbose = verbose)

  # # ................................ ###

  ## Current species info
  ecokit::cat_time("Current species info", level = 1L, verbose = verbose)

  Species2 <- ecokit::replace_space(species)
  IAS_ID <- readr::read_tsv(
    file = Path_TaxaInfo, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(Species_name2 == Species2)
  Sp_File <- unique(IAS_ID$Species_File)
  IAS_ID <- dplyr::pull(IAS_ID, "IAS_ID") %>%
    unique() %>%
    stringr::str_pad(width = 4, pad = "0")
  GBIF_Keys <- ecokit::load_as(Path_TaxaInfo_RData) %>%
    dplyr::filter(Species_name == species) %>%
    dplyr::pull("speciesKey")

  # # ................................ ###

  ## Check directories ----
  ecokit::cat_time("Check directories", level = 1L, verbose = verbose)

  ecokit::cat_time(
    "Check and create directories", level = 2L, verbose = verbose)
  Path_PA_Summary <- fs::path(Path_PA, "PA_summary")
  Path_PA_tif <- fs::path(Path_PA, "PA_tif")
  Path_PA_RData <- fs::path(Path_PA, "PA_RData")
  Path_PA_JPEG <- fs::path(Path_PA, "Distribution_JPEG")

  c(Path_PA_Summary, Path_PA_tif, Path_PA_RData, Path_PA_JPEG) %>%
    purrr::walk(fs::dir_create)

  Out_JPEG <- fs::path(Path_PA_JPEG, paste0(Sp_File, ".jpeg"))
  out_qs2 <- fs::path(Path_PA_Summary, paste0(Sp_File, ".qs2"))

  # # ................................ ###

  ecokit::cat_time(
    "Check path for EASIN and GBIF data", level = 2L, verbose = verbose)
  Path_GBIF_DT <- fs::path(Path_GBIF, "Sp_Data")
  Path_EASIN <- fs::path(Path_EASIN, "Sp_DT")

  if (!dir.exists(Path_GBIF_DT)) {
    ecokit::stop_ctx(
      "Required path for GBIF data do not exist", Path_GBIF_DT = Path_GBIF_DT,
      include_backtrace = TRUE)
  }

  if (!dir.exists(Path_EASIN)) {
    ecokit::stop_ctx(
      "Required path for EASIN data do not exist", Path_EASIN = Path_EASIN,
      include_backtrace = TRUE)
  }

  # # ................................ ###

  ecokit::cat_time("eLTER data", level = 1L, verbose = verbose)
  eLTER_DT <- ecokit::load_as(Path_eLTER) %>%
    dplyr::filter(Species_name == species)

  # # ................................ ###

  # Loading reference grids ----
  ecokit::cat_time("Loading reference grids", level = 1L, verbose = verbose)
  grid_100_file <- fs::path(Path_Grid_Ref, "Grid_100_sf.RData")
  if (!file.exists(grid_100_file)) {
    ecokit::stop_ctx(
      "grid file does not exist: Grid_100_sf.RData",
      grid_100_file = grid_100_file, include_backtrace = TRUE)
  }

  GridsPath <- fs::path(
    Path_Grid,
    c(
      "Grid_10_Land_Crop_sf_Country.RData", "Grid_10_Land_Crop.RData",
      "Grid_10_Land_sf.RData"))

  if (!all(file.exists(GridsPath))) {
    ecokit::stop_ctx(
      paste0(
        "grid files do not exist: \n  >>> ",
        paste(GridsPath[!file.exists(GridsPath)], collapse = "\n  >>> ")),
      GridsPath = GridsPath,
      missing_GridsPath = GridsPath[!file.exists(GridsPath)],
      include_backtrace = TRUE)
  }

  ### sf - 100 km ----
  ecokit::cat_time("sf - 100 km", level = 2L, verbose = verbose)
  Grid_100_sf <- ecokit::load_as(grid_100_file) %>%
    magrittr::extract2("Grid_100_sf_s")

  ### sf - 10 km ----
  ecokit::cat_time("sf - 10 km - CNT", level = 2L, verbose = verbose)
  Grid_10_CNT <- fs::path(
    Path_Grid, "Grid_10_Land_Crop_sf_Country.RData") %>%
    ecokit::load_as()

  ### raster - 10 km ----
  ecokit::cat_time("raster - 10 km", level = 2L, verbose = verbose)
  RefGrid <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData") %>%
    ecokit::load_as(unwrap_r = TRUE)

  ### raster - 100 km ----
  ecokit::cat_time("raster - 100 km", level = 2L, verbose = verbose)
  Grid_100_Land <- fs::path(Path_Grid, "Grid_10_Land_sf.RData") %>%
    ecokit::load_as() %>%
    sf::st_geometry() %>%
    sf::st_centroid() %>%
    sf::st_filter(x = Grid_100_sf, join = sf::st_within)

  Grid100Empty <- dplyr::slice(Grid_100_sf, 0)

  rm(Grid_100_sf, envir = environment())
  invisible(gc())

  # # ................................ ###

  ## `PA_100km` function - species distribution at 100 km ----
  ecokit::cat_time(
    "`PA_100km` function - species distribution at 100 km", level = 1L,
    verbose = verbose)
  PA_100km <- function(Sp_Sf, Grid100 = Grid_100_Land) {
    sf::st_filter(x = Grid100, y = Sp_Sf, join = sf::st_within) %>%
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
  Countries2Exclude <- readxl::read_xlsx(path = Path_TaxaCNT, sheet = 1) %>%
    dplyr::rename(status_decision = `status-decision`) %>%
    # filter data on the current species and only keep countries with `cult+cas`
    # or `delete` status
    dplyr::filter(
      gbif_key %in% GBIF_Keys, status_decision %in% c("cult+cas", "delete")
    ) %>%
    # rename countries to match country names of the countries boundaries
    dplyr::mutate(
      country = dplyr::case_when(
        country == "United_Kingdom" ~ "United Kingdom",
        country == "North_Macedonia" ~ "North Macedonia",
        country == "Isle_of_Man" ~ "Isle of Man",
        country == "Bosnia_and_Herzegovina" ~ "Bosnia and Herzegovina",
        .default = country)
    ) %>%
    # Exclude Turkey, as it is currently not in the current study area
    dplyr::filter(country != "Turkey") %>%
    dplyr::pull("country")

  ecokit::cat_time(
    paste0("There are ", length(Countries2Exclude), " countries to exclude:"),
    level = 2L, cat_timestamp = FALSE, verbose = verbose)

  # Mask grid to exclude countries - `TRUE` for grid cells to be considered as
  # presence if present in any of the data source; `FALSE` for grid cells need
  # to be masked as 0 in species distribution maps (1 becomes 0)
  if (length(Countries2Exclude) > 0) {

    ecokit::cat_time(
      paste(sort(Countries2Exclude), collapse = " + "),
      level = 3L, cat_timestamp = FALSE, verbose = verbose)

    Mask_Keep <- Grid_10_CNT %>%
      dplyr::mutate(Keep = !(Country %in% Countries2Exclude)) %>%
      dplyr::select("Keep") %>%
      terra::rasterize(y = RefGrid, field = "Keep") %>%
      terra::as.bool() %>%
      terra::toMemory() %>%
      stats::setNames("Mask_Keep")
  } else {
    Mask_Keep <- terra::as.bool(RefGrid) %>%
      terra::toMemory() %>%
      stats::setNames("Mask_Keep")
  }

  rm(Grid_10_CNT, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Preparing species data from the 3 sources -----
  ecokit::cat_time(
    "Preparing species data from the 3 sources", verbose = verbose)

  ## 1. GBIF -----
  ecokit::cat_time("1. GBIF", level = 1L, verbose = verbose)

  GBIF_Keys <- paste(GBIF_Keys, collapse = "_")

  # Path for the current species GBIF data
  Path_GBIF_D <- fs::path(Path_GBIF_DT, paste0(Sp_File, ".RData"))
  Path_GBIF_R <- fs::path(
    Path_GBIF, "Sp_Raster", paste0(Sp_File, "_Raster.RData"))

  if (all(file.exists(Path_GBIF_D, Path_GBIF_R))) {
    # GBIF data as sf object
    GBIF_DT <- ecokit::load_as(Path_GBIF_D)

    # GBIF data as raster map
    GBIF_R <- ecokit::load_as(Path_GBIF_R) %>%
      terra::unwrap() %>%
      # convert raster map into binary (1/0)
      ecokit::raster_to_pres_abs() %>%
      terra::mask(RefGrid) %>%
      stats::setNames("GBIF")

    # presence 100 km grid
    GBIF_Gr100 <- PA_100km(GBIF_DT)

    rm(GBIF_DT, envir = environment())
    invisible(gc())

    # Map to export out of this function
    GBIF_R_Out <- stats::setNames(GBIF_R, Sp_File) %>%
      terra::wrap() %>%
      list()

    GBIF_Path <- Path_GBIF_D
  } else {
    GBIF_R <- terra::classify(RefGrid, cbind(1, 0)) %>%
      stats::setNames("GBIF")

    GBIF_R_Out <- stats::setNames(GBIF_R, Sp_File) %>%
      terra::wrap() %>%
      list()

    GBIF_Gr100 <- Grid100Empty
    GBIF_Path <- NA_character_
  }

  # # .................................... ###

  ## 2. EASIN -----
  ecokit::cat_time("2. EASIN", level = 1L, verbose = verbose)

  Path_EASIN_D <- fs::path(Path_EASIN, paste0(Sp_File, "_DT.RData"))

  if (file.exists(Path_EASIN_D)) {
    EASIN_DT <- ecokit::load_as(Path_EASIN_D)

    # raster map
    EASIN_R <- dplyr::select(EASIN_DT, "Species_name") %>%
      terra::rasterize(RefGrid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(RefGrid) %>%
      stats::setNames("EASIN")

    EASIN_Gr100 <- PA_100km(EASIN_DT)

    rm(EASIN_DT, envir = environment())
    invisible(gc())

    EASIN_R_Out <- stats::setNames(EASIN_R, Sp_File) %>%
      terra::wrap() %>%
      list()

    EASIN_Path <- Path_EASIN_D
  } else {
    EASIN_R <- terra::classify(RefGrid, cbind(1, 0)) %>%
      stats::setNames("EASIN")

    EASIN_R_Out <- stats::setNames(EASIN_R, Sp_File) %>%
      terra::wrap() %>%
      list()

    EASIN_Gr100 <- Grid100Empty
    EASIN_Path <- NA_character_
  }

  # # .................................... ###

  ## 3. eLTER -----
  ecokit::cat_time("3. eLTER", level = 1L, verbose = verbose)

  if (nrow(eLTER_DT) > 0) {
    eLTER_R <- dplyr::select(eLTER_DT, "Species_name") %>%
      terra::rasterize(RefGrid) %>%
      ecokit::raster_to_pres_abs() %>%
      terra::mask(RefGrid) %>%
      stats::setNames("eLTER")

    eLTER_Gr100 <- PA_100km(eLTER_DT)

    eLTER_R_Out <- stats::setNames(eLTER_R, Sp_File) %>%
      terra::wrap() %>%
      list()
  } else {
    eLTER_R <- terra::classify(RefGrid, cbind(1, 0)) %>%
      stats::setNames("eLTER")

    eLTER_R_Out <- stats::setNames(eLTER_R, Sp_File) %>%
      terra::wrap() %>%
      list()

    eLTER_Gr100 <- Grid100Empty
  }

  rm(eLTER_DT, Grid_100_Land, Grid100Empty, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ## 4. Merging data from the 3 data sources -----
  ecokit::cat_time(
    "Merging data from the 3 data sources", level = 1L, verbose = verbose)

  Sp_PA <- sum(GBIF_R, EASIN_R, eLTER_R) %>%
    ecokit::raster_to_pres_abs() %>%
    stats::setNames("PA") %>%
    terra::mask(RefGrid)
  Sp_PA$PA_Masked <- (Sp_PA$PA * Mask_Keep)
  Sp_PA <- c(GBIF_R, EASIN_R, eLTER_R, Sp_PA, Mask_Keep) %>%
    # Ensure that values are read from memory
    terra::toMemory()

  rm(RefGrid, Mask_Keep, envir = environment())
  invisible(gc())

  # number of cells with values
  PA_NCells_All <- terra::global(Sp_PA$PA == 1, sum, na.rm = TRUE) %>%
    as.integer()
  PA_NCells_Naturalized <- terra::global(
    x = Sp_PA$PA_Masked == 1, sum, na.rm = TRUE) %>%
    as.integer()

  # # ..................................................................... ###

  # Processing species data ----
  ecokit::cat_time("Processing species data", verbose = verbose)

  # If there is no presence grid cell for the current species, return NA early
  if (PA_NCells_All == 0) {
    return(tibble::tibble(species = species, PA_summary = NA_character_))
  }

  # # .................................... ###

  ## Save maps ------

  ecokit::cat_time("Save species data", level = 1L, verbose = verbose)

  ### RData -----
  ecokit::cat_time("`.RData`", level = 2L, verbose = verbose)
  path_RData <- fs::path(Path_PA_RData, paste0(Sp_File, "_PA.RData"))
  ecokit::save_as(
    object = terra::wrap(Sp_PA), object_name = paste0(Sp_File, "_PA"),
    out_path = path_RData)

  ### tif - all presence grid cells -----
  ecokit::cat_time(
    "`.tif` - all presence grid cells", level = 2L, verbose = verbose)
  terra::writeRaster(
    x = Sp_PA$PA, overwrite = TRUE,
    filename = fs::path(Path_PA_tif, paste0(Sp_File, "_All.tif")))

  ### tif - Excluding cultivated or casual observations -----
  ecokit::cat_time(
    "`.tif` - Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)
  terra::writeRaster(
    x = Sp_PA$PA_Masked, overwrite = TRUE,
    filename = fs::path(Path_PA_tif, paste0(Sp_File, "_Masked.tif")))

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

  BioReg_R <- fs::path(Path_BioReg, "BioReg_R.RData")
  if (isFALSE(ecokit::check_rdata(BioReg_R))) {
    ecokit::stop_ctx(
      "Required file for biogeographical regions does not exist",
      BioReg_R = BioReg_R, include_backtrace = TRUE)
  }
  BioReg_R <- ecokit::load_as(BioReg_R, unwrap_r = TRUE)

  # name of Biogeographical regions
  BioReg_Names <- c(
    "alpine", "arctic", "atlantic", "blackSea", "boreal",
    "continental", "mediterranean", "pannonian", "steppic") %>%
    paste0("BioReg_", .)

  Sp_BiogeoRegions <- terra::classify(Sp_PA$PA, cbind(0, NA)) %>%
    terra::mask(mask = ., x = BioReg_R) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(BiogReg = "short_name") %>%
    dplyr::count(BiogReg) %>%
    tidyr::pivot_wider(names_from = BiogReg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("BioReg_", .x)) %>%
    ecokit::add_missing_columns(0L, BioReg_Names) %>%
    # Sort columns
    dplyr::select(tidyselect::all_of(BioReg_Names))

  BioRegsSumm_Min <- min(Sp_BiogeoRegions, na.rm = TRUE)
  BioRegsSumm_Max <- max(Sp_BiogeoRegions, na.rm = TRUE)
  BioRegsSumm_Mean <- unlist(Sp_BiogeoRegions) %>%
    mean(na.rm = TRUE) %>%
    round(1)
  BioRegsSumm_N <- sum(Sp_BiogeoRegions > 0)

  ### Masked -----
  ecokit::cat_time(
    "Excluding cultivated or casual observations",
    level = 2L, verbose = verbose)

  BioReg_Names2 <- stringr::str_replace(
    string = BioReg_Names, pattern = "BioReg_", replacement = "BioReg_Masked_")

  Sp_BiogeoRegions_Masked <- terra::classify(Sp_PA$PA_Masked, cbind(0, NA)) %>%
    terra::mask(mask = ., x = BioReg_R) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    dplyr::rename(BiogReg = "short_name") %>%
    dplyr::count(BiogReg) %>%
    tidyr::pivot_wider(names_from = BiogReg, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("BioReg_Masked_", .x)) %>%
    ecokit::add_missing_columns(0L, BioReg_Names2) %>%
    dplyr::select(tidyselect::all_of(BioReg_Names2))


  rm(BioReg_R, envir = environment())
  invisible(gc())

  if (nrow(Sp_BiogeoRegions_Masked) > 0) {
    BioRegsMaskSumm_Min <- min(Sp_BiogeoRegions_Masked, na.rm = TRUE)
    BioRegsMaskSumm_Max <- max(Sp_BiogeoRegions_Masked, na.rm = TRUE)
    BioRegsMaskSumm_Mean <- unlist(Sp_BiogeoRegions_Masked) %>%
      mean(na.rm = TRUE) %>%
      round(1)
    BioRegsMaskSumm_N <- sum(Sp_BiogeoRegions_Masked > 0)
  } else {
    BioRegsMaskSumm_Min <- 0
    BioRegsMaskSumm_Max <- 0
    BioRegsMaskSumm_Mean <- 0
    BioRegsMaskSumm_N <- 0
  }

  # # .................................... ###

  ## Number of presence grid cells per data provider -----
  ecokit::cat_time(
    "# presence grid cells per data provider", level = 1L, verbose = verbose)

  # presence grid cells per data type
  RVals <- c(GBIF_R, EASIN_R, eLTER_R) %>%
    as.data.frame(na.rm = TRUE) %>%
    tibble::tibble() %>%
    # Keep only grid cells present for any data provider
    dplyr::filter_all(dplyr::any_vars(. > 0))

  RVals_Masked <- c(GBIF_R, EASIN_R, eLTER_R) %>%
    magrittr::multiply_by(Sp_PA$PA_Masked) %>%
    as.data.frame(na.rm = TRUE) %>%
    tibble::tibble() %>%
    dplyr::filter_all(dplyr::any_vars(. > 0))

  N_GBIF <- as.integer(terra::global(GBIF_R, sum, na.rm = TRUE))
  N_GBIF_Unique <- dplyr::filter(RVals, GBIF == 1, EASIN == 0, eLTER == 0) %>%
    dplyr::pull("GBIF") %>%
    sum()

  N_GBIF_Masked <- terra::global(
    (Sp_PA$GBIF * Sp_PA$Mask_Keep), sum, na.rm = TRUE) %>%
    as.integer()
  N_GBIF_Unique_Masked <- dplyr::filter(
    RVals_Masked, GBIF == 1, EASIN == 0, eLTER == 0) %>%
    dplyr::pull("GBIF") %>%
    sum()

  N_EASIN <- as.integer(terra::global(EASIN_R, sum, na.rm = TRUE))
  N_EASIN_Unique <- dplyr::filter(RVals, GBIF == 0, EASIN == 1, eLTER == 0) %>%
    dplyr::pull("EASIN") %>%
    sum()

  N_EASIN_Masked <- terra::global(
    (Sp_PA$EASIN * Sp_PA$Mask_Keep), sum, na.rm = TRUE) %>%
    as.integer()
  N_EASIN_Unique_Masked <- dplyr::filter(
    RVals_Masked, GBIF == 0, EASIN == 1, eLTER == 0) %>%
    dplyr::pull("EASIN") %>%
    sum()

  N_eLTER <- as.integer(terra::global(eLTER_R, sum, na.rm = TRUE))
  N_eLTER_Unique <- dplyr::filter(RVals, GBIF == 0, EASIN == 0, eLTER == 1) %>%
    dplyr::pull("eLTER") %>%
    sum()

  N_eLTER_Masked <- terra::global(
    (Sp_PA$eLTER * Sp_PA$Mask_Keep), sum, na.rm = TRUE) %>%
    as.integer()
  N_eLTER_Unique_Masked <- dplyr::filter(
    RVals_Masked, GBIF == 0, EASIN == 0, eLTER == 1) %>%
    dplyr::pull("eLTER") %>%
    sum()

  rm(RVals, RVals_Masked, GBIF_R, EASIN_R, eLTER_R, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## Number of presence grid cells per country -----
  ecokit::cat_time(
    "# presence grid cells per country", level = 1L, verbose = verbose)

  CountryList <- c(
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

  Grid_CNT <- fs::path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData") %>%
    ecokit::load_as() %>%
    dplyr::select("Country")

  SpCountry <- terra::as.points(Sp_PA$PA) %>%
    sf::st_as_sf() %>%
    dplyr::filter(PA == 1) %>%
    dplyr::select(-"PA") %>%
    sf::st_join(Grid_CNT) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble() %>%
    dplyr::mutate(Country = ecokit::replace_space(Country)) %>%
    dplyr::count(Country) %>%
    tidyr::pivot_wider(names_from = Country, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("CNT_", .x)) %>%
    ecokit::add_missing_columns(0L, CountryList) %>%
    dplyr::select(tidyselect::all_of(CountryList))

  rm(Grid_CNT, envir = environment())
  invisible(gc())

  # # .................................... ###

  ## Number of unique iNaturalist grid cells -----
  ecokit::cat_time(
    "# unique iNaturalist grid cells", level = 1L, verbose = verbose)

  iNatur_DT <- fs::path(Path_GBIF, "iNaturalist_Count.RData") %>%
    ecokit::load_as() %>%
    dplyr::filter(species == !!species)

  if (nrow(iNatur_DT) > 0) {
    iNatur_Unique <- iNatur_DT$iNaturalist_Unique
    iNatur_Perc <- round(100 * (iNatur_DT$iNaturalist_Unique / N_GBIF), 1)
  } else {
    iNatur_Unique <- iNatur_Perc <- 0
  }

  rm(iNatur_DT, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Prepare and export species summary info -------
  ecokit::cat_time("Prepare and export species summary info", verbose = verbose)

  Binary_R_Out <- stats::setNames(Sp_PA$PA, Sp_File) %>%
    terra::wrap() %>%
    list()

  Binary_R_Masked_Out <- stats::setNames(Sp_PA$PA_Masked, Sp_File) %>%
    terra::wrap() %>%
    list()

  Binary_R_Mask_Keep <- stats::setNames(Sp_PA$Mask_Keep, Sp_File) %>%
    terra::wrap() %>%
    list()

  rm(Sp_PA, envir = environment())
  invisible(gc())

  IntegerCols <- c(
    "NCells_All", "NCells_Naturalized",
    "GBIF", "GBIF_Unique", "GBIF_Masked", "GBIF_Masked_Unique",
    "EASIN", "EASIN_Unique", "EASIN_Masked", "EASIN_Masked_Unique",
    "eLTER", "eLTER_Unique", "eLTER_Masked", "eLTER_Masked_Unique",
    "BioRegsSumm_Min", "BioRegsSumm_Max", "BioRegsSumm_N",
    "BioRegsMaskSumm_Min", "BioRegsMaskSumm_Max", "BioRegsMaskSumm_N",
    "iNaturalist_Unique")

  Results <- tibble::tibble(
    Species = species,
    GBIF_Keys = GBIF_Keys,
    species_ID = IAS_ID,
    NCells_All = PA_NCells_All,
    NCells_Naturalized = PA_NCells_Naturalized,

    GBIF = N_GBIF,
    GBIF_Unique = N_GBIF_Unique,
    GBIF_Masked = N_GBIF_Masked,
    GBIF_Masked_Unique = N_GBIF_Unique_Masked,
    GBIF_Path = GBIF_Path,
    GBIF_Gr100 = list(GBIF_Gr100),
    GBIF_R = GBIF_R_Out,

    EASIN = N_EASIN,
    EASIN_Unique = N_EASIN_Unique,
    EASIN_Masked = N_EASIN_Masked,
    EASIN_Masked_Unique = N_EASIN_Unique_Masked,
    EASIN_Path = EASIN_Path,
    EASIN_Gr100 = list(EASIN_Gr100),
    EASIN_R = EASIN_R_Out,

    eLTER = N_eLTER,
    eLTER_Unique = N_eLTER_Unique,
    eLTER_Masked = N_eLTER_Masked,
    eLTER_Masked_Unique = N_eLTER_Unique_Masked,
    eLTER_Gr100 = list(eLTER_Gr100),
    eLTER_R = eLTER_R_Out,

    PA_Map = Binary_R_Out,
    PA_Masked_Map = Binary_R_Masked_Out,
    Mask_Keep = Binary_R_Mask_Keep,

    SpCountry = list(SpCountry),
    BioRegs_DT = list(Sp_BiogeoRegions),
    BioRegsSumm_Min = BioRegsSumm_Min,
    BioRegsSumm_Max = BioRegsSumm_Max,
    BioRegsSumm_Mean = BioRegsSumm_Mean,
    BioRegsSumm_N = BioRegsSumm_N,
    BioRegsMask_DT = list(Sp_BiogeoRegions_Masked),
    BioRegsMaskSumm_Min = BioRegsMaskSumm_Min,
    BioRegsMaskSumm_Max = BioRegsMaskSumm_Max,
    BioRegsMaskSumm_Mean = BioRegsMaskSumm_Mean,
    BioRegsMaskSumm_N = BioRegsMaskSumm_N,
    iNaturalist_Unique = iNatur_Unique,
    iNatur_Perc = iNatur_Perc,
    Countries2Exclude = list(Countries2Exclude),
    path_JPEG = Out_JPEG) %>%
    dplyr::mutate_at(IntegerCols, ~ as.integer(.))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing species data was finished in ", ... = "\n",
    verbose = verbose)

  # save species data
  ecokit::save_as(object = Results, out_path = out_qs2)

  return(tibble::tibble(species = species, PA_summary = out_qs2))

}
