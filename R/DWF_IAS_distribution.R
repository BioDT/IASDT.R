# # |------------------------------------------------------------------------| #
# IAS_distribution ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name IAS_data
#' @rdname IAS_data
#' @order 2

IAS_distribution <- function(
    species = NULL, env_file = ".env", verbose = FALSE, overwrite = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(species) || is.na(species) || species == "") {
    stop("`species` cannot be empty", call. = FALSE)
  }

  if (isFALSE(verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  IASDT.R::info_chunk(species)

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::cat_time("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("species", "env_file"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("verbose", "overwrite"))

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
  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
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
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Preparing input data ----
  IASDT.R::cat_time("Preparing input data")

  # # ................................ ###

  ## Current species info
  IASDT.R::cat_time("Current species info", level = 1)

  Species2 <- IASDT.R::replace_space(species)
  IAS_ID <- readr::read_tsv(
    file = Path_TaxaInfo, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(Species_name2 == Species2)
  Sp_File <- unique(IAS_ID$Species_File)
  IAS_ID <- dplyr::pull(IAS_ID, "IAS_ID") %>%
    unique() %>%
    stringr::str_pad(width = 4, pad = "0")
  GBIF_Keys <- IASDT.R::load_as(Path_TaxaInfo_RData) %>%
    dplyr::filter(Species_name == species) %>%
    dplyr::pull("speciesKey")

  # # ................................ ###

  ## Check directories ----
  IASDT.R::cat_time("Check directories", level = 1)

  IASDT.R::cat_time("Check/create directories", level = 2)
  Path_PA_Summary <- IASDT.R::path(Path_PA, "SpSummary")
  Path_PA_tif <- IASDT.R::path(Path_PA, "tif")
  Path_PA_RData <- IASDT.R::path(Path_PA, "RData")
  Path_PA_JPEG <- IASDT.R::path(Path_PA, "JPEG_Maps")
  Path_PA_All <- c(
    Path_PA_Summary, Path_PA_tif, Path_PA_RData, Path_PA_JPEG)
  Path_PA_Missing <- !all(dir.exists(Path_PA_All))
  if (Path_PA_Missing) {
    Missing <- which(!dir.exists(Path_PA_All))
    fs::dir_create(Path_PA_All[Missing])
  }

  Out_PA <- IASDT.R::path(Path_PA_RData, paste0(Sp_File, "_PA.RData"))
  Out_Summary <- IASDT.R::path(
    Path_PA_Summary, paste0(Sp_File, "_Summary.RData"))
  Out_tif_All <- IASDT.R::path(Path_PA_tif, paste0(Sp_File, "_All.tif"))
  Out_tif_Masked <- IASDT.R::path(Path_PA_tif, paste0(Sp_File, "_Masked.tif"))
  Out_JPEG <- IASDT.R::path(Path_PA_JPEG, paste0(Sp_File, ".jpeg"))
  Out_Exists <- c(
    Out_PA, Out_Summary, Out_tif_All, Out_tif_Masked, Out_JPEG) %>%
    file.exists() %>%
    all()
  if (isFALSE(overwrite) && Out_Exists) {
    return(NULL)
  }

  # # ................................ ###

  IASDT.R::cat_time("Check path for EASIN and GBIF data", level = 2)
  Path_GBIF_DT <- IASDT.R::path(Path_GBIF, "Sp_Data")
  Path_EASIN <- IASDT.R::path(Path_EASIN, "Sp_DT")

  if (!dir.exists(Path_GBIF_DT)) {
    stop(
      "Required path for GBIF data do not exist: ", Path_GBIF_DT, call. = FALSE)
  }

  if (!dir.exists(Path_EASIN)) {
    stop(
      "Required path for EASIN data do not exist: ", Path_EASIN, call. = FALSE)
  }

  # # ................................ ###

  IASDT.R::cat_time("eLTER data", level = 1)
  eLTER_DT <- IASDT.R::load_as(Path_eLTER) %>%
    dplyr::filter(Species_name == species)

  # # ................................ ###

  # Loading reference grids ----
  IASDT.R::cat_time("Loading reference grids", level = 1)
  if (!file.exists(IASDT.R::path(Path_Grid_Ref, "Grid_100_sf.RData"))) {
    stop(
      "The following grid file does not exist: Grid_100_sf.RData",
      call. = FALSE)
  }

  GridsPath <- IASDT.R::path(
    Path_Grid,
    c(
      "Grid_10_Land_Crop_sf_Country.RData", "Grid_10_Land_Crop.RData",
      "Grid_10_Land_sf.RData"))

  if (!all(file.exists(GridsPath))) {
    stop(
      "The following grid files do not exist: \n  >>> ",
      paste(GridsPath[!file.exists(GridsPath)], collapse = "\n  >>> "),
      call. = FALSE)
  }

  ### sf - 100 km ----
  IASDT.R::cat_time("sf - 100 km", level = 2)
  Grid_100_sf <- IASDT.R::path(Path_Grid_Ref, "Grid_100_sf.RData") %>%
    IASDT.R::load_as() %>%
    magrittr::extract2("Grid_100_sf_s")

  ### sf - 10 km ----
  IASDT.R::cat_time("sf - 10 km - CNT", level = 2)
  Grid_10_CNT <- IASDT.R::path(
    Path_Grid, "Grid_10_Land_Crop_sf_Country.RData") %>%
    IASDT.R::load_as()

  ### raster - 10 km ----
  IASDT.R::cat_time("raster - 10 km", level = 2)
  RefGrid <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData") %>%
    IASDT.R::load_as() %>%
    terra::unwrap()

  ### raster - 100 km ----
  IASDT.R::cat_time("raster - 100 km", level = 2)
  Grid_100_Land <- IASDT.R::path(Path_Grid, "Grid_10_Land_sf.RData") %>%
    IASDT.R::load_as() %>%
    sf::st_geometry() %>%
    sf::st_centroid() %>%
    sf::st_filter(x = Grid_100_sf, join = sf::st_within)

  Grid100Empty <- dplyr::slice(Grid_100_sf, 0)

  # # ................................ ###

  ## `PA_100km` function - species distribution at 100 km ----
  IASDT.R::cat_time(
    "`PA_100km` function - species distribution at 100 km", level = 1)
  PA_100km <- function(Sp_Sf, Grid100 = Grid_100_Land) {
    sf::st_filter(x = Grid100, y = Sp_Sf, join = sf::st_within) %>%
      dplyr::select("geometry") %>%
      tibble::tibble() %>%
      dplyr::distinct() %>%
      sf::st_as_sf()
  }

  # # ................................ ###

  ## Exclude countries with only cultivated or casual observations ----
  IASDT.R::cat_time(
    "Exclude countries with only cultivated or casual observations",
    level = 1)

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

  IASDT.R::cat_time(
    paste0("There are ", length(Countries2Exclude), " countries to exclude:"),
    level = 2, time = FALSE)

  # Mask grid to exclude countries - `TRUE` for grid cells to be considered as
  # presence if present in any of the data source; `FALSE` for grid cells need
  # to be masked as 0 in species distribution maps (1 becomes 0)
  if (length(Countries2Exclude) > 0) {

    IASDT.R::cat_time(
      paste(sort(Countries2Exclude), collapse = " + "), level = 3, time = FALSE)

    Mask_Keep <- Grid_10_CNT %>%
      dplyr::mutate(Keep = !(Country %in% Countries2Exclude)) %>%
      dplyr::select("Keep") %>%
      terra::rasterize(y = RefGrid, field = "Keep") %>%
      terra::as.bool() %>%
      IASDT.R::set_raster_values() %>%
      stats::setNames("Mask_Keep")
  } else {
    Mask_Keep <- terra::as.bool(RefGrid) %>%
      IASDT.R::set_raster_values() %>%
      stats::setNames("Mask_Keep")
  }

  rm(Grid_10_CNT, envir = environment())

  GBIF_Keys <- paste(GBIF_Keys, collapse = "_")

  # # ..................................................................... ###

  # Preparing species data from the 3 sources -----
  IASDT.R::cat_time("Preparing species data from the 3 sources")

  ## 1. GBIF -----
  IASDT.R::cat_time("1. GBIF", level = 1)

  # Path for the current species GBIF data
  Path_GBIF_D <- IASDT.R::path(Path_GBIF_DT, paste0(Sp_File, ".RData"))
  Path_GBIF_R <- IASDT.R::path(
    Path_GBIF, "Sp_Raster", paste0(Sp_File, "_Raster.RData"))

  if (all(file.exists(Path_GBIF_D, Path_GBIF_R))) {
    # GBIF data as sf object
    GBIF_DT <- IASDT.R::load_as(Path_GBIF_D)

    # GBIF data as raster map
    GBIF_R <- IASDT.R::load_as(Path_GBIF_R) %>%
      terra::unwrap() %>%
      # convert raster map into binary (1/0)
      IASDT.R::raster_to_PA() %>%
      terra::mask(RefGrid) %>%
      stats::setNames("GBIF")

    # presence 100 km grid
    GBIF_Gr100 <- PA_100km(GBIF_DT)

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
  IASDT.R::cat_time("2. EASIN", level = 1)

  Path_EASIN_D <- IASDT.R::path(Path_EASIN, paste0(Sp_File, "_DT.RData"))

  if (file.exists(Path_EASIN_D)) {
    EASIN_DT <- IASDT.R::load_as(Path_EASIN_D)

    # raster map
    EASIN_R <- dplyr::select(EASIN_DT, "Species_name") %>%
      terra::rasterize(RefGrid) %>%
      IASDT.R::raster_to_PA() %>%
      terra::mask(RefGrid) %>%
      stats::setNames("EASIN")

    EASIN_Gr100 <- PA_100km(EASIN_DT)

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
  IASDT.R::cat_time("3. eLTER", level = 1)

  if (nrow(eLTER_DT) > 0) {
    eLTER_R <- dplyr::select(eLTER_DT, "Species_name") %>%
      terra::rasterize(RefGrid) %>%
      IASDT.R::raster_to_PA() %>%
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

  # # ..................................................................... ###

  ## 4. Merging data from the 3 data sources -----
  IASDT.R::cat_time("Merging data from the 3 data sources", level = 1)

  Sp_PA <- sum(GBIF_R, EASIN_R, eLTER_R) %>%
    IASDT.R::raster_to_PA() %>%
    stats::setNames("PA") %>%
    terra::mask(RefGrid)
  Sp_PA$PA_Masked <- (Sp_PA$PA * Mask_Keep)
  Sp_PA <- c(GBIF_R, EASIN_R, eLTER_R, Sp_PA, Mask_Keep) %>%
    # Ensure that values are read from memory
    IASDT.R::set_raster_values()

  # number of cells with values
  PA_NCells_All <- terra::global(Sp_PA$PA == 1, sum, na.rm = TRUE) %>%
    as.integer()
  PA_NCells_Naturalized <- terra::global(
    x = Sp_PA$PA_Masked == 1, sum, na.rm = TRUE) %>%
    as.integer()

  # # ..................................................................... ###

  # Processing species data ----
  IASDT.R::cat_time("Processing species data")

  # If there is no presence grid cell for the current species, return NULL early
  if (PA_NCells_All == 0) {
    return(invisible(NULL))
  }

  # # .................................... ###

  ## Save maps ------

  IASDT.R::cat_time("Save species data", level = 1)

  ### RData -----
  IASDT.R::cat_time("`.RData`", level = 2)
  path_RData <- IASDT.R::path(Path_PA_RData, paste0(Sp_File, "_PA.RData"))
  IASDT.R::save_as(
    object = terra::wrap(Sp_PA), object_name = paste0(Sp_File, "_PA"),
    out_path = path_RData)

  ### tif - all presence grid cells -----
  IASDT.R::cat_time("`.tif` - all presence grid cells", level = 2)
  terra::writeRaster(
    x = Sp_PA$PA, overwrite = TRUE,
    filename = IASDT.R::path(Path_PA_tif, paste0(Sp_File, "_All.tif"))
  )

  ### tif - Excluding cultivated or casual observations -----
  IASDT.R::cat_time(
    "`.tif` - Excluding cultivated or casual observations", level = 2)
  terra::writeRaster(
    x = Sp_PA$PA_Masked, overwrite = TRUE,
    filename = IASDT.R::path(Path_PA_tif, paste0(Sp_File, "_Masked.tif")))

  ### NetCDF -----
  # NOT IMPLEMENTED YET

  # # .................................... ###

  ## Biogeographical regions ----
  IASDT.R::cat_time("Analysis per biogeographical regions", level = 1)

  ### All -----
  IASDT.R::cat_time("all data", level = 2)

  # Number of grid per each biogeographical region
  #
  # number of biogeographical regions per species and minimum / maximum / mean
  # number of grid cells per biogeographical regions

  BioReg_R <- IASDT.R::path(Path_BioReg, "BioReg_R.RData")
  if (isFALSE(IASDT.R::check_RData(BioReg_R))) {
    stop(
      "Required file for biogeographical regions does not exist: ", BioReg_R,
      call. = FALSE)
  }

  BioReg_R <- terra::unwrap(IASDT.R::load_as(BioReg_R))

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
    IASDT.R::add_missing_columns(0L, BioReg_Names) %>%
    # Sort columns
    dplyr::select(tidyselect::all_of(BioReg_Names))

  BioRegsSumm_Min <- min(Sp_BiogeoRegions, na.rm = TRUE)
  BioRegsSumm_Max <- max(Sp_BiogeoRegions, na.rm = TRUE)
  BioRegsSumm_Mean <- unlist(Sp_BiogeoRegions) %>%
    mean(na.rm = TRUE) %>%
    round(1)
  BioRegsSumm_N <- sum(Sp_BiogeoRegions > 0)

  ### Masked -----
  IASDT.R::cat_time("Excluding cultivated or casual observations", level = 2)

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
    IASDT.R::add_missing_columns(0L, BioReg_Names2) %>%
    dplyr::select(tidyselect::all_of(BioReg_Names2))

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
  IASDT.R::cat_time("# presence grid cells per data provider", level = 1)

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

  # # .................................... ###

  ## Number of presence grid cells per country -----
  IASDT.R::cat_time("# presence grid cells per country", level = 1)

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
    IASDT.R::replace_space() %>%
    paste0("CNT_", .)

  Grid_CNT <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData") %>%
    IASDT.R::load_as() %>%
    dplyr::select("Country")

  SpCountry <- terra::as.points(Sp_PA$PA) %>%
    sf::st_as_sf() %>%
    dplyr::filter(PA == 1) %>%
    dplyr::select(-"PA") %>%
    sf::st_join(Grid_CNT) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble() %>%
    dplyr::mutate(Country = IASDT.R::replace_space(Country)) %>%
    dplyr::count(Country) %>%
    tidyr::pivot_wider(names_from = Country, values_from = n) %>%
    dplyr::rename_all(~ stringr::str_c("CNT_", .x)) %>%
    IASDT.R::add_missing_columns(0L, CountryList) %>%
    dplyr::select(tidyselect::all_of(CountryList))

  # # .................................... ###

  ## Number of unique iNaturalist grid cells -----
  IASDT.R::cat_time("# unique iNaturalist grid cells", level = 1)

  iNatur_DT <- IASDT.R::path(Path_GBIF, "iNaturalist_Count.RData") %>%
    IASDT.R::load_as() %>%
    dplyr::filter(species == !!species)

  if (nrow(iNatur_DT) > 0) {
    iNatur_Unique <- iNatur_DT$iNaturalist_Unique
    iNatur_Perc <- round(100 * (iNatur_DT$iNaturalist_Unique / N_GBIF), 1)
  } else {
    iNatur_Unique <- iNatur_Perc <- 0
  }

  # # ..................................................................... ###

  # Prepare/export species summary info -------
  IASDT.R::cat_time("Prepare/export species summary info")

  Binary_R_Out <- stats::setNames(Sp_PA$PA, Sp_File) %>%
    terra::wrap() %>%
    list()

  Binary_R_Masked_Out <- stats::setNames(Sp_PA$PA_Masked, Sp_File) %>%
    terra::wrap() %>%
    list()

  Binary_R_Mask_Keep <- stats::setNames(Sp_PA$Mask_Keep, Sp_File) %>%
    terra::wrap() %>%
    list()

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

  dplyr::select(Results, -GBIF_R, -EASIN_R, -eLTER_R) %>%
    IASDT.R::save_as(
      object = ., object_name = paste0(Sp_File, "_Summary"),
      out_path = IASDT.R::path(
        Path_PA_Summary, paste0(Sp_File, "_Summary.RData")))

  # # ..................................................................... ###

  # Plotting species distribution -----
  IASDT.R::cat_time("Plotting species distribution")

  IASDT.R::IAS_plot(
    species = species, env_file = env_file, overwrite = overwrite)
  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTime,
    prefix = "\nProcessing species data was finished in ", ... = "\n")

  return(dplyr::select(Results, -GBIF_Gr100, -EASIN_Gr100, -eLTER_Gr100))
}
