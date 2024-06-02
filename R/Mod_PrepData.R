## |------------------------------------------------------------------------| #
# Mod_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare habitat-specific data for the models
#'
#' Prepare habitat-specific data for the models
#'
#' @param Hab_Abb String. Habitat type. This has to be one of the following: c("0", "1", "2", "3", "4a", "4b", "5", "6", "8", "10", "12a", "12b"). "0" means prepare data irrespective of the habitat type
#' @param MinPresGrids Integer. Minimum number of presence grid cells per species. Only species with ≥ this number will be considered
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param BioVars String. The bioclimatic variables to get from CHELSA. Default value: `c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18")`
#' @param ReturnData Logical. Should the resulted data be returned as an R object? Default: `FALSE`
#' @param OutputPath String. Path to save the output file.
#' @param VerboseProgress Logical. Show messages for the progress of creating files
#' @name Mod_PrepData
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_PrepData <- function(
    Hab_Abb = NULL, MinPresGrids = 50, EnvFile = ".env",
    BioVars = c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18"),
    ReturnData = FALSE, OutputPath = NULL, VerboseProgress = FALSE) {


  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- SpeciesID <- Species_name <- Species_File <- PA <-
    NAME_ENGL <- cell <- geometry <- Country <- Country2 <-
    x <- y <- NULL


  if (magrittr::not(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  Hab_Abb <- as.character(Hab_Abb)

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c("EnvFile", "Hab_Abb", "OutputPath")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Input data paths - these are read from the .env file

  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_Grid <- Sys.getenv("DP_R_Mod_Path_Grid")
    Path_Bound <- Sys.getenv("DP_R_Mod_Path_Bound")
    Path_PA <- Sys.getenv("DP_R_Mod_Path_PA")
    Path_CLC_Summ <- Sys.getenv("DP_R_Mod_Path_CLC_Summ")
    Path_Chelsa_Time_CC <- Sys.getenv("DP_R_Mod_Path_Chelsa_Time_CC")
    Path_Roads <- Sys.getenv("DP_R_Mod_Path_Roads")
    Path_Rail <- Sys.getenv("DP_R_Mod_Path_Rail")
    Path_Bias <- Sys.getenv("DP_R_Mod_Path_Bias")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 ## Paths checking ----

  IASDT.R::CatTime("Checking paths")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c(
    "Path_Grid", "Path_Bound", "Path_PA", "Path_CLC_Summ",
    "Path_Chelsa_Time_CC", "Path_Roads", "Path_Rail", "Path_Bias")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("ReturnData", "VerboseProgress"),
    Type = "logical")

  Path_List <- list(
    Path_Grid = Path_Grid, Path_Bound = Path_Bound, Path_PA = Path_PA,
    Path_CLC_Summ = Path_CLC_Summ, Path_Chelsa_Time_CC = Path_Chelsa_Time_CC,
    Path_Roads = Path_Roads, Path_Rail = Path_Rail, Path_Bias = Path_Bias)

  MissingPaths <- Path_List %>%
    purrr::map(fs::dir_exists) %>%
    purrr::discard(.p = isTRUE) %>%
    names() %>%
    sort()
  if (length(MissingPaths) > 0) {
    stop(
      paste0("The following path(s) does not exist.",
             " Please check provided values in the .env file\n >> ",
             paste0(MissingPaths, collapse = " | ")))
  }
  rm(CharArgs)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking arguments")
  ValidHabAbbs <- c(0:3, "4a", "4b", 5, 6, 8, 10, "12a", "12b")
  if (magrittr::not(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    stop(
      paste0("Hab_Abb has to be one of the following:\n >> ",
             paste0(ValidHabAbbs, collapse = " | ")))
  }
  Hab_Abb <- as.character(Hab_Abb)

  fs::dir_create(OutputPath)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load data")

  IASDT.R::CatTime(">> Load species data summary")
  R_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData") %>%
    IASDT.R::LoadAs()

  if (Hab_Abb == "0") {
    Hab_column <- NULL
  } else {
    Hab_column <- c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_5_Sandy", "Hab_6_Rocky", "Hab_8_Saline", "Hab_10_Wetland",
      "Hab_12a_Ruderal_habitats", "Hab_12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("_", as.character(Hab_Abb), "_"))

    R_Sp <- dplyr::filter(R_Sp, !!as.symbol(Hab_column))
  }

  IASDT.R::CatTime(">> Load species PA data")
  R_Sp <- R_Sp %>%
    dplyr::filter(NCells >= MinPresGrids) %>%
    dplyr::select(SpeciesID, Species_name, Species_File) %>%
    dplyr::mutate(
      PA = purrr::map2(
        .x = Species_File, .y = SpeciesID,
        .f = ~{
          R <- file.path(Path_PA, "RData", stringr::str_c(.x, "_PA.RData")) %>%
            IASDT.R::LoadAs() %>%
            terra::unwrap() %>%
            magrittr::extract2("PA") %>%
            stats::setNames(paste0("Sp_", .y))
        })) %>%
    dplyr::pull(PA) %>%
    terra::rast()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## CHELSA -----

  IASDT.R::CatTime(">> Load CHELSA data")
  R_Chelsa <- Path_Chelsa_Time_CC %>%
    file.path("St_1981_2010.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    stats::setNames(stringr::str_remove(names(.), "_1981_2010")) %>%
    terra::subset(BioVars)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  IASDT.R::CatTime(">> Load Habitat coverage")
  if (Hab_Abb == "0") {
    R_Hab <- file.path(Path_CLC_Summ, "PercCov_SynHab_Crop.RData") %>%
      IASDT.R::LoadAs() %>%
      terra::unwrap() %>%
      magrittr::extract2(1) %>%
      stats::setNames("Hab")
    terra::set.values(
      x = R_Hab, cells = seq_len(terra::ncell(R_Hab)), values = 1)
  } else {
    R_Hab <- file.path(Path_CLC_Summ, "PercCov_SynHab_Crop.RData") %>%
      IASDT.R::LoadAs() %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
      stats::setNames("Hab")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road ----
  IASDT.R::CatTime(">> Load Road intensity")
  R_RoadInt <- file.path(Path_Roads, "Road_Length.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt")
  R_RoadIntLog <- R_RoadInt %>%
    magrittr::add(0.1) %>%
    log10() %>%
    stats::setNames("RoadIntLog")
  R_RoadDist <- file.path(Path_Roads, "Dist2Road.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    stats::setNames("RoadDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Rail ----

  IASDT.R::CatTime(">> Load railway intensity")
  R_RailInt <- file.path(Path_Rail, "Railway_Length.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    magrittr::extract2("rail") %>%
    stats::setNames("RailInt")
  R_RailIntLog <- R_RailInt %>%
    magrittr::add(0.1) %>%
    log10() %>%
    stats::setNames("RailIntLog")
  R_RailDist <- file.path(Path_Rail, "Dist2Rail.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    stats::setNames("RailDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road + rail ----

  IASDT.R::CatTime(">> Load railway + road intensity")
  R_RoadRail <- (R_RoadInt + R_RailInt) %>%
    stats::setNames("RoadRail")
  R_RoadRailLog <- log10(R_RoadRail + 0.1) %>%
    stats::setNames("RoadRailLog")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Sampling intensity ----
  IASDT.R::CatTime(">> Load sampling intensity")
  R_Bias <- file.path(Path_Bias, "Bias_GBIF_SummaryR.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    # NObs or NObs_Native
    magrittr::extract2("NObs") %>%
    stats::setNames("Bias")
  R_BiasLog <- log10(R_Bias + 0.1) %>%
    stats::setNames("BiasLog")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Reference grid -----

  IASDT.R::CatTime(">> Load reference grid")
  EU_Grid <- Path_Grid %>%
    file.path("Grid_10_sf.RData") %>%
    IASDT.R::LoadAs() %>%
    magrittr::extract2("Grid_10_sf_s")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Country boundary -----

  IASDT.R::CatTime(">> Load country boundaries")
  EU_Bound <- Path_Bound %>%
    file.path("Bound_sf_Eur.RData") %>%
    IASDT.R::LoadAs() %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_01") %>%
    dplyr::select(NAME_ENGL)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Merging data together -----
  IASDT.R::CatTime("Merging data together")

  ### Combine maps and convert to tibble -----
  IASDT.R::CatTime(">> Combine maps and convert to tibble")
  DT_All <- c(
    R_Chelsa, R_Hab, R_RoadInt, R_RoadIntLog, R_RoadDist,
    R_RailInt, R_RailIntLog, R_RailDist, R_RoadRail, R_RoadRailLog,
    R_Bias, R_BiasLog, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # get country and grid cell ID
  IASDT.R::CatTime(">> get country and grid cell ID")

  # List of countries in the boundaries shapefile
  CountryNames <- dplyr::pull(EU_Bound, NAME_ENGL)

  # Find spatially matching countries
  IASDT.R::CatTime(">> find spatially matching countries")
  DT_Country <- DT_All %>%
    dplyr::select(cell, x, y) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  IASDT.R::CatTime(">> Add grid cell code")
  DT_Country <- sf::st_join(DT_Country, EU_Grid)

  IASDT.R::CatTime(">> Add country name")
  DT_Country <- sf::st_join(DT_Country, EU_Bound) %>%
    dplyr::rename(Country = NAME_ENGL)

  # find nearest countries for unmatched grid cells
  IASDT.R::CatTime(">> find nearest countries for unmatched grid cells")
  MissingCountries <- DT_Country %>%
    dplyr::filter(is.na(Country)) %>%
    dplyr::mutate(
      Country2 = CountryNames[sf::st_nearest_feature(geometry, EU_Bound)]) %>%
    dplyr::select(cell, Country2) %>%
    sf::st_drop_geometry()

  DT_Country <- DT_Country %>%
    sf::st_drop_geometry() %>%
    dplyr::left_join(MissingCountries, by = "cell") %>%
    dplyr::mutate(Country = dplyr::coalesce(Country, Country2)) %>%
    dplyr::select(-Country2)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  DT_All <- dplyr::full_join(DT_Country, DT_All, by = "cell")

  invisible(gc())

  if (Hab_Abb == "0") {
    DT_All$Hab <- NA_real_
    OutObjName <- paste0("ModelData_", MinPresGrids, "Grids_0_All")
  } else {
    OutObjName <- paste0("ModelData_", MinPresGrids, "Grids_",
                         stringr::str_remove(Hab_column, "Hab_"))
  }

  IASDT.R::SaveAs(
    InObj = DT_All, OutObj = OutObjName,
    OutPath = file.path(OutputPath, paste0(OutObjName, ".RData")))

  if (ReturnData) {
    return(DT_All)
  } else {
    return(invisible(NULL))
  }
}
