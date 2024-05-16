# |---------------------------------------------------| #
# PrepModData ----
# |---------------------------------------------------| #

#' Prepare habitat-specific data for the models
#'
#' Prepare habitat-specific data for the models
#' @param Hab_Abb Habitat type
#' @param MinPresGrids Minimum number of presence grid cells per species. Only species with >= this number will be considered
#' @param Path_EnvFile Path to read the environment variables
#' @param ReturnData return data object
#' @param OutputPath Dir path for the output file
#' @name PrepModData
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PrepModData <- function(
    Hab_Abb = NULL, MinPresGrids = 50, Path_EnvFile = ".env",
    ReturnData = FALSE, OutputPath = "Data/Models") {

  MissingArgs <- list(
    Path_EnvFile = Path_EnvFile, Hab_Abb = Hab_Abb, OutputPath = OutputPath) %>%
    purrr::map(~inherits(.x, "character") && nchar(.x) > 0) %>%
    purrr::discard(.p = isTRUE) %>%
    names() %>%
    sort()
  if (length(MissingArgs) > 0) {
    stop(
      paste0("The following argument(s) must be provided as an string with nchar > 0:\n >> ",
             paste0(MissingArgs, collapse = " | ")))
  }

  readRenviron(Path_EnvFile)

  Path_Grid <- Sys.getenv("DP_R_Mod_Path_Grid")
  Path_Bound <- Sys.getenv("DP_R_Mod_Path_Bound")
  Path_PA <- Sys.getenv("DP_R_Mod_Path_PA")
  Path_CLC_Summ <- Sys.getenv("DP_R_Mod_Path_CLC_Summ")
  Path_Chelsa_Time_CC <- Sys.getenv("DP_R_Mod_Path_Chelsa_Time_CC")
  Path_Roads <- Sys.getenv("DP_R_Mod_Path_Roads")
  Path_Rail <- Sys.getenv("DP_R_Mod_Path_Rail")
  Path_Bias <- Sys.getenv("DP_R_Mod_Path_Bias")

  MissingPaths <- list(
    Path_Grid = Path_Grid, Path_Bound = Path_Bound, Path_PA = Path_PA,
    Path_CLC_Summ = Path_CLC_Summ, Path_Chelsa_Time_CC = Path_Chelsa_Time_CC,
    Path_Roads = Path_Roads, Path_Rail = Path_Rail, Path_Bias = Path_Bias) %>%
    purrr::map(~inherits(.x, "character") && nchar(.x) > 0) %>%
    purrr::discard(.p = isTRUE) %>%
    names() %>%
    sort()
  if (length(MissingPaths) > 0) {
    stop(
      paste0("The following paths must be provided in the .env file > 0:\n >> ",
             paste0(MissingPaths, collapse = " | ")))
  }

  ValidHabAbbs <- c(1:3, "4a", "4b", 5, 6, 8, 10, "12a", "12b")
  if (magrittr::not(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    stop(
      paste0("Hab_Abb has to be one of the following:\n >> ",
             paste0(ValidHabAbbs, collapse = " | ")))
  }

  fs::dir_create(OutputPath)

  Hab_column <- c(
    "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
    "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
    "Hab_5_Sandy", "Hab_6_Rocky", "Hab_8_Saline", "Hab_10_Wetland",
    "Hab_12a_Ruderal_habitats", "Hab_12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("_", as.character(Hab_Abb), "_"))

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # TaxaList <- IASDT.R::LoadAs(Path_TaxaList)

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  R_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::filter(!!as.symbol(Hab_column), NCells >= MinPresGrids) %>%
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

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## CHELSA -----

  BioVars <- c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18")
  R_Chelsa <- Path_Chelsa_Time_CC %>%
    file.path("St_1981_2010.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    stats::setNames(stringr::str_remove(names(.), "_1981_2010")) %>%
    terra::subset(BioVars)

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## SynHab -----

  R_Hab <- file.path(Path_CLC_Summ, "PercCov_SynHab_Crop.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
    stats::setNames("Hab")

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road ----

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

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Rail ----

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

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road + rail ----

  R_RoadRail <- (R_RoadInt + R_RailInt) %>%
    stats::setNames("RoadRail")
  R_RoadRailLog <- log10(R_RoadRail + 0.1) %>%
    stats::setNames("RoadRailLog")

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Sampling intensity ----

  R_Bias <- file.path(Path_Bias, "Bias_GBIF_SummaryR.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap() %>%
    magrittr::extract2("NObs") %>%      # NObs_Native
    stats::setNames("Bias")
  R_BiasLog <- log10(R_Bias + 0.1) %>%
    stats::setNames("BiasLog")

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Reference grid -----
  EU_Grid <- Path_Grid %>%
    file.path("Grid_10_sf.RData") %>%
    IASDT.R::LoadAs() %>%
    magrittr::extract2("Grid_10_sf_s")

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Country boundary -----
  EU_Bound <- Path_Bound %>%
    file.path("Bound_sf_Eur.RData") %>%
    IASDT.R::LoadAs() %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_01") %>%
    dplyr::select(NAME_ENGL)

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Merging data together -----

  ### Combine maps and convert to tibble -----
  DT_All <- c(
    R_Chelsa, R_Hab, R_RoadInt, R_RoadIntLog, R_RoadDist, R_RailInt,
    R_RailIntLog, R_RailDist, R_RoadRail, R_RoadRailLog, R_Bias, R_BiasLog, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble()

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # get country and grid cell ID
  DT_Country <- DT_All %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    dplyr::select(cell, geometry) %>%
    sf::st_join(EU_Grid) %>%
    sf::st_join(EU_Bound) %>%
    dplyr::rename(Country = NAME_ENGL) %>%
    sf::st_drop_geometry()

  ## # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  DT_All <- dplyr::full_join(DT_Country, DT_All, by = "cell")

  invisible(gc())

  OutObjName <- paste0("ModelDT_", stringr::str_remove(Hab_column, "Hab_"))
  OutputPath2 <- file.path(OutputPath, paste0(OutObjName, ".RData"))

  IASDT.R::SaveAs(InObj = DT_All, OutObj = OutObjName, OutPath = OutputPath2)

  # cat(paste0(
  #   "Data for ", crayon::bold(Hab_column), " was saved to: ", OutputPath2))

  if (ReturnData) {
    return(DT_All)
  } else {
    return(invisible(NULL))
  }
}
