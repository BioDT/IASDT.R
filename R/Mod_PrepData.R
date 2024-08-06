## |------------------------------------------------------------------------| #
# Mod_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare habitat-specific data for Hmsc models
#'
#' This function prepares habitat-specific data for Hmsc models by processing
#' environmental and species presence data. It checks input arguments, reads
#' environment variables from a file, verifies paths, loads and filters species
#' data based on habitat type and minimum presence grid cells per species, and
#' merges various environmental layers (e.g., CHELSA bioclimatic variables,
#' habitat coverage, road and railway intensity, sampling intensity) into a
#' single dataset. The processed data can be saved to disk and/or returned as an
#' R object.
#' @param Hab_Abb Character. Habitat abbreviation indicating the specific
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type to
#'   prepare data for. Valid values include "0", "1", "2", "3", "4a", "4b",
#'   "10", "12a", "12b". If `Hab_Abb` = "0", data is prepared irrespective of
#'   the habitat type.
#' @param MinPresGrids Integer. Minimum number of presence grid cells required
#'   for a species to be included in the analysis. Defaults to 50.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param BioVars Character vector. Specifies the bioclimatic variables to be
#'   included from the CHELSA dataset. Defaults to a selection of variables:
#'   "bio4", "bio6", "bio8", "bio12", "bio15", and "bio18".
#' @param ReturnData Logical. Indicates whether the processed data should be
#'   returned as an R object. Defaults to `FALSE`, meaning the data will not be
#'   returned but only saved to disk.
#' @param OutputPath Character. Path where the output file should be saved.
#' @param VerboseProgress Logical. Indicates whether progress messages should be
#'   displayed. Defaults to `FALSE`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @name Mod_PrepData
#' @author Ahmed El-Gabbas
#' @return If `ReturnData` is `TRUE`, returns a data frame containing the
#'   modelling data. Otherwise, invisibly returns `NULL`.
#' @details The function reads the following environment variables:
#'    - **`DP_R_Grid`** (if `FromHPC` = `TRUE`) or **`DP_R_Grid_Local`** (if
#'   `FromHPC` = `FALSE`). The function reads the content of the following files
#'   in this path: `Grid_10_sf.RData`, `Grid_10_Land_Crop.RData`, and
#'   `Grid_10_Land_Crop_sf_Country.RData`
#'    - **`DP_R_PA`** (if `FromHPC` = `TRUE`) or **`DP_R_PA_Local`** (if
#'    `FromHPC` = `FALSE`). The function reads the content of the
#'     `Sp_PA_Summary_DF.RData` file from this path.
#'    - **`DP_R_CLC_Summary`** / **`DP_R_CLC_Summary_Local`**: Path containing
#'    the `PercCov_SynHab_Crop.RData` file. This file contains maps for the
#'   percentage coverage of each SynHab habitat type per grid cell.
#'    - **`DP_R_CHELSA_Time_CC`** / **`DP_R_CHELSA_Time_CC_Local`**: Path
#'    containing the `St_1981_2010.RData` file. This file contains processed
#'    `CHELSA` data for the current climate.
#'    - **`DP_R_Roads`** / **`DP_R_Roads_Local`**: Path for processed road data.
#'    The function reads the contents of: `Road_Length.RData` for the total
#'    length of any road type per grid cell and `Dist2Road.RData` for the
#'    minimum distance between the centroid of each grid cell and the nearest
#'    road.
#'    - **`DP_R_Railway`** / **`DP_R_Railway_Local`**: Path for processed
#'    railway data. The function reads the contents of: `Railway_Length.RData`
#'    for the total length of any railway type per grid cell and
#'    `Dist2Rail.RData` for the minimum distance between the centroid of each
#'    grid cell and the nearest railway.
#'    - **`DP_R_Bias`** / **`DP_R_Bias_Local`**: Path for processed sampling
#'    efforts analysis. The function reads the content of
#'    `Bias_GBIF_SummaryR.RData` file containing the total number of GBIF
#'    vascular plant observations per grid cell.
#' @export

Mod_PrepData <- function(
    Hab_Abb = NULL, MinPresGrids = 50, EnvFile = ".env",
    BioVars = c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18"),
    ReturnData = FALSE, OutputPath = NULL, VerboseProgress = FALSE,
    FromHPC = TRUE) {

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 2023
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # We decided to exclude the following SynHab habitat types from our analysis:
  # # 4. Grasslands: use 4a and 4b separately
  # # 7. Dryland
  # # 9. Riparian
  # # 11. Aquatic
  # # 12. Man-made: use 12a and 12b separately

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # July 2024
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # We decided to further exclude three SynHab habitat types
  # # 5. Sandy
  # # 6. Rocky
  # # 8. Saline

  if (is.null(Hab_Abb) || is.null(OutputPath) || is.null(EnvFile)) {
    stop("Hab_Abb, OutputPath, and EnvFile cannot be NULL")
  }
  Hab_Abb <- as.character(Hab_Abb)

  if (magrittr::not(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- SpeciesID <- Species_name <- Species_File <- PA <-
    cell <- x <- Path_PA <- Path_Grid <- Path_CLC_Summ <- Path_Roads <-
    Path_Rail <- Path_Bias <- Path_Chelsa_Time_CC <- NULL

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

  if (magrittr::not(file.exists(EnvFile))) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_PA", "DP_R_PA", TRUE, FALSE,
      "Path_CLC_Summ", "DP_R_CLC_Summary", TRUE, FALSE,
      "Path_Chelsa_Time_CC", "DP_R_CHELSA_Time_CC", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads", TRUE, FALSE,
      "Path_Rail", "DP_R_Railway", TRUE, FALSE,
      "Path_Bias", "DP_R_Bias", TRUE, FALSE)

  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE,
      "Path_CLC_Summ", "DP_R_CLC_Summary_Local", TRUE, FALSE,
      "Path_Chelsa_Time_CC", "DP_R_CHELSA_Time_CC_Local", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads_Local", TRUE, FALSE,
      "Path_Rail", "DP_R_Railway_Local", TRUE, FALSE,
      "Path_Bias", "DP_R_Bias_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking arguments")

  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")

  if (magrittr::not(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0("Hab_Abb has to be one of the following:\n >> ",
             paste0(ValidHabAbbs, collapse = ", ")))
  }

  fs::dir_create(OutputPath)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading data")

  IASDT.R::CatTime("   >>>   Species data summary")

  R_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")

  if (magrittr::not(file.exists(R_Sp))) {
    stop(paste0(R_Sp, " file does not exist"))
  }

  R_Sp <- IASDT.R::LoadAs(R_Sp)

  if (Hab_Abb == "0") {
    Hab_column <- NULL
  } else {
    Hab_column <- c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("_", as.character(Hab_Abb), "_"))

    R_Sp <- dplyr::filter(R_Sp, !!as.symbol(Hab_column))
  }

  IASDT.R::CatTime("   >>>   Species Presence-absence data")
  R_Sp <- dplyr::filter(R_Sp, NCells >= MinPresGrids) %>%
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

  IASDT.R::CatTime("   >>>   CHELSA data")
  R_Chelsa <- file.path(Path_Chelsa_Time_CC, "St_1981_2010.RData")
  if (magrittr::not(file.exists(R_Chelsa))) {
    stop(paste0(R_Chelsa, " file does not exist"))
  }
  R_Chelsa <- IASDT.R::LoadAs(R_Chelsa) %>%
    terra::unwrap() %>%
    stats::setNames(stringr::str_remove(names(.), "_1981_2010")) %>%
    terra::subset(BioVars)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Reference grid -----

  IASDT.R::CatTime("   >>>   Reference grid")

  IASDT.R::CatTime("   >>>   Reference grid   >>>   sf")
  # Reference grid as sf
  Grid_SF <- file.path(Path_Grid, "Grid_10_sf.RData")
  if (magrittr::not(file.exists(Grid_SF))) {
    stop(paste0(Grid_SF, " file does not exist"))
  }
  Grid_SF <- IASDT.R::LoadAs(Grid_SF) %>%
    magrittr::extract2("Grid_10_sf_s")

  # Reference grid as raster
  IASDT.R::CatTime("   >>>   Reference grid   >>>   raster")
  Grid_R <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (magrittr::not(file.exists(Grid_R))) {
    stop(paste0(Grid_R, " file does not exist"))
  }
  Grid_R <- IASDT.R::LoadAs(Grid_R) %>%
    terra::unwrap()

  # Reference grid as sf - country names
  IASDT.R::CatTime("   >>>   Reference grid   >>>   country names")
  Grid_CNT <- file.path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData")
  if (magrittr::not(file.exists(Grid_CNT))) {
    stop(paste0(Grid_CNT, " file does not exist"))
  }
  Grid_CNT <- IASDT.R::LoadAs(Grid_CNT)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  IASDT.R::CatTime("   >>>   Habitat coverage")

  Path_Hab <- file.path(Path_CLC_Summ, "PercCov_SynHab_Crop.RData")
  if (magrittr::not(file.exists(Path_Hab))) {
    stop("Path_Hab file: ", Path_Hab, " does not exist")
  }

  if (Hab_Abb == "0") {
    # Use dummy habitat values
    R_Hab <- stats::setNames(Grid_R, "Hab")
    R_HabLog <- stats::setNames(Grid_R, "HabLog")
  } else {
    R_Hab <- IASDT.R::LoadAs(Path_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
      terra::crop(Grid_R) %>%
      terra::mask(Grid_R) %>%
      stats::setNames("Hab")

    R_HabLog <- log10(R_Hab + 0.1) %>%
      stats::setNames("HabLog")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road ----
  IASDT.R::CatTime("   >>>   Road intensity")

  # road intensity of any road type
  R_RoadInt <- file.path(Path_Roads, "Road_Length.RData")
  if (magrittr::not(file.exists(R_RoadInt))) {
    stop(paste0(R_RoadInt, " file does not exist"))
  }
  R_RoadInt <- IASDT.R::LoadAs(R_RoadInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt")

  # Log of road intensity
  R_RoadIntLog <- log10(R_RoadInt + 0.1) %>%
    stats::setNames("RoadIntLog")

  # Distance to roads
  R_RoadDist <- file.path(Path_Roads, "Dist2Road.RData")
  if (magrittr::not(file.exists(R_RoadDist))) {
    stop(paste0(R_RoadDist, " file does not exist"))
  }
  R_RoadDist <- IASDT.R::LoadAs(R_RoadDist) %>%
    terra::unwrap() %>%
    stats::setNames("RoadDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Rail ----

  IASDT.R::CatTime("   >>>   Railway intensity")

  # Railway intensity
  R_RailInt <- file.path(Path_Rail, "Railway_Length.RData")
  if (magrittr::not(file.exists(R_RailInt))) {
    stop(paste0(R_RailInt, " file does not exist"))
  }
  R_RailInt <- IASDT.R::LoadAs(R_RailInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("rail") %>%
    stats::setNames("RailInt")

  # Log of railway intensity
  R_RailIntLog <- log10(R_RailInt + 0.1) %>%
    stats::setNames("RailIntLog")

  # Distance to nearest rail
  R_RailDist <- file.path(Path_Rail, "Dist2Rail.RData")
  if (magrittr::not(file.exists(R_RailDist))) {
    stop(paste0(R_RailDist, " file does not exist"))
  }
  R_RailDist <- IASDT.R::LoadAs(R_RailDist) %>%
    terra::unwrap() %>%
    stats::setNames("RailDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road + rail ----

  IASDT.R::CatTime("   >>>   Railway + road intensity")
  R_RoadRail <- (R_RoadInt + R_RailInt) %>%
    stats::setNames("RoadRail")
  R_RoadRailLog <- log10(R_RoadRail + 0.1) %>%
    stats::setNames("RoadRailLog")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Sampling intensity / efforts ----

  IASDT.R::CatTime("   >>>   Sampling intensity")
  R_Bias <- file.path(Path_Bias, "Bias_GBIF_SummaryR.RData")
  if (magrittr::not(file.exists(R_Bias))) {
    stop(paste0(R_Bias, " file does not exist"))
  }
  R_Bias <- IASDT.R::LoadAs(R_Bias) %>%
    terra::unwrap() %>%
    # NObs or NObs_Native
    magrittr::extract2("NObs") %>%
    stats::setNames("Bias")

  # log sampling efforts
  R_BiasLog <- log10(R_Bias + 0.1) %>%
    stats::setNames("BiasLog")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Merging data together -----
  IASDT.R::CatTime("Merging data together")

  ### Combine maps and convert to tibble -----
  IASDT.R::CatTime("   >>>   Combine maps and convert to tibble")
  DT_All <- c(
    R_Chelsa, R_Hab, R_HabLog, R_RoadInt, R_RoadIntLog, R_RoadDist,
    R_RailInt, R_RailIntLog, R_RailDist, R_RoadRail, R_RoadRailLog,
    R_Bias, R_BiasLog, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Add country name
  ColumnsFirst <- c("CellNum", "CellCode", "Country", "Country_Nearest")

  DT_All <- dplyr::select(DT_All, "x", "y") %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE) %>%
    sf::st_join(Grid_CNT) %>%
    sf::st_drop_geometry() %>%
    dplyr::right_join(DT_All, by = c("x", "y")) %>%
    dplyr::rename(Country_Nearest = "Nearest", CellNum = cell) %>%
    dplyr::select(tidyselect::all_of(ColumnsFirst), tidyselect::everything())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save model data to disk")

  if (Hab_Abb == "0") {
    DT_All$Hab <- NA_real_
    OutObjName <- paste0("ModDT_", MinPresGrids, "Grids_0_All")
  } else {
    OutObjName <- paste0(
      "ModDT_", MinPresGrids, "Grids_", stringr::str_remove(Hab_column, "Hab_"))
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
