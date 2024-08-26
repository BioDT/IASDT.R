## |------------------------------------------------------------------------| #
# Mod_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare habitat-specific data for Hmsc models
#'
#' This function prepares habitat-specific data for Hmsc models by processing
#' environmental and species presence data. It checks input arguments, reads
#' environment variables from a file, verifies paths, loads and filters species
#' data based on habitat type and minimum presence grid cells per species, and
#' merges various environmental layers (e.g., CHELSA Bioclimatic variables,
#' habitat coverage, road and railway intensity, sampling intensity) into a
#' single dataset. The processed data can be saved to disk as `*.RData` file .
#' @param Hab_Abb Character. Habitat abbreviation indicating the specific
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type to
#'   prepare data for. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`,
#'   `12a`, `12b`. If `Hab_Abb` = `0`, data is prepared irrespective of the
#'   habitat type. For more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param MinEffortsSp Minimum number of vascular plant species per grid cell in
#'   GBIF database for a grid cell to be included in the models. This to exclude
#'   grid cells with very little sampling efforts. Defaults to `100`.
#' @param NVars Integer. The number of variables used in the model. This
#'   argument has to be provided and can not be set to `NULL` (default).
#' @param PresPerVar Integer. Number of presence grid cells per predictor
#'   required for a species to be included in the analysis. The number of
#'   presence grid cells per species is calculated as the product of this factor
#'   and the number of variables used in the models `NVars` and is calculated
#'   after discarding grid cells with low sampling efforts (`MinEffortsSp`).
#'   Defaults to `10`.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param BioVars Character vector. Specifies the Bioclimatic variables to be
#'   included from the CHELSA dataset. Defaults to a selection of variables:
#'   `bio4`, `bio6`, `bio8`, `bio12`, `bio15`, and `bio18`.
#' @param Path_Model Character. Path where the output file should be saved.
#' @param VerboseProgress Logical. Indicates whether progress messages should be
#'   displayed. Defaults to `FALSE`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param SaveData Logical. Indicates whether the processed data should be saved
#'   as RData file. Defaults to `FALSE`.
#' @name Mod_PrepData
#' @author Ahmed El-Gabbas
#' @return a tibble containing modelling data.
#' @details The function reads the following environment variables:
#'    - **`DP_R_Grid`** (if `FromHPC` = `TRUE`) or **`DP_R_Grid_Local`** (if
#'   `FromHPC` = `FALSE`). The function reads the content of the
#'   `Grid_10_Land_Crop.RData` and `Grid_10_Land_Crop_sf_Country.RData` files
#'    - **`DP_R_Grid_Ref`** or **`DP_R_Grid_Ref_Local`**: The function reads the
#'   content of `Grid_10_sf.RData` file from this path.
#'    - **`DP_R_PA`** or **`DP_R_PA_Local`**: The function reads the contents of the
#'   `Sp_PA_Summary_DF.RData` file from this path.
#'    - **`DP_R_CLC_Summary`** / **`DP_R_CLC_Summary_Local`**: Path containing
#'   the `PercCov_SynHab_Crop.RData` file. This file contains maps for the
#'   percentage coverage of each SynHab habitat type per grid cell.
#'    - **`DP_R_CHELSA_Time_CC`** / **`DP_R_CHELSA_Time_CC_Local`**: Path
#'   containing the `St_1981_2010.RData` file. This file contains processed
#'   `CHELSA` data for the current climate.
#'    - **`DP_R_Roads`** / **`DP_R_Roads_Local`**: Path for processed road data.
#'   The function reads the contents of: `Road_Length.RData` for the total
#'   length of any road type per grid cell and `Dist2Road.RData` for the minimum
#'   distance between the centroid of each grid cell and the nearest road.
#'    - **`DP_R_Railway`** / **`DP_R_Railway_Local`**: Path for processed
#'   railway data. The function reads the contents of: `Railway_Length.RData`
#'   for the total length of any railway type per grid cell and
#'   `Dist2Rail.RData` for the minimum distance between the centroid of each
#'   grid cell and the nearest railway.
#'    - **`DP_R_Efforts`** / **`DP_R_Efforts_Local`**: Path for processed sampling
#'   efforts analysis. The function reads the content of
#'   `Bias_GBIF_SummaryR.RData` file containing the total number of GBIF
#'   vascular plant observations per grid cell.
#' @export

Mod_PrepData <- function(
    Hab_Abb = NULL, MinEffortsSp = 100L, NVars = NULL, PresPerVar = 10L,
    EnvFile = ".env",
    BioVars = c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18"),
    Path_Model = NULL, VerboseProgress = FALSE, FromHPC = TRUE,
    SaveData = FALSE) {

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

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # Check input parameters ----
  # # |||||||||||||||||||||||||||||||||||

  CheckNULL <- c("Hab_Abb", "Path_Model", "EnvFile", "NVars")
  IsNull <- purrr::map_lgl(CheckNULL, ~is.null(get(.x)))
  if (any(IsNull)) {
    stop(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"), call. = FALSE)
  }

  Hab_Abb <- as.character(Hab_Abb)

  if (!VerboseProgress) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpeciesID <- Species_name <- Species_File <- PA <-
    cell <- x <- Path_PA <- Path_Grid <- Path_Grid_Ref <- Path_CLC_Summ <-
    Path_Roads <- Path_Rail <- Path_Bias <- Path_Chelsa_Time_CC <-
    NGrids <- NSp <- EU_Bound <- NCells <- SpPA <- NPres <- Grid_R <- NULL

  IASDT.R::CatTime("Checking input arguments")
  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c("EnvFile", "Hab_Abb", "Path_Model")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("MinEffortsSp", "PresPerVar", "NVars"),
    Type = "numeric")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # Reading/checking environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  # Input data paths - these are read from the .env file

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_PA", "DP_R_PA", TRUE, FALSE,
      "Path_CLC_Summ", "DP_R_CLC_Summary", TRUE, FALSE,
      "Path_Chelsa_Time_CC", "DP_R_CHELSA_Time_CC", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads", TRUE, FALSE,
      "Path_Rail", "DP_R_Railway", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE,
      "Path_CLC_Summ", "DP_R_CLC_Summary_Local", TRUE, FALSE,
      "Path_Chelsa_Time_CC", "DP_R_CHELSA_Time_CC_Local", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads_Local", TRUE, FALSE,
      "Path_Rail", "DP_R_Railway_Local", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts_Local", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # Checking arguments ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking arguments")

  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")

  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0("Hab_Abb has to be one of the following:\n >> ",
             paste0(ValidHabAbbs, collapse = ", ")), call. = FALSE)
  }

  fs::dir_create(Path_Model)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # Loading data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading data")

  ## Sampling intensity / efforts ----
  IASDT.R::CatTime("Sampling intensity/efforts", Level = 1)
  R_Bias <- file.path(Path_Bias, "Bias_GBIF_SummaryR.RData")
  if (!file.exists(R_Bias)) {
    stop(paste0(R_Bias, " file does not exist"), call. = FALSE)
  }

  EffortsMask <- IASDT.R::LoadAs(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NSp") %>%
    terra::classify(
      rcl = matrix(
        c(0, MinEffortsSp, NA, MinEffortsSp, Inf, 1),
        byrow = TRUE, ncol = 3),
      include.lowest = TRUE, right = FALSE)

  ## Species data summary ----

  IASDT.R::CatTime("Species data summary", Level = 1)

  DT_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")

  if (!file.exists(DT_Sp)) {
    stop(paste0(DT_Sp, " file does not exist"), call. = FALSE)
  }

  DT_Sp <- IASDT.R::LoadAs(DT_Sp)

  if (Hab_Abb == "0") {
    Hab_column <- NULL
  } else {
    Hab_column <- c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("_", as.character(Hab_Abb), "_"))

    DT_Sp <- dplyr::filter(DT_Sp, !!as.symbol(Hab_column))
  }

  ## Species Presence-absence data ----

  IASDT.R::CatTime("Species Presence-absence data", Level = 1)

  # minimum number of presence grids per species
  MinPresGrids <- PresPerVar * NVars

  R_Sp <- DT_Sp %>%
    # Exclude species with too few presence grid cells. There will be further
    # exclusion of species with few grid cells in this pipeline, but excluding
    # this first may help to reduce processing time
    dplyr::filter(NCells >= MinPresGrids) %>%
    dplyr::select(SpeciesID, Species_name, Species_File) %>%
    # Mask each species map with the filtered grid cells
    dplyr::mutate(
      PA = purrr::map2(
        .x = Species_File, .y = SpeciesID,
        .f = ~{
          # Masked raster map
          SpPA <- file.path(
            Path_PA, "RData", stringr::str_c(.x, "_PA.RData")) %>%
            IASDT.R::LoadAs() %>%
            terra::unwrap() %>%
            magrittr::extract2("PA") %>%
            terra::mask(EffortsMask) %>%
            stats::setNames(paste0("Sp_", .y))
          # Number of presence grid cells after masking
          NPres <- as.integer(terra::global(SpPA, "sum", na.rm = TRUE))
          return(tibble::tibble(SpPA = list(SpPA), NPres = NPres))
        })) %>%
    tidyr::unnest_wider(PA) %>%
    dplyr::mutate(SpPA = unlist(SpPA)) %>%
    # filter species with too few observations (after masking)
    dplyr::filter(NPres >= MinPresGrids) %>%
    dplyr::pull(SpPA) %>%
    terra::rast()

  # Change the efforts mask to exclude grid cells with no species
  ZeroSpGrids <- (sum(R_Sp) == 0)
  EffortsMask[ZeroSpGrids] <- NA
  R_Sp[ZeroSpGrids] <- NA

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plotting number of grid cells per species", Level = 2)

  EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_03")
  R_Sp_sum <- sum(R_Sp, na.rm = TRUE)
  R_Sp_sumP <- terra::classify(R_Sp_sum, cbind(0, NA))
  NGridsWzSpecies <- terra::classify(R_Sp_sumP, cbind(1, Inf, 1)) %>%
    terra::global("sum", na.rm = TRUE) %>%
    as.integer() %>%
    format(big.mark = ",")
  Limits <- terra::trim(R_Sp_sumP) %>%
    terra::ext() %>%
    as.vector()
  NSpPerGrid <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = R_Sp_sumP) +
    tidyterra::scale_fill_whitebox_c(
      na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
    ggplot2::geom_sf(
      data = EU_Bound, fill = "transparent", colour = "black",
      linewidth = 0.15) +
    ggplot2::labs(
      title = paste0(
        '<span style="color:blue; font-size:20px;"><b>',
        "Number of IAS per grid cell to be used in the models</b></span>",
        '<span style="color:black; font-size:16px;"> (',
        stringr::str_remove(Hab_column, "Hab_"), ")</span>"),
      caption  = paste0(
        "Only grid cells with &#8805;", MinEffortsSp,
        " vascular plant species in GBIF (", NGridsWzSpecies,
        ") and species with ",
        "&#8805;", MinPresGrids, " presence grid cells (", terra::nlyr(R_Sp),
        ") are considered")) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = Limits[3:4]) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = Limits[1:2]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.25, 0, 0.25, 0, "cm"),
      plot.title = ggtext::element_markdown(
        size = 16, color = "blue",
        margin = ggplot2::margin(0, 0, 0.1, 0, "cm")),
      plot.caption = ggtext::element_markdown(
        size = 12, colour = "grey40", hjust = 0.3),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.9),
      legend.key.size = grid::unit(0.8, "cm"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank())

  ggplot2::ggsave(
    plot = NSpPerGrid, width = 25, height = 27, units = "cm", dpi = 600,
    filename = file.path(Path_Model, "NSpPerGrid.jpeg"))

  rm(Limits, NSpPerGrid, R_Sp_sum, R_Sp_sumP, EU_Bound)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## CHELSA -----

  IASDT.R::CatTime("CHELSA", Level = 1)
  R_Chelsa <- file.path(Path_Chelsa_Time_CC, "St_1981_2010.RData")
  if (!file.exists(R_Chelsa)) {
    stop(paste0(R_Chelsa, " file does not exist"), call. = FALSE)
  }
  R_Chelsa <- IASDT.R::LoadAs(R_Chelsa) %>%
    terra::unwrap() %>%
    stats::setNames(stringr::str_remove(names(.), "_1981_2010")) %>%
    terra::subset(BioVars) %>%
    terra::mask(EffortsMask)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Reference grid -----

  IASDT.R::CatTime("Reference grid - sf", Level = 1)
  # Reference grid as sf
  Grid_SF <- file.path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(Grid_SF)) {
    stop(paste0(Grid_SF, " file does not exist"), call. = FALSE)
  }
  Grid_SF <- IASDT.R::LoadAs(Grid_SF) %>%
    magrittr::extract2("Grid_10_sf_s")


  # Reference grid as sf - country names
  IASDT.R::CatTime("Reference grid - country names", Level = 1)
  Grid_CNT <- file.path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData")
  if (!file.exists(Grid_CNT)) {
    stop(paste0(Grid_CNT, " file does not exist"), call. = FALSE)
  }
  Grid_CNT <- IASDT.R::LoadAs(Grid_CNT) %>%
    dplyr::mutate(x = sf::st_coordinates(.)[, 1],
                  y = sf::st_coordinates(.)[, 2]) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  IASDT.R::CatTime("Habitat coverage", Level = 1)

  Path_Hab <- file.path(Path_CLC_Summ, "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    stop(
      paste0("Path_Hab file: ", Path_Hab, " does not exist"), call. = FALSE)
  }

  if (Hab_Abb == "0") {
    # Use dummy habitat values
    R_Hab <- stats::setNames(Grid_R, "Hab")
    R_HabLog <- stats::setNames(Grid_R, "HabLog")
  } else {
    R_Hab <- IASDT.R::LoadAs(Path_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
      terra::mask(EffortsMask) %>%
      stats::setNames("Hab")

    R_HabLog <- log10(R_Hab + 0.1) %>%
      stats::setNames("HabLog")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road ----
  IASDT.R::CatTime("Road intensity", Level = 1)

  # road intensity of any road type
  R_RoadInt <- file.path(Path_Roads, "Road_Length.RData")
  if (!file.exists(R_RoadInt)) {
    stop(paste0(R_RoadInt, " file does not exist"), call. = FALSE)
  }
  R_RoadInt <- IASDT.R::LoadAs(R_RoadInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt") %>%
    terra::mask(EffortsMask)

  # Log of road intensity
  R_RoadIntLog <- log10(R_RoadInt + 0.1) %>%
    stats::setNames("RoadIntLog")

  # Distance to roads
  R_RoadDist <- file.path(Path_Roads, "Dist2Road.RData")
  if (!file.exists(R_RoadDist)) {
    stop(paste0(R_RoadDist, " file does not exist"), call. = FALSE)
  }
  R_RoadDist <- IASDT.R::LoadAs(R_RoadDist) %>%
    terra::unwrap() %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("RoadDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Rail ----

  IASDT.R::CatTime("Railway intensity", Level = 1)

  # Railway intensity
  R_RailInt <- file.path(Path_Rail, "Railways_Length.RData")
  if (!file.exists(R_RailInt)) {
    stop(paste0(R_RailInt, " file does not exist"), call. = FALSE)
  }
  R_RailInt <- IASDT.R::LoadAs(R_RailInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("rail") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("RailInt")

  # Log of railway intensity
  R_RailIntLog <- log10(R_RailInt + 0.1) %>%
    stats::setNames("RailIntLog")

  # Distance to nearest rail
  R_RailDist <- file.path(Path_Rail, "Dist2Rail.RData")
  if (!file.exists(R_RailDist)) {
    stop(paste0(R_RailDist, " file does not exist"), call. = FALSE)
  }
  R_RailDist <- IASDT.R::LoadAs(R_RailDist) %>%
    terra::unwrap() %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("RailDist")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Road + rail ----

  IASDT.R::CatTime("Railway + road intensity", Level = 1)
  R_RoadRail <- (R_RoadInt + R_RailInt) %>%
    stats::setNames("RoadRail")
  R_RoadRailLog <- log10(R_RoadRail + 0.1) %>%
    stats::setNames("RoadRailLog")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Merging data together -----
  IASDT.R::CatTime("Merging data together")

  ColumnsFirst <- c("CellNum", "CellCode", "Country", "Country_Nearest")
  DT_All <- c(
    R_Chelsa, R_Hab, R_HabLog, R_RoadInt, R_RoadIntLog, R_RoadDist,
    R_RailInt, R_RailIntLog, R_RailDist, R_RoadRail, R_RoadRailLog, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble() %>%
    # Add country name
    dplyr::left_join(Grid_CNT, by = c("x", "y")) %>%
    dplyr::rename(Country_Nearest = "Nearest", CellNum = cell) %>%
    dplyr::select(tidyselect::all_of(ColumnsFirst), tidyselect::everything())

  if (Hab_Abb == "0") {
    DT_All$Hab <- NA_real_
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model data to disk -----

  if (SaveData) {
    IASDT.R::CatTime("Save model data to disk")
    if (Hab_Abb == "0") {
      OutObjName <- "ModDT_0_All"
    } else {
      OutObjName <- paste0("ModDT_", stringr::str_remove(Hab_column, "Hab_"))
    }
    IASDT.R::SaveAs(
      InObj = DT_All, OutObj = OutObjName,
      OutPath = file.path(Path_Model, paste0(OutObjName, ".RData")))
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  return(DT_All)
}
