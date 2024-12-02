## |------------------------------------------------------------------------| #
# Mod_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare habitat-specific data for Hmsc models
#'
#' This function processes environmental and species presence data to prepare
#' habitat-specific datasets for use in Hmsc models. It checks input arguments,
#' reads environment variables from a file, verifies paths, loads and filters
#' species data based on habitat type and minimum presence grid cells per
#' species, and merges various environmental layers (e.g., CHELSA Bioclimatic
#' variables, habitat coverage, road and railway intensity, sampling efforts)
#' into a single dataset. The processed data can be saved to disk as an
#' `*.RData` file.
#'
#' @param Hab_Abb Character. Abbreviation for the habitat type (based on
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548)) for which to prepare
#'   data. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`.
#'   If `Hab_Abb` = `0`, data is prepared irrespective of the habitat type. For
#'   more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param MinEffortsSp Integer specifying the minimum number of vascular plant
#'   species per grid cell (from GBIF data) required for inclusion in the
#'   models. This is to exclude grid cells with very little sampling efforts.
#'   Defaults to `100`.
#' @param PresPerSpecies Integer. The minimum number of presence grid cells for
#'   a species to be included in the analysis. The number of presence grid cells
#'   per species is calculated after discarding grid cells with low sampling
#'   efforts (`MinEffortsSp`). Defaults to `80`.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Path_Model Character. Path where the output file should be saved.
#' @param VerboseProgress Logical. Indicates whether progress messages should be
#'   displayed. Defaults to `TRUE`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param SaveData Logical. Indicates whether the processed data should be saved
#'   as RData file. Defaults to `TRUE`.
#' @param ExcludeCult Logical. Indicates whether to exclude countries with
#'   cultivated or casual observations per species. Defaults to `TRUE`.
#' @name Mod_PrepData
#' @author Ahmed El-Gabbas
#' @return a tibble containing modelling data.
#' @importFrom rlang .data
#' @details The function reads the following environment variables:
#'    - **`DP_R_Grid`** (if `FromHPC` = `TRUE`) or **`DP_R_Grid_Local`** (if
#'   `FromHPC` = `FALSE`). The function reads the content of the
#'   `Grid_10_Land_Crop.RData` and `Grid_10_Land_Crop_sf_Country.RData` files
#'    - **`DP_R_Grid_Ref`** or **`DP_R_Grid_Ref_Local`**: The function reads the
#'   content of `Grid_10_sf.RData` file from this path.
#'    - **`DP_R_PA`** or **`DP_R_PA_Local`**: The function reads the contents
#'   of the `Sp_PA_Summary_DF.RData` file from this path.
#'    - **`DP_R_CLC_Summary`** / **`DP_R_CLC_Summary_Local`**: Path containing
#'   the `PercCov_SynHab_Crop.RData` file. This file contains maps for the
#'   percentage coverage of each SynHab habitat type per grid cell.
#'    - **`DP_R_CHELSA_Output`** / **`DP_R_CHELSA_Output_Local`**: Path
#'   for processed CHELSA data.
#'    - **`DP_R_Roads`** / **`DP_R_Roads_Local`**: Path for processed road data.
#'   The function reads the contents of: `Road_Length.RData` for the total
#'   length of any road type per grid cell.
#'    - **`DP_R_Railway`** / **`DP_R_Railway_Local`**: Path for processed
#'   railway data. The function reads the contents of: `Railway_Length.RData`
#'   for the total length of any railway type per grid cell.
#'    - **`DP_R_Efforts`** / **`DP_R_Efforts_Local`**: Path for processed
#'   sampling efforts analysis. The function reads the content of
#'   `Bias_GBIF_SummaryR.RData` file containing the total number of GBIF
#'   vascular plant observations per grid cell.
#' @export

Mod_PrepData <- function(
    Hab_Abb = NULL, MinEffortsSp = 100L, PresPerSpecies = 80L, EnvFile = ".env",
    Path_Model = NULL, VerboseProgress = TRUE, FromHPC = TRUE,
    SaveData = TRUE, ExcludeCult = TRUE) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 2023
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||
  # We decided to exclude the following SynHab habitat types in our analysis:
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

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Check input parameters ----
  # # |||||||||||||||||||||||||||||||||||

  CheckNULL <- c("Hab_Abb", "Path_Model", "EnvFile")
  IsNull <- purrr::map_lgl(CheckNULL, ~ is.null(get(.x)))
  if (any(IsNull)) {
    stop(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"),
      call. = FALSE)
  }

  Hab_Abb <- as.character(Hab_Abb)

  if (isFALSE(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpeciesID <- Species_name <- Species_File <- PA <-
    cell <- Path_PA <- Path_Grid <- Path_Grid_Ref <- Path_CLC <-
    Path_Roads <- Path_Rail <- Path_Bias <- Path_CHELSA <-
    EU_Bound <- SpPA <- NPres <- Grid_R <- IAS_ID <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime("Checking input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c("EnvFile", "Hab_Abb", "Path_Model")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("MinEffortsSp", "PresPerSpecies"),
    Type = "numeric")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Reading/checking environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  # Input data paths - these are read from the .env file

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_PA", "DP_R_PA", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC", TRUE, FALSE,
      "Path_CHELSA", "DP_R_CHELSA_Output", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads", TRUE, FALSE,
      "Path_Rail", "DP_R_Railways", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC_Local", TRUE, FALSE,
      "Path_CHELSA", "DP_R_CHELSA_Output_Local", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads_Local", TRUE, FALSE,
      "Path_Rail", "DP_R_Railways_Local", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts_Local", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Checking arguments ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking arguments")

  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0(
        "Hab_Abb has to be one of the following:\n >> ",
        paste0(ValidHabAbbs, collapse = ", ")),
      call. = FALSE)
  }

  fs::dir_create(Path_Model)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Loading data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading data")

  ## Sampling efforts ----
  IASDT.R::CatTime("Sampling efforts", Level = 1)
  R_Bias <- file.path(Path_Bias, "Efforts_SummaryR.RData")
  if (!file.exists(R_Bias)) {
    stop(paste0(R_Bias, " file does not exist"), call. = FALSE)
  }

  # This mask layer represents grid cells with minimum accepted efforts
  # (`MinEffortsSp`). All subsequent maps will be masked to this layer
  EffortsMask <- IASDT.R::LoadAs(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NSp") %>%
    terra::classify(
      rcl = matrix(
        c(0, MinEffortsSp, NA, MinEffortsSp, Inf, 1), byrow = TRUE, ncol = 3),
      include.lowest = TRUE, right = FALSE)

  # # ..................................................................... ###

  ## Species data summary ----
  IASDT.R::CatTime("Species data summary", Level = 1)

  # Extract the list of species for the current habitat type
  DT_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (!file.exists(DT_Sp)) {
    stop(paste0(DT_Sp, " file does not exist"), call. = FALSE)
  }
  DT_Sp <- IASDT.R::LoadAs(DT_Sp) %>%
    dplyr::arrange(IAS_ID)

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

  # # ..................................................................... ###

  ## Species Presence-absence data ----
  IASDT.R::CatTime("Species Presence-absence data", Level = 1)

  # minimum number of presence grids per species
  NCellsCol <- dplyr::if_else(ExcludeCult, "NCells_Naturalized", "NCells_All")

  R_Sp <- DT_Sp %>%
    # Exclude species with too few presence grid cells. There will be further
    # exclusion of species with few grid cells in this pipeline, but excluding
    # this first may help to reduce processing time
    dplyr::filter(.data[[NCellsCol]] >= PresPerSpecies) %>%
    dplyr::select(SpeciesID, Species_name, Species_File) %>%
    # Mask each species map with the filtered grid cells
    dplyr::mutate(
      PA = purrr::map2(
        .x = Species_File, .y = SpeciesID,
        .f = ~{

          PA_Layer <- dplyr::if_else(ExcludeCult, "PA_Masked", "PA")

          # Masked raster map
          SpPA <- file.path(
            Path_PA, "RData", stringr::str_c(.x, "_PA.RData")) %>%
            IASDT.R::LoadAs() %>%
            terra::unwrap() %>%
            magrittr::extract2(PA_Layer) %>%
            terra::mask(EffortsMask) %>%
            stats::setNames(paste0("Sp_", .y))

          # Number of presence grid cells after masking
          NPres <- as.integer(terra::global(SpPA, "sum", na.rm = TRUE))

          return(tibble::tibble(SpPA = list(SpPA), NPres = NPres))
        })) %>%
    tidyr::unnest_wider(PA) %>%
    dplyr::mutate(SpPA = unlist(SpPA)) %>%
    # filter species with too few observations (after masking)
    dplyr::filter(NPres >= PresPerSpecies) %>%
    dplyr::pull(SpPA) %>%
    terra::rast()

  # # ................................... ###

  # Change the efforts mask to exclude grid cells with no species
  #
  # This is not necessary anymore. By excluding less sampled grid cells, there
  # is no need to further exclude grid cells with no species as this could be
  # for an ecological reason
  #
  # ZeroSpGrids <- (sum(R_Sp, na.rm = TRUE) == 0)
  # EffortsMask[ZeroSpGrids] <- NA
  # R_Sp[ZeroSpGrids] <- NA

  # # ................................... ###

  ### Plotting number of IAS per grid cell -----
  IASDT.R::CatTime("Plotting number of IAS per grid cell", Level = 2)

  EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_03")
  R_Sp_sum <- sum(R_Sp, na.rm = TRUE)

  NGridsWzSpecies <- terra::global(R_Sp_sum, fun = "notNA") %>%
    as.integer() %>%
    format(big.mark = ",")
  Limits <- terra::trim(R_Sp_sum) %>%
    terra::ext() %>%
    as.vector()

  NSpPerGrid_gg <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = R_Sp_sum) +
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
        " vascular plant species in GBIF and IAS with ",
        "&#8805;", PresPerSpecies, " presence grid cells are considered (",
        NGridsWzSpecies, " grid cells & ", terra::nlyr(R_Sp), " IAS)")) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = Limits[c(3, 4)]) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = Limits[c(1, 2)]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.25, 0, 0.25, 0, "cm"),
      plot.title = ggtext::element_markdown(
        size = 16, color = "blue",
        margin = ggplot2::margin(0, 0, 0.1, 0, "cm")),
      plot.caption = ggtext::element_markdown(
        size = 11, colour = "grey40", hjust = 0.3),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.9),
      legend.key.size = grid::unit(0.8, "cm"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank())

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_Model, "NSpPerGrid.jpeg"),
    width = 25, height = 27, units = "cm", quality = 100, res = 600)
  print(NSpPerGrid_gg)
  grDevices::dev.off()

  rm(Limits, NSpPerGrid_gg, R_Sp_sum, EU_Bound, envir = environment())

  # # ..................................................................... ###

  ## CHELSA -----

  IASDT.R::CatTime("CHELSA", Level = 1)
  R_CHELSA <- file.path(Path_CHELSA, "Processed", "R_Current.RData")
  if (!file.exists(R_CHELSA)) {
    stop(paste0(R_CHELSA, " file does not exist"), call. = FALSE)
  }
  R_CHELSA <- IASDT.R::LoadAs(R_CHELSA) %>%
    terra::unwrap() %>%
    terra::mask(EffortsMask)

  # # ..................................................................... ###

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
    dplyr::mutate(
      x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble()

  # # ..................................................................... ###

  ## Habitat coverage -----

  IASDT.R::CatTime("Habitat coverage", Level = 1)

  Path_Hab <- file.path(Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    stop(
      paste0("Path_Hab file: ", Path_Hab, " does not exist"),
      call. = FALSE)
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

  # # ..................................................................... ###

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

  # # ..................................................................... ###

  ## Railways ----
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

  # # ..................................................................... ###

  ## Road + rail ----

  IASDT.R::CatTime("Railway + road intensity", Level = 1)
  R_RoadRail <- (R_RoadInt + R_RailInt) %>%
    stats::setNames("RoadRail")
  R_RoadRailLog <- log10(R_RoadRail + 0.1) %>%
    stats::setNames("RoadRailLog")

  # # ..................................................................... ###

  ## Sampling effort ----
  IASDT.R::CatTime("Sampling effort", Level = 1)

  R_Efforts <- IASDT.R::LoadAs(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NObs") %>%
    stats::setNames("Efforts")
  R_EffortsLog <- log10(R_Efforts + 0.1) %>%
    stats::setNames("EffortsLog")

  # # ..................................................................... ###

  ## Merging data together -----
  IASDT.R::CatTime("Merging data together")

  ColumnsFirst <- c("CellNum", "CellCode", "Country", "Country_Nearest")
  DT_All <- c(
    R_CHELSA, R_Hab, R_HabLog, R_RoadInt, R_RoadIntLog,
    R_RailInt, R_RailIntLog, R_RoadRail, R_RoadRailLog,
    R_Efforts, R_EffortsLog, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble() %>%
    # Add country name
    dplyr::left_join(Grid_CNT, by = c("x", "y")) %>%
    dplyr::rename(Country_Nearest = "Nearest", CellNum = cell) %>%
    dplyr::select(tidyselect::all_of(ColumnsFirst), tidyselect::everything())

  if (Hab_Abb == "0") {
    DT_All$Hab <- NA_real_
  }

  # # ..................................................................... ###

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

  # # ..................................................................... ###

  return(DT_All)
}
