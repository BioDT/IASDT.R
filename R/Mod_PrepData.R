## |------------------------------------------------------------------------| #
# Mod_PrepData ----
## |------------------------------------------------------------------------| #

#' @export
#' @name Mod_inputs
#' @rdname Mod_inputs
#' @order 1
#' @author Ahmed El-Gabbas

Mod_PrepData <- function(
    Hab_Abb = NULL, DirName = NULL, MinEffortsSp = 100L, ExcludeCult = TRUE,
    ExcludeZeroHabitat = TRUE, PresPerSpecies = 80L, EnvFile = ".env",
    VerboseProgress = TRUE) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Check input parameters ----
  # # |||||||||||||||||||||||||||||||||||

  CheckNULL <- c("Hab_Abb", "DirName", "EnvFile")
  IsNull <- purrr::map_lgl(CheckNULL, ~ is.null(get(.x)))
  if (any(IsNull)) {
    stop(
      paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
      " can not be empty", call. = FALSE)
  }

  Hab_Abb <- as.character(Hab_Abb)

  if (isFALSE(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpeciesID <- Species_name <- Species_File <- PA <- Path_Rivers <-
    cell <- Path_PA <- Path_Grid <- Path_Grid_Ref <- Path_CLC <-
    Path_Roads <- Path_Rail <- Path_Bias <- Path_CHELSA <- Path_Model <-
    EU_Bound <- SpPA <- NPres <- Grid_R <- IAS_ID <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime("Checking input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c("EnvFile", "Hab_Abb", "DirName")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("MinEffortsSp", "PresPerSpecies"),
    Type = "numeric")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("ExcludeCult", "ExcludeZeroHabitat", "VerboseProgress"))

  # Valid habitat type values
  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      "Hab_Abb has to be one of the following:\n >> ", toString(ValidHabAbbs),
      call. = FALSE)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Reading/checking environment variables ----
  IASDT.R::CatTime("Reading/checking environment variables")

  # # |||||||||||||||||||||||||||||||||||

  # Input data paths - these are read from the .env file

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", FALSE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_PA", "DP_R_PA", TRUE, FALSE,
    "Path_CLC", "DP_R_CLC_processed", TRUE, FALSE,
    "Path_CHELSA", "DP_R_CHELSA_processed", TRUE, FALSE,
    "Path_Roads", "DP_R_Roads_processed", TRUE, FALSE,
    "Path_Rail", "DP_R_Railways_processed", TRUE, FALSE,
    "Path_Bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "Path_Rivers", "DP_R_Rivers_processed", FALSE, TRUE,
    "Path_Model", "DP_R_Model_path", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Path_Model <- IASDT.R::Path(Path_Model, DirName)
  if (fs::dir_exists(Path_Model)) {
    stop("Model directory already exists: ", Path_Model, call. = FALSE)
  }
  fs::dir_create(Path_Model)

  IASDT.R::RecordArgs(
    ExportPath = IASDT.R::Path(Path_Model, "Args_Mod_PrepData.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Loading data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading data")

  ## Sampling efforts ----
  IASDT.R::CatTime("Sampling efforts", Level = 1)
  R_Bias <- IASDT.R::Path(Path_Bias, "Efforts_SummaryR.RData")
  if (!file.exists(R_Bias)) {
    stop(R_Bias, " file does not exist", call. = FALSE)
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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  IASDT.R::CatTime("Habitat coverage", Level = 1)

  Path_Hab <- IASDT.R::Path(
    Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    stop("Path_Hab file: ", Path_Hab, " does not exist", call. = FALSE)
  }

  if (Hab_Abb == "0") {
    # Use dummy habitat values
    R_Hab <- stats::setNames(Grid_R, "Hab")
    R_HabLog <- stats::setNames(Grid_R, "HabLog")
  } else {
    # Load habitat coverage data and mask by the efforts mask
    R_Hab <- IASDT.R::LoadAs(Path_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
      terra::mask(EffortsMask) %>%
      stats::setNames("Hab")

    # Exclude grid cells with zero habitat coverage
    if (ExcludeZeroHabitat) {
      IASDT.R::CatTime(
        "Exclude grid cells with zero habitat coverage", Level = 2)
      ZeroHabGrids <- (R_Hab == 0)
      # Update the efforts mask to exclude grid cells with zero habitat coverage
      EffortsMask[ZeroHabGrids] <- NA
      R_Hab[ZeroHabGrids] <- NA
    }

    # Log of habitat coverage
    R_HabLog <- log10(R_Hab + 0.1) %>%
      stats::setNames("HabLog")
  }

  # # ..................................................................... ###

  ## Species data summary ----
  IASDT.R::CatTime("Species data summary", Level = 1)

  # Extract the list of species for the current habitat type
  DT_Sp <- IASDT.R::Path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (!file.exists(DT_Sp)) {
    stop(DT_Sp, " file does not exist", call. = FALSE)
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

  ## Species presence-absence data ----
  IASDT.R::CatTime("Species presence-absence data", Level = 1)

  # Minimum number of presence grids per species
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
          SpPA <- IASDT.R::Path(
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

  Xlim <- c(2600000, 6550000)
  Ylim <- c(1450000, 5420000)

  Caption <- stringr::str_glue(
    "<span style='color:red; font-size:18px;'>\\
    **{NGridsWzSpecies} grid cells --- {terra::nlyr(R_Sp)} IAS**</span><br/>\\
    - Excluding grid cells with < {MinEffortsSp} vascular plant species \\
    in GBIF",
    dplyr::if_else(
      ExcludeZeroHabitat,
      " or with 0% habitat coverage<br/>", "<br/>"),
    "- Considering only IAS with &#8805; {PresPerSpecies} presence grid cells")

  NSpPerGrid_gg <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = EU_Bound, fill = "gray95", colour = "black", linewidth = 0.15) +
    tidyterra::geom_spatraster(data = R_Sp_sum) +
    tidyterra::scale_fill_whitebox_c(
      na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
    ggplot2::labs(
      title = paste0(
        '<span style="color:blue; font-size:20px;"><b>',
        "Number of IAS per grid cell to be used in the model</b></span>",
        '<span style="color:black; font-size:16px;"> (',
        stringr::str_remove(Hab_column, "Hab_"), ")</span>"),
      caption = Caption) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = Ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = Xlim) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.25, 0, 0.125, 0, "cm"),
      plot.title = ggtext::element_markdown(
        size = 16, color = "blue",
        margin = ggplot2::margin(0, 0, 0.1, 0.25, "cm")),
      plot.caption = ggtext::element_markdown(
        size = 11, colour = "grey40", hjust = 0, lineheight = 1.3,
        margin = ggplot2::margin(0.1, 0, 0, 0.25, "cm")),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.9),
      legend.key.size = grid::unit(0.8, "cm"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank())

  ragg::agg_jpeg(
    filename = IASDT.R::Path(Path_Model, "NSpPerGrid.jpeg"),
    width = 25.5, height = 28, res = 600, quality = 100, units = "cm")
  print(NSpPerGrid_gg)
  grDevices::dev.off()

  rm(NSpPerGrid_gg, R_Sp_sum, EU_Bound, envir = environment())

  # # ..................................................................... ###

  ## CHELSA -----

  IASDT.R::CatTime("CHELSA", Level = 1)
  R_CHELSA <- IASDT.R::Path(Path_CHELSA, "Processed", "R_Current.RData")
  if (!file.exists(R_CHELSA)) {
    stop(R_CHELSA, " file does not exist", call. = FALSE)
  }
  R_CHELSA <- IASDT.R::LoadAs(R_CHELSA) %>%
    terra::unwrap() %>%
    terra::mask(EffortsMask)

  # # ..................................................................... ###

  ## Reference grid -----
  IASDT.R::CatTime("Reference grid", Level = 1)

  IASDT.R::CatTime("Reference grid - sf", Level = 2)
  # Reference grid as sf
  Grid_SF <- IASDT.R::Path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(Grid_SF)) {
    stop(Grid_SF, " file does not exist", call. = FALSE)
  }
  Grid_SF <- IASDT.R::LoadAs(Grid_SF) %>%
    magrittr::extract2("Grid_10_sf_s")

  # # ||||||||||||||||||||||||||||||||||||||||||

  # Reference grid as sf - country names
  IASDT.R::CatTime("Reference grid - country names", Level = 2)
  Grid_CNT <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData")
  if (!file.exists(Grid_CNT)) {
    stop(Grid_CNT, " file does not exist", call. = FALSE)
  }
  Grid_CNT <- IASDT.R::LoadAs(Grid_CNT) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble()

  # # ..................................................................... ###

  ## Railway + road intensity ----
  IASDT.R::CatTime("Railway + road intensity", Level = 1)

  ### Road ----
  IASDT.R::CatTime("Road intensity", Level = 2)

  # road intensity of any road type
  R_RoadInt <- IASDT.R::Path(Path_Roads, "Road_Length.RData")
  if (!file.exists(R_RoadInt)) {
    stop(R_RoadInt, " file does not exist", call. = FALSE)
  }
  R_RoadInt <- IASDT.R::LoadAs(R_RoadInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt") %>%
    terra::mask(EffortsMask)

  # Log of road intensity
  R_RoadIntLog <- log10(R_RoadInt + 0.1) %>%
    stats::setNames("RoadIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Railways ----
  IASDT.R::CatTime("Railway intensity", Level = 2)

  # Railway intensity
  R_RailInt <- IASDT.R::Path(Path_Rail, "Railways_Length.RData")
  if (!file.exists(R_RailInt)) {
    stop(R_RailInt, " file does not exist", call. = FALSE)
  }
  R_RailInt <- IASDT.R::LoadAs(R_RailInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("rail") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("RailInt")

  # Log of railway intensity
  R_RailIntLog <- log10(R_RailInt + 0.1) %>%
    stats::setNames("RailIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Merging Road + rail ----
  IASDT.R::CatTime("Merging Railway + road intensity", Level = 2)
  R_RoadRail <- stats::setNames((R_RoadInt + R_RailInt), "RoadRail")
  R_RoadRailLog <- stats::setNames(log10(R_RoadRail + 0.1), "RoadRailLog")

  # # ..................................................................... ###

  ## Sampling effort ----
  IASDT.R::CatTime("Sampling effort", Level = 1)

  R_Efforts <- IASDT.R::LoadAs(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NObs") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("Efforts")
  R_EffortsLog <- log10(R_Efforts + 0.1) %>%
    stats::setNames("EffortsLog")

  # # ..................................................................... ###

  ## River length ----
  IASDT.R::CatTime("River length", Level = 1)

  R_Rivers <- IASDT.R::Path(Path_Rivers, "River_Lengths.RData")
  if (!file.exists(R_Rivers)) {
    stop(R_Rivers, " file does not exist", call. = FALSE)
  }
  R_Rivers <- IASDT.R::LoadAs(R_Rivers) %>%
    terra::unwrap() %>%
    magrittr::extract2("STRAHLER_5") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("Rivers")
  R_RiversLog <- log10(R_Rivers + 0.1) %>%
    stats::setNames("RiversLog")

  # # ..................................................................... ###

  ## Merging data together -----
  IASDT.R::CatTime("Merging data together")

  ColumnsFirst <- c("CellNum", "CellCode", "Country", "Country_Nearest")
  DT_All <- c(
    R_CHELSA, R_Hab, R_HabLog, R_RoadInt, R_RoadIntLog,
    R_RailInt, R_RailIntLog, R_RoadRail, R_RoadRailLog,
    R_Efforts, R_EffortsLog, R_Rivers, R_RiversLog, R_Sp) %>%
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

  IASDT.R::CatTime("Save model data to disk")
  IASDT.R::SaveAs(
    InObj = DT_All, OutObj = "ModDT",
    OutPath = IASDT.R::Path(Path_Model, "ModDT.RData"))

  # # ..................................................................... ###

  return(DT_All)
}
