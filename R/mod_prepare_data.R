## |------------------------------------------------------------------------| #
# mod_prepare_data ----
## |------------------------------------------------------------------------| #

#' @export
#' @name mod_inputs
#' @rdname mod_inputs
#' @order 1
#' @author Ahmed El-Gabbas

mod_prepare_data <- function(
    hab_abb = NULL, directory_name = NULL, min_efforts_n_species = 100L,
    exclude_cultivated = TRUE, exclude_0_habitat = TRUE,
    n_pres_per_species = 80L, env_file = ".env", verbose_progress = TRUE) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Check input parameters ----
  # # |||||||||||||||||||||||||||||||||||

  CheckNULL <- c("hab_abb", "directory_name", "env_file")
  IsNull <- purrr::map_lgl(CheckNULL, ~ is.null(get(.x)))
  if (any(IsNull)) {
    IASDT.R::stop_ctx(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"),
      hab_abb = hab_abb, directory_name = directory_name, env_file = env_file)
  }

  hab_abb <- as.character(hab_abb)

  if (isFALSE(verbose_progress)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  species_ID <- Species_name <- Species_File <- PA <- Path_Rivers <-
    cell <- Path_PA <- Path_Grid <- Path_Grid_Ref <- Path_CLC <-
    Path_Roads <- Path_Rail <- Path_Bias <- Path_CHELSA <- path_model <-
    EU_Bound <- SpPA <- NPres <- Grid_R <- IAS_ID <- NULL

  # # ..................................................................... ###

  IASDT.R::cat_time("Checking input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  CharArgs <- c("env_file", "hab_abb", "directory_name")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = CharArgs, args_type = "character")
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("min_efforts_n_species", "n_pres_per_species"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "exclude_cultivated", "exclude_0_habitat", "verbose_progress"))

  # Valid habitat type values
  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (!(hab_abb %in% ValidHabAbbs)) {
    IASDT.R::stop_ctx(
      paste0(
        "`hab_abb` has to be one of the following:\n >> ",
        toString(ValidHabAbbs)),
      hab_abb = hab_abb)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Reading/checking environment variables ----
  IASDT.R::cat_time("Reading/checking environment variables")

  # # |||||||||||||||||||||||||||||||||||

  # Input data paths - these are read from the .env file

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", FALSE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_PA", "DP_R_PA", TRUE, FALSE,
    "Path_CLC", "DP_R_CLC_processed", TRUE, FALSE,
    "Path_CHELSA", "DP_R_CHELSA_processed", TRUE, FALSE,
    "Path_Roads", "DP_R_Roads_processed", TRUE, FALSE,
    "Path_Rail", "DP_R_Railways_processed", TRUE, FALSE,
    "Path_Bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "Path_Rivers", "DP_R_Rivers_processed", FALSE, TRUE,
    "path_model", "DP_R_Model_path", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  path_model <- fs::path(path_model, directory_name)
  fs::dir_create(path_model)

  IASDT.R::record_arguments(
    out_path = fs::path(path_model, "Args_prepare_data.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Loading data ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Loading data")

  ## Sampling efforts ----
  IASDT.R::cat_time("Sampling efforts", level = 1L)
  R_Bias <- fs::path(Path_Bias, "Efforts_SummaryR.RData")
  if (!file.exists(R_Bias)) {
    IASDT.R::stop_ctx("R_Bias file does not exist", R_Bias = R_Bias)
  }

  # This mask layer represents grid cells with minimum accepted efforts
  # (`min_efforts_n_species`). All subsequent maps will be masked to this layer
  EffortsMask <- IASDT.R::load_as(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NSp") %>%
    terra::classify(
      rcl = matrix(
        c(0, min_efforts_n_species, NA, min_efforts_n_species, Inf, 1),
        byrow = TRUE, ncol = 3),
      include.lowest = TRUE, right = FALSE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  IASDT.R::cat_time("Habitat coverage", level = 1L)

  Path_Hab <- fs::path(Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    IASDT.R::stop_ctx("Path_Hab file does not exist", Path_Hab = Path_Hab)
  }

  if (hab_abb == "0") {
    # Use dummy habitat values
    R_Hab <- stats::setNames(Grid_R, "Hab")
    R_HabLog <- stats::setNames(Grid_R, "HabLog")
  } else {
    # Load habitat coverage data and mask by the efforts mask
    R_Hab <- IASDT.R::load_as(Path_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", hab_abb)) %>%
      terra::mask(EffortsMask) %>%
      stats::setNames("Hab")

    # Exclude grid cells with zero habitat coverage
    if (exclude_0_habitat) {
      IASDT.R::cat_time(
        "Exclude grid cells with zero habitat coverage", level = 2L)
      zero_hab_grids <- (R_Hab == 0)
      # Update the efforts mask to exclude grid cells with zero habitat coverage
      EffortsMask[zero_hab_grids] <- NA
      R_Hab[zero_hab_grids] <- NA
    }

    # Log of habitat coverage
    R_HabLog <- log10(R_Hab + 0.1) %>%
      stats::setNames("HabLog")
  }

  # # ..................................................................... ###

  ## Species data summary ----
  IASDT.R::cat_time("Species data summary", level = 1L)

  # Extract the list of species for the current habitat type
  DT_Sp <- fs::path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (!file.exists(DT_Sp)) {
    IASDT.R::stop_ctx("DT_Sp file does not exist", DT_Sp = DT_Sp)
  }
  DT_Sp <- IASDT.R::load_as(DT_Sp) %>%
    dplyr::arrange(IAS_ID)

  if (hab_abb == "0") {
    Hab_column <- NULL
  } else {
    Hab_column <- c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("_", as.character(hab_abb), "_"))

    DT_Sp <- dplyr::filter(DT_Sp, !!as.symbol(Hab_column))
  }

  # # ..................................................................... ###

  ## Species presence-absence data ----
  IASDT.R::cat_time("Species presence-absence data", level = 1L)

  # Minimum number of presence grids per species
  NCellsCol <- dplyr::if_else(
    exclude_cultivated, "NCells_Naturalized", "NCells_All")

  R_Sp <- DT_Sp %>%
    # Exclude species with too few presence grid cells. There will be further
    # exclusion of species with few grid cells in this pipeline, but excluding
    # this first may help to reduce processing time
    dplyr::filter(.data[[NCellsCol]] >= n_pres_per_species) %>%
    dplyr::select(species_ID, Species_name, Species_File) %>%
    # Mask each species map with the filtered grid cells
    dplyr::mutate(
      PA = purrr::map2(
        .x = Species_File, .y = species_ID,
        .f = ~{

          PA_Layer <- dplyr::if_else(exclude_cultivated, "PA_Masked", "PA")

          # Masked raster map
          SpPA <- fs::path(
            Path_PA, "RData", stringr::str_c(.x, "_PA.RData")) %>%
            IASDT.R::load_as() %>%
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
    dplyr::filter(NPres >= n_pres_per_species) %>%
    dplyr::pull(SpPA) %>%
    terra::rast()

  # # ................................... ###

  # Change the efforts mask to exclude grid cells with no species
  #
  # This is not necessary anymore. By excluding less sampled grid cells, there
  # is no need to further exclude grid cells with no species as this could be
  # for an ecological reason
  #
  # zero_sp_grids <- (sum(R_Sp, na.rm = TRUE) == 0)
  # EffortsMask[zero_sp_grids] <- NA
  # R_Sp[zero_sp_grids] <- NA

  # # ................................... ###

  ### Plotting number of IAS per grid cell -----
  IASDT.R::cat_time("Plotting number of IAS per grid cell", level = 2L)

  EU_Bound <- IASDT.R::load_as(EU_Bound) %>%
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
    - Excluding grid cells with < {min_efforts_n_species} vascular plant \\
    species in GBIF",
    dplyr::if_else(
      exclude_0_habitat,
      " or with 0% habitat coverage<br/>", "<br/>"),
    paste0(
      "- Considering only IAS with &#8805; {n_pres_per_species}",
      " presence grid cells"))

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
    filename = fs::path(path_model, "NSpPerGrid.jpeg"),
    width = 25.5, height = 28, res = 600, quality = 100, units = "cm")
  print(NSpPerGrid_gg)
  grDevices::dev.off()

  rm(NSpPerGrid_gg, R_Sp_sum, EU_Bound, envir = environment())

  # # ..................................................................... ###

  ## CHELSA -----

  IASDT.R::cat_time("CHELSA", level = 1L)
  R_CHELSA <- fs::path(Path_CHELSA, "Processed", "R_Current.RData")
  if (!file.exists(R_CHELSA)) {
    IASDT.R::stop_ctx("R_CHELSA file does not exist", R_CHELSA = R_CHELSA)
  }
  R_CHELSA <- IASDT.R::load_as(R_CHELSA) %>%
    terra::unwrap() %>%
    terra::mask(EffortsMask)

  # # ..................................................................... ###

  ## Reference grid -----
  IASDT.R::cat_time("Reference grid", level = 1L)

  IASDT.R::cat_time("Reference grid - sf", level = 2L)
  # Reference grid as sf
  Grid_SF <- fs::path(Path_Grid_Ref, "Grid_10_sf.RData")
  if (!file.exists(Grid_SF)) {
    IASDT.R::stop_ctx("Grid_SF file does not exist", Grid_SF = Grid_SF)
  }
  Grid_SF <- IASDT.R::load_as(Grid_SF) %>%
    magrittr::extract2("Grid_10_sf_s")

  # # ||||||||||||||||||||||||||||||||||||||||||

  # Reference grid as sf - country names
  IASDT.R::cat_time("Reference grid - country names", level = 2L)
  Grid_CNT <- fs::path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData")
  if (!file.exists(Grid_CNT)) {
    IASDT.R::stop_ctx("Grid_CNT file does not exist", Grid_CNT = Grid_CNT)
  }
  Grid_CNT <- IASDT.R::load_as(Grid_CNT) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble()

  # # ..................................................................... ###

  ## Railway + road intensity ----
  IASDT.R::cat_time("Railway + road intensity", level = 1L)

  ### Road ----
  IASDT.R::cat_time("Road intensity", level = 2L)

  # road intensity of any road type
  R_RoadInt <- fs::path(Path_Roads, "Road_Length.RData")
  if (!file.exists(R_RoadInt)) {
    IASDT.R::stop_ctx("R_RoadInt file does not exist", R_RoadInt = R_RoadInt)
  }
  R_RoadInt <- IASDT.R::load_as(R_RoadInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt") %>%
    terra::mask(EffortsMask)

  # Log of road intensity
  R_RoadIntLog <- log10(R_RoadInt + 0.1) %>%
    stats::setNames("RoadIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Railways ----
  IASDT.R::cat_time("Railway intensity", level = 2L)

  # Railway intensity
  R_RailInt <- fs::path(Path_Rail, "Railways_Length.RData")
  if (!file.exists(R_RailInt)) {
    IASDT.R::stop_ctx("R_RailInt file does not exist", R_RailInt = R_RailInt)
  }
  R_RailInt <- IASDT.R::load_as(R_RailInt) %>%
    terra::unwrap() %>%
    magrittr::extract2("rail") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("RailInt")

  # Log of railway intensity
  R_RailIntLog <- log10(R_RailInt + 0.1) %>%
    stats::setNames("RailIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Merging Road + rail ----
  IASDT.R::cat_time("Merging Railway + road intensity", level = 2L)
  R_RoadRail <- stats::setNames((R_RoadInt + R_RailInt), "RoadRail")
  R_RoadRailLog <- stats::setNames(log10(R_RoadRail + 0.1), "RoadRailLog")

  # # ..................................................................... ###

  ## Sampling effort ----
  IASDT.R::cat_time("Sampling effort", level = 1L)

  R_Efforts <- IASDT.R::load_as(R_Bias) %>%
    terra::unwrap() %>%
    magrittr::extract2("NObs") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("Efforts")
  R_EffortsLog <- log10(R_Efforts + 0.1) %>%
    stats::setNames("EffortsLog")

  # # ..................................................................... ###

  ## River length ----
  IASDT.R::cat_time("River length", level = 1L)

  R_Rivers <- fs::path(Path_Rivers, "River_Lengths.RData")
  if (!file.exists(R_Rivers)) {
    IASDT.R::stop_ctx("R_Rivers file does not exist", R_Rivers = R_Rivers)
  }
  R_Rivers <- IASDT.R::load_as(R_Rivers) %>%
    terra::unwrap() %>%
    magrittr::extract2("STRAHLER_5") %>%
    terra::mask(EffortsMask) %>%
    stats::setNames("Rivers")
  R_RiversLog <- log10(R_Rivers + 0.1) %>%
    stats::setNames("RiversLog")

  # # ..................................................................... ###

  ## Merging data together -----
  IASDT.R::cat_time("Merging data together")

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

  if (hab_abb == "0") {
    DT_All$Hab <- NA_real_
  }

  # # ..................................................................... ###

  # Save model data to disk -----

  IASDT.R::cat_time("Save model data to disk")
  IASDT.R::save_as(
    object = DT_All, object_name = "ModDT",
    out_path = fs::path(path_model, "ModDT.RData"))

  # # ..................................................................... ###

  return(DT_All)
}
