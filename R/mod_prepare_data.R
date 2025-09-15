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

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  ecokit::check_args(args_to_check = "directory_name", args_type = "character")
  ecokit::check_args(
    args_to_check = c("min_efforts_n_species", "n_pres_per_species"),
    args_type = "numeric")
  ecokit::check_args(
    args_to_check = c(
      "exclude_cultivated", "exclude_0_habitat", "verbose_progress"),
    args_type = "logical")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  species_id <- species_file <- PA <- path_river <- cell <- path_PA <-
    path_grid <- path_grid_ref <- path_clc <- path_roads <- path_railway <-
    path_efforts <- path_chelsa <- path_model <- EU_boundaries <- sp_PA <-
    n_pres <- grid_r <- ias_id <- path_soil <- path_wetness <- NULL

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Reading/checking environment variables ----
  ecokit::cat_time(
    "Reading/checking environment variables", verbose = verbose_progress)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  # Input data paths - these are read from the .env file
  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", FALSE, FALSE,
    "path_grid_ref", "DP_R_grid_raw", TRUE, FALSE,
    "path_PA", "DP_R_PA", TRUE, FALSE,
    "path_clc", "DP_R_clc_processed", TRUE, FALSE,
    "path_chelsa", "DP_R_chelsa_processed", TRUE, FALSE,
    "path_roads", "DP_R_road_processed", TRUE, FALSE,
    "path_railway", "DP_R_railway_processed", TRUE, FALSE,
    "path_efforts", "DP_R_efforts_processed", TRUE, FALSE,
    "path_river", "DP_R_river_processed", FALSE, TRUE,
    "path_wetness", "DP_R_wetness_processed", FALSE, TRUE,
    "path_soil", "DP_R_soil_density", FALSE, TRUE,
    "path_model", "DP_R_model_root_path", TRUE, FALSE,
    "EU_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  path_model <- fs::path(path_model, directory_name)
  fs::dir_create(path_model)

  ecokit::record_arguments(
    out_path = fs::path(path_model, "Args_prepare_data.RData"))

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Loading data ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Loading data", verbose = verbose_progress)

  ## Sampling efforts ----
  ecokit::cat_time("Sampling efforts", level = 1L, verbose = verbose_progress)
  R_Bias <- fs::path(path_efforts, "efforts_summary_r.RData")
  if (!fs::file_exists(R_Bias)) {
    ecokit::stop_ctx(
      "R_Bias file does not exist", R_Bias = R_Bias, include_backtrace = TRUE)
  }

  # This mask layer represents grid cells with minimum accepted efforts
  # (`min_efforts_n_species`). All subsequent maps will be masked to this layer
  efforts_mask <- ecokit::load_as(R_Bias, unwrap_r = TRUE) %>%
    magrittr::extract2("NSp") %>%
    terra::classify(
      rcl = matrix(
        c(0, min_efforts_n_species, NA, min_efforts_n_species, Inf, 1),
        byrow = TRUE, ncol = 3),
      include.lowest = TRUE, right = FALSE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Habitat coverage -----

  ecokit::cat_time("Habitat coverage", level = 1L, verbose = verbose_progress)

  path_habitat <- fs::path(
    path_clc, "summary_rdata", "perc_cover_synhab_crop.RData")
  if (!fs::file_exists(path_habitat)) {
    ecokit::stop_ctx(
      "path_habitat file does not exist", path_habitat = path_habitat,
      include_backtrace = TRUE)
  }

  if (hab_abb == "0") {
    # Use dummy habitat values
    r_habitat <- stats::setNames(grid_r, "Hab")
    r_habitat_log <- stats::setNames(grid_r, "HabLog")
  } else {
    # Load habitat coverage data and mask by the efforts mask
    r_habitat <- ecokit::load_as(path_habitat, unwrap_r = TRUE) %>%
      magrittr::extract2(paste0("SynHab_", hab_abb)) %>%
      terra::mask(efforts_mask) %>%
      stats::setNames("Hab")

    # Exclude grid cells with zero habitat coverage
    if (exclude_0_habitat) {
      ecokit::cat_time(
        "Exclude grid cells with zero habitat coverage",
        level = 2L, verbose = verbose_progress)
      zero_hab_grids <- (r_habitat == 0)
      # Update the efforts mask to exclude grid cells with zero habitat coverage
      efforts_mask[zero_hab_grids] <- NA
      r_habitat[zero_hab_grids] <- NA
    }

    # Log of habitat coverage
    r_habitat_log <- stats::setNames(log10(r_habitat + 0.1), "HabLog")
  }

  # # ..................................................................... ###

  ## Species data summary ----
  ecokit::cat_time(
    "Species data summary", level = 1L, verbose = verbose_progress)

  # Extract the list of species for the current habitat type
  DT_Sp <- fs::path(path_PA, "sp_pa_summary_df.RData")
  if (!fs::file_exists(DT_Sp)) {
    ecokit::stop_ctx(
      "DT_Sp file does not exist", DT_Sp = DT_Sp, include_backtrace = TRUE)
  }
  DT_Sp <- ecokit::load_as(DT_Sp) %>%
    dplyr::arrange(ias_id)

  if (hab_abb == "0") {
    Hab_column <- NULL
  } else {
    Hab_column <- c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("_", hab_abb, "_"))

    DT_Sp <- dplyr::filter(DT_Sp, !!as.symbol(Hab_column))
  }

  # # ..................................................................... ###

  ## Species presence-absence data ----
  ecokit::cat_time(
    "Species presence-absence data", level = 1L, verbose = verbose_progress)

  # Minimum number of presence grids per species
  n_cells_col <- dplyr::if_else(
    exclude_cultivated, "n_cells_naturalized", "n_cells_all")

  R_Sp <- DT_Sp %>%
    # Exclude species with too few presence grid cells. There will be further
    # exclusion of species with few grid cells in this pipeline, but excluding
    # this first may help to reduce processing time
    dplyr::filter(.data[[n_cells_col]] >= n_pres_per_species) %>%
    dplyr::select(
      tidyselect::all_of(c("species_id", "species_name", "species_file"))) %>%
    # Mask each species map with the filtered grid cells
    dplyr::mutate(
      PA = purrr::map2(
        .x = species_file, .y = species_id,
        .f = ~{

          PA_Layer <- dplyr::if_else(exclude_cultivated, "PA_masked", "PA")

          # Masked raster map
          sp_PA <- fs::path(
            path_PA, "PA_rdata", stringr::str_c(.x, "_PA.RData")) %>%
            ecokit::load_as(unwrap_r = TRUE) %>%
            magrittr::extract2(PA_Layer) %>%
            terra::mask(efforts_mask) %>%
            stats::setNames(paste0("Sp_", .y))

          # Number of presence grid cells after masking
          n_pres <- as.integer(terra::global(sp_PA, "sum", na.rm = TRUE))

          return(tibble::tibble(sp_PA = list(sp_PA), n_pres = n_pres))
        })) %>%
    tidyr::unnest_wider(PA) %>%
    dplyr::mutate(sp_PA = unlist(sp_PA)) %>%
    # filter species with too few observations (after masking)
    dplyr::filter(n_pres >= n_pres_per_species) %>%
    dplyr::pull(sp_PA) %>%
    terra::rast()

  # # ................................... ###

  # Change the efforts mask to exclude grid cells with no species
  #
  # This is not necessary any more. By excluding less sampled grid cells, there
  # is no need to further exclude grid cells with no species as this could be
  # for an ecological reason
  #
  # zero_sp_grids <- (sum(R_Sp, na.rm = TRUE) == 0)
  # efforts_mask[zero_sp_grids] <- NA
  # R_Sp[zero_sp_grids] <- NA

  # # ................................... ###

  ### Plotting number of IAS per grid cell -----
  ecokit::cat_time(
    "Plotting number of IAS per grid cell",
    level = 2L, verbose = verbose_progress)

  EU_boundaries <- ecokit::load_as(EU_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_03")
  R_Sp_sum <- sum(R_Sp, na.rm = TRUE)

  NGridsWzSpecies <- terra::global(R_Sp_sum, fun = "notNA") %>%      # nolint: object_name_linter
    as.integer() %>%
    format(big.mark = ",")

  Xlim <- c(2600000, 6550000)
  Ylim <- c(1450000, 5420000)

  plot_caption <- stringr::str_glue(
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

  n_sp_per_grid_gg <- ggplot2::ggplot(environment = emptyenv()) +
    ggplot2::geom_sf(
      data = EU_boundaries, fill = "gray95",
      colour = "black", linewidth = 0.15) +
    tidyterra::geom_spatraster(data = R_Sp_sum) +
    tidyterra::scale_fill_whitebox_c(
      na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
    ggplot2::labs(
      title = paste0(
        '<span style="color:blue; font-size:20px;"><b>',
        "Number of IAS per grid cell to be used in the model</b></span>",
        '<span style="color:black; font-size:16px;"> (',
        stringr::str_remove(Hab_column, "Hab_"), ")</span>"),
      caption = plot_caption) +
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
  print(n_sp_per_grid_gg)
  grDevices::dev.off()

  rm(n_sp_per_grid_gg, R_Sp_sum, EU_boundaries, envir = environment())

  # # ..................................................................... ###

  ## CHELSA -----

  ecokit::cat_time("CHELSA", level = 1L, verbose = verbose_progress)
  r_chelsa <- fs::path(path_chelsa, "Processed", "R_Current.RData")
  if (!fs::file_exists(r_chelsa)) {
    ecokit::stop_ctx(
      "r_chelsa file does not exist", r_chelsa = r_chelsa,
      include_backtrace = TRUE)
  }
  r_chelsa <- ecokit::load_as(r_chelsa, unwrap_r = TRUE) %>%
    terra::mask(efforts_mask)

  # # ..................................................................... ###

  ## Reference grid -----
  ecokit::cat_time("Reference grid", level = 1L, verbose = verbose_progress)

  ecokit::cat_time(
    "Reference grid - sf", level = 2L, verbose = verbose_progress)
  # Reference grid as sf
  grid_sf <- fs::path(path_grid_ref, "Grid_10_sf.RData")
  if (!fs::file_exists(grid_sf)) {
    ecokit::stop_ctx(
      "grid_sf file does not exist", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }
  grid_sf <- ecokit::load_as(grid_sf) %>%
    magrittr::extract2("Grid_10_sf_s")

  # # ||||||||||||||||||||||||||||||||||||||||||

  # Reference grid as sf - country names
  ecokit::cat_time(
    "Reference grid - country names", level = 2L, verbose = verbose_progress)
  grid_country <- fs::path(path_grid, "grid_10_land_crop_sf_country.RData")
  if (!fs::file_exists(grid_country)) {
    ecokit::stop_ctx(
      "grid_country file does not exist", grid_country = grid_country,
      include_backtrace = TRUE)
  }
  grid_country <- ecokit::load_as(grid_country) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) %>%
    sf::st_drop_geometry() %>%
    tibble::tibble()

  # # ..................................................................... ###

  ## Railway + road intensity ----
  ecokit::cat_time(
    "Railway + road intensity", level = 1L, verbose = verbose_progress)

  ### Road ----
  ecokit::cat_time("Road intensity", level = 2L, verbose = verbose_progress)

  # road intensity of any road type
  r_road_int <- fs::path(path_roads, "poad_length.RData")
  if (!fs::file_exists(r_road_int)) {
    ecokit::stop_ctx(
      "r_road_int file does not exist", r_road_int = r_road_int,
      include_backtrace = TRUE)
  }
  r_road_int <- ecokit::load_as(r_road_int, unwrap_r = TRUE) %>%
    magrittr::extract2("All") %>%
    stats::setNames("RoadInt") %>%
    terra::mask(efforts_mask)

  # Log of road intensity
  # add 1 (older versions 0.1) to get log for 0 values
  # [only for rivers/roads/efforts, not hab/rivers]
  r_road_int_log <- stats::setNames(log10(r_road_int + 1), "RoadIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Railways ----
  ecokit::cat_time("Railway intensity", level = 2L, verbose = verbose_progress)

  # Railway intensity
  r_rail_int <- fs::path(path_railway, "railway_length.RData")
  if (!fs::file_exists(r_rail_int)) {
    ecokit::stop_ctx(
      "r_rail_int file does not exist", r_rail_int = r_rail_int,
      include_backtrace = TRUE)
  }
  r_rail_int <- ecokit::load_as(r_rail_int, unwrap_r = TRUE) %>%
    magrittr::extract2("rail") %>%
    terra::mask(efforts_mask) %>%
    stats::setNames("RailInt")

  # Log of railway intensity
  # add 1 (older versions 0.1) to get log for 0 values
  # [only for rivers/roads/efforts, not hab/rivers]
  r_rail_int_log <- stats::setNames(log10(r_rail_int + 1), "RailIntLog")

  # # ||||||||||||||||||||||||||||||||||||||||||

  ### Merging Road + rail ----
  ecokit::cat_time(
    "Merging Railway + road intensity", level = 2L, verbose = verbose_progress)
  r_road_rail <- stats::setNames((r_road_int + r_rail_int), "RoadRail")
  # add 1 (older versions 0.1) to get log for 0 values
  # [only for rivers/roads/efforts, not hab/rivers]
  r_road_rail_log <- stats::setNames(log10(r_road_rail + 1), "RoadRailLog")

  # # ..................................................................... ###

  ## Sampling effort ----
  ecokit::cat_time("Sampling effort", level = 1L, verbose = verbose_progress)

  r_efforts <- ecokit::load_as(R_Bias, unwrap_r = TRUE) %>%
    magrittr::extract2("NObs") %>%
    terra::mask(efforts_mask) %>%
    stats::setNames("Efforts")
  # add 1 (older versions 0.1) to get log for 0 values
  # [only for rivers/roads/efforts, not hab/rivers]
  r_efforts_log <- stats::setNames(log10(r_efforts + 1), "EffortsLog")

  # # ..................................................................... ###

  ## River length ----
  ecokit::cat_time("River length", level = 1L, verbose = verbose_progress)

  r_rivers <- fs::path(path_river, "river_lengths.RData")
  if (!fs::file_exists(r_rivers)) {
    ecokit::stop_ctx(
      "r_rivers file does not exist", r_rivers = r_rivers,
      include_backtrace = TRUE)
  }
  r_rivers <- ecokit::load_as(r_rivers, unwrap_r = TRUE) %>%
    magrittr::extract2("STRAHLER_5") %>%
    terra::mask(efforts_mask) %>%
    stats::setNames("Rivers")
  r_rivers_log <- stats::setNames(log10(r_rivers + 0.1), "RiversLog")

  # # ..................................................................... ###

  ## soil bulk density ----
  ecokit::cat_time("soil bulk density", level = 1L, verbose = verbose_progress)

  r_soil <- fs::path(path_soil, "soil_density.RData")
  if (!fs::file_exists(r_soil)) {
    ecokit::stop_ctx(
      "r_soil file does not exist", r_soil = r_soil,
      include_backtrace = TRUE)
  }
  r_soil <- ecokit::load_as(r_soil, unwrap_r = TRUE) %>%
    magrittr::extract2("bdod_5_15_mean") %>%
    terra::mask(efforts_mask) %>%
    stats::setNames("soil")

  # # ..................................................................... ###

  ## Topographic Wetness Index ----
  ecokit::cat_time(
    "Topographic Wetness Index", level = 1L, verbose = verbose_progress)

  r_wetness <- fs::path(path_wetness, "wetness_index.RData")
  if (!fs::file_exists(r_wetness)) {
    ecokit::stop_ctx(
      "r_wetness file does not exist", r_wetness = r_wetness,
      include_backtrace = TRUE)
  }
  r_wetness <- ecokit::load_as(r_wetness, unwrap_r = TRUE) %>%
    terra::mask(efforts_mask) %>%
    stats::setNames("wetness")

  # # ..................................................................... ###

  ## Merging data together -----
  ecokit::cat_time("Merging data together", verbose = verbose_progress)

  columns_first <- c("CellNum", "CellCode", "country", "country_nearest")
  data_all <- c(
    r_chelsa, r_habitat, r_habitat_log, r_road_int, r_road_int_log,
    r_rail_int, r_rail_int_log, r_road_rail, r_road_rail_log,
    r_efforts, r_efforts_log, r_rivers, r_rivers_log,
    r_soil, r_wetness, R_Sp) %>%
    as.data.frame(na.rm = TRUE, xy = TRUE, cells = TRUE) %>%
    tibble::tibble() %>%
    # Add country name
    dplyr::left_join(grid_country, by = c("x", "y")) %>%
    dplyr::rename(country_nearest = "Nearest", CellNum = cell) %>%
    dplyr::select(tidyselect::all_of(columns_first), tidyselect::everything())

  if (hab_abb == "0") {
    data_all$Hab <- NA_real_
  }

  # # ..................................................................... ###

  # Save model data to disk -----

  ecokit::cat_time("Save model data to disk", verbose = verbose_progress)
  ecokit::save_as(
    object = data_all, object_name = "ModDT",
    out_path = fs::path(path_model, "ModDT.RData"))

  # # ..................................................................... ###

  return(data_all)
}
