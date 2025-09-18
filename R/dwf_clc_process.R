## |------------------------------------------------------------------------| #
# clc_process ------
## |------------------------------------------------------------------------| #

#' Process Corine Land Cover (CLC) data for the `IASDT`
#'
#' Processes [Corine Land Cover
#' (CLC)](https://land.copernicus.eu/pan-european/corine-land-cover/clc2018)
#' data for the Invasive Alien Species Digital Twin (`IASDT`). Calculates
#' percentage coverage and most common classes per grid cell at three CLC
#' levels, plus `eunis19` and `SynHab` habitat types. Prepares a reference grid
#' and optionally generates percentage coverage maps as JPEG.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param min_land_percent Numeric. Minimum land percentage per grid cell for
#'   the reference grid. Default: `15`.
#' @param plot_clc Logical. If `TRUE`, plots percentage coverage for CLC levels
#'   and habitat types. Default: `TRUE`.
#' @param n_cores_plot Integer. Number of CPU cores to use for parallel
#'   plotting. Default: 5L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @return Returns `invisible(NULL)`; saves processed data and optional plots to
#'   disk.
#' @name clc_process
#'
#' @author Ahmed El-Gabbas
#' @export
#' @references
#' - Data source:
#' <https://land.copernicus.eu/pan-european/corine-land-cover/clc2018>
#' - Data citation:
#' <https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac>

clc_process <- function(
    env_file = ".env", min_land_percent = 15L, plot_clc = TRUE,
    n_cores_plot = 5L, strategy = "multisession") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid warning while reading CLC data
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CNTR_ID <- geometry <- name <- path_clc <- path_grid <-
    path_grid_ref <- km <- majority <- path_clc_tif <- path_clc_crosswalk <-
    eu_boundaries <- value <- country <- country_2 <- map <- Value <- NULL

  if (!is.numeric(min_land_percent) ||
      !dplyr::between(min_land_percent, 0, 100)) {
    ecokit::stop_ctx(
      "`min_land_percent` must be a numeric value between 0 and 100.",
      min_land_percent = min_land_percent, include_backtrace = TRUE)
  }

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  ## ---------------------------------------------------
  # Loading data ------
  ## ---------------------------------------------------

  ecokit::cat_time("Loading data")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Environment variables
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", FALSE, FALSE,
    "path_grid_ref", "DP_R_grid_raw", TRUE, FALSE,
    "path_clc", "DP_R_clc_processed", FALSE, FALSE,
    "path_clc_tif", "DP_R_clc_tif", FALSE, TRUE,
    "path_clc_crosswalk", "DP_R_clc_crosswalk", FALSE, TRUE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Check files/directories
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Check files and directories", level = 1L)
  fs::dir_create(path_clc)

  path_grid_sf <- fs::path(path_grid_ref, "Grid_10_sf.RData")
  path_grid_rast <- fs::path(path_grid_ref, "Grid_10_Raster.RData")

  required_paths <- c(path_grid_sf, path_clc_crosswalk, path_grid_rast)

  purrr::walk(
    .x = required_paths,
    .f = function(path) {
      if (!file.exists(path)) {
        ecokit::stop_ctx(
          "Required path does not exist", path = path,
          include_backtrace = TRUE)
      }})

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Cross-walk
  # description of CLC values and custom crosswalks
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Loading crosswalks", level = 1L)
  clc_crosswalk <- readr::read_delim(
    file = path_clc_crosswalk, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::rename_all(~stringr::str_replace(.x, "CLC", "clc")) %>%
    dplyr::rename_all(~stringr::str_replace(.x, "_L", "_l")) %>%
    dplyr::rename_all(~stringr::str_replace(.x, "EUNIS_2019", "eunis2019")) %>%
    dplyr::rename_all(~stringr::str_replace(.x, "SynHab", "synhab")) %>%
    dplyr::rename(value = Value)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # reference grid
  # # ||||||||||||||||||||||||||||||||||||||||||||

  # Path of reference grid
  ecokit::cat_time("Loading reference grid", level = 1L)

  # sf
  ecokit::cat_time("sf", level = 2L)

  grid_sf <- ecokit::load_as(path_grid_sf) %>%
    magrittr::extract2("Grid_10_sf_s")

  # raster
  ecokit::cat_time("raster", level = 2L)
  grid_r <- ecokit::load_as(path_grid_rast)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # country boundaries
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Loading country boundaries", level = 1L)
  eu_boundaries_sf <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Creating folders
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Creating folders", level = 1L)

  # sub-folders to store tif files (no masking to final reference grid)
  path_clc_summary_tif <- fs::path(path_clc, "summary_tif")

  # sub-folders to store tif files (masked to final reference grid)
  path_clc_summary_tif_crop <- fs::path(path_clc, "summary_tif_crop")

  # sub-folders to store output RData files
  path_clc_summary_rdata <- fs::path(path_clc, "summary_rdata")

  # Create folders when necessary
  c(
    path_clc, path_clc_summary_tif, path_clc_summary_tif_crop,
    path_clc_summary_rdata, path_grid) %>%
    purrr::walk(fs::dir_create)

  if (plot_clc) {
    # sub-folders to store JPEG files
    path_clc_summary_jpeg <- fs::path(path_clc, "summary_jpeg")
    fs::dir_create(path_clc_summary_jpeg)
  }

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Calculate fraction of each CLC class at each grid cell ----
  # ---------------------------------------------------

  ecokit::cat_time("Calculate fraction of each CLC class at each grid cell")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Loading CLC tif file
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Loading CLC tif file", level = 1L)
  clc_rast <- terra::rast(path_clc_tif)
  terra::NAflag(clc_rast) <- 128

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Calculate fraction for each CLC value
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Processing using `exactextractr::exact_extract` function",
    level = 1L)

  clc_fracs <- grid_sf %>%
    # Ensure that the projection of x and y parameters of
    # `exactextractr::exact_extract` suppress warning: Polygons transformed to
    # raster CRS (EPSG:3035)
    # https://github.com/isciences/exactextractr/issues/103
    sf::st_transform(sf::st_crs(clc_rast)) %>%
    dplyr::mutate(
      exactextractr::exact_extract(
        x = clc_rast, y = ., fun = "frac",
        force_df = TRUE, default_value = 44, progress = FALSE)) %>%
    sf::st_transform(crs = 3035) %>%
    tibble::as_tibble() %>%
    sf::st_as_sf()

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Save fraction results
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save fraction results", level = 1L)
  save(
    clc_fracs,
    file = fs::path(path_clc_summary_rdata, "clc_fracs.RData"))

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Convert fractions to raster
  # # ||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Convert fractions to raster", level = 1L)
  # convert to SpatVector objects for faster rasterization
  clc_fracs_vect <- terra::vect(clc_fracs)

  clc_fracs_r <- names(clc_fracs) %>%
    # Exclude processing CLC values of NODATA "frac_48"
    setdiff(y = c("CellCode", "geometry", "frac_48")) %>%
    purrr::map(
      .f = ~ {
        map <- terra::rasterize(
          x = clc_fracs_vect, y = terra::rast(grid_r), field = .x) %>%
          magrittr::multiply_by(100) %>%
          stats::setNames(.x) %>%
          ecokit::set_raster_crs(crs = "epsg:3035")

        # For `frac_44` (CLC 5.2.3 Sea and ocean), the original values do not
        # cover the whole study area, resulting in values of zero in the sea.
        # The following force seas and oceans outside mainland Europe to have
        # 100 value
        if (.x == "frac_44") {
          map <- terra::trim(map, value = 0) %>%
            terra::extend(terra::rast(grid_r), fill = 100)
        }
        return(map)
      }
    ) %>%
    terra::rast()

  rm(clc_fracs_vect, clc_fracs, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Calculate % coverage of different crosswalks per grid cell ----
  # ---------------------------------------------------

  ecokit::cat_time(
    "Calculate % coverage of different crosswalks per grid cell")

  perc_cover_maps <- purrr::map_dfr(
    .x = c("synhab", "clc_l1", "clc_l2", "clc_l3", "eunis2019"),
    .f = clc_get_percentage,
    clc_crosswalk = clc_crosswalk, clc_fracs_r = clc_fracs_r,
    path_tif = path_clc_summary_tif, path_rdata = path_clc_summary_rdata)

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Prepare reference grid --- Exclude areas from the study area ----
  # ---------------------------------------------------

  ecokit::cat_time("Reference grid --- Exclude areas from the study area")

  ### |||||||||||||||||||||||||||||||||||||||
  # Islands to exclude
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Islands to exclude", level = 1L)
  exclude_islands <- list(
    melilla = c(2926000, 2944000, 1557000, 1574000),
    ceuta = c(3098000, 3219000, 1411000, 1493000),
    atlantic_islands = c(342000, 2419000, 687000, 2990000)) %>%
    purrr::map(
      .f = ~ {
        raster::extent(.x) %>%
          methods::as("SpatialPolygons") %>%
          sf::st_as_sf()
      }
    ) %>%
    dplyr::bind_rows() %>%
    sf::st_set_crs(3035)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Turkey --- boundaries / extent
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Turkey --- boundaries and extent", level = 1L)
  turkey_boundaries <- eu_boundaries_sf$L_01 %>%
    dplyr::filter(CNTR_ID == "TR") %>%
    dplyr::select(-CNTR_ID)

  # extent to exclude some left-over cells in the east of Turkey
  extent_turkey <- raster::extent(6604000, 7482000, 1707000, 2661000) %>%
    methods::as("SpatialPolygons") %>%
    sf::st_as_sf() %>%
    sf::st_set_crs(3035)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Combine areas to be excluded
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Combine areas to be excluded", level = 1L)
  exclude_area <- grid_sf %>%
    dplyr::mutate(
      turkey_boundaries = as.integer(
        !sf::st_intersects(geometry, turkey_boundaries))) %>%
    dplyr::filter(is.na(turkey_boundaries)) %>%
    sf::st_union(extent_turkey) %>%
    sf::st_union(exclude_islands) %>%
    sf::st_union() %>%
    smoothr::fill_holes(units::set_units(100000, km^2)) %>%
    terra::vect() %>%
    terra::rasterize(terra::rast(grid_r)) %>%
    # Suppress the warning: attribute variables are assumed to be spatially
    # constant throughout all geometries
    suppressWarnings()

  rm(grid_r, envir = environment())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Percentage of water habitats per grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Calculate the % of water per grid cell", level = 1L)
  grid_10_land <- perc_cover_maps %>%
    dplyr::filter(name == "perc_cover_clc_l3") %>%
    dplyr::pull(map) %>%
    magrittr::extract2(1) %>%
    terra::unwrap() %>%
    magrittr::extract2(
      c(
        "clc_l3_423", "clc_l3_511", "clc_l3_512",
        "clc_l3_521", "clc_l3_522", "clc_l3_523")) %>%
    # calculate the sum of these classes
    sum() %>%
    stats::setNames("grid_10_land") %>%
    ecokit::set_raster_crs(crs = "epsg:3035")

  ## ||||||||||||||||||||||||||||||||||||||||
  # Reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Reference grid", level = 1L)

  ecokit::cat_time("Reference grid --- land only", level = 2L)
  grid_10_land[exclude_area == 1] <- NA
  grid_10_land[grid_10_land > (100 - min_land_percent)] <- NA
  grid_10_land[!is.na(grid_10_land)] <- 1

  ecokit::cat_time("Reference grid --- cropped", level = 2L)
  grid_10_land_crop <- terra::trim(grid_10_land) %>%
    stats::setNames("grid_10_land_crop") %>%
    ecokit::set_raster_crs(crs = "epsg:3035")

  ecokit::cat_time("Reference grid --- sf object", level = 2L)
  grid_10_land_sf <- terra::as.points(grid_10_land) %>%
    sf::st_as_sf() %>%
    sf::st_join(x = grid_sf, y = .) %>%
    dplyr::filter(grid_10_land == 1) %>%
    dplyr::select(-grid_10_land)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save reference grid --- RData", level = 1L)
  grid_10_land <- terra::wrap(grid_10_land)
  grid_10_land_crop <- terra::wrap(grid_10_land_crop)

  save(grid_10_land, file = fs::path(path_grid, "grid_10_land.RData"))
  save(
    grid_10_land_crop,
    file = fs::path(path_grid, "grid_10_land_crop.RData"))
  save(
    grid_10_land_sf, file = fs::path(path_grid, "grid_10_land_sf.RData"))

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save calculated % coverage
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save calculated % coverage", level = 1L)
  clc_fracs_r <- terra::wrap(clc_fracs_r)
  save(
    clc_fracs_r,
    file = fs::path(path_clc_summary_rdata, "clc_fracs_r.RData"))

  rm(clc_fracs_r, envir = environment())
  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid - tif
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save reference grid --- tif", level = 1L)
  terra::writeRaster(
    terra::unwrap(grid_10_land), overwrite = TRUE,
    filename = fs::path(path_grid, "grid_10_land.tif"))

  terra::writeRaster(
    terra::unwrap(grid_10_land_crop), overwrite = TRUE,
    filename = fs::path(path_grid, "grid_10_land_crop.tif"))

  # # ..................................................................... ###

  ## ||||||||||||||||||||||||||||||||||||||||
  # Add country name to grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  # country boundaries
  eu_boundaries_countries <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_01") %>%
    dplyr::select(country = "NAME_ENGL") %>%
    dplyr::filter(!(country %in% c("Russian Federation", "Belarus", "Ukraine")))

  # Reference grid as points, with spatially matched country names
  grid_country <- terra::unwrap(grid_10_land_crop) %>%
    terra::as.points() %>%
    sf::st_as_sf() %>%
    sf::st_join(eu_boundaries_countries) %>%
    dplyr::select(-"grid_10_land_crop")

  # Get the nearest country names for some grid cells that their centroid do not
  # overlap with country boundaries
  grid_country_near <- dplyr::filter(grid_country, is.na(country)) %>%
    dplyr::select(-"country") %>%
    sf::st_join(y = eu_boundaries_countries, join = sf::st_nearest_feature) %>%
    dplyr::rename(country_2 = "country")

  # Merge data and add a new column representing whether the country information
  # was retrieved by spatial joining of the grid centroid and country boundaries
  # or estimated as the nearest country
  grid_country <- sf::st_join(x = grid_country, y = grid_country_near) %>%
    sf::st_join(grid_10_land_sf) %>%
    dplyr::mutate(
      Nearest = is.na(country),
      country = dplyr::coalesce(country, country_2)) %>%
    dplyr::select(-"country_2")

  ecokit::save_as(
    object = grid_country, object_name = "grid_10_land_sf_country",
    out_path = fs::path(path_grid, "grid_10_land_sf_country.RData"))

  rm(
    exclude_area, turkey_boundaries, grid_country, grid_country_near,
    envir = environment())

  invisible(gc())

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Crop % coverage results ----
  # ---------------------------------------------------

  ecokit::cat_time("Crop % coverage results")

  perc_cover_maps <- dplyr::mutate(
    perc_cover_maps,
    map_crop = purrr::map2(
      .x = name,
      .y = map,
      .f = ~ {
        clc_type <- stringr::str_remove(.x, "perc_cover_")
        ecokit::cat_time(clc_type, level = 2L)

        map_cropped <- terra::crop(
          terra::unwrap(.y), terra::unwrap(grid_10_land_crop)) %>%
          terra::mask(terra::unwrap(grid_10_land_crop)) %>%
          ecokit::set_raster_crs(crs = "epsg:3035")

        terra::writeRaster(
          x = map_cropped, overwrite = TRUE,
          filename = fs::path(
            path_clc_summary_tif_crop, paste0("perc_cover_", names(map_cropped),
                                              ".tif")))

        out_obj_name <- paste0("perc_cover_", clc_type, "_crop")
        ecokit::save_as(
          object = terra::wrap(map_cropped), object_name = out_obj_name,
          out_path = fs::path(
            path_clc_summary_rdata, paste0(out_obj_name, ".RData")))

        return(terra::wrap(map_cropped))
      }
    ))

  ecokit::save_as(
    object = perc_cover_maps, object_name = "perc_cover_maps",
    out_path = fs::path(path_clc_summary_rdata, "perc_cover_maps.RData"))

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Majority per grid cell ----
  # ---------------------------------------------------

  ecokit::cat_time("Identify major CLC class per per grid cell")

  ## ||||||||||||||||||||||||||||||||||||||||
  # Processing using `exactextractr::exact_extract`
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Processing using `exactextractr::exact_extract`", level = 1L)

  # Ensure that the projection of x and y parameters of
  # `exactextractr::exact_extract` suppress warning: Polygons transformed to
  # raster CRS (EPSG:3035) https://github.com/isciences/exactextractr/issues/103

  clc_majority <- sf::st_transform(grid_10_land_sf, sf::st_crs(clc_rast)) %>%
    dplyr::mutate(
      majority = exactextractr::exact_extract(
        x = clc_rast, y = ., fun = "majority",
        default_value = 44, progress = FALSE)) %>%
    dplyr::left_join(clc_crosswalk, by = dplyr::join_by(majority == value)) %>%
    sf::st_transform(crs = 3035) %>%
    dplyr::mutate(majority = factor(majority)) %>%
    tibble::tibble() %>%
    sf::st_as_sf()

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save majority results", level = 1L)
  save(
    clc_majority,
    file = fs::path(path_clc_summary_rdata, "clc_majority.RData"))

  rm(clc_rast, grid_sf, grid_10_land_sf, envir = environment())
  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # post-processing majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Post-processing majority results", level = 1L)

  majority_maps <- purrr::map_dfr(
    .x = c("synhab", "clc_l1", "clc_l2", "clc_l3", "eunis2019"),
    .f = clc_get_majority,
    clc_majority = clc_majority, path_tif = path_clc_summary_tif,
    path_tif_crop = path_clc_summary_tif_crop,
    path_rdata = path_clc_summary_rdata,
    grid_10_land = grid_10_land, grid_10_land_crop = grid_10_land_crop)

  rm(
    clc_majority, grid_10_land, grid_10_land_crop, majority_maps,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # ---------------------------------------------------
  # Plotting ----
  # ---------------------------------------------------

  if (plot_clc) {

    ecokit::cat_time("Plotting")

    if (n_cores_plot == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores_plot, level = 1L, strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    eu_map <- eu_boundaries_sf$L_03

    # packages to be loaded in parallel
    pkg_to_export <- ecokit::load_packages_future(
      packages = c(
        "cowplot", "dplyr", "ecokit", "fs", "ggplot2", "ggtext", "grid",
        "magrittr", "paletteer", "purrr", "ragg", "scales", "stringi",
        "stringr", "terra", "tibble", "tidyselect", "tidyterra", "viridis"),
      strategy = strategy)

    plot_list <- future.apply::future_lapply(
      X = c(
        "perc_cover_synhab", "perc_cover_clc_l1", "perc_cover_clc_l2",
        "perc_cover_clc_l3", "perc_cover_eunis2019"),
      FUN = function(x) {
        clc_plot(
          clc_name = x, clc_map = perc_cover_maps, eu_map = eu_map,
          crosswalk = clc_crosswalk, path_jpeg = path_clc_summary_jpeg)
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c(
        "clc_plot", "perc_cover_maps", "eu_map", "clc_crosswalk",
        "path_clc_summary_jpeg"))

    rm(plot_list, envir = environment())
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nProcessing CLC data was finished in ")

  return(invisible(NULL))
}

# --------------------------------------------------- ####
# --------------------------------------------------- ####

## |------------------------------------------------------------------------| #
# clc_get_percentage ------
## |------------------------------------------------------------------------| #

#' clc_get_percentage
#'
#' This function calculate % coverage of different crosswalks per grid cell.
#' The function outputs raster and `RData` files containing the percentage cover
#' information.
#' @param clc_type Character. The crosswalk type to be processed. This has to
#'   be one of `synhab`, `clc_l1`, `clc_l2`, `clc_l3`, and `eunis2019`.
#' @param clc_crosswalk `data.frame`. A data frame containing the crosswalk data
#'   between values and habitat types.
#' @param clc_fracs_r A list of rasters for percent coverage results.
#' @param path_tif,path_rdata Character. Path where the output `TIFF` and
#'   `RData` files, respectively, will be saved.
#' @noRd
#' @keywords internal
#' @return A tibble containing the processed raster object.
#' @author Ahmed El-Gabbas

clc_get_percentage <- function(
    clc_type, clc_crosswalk, clc_fracs_r, path_tif, path_rdata) {
  # # ..................................................................... ###

  if (is.null(clc_type) || is.null(clc_crosswalk) || is.null(clc_fracs_r) ||
      is.null(path_tif) || is.null(path_rdata)) {
    ecokit::stop_ctx(
      "None of the input parameters can be empty",
      clc_type = clc_type, clc_crosswalk = clc_crosswalk,
      clc_fracs_r = clc_fracs_r, path_tif = path_tif, path_rdata = path_rdata,
      include_backtrace = TRUE)
  }

  if (!(
    clc_type %in% c("synhab", "clc_l1", "clc_l2", "clc_l3", "eunis2019"))) {
    ecokit::stop_ctx(
      paste0(
        "clc_type has to be one of synhab, clc_l1, clc_l2, ",
        "clc_l3, and eunis2019"),
      clc_type = clc_type, include_backtrace = TRUE)
  }
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  fracs <- cw_class <- hab_percent <- NULL

  # # ..................................................................... ###

  ecokit::cat_time(clc_type, level = 1L)
  out_obj_name <- paste0("perc_cover_", clc_type)

  map <- dplyr::select(
    clc_crosswalk, "value", tidyselect::all_of(clc_type)) %>%
    stats::setNames(c("fracs", "cw_class")) %>%
    tidyr::nest(fracs = -cw_class) %>%
    dplyr::slice(gtools::mixedorder(cw_class)) %>%
    dplyr::mutate(
      fracs = purrr::map(.x = fracs, .f = ~ as.vector(unlist(.x))),
      hab_percent = purrr::map2(
        .x = fracs, .y = cw_class,
        .f = ~ {
          r_name <- paste0(clc_type, "_", .y) %>%
            stringr::str_replace_all("\\.", "") %>%
            stringr::str_trim()

          clc_fracs_r[[paste0("frac_", .x)]] %>%
            sum() %>%
            stats::setNames(r_name)
        }
      )
    ) %>%
    dplyr::pull(hab_percent) %>%
    terra::rast()

  terra::writeRaster(
    map, overwrite = TRUE,
    filename = fs::path(path_tif, paste0("perc_cover_", names(map), ".tif")))

  ecokit::save_as(
    object = terra::wrap(map), object_name = out_obj_name,
    out_path = fs::path(path_rdata, paste0(out_obj_name, ".RData")))

  tibble::tibble(name = out_obj_name, map = list(terra::wrap(map)))
}

# --------------------------------------------------- ####
# --------------------------------------------------- ####

## |------------------------------------------------------------------------| #
# clc_get_majority ------
## |------------------------------------------------------------------------| #

#' clc_get_majority
#'
#' This function processes CLC majority data for a specified type and generates
#' corresponding raster and vector data. It includes functionality to crop and
#' mask the data to a specific region.
#' @param clc_type Character. The type of data to process. Must be one of
#'   "synhab", "clc_l1", "clc_l2", "clc_l3", or "eunis2019".
#' @param clc_majority A data frame containing the majority class data.
#' @param path_tif Character. The path where the output TIFF files will
#'   be saved.
#' @param path_tif_crop Character. The path where the cropped output
#'   TIFF files will be saved.
#' @param path_rdata Character. he path where the output RData files
#'   will be saved.
#' @param grid_10_land A `SpatRaster` object representing the land grid for
#'   rasterization.
#' @param grid_10_land_crop A `SpatRaster` object representing the cropped land
#'   grid for masking.
#' @return A tibble containing the type and the processed raster and cropped
#'   raster objects.
#' @noRd
#' @keywords internal
#' @author Ahmed El-Gabbas

clc_get_majority <- function(
    clc_type, clc_majority, path_tif, path_tif_crop, path_rdata,
    grid_10_land, grid_10_land_crop) {

  # # ..................................................................... ###

  if (is.null(clc_type) || is.null(clc_majority) || is.null(path_tif) ||
      is.null(path_tif_crop) || is.null(path_rdata) || is.null(grid_10_land) ||
      is.null(grid_10_land_crop)) {
    ecokit::stop_ctx(
      "None of the input parameters can be empty",
      clc_type = clc_type, clc_majority = clc_majority,
      path_tif = path_tif, path_tif_crop = path_tif_crop,
      path_rdata = path_rdata, grid_10_land = grid_10_land,
      grid_10_land_crop = grid_10_land_crop, include_backtrace = TRUE)
  }

  if (!(
    clc_type %in% c("synhab", "clc_l1", "clc_l2", "clc_l3", "eunis2019"))) {
    ecokit::stop_ctx(
      paste0(
        "clc_type has to be one of synhab, clc_l1, clc_l2, ",
        "clc_l3, and eunis2019"),
      clc_type = clc_type, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  cw_label <- cw_class <- NULL

  # # ..................................................................... ###

  ecokit::cat_time(clc_type, level = 2L)
  out_obj_name <- paste0("majority_", clc_type)
  out_obj_name_crop <- paste0(out_obj_name, "_crop")

  na_classes <- c(
    "Marine habitats", "5_Water bodies",
    "5_2_Marine waters", "5_2_3_Sea and ocean", "A_Marine_Marine habitats")

  map <- clc_majority %>%
    dplyr::select(
      tidyselect::starts_with(clc_type), -tidyselect::ends_with("_desc")) %>%
    stats::setNames(c("cw_class", "cw_label", "geometry")) %>%
    dplyr::mutate(
      cw_class = paste0(cw_class, "_", cw_label),
      cw_class = stringr::str_replace_all(cw_class, "\\.|\\._", "_"),
      cw_class = stringr::str_replace_all(cw_class, "__", "_"),
      cw_class = stringr::str_replace_all(cw_class, "NA_NA", NA_character_),
      cw_label = NULL) %>%
    dplyr::filter(!is.na(cw_class), !cw_class %in% na_classes) %>%
    terra::rasterize(terra::unwrap(grid_10_land), field = "cw_class") %>%
    stats::setNames(out_obj_name)

  # Remap raster cell values
  map_levels <- magrittr::extract2(terra::levels(map), 1) %>%
    ecokit::arrange_alphanum("cw_class") %>%
    dplyr::mutate(cw_id = seq_len(dplyr::n()))
  # Create a named vector to map old IDs to new sequential IDs
  id_map <- stats::setNames(map_levels$cw_id, map_levels$ID)
  # Create a lookup for all possible raster values
  lookup <- rep(NA_integer_, max(map_levels$ID) + 1)
  lookup[as.integer(names(id_map)) + 1] <- id_map
  map <- terra::classify(map, rcl = cbind(map_levels$ID, map_levels$cw_id))
  # Now set levels with sequential IDs and your desired order
  levels(map)[[1]] <- map_levels[, c("cw_id", "cw_class")]

  terra::writeRaster(
    x = map, overwrite = TRUE,
    filename = fs::path(path_tif, paste0(out_obj_name, ".tif")))
  terra::levels(map) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = cw_class) %>%
    foreign::write.dbf(
      file = fs::path(path_tif, paste0(out_obj_name, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)
  ecokit::save_as(
    object = terra::wrap(map), object_name = out_obj_name,
    out_path = fs::path(path_rdata, paste0(out_obj_name, ".RData")))

  # cropping
  map_crop <- terra::crop(x = map, y = terra::unwrap(grid_10_land_crop))

  terra::writeRaster(
    x = map_crop, overwrite = TRUE,
    filename = fs::path(path_tif_crop, paste0(out_obj_name_crop, ".tif")))
  terra::levels(map_crop) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = cw_class) %>%
    foreign::write.dbf(
      file = fs::path(
        path_tif_crop, paste0(out_obj_name_crop, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)
  ecokit::save_as(
    object = terra::wrap(map_crop), object_name = out_obj_name_crop,
    out_path = fs::path(path_rdata, paste0(out_obj_name_crop, ".RData")))

  tibble::tibble(
    clc_type = clc_type,
    map = list(terra::wrap(map)),
    map_crop = list(terra::wrap(map_crop)))
}
