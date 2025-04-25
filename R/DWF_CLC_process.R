## |------------------------------------------------------------------------| #
# CLC_process ------
## |------------------------------------------------------------------------| #

#' Process Corine Land Cover (CLC) data for the `IAS-pDT`
#'
#' Processes [Corine Land Cover
#' (CLC)](https://land.copernicus.eu/pan-european/corine-land-cover/clc2018)
#' data for the Invasive Alien Species prototype Digital Twin (`IAS-pDT`).
#' Calculates percentage coverage and most common classes per grid cell at three
#' CLC levels, plus `EUNIS_19` and `SynHab` habitat types. Prepares a reference
#' grid and optionally generates percentage coverage maps as JPEG.
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param min_land_percent Numeric. Minimum land percentage per grid cell for
#'   the reference grid. Default: `15`.
#' @param plot_CLC Logical. If `TRUE`, plots percentage coverage for CLC levels
#'   and habitat types. Default: `TRUE`.
#' @return Returns `invisible(NULL)`; saves processed data and optional plots to
#'   disk.
#' @name CLC_process
#'
#' @author Ahmed El-Gabbas
#' @export
#' @references
#' - Data source:
#' <https://land.copernicus.eu/pan-european/corine-land-cover/clc2018>
#' - Data citation:
#' <https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac>

CLC_process <- function(
    env_file = ".env", min_land_percent = 15, plot_CLC = TRUE) {
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid warning while reading CLC data
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SynHab_desc <- CNTR_ID <- geometry <- Name <- Path_CLC <- Path_Grid <-
    Path_Grid_Ref <- km <- Majority <- Path_CLC_tif <- Path_CLC_CW <-
    EU_Bound <- Value <- Country <- Country2 <- NULL

  if (is.null(env_file)) {
    IASDT.R::stop_ctx("env_file can not be empty", env_file = env_file)
  }

  if (!is.numeric(min_land_percent) ||
      !dplyr::between(min_land_percent, 0, 100)) {
    IASDT.R::stop_ctx(
      "`min_land_percent` must be a numeric value between 0 and 100.",
      min_land_percent = min_land_percent)
  }

  if (!file.exists(env_file)) {
    IASDT.R::stop_ctx(
      "Path to environment variables was not found", env_file = env_file)
  }

  # # ..................................................................... ###

  ## ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Loading data ------
  ## ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time("Loading data")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Environment variables
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", FALSE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "Path_CLC", "DP_R_CLC_processed", FALSE, FALSE,
    "Path_CLC_tif", "DP_R_CLC_tif", FALSE, TRUE,
    "Path_CLC_CW", "DP_R_CLC_crosswalk", FALSE, TRUE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Check files/directories
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Check files/directories", level = 1)
  fs::dir_create(Path_CLC)

  Path_Grid_sf <- IASDT.R::path(Path_Grid_Ref, "Grid_10_sf.RData")
  Path_Grid_rast <- IASDT.R::path(Path_Grid_Ref, "Grid_10_Raster.RData")

  requiredPaths <- c(Path_Grid_sf, Path_CLC_CW, Path_Grid_rast)

  purrr::walk(
    .x = requiredPaths,
    .f = function(path) {
      if (!file.exists(path)) {
        IASDT.R::stop_ctx("Required path does not exist", path = path)
      }
    }
  )

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Cross-walk
  # description of CLC values and custom cross-walks
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Loading cross-walk", level = 1)
  CLC_crossWalk <- readr::read_delim(Path_CLC_CW, show_col_types = FALSE) %>%
    dplyr::select(-SynHab_desc)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # reference grid
  # # ||||||||||||||||||||||||||||||||||||||||||||

  # Path of reference grid
  IASDT.R::cat_time("Loading reference grid", level = 1)

  # sf
  IASDT.R::cat_time("sf", level = 2)

  Grid_sf <- IASDT.R::load_as(Path_Grid_sf) %>%
    magrittr::extract2("Grid_10_sf_s")

  # raster
  IASDT.R::cat_time("raster", level = 2)
  Grid_R <- IASDT.R::load_as(Path_Grid_rast)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # country boundaries
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Loading country boundaries", level = 1)
  EUBound_sf <- IASDT.R::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Creating folders
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Creating folders", level = 1)

  # sub-folders to store tif files (no masking to final reference grid)
  Path_CLC_Summary_Tif <- IASDT.R::path(Path_CLC, "Summary_Tif")

  # sub-folders to store tif files (masked to final reference grid)
  Path_CLC_Summary_Tif_Crop <- IASDT.R::path(Path_CLC, "Summary_Tif_Crop")

  # sub-folders to store output RData files
  Path_CLC_Summary_RData <- IASDT.R::path(Path_CLC, "Summary_RData")

  # Create folders when necessary
  c(
    Path_CLC, Path_CLC_Summary_Tif, Path_CLC_Summary_Tif_Crop,
    Path_CLC_Summary_RData) %>%
    purrr::walk(fs::dir_create)


  if (plot_CLC) {
    # sub-folders to store JPEG files
    Path_CLC_Summary_JPEG <- IASDT.R::path(Path_CLC, "Summary_JPEG")
    Path_CLC_Summary_JPEG_Free <- IASDT.R::path(
      Path_CLC_Summary_JPEG, "FreeLegend")

    fs::dir_create(c(Path_CLC_Summary_JPEG, Path_CLC_Summary_JPEG_Free))
  }

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Calculate fraction of each CLC class at each grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time("Calculate fraction of each CLC class at each grid cell")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Loading CLC tif file
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Loading CLC tif file", level = 1)
  CLC_Rast <- terra::rast(Path_CLC_tif)
  terra::NAflag(CLC_Rast) <- 128

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Calculate fraction for each CLC value
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time(
    "Processing using `exactextractr::exact_extract` function",
    level = 1)

  CLC_Fracs <- Grid_sf %>%
    # Ensure that the projection of x and y parameters of
    # `exactextractr::exact_extract` suppress warning: Polygons transformed to
    # raster CRS (EPSG:3035)
    # https://github.com/isciences/exactextractr/issues/103
    sf::st_transform(sf::st_crs(CLC_Rast)) %>%
    dplyr::mutate(
      exactextractr::exact_extract(
        x = CLC_Rast, y = ., fun = "frac",
        force_df = TRUE, default_value = 44, progress = FALSE)) %>%
    sf::st_transform(crs = 3035)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Save fraction results
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Save fraction results", level = 1)
  save(
    CLC_Fracs,
    file = IASDT.R::path(Path_CLC_Summary_RData, "CLC_Fracs.RData"))

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Convert fractions to raster
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Convert fractions to raster", level = 1)
  # convert to SpatVector objects for faster rasterization
  CLC_Fracs_vect <- terra::vect(CLC_Fracs)

  CLC_fracs_r <- names(CLC_Fracs) %>%
    # Exclude processing CLC values of NODATA "frac_48"
    setdiff(y = c("CellCode", "geometry", "frac_48")) %>%
    purrr::map(
      .f = ~ {
        Map <- terra::rasterize(
          x = CLC_Fracs_vect, y = terra::rast(Grid_R), field = .x) %>%
          magrittr::multiply_by(100) %>%
          stats::setNames(.x) %>%
          IASDT.R::set_raster_CRS()

        # For `frac_44` (CLC 5.2.3 Sea and ocean), the original values do not
        # cover the whole study area, resulting in values of zero in the sea.
        # The following force seas and oceans outside mainland Europe to have
        # 100 value
        if (.x == "frac_44") {
          Map <- terra::trim(Map, value = 0) %>%
            terra::extend(terra::rast(Grid_R), fill = 100)
        }
        return(Map)
      }
    ) %>%
    terra::rast()

  rm(CLC_Fracs_vect, CLC_Fracs, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Calculate % coverage of different cross-walks per grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time(
    "Calculate % coverage of different cross-walks per grid cell")

  PercCovMaps <- purrr::map_dfr(
    .x = c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"),
    .f = CLC_get_percentage,
    CLC_crossWalk = CLC_crossWalk, CLC_fracs_r = CLC_fracs_r,
    path_tif = Path_CLC_Summary_Tif, path_RData = Path_CLC_Summary_RData)

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Prepare reference grid --- Exclude areas from the study area ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time("Reference grid --- Exclude areas from the study area")

  ### |||||||||||||||||||||||||||||||||||||||
  # Islands to exclude
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Islands to exclude", level = 1)
  Exclude_Islands <- list(
    melilla = c(2926000, 2944000, 1557000, 1574000),
    ceuta = c(3098000, 3219000, 1411000, 1493000),
    AtlanticIslands = c(342000, 2419000, 687000, 2990000)) %>%
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

  IASDT.R::cat_time("Turkey --- boundaries / extent", level = 1)
  TR <- EUBound_sf$L_01 %>%
    dplyr::filter(CNTR_ID == "TR") %>%
    dplyr::select(-CNTR_ID)

  # extent to exclude some left-over cells in the east of Turkey
  Extent_TR <- raster::extent(6604000, 7482000, 1707000, 2661000) %>%
    methods::as("SpatialPolygons") %>%
    sf::st_as_sf() %>%
    sf::st_set_crs(3035)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Combine areas to be excluded
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Combine areas to be excluded", level = 1)
  Exclude_Area <- Grid_sf %>%
    dplyr::mutate(TR = as.integer(!sf::st_intersects(geometry, TR))) %>%
    dplyr::filter(is.na(TR)) %>%
    sf::st_union(Extent_TR) %>%
    sf::st_union(Exclude_Islands) %>%
    sf::st_union() %>%
    smoothr::fill_holes(units::set_units(100000, km^2)) %>%
    terra::vect() %>%
    terra::rasterize(terra::rast(Grid_R)) %>%
    # Suppress the warning: attribute variables are assumed to be spatially
    # constant throughout all geometries
    suppressWarnings()

  rm(Grid_R, envir = environment())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Percentage of water habitats per grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Calculate the % of water per grid cell", level = 1)
  Grid_10_Land <- dplyr::filter(PercCovMaps, Name == "PercCov_CLC_L3") %>%
    dplyr::pull(Map) %>%
    magrittr::extract2(1) %>%
    magrittr::extract2(
      c(
        "CLC_L3_423", "CLC_L3_511", "CLC_L3_512",
        "CLC_L3_521", "CLC_L3_522", "CLC_L3_523")) %>%
    # calculate the sum of these classes
    sum() %>%
    stats::setNames("Grid_10_Land") %>%
    IASDT.R::set_raster_CRS()

  ## ||||||||||||||||||||||||||||||||||||||||
  # Reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Reference grid", level = 1)

  IASDT.R::cat_time("Reference grid --- land only", level = 2)
  Grid_10_Land[Exclude_Area == 1] <- NA
  Grid_10_Land[Grid_10_Land > (100 - min_land_percent)] <- NA
  Grid_10_Land[!is.na(Grid_10_Land)] <- 1

  IASDT.R::cat_time("Reference grid --- cropped", level = 2)
  Grid_10_Land_Crop <- terra::trim(Grid_10_Land) %>%
    stats::setNames("Grid_10_Land_Crop") %>%
    IASDT.R::set_raster_CRS()

  IASDT.R::cat_time("Reference grid --- sf object", level = 2)
  Grid_10_Land_sf <- terra::as.points(Grid_10_Land) %>%
    sf::st_as_sf() %>%
    sf::st_join(x = Grid_sf, y = .) %>%
    dplyr::filter(Grid_10_Land == 1) %>%
    dplyr::select(-Grid_10_Land)

  IASDT.R::cat_time("Reference grid --- sf object - cropped", level = 2)
  Grid_10_Land_Crop_sf <- terra::as.points(Grid_10_Land_Crop) %>%
    sf::st_as_sf() %>%
    sf::st_join(x = Grid_sf, y = .) %>%
    dplyr::filter(Grid_10_Land_Crop == 1) %>%
    dplyr::select(-Grid_10_Land_Crop)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Save reference grid --- RData", level = 1)
  Grid_10_Land <- terra::wrap(Grid_10_Land)
  Grid_10_Land_Crop <- terra::wrap(Grid_10_Land_Crop)

  fs::dir_create(Path_Grid)
  save(Grid_10_Land, file = IASDT.R::path(Path_Grid, "Grid_10_Land.RData"))
  save(
    Grid_10_Land_Crop,
    file = IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData"))
  save(
    Grid_10_Land_sf,  file = IASDT.R::path(Path_Grid, "Grid_10_Land_sf.RData"))
  save(
    Grid_10_Land_Crop_sf,
    file = IASDT.R::path(Path_Grid, "Grid_10_Land_Crop_sf.RData"))

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save calculated % coverage
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Save calculated % coverage", level = 1)
  CLC_fracs_r <- terra::wrap(CLC_fracs_r)
  save(
    CLC_fracs_r,
    file = IASDT.R::path(Path_CLC_Summary_RData, "CLC_fracs_r.RData"))

  rm(CLC_fracs_r, envir = environment())
  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid - tif
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Save reference grid --- tif", level = 1)
  terra::writeRaster(
    terra::unwrap(Grid_10_Land), overwrite = TRUE,
    filename = IASDT.R::path(Path_Grid, "Grid_10_Land.tif"))

  terra::writeRaster(
    terra::unwrap(Grid_10_Land_Crop), overwrite = TRUE,
    filename = IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.tif"))

  # # ..................................................................... ###

  ## ||||||||||||||||||||||||||||||||||||||||
  # Add country name to grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  # Country boundaries
  EU_BoundCNT <- IASDT.R::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_01") %>%
    dplyr::select(Country = "NAME_ENGL") %>%
    dplyr::filter(!(Country %in% c("Russian Federation", "Belarus", "Ukraine")))

  # Reference grid as points, with spatially matched Country names
  Grid_CNT <- terra::unwrap(Grid_10_Land_Crop) %>%
    terra::as.points() %>%
    sf::st_as_sf() %>%
    sf::st_join(EU_BoundCNT) %>%
    dplyr::select(-"Grid_10_Land_Crop")

  # Get the nearest country names for some grid cells that their centroid do not
  # overlap with country boundaries
  Grid_CNT_Near <- dplyr::filter(Grid_CNT, is.na(Country)) %>%
    dplyr::select(-"Country") %>%
    sf::st_join(y = EU_BoundCNT, join = sf::st_nearest_feature) %>%
    dplyr::rename(Country2 = "Country")

  # Merge data and add a new column representing whether the country information
  # was retrieved by spatial joining of the grid centroid and country boundaries
  # or estimated as the nearest country
  Grid_CNT <- sf::st_join(x = Grid_CNT, y = Grid_CNT_Near) %>%
    sf::st_join(Grid_10_Land_Crop_sf) %>%
    dplyr::mutate(
      Nearest = is.na(Country),
      Country = dplyr::coalesce(Country, Country2)) %>%
    dplyr::select(-"Country2")

  IASDT.R::save_as(
    object = Grid_CNT, object_name = "Grid_10_Land_Crop_sf_Country",
    out_path = IASDT.R::path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData"))

  rm(
    Grid_10_Land_Crop_sf, Grid_10_Land_sf, Exclude_Area,
    TR, Grid_CNT, Grid_CNT_Near, envir = environment())

  invisible(gc())

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Crop % coverage results ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time("Crop % coverage results")

  PercCovMaps <- dplyr::mutate(
    PercCovMaps,
    Map_Crop = purrr::map2(
      .x = Name, .y = Map,
      .f = ~ {
        CLC_type <- stringr::str_remove(.x, "PercCov_")
        IASDT.R::cat_time(CLC_type, level = 2)

        Map <- terra::crop(.y, terra::unwrap(Grid_10_Land_Crop)) %>%
          terra::mask(terra::unwrap(Grid_10_Land_Crop)) %>%
          IASDT.R::set_raster_CRS()

        terra::writeRaster(
          x = Map, overwrite = TRUE,
          filename = IASDT.R::path(
            Path_CLC_Summary_Tif_Crop,
            paste0("PercCov_", names(Map), ".tif")))

        OutObjName <- paste0("PercCov_", CLC_type, "_Crop")
        IASDT.R::save_as(
          object = terra::wrap(Map), object_name = OutObjName,
          out_path = IASDT.R::path(
            Path_CLC_Summary_RData, paste0(OutObjName, ".RData")))

        return(Map)
      }
    )
  )

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Majority per grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::cat_time("Identify major CLC class per per grid cell")

  ## ||||||||||||||||||||||||||||||||||||||||
  # Processing using `exactextractr::exact_extract`
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time(
    "Processing using `exactextractr::exact_extract`", level = 1)

  CLC_majority <- Grid_sf %>%
    # Ensure that the projection of x and y parameters of
    # `exactextractr::exact_extract` suppress warning: Polygons transformed to
    # raster CRS (EPSG:3035)
    # https://github.com/isciences/exactextractr/issues/103
    sf::st_transform(sf::st_crs(CLC_Rast)) %>%
    dplyr::mutate(
      Majority = exactextractr::exact_extract(
        x = CLC_Rast, y = ., fun = "majority",
        default_value = 44, progress = FALSE)) %>%
    dplyr::left_join(CLC_crossWalk, by = dplyr::join_by(Majority == Value)) %>%
    tibble::tibble() %>%
    sf::st_as_sf() %>%
    sf::st_transform(crs = 3035)

  rm(CLC_Rast, Grid_sf, envir = environment())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Save majority results", level = 1)
  save(
    CLC_majority,
    file = IASDT.R::path(Path_CLC_Summary_RData, "CLC_majority.RData"))

  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # post-processing majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_time("Post-processing majority results", level = 1)

  MajorityMaps <- purrr::map_dfr(
    .x = c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"),
    .f = CLC_get_majority,
    CLC_majority = CLC_majority, path_tif = Path_CLC_Summary_Tif,
    path_tif_crop = Path_CLC_Summary_Tif_Crop,
    path_RData = Path_CLC_Summary_RData,
    Grid_10_Land = Grid_10_Land, Grid_10_Land_Crop = Grid_10_Land_Crop)

  rm(
    CLC_majority, Grid_10_Land, Grid_10_Land_Crop, MajorityMaps,
    envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Plotting ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  if (plot_CLC) {
    IASDT.R::cat_time("Plotting")
    c(
      "PercCov_SynHab", "PercCov_CLC_L1", "PercCov_CLC_L2",
      "PercCov_CLC_L3", "PercCov_EUNIS_2019") %>%
      purrr::walk(
        .f = CLC_plot,
        CLC_map = PercCovMaps, EU_map = EUBound_sf$L_03,
        crosswalk = CLC_crossWalk, path_JPEG = Path_CLC_Summary_JPEG,
        path_JPEG_free = Path_CLC_Summary_JPEG_Free)
  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nProcessing CLC data was finished in ")

  return(invisible(NULL))
}

# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####
# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####

## |------------------------------------------------------------------------| #
# CLC_get_percentage ------
## |------------------------------------------------------------------------| #

#' CLC_get_percentage
#'
#' This function Calculate % coverage of different cross-walks per grid cell.
#' The function outputs raster and `RData` files containing the percentage cover
#' information.
#' @param CLC_type Character. The cross-walk type to be processed. This has to
#'   be one of `SynHab`, `CLC_L1`, `CLC_L2`, `CLC_L3`, and `EUNIS_2019`.
#' @param CLC_crossWalk `data.frame`. A data frame containing the crosswalk data
#'   between values and habitat types.
#' @param CLC_fracs_r A list of rasters for percent coverage results.
#' @param path_tif,path_RData Character. Path where the output `TIFF` and
#'   `RData` files, respectively, will be saved.
#' @noRd
#' @keywords internal
#' @return A tibble containing the processed raster object.
#' @author Ahmed El-Gabbas

CLC_get_percentage <- function(
    CLC_type, CLC_crossWalk, CLC_fracs_r, path_tif, path_RData) {
  # # ..................................................................... ###

  if (is.null(CLC_type) || is.null(CLC_crossWalk) || is.null(CLC_fracs_r) ||
      is.null(path_tif) || is.null(path_RData)) {
    IASDT.R::stop_ctx(
      "None of the input parameters can be empty",
      CLC_type = CLC_type, CLC_crossWalk = CLC_crossWalk,
      CLC_fracs_r = CLC_fracs_r, path_tif = path_tif, path_RData = path_RData)
  }

  if (!(
    CLC_type %in% c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"))) {
    IASDT.R::stop_ctx(
      paste0(
        "CLC_type has to be one of SynHab, CLC_L1, CLC_L2, ",
        "CLC_L3, and EUNIS_2019"),
      CLC_type = CLC_type)
  }
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Fracs <- Class <- HabPerc <- NULL

  # # ..................................................................... ###

  IASDT.R::cat_time(CLC_type, level = 1)
  OutObjName <- paste0("PercCov_", CLC_type)

  Map <- dplyr::select(
    CLC_crossWalk, "Value", tidyselect::all_of(CLC_type)) %>%
    stats::setNames(c("Fracs", "Class")) %>%
    tidyr::nest(Fracs = -Class) %>%
    dplyr::slice(gtools::mixedorder(Class)) %>%
    dplyr::mutate(
      Fracs = purrr::map(.x = Fracs, .f = ~ as.vector(unlist(.x))),
      HabPerc = purrr::map2(
        .x = Fracs, .y = Class,
        .f = ~ {
          RName <- paste0(CLC_type, "_", .y) %>%
            stringr::str_replace_all("\\.", "") %>%
            stringr::str_trim()

          CLC_fracs_r[[paste0("frac_", .x)]] %>%
            sum() %>%
            stats::setNames(RName)
        }
      )
    ) %>%
    dplyr::pull(HabPerc) %>%
    terra::rast()

  terra::writeRaster(
    Map, overwrite = TRUE,
    filename = IASDT.R::path(path_tif, paste0("PercCov_", names(Map), ".tif")))

  IASDT.R::save_as(
    object = terra::wrap(Map), object_name = OutObjName,
    out_path = IASDT.R::path(path_RData, paste0(OutObjName, ".RData")))

  return(tibble::tibble(Name = OutObjName, Map = list(Map)))
}

# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####
# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####

## |------------------------------------------------------------------------| #
# CLC_get_majority ------
## |------------------------------------------------------------------------| #

#' CLC_get_majority
#'
#' This function processes CLC majority data for a specified type and generates
#' corresponding raster and vector data. It includes functionality to crop and
#' mask the data to a specific region.
#' @param CLC_type Character. The type of data to process. Must be one of
#'   "SynHab", "CLC_L1", "CLC_L2", "CLC_L3", or "EUNIS_2019".
#' @param CLC_majority A data frame containing the majority class data.
#' @param path_tif Character. The path where the output TIFF files will
#'   be saved.
#' @param path_tif_crop Character. The path where the cropped output
#'   TIFF files will be saved.
#' @param path_RData Character. he path where the output RData files
#'   will be saved.
#' @param Grid_10_Land A `SpatRaster` object representing the land grid for
#'   rasterization.
#' @param Grid_10_Land_Crop A `SpatRaster` object representing the cropped land
#'   grid for masking.
#' @return A tibble containing the type and the processed raster and cropped
#'   raster objects.
#' @noRd
#' @keywords internal
#' @author Ahmed El-Gabbas

CLC_get_majority <- function(
    CLC_type, CLC_majority, path_tif, path_tif_crop, path_RData,
    Grid_10_Land, Grid_10_Land_Crop) {
  # # ..................................................................... ###

  if (is.null(CLC_type) || is.null(CLC_majority) || is.null(path_tif) ||
      is.null(path_tif_crop) || is.null(path_RData) || is.null(Grid_10_Land) ||
      is.null(Grid_10_Land_Crop)) {
    IASDT.R::stop_ctx(
      "None of the input parameters can be empty",
      CLC_type = CLC_type, CLC_majority = CLC_majority,
      path_tif = path_tif, path_tif_crop = path_tif_crop,
      path_RData = path_RData, Grid_10_Land = Grid_10_Land,
      Grid_10_Land_Crop = Grid_10_Land_Crop)
  }

  if (!(
    CLC_type %in% c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"))) {
    IASDT.R::stop_ctx(
      paste0(
        "CLC_type has to be one of SynHab, CLC_L1, CLC_L2, ",
        "CLC_L3, and EUNIS_2019"),
      CLC_type = CLC_type)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/

  Label <- Class <- ID <- NULL

  # # ..................................................................... ###

  IASDT.R::cat_time(CLC_type, level = 2)
  OutObjName <- paste0("Majority_", CLC_type)
  OutObjName_Cr <- paste0(OutObjName, "_Crop")

  Map <- dplyr::filter(CLC_majority, !is.na(CLC_type)) %>%
    dplyr::select(
      tidyselect::starts_with(CLC_type), -tidyselect::ends_with("_desc")) %>%
    stats::setNames(c("ID", "Label", "geometry")) %>%
    dplyr::mutate(
      ID = paste0(ID, "_", Label),
      ID = stringr::str_replace_all(ID, "\\.|\\._", "_"),
      ID = stringr::str_replace_all(ID, "__", "_"),
      ID = stringr::str_replace_all(ID, "NA_NA+", NA_character_),
      Label = NULL) %>%
    # https://stackoverflow.com/questions/43487773/
    dplyr::rename(!!OutObjName := ID) %>%
    terra::rasterize(terra::unwrap(Grid_10_Land), field = OutObjName)

  MapLevels <- magrittr::extract2(terra::levels(Map), 1)
  NA_Flag <- dplyr::pull(dplyr::filter(MapLevels, is.na(get(OutObjName))), 1)
  terra::NAflag(Map) <- (NA_Flag)
  Map <- terra::droplevels(Map)

  MapLevelsNew <- MapLevels %>%
    dplyr::slice(gtools::mixedorder(.[, 2])) %>%
    dplyr::mutate(ID = seq_len(dplyr::n()))

  MapLevelsM <- dplyr::left_join(
    x = MapLevels, y = MapLevelsNew, by = names(MapLevels)[2]) %>%
    dplyr::select(
      tidyselect::all_of(OutObjName), tidyselect::everything())

  Map <- terra::classify(Map, MapLevelsM[, -1])
  levels(Map) <- list(MapLevelsNew)

  NAClasses <- c(
    "Marine_Marine habitats", "5_Water bodies",
    "5_2_Marine waters", "5_2_3_Sea and ocean", "A_Marine habitats")

  VV <- stats::setNames(MapLevelsNew, c("ID", "Class")) %>%
    dplyr::filter(Class %in% NAClasses) %>%
    dplyr::pull(ID)

  Map <- IASDT.R::set_raster_CRS(terra::classify(Map, cbind(NA, VV)))
  levels(Map) <- list(MapLevelsNew)

  terra::writeRaster(
    x = Map, overwrite = TRUE,
    filename = IASDT.R::path(path_tif, paste0(OutObjName, ".tif")))

  terra::levels(Map) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = IASDT.R::path(path_tif, paste0(OutObjName, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)

  IASDT.R::save_as(
    object = terra::wrap(Map), object_name = OutObjName,
    out_path = IASDT.R::path(path_RData, paste0(OutObjName, ".RData")))

  # CROPPING
  Map_Cr <- terra::crop(x = Map, y = terra::unwrap(Grid_10_Land_Crop)) %>%
    terra::mask(mask = terra::unwrap(Grid_10_Land_Crop)) %>%
    IASDT.R::set_raster_CRS()

  terra::writeRaster(
    x = Map_Cr, overwrite = TRUE,
    filename = IASDT.R::path(path_tif_crop, paste0(OutObjName_Cr, ".tif")))

  terra::levels(Map_Cr) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = IASDT.R::path(
        path_tif_crop, paste0(OutObjName_Cr, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)

  IASDT.R::save_as(
    object = terra::wrap(Map_Cr), object_name = OutObjName_Cr,
    out_path = IASDT.R::path(path_RData, paste0(OutObjName_Cr, ".RData")))

  return(
    tibble::tibble(Type = CLC_type, Map = list(Map), Map_Cr = list(Map_Cr)))
}
