## |------------------------------------------------------------------------| #
# CLC_Process ------
## |------------------------------------------------------------------------| #

#' Process Corine Land Cover (CLC) data
#'
#' Processes Corine Land Cover (CLC) data for environmental modeling purposes.
#' The function calculates the percentage coverage of CLC per grid cell. It
#' estimate the percent coverage at the 3 CLC levels, and cross-walk for
#' `EUNIS_19` and `SynHab` habitat types. Similarly, the function estimates the
#' most common class per grid cell (3 CLC levels, EUNIS, and SynHab habitat
#' types) and prepares reference grid for models. The function optionally plots
#' % coverage maps.
#' @param EnvFile String specifying the path to the environment file containing
#'   necessary paths and configurations.
#' @param FromHPC Logical indicating whether the processing is being done on a
#'   High-Performance Computing (HPC) environment. Defaults to `TRUE`.
#' @param MinLandPerc A numeric value indicating the minimum percentage of land
#'   per grid cell to be used in the reference grid cell. Defaults to 15.
#' @param PlotCLC Logical indicating whether to plot the percentage coverage of
#'   different levels of CLC and custom habitat types. Defaults to `TRUE`.
#' @return The function does not return a value but produces side effects such
#'   as saving processed data and plots to specified directories.
#' @name CLC_Process
#' @author Ahmed El-Gabbas
#' @export
#' @details
#' - [Data source](https://land.copernicus.eu/pan-european/corine-land-cover/clc2018)
#' - Citation](https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac)
#'
#' The function reads the following environment variable:
#'    - **`DP_R_Grid`** (if `FromHPC = TRUE`) or
#'    **`DP_R_Grid_Local`** (if `FromHPC = FALSE`): Path for saving the processed CLC data.
#'    - **`DP_R_Grid_Ref`** (if `FromHPC = TRUE`) or
#'    **`DP_R_Grid_Ref_Local`** (if `FromHPC = FALSE`). The function reads
#' the content of `Grid_10_sf.RData` and `Grid_10_Raster.RData` files from this
#' path.
#'    - **`DP_R_EUBound_sf`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_EUBound_sf_Local`** (if `FromHPC` = `FALSE`): path for the
#' `RData` file containing the country boundaries (`sf` object).
#'   - **`DP_R_CLC`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_CLC_Local`** (if `FromHPC` = `FALSE`): directory where the outputs of CLC data processing are processed.
#'   - **`DP_R_CLC_CW`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_CLC_CW_Local`** (if `FromHPC` = `FALSE`): path for the `crossWalk.txt` file containing custom cross-walk between CLC values and their corresponding values for three levels of CLC and EUNIS_19 & SynHab habitat types.
#'   - **`DP_R_CLC_tif`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_CLC_tif_Local`** (if `FromHPC` = `FALSE`): path for the input CLC `.tif` file.


CLC_Process <- function(
    EnvFile = ".env", FromHPC = TRUE, MinLandPerc = 15, PlotCLC = TRUE) {

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid warning while reading CLC data
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SynHab_desc <- CNTR_ID <- geometry <- Name <- Path_CLC <- Path_Grid <-
    Path_Grid_Ref <- km <- Majority <- Path_CLC_tif <- Path_CLC_CW <-
    EU_Bound <- Value <- Country <- Country2 <- NULL

  if (is.null(EnvFile)) {
    stop("EnvFile can not be empty", .call = FALSE)
  }

  if (!is.numeric(MinLandPerc) || !dplyr::between(MinLandPerc, 0, 100)) {
    stop("MinLandPerc must be a numeric value between 0 and 100.", .call = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(paste0("Path for environment variables (`EnvFile`): ",
                EnvFile, " was not found"), .call = FALSE)
  }

  ## ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Loading data ------
  ## ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk("Loading data")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Loading environment variables
  # # ||||||||||||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("Loading and checking environment variables", Level = 1)
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC", FALSE, FALSE,
      "Path_CLC_tif", "DP_R_CLC_tif", FALSE, TRUE,
      "Path_CLC_CW", "DP_R_CLC_CW", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", FALSE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC_Local", FALSE, FALSE,
      "Path_CLC_tif", "DP_R_CLC_tif_Local", FALSE, TRUE,
      "Path_CLC_CW", "DP_R_CLC_CW_Local", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Check files/directories
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check files/directories", Level = 1)
  fs::dir_create(Path_CLC)

  Path_Grid_sf <- file.path(Path_Grid_Ref, "Grid_10_sf.RData")
  Path_Grid_rast <- file.path(Path_Grid_Ref, "Grid_10_Raster.RData")

  requiredPaths <- c(Path_Grid_sf, Path_CLC_CW, Path_Grid_rast)

  purrr::walk(
    .x = requiredPaths,
    .f = function(path) {
      if (!file.exists(path)) {
        stop(paste0("Required path does not exist: ", path), .call = FALSE)
      }
    })

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Cross-walk
  # description of CLC values and custom cross-walks
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading cross-walk", Level = 1)
  CLC_CrossWalk <- readr::read_delim(Path_CLC_CW, show_col_types = FALSE) %>%
    dplyr::select(-SynHab_desc)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # reference grid
  # # ||||||||||||||||||||||||||||||||||||||||||||

  # Path of reference grid
  IASDT.R::CatTime("Loading reference grid", Level = 1)

  # sf
  IASDT.R::CatTime("sf", Level = 2)

  Grid_sf <- IASDT.R::LoadAs(Path_Grid_sf) %>%
    magrittr::extract2("Grid_10_sf_s")

  # raster
  IASDT.R::CatTime("raster", Level = 2)
  Grid_R <- IASDT.R::LoadAs(Path_Grid_rast)

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # country boundaries
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading country boundaries", Level = 1)
  EUBound_sf <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Creating folders
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Creating folders", Level = 1)

  # sub-folders to store tif files (no masking to final reference grid)
  Path_CLC_Summary_Tif <- file.path(Path_CLC, "Summary_Tif")

  # sub-folders to store tif files (masked to final reference grid)
  Path_CLC_Summary_Tif_Crop <- file.path(Path_CLC, "Summary_Tif_Crop")

  # sub-folders to store output RData files
  Path_CLC_Summary_RData <- file.path(Path_CLC, "Summary_RData")

  # Create folders when necessary
  c(Path_CLC, Path_CLC_Summary_Tif, Path_CLC_Summary_Tif_Crop,
    Path_CLC_Summary_RData) %>%
    purrr::walk(fs::dir_create)


  if (PlotCLC) {
    # sub-folders to store JPEG files
    Path_CLC_Summary_JPEG <- file.path(Path_CLC, "Summary_JPEG")
    Path_CLC_Summary_JPEG_Free <- file.path(
      Path_CLC_Summary_JPEG, "FreeLegend")

    c(Path_CLC_Summary_JPEG, Path_CLC_Summary_JPEG_Free) %>%
      purrr::walk(fs::dir_create)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Calculate fraction of each CLC class at each grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk(
    Message = "Calculate fraction of each CLC class at each grid cell")

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Loading CLC tif file
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading CLC tif file", Level = 1)
  CLC_Rast <- terra::rast(Path_CLC_tif)
  terra::NAflag(CLC_Rast) <- 128

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Calculate fraction for each CLC value
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime(
    "Processing using exactextractr::exact_extract function", Level = 1)

  CLC_Fracs <- Grid_sf %>%
    # Ensure that the projection of x and y parameters of exactextractr::exact_extract
    # suppress warning: Polygons transformed to raster CRS (EPSG:3035)
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

  IASDT.R::CatTime("Save fraction results", Level = 1)
  save(CLC_Fracs,
       file = file.path(Path_CLC_Summary_RData, "CLC_Fracs.RData"))

  # # ||||||||||||||||||||||||||||||||||||||||||||
  # Convert fractions to raster
  # # ||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Convert fractions to raster", Level = 1)
  # convert to SpatVector objects for faster rasterization
  CLC_Fracs_vect <- terra::vect(CLC_Fracs)

  CLC_FracsR <- names(CLC_Fracs) %>%
    # Exclude processing CLC values of NODATA "frac_48"
    setdiff(y = c("CellCode", "geometry", "frac_48")) %>%
    purrr::map(
      .f = ~{
        Map <- terra::rasterize(
          x = CLC_Fracs_vect, y = terra::rast(Grid_R), field = .x) %>%
          magrittr::multiply_by(100) %>%
          stats::setNames(.x) %>%
          IASDT.R::setRastCRS()

        # For `frac_44` (CLC 5.2.3 Sea and ocean), the original values do not
        # cover the whole study area, resulting in values of zero in the sea.
        # The following force seas and oceans outside mainland Europe to have
        # 100 value
        if (.x == "frac_44") {
          Map <- terra::trim(Map, value = 0) %>%
            terra::extend(terra::rast(Grid_R), fill = 100)
        }
        return(Map)
      }) %>%
    terra::rast()

  rm(CLC_Fracs_vect, CLC_Fracs)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Calculate % coverage of different cross-walks per grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk(
    Message = "Calculate % coverage of different cross-walks per grid cell")

  PercCovMaps <- purrr::map_dfr(
    .x = c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"),
    .f = CLC_GetPerc,
    CLC_CrossWalk = CLC_CrossWalk, CLC_FracsR = CLC_FracsR,
    Path_Tif = Path_CLC_Summary_Tif, Path_RData = Path_CLC_Summary_RData)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Prepare reference grid --- Exclude areas from the study area ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk("Reference grid --- Exclude areas from the study area")

  ### |||||||||||||||||||||||||||||||||||||||
  # Islands to exclude
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Islands to exclude", Level = 1)
  Exclude_Islands <- list(
    melilla = c(2926000, 2944000, 1557000, 1574000),
    ceuta = c(3098000, 3219000, 1411000, 1493000),
    AtlanticIslands = c(342000, 2419000, 687000, 2990000)) %>%
    purrr::map(
      .f = ~{
        raster::extent(.x) %>%
          methods::as("SpatialPolygons") %>%
          sf::st_as_sf()
      }) %>%
    dplyr::bind_rows() %>%
    sf::st_set_crs(3035)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Turkey --- boundaries
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Turkey --- boundaries", Level = 1)
  TR <- EUBound_sf$L_01 %>%
    dplyr::filter(CNTR_ID == "TR") %>%
    dplyr::select(-CNTR_ID)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Turkey --- extent to exclude some left-over cells in the east of Turkey
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Turkey --- extent", Level = 1)
  Extent_TR <- raster::extent(6604000, 7482000, 1707000, 2661000) %>%
    methods::as("SpatialPolygons") %>%
    sf::st_as_sf() %>%
    sf::st_set_crs(3035)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Combine areas to be excluded
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Combine areas to be excluded", Level = 1)
  Exclude_Area <- Grid_sf %>%
    dplyr::mutate(TR = as.integer(!sf::st_intersects(geometry, TR))) %>%
    dplyr::filter(is.na(TR)) %>%
    sf::st_union(Extent_TR) %>%
    sf::st_union(Exclude_Islands) %>%
    sf::st_union() %>%
    smoothr::fill_holes(units::set_units(100000, km^2)) %>%
    terra::vect() %>%
    terra::rasterize(terra::rast(Grid_R)) %>%
    #   # Suppress the warning: attribute variables are assumed to be spatially
    #   # constant throughout all geometries
    suppressWarnings()

  rm(Grid_R)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Percentage of water habitats per grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Calculate the % of water per grid cell", Level = 1)
  Grid_10_Land <- dplyr::filter(PercCovMaps, Name == "PercCov_CLC_L3") %>%
    dplyr::pull(Map) %>%
    magrittr::extract2(1) %>%
    magrittr::extract2(
      c("CLC_L3_423", "CLC_L3_511", "CLC_L3_512",
        "CLC_L3_521", "CLC_L3_522", "CLC_L3_523")) %>%
    # calculate the sum of these classes
    sum() %>%
    stats::setNames("Grid_10_Land") %>%
    IASDT.R::setRastCRS()

  ## ||||||||||||||||||||||||||||||||||||||||
  # Reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Reference grid", Level = 1)

  IASDT.R::CatTime("Reference grid --- land only", Level = 2)
  Grid_10_Land[Exclude_Area == 1] <- NA
  Grid_10_Land[Grid_10_Land > (100 - MinLandPerc)] <- NA
  Grid_10_Land[!is.na(Grid_10_Land)] <- 1

  IASDT.R::CatTime("Reference grid --- cropped", Level = 2)
  Grid_10_Land_Crop <- terra::trim(Grid_10_Land) %>%
    stats::setNames("Grid_10_Land_Crop") %>%
    IASDT.R::setRastCRS()

  IASDT.R::CatTime("Reference grid --- sf object", Level = 2)
  Grid_10_Land_sf <- terra::as.points(Grid_10_Land) %>%
    sf::st_as_sf() %>%
    sf::st_join(x = Grid_sf, y = .) %>%
    dplyr::filter(Grid_10_Land == 1) %>%
    dplyr::select(-Grid_10_Land)

  IASDT.R::CatTime("Reference grid --- sf object - cropped", Level = 2)
  Grid_10_Land_Crop_sf <- terra::as.points(Grid_10_Land_Crop) %>%
    sf::st_as_sf() %>%
    sf::st_join(x = Grid_sf, y = .) %>%
    dplyr::filter(Grid_10_Land_Crop == 1) %>%
    dplyr::select(-Grid_10_Land_Crop)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save reference grid --- RData", Level = 1)
  Grid_10_Land <- terra::wrap(Grid_10_Land)
  Grid_10_Land_Crop <- terra::wrap(Grid_10_Land_Crop)

  fs::dir_create(Path_Grid)
  save(Grid_10_Land, file = file.path(Path_Grid, "Grid_10_Land.RData"))
  save(Grid_10_Land_Crop,
       file = file.path(Path_Grid, "Grid_10_Land_Crop.RData"))
  save(Grid_10_Land_sf, file = file.path(Path_Grid, "Grid_10_Land_sf.RData"))
  save(Grid_10_Land_Crop_sf,
       file = file.path(Path_Grid, "Grid_10_Land_Crop_sf.RData"))

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save calculated % coverage
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save calculated % coverage", Level = 1)
  CLC_FracsR <- terra::wrap(CLC_FracsR)
  save(CLC_FracsR,
       file = file.path(Path_CLC_Summary_RData, "CLC_FracsR.RData"))

  rm(CLC_FracsR)
  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save reference grid - tif
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save reference grid --- tif", Level = 1)
  terra::writeRaster(
    terra::unwrap(Grid_10_Land), overwrite = TRUE,
    filename = file.path(Path_Grid, "Grid_10_Land.tif"))
  terra::writeRaster(
    terra::unwrap(Grid_10_Land_Crop), overwrite = TRUE,
    filename = file.path(Path_Grid, "Grid_10_Land_Crop.tif"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## ||||||||||||||||||||||||||||||||||||||||
  # Add country name to grid cells
  ## ||||||||||||||||||||||||||||||||||||||||

  # Country boundaries
  EU_BoundCNT <- IASDT.R::LoadAs(EU_Bound) %>%
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

  # Get the nearest country names for some grid cells that their centroid do not overlap with country boundaries
  Grid_CNT_Near <- dplyr::filter(Grid_CNT, is.na(Country)) %>%
    dplyr::select(-"Country") %>%
    sf::st_join(y = EU_BoundCNT, join = sf::st_nearest_feature) %>%
    dplyr::rename(Country2 = "Country")

  # Merge data and add a new column representing whether the country information was retrieved by spatial joining of the grid centroid and country boundaries or estimated as the nearest country
  Grid_CNT <- sf::st_join(x = Grid_CNT, y = Grid_CNT_Near) %>%
    sf::st_join(Grid_10_Land_Crop_sf) %>%
    dplyr::mutate(
      Nearest = is.na(Country),
      Country = dplyr::coalesce(Country, Country2)) %>%
    dplyr::select(-"Country2")

  IASDT.R::SaveAs(
    InObj = Grid_CNT, OutObj = "Grid_10_Land_Crop_sf_Country",
    OutPath = file.path(Path_Grid, "Grid_10_Land_Crop_sf_Country.RData"))

  rm(Grid_10_Land_Crop_sf, Grid_10_Land_sf, Exclude_Area,
     TR, Grid_CNT, Grid_CNT_Near)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Crop % coverage results ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk("Crop % coverage results")

  PercCovMaps <- dplyr::mutate(
    PercCovMaps,
    Map_Crop = purrr::map2(
      .x = Name, .y = Map,
      .f = ~{
        Type <- stringr::str_remove(.x, "PercCov_")
        IASDT.R::CatTime(Type, Level = 2)

        Map <- terra::crop(.y, terra::unwrap(Grid_10_Land_Crop)) %>%
          terra::mask(terra::unwrap(Grid_10_Land_Crop)) %>%
          IASDT.R::setRastCRS()

        terra::writeRaster(
          x = Map, overwrite = TRUE,
          filename = file.path(
            Path_CLC_Summary_Tif_Crop, paste0("PercCov_", names(Map), ".tif")))

        OutObjName <- paste0(Type, "_Crop")
        IASDT.R::SaveAs(
          InObj = terra::wrap(Map), OutObj = OutObjName,
          OutPath = file.path(
            Path_CLC_Summary_RData, paste0(OutObjName, ".RData")))

        return(Map)
      }))

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Majority per grid cell ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  IASDT.R::InfoChunk("Identify major CLC class per per grid cell")

  ## ||||||||||||||||||||||||||||||||||||||||
  # Processing using exactextractr::exact_extract
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Processing using exactextractr::exact_extract", Level = 1)

  CLC_Majority <- Grid_sf %>%
    # Ensure that the projection of x and y parameters of exactextractr::exact_extract
    # suppress warning: Polygons transformed to raster CRS (EPSG:3035)
    # https://github.com/isciences/exactextractr/issues/103
    sf::st_transform(sf::st_crs(CLC_Rast)) %>%
    dplyr::mutate(
      Majority = exactextractr::exact_extract(
        x = CLC_Rast, y = ., fun = "majority",
        default_value = 44, progress = FALSE)) %>%
    dplyr::left_join(CLC_CrossWalk, by = dplyr::join_by(Majority == Value)) %>%
    tibble::tibble() %>%
    sf::st_as_sf() %>%
    sf::st_transform(crs = 3035)

  rm(CLC_Rast, Grid_sf)

  ## ||||||||||||||||||||||||||||||||||||||||
  # Save majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save majority results", Level = 1)
  save(CLC_Majority,
       file = file.path(Path_CLC_Summary_RData, "CLC_Majority.RData"))

  invisible(gc())

  ## ||||||||||||||||||||||||||||||||||||||||
  # post-processing majority results
  ## ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Post-processing majority results", Level = 1)

  MajorityMaps <- purrr::map_dfr(
    .x = c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"),
    .f = CLC_ProcessMajority,
    CLC_Majority = CLC_Majority, Path_Tif = Path_CLC_Summary_Tif,
    Path_Tif_Crop = Path_CLC_Summary_Tif_Crop,
    Path_RData = Path_CLC_Summary_RData,
    Grid_10_Land = Grid_10_Land, Grid_10_Land_Crop = Grid_10_Land_Crop)

  rm(CLC_Majority, Grid_10_Land, Grid_10_Land_Crop, MajorityMaps)
  invisible(gc())

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪
  # Plotting ----
  # ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪

  if (PlotCLC) {
    IASDT.R::InfoChunk("Plotting")
    c("PercCov_SynHab", "PercCov_CLC_L1", "PercCov_CLC_L2",
      "PercCov_CLC_L3", "PercCov_EUNIS_2019") %>%
      purrr::walk(
        .f = CLC_Plot,
        CLC_Map = PercCovMaps, EU_Map = EUBound_sf$L_03,
        CrossWalk = CLC_CrossWalk,
        Path_JPEG = Path_CLC_Summary_JPEG,
        Path_JPEG_Free = Path_CLC_Summary_JPEG_Free)
  }

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  return(invisible(NULL))

}

# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####
# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####

## |------------------------------------------------------------------------| #
# CLC_GetPerc ------
## |------------------------------------------------------------------------| #

#' CLC_GetPerc
#'
#' This function Calculate % coverage of different cross-walks per grid cell.
#' The function outputs raster and RData files containing the percentage cover
#' information.
#' @param Type A string indicating the type of habitat to process. This has to
#'   be one of one of SynHab, CLC_L1, CLC_L2, CLC_L3, and EUNIS_2019.
#' @param CLC_CrossWalk A data frame containing the crosswalk data between
#'   values and habitat types.
#' @param CLC_FracsR A list of rasters for percent coverage results.
#' @param Path_Tif A string specifying the path where the output TIFF files will
#'   be saved.
#' @param Path_RData A string specifying the path where the output RData files
#'   will be saved.
#' @noRd
#' @keywords internal
#' @return A tibble containing the processed raster object.
#' @author Ahmed El-Gabbas

CLC_GetPerc <- function(Type, CLC_CrossWalk, CLC_FracsR, Path_Tif, Path_RData) {

  if (is.null(Type) || is.null(CLC_CrossWalk) || is.null(CLC_FracsR) ||
      is.null(Path_Tif) || is.null(Path_RData)) {
    stop("None of the input parameters can be empty", .call = FALSE)
  }

  if (!(
    Type %in% c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"))) {
    stop("Type has to be one of SynHab, CLC_L1, CLC_L2, CLC_L3, and EUNIS_2019", .call = FALSE)
  }

  Fracs <- Class <- HabPerc <- NULL

  IASDT.R::CatTime(Type, Level = 1)
  OutObjName <- paste0("PercCov_", Type)

  Map <- dplyr::select(
    CLC_CrossWalk, "Value", tidyselect::all_of(Type)) %>%
    stats::setNames(c("Fracs", "Class")) %>%
    tidyr::nest(Fracs = -Class) %>%
    dplyr::slice(gtools::mixedorder(Class)) %>%
    dplyr::mutate(
      Fracs = purrr::map(.x = Fracs, .f = ~as.vector(unlist(.x))),
      HabPerc = purrr::map2(
        .x = Fracs, .y = Class,
        .f = ~{
          RName <- paste0(Type, "_", .y) %>%
            stringr::str_replace_all("\\.", "") %>%
            stringr::str_trim()

          CLC_FracsR[[paste0("frac_", .x)]] %>%
            sum() %>%
            stats::setNames(RName)
        })) %>%
    dplyr::pull(HabPerc) %>%
    terra::rast()

  terra::writeRaster(
    Map, overwrite = TRUE,
    filename = file.path(Path_Tif, paste0("PercCov_", names(Map), ".tif")))

  IASDT.R::SaveAs(
    InObj = terra::wrap(Map), OutObj = OutObjName,
    OutPath = file.path(Path_RData, paste0(OutObjName, ".RData")))

  return(tibble::tibble(Name = OutObjName, Map = list(Map)))
}

# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####
# ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ####

## |------------------------------------------------------------------------| #
# CLC_ProcessMajority ------
## |------------------------------------------------------------------------| #

#' CLC_ProcessMajority
#'
#' This function processes CLC majority data for a specified type and generates
#' corresponding raster and vector data. It includes functionality to crop and
#' mask the data to a specific region.
#' @param Type A string indicating the type of data to process. Must be one of
#'   "SynHab", "CLC_L1", "CLC_L2", "CLC_L3", or "EUNIS_2019".
#' @param CLC_Majority A data frame containing the majority class data.
#' @param Path_Tif A string specifying the path where the output TIFF files will
#'   be saved.
#' @param Path_Tif_Crop A string specifying the path where the cropped output
#'   TIFF files will be saved.
#' @param Path_RData A string specifying the path where the output RData files
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

CLC_ProcessMajority <- function(
    Type, CLC_Majority, Path_Tif, Path_Tif_Crop, Path_RData,
    Grid_10_Land, Grid_10_Land_Crop) {

  if (is.null(Type) || is.null(CLC_Majority) || is.null(Path_Tif) ||
      is.null(Path_Tif_Crop) || is.null(Path_RData) || is.null(Grid_10_Land) ||
      is.null(Grid_10_Land_Crop)) {
    stop("None of the input parameters can be empty", .call = FALSE)
  }

  if (!(
    Type %in% c("SynHab", "CLC_L1", "CLC_L2", "CLC_L3", "EUNIS_2019"))) {
    stop("Type has to be one of SynHab, CLC_L1, CLC_L2, CLC_L3, and EUNIS_2019", .call = FALSE)
  }

  Label <- Class <- ID <- NULL

  IASDT.R::CatTime(Type, Level = 2)
  OutObjName <- paste0("Majority_", Type)
  OutObjName_Cr <- paste0(OutObjName, "_Crop")

  Map <- dplyr::filter(CLC_Majority, !is.na(Type)) %>%
    dplyr::select(
      tidyselect::starts_with(Type), -tidyselect::ends_with("_desc")) %>%
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
      tidyselect::all_of(OutObjName),
      tidyselect::everything())

  Map <- terra::classify(Map, MapLevelsM[, -1])
  levels(Map) <- list(MapLevelsNew)

  NAClasses <- c(
    "Marine_Marine habitats", "5_Water bodies",
    "5_2_Marine waters", "5_2_3_Sea and ocean", "A_Marine habitats")

  VV <- stats::setNames(MapLevelsNew, c("ID", "Class")) %>%
    dplyr::filter(Class %in% NAClasses) %>%
    dplyr::pull(ID)

  Map <- IASDT.R::setRastCRS(terra::classify(Map, cbind(NA, VV)))
  levels(Map) <- list(MapLevelsNew)

  terra::writeRaster(
    x = Map, overwrite = TRUE,
    filename = file.path(Path_Tif, paste0(OutObjName, ".tif")))

  terra::levels(Map) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = file.path(Path_Tif, paste0(OutObjName, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)

  IASDT.R::SaveAs(
    InObj = terra::wrap(Map), OutObj = OutObjName,
    OutPath = file.path(Path_RData, paste0(OutObjName, ".RData")))

  # CROPPING
  Map_Cr <- terra::crop(x = Map, y = terra::unwrap(Grid_10_Land_Crop)) %>%
    terra::mask(mask = terra::unwrap(Grid_10_Land_Crop)) %>%
    IASDT.R::setRastCRS()

  terra::writeRaster(
    x = Map_Cr, overwrite = TRUE,
    filename = file.path(Path_Tif_Crop, paste0(OutObjName_Cr, ".tif")))

  terra::levels(Map_Cr) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(VALUE = ID) %>%
    foreign::write.dbf(
      file = file.path(Path_Tif_Crop, paste0(OutObjName_Cr, ".tif.vat.dbf")),
      factor2char = TRUE, max_nchar = 254)

  IASDT.R::SaveAs(
    InObj = terra::wrap(Map_Cr), OutObj = OutObjName_Cr,
    OutPath = file.path(Path_RData, paste0(OutObjName_Cr, ".RData")))

  return(tibble::tibble(Type = Type, Map = list(Map), Map_Cr = list(Map_Cr)))
}
