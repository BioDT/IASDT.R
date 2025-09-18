# # |------------------------------------------------------------------------| #
# gbif_read_chunk ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name gbif_data
#' @rdname gbif_data
#' @order 3

gbif_read_chunk <- function(
    chunk_file, env_file = ".env", max_uncertainty = 10L, start_year = 1981L,
    save_rdata = TRUE, return_data = FALSE, overwrite = FALSE) {

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::check_args(args_to_check = "chunk_file", args_type = "character")
  ecokit::check_args(
    args_to_check = c("save_rdata", "return_data", "overwrite"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c("max_uncertainty", "start_year"), args_type = "numeric")

  if (isFALSE(save_rdata) && isFALSE(return_data)) {
    ecokit::stop_ctx(
      "At least one of `save_rdata` and `return_data` has to be `TRUE`",
      save_rdata = save_rdata, return_data = return_data,
      include_backtrace = TRUE)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  chunk_out_path <- stringr::str_replace(chunk_file, ".txt$", ".RData")

  if (isFALSE(overwrite) &&
      ecokit::check_data(chunk_out_path, warning = FALSE)) {
    if (return_data) {
      return(ecokit::load_as(chunk_out_path))
    } else {
      return(invisible(NULL))
    }
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  clc_tif <- synhab_desc <- clc_crosswalk <- Longitude <- Latitude <-
    uncertain_km <- country_codes <- countryName <- hasCoordinate <-
    hasGeospatialIssues <- phylum <- phylumKey <- path_grid <-
    occurrenceStatus <- value <- path_gbif <- selected_columns <- int_cols <-
    lgl_cols <- dbl_cols <- int64_cols <- coordinatePrecision <- n_dec_long <-
    n_dec_lat <- year <- taxonRank <- sort_columns <- CellCode <- NULL

  # # ..................................................................... ###

  # Environment variables

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "clc_tif", "DP_R_clc_tif", FALSE, TRUE,
    "clc_crosswalk", "DP_R_clc_crosswalk", FALSE, TRUE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_gbif", "DP_R_gbif_processed", FALSE, FALSE,
    "path_gbif_Interim", "DP_R_gbif_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  load(fs::path(path_gbif, "selected_columns.RData"))

  # # ..................................................................... ###

  # grid_10_land_crop
  grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_r)) {
    ecokit::stop_ctx(
      "Reference grid file not found", grid_r = grid_r,
      include_backtrace = TRUE)
  }
  grid_r <- ecokit::load_as(grid_r, unwrap_r = TRUE)

  # # grid_10_land_crop_sf
  grid_sf <- fs::path(path_grid, "grid_10_land_crop_sf.RData")
  if (!file.exists(grid_sf)) {
    ecokit::stop_ctx(
      "Reference grid file (sf) not found", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }
  grid_sf <- ecokit::load_as(grid_sf)

  # CLC - tif
  corine_r <- terra::rast(clc_tif)
  terra::activeCat(corine_r) <- 0

  # CLC cross-walk to match observations
  clc_levels <- readr::read_delim(
    file = clc_crosswalk, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(-synhab_desc)

  # publishers for citizen science data
  cz_publishers <- c(
    "Pl@ntNet", "iNaturalist.org", "Observation.org", "naturgucker.de",
    "iSpot", "BioDiversity4All", "Questagame")

  chunk_data <- readr::read_lines(chunk_file, progress = FALSE) %>%
    # read the data by lines and convert to tibble
    dplyr::tibble(Data = .) %>%
    # split into columns and assign column names
    tidyr::separate_wider_delim(
      cols = "Data", delim = "\t", names = selected_columns$Col) %>%
    # Convert classes for some columns
    dplyr::mutate(
      # convert empty strings to NA
      dplyr::across(tidyselect::everything(), ~ dplyr::na_if(., "")),
      # number of decimal places for longitude / latitude
      n_dec_long = purrr::map_int(Longitude, ecokit::n_decimals),
      n_dec_lat = purrr::map_int(Latitude, ecokit::n_decimals),
      # change column classes
      dplyr::across(tidyselect::all_of(int_cols), as.integer),
      dplyr::across(tidyselect::all_of(lgl_cols), as.logical),
      dplyr::across(tidyselect::all_of(dbl_cols), as.double),
      dplyr::across(tidyselect::all_of(int64_cols), bit64::as.integer64),
      # convert uncertainty to kilometre
      uncertain_km = uncertain_km / 1000,
      # publisher type: citizen science vs others
      publisher_type = dplyr::case_when(
        publisher %in% cz_publishers ~ "citizen_science",
        .default = "others")) %>%
    # filtering data
    dplyr::filter(
      # exclude occurrences with empty coordinates
      !is.na(Longitude) | !is.na(Latitude),
      # only occurrences with coordinates and without spatial issues
      hasCoordinate, !hasGeospatialIssues,
      # exclude high spatial uncertainty (keep empty uncertainty values)
      uncertain_km <= max_uncertainty | is.na(uncertain_km),
      # exclude occurrences with less precision (keep empty precision values)
      coordinatePrecision <= 0.05 | is.na(coordinatePrecision),
      # exclude occurrences if either latitude/longitude has no decimals
      # (integer)
      (n_dec_long > 0 & n_dec_lat > 0),
      # keep only occurrences recorder after specific year
      year >= start_year,
      # only "PRESENT" data (i.e. exclude ABSENT)
      occurrenceStatus == "PRESENT" | is.na(occurrenceStatus),
      # only vascular plants (not necessary as all requested data are for
      # vascular plants)
      phylum == "Tracheophyta", phylumKey == 7707728,
      # only accepted taxonomy ranks (species level and below)
      taxonRank %in% c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")) %>%
    # re-order columns
    dplyr::select(
      tidyselect::all_of(sort_columns), tidyselect::everything()) %>%
    # convert to sf object (keep original coordinate columns)
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    # project into EPSG:3035
    sf::st_transform(crs = 3035) %>%
    # Extract coordinates in the new projection
    ecokit::sf_add_coords("Longitude_3035", "Latitude_3035") %>%
    # add country name (match original data iso name for countries)
    dplyr::left_join(y = country_codes, by = "countryCode") %>%
    dplyr::relocate(countryName, .after = "countryCode") %>%
    # add CellCode
    sf::st_join(grid_sf) %>%
    # remove points not overlapping with land grid cells
    dplyr::filter(!is.na(CellCode)) %>%
    dplyr::select(
      -hasCoordinate, -hasGeospatialIssues, -phylum, -phylumKey,
      -occurrenceStatus)

  rm(cz_publishers, envir = environment())

  if (nrow(chunk_data) > 0) {
    # Extract CLC data for occurrences (original data at resolution of 100 m)
    # extract coordinates
    chunk_data <- sf::st_coordinates(chunk_data) %>%
      # extract matching CLC class
      terra::extract(x = corine_r, y = .) %>%
      dplyr::pull(value) %>%
      # convert to tibble (integer)
      dplyr::tibble(value = .) %>%
      dplyr::mutate(
        # replace 999 with NA
        value = dplyr::if_else(
          as.integer(value) == 999, NA_integer_, as.integer(value))) %>%
      # add information on CLC (L1//L2/L3) classes
      dplyr::left_join(clc_levels, by = "value") %>%
      dplyr::select(-value) %>%
      # merge with cleaned dataset
      dplyr::bind_cols(chunk_data, .)

    if (save_rdata) {
      chunk_out_name <- stringr::str_remove_all(basename(chunk_file), ".txt$")
      ecokit::save_as(
        object = chunk_data, object_name = chunk_out_name,
        out_path = chunk_out_path)
    }

    if (return_data) {
      return(chunk_data)
    } else {
      return(invisible(NULL))
    }
  } else {
    # If there are no observations after filtering and SaveData = TRUE, return
    # an empty tibble. This is useful to indicate that this chunk was processed
    # to avoid errors from gbif_process function

    if (save_rdata) {
      chunk_out_name <- stringr::str_remove_all(basename(chunk_file), ".txt$")
      ecokit::save_as(
        object = tibble::tibble(),
        object_name = chunk_out_name, out_path = chunk_out_path)
    }

    return(invisible(NULL))
  }
}
