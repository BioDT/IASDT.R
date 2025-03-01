# # |------------------------------------------------------------------------| #
# GBIF_ReadChunk ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 4

GBIF_ReadChunk <- function(
    ChunkFile, EnvFile = ".env", MaxUncert = 10L, StartYear = 1981L,
    SaveRData = TRUE, ReturnData = FALSE, Overwrite = FALSE) {

  # # ..................................................................... ###

  if (isFALSE(SaveRData) && isFALSE(ReturnData)) {
    stop(
      "At least one of SaveRData and ReturnData has to be `TRUE`",
      call. = FALSE)
  }

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
  # official parameters (overriding the ones from GeoTIFF keys)
  # see: https://stackoverflow.com/questions/78007307
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # # ..................................................................... ###

  # Checking arguments ----
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("ChunkFile", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("SaveRData", "ReturnData", "Overwrite"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("MaxUncert", "StartYear"))

  # # ..................................................................... ###

  ChunkOutPath <- stringr::str_replace(ChunkFile, ".txt$", ".RData")

  if (isFALSE(Overwrite) && IASDT.R::CheckData(ChunkOutPath, warning = FALSE)) {
    if (ReturnData) {
      return(IASDT.R::LoadAs(ChunkOutPath))
    } else {
      return(invisible(NULL))
    }
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CLC_Tif <- SynHab_desc <- CLC_CW <- Longitude <- Latitude <- UncertainKm <-
    CountryCodes <- countryName <- hasCoordinate <- hasGeospatialIssues <-
    phylum <- phylumKey <- Path_Grid <- occurrenceStatus <- Value <-
    Path_GBIF <- SelectedCols <- Int_cols <- lgl_cols <- Dbl_cols <-
    Int64_cols <- coordinatePrecision <- NDecLong <- NDecLat <- year <-
    taxonRank <- SortCols <- CellCode <- NULL

  # # ..................................................................... ###

  # Environment variables
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "CLC_Tif", "DP_R_CLC_tif", FALSE, TRUE,
    "CLC_CW", "DP_R_CLC_crosswalk", FALSE, TRUE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "Path_GBIF_Interim", "DP_R_GBIF_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  load(IASDT.R::Path(Path_GBIF, "SelectedCols.RData"))

  # # ..................................................................... ###

  # Grid_10_Land_Crop
  GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop("Reference grid file not found at: ", GridR, call. = FALSE)
  }
  GridR <- terra::unwrap(IASDT.R::LoadAs(GridR))

  # # Grid_10_Land_Crop_sf
  GridSf <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    stop("Reference grid file (sf) not found at: ", GridSf, call. = FALSE)
  }
  GridSf <- IASDT.R::LoadAs(GridSf)

  # CLC - tif
  Corine <- terra::rast(CLC_Tif)
  terra::activeCat(Corine) <- 0

  # CLC cross-walk to match observations
  CLC_Levels <- readr::read_delim(file = CLC_CW, show_col_types = FALSE) %>%
    dplyr::select(-SynHab_desc)

  ChunkData <- readr::read_lines(ChunkFile, progress = FALSE) %>%
    # read the data by lines and convert to tibble
    dplyr::tibble(Data = .) %>%
    # split into columns and assign column names
    tidyr::separate_wider_delim(
      cols = "Data", delim = "\t", names = SelectedCols$Col) %>%
    # Convert classes for some columns
    dplyr::mutate(
      # convert empty strings to NA
      dplyr::across(tidyselect::everything(), ~ dplyr::na_if(., "")),
      # number of decimal places for longitude / latitude
      NDecLong = purrr::map_int(Longitude, IASDT.R::NDecimals),
      NDecLat = purrr::map_int(Latitude, IASDT.R::NDecimals),
      # change column classes
      dplyr::across(tidyselect::all_of(Int_cols), as.integer),
      dplyr::across(tidyselect::all_of(lgl_cols), as.logical),
      dplyr::across(tidyselect::all_of(Dbl_cols), as.double),
      dplyr::across(tidyselect::all_of(Int64_cols), bit64::as.integer64),
      # convert uncertainty to kilometer
      UncertainKm = UncertainKm / 1000) %>%
    # filtering data
    dplyr::filter(
      # exclude occurrences with empty coordinates
      !is.na(Longitude) | !is.na(Latitude),
      # only occurrences with coordinates and without spatial issues
      hasCoordinate, !hasGeospatialIssues,
      # exclude high spatial uncertainty (keep empty uncertainty values)
      UncertainKm <= MaxUncert | is.na(UncertainKm),
      # exclude occurrences with less precision (keep empty precision values)
      coordinatePrecision <= 0.05 | is.na(coordinatePrecision),
      # exclude occurrences if either latitude/longitude has no decimals
      # (integer)
      (NDecLong > 0 & NDecLat > 0),
      # keep only occurrences recorder after specific year
      year >= StartYear,
      # only "PRESENT" data (i.e. exclude ABSENT)
      occurrenceStatus == "PRESENT" | is.na(occurrenceStatus),
      # only vascular plants (not necessary as all requested data are for
      # vascular plants)
      phylum == "Tracheophyta", phylumKey == 7707728,
      # only accepted taxonomy ranks (species level and below)
      taxonRank %in% c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")) %>%
    # re-order columns
    dplyr::select(tidyselect::all_of(SortCols), tidyselect::everything()) %>%
    # convert to sf object (keep original coordinate columns)
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    # project into EPSG:3035
    sf::st_transform(crs = 3035) %>%
    # Extract coordinates in the new projection
    IASDT.R::sf_add_coords("Longitude_3035", "Latitude_3035") %>%
    # add country name (match original data iso name for countries)
    dplyr::left_join(y = CountryCodes, by = "countryCode") %>%
    dplyr::relocate(countryName, .after = "countryCode") %>%
    # add CellCode
    sf::st_join(GridSf) %>%
    # remove points not overlapping with land grid cells
    dplyr::filter(!is.na(CellCode)) %>%
    dplyr::select(
      -hasCoordinate, -hasGeospatialIssues, -phylum, -phylumKey,
      -occurrenceStatus)

  if (nrow(ChunkData) > 0) {
    # Extract CLC data for occurrences (original data at resolution of 100 m)
    # extract coordinates
    ChunkData <- sf::st_coordinates(ChunkData) %>%
      # extract matching CLC class
      terra::extract(x = Corine, y = .) %>%
      dplyr::pull(Value) %>%
      # convert to tibble (integer)
      dplyr::tibble(Value = .) %>%
      dplyr::mutate(
        # replace 999 with NA
        Value = dplyr::if_else(
          as.integer(Value) == 999, NA_integer_, as.integer(Value))) %>%
      # add information on CLC (L1//L2/L3) classes
      dplyr::left_join(CLC_Levels, by = "Value") %>%
      dplyr::select(-Value) %>%
      # merge with cleaned dataset
      dplyr::bind_cols(ChunkData, .)


    if (SaveRData) {
      ChunkOutName <- stringr::str_remove_all(basename(ChunkFile), ".txt$")
      IASDT.R::SaveAs(
        InObj = ChunkData, OutObj = ChunkOutName, OutPath = ChunkOutPath)
    }

    if (ReturnData) {
      return(ChunkData)
    } else {
      return(invisible(NULL))
    }
  } else {
    # If there are no observations after filtering and SaveData = TRUE, return
    # an empty tibble. This is useful to indicate that this chunk was processed
    # to avoid errors from GBIF_Process function

    if (SaveRData) {
      ChunkOutName <- stringr::str_remove_all(basename(ChunkFile), ".txt$")
      IASDT.R::SaveAs(
        InObj = tibble::tibble(),
        OutObj = ChunkOutName, OutPath = ChunkOutPath)
    }

    return(invisible(NULL))
  }
}
