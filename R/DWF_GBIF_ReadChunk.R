# # |------------------------------------------------------------------------| #
# GBIF_ReadChunk ----
## |------------------------------------------------------------------------| #

#' Read GBIF chunk data
#'
#' This function reads and processes a chunk of GBIF (Global Biodiversity
#' Information Facility) data. It filters the data based on several criteria
#' including spatial uncertainty, collection year, coordinate precision, and
#' taxonomic rank and returns only selected columns.
#'
#' @param ChunkFile A string specifying the path of the chunk file to be
#'   processed.
#' @param FromHPC A logical value indicating whether the environment variables
#'   should be read from an HPC (High-Performance Computing) environment.
#'   Default is TRUE.
#' @param EnvFile A string specifying the path of the environment file. Default
#'   is ".env".
#' @param MaxUncert A numeric value specifying the maximum accepted spatial
#'   uncertainty in kilometers. Default is 10 km.
#' @param MinYear An integer specifying the earliest collection year to be
#'   included. Default is 1981.
#' @param SaveRData Logical; indicating whether to save the cleaned data for the
#'   current chunk as `*.RData` file.
#' @param ReturnData Logical; indicating whether to return the cleaned data for
#'   the current chunk. Defaults to `FALSE`.
#' @param Overwrite Logical; indicating whether to process the current chunk
#'   file if it has already processed and saved as `*.RData` file. This helps to
#'   continue working on previously processed chunks if the previous try failed,
#'   e.g. due to memory issue.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [GBIF_Process] function.
#' @return The output of the function is a tibble (sf) object containing the
#'   processed and filtered GBIF data for the input chunk file. Whether the
#'   tibble is exported or saved as `*.RData` file depends on the values of the
#'   `SaveRData`, `ReturnData`, and `Overwrite` parameters.
#'
#'   If `SaveRData = TRUE` (default), the processed data will be saved as
#'   `RData` file with the same base name of the chunk file and at the same
#'   directory. By default, the function does not return any value, unless
#'   `ReturnData` is set to `TRUE`. `SaveRData` and `ReturnData` can not both
#'   set to `FALSE`. The function can optionally skip processing the current
#'   chunk file if the `*.RData` file for this chunk already exist and
#'   `Overwrite` is set as `FALSE` (default). If `Overwrite` is set to `TRUE`,
#'   the data for the current chunk will be re-processed irrespective of the
#'   existence of the RData file. This helps to continue working on previously
#'   processed chunks if a previous run of [GBIF_Process] failed, e.g. due to
#'   memory issue.
#' @name GBIF_ReadChunk
#' @author Ahmed El-Gabbas
#' @export

GBIF_ReadChunk <- function(
    ChunkFile, EnvFile = ".env", FromHPC = TRUE, MaxUncert = 10,
    MinYear = 1981, SaveRData = TRUE, ReturnData = FALSE, Overwrite = FALSE) {

  if (!SaveRData && !ReturnData) {
    stop(
      "At least one of SaveRData and ReturnData has to be `TRUE`",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("ChunkFile", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("FromHPC", "SaveRData", "ReturnData", "Overwrite"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("MaxUncert", "MinYear"))

  # # ..................................................................... ###

  ChunkOutPath <- stringr::str_replace(ChunkFile, ".txt$", ".RData")

  if (!Overwrite && file.exists(ChunkOutPath)) {
    if (IASDT.R::CheckRData(ChunkOutPath)) {
      if (ReturnData) {
        return(IASDT.R::LoadAs(ChunkOutPath))
      } else {
        return(invisible(NULL))
      }
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
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "CLC_Tif", "DP_R_CLC_tif", FALSE, TRUE,
      "CLC_CW", "DP_R_CLC_CW", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_GBIF", "DP_R_GBIF", FALSE, FALSE,
      "Path_GBIF_Interim", "DP_R_GBIF_Interim", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "CLC_Tif", "DP_R_CLC_tif_Local", FALSE, TRUE,
      "CLC_CW", "DP_R_CLC_CW_Local", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_GBIF", "DP_R_GBIF_Local", FALSE, FALSE,
      "Path_GBIF_Interim", "DP_R_GBIF_Interim_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  load(file.path(Path_GBIF, "SelectedCols.RData"))

  # # ..................................................................... ###

  # Grid_10_Land_Crop
  GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(GridR)) {
    stop(paste0("Reference grid file not found at: ", GridR), call. = FALSE)
  }
  GridR <- terra::unwrap(IASDT.R::LoadAs(GridR))

  # # Grid_10_Land_Crop_sf
  GridSf <- file.path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(GridSf)) {
    stop(
      paste0("Reference grid file (sf) not found at: ", GridSf), call. = FALSE)
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
      dplyr::across(tidyselect::everything(), ~dplyr::na_if(., "")),
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
      # exclude occurrences if either latitude/longitude has no decimals (integer)
      (NDecLong > 0 & NDecLat > 0),
      # keep only occurrences recorder after specific year
      year >= MinYear,
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
        InObj = tibble::tibble(), OutObj = ChunkOutName, OutPath = ChunkOutPath)
    }

    return(invisible(NULL))
  }
}
