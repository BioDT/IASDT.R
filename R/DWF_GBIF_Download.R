# # |------------------------------------------------------------------------| #
# GBIF_Download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 3

GBIF_Download <- function(
    FromHPC = TRUE, EnvFile = ".env", Renviron = ".Renviron",
    Request = TRUE, Download = TRUE, SplitChunks = TRUE,
    ChunkSize = 50000L, Boundaries = c(-30, 50, 25, 75), StartYear = 1981L) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Renviron", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("FromHPC", "Request", "Download", "SplitChunks"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("ChunkSize", "Boundaries", "StartYear"))

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SortID <- ID <- Col <- Path_GBIF_Interim <- Path_GBIF_Raw <- Path_GBIF <-
    CountryCodes <- TaxaInfo <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime(
    "Ensure that GBIF access information is available",
    Level = 1)
  IASDT.R::GBIF_Check(Renviron = Renviron)

  # # ..................................................................... ###

  IASDT.R::CatTime("Check system commands")
  Commands <- c("unzip", "nl", "head", "cut", "sed", "split")
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::CheckCommands)
  if (!all(CommandsAvail)) {
    Missing <- paste(Commands[!CommandsAvail], collapse = " + ")
    stop("The following command(s) are not available: ", Missing, call. = FALSE)
  }

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_GBIF", "DP_R_GBIF", FALSE, FALSE,
      "Path_GBIF_Raw", "DP_R_GBIF_Raw", FALSE, FALSE,
      "Path_GBIF_Interim", "DP_R_GBIF_Interim", FALSE, FALSE,
      "CountryCodes", "DP_R_CountryCodes", FALSE, TRUE,
      "TaxaInfo", "DP_R_TaxaInfo_RData", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_GBIF", "DP_R_GBIF_Local", FALSE, FALSE,
      "Path_GBIF_Raw", "DP_R_GBIF_Raw_Local", FALSE, FALSE,
      "Path_GBIF_Interim", "DP_R_GBIF_Interim_Local", FALSE, FALSE,
      "CountryCodes", "DP_R_CountryCodes_Local", FALSE, TRUE,
      "TaxaInfo", "DP_R_TaxaInfo_RData_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # Input data ------
  IASDT.R::CatTime("Loading input data")

  ## Create paths -----
  IASDT.R::CatTime("Create paths", Level = 1)
  fs::dir_create(c(Path_GBIF, Path_GBIF_Raw, Path_GBIF_Interim))

  ## Country codes -----
  IASDT.R::CatTime("Country codes", Level = 1)
  CountryCodes <- readr::read_csv(
    file = CountryCodes, col_types = list(readr::col_character()),
    progress = FALSE, col_select = c("countryName", "countryCode"))

  ## Species list -----
  IASDT.R::CatTime("Species list", Level = 1)
  TaxaList <- IASDT.R::LoadAs(TaxaInfo)

  # # ..................................................................... ###

  # Request GBIF data ------

  IASDT.R::CatTime("Request GBIF data")

  if (Request) {
    # This can take 1-3 hours for the data to be ready
    .StartTimeRequest <- lubridate::now(tzone = "CET")

    ## Request GBIF data -----
    IASDT.R::CatTime("Requesting GBIF data", Level = 1)

    # a new DOI will be created; a couple of hours waiting time is expected

    GBIF_Request <- rgbif::occ_download(
      # list of species keys
      rgbif::pred_in("taxonKey", TaxaList$speciesKey),
      # Only with coordinates & no spatial issues
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      # Only after (>=) a certain year
      rgbif::pred_gte("year", StartYear),
      # Only within specific boundaries
      rgbif::pred_within(
        IASDT.R::DownBoundary(
          Left = Boundaries[1], Right = Boundaries[2],
          Bottom = Boundaries[3], Top = Boundaries[4])
      )
    )

    IASDT.R::CatTime("Save data request", Level = 1)
    save(GBIF_Request, file = IASDT.R::Path(Path_GBIF, "GBIF_Request.RData"))

    # Waiting for data to be ready ------
    IASDT.R::CatTime("Waiting for data to be ready", Level = 1)
    StatusDetailed <- rgbif::occ_download_wait(GBIF_Request)

    IASDT.R::CatDiff(
      InitTime = .StartTimeRequest,
      Prefix = "Requesting GBIF data took ", Level = 1)

    IASDT.R::CatTime("Save status details", Level = 1)
    save(
      StatusDetailed,
      file = IASDT.R::Path(Path_GBIF, "StatusDetailed.RData"))

    IASDT.R::CatTime("Data is ready - status summary:", ... = "\n", Level = 1)
    print(rgbif::occ_download_meta(key = StatusDetailed$key))
  } else {
    Path_Request <- IASDT.R::Path(Path_GBIF, "GBIF_Request.RData")

    if (!file.exists(Path_Request)) {
      stop(
        "Path to previously requested data does not exist:", Path_Request,
        call. = FALSE)
    }

    IASDT.R::CatTime("Loading previous GBIF request", Level = 1)
    GBIF_Request <- IASDT.R::LoadAs(Path_Request)

    Path_Status <- IASDT.R::Path(Path_GBIF, "StatusDetailed.RData")

    if (!file.exists(Path_Status)) {
      stop("Path to status info does not exist:", Path_Status, call. = FALSE)
    }

    IASDT.R::CatTime(
      "Loading `status` information from a previous data request",
      Level = 1)
    StatusDetailed <- IASDT.R::LoadAs(Path_Status)
  }

  # # ..................................................................... ###

  # Download the data to disk ----

  IASDT.R::CatTime("Download GBIF data")

  if (Download) {
    IASDT.R::CatTime("Downloading GBIF data", Level = 1)

    .StartTimeDownload <- lubridate::now(tzone = "CET")

    Dwn <- rgbif::occ_download_get(
      GBIF_Request, path = Path_GBIF_Raw, overwrite = FALSE) %>%
      suppressMessages()

    IASDT.R::CatTime(
      "GBIF data was downloaded at the following path:", Level = 3)
    IASDT.R::CatTime(as.character(Dwn), Level = 2)
    IASDT.R::CatDiff(
      InitTime = .StartTimeDownload,
      Prefix = "Downloading GBIF data took ", Level = 2)

    # Extract/save metadata info
    IASDT.R::CatTime("Extract/save metadata info", Level = 1)
    GBIF_Metadata <- list(
      GBIF_Request = GBIF_Request,
      StatusDetailed = StatusDetailed,
      Citation = attr(GBIF_Request, "citation"),
      DwnKey = StatusDetailed$key,
      DOI = StatusDetailed$doi,
      CreatedTime = StatusDetailed$created,
      ModTime = StatusDetailed$modified,
      DownLink = StatusDetailed$downloadLink,
      FileSizeM = StatusDetailed$size / (1024 * 1024),
      NRecordsM = StatusDetailed$totalRecords / 1e6,
      NDatasets = StatusDetailed$numberDatasets,
      Status = StatusDetailed$status,
      DwnPath = as.character(Dwn))

    save(GBIF_Metadata, file = IASDT.R::Path(Path_GBIF, "GBIF_Metadata.RData"))
  } else {
    IASDT.R::CatTime("Data was NOT downloaded", Level = 1)

    GBIF_Metadata <- IASDT.R::Path(Path_GBIF, "GBIF_Metadata.RData")
    if (!file.exists(GBIF_Metadata)) {
      stop(
        "GBIF metadata file does not exist: ", GBIF_Metadata, call. = FALSE)
    }

    GBIF_Metadata <- IASDT.R::LoadAs(GBIF_Metadata)
  }

  # # ..................................................................... ###

  # Selected columns and data types ----

  IASDT.R::CatTime("Prepare selected columns and data types")

  # Class types of some selected columns
  # integer columns
  Int_cols <- c(
    "LineNum", "day", "month", "year", "speciesKey", "acceptedNameUsageID",
    "taxonKey", "acceptedTaxonKey", "phylumKey", "classKey", "orderKey",
    "familyKey", "genusKey")

  # integer64 columns
  Int64_cols <- c("gbifID", "catalogNumber")

  # double columns
  Dbl_cols <- c("UncertainKm", "Longitude", "Latitude", "coordinatePrecision")

  # logical columns
  lgl_cols <- c("hasCoordinate", "hasGeospatialIssues")

  # Names of selected columns
  SelectedCols <- c(
    # data / source ID
    "gbifID", "datasetKey", "catalogNumber", "occurrenceID", "basisOfRecord",
    "publisher", "ownerInstitutionCode", "institutionCode", "datasetName",
    "datasetID",

    # spatial info
    "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters",
    "coordinatePrecision", "georeferenceSources",
    "georeferenceVerificationStatus", "georeferenceRemarks", "hasCoordinate",
    "hasGeospatialIssues", "georeferenceProtocol", "footprintWKT",

    # Location info
    "countryCode", "publishingCountry", "continent", "verbatimLocality",
    "locality", "stateProvince", "county", "municipality",

    # time info
    "day", "month", "year",

    # ecological info
    "issue", "occurrenceStatus", "habitat", "lifeStage", "samplingProtocol",
    "type",

    # taxonomy
    "phylum", "phylumKey", "class", "classKey", "order", "orderKey", "family",
    "familyKey", "genus", "genusKey", "species", "speciesKey", "vernacularName",
    "scientificName", "acceptedScientificName", "taxonomicStatus", "taxonKey",
    "taxonID", "acceptedTaxonKey", "taxonRank", "verbatimScientificName",
    "specificEpithet", "infraspecificEpithet", "originalNameUsage",
    "acceptedNameUsage", "acceptedNameUsageID", "previousIdentifications") %>%
    dplyr::tibble(Col = .) %>%
    # SortID: an ID to sort columns in the final data
    dplyr::mutate(SortID = seq_len(dplyr::n()))


  SelectedCols <-
    # extract column names and their numbers from the zipped file without
    # extraction read first line
    "unzip -p {GBIF_Metadata$DwnPath} occurrence.txt | head -n 1" %>%
    stringr::str_glue() %>%
    IASDT.R::System() %>%
    # Split the first row into column names. Data is tab-separated
    stringr::str_split("\t") %>%
    magrittr::extract2(1) %>%
    dplyr::tibble(Col = .) %>%
    # column number in the original data
    dplyr::mutate(ID = seq_len(dplyr::n())) %>%
    # only keep information for selected columns
    dplyr::right_join(SelectedCols, by = "Col") %>%
    dplyr::mutate(
      # add information on column classes
      Class = dplyr::case_when(
        Col %in% Int_cols ~ "integer",
        Col %in% Int64_cols ~ "integer64",
        Col %in% Dbl_cols ~ "double",
        Col %in% lgl_cols ~ "logical",
        .default = "character"),
      # shorten some columns
      Col = dplyr::case_when(
        Col == "decimalLongitude" ~ "Longitude",
        Col == "decimalLatitude" ~ "Latitude",
        Col == "coordinateUncertaintyInMeters" ~ "UncertainKm",
        .default = Col)
    ) %>%
    # arrange by column number in the original data
    dplyr::arrange(ID) %>%
    # add a new column representing the line number of original data
    dplyr::mutate(ID = as.integer(ID + 1), SortID = as.integer(SortID + 1)) %>%
    dplyr::bind_rows(
      tibble::tibble(
        Col = "LineNum", ID = 1L, SortID = 1L, Class = "integer"),
      .)

  SortCols <- dplyr::pull(dplyr::arrange(SelectedCols, SortID), Col)

  save(
    SelectedCols, Int_cols, Int64_cols, Dbl_cols, lgl_cols,
    SortCols, CountryCodes,
    file = IASDT.R::Path(Path_GBIF, "SelectedCols.RData"))

  # # ..................................................................... ###

  # Split data into chunks using bash ------

  if (SplitChunks) {
    IASDT.R::CatTime("Split data into chunks using bash")

    # This bash command implements the following:
    ## unzip: read the content of occurrences without extracting its content
    ## nl: add line number to each line, followed by a tab
    ## cut: select only specific columns
    ## sed: exclude first row containing header (column names)
    ## split: split into smaller chunks

    # ensure that ChunkSize is not formatted in scientific notation
    ChunkSize <- format(ChunkSize, scientific = FALSE)

    paste0(
      "unzip -p {GBIF_Metadata$DwnPath} occurrence.txt |",
      ' nl -w1 -n "ln" -s "\t" | cut -f{stringr::str_c(SelectedCols$ID, ',
      'collapse = ",")} -d "\t" | sed -n "1!p" | split -l {ChunkSize} ',
      '-a 3 -d - "{Path_GBIF_Interim}/Chunk_" --additional-suffix=.txt') %>%
      stringr::str_glue() %>%
      IASDT.R::System(RObj = FALSE) %>%
      invisible()
  } else {
    IASDT.R::CatTime("`SplitChunks = FALSE`; no data split was made")
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = paste0(
      "Requesting/downloading GBIF data and split data into chunks took "))

  return(invisible(NULL))
}
