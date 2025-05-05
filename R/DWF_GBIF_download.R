# # |------------------------------------------------------------------------| #
# GBIF_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 3

GBIF_download <- function(
    env_file = ".env", r_environ = ".Renviron", request = TRUE, download = TRUE,
    split_chunks = TRUE, chunk_size = 50000L, boundaries = c(-30, 50, 25, 75),
    start_year = 1981L) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::cat_time("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("r_environ", "env_file"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("request", "download", "split_chunks"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("chunk_size", "boundaries", "start_year"))

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SortID <- ID <- Col <- Path_GBIF_Interim <- Path_GBIF_Raw <- Path_GBIF <-
    CountryCodes <- TaxaInfo <- NULL

  # # ..................................................................... ###

  IASDT.R::cat_time(
    "Ensure that GBIF access information is available", level = 1L)
  IASDT.R::GBIF_check(r_environ = r_environ)

  # # ..................................................................... ###

  IASDT.R::cat_time("Check system commands")
  Commands <- c("unzip", "nl", "head", "cut", "sed", "split")
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::check_system_command)
  if (!all(CommandsAvail)) {
    Missing <- paste(Commands[!CommandsAvail], collapse = " + ")
    IASDT.R::stop_ctx("Missing commands", missing_commands = Missing)
  }

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_GBIF", "DP_R_GBIF_processed", FALSE, FALSE,
    "Path_GBIF_Raw", "DP_R_GBIF_raw", FALSE, FALSE,
    "Path_GBIF_Interim", "DP_R_GBIF_interim", FALSE, FALSE,
    "CountryCodes", "DP_R_Countrycodes", FALSE, TRUE,
    "TaxaInfo", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Input data ------
  IASDT.R::cat_time("Loading input data")

  ## Create paths -----
  IASDT.R::cat_time("Create paths", level = 1L)
  fs::dir_create(c(Path_GBIF, Path_GBIF_Raw, Path_GBIF_Interim))

  ## Country codes -----
  IASDT.R::cat_time("Country codes", level = 1L)
  CountryCodes <- readr::read_csv(
    file = CountryCodes, col_types = list(readr::col_character()),
    progress = FALSE, col_select = c("countryName", "countryCode"))

  ## Species list -----
  IASDT.R::cat_time("Species list", level = 1L)
  TaxaList <- IASDT.R::load_as(TaxaInfo)

  # # ..................................................................... ###

  # request GBIF data ------

  IASDT.R::cat_time("Request GBIF data")

  if (request) {
    # This can take 1-3 hours for the data to be ready
    .StartTimeRequest <- lubridate::now(tzone = "CET")

    ## Request GBIF data -----
    IASDT.R::cat_time("Requesting GBIF data", level = 1L)

    # a new DOI will be created; a couple of hours waiting time is expected

    GBIF_Request <- rgbif::occ_download(
      # list of species keys
      rgbif::pred_in("taxonKey", TaxaList$speciesKey),
      # Only with coordinates & no spatial issues
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      # Only after (>=) a certain year
      rgbif::pred_gte("year", start_year),
      # Only within specific boundaries
      rgbif::pred_within(
        IASDT.R::boundary_to_WKT(
          left = boundaries[1], right = boundaries[2],
          bottom = boundaries[3], top = boundaries[4])
      )
    )

    IASDT.R::cat_time("Save data request", level = 1L)
    save(GBIF_Request, file = fs::path(Path_GBIF, "GBIF_Request.RData"))

    # Waiting for data to be ready ------
    IASDT.R::cat_time("Waiting for data to be ready", level = 1L)
    StatusDetailed <- rgbif::occ_download_wait(GBIF_Request)

    IASDT.R::cat_diff(
      init_time = .StartTimeRequest,
      prefix = "Requesting GBIF data took ", level = 1L)

    IASDT.R::cat_time("Save status details", level = 1L)
    save(StatusDetailed, file = fs::path(Path_GBIF, "StatusDetailed.RData"))

    IASDT.R::cat_time("Data is ready - status summary:", ... = "\n", level = 1L)
    print(rgbif::occ_download_meta(key = StatusDetailed$key))
  } else {
    Path_Request <- fs::path(Path_GBIF, "GBIF_Request.RData")

    if (!file.exists(Path_Request)) {
      IASDT.R::stop_ctx(
        "Path to previously requested data does not exist",
        Path_Request = Path_Request)
    }

    IASDT.R::cat_time("Loading previous GBIF request", level = 1L)
    GBIF_Request <- IASDT.R::load_as(Path_Request)

    Path_Status <- fs::path(Path_GBIF, "StatusDetailed.RData")

    if (!file.exists(Path_Status)) {
      IASDT.R::stop_ctx(
        "Path to status info does not exist", Path_Status = Path_Status)
    }

    IASDT.R::cat_time(
      "Loading `status` information from a previous data request",
      level = 1L)
    StatusDetailed <- IASDT.R::load_as(Path_Status)
  }

  # # ..................................................................... ###

  # download the data to disk ----

  IASDT.R::cat_time("Download GBIF data")

  if (download) {
    IASDT.R::cat_time("Downloading GBIF data", level = 1L)

    .StartTimeDownload <- lubridate::now(tzone = "CET")

    Dwn <- rgbif::occ_download_get(
      GBIF_Request, path = Path_GBIF_Raw, overwrite = FALSE) %>%
      suppressMessages()

    IASDT.R::cat_time(
      "GBIF data was downloaded at the following path:", level = 3L)
    IASDT.R::cat_time(as.character(Dwn), level = 2L)
    IASDT.R::cat_diff(
      init_time = .StartTimeDownload,
      prefix = "Downloading GBIF data took ", level = 2L)

    # Extract/save metadata info
    IASDT.R::cat_time("Extract and save metadata info", level = 1L)
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

    save(GBIF_Metadata, file = fs::path(Path_GBIF, "GBIF_Metadata.RData"))
  } else {
    IASDT.R::cat_time("Data was NOT downloaded", level = 1L)

    GBIF_Metadata <- fs::path(Path_GBIF, "GBIF_Metadata.RData")
    if (!file.exists(GBIF_Metadata)) {
      IASDT.R::stop_ctx(
        "GBIF metadata file does not exist", GBIF_Metadata = GBIF_Metadata)
    }

    GBIF_Metadata <- IASDT.R::load_as(GBIF_Metadata)
  }

  # # ..................................................................... ###

  # Selected columns and data types ----

  IASDT.R::cat_time("Prepare selected columns and data types")

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
    IASDT.R::system_command() %>%
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
    file = fs::path(Path_GBIF, "SelectedCols.RData"))

  # # ..................................................................... ###

  # Split data into chunks using bash ------

  if (split_chunks) {
    IASDT.R::cat_time("Split data into chunks using bash")

    # This bash command implements the following:
    ## unzip: read the content of occurrences without extracting its content
    ## nl: add line number to each line, followed by a tab
    ## cut: select only specific columns
    ## sed: exclude first row containing header (column names)
    ## split: split into smaller chunks

    # ensure that chunk_size is not formatted in scientific notation
    chunk_size <- format(chunk_size, scientific = FALSE)

    paste0(
      "unzip -p {GBIF_Metadata$DwnPath} occurrence.txt |",
      ' nl -w1 -n "ln" -s "\t" | cut -f{stringr::str_c(SelectedCols$ID, ',
      'collapse = ",")} -d "\t" | sed -n "1!p" | split -l {chunk_size} ',
      '-a 3 -d - "{Path_GBIF_Interim}/Chunk_" --additional-suffix=.txt') %>%
      stringr::str_glue() %>%
      IASDT.R::system_command(r_object = FALSE) %>%
      invisible()
  } else {
    IASDT.R::cat_time("`split_chunks = FALSE`; no data split was made")
  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .start_time,
    prefix = paste0(
      "Requesting or downloading GBIF data and split data into chunks took "))

  return(invisible(NULL))
}
