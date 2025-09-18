# # |------------------------------------------------------------------------| #
# gbif_download ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name gbif_data
#' @rdname gbif_data
#' @order 2

gbif_download <- function(
    env_file = ".env", r_environ = ".Renviron", request = TRUE, download = TRUE,
    split_chunks = TRUE, chunk_size = 50000L,
    boundaries = c(-30L, 50L, 25L, 75L), start_year = 1981L) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")

  ecokit::check_args("r_environ", args_type = "character")
  ecokit::check_args(
    args_to_check = c("request", "download", "split_chunks"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c("chunk_size", "start_year", "boundaries"),
    args_type = "numeric", arg_length = c(1L, 1L, 4L))

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  sort_id <- id <- column <- path_gbif_interim <- path_gbif_raw <- path_gbif <-
    country_codes <- taxa_info <- NULL

  # # ..................................................................... ###

  ecokit::cat_time(
    "Ensure that GBIF access information is available", level = 1L)
  ecokit::check_gbif(r_environ = r_environ)

  # # ..................................................................... ###

  ecokit::cat_time("Check system commands")
  commands_to_check <- c("unzip", "nl", "head", "cut", "sed", "split")
  commands_avail <- purrr::map_lgl(
    commands_to_check, ecokit::check_system_command)
  if (!all(commands_avail)) {
    missing_commands <- paste(
      commands_to_check[!commands_avail], collapse = " + ")
    ecokit::stop_ctx(
      "Missing commands", missing_commands = missing_commands,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_gbif", "DP_R_gbif_processed", FALSE, FALSE,
    "path_gbif_raw", "DP_R_gbif_raw", FALSE, FALSE,
    "path_gbif_interim", "DP_R_gbif_interim", FALSE, FALSE,
    "country_codes", "DP_R_country_codes", FALSE, TRUE,
    "taxa_info", "DP_R_taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # Input data ------
  ecokit::cat_time("Loading input data")

  ## Create paths -----
  ecokit::cat_time("Create paths", level = 1L)
  fs::dir_create(c(path_gbif, path_gbif_raw, path_gbif_interim))

  ## Country codes -----
  ecokit::cat_time("Country codes", level = 1L)
  country_codes <- readr::read_csv(
    file = country_codes, col_types = list(readr::col_character()),
    progress = FALSE, col_select = c("countryName", "countryCode"))

  ## Species list -----
  ecokit::cat_time("Species list", level = 1L)
  taxa_list <- ecokit::load_as(taxa_info)

  # # ..................................................................... ###

  # request GBIF data ------

  ecokit::cat_time("Request GBIF data")

  if (request) {
    # This can take 1-3 hours for the data to be ready
    .start_time_request <- lubridate::now(tzone = "CET")

    ## Request GBIF data -----
    ecokit::cat_time("Requesting GBIF data", level = 1L)

    # a new DOI will be created; a couple of hours waiting time is expected

    gbif_request <- rgbif::occ_download(
      # list of species keys
      rgbif::pred_in("taxonKey", taxa_list$speciesKey),
      # Only with coordinates & no spatial issues
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      # Only after (>=) a certain year
      rgbif::pred_gte("year", start_year),
      # Only within specific boundaries
      rgbif::pred_within(
        ecokit::boundary_to_wkt(
          left = boundaries[1], right = boundaries[2],
          bottom = boundaries[3], top = boundaries[4])
      )
    )

    ecokit::cat_time("Save data request", level = 1L)
    save(gbif_request, file = fs::path(path_gbif, "gbif_request.RData"))

    # Waiting for data to be ready ------
    ecokit::cat_time("Waiting for data to be ready", level = 1L)
    gbif_status <- rgbif::occ_download_wait(gbif_request)

    ecokit::cat_diff(
      init_time = .start_time_request,
      prefix = "Requesting GBIF data took ", level = 1L)

    ecokit::cat_time("Save status details", level = 1L)
    save(gbif_status, file = fs::path(path_gbif, "gbif_status.RData"))

    ecokit::cat_time("Data is ready - status summary:", ... = "\n", level = 1L)
    print(rgbif::occ_download_meta(key = gbif_status$key))
  } else {
    path_request <- fs::path(path_gbif, "gbif_request.RData")

    if (!file.exists(path_request)) {
      ecokit::stop_ctx(
        "Path to previously requested data does not exist",
        path_request = path_request, include_backtrace = TRUE)
    }

    ecokit::cat_time("Loading previous GBIF request", level = 1L)
    gbif_request <- ecokit::load_as(path_request)

    path_status <- fs::path(path_gbif, "gbif_status.RData")

    if (!file.exists(path_status)) {
      ecokit::stop_ctx(
        "Path to status info does not exist", path_status = path_status,
        include_backtrace = TRUE)
    }

    ecokit::cat_time(
      "Loading `status` information from a previous data request",
      level = 1L)
    gbif_status <- ecokit::load_as(path_status)
  }

  # # ..................................................................... ###

  # download the data to disk ----

  ecokit::cat_time("Download GBIF data")

  if (download) {
    ecokit::cat_time("Downloading GBIF data", level = 1L)

    .start_time_download <- lubridate::now(tzone = "CET")

    dwnload <- rgbif::occ_download_get(
      gbif_request, path = path_gbif_raw, overwrite = FALSE) %>%
      suppressMessages()

    ecokit::cat_time(
      "GBIF data was downloaded at the following path:", level = 2L,
      cat_timestamp = FALSE)
    ecokit::cat_time(as.character(dwnload), level = 2L, cat_timestamp = FALSE)
    ecokit::cat_diff(
      init_time = .start_time_download,
      prefix = "Downloading GBIF data took ", level = 2L)

    # Extract/save metadata info
    ecokit::cat_time("Extract and save metadata info", level = 1L)
    gbif_metadata <- list(
      gbif_request = gbif_request,
      gbif_status = gbif_status,
      citation = attr(gbif_request, "citation"),
      download_key = gbif_status$key,
      DOI = gbif_status$doi,
      created_time = gbif_status$created,
      modified_time = gbif_status$modified,
      download_link = gbif_status$downloadLink,
      file_size_mb = gbif_status$size / (1024 * 1024),
      n_records_millions = gbif_status$totalRecords / 1e6,
      n_datasets = gbif_status$numberDatasets,
      status = gbif_status$status,
      download_path = as.character(dwnload))

    save(gbif_metadata, file = fs::path(path_gbif, "gbif_metadata.RData"))
  } else {
    ecokit::cat_time("Data was NOT downloaded", level = 1L)

    gbif_metadata <- fs::path(path_gbif, "gbif_metadata.RData")
    if (!file.exists(gbif_metadata)) {
      ecokit::stop_ctx(
        "GBIF metadata file does not exist", gbif_metadata = gbif_metadata,
        include_backtrace = TRUE)
    }

    gbif_metadata <- ecokit::load_as(gbif_metadata)
  }

  # # ..................................................................... ###

  # Selected columns and data types ----

  ecokit::cat_time("Prepare selected columns and data types")

  # class types of some selected columns
  # integer columns
  int_cols <- c(
    "line_number", "day", "month", "year", "speciesKey", "acceptedNameUsageID",
    "taxonKey", "acceptedTaxonKey", "phylumKey", "classKey", "orderKey",
    "familyKey", "genusKey")

  # integer64 columns
  int64_cols <- c("gbifID", "catalogNumber")

  # double columns
  dbl_cols <- c("uncertain_km", "Longitude", "Latitude", "coordinatePrecision")

  # logical columns
  lgl_cols <- c("hasCoordinate", "hasGeospatialIssues")

  # Names of selected columns
  selected_columns <- c(
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
    dplyr::tibble(column = .) %>%
    # sort_id: an ID to sort columns in the final data
    dplyr::mutate(sort_id = seq_len(dplyr::n()))


  selected_columns <-
    # extract column names and their numbers from the zipped file without
    # extraction read first line
    "unzip -p {gbif_metadata$download_path} occurrence.txt | head -n 1" %>%
    stringr::str_glue() %>%
    ecokit::system_command() %>%
    # Split the first row into column names. Data is tab-separated
    stringr::str_split("\t") %>%
    magrittr::extract2(1) %>%
    dplyr::tibble(column = .) %>%
    # column number in the original data
    dplyr::mutate(id = seq_len(dplyr::n())) %>%
    # only keep information for selected columns
    dplyr::right_join(selected_columns, by = "column") %>%
    dplyr::mutate(
      # add information on column classes
      class = dplyr::case_when(
        column %in% int_cols ~ "integer",
        column %in% int64_cols ~ "integer64",
        column %in% dbl_cols ~ "double",
        column %in% lgl_cols ~ "logical",
        .default = "character"),
      # shorten some columns
      column = dplyr::case_when(
        column == "decimalLongitude" ~ "Longitude",
        column == "decimalLatitude" ~ "Latitude",
        column == "coordinateUncertaintyInMeters" ~ "uncertain_km",
        .default = column)
    ) %>%
    # arrange by column number in the original data
    dplyr::arrange(id) %>%
    # add a new column representing the line number of original data
    dplyr::mutate(
      id = as.integer(id + 1), sort_id = as.integer(sort_id + 1)) %>%
    dplyr::bind_rows(
      tibble::tibble(
        column = "line_number", id = 1L, sort_id = 1L, class = "integer"),
      .)

  sort_columns <- dplyr::pull(dplyr::arrange(selected_columns, sort_id), column)

  save(
    selected_columns, int_cols, int64_cols, dbl_cols, lgl_cols,
    sort_columns, country_codes,
    file = fs::path(path_gbif, "selected_columns.RData"))

  # # ..................................................................... ###

  # Split data into chunks using bash ------

  if (split_chunks) {
    ecokit::cat_time("Split data into chunks using bash")

    # This bash command implements the following:
    ## unzip: read the content of occurrences without extracting its content
    ## nl: add line number to each line, followed by a tab
    ## cut: select only specific columns
    ## sed: exclude first row containing header (column names)
    ## split: split into smaller chunks

    # ensure that chunk_size is not formatted in scientific notation
    chunk_size <- format(chunk_size, scientific = FALSE)

    paste0(
      "unzip -p {gbif_metadata$download_path} occurrence.txt |",
      ' nl -w1 -n "ln" -s "\t" | cut -f{stringr::str_c(selected_columns$id, ',
      'collapse = ",")} -d "\t" | sed -n "1!p" | split -l {chunk_size} ',
      '-a 3 -d - "{path_gbif_interim}/Chunk_" --additional-suffix=.txt') %>%
      stringr::str_glue() %>%
      ecokit::system_command(r_object = FALSE) %>%
      invisible()
  } else {
    ecokit::cat_time("`split_chunks = FALSE`; no data split was made")
  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = paste0(
      "Requesting or downloading GBIF data and split data into chunks took "))

  return(invisible(NULL))
}
