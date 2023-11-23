# |---------------------------------------------------| #
# Match_to_GBIF ----
# |---------------------------------------------------| #

#' Match taxonomy with GBIF; may return >1 match
#'
#' Match taxonomy with GBIF; may return >1 match
#' @param taxon_name taxonomy name
#' @param taxon_id taxonomy ID
#' @param include_genus include matches at genus level; default: `FALSE`
#' @param Parallel logical; whether to implement standardization on parallel; default: `FALSE`
#' @param Progress logical; whether to print progress bar; default: `FALSE`
#' @name Match_to_GBIF
#' @importFrom rlang .data
#' @author Marina Golivets
#' @return a tibble for the standardization results
#' @export
#' @details
#' as input, provide a vector of verbatim taxon names (preferably with authorship) and a vector of existing local identifiers for those names

Match_to_GBIF <- function(
    taxon_name, taxon_id = NULL, include_genus = FALSE,
    Parallel = FALSE, Progress = FALSE) {

  if (is.null(taxon_id)) taxon_id <- seq_along(taxon_name)

  if (Parallel) {
    all_matches <- furrr::future_map(
      taxon_name, rgbif::name_backbone_verbose,
      kingdom = "plants", strict = TRUE,
      .progress = Progress, .options = furrr::furrr_options(seed = TRUE))
  } else {
    if (Progress) {
      ProgrOptns <- list(
        type = "iterator", clear = TRUE,
        format = "{cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]")
    } else {
      ProgrOptns <- FALSE
    }

    all_matches <- purrr::map(
      taxon_name, rgbif::name_backbone_verbose,
      kingdom = "plants", strict = TRUE, .progress = ProgrOptns)
  }

  # retrieve alternative matches
  alternative_matches <- lapply(
    all_matches,
    function(x) {
      y <- x$alternatives
      if (nrow(y) == 0) {
        y[1, 1] <- NA
        colnames(y) <- "usageKey"
      } else {
        y <- y
      }
      return(y)
    }
  ) %>%
    mapply(
      cbind, .,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::filter(!is.na(.data$usageKey)) %>%
    dplyr::distinct() %>%
    # filter only if phylum column exists
    dplyr::filter({if ("phylum" %in% names(.)) phylum else NULL} == "Tracheophyta")

  # retrieve best matches
  best_matches <- lapply(all_matches, function(x) x$data) %>%
    mapply(
      cbind, .,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::distinct() %>%
    dplyr::filter({if ("phylum" %in% names(.)) phylum else NULL} == "Tracheophyta")

  matched <- best_matches %>%
    dplyr::filter(!(.data$matchType %in% c("NONE", "HIGHERRANK")))

  matched_alternative <- try(
    alternative_matches %>%
      # use only vascular plants
      # filter only if phylum column exists
      dplyr::filter({if ("phylum" %in% names(.)) phylum else NULL} == "Tracheophyta") %>%
      dplyr::filter(.data$confidence >= 0) %>%
      dplyr::filter(!taxon_id %in% matched$taxon_id)
  )
  if (class(matched_alternative)[1] == "try-error") {
    taxon_list <- matched
  } else {
    taxon_list <- dplyr::bind_rows(matched, matched_alternative)
  }

  if (include_genus == FALSE) {
    taxon_list <- taxon_list %>%
      dplyr::filter(rank != "GENUS")
  }


  # get names that were matched as accepted
  accepted <- taxon_list %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$status == "ACCEPTED")
  if (nrow(accepted) > 0) {
    accepted <- accepted %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    accepted <- dplyr::ungroup(accepted)
  }

  # get names that were matched as synonyms only
  synonyms <- taxon_list %>%
    dplyr::group_by(.data$taxon_id) %>%
    dplyr::summarise(has_accepted = dplyr::n_distinct(.data$status == "ACCEPTED") > 1) %>%
    dplyr::full_join(taxon_list, by = dplyr::join_by(taxon_id)) %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::filter(.data$status == "SYNONYM")
  if (nrow(synonyms) > 0) {
    synonyms <- synonyms %>%
      dplyr::group_by(taxon_id) %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    synonyms <- dplyr::ungroup(synonyms)
  }

  # get names that were matched as doubtful only
  doubtful <- taxon_list %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(has_accepted = dplyr::n_distinct(.data$status == "ACCEPTED") > 1) %>%
    dplyr::full_join(taxon_list, by = dplyr::join_by(taxon_id)) %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$status == "DOUBTFUL")
  if (nrow(doubtful) > 0) {
    doubtful <- doubtful %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    doubtful <- dplyr::ungroup(doubtful)
  }

  # combine all names
  taxon_list_final <- dplyr::bind_rows(accepted, synonyms, doubtful) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$confidence == max(.data$confidence)) %>%
    dplyr::filter(.data$status != "NONE") %>% # exclude non-matched names
    dplyr::select(-"has_accepted") %>%
    dplyr::ungroup() %>%
    dplyr::relocate(taxon_name, taxon_id)

  return(taxon_list_final)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Get_EASIN_Data ----
# |---------------------------------------------------| #

#' Extract EASIN data using EASIN ID
#'
#' Extract EASIN data using EASIN ID
#' @param SpKey EASIN ID
#' @param NSearch number of grid cells to extract each time
#' @name Get_EASIN_Data
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @details
#' A function to extract EASIN data for a given EASIN_ID
#' @examples
#' \dontrun{
#' c("R00006", "R18689", "R00008") %>%
#'   purrr::map(Get_EASIN_Data)
#'}

Get_EASIN_Data <- function(SpKey, NSearch = 500) {

  withr::local_options(list(scipen = 999))

  Skip <- 0         # Skip = 0; start at the first presence grid
  ID <- 0           # iteration ID
  DT <- list()      # List object to save the data

  # Looping over data chunks
  while (TRUE) {
    ID <- ID + 1
    URL <- stringr::str_glue("https://alien.jrc.ec.europa.eu/apixg/geoxg/speciesid/{SpKey}/layertype/grid/skip/{Skip}/take/{NSearch}")

    # Extract data from JSON
    Data <- jsonlite::fromJSON(URL)

    if (all(names(Data) == "Empty")) {
      # If there is no data for this EASIN ID, quit the loop and return a value of NULL
      DT[[ID]] <- NULL
      break()
    } else {
      # Add data on the current chunk to the list
      DT[[ID]] <- tibble::tibble(Data)
    }

    # if the number of grids in the current chunk < number of grids per chunks break the loop
    if (nrow(Data) < NSearch) {
      break()
    } else {
      Skip <- Skip + NSearch
    }
  }

  # merge the list items together
  if (length(DT) > 0) {
    return(dplyr::bind_rows(DT))
  } else {
    return(NULL)
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# ChangeClass ----
# |---------------------------------------------------| #

#' Change the class of some of GBIF columns
#'
#' Change the class of some of GBIF columns
#' @param DF Data frame
#' @name ChangeClass
#' @author Ahmed El-Gabbas
#' @details
#' A function to change the class of some of GBIF columns
#' @return NULL
#' @keywords internal
#' @noRd

ChangeClass <- function(DF) {
  VarsInt <- c(
    "year", "month", "day", "coordinateUncertaintyInMeters", "acceptedNameUsageID",
    "taxonKey", "acceptedTaxonKey", "kingdomKey", "phylumKey", "classKey",
    "orderKey", "familyKey", "genusKey", "speciesKey")
  VarsInt64 <- c("gbifID", "catalogNumber", "recordNumber", "taxonID", "identificationID")
  VarsDbl <- c("decimalLatitude", "decimalLongitude")
  VarsLgl <- c("hasCoordinate", "hasGeospatialIssues", "repatriated")
  VarsDte <- c("modified", "eventDate", "dateIdentified", "lastInterpreted", "lastParsed", "lastCrawled")

  DF %>%
    dplyr::mutate_at(VarsInt, as.integer) %>%
    dplyr::mutate_at(VarsInt64, bit64::as.integer64) %>%
    dplyr::mutate_at(VarsDbl, as.double)  %>%
    dplyr::mutate_at(VarsLgl, as.logical)  %>%
    dplyr::mutate_at(VarsDte, lubridate::as_date) %>%
    return()
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# GetAcceptedName ----
# |---------------------------------------------------| #

#' Get accepted name of a taxa
#'
#' Get accepted name of a taxa
#' @name GetAcceptedName
#' @param ID GBIF id
#' @author Ahmed El-Gabbas
#' @return Scientific name of the accepted taxa
#' @export
#' @details
#' Get accepted name of a taxa
#' @examples
#' # https://www.gbif.org/species/5372559
#' GetAcceptedName(5372559)

GetAcceptedName <- function(ID) {
  if (!is.na(ID)) {
    rgbif::name_usage(ID)$data$scientificName
  } else {
    NA_character_
  }
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Extract_BB ----
# |---------------------------------------------------| #

#' Extract a specific column from the output of GBIF standardization
#'
#' Extract a specific column from the output of GBIF standardization
#' @name Extract_BB
#' @param x GBIF results of GBIF standardization `rgbif::name_backbone`
#' @param var column name to extract
#' @author Ahmed El-Gabbas
#' @return Scientific name of the accepted taxa
#' @export
#' @details
#' Extract a specific column from the output of GBIF
#' @examples
#' rgbif::name_backbone(name = "Helianthus annuus", kingdom = "plants") %>%
#'    Extract_BB(status)
#'
#' # ------------------------
#'
#' c("Helianthus annuus", "Tagetes patula L.") %>%
#'    tibble::tibble(Taxa = .) %>%
#'    dplyr::mutate(
#'       BB = purrr::map(Taxa, rgbif::name_backbone),
#'       status = purrr::map_chr(BB, Extract_BB, status),
#'       SpKey = purrr::map_int(BB, Extract_BB, speciesKey))

Extract_BB <- function(x = ., var) {
  var <- rlang::ensyms(var) %>%
    as.character()
  if (var %in% names(x)) {
    x[[var]]
  } else {
    NA
  }
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Extract_Matching ----
# |---------------------------------------------------| #

#' Extract time / climate models & scenarios from Chelsa URL
#'
#' A function to extract which climate model and scenario each file represent
#' @name Chelsa_Extract_Matching
#' @param String URL string
#' @param Time Time period to match
#' @param Matches further string to match
#' @author Ahmed El-Gabbas
#' @export
#' @details
#' A function to extract which climate model and scenario each file represent

Chelsa_Extract_Matching <- function(String, Time, Matches) {
  # assign "Current" for files represent current climates
  if (Time == "1981-2010") return("Current")

  # Index of matched text
  Which <- sapply(
    Matches,
    function(x) {
      grepl(x, String, ignore.case = "True")
    }) %>%
    which()

  if (length(Which) > 0) {
    return(Matches[Which])
  } else {
    return(NA_character_)
  }
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Prepare_List ----
# |---------------------------------------------------| #

#' Prepare list of potential variables
#'
#' Prepare list of potential variables
#' @name Chelsa_Prepare_List
#' @param Down Should the Chelsa files be downloaed
#' @param DownParallel if `Down` was set as `TRUE`, should the download be on parallel
#' @param DwnPath Path for download
#' @author Ahmed El-Gabbas
#' @export
#' @details
#' list of variables exist under current and future climates
#' 46 variables available at 46 options (current and 45 future scenarios)

Chelsa_Prepare_List <- function(Down = FALSE, DownParallel = TRUE, DwnPath = NULL) {

  IASDT.R::CatTime("Preparing Chelsa climatology data")

  BaseURL <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/"
  ClimateModels <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  ClimateScenarios <- c("ssp126", "ssp370", "ssp585")

  # Variables to exclude
  Exclude <- c(
    "ai", "hurs", "clt", "sfcWind", "vpd", "rsds", "pet", "cmi", "swb",
    "pr_", "tasmax_", "tasmin_", "tas_") %>%
    stringr::str_c(collapse = "|")

  # List of Chelsa Vars we are interested in (available in current/future)
  Chelsa_Vars <- file.path(Path_Chelsa, "Chelsa_Vars.txt") %>%
    read_delim(delim = "\t", show_col_types = FALSE) %>%
    dplyr::select(-scale, -offset)

  ChelsaClimData <- file.path(Path_Chelsa, "DwnLinks") %>%
    # files containing download links for climatology data
    list.files(
      pattern = "DwnLinks_Climatologies_.+txt$", recursive = TRUE,
      full.names = TRUE) %>%
    dplyr::tibble(URL_File = .) %>%
    dplyr::group_by(URL_File) %>%

    # Add download links
    dplyr::mutate(
      URL = purrr::map(
        .x = URL_File,
        .f = ~{
          .x %>%
            readr::read_lines() %>%
            trimws() # remove trailing spaces
        }),
      URL_File = basename(URL_File)) %>%
    tidyr::unnest(cols = c(URL)) %>%
    dplyr::ungroup() %>%

    dplyr::mutate(
      # The name of the downloaded file and folder
      Folder = stringr::str_remove_all(string = URL, pattern = BaseURL),
      File = basename(Folder),
      Folder = dirname(Folder),

      # Extract time period
      TimePeriod = stringr::str_remove_all(
        string = URL_File, pattern = "DwnLinks_Climatologies_|.txt"),

      # File extension
      Ext = tools::file_ext(URL),

      # which climate model
      ClimModel = purrr::map2_chr(
        .x = Folder, .y = TimePeriod, .f = Chelsa_Extract_Matching,
        Matches = ClimateModels),

      # which climate scenario
      ClimScenario = purrr::map2_chr(
        .x = Folder, .y = TimePeriod, .f = Chelsa_Extract_Matching,
        Matches = ClimateScenarios),

      URL_File = NULL) %>%

    # Extract variable name from the file name
    dplyr::rowwise() %>%

    dplyr::mutate(
      Var = purrr::map_chr(
        .x = File,
        .f = stringr::str_remove_all,
        pattern = glue("_r1i1p1f1_w5e5_|_norm|CHELSA_|V.2.1|_V\\.2\\.1|{TimePeriod}|.{Ext}|{ClimScenario}")),

      Var = purrr::map2_chr(
        .x = Var, .y = ClimModel,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue("{.y}|{tolower(.y)}"))
        }),

      Var = purrr::map2_chr(
        .x = Var, .y = TimePeriod,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue('{.y}|{stringr::str_replace(.y, "-", "_")}'))
        }),
      Var = purrr::map_chr(
        .x = Var, .f = stringr::str_remove_all, pattern = "__|___"),
      Var = purrr::map_chr(
        .x = Var, .f = stringr::str_remove_all, pattern = "^_|_$"),

      DownPath = file.path(DwnPath, File),
      DownCommand = stringr::str_glue('curl "{URL}" -o "{DownPath}" --silent'),

      # Unique name for variable / time combination
      OutName = paste0(
        Var, "_", TimePeriod, "_", ClimModel, "_", ClimScenario),

      OutName = stringr::str_replace(
        string = OutName,
        pattern = "1981-2010_Current_Current",
        replacement = "1981-2010_Current")) %>%

    dplyr::ungroup() %>%
    dplyr::filter(
      # Only tif files
      Ext == "tif",

      # Exclude previously determined list of variables
      stringr::str_detect(string = Var, pattern = Exclude, negate = TRUE),

      # Exclude duplicated files on the Chelsa server
      Folder != "climatologies/2011-2040/UKESM1-0-LL/ssp126") %>%
    dplyr::select(-Folder) %>%
    dplyr::left_join(Chelsa_Vars, by = dplyr::join_by(Var))


  # Download in parallel
  if (Down && DownParallel) {
    ChelsaClimData %>%
      dplyr::pull(DownCommand) %>%
      furrr::future_walk(
        IASDT.R::System, RObj = FALSE,
        .options = furrr::furrr_options(seed = TRUE), .progress = FALSE)
  }

  # Download sequentially
  if (Down && !DownParallel) {
    ChelsaClimData %>%
      dplyr::pull(DownCommand) %>%
      purrr::walk(IASDT.R::System, RObj = FALSE, .progress = FALSE)
  }

  # Save to disk
  save(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.RData"))
  readr::write_csv(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.csv"))
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Process ----
# |---------------------------------------------------| #

#' Processing Chelsa data
#'
#' Processing Chelsa data
#' @name Chelsa_Process
#' @param InputPath Input path (all `*.tif` files in this path will be processed)
#' @param OutPath Output path. Processed objects will be saved to this path (with the same original file name)
#' @param GridFile Path for the `*.RData` file containing the reference grid. This grid will be used as reference grid for projection and the resulted file will be masked to it
#' @param Verbose should the name of the processed file be printed to the console
#' @author Ahmed El-Gabbas
#' @export

Chelsa_Process <- function(InputPath, OutPath, GridFile, Verbose = FALSE) {
  # reference grid layer
  IASDT.R::CatTime("  >>> Loading reference gird")
  GridR <- GridFile %>%
    IASDT.R::LoadAs() %>%
    terra::rast()

  # List of input and output files
  InOut <- InputPath %>%
    list.files(pattern = ".tif$", full.names = TRUE) %>%
    tibble::tibble(In = .) %>%
    dplyr::mutate(
      Out = purrr::map_chr(
        .x = In, stringr::str_replace, pattern = InputPath, replacement = OutPath))

  if (Verbose) IASDT.R::CatTime("  >>> Processing Chelsa files sequentially")

  InOut %>%
    dplyr::mutate(
      Down = purrr::walk2(
        .x = In, .y = Out,
        .f = purrr::safely(~{

          if (Verbose)  CatTime(glue::glue("     --> {basename(.x)}"))

          IASDT.R::LoadPackages(dplyr, raster, terra)

          Rstr <- .x %>%
            # read tif file as terra rast object
            terra::rast() %>%
            # project to reference grid
            terra::project(GridR, method = "average", threads = TRUE) %>%
            # convert back to raster object
            raster::raster() %>%
            raster::mask(IASDT.R::LoadAs(GridFile))

          # Ensure that the object is located in memory, not reading from temporary file
          # This may not be necessary as we save the file as .tif file not .RData
          if (raster::fromDisk(Rstr)) {
            Rstr <- raster::readAll(Rstr)
          }

          invisible(gc())

          # Write file to disk
          terra::writeRaster(x = terra::rast(Rstr), filename = .y, overwrite = TRUE)
        }), .progress = FALSE))

  # Check if all files were processed successfully

  OutFilesExist <- all(file.exists(InOut$Out))
  OutFilesOkay <- InOut$Out %>%
    purrr::map_lgl(CheckTiff) %>%
    all()

  if (all(OutFilesExist, OutFilesOkay)) {
    IASDT.R::CatTime(stringr::str_glue("  >>>  All CHELSA files ({nrow(InOut)} files) were processed successfully"))
  }

  if (!OutFilesExist) {
    MissingFiles <- InOut$Out %>%
      purrr::discard(~file.exists(.x)) %>%
      basename()
    IASDT.R::CatTime(stringr::str_glue("  >>>  {length(MissingFiles)} files were not processed"))
    cat(stringr::str_glue("     --> {MissingFiles}"), sep = "\n")
  }

  if (!OutFilesOkay) {
    CorruptFiles <- InOut$Out %>%
      purrr::map_lgl(CheckTiff) %>%
      which() %>%
      magrittr::not()
    CorruptFiles <- InOut$Out[CorruptFiles] %>%
      basename()
    IASDT.R::CatTime(stringr::str_glue("  >>>  {length(CorruptFiles)} files are corrupted"))
    cat(stringr::str_glue("     --> {CorruptFiles}"), sep = "\n")
  }

  return(invisible(NULL))
}
