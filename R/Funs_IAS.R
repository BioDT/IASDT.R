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
    dplyr::filter((if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta")

  # retrieve best matches
  best_matches <- lapply(all_matches, function(x) x$data) %>%
    mapply(
      cbind, .,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::distinct() %>%
    dplyr::filter((if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta")

  matched <- best_matches %>%
    dplyr::filter(!(.data$matchType %in% c("NONE", "HIGHERRANK")))

  matched_alternative <- try(
    alternative_matches %>%
      # use only vascular plants
      # filter only if phylum column exists
      dplyr::filter((if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta") %>%
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
    URL <- stringr::str_glue("https://easin.jrc.ec.europa.eu/apixg/geoxg/speciesid/{SpKey}/layertype/grid/skip/{Skip}/take/{NSearch}")

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
#' @param Down logical; Download Chelsa files?
#' @param DownParallel logical; Download input files on parallel (if `Down` = `TRUE`)
#' @param DwnPath Download path
#' @param OutPath Output path
#' @param UpdateExisting logical; Reprocess existing file?
#' @param Path_Chelsa Path for Chelsa analyses (including `Chelsa_Vars.txt` file)
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
#' @details
#' list of variables exist under current and future climates.
#'
#' 46 variables available at 46 options (current and 45 future scenarios)

Chelsa_Prepare_List <- function(
    Down = FALSE, DownParallel = TRUE, DwnPath = NULL, OutPath = NULL,
    UpdateExisting = FALSE, Path_Chelsa = "Data/Chelsa") {

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
    readr::read_delim(delim = "\t", show_col_types = FALSE) %>%
    dplyr::select(-"scale", -"offset")

  ChelsaClimData <- file.path(Path_Chelsa, "DwnLinks") %>%
    # files containing download links for climatology data
    list.files(
      pattern = "DwnLinks_Climatologies_.+txt$", recursive = TRUE,
      full.names = TRUE) %>%
    dplyr::tibble(URL_File = .) %>%
    dplyr::group_by(.data$URL_File) %>%

    # Add download links
    dplyr::mutate(
      URL = purrr::map(
        .x = .data$URL_File,
        .f = ~{
          .x %>%
            readr::read_lines() %>%
            trimws() # remove trailing spaces
        }),
      URL_File = basename(.data$URL_File)) %>%
    tidyr::unnest(cols = "URL") %>%
    dplyr::ungroup() %>%

    dplyr::mutate(
      # The name of the downloaded file and folder
      Folder = stringr::str_remove_all(string = .data$URL, pattern = BaseURL),
      File = basename(.data$Folder),
      Folder = dirname(.data$Folder),

      # Extract time period
      TimePeriod = stringr::str_remove_all(
        string = .data$URL_File, pattern = "DwnLinks_Climatologies_|.txt"),

      # File extension
      Ext = tools::file_ext(.data$URL),

      # which climate model
      ClimModel = purrr::map2_chr(
        .x = .data$Folder, .y = .data$TimePeriod,
        .f = Chelsa_Extract_Matching, Matches = ClimateModels),

      # which climate scenario
      ClimScenario = purrr::map2_chr(
        .x = .data$Folder, .y = .data$TimePeriod,
        .f = Chelsa_Extract_Matching, Matches = ClimateScenarios),

      URL_File = NULL) %>%

    # Extract variable name from the file name
    dplyr::rowwise() %>%

    dplyr::mutate(
      Var = purrr::map_chr(
        .x = .data$File,
        .f = stringr::str_remove_all,
        pattern = stringr::str_glue("_r1i1p1f1_w5e5_|_norm|CHELSA_|V.2.1|_V\\.2\\.1|{TimePeriod}|.{Ext}|{ClimScenario}")),

      Var = purrr::map2_chr(
        .x = .data$Var, .y = .data$ClimModel,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue("{.y}|{tolower(.y)}"))
        }),

      Var = purrr::map2_chr(
        .x = .data$Var, .y = .data$TimePeriod,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue('{.y}|{stringr::str_replace(.y, "-", "_")}'))
        }),
      Var = purrr::map_chr(
        .x = .data$Var, .f = stringr::str_remove_all, pattern = "__|___"),
      Var = purrr::map_chr(
        .x = .data$Var, .f = stringr::str_remove_all, pattern = "^_|_$"),

      DownFile = file.path(DwnPath, .data$File),

      OutFile = purrr::map_chr(
        .x = .data$DownFile, stringr::str_replace,
        pattern = DwnPath, replacement = OutPath),

      DownCommand = stringr::str_glue('curl -k -L "{URL}" -o "{DownFile}" --silent'),

      # Unique name for variable / time combination
      OutName = paste0(
        .data$Var, "_", .data$TimePeriod, "_",
        .data$ClimModel, "_", .data$ClimScenario),

      OutName = stringr::str_replace(
        string = .data$OutName,
        pattern = "1981-2010_Current_Current",
        replacement = "1981-2010_Current")) %>%

    dplyr::ungroup() %>%
    dplyr::filter(
      # Only tif files
      .data$Ext == "tif",

      # Exclude previously determined list of variables
      stringr::str_detect(string = .data$Var, pattern = Exclude, negate = TRUE),

      # Exclude duplicated files on the Chelsa server
      .data$Folder != "climatologies/2011-2040/UKESM1-0-LL/ssp126") %>%
    dplyr::select(-"Folder") %>%
    dplyr::left_join(Chelsa_Vars, by = dplyr::join_by("Var"))


  if (Down) {

    if (UpdateExisting) {
      Data2Down <- ChelsaClimData
    } else {
      Data2Down <- dplyr::filter(ChelsaClimData, magrittr::not(file.exists(DownFile)))
    }

    # Download in parallel
    if (DownParallel && nrow(Data2Down) > 0) {
      Data2Down %>%
        dplyr::pull(.data$DownCommand) %>%
        furrr::future_walk(
          IASDT.R::System, RObj = FALSE,
          .options = furrr::furrr_options(seed = TRUE), .progress = FALSE)
    }

    # Download sequentially
    if (!DownParallel && nrow(Data2Down) > 0) {
      Data2Down %>%
        dplyr::pull(.data$DownCommand) %>%
        purrr::walk(IASDT.R::System, RObj = FALSE, .progress = FALSE)
    }

  }

  # Save to disk
  save(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.RData"))
  readr::write_csv(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.csv"))

  return(ChelsaClimData)
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Vars ----
# |---------------------------------------------------| #

#' Detailed information on Chelsa data
#'
#' Detailed information on Chelsa data
#' @name Chelsa_Vars
#' @author Ahmed El-Gabbas
#' @examples
#' print(Chelsa_Vars(), n = Inf)
#' @export
#' @details
#' https://chelsa-climate.org/bioclim/

Chelsa_Vars <- function() {
  tibble::tribble(
    ~Variable, ~Long_name, ~unit, ~scale, ~offset, ~explanation,
    "bio1", "mean annual air temperature", "deg C", 0.1, -273.15, "mean annual daily mean air temperatures averaged over 1 year",
    "bio2", "mean diurnal air temperature range", "deg C", 0.1, 0, "mean diurnal range of temperatures averaged over 1 year",
    "bio3", "isothermality", "deg C", 0.1, 0, "ratio of diurnal variation to annual variation in temperatures",
    "bio4", "temperature seasonality", "deg C/100", 0.1, 0, "standard deviation of the monthly mean temperatures",
    "bio5", "mean daily maximum air temperature of the warmest month", "deg C", 0.1, -273.15, "The highest temperature of any monthly daily mean maximum temperature",
    "bio6", "mean daily minimum air temperature of the coldest month", "deg C", 0.1, -273.15, "The lowest temperature of any monthly daily mean maximum temperature",
    "bio7", "annual range of air temperature", "deg C", 0.1, 0, "Difference between the Maximum Temperature of Warmest month and the Minimum Temperature of Coldest month (bio5-bio6)",
    "bio8", "mean daily mean air temperatures of the wettest quarter", "deg C", 0.1, -273.15, "The wettest quarter of the year is determined (to the nearest month)",
    "bio9", "mean daily mean air temperatures of the driest quarter", "deg C", 0.1, -273.15, "The driest quarter of the year is determined (to the nearest month)",
    "bio10", "mean daily mean air temperatures of the warmest quarter", "deg C", 0.1, -273.15, "The warmest quarter of the year is determined (to the nearest month)",
    "bio11", "mean daily mean air temperatures of the coldest quarter", "deg C", 0.1, -273.15, "The coldest quarter of the year is determined (to the nearest month)",
    "bio12", "annual precipitation amount", "kg m-2", 0.1, 0, "Accumulated precipitation amount over 1 year",
    "bio13", "precipitation amount of the wettest month", "kg m-2", 0.1, 0, "The precipitation of the wettest month.",
    "bio14", "precipitation amount of the driest month", "kg m-2", 0.1, 0, "The precipitation of the driest month.",
    "bio15", "precipitation seasonality", "kg m-2", 0.1, 0, "The Coefficient of Variation = 100*standard deviation of the monthly precipitation / mean ",
    "bio16", "mean monthly precipitation amount of the wettest quarter", "kg m-2", 0.1, 0, "The wettest quarter of the year is determined (to the nearest month)",
    "bio17", "mean monthly precipitation amount of the driest quarter", "kg m-2", 0.1, 0, "The driest quarter of the year is determined (to the nearest month)",
    "bio18", "mean monthly precipitation amount of the warmest quarter", "kg m-2", 0.1, 0, "The warmest quarter of the year is determined (to the nearest month)",
    "bio19", "mean monthly precipitation amount of the coldest quarter", "kg m-2", 0.1, 0, "The coldest quarter of the year is determined (to the nearest month)",
    "gdgfgd0", "First growing degree day above 0 deg C", "julian day", 1, 0, "First day of the year above 0 deg C",
    "gdgfgd5", "First growing degree day above 5 deg C", "julian day", 1, 0, "First day of the year above 5 deg C",
    "gdgfgd10", "First growing degree day above 10 deg C", "julian day", 1, 0, "First day of the year above 10 deg C",
    "gddlgd0", "Last growing degree day above 0 deg C", "julian day", 1, 0, "Last day of the year above 0 deg C",
    "gddlgd5", "Last growing degree day above 5 deg C", "julian day", 1, 0, "Last day of the year above 5 deg C",
    "gddlgd10", "Last growing degree day above 10 deg C", "julian day", 1, 0, "Last day of the year above 10 deg C",
    "gdd0", "Growing degree days heat sum above 0 deg C", "deg C", 0.1, 0, "heat sum of all days above the 0 deg C temperature accumulated over 1 year.",
    "gdd5", "Growing degree days heat sum above 5 deg C", "deg C", 0.1, 0, "heat sum of all days above the 5 deg C temperature accumulated over 1 year.",
    "gdd10", "Growing degree days heat sum above 10 deg C", "deg C", 0.1, 0, "heat sum of all days above the 10 deg C temperature accumulated over 1 year.",
    "ngd0", "Number of growing degree days", "# days", 1, 0, "Number of days at which tas > 0 deg C",
    "ngd5", "Number of growing degree days", "# days", 1, 0, "Number of days at which tas > 5 deg C",
    "ngd10", "Number of growing degree days", "# days", 1, 0, "Number of days at which tas > 10 deg C",
    "gsl", "growing season length TREELIM", "# days", 1, 0, "Length of the growing season",
    "gst", "Mean temperature of the growing season TREELIM", "deg C", 0.1, -273.15, "Mean temperature of all growing season days based on TREELIM",
    "lgd", "last day of the growing season TREELIM", "julian day", 1, 0, "Last day of the growing season according to TREELIM",
    "fgd", "first day of the growing season TREELIM", "julian day", 1, 0, "first day of the growing season according to TREELIM",
    "gsp", "Accumulated precipiation amount on growing season days TREELIM", "kg m-2", 0.1, 0, "precipitation sum accumulated on all days during the growing season based on TREELIM",
    "fcf", "Frost change frequency", "count", 1, 0, "Number of events in which tmin or tmax go above, or below 0 deg C",
    "scd", "Snow cover days", "count", 1, 0, "Number of days with snowcover calculated using the snowpack model implementation in from TREELIM",
    "swe", "Snow water equivalent", "kg m-2", 0.1, 0, "Amount of luquid water if snow is melted",
    "kg0", "Koeppen-Geiger climate classification", "category", 1, 0, "Koeppen Geiger - Koeppen&Geiger (1936)",
    "kg1", "Koeppen-Geiger climate classification", "category", 1, 0, "Koeppen Geiger without As/Aw differentiation - Koeppen&Geiger (1936)",
    "kg2", "Koeppen-Geiger climate classification", "category", 1, 0, "Koeppen Geiger after Peel et al. 2007",
    "kg3", "Koeppen-Geiger climate classification", "category", 1, 0, "Wissmann 1939",
    "kg4", "Koeppen-Geiger climate classification", "category", 1, 0, "Thornthwaite 1931",
    "kg5", "Koeppen-Geiger climate classification", "category", 1, 0, "Troll-Pfaffen - Troll&Paffen (1964)",
    "npp", "Net primary productivity", "g C m-2 yr-1", 0.1, 0, "Calculated based on the `Miami model`, Lieth, H., 1972."
  )
}


# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Info ----
# |---------------------------------------------------| #

#' Extract info from Chelsa file name
#'
#' Extract info from Chelsa file name
#' @name Chelsa_Info
#' @param FileName character; URL, file path or file name
#' @author Ahmed El-Gabbas
#' @examples
#' LFiles <- c(
#'   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2041-2070/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio14_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1.tif",
#'   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2011-2040/MRI-ESM2-0/ssp126/bio/CHELSA_bio1_2011-2040_mri-esm2-0_ssp126_V.2.1.tif",
#'   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2011-2040/IPSL-CM6A-LR/ssp370/bio/CHELSA_lgd_2011-2040_ipsl-cm6a-lr_ssp370_V.2.1.tif",
#'   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2041-2070/MRI-ESM2-0/ssp370/bio/CHELSA_bio5_2041-2070_mri-esm2-0_ssp370_V.2.1.tif",
#'   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2011-2040/MRI-ESM2-0/ssp585/bio/CHELSA_scd_2011-2040_mri-esm2-0_ssp585_V.2.1.tif")
#'
#' Chelsa_Info(LFiles[1])
#'
#' Chelsa_Info(LFiles)
#' @export

Chelsa_Info <- function(FileName) {

  FileName %>%
    basename() %>%
    purrr::map_dfr(
      .f = ~{
        TimePeriod <- dplyr::case_when(
          stringr::str_detect( .x, "1981-2010") ~ "Current",
          stringr::str_detect( .x, "2011-2040") ~ "2011_2040",
          stringr::str_detect( .x, "2041-2070") ~ "2041_2070",
          stringr::str_detect( .x, "2071-2100") ~ "2071_2100",
          .default = NA_character_)

        ClimModel <- dplyr::case_when(
          stringr::str_detect( .x, "1981-2010") ~ "Current",
          # National Oceanic and Atmospheric Administration, Geophysical Fluid Dynamics Laboratory, Princeton, NJ 08540, USA
          stringr::str_detect( .x, "GFDL-ESM4|gfdl-esm4") ~ "GFDL_ESM4",
          # Institut Pierre Simon Laplace, Paris 75252, France
          stringr::str_detect( .x, "IPSL-CM6A-LR|ipsl-cm6a-lr") ~ "IPSL_CM6A",
          # Max Planck Institute for Meteorology, Hamburg 20146, Germany
          stringr::str_detect( .x, "MPI-ESM1-2-HR|mpi-esm1-2-hr") ~ "MPI_ESM",
          # Met Office Hadley Centre, Fitzroy Road, Exeter, Devon, EX1 3PB, UK
          stringr::str_detect( .x, "UKESM1-0-LL|ukesm1-0-ll") ~ "UKESM",
          # Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan
          stringr::str_detect( .x, "MRI-ESM2-0|mri-esm2-0") ~ "MRI_ESM2",
          .default = NA_character_)

        ClimScenario <- dplyr::case_when(
          stringr::str_detect( .x, "1981-2010") ~ "Current",
          # SSP1-RCP2.6 climate as simulated by the GCMs
          stringr::str_detect( .x, "ssp126") ~ "ssp126",
          # ssp370 SSP3-RCP7 climate as simulated by the GCMs
          stringr::str_detect( .x, "ssp370") ~ "ssp370",
          # ssp585 SSP5-RCP8.5 climate as simulated by the GCMs
          stringr::str_detect( .x, "ssp585") ~ "ssp585",
          .default = NA_character_)

        CurrVar <-  .x %>%
          stringr::str_remove_all(
            c("1981-2010", "2011-2040", "2041-2070", "2071-2100", "GFDL-ESM4|gfdl-esm4",
              "IPSL-CM6A-LR|ipsl-cm6a-lr", "MPI-ESM1-2-HR|mpi-esm1-2-hr",
              "UKESM1-0-LL|ukesm1-0-ll", "MRI-ESM2-0|mri-esm2-0",
              "ssp126", "ssp370", "ssp585", "CHELSA") %>%
              paste0(collapse = "|", sep = "_") %>%
              paste0("|V.2.1.tif")) %>%
          stringr::str_remove_all("_")

        tibble::tibble(
          FileName =  .x, Variable = CurrVar, TimePeriod = TimePeriod,
          ClimModel = ClimModel, ClimScenario = ClimScenario)
      }) %>%
    dplyr::left_join(IASDT.R::Chelsa_Vars(), by = "Variable")
}

# ****************************************************
# ****************************************************

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Chelsa_Project ----
# |---------------------------------------------------| #

#' Project Chelsa data to the study area
#'
#' Project Chelsa data to the study area
#' @name Chelsa_Project
#' @param InputFile Path or URL for input tif file
#' @param OutFile Path for output tif file
#' @param GridFile `raster` or `SpatRaster` for the reference grid. This grid will be used as reference grid for projection and the resulted file will be masked to it
#' @param ReturnMap logical; should the processed map be returned by the end of the function?
#' @param DownPath Where to save downloaded files
#' @param KeepDownloaded if URL is provided as input file, the file will be downloaded to disk first before processing. Should the downloaded file be kept in disk.
#' @param SaveTiff Also save output map as *.tif file
#' @returns if `ReturnMap = TRUE`, the a wrapped SpatRaster object is returned (`PackedSpatRaster`); otherwise nothing is returned. By default the function exports a NetCDF file. If `SaveTiff` is set as `TRUE`, additional tiff file will be saved to disk.
#' @author Ahmed El-Gabbas
#' @export

Chelsa_Project <- function(
    InputFile = NULL, OutFile = NULL, GridFile = NULL, ReturnMap = FALSE,
    DownPath = NULL, KeepDownloaded = TRUE, SaveTiff = FALSE) {

  # Ensure that the reference grid is not null
  if (is.null(GridFile)) stop("GridFile can not be empty")

  # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use official parameters (overriding the ones from GeoTIFF keys)
  terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

  # Input file name
  InputName <- basename(InputFile)

  # Variable description
  VarDesc <- IASDT.R::Chelsa_Info(InputFile)
  # Scale and offset information
  VarScale <- VarDesc$scale
  VarOffset <- VarDesc$offset

  # ||||||||||||||||||||||||||||||||||||||||
  # Ensure that some packages are loaded
  # ||||||||||||||||||||||||||||||||||||||||

  IASDT.R::LoadPackages("dplyr", "raster", "terra") %>%
    suppressMessages() %>%
    suppressWarnings()

  # ||||||||||||||||||||||||||||||||||||||||
  # Loading reference grid
  # ||||||||||||||||||||||||||||||||||||||||

  if (inherits(GridFile, "RasterLayer")) {
    GridR <- terra::rast(GridFile)
  } else {
    if (inherits(GridFile, "PackedSpatRaster")) {
      GridR <- terra::unwrap(GridFile)
    } else {
      GridR <- GridFile
    }
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Remote or local
  # ||||||||||||||||||||||||||||||||||||||||

  Remote <- dplyr::if_else(stringr::str_detect(InputFile, "^http"), TRUE, FALSE)

  if (Remote) {

    if (is.null(OutFile)) stop("OutFile can not be empty if the input file is an URL")

    # Folder to download the file
    if (is.null(DownPath)) {
      # download as temporary file
      DownPath <- tempfile(fileext = ".tif")
      # delete the temporary file after processing
      KeepDownloaded <- FALSE
    } else {
      DownPath <- file.path(DownPath, InputName)
    }

    # Ensure that output dir exists
    IASDT.R::DirCreate(dirname(DownPath), Verbose = FALSE)

    # increase time out to allow more time to download input file
    # This may not be necessary for LUMI HPC
    # EVE has a download limit. This allows more time for download
    options(timeout = max(300, getOption("timeout")))

    # Download file to disk
    utils::download.file(
      url = InputFile, destfile = DownPath, quiet = TRUE, mode = "wb")

  } else {
    # if input file is located locally,
    DownPath <- InputFile
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Land mask
  # ||||||||||||||||||||||||||||||||||||||||

  # Extent to crop the maps prior to processing.
  # This ensures that the object reads from the memory. See below
  CropExtent <- terra::ext(-26, 37.5, 34, 72)

  # source: CHELSA-W5E5 v1.0: W5E5 v1.0 downscaled with CHELSA v2.0
  # https://data.isimip.org/10.48364/ISIMIP.836809.3
  # Version: 1.0.3
  # This file was copied to the package data to make it easier to use it
  LandMaskL <- system.file(
    "extdata", "LandMask.nc", package = "IASDT.R", mustWork = TRUE) %>%
    terra::rast() %>%
    terra::crop(CropExtent) %>%
    suppressWarnings() %>%  # suppress warning on LUMI while cropping
    terra::classify(cbind(0, NA))

  # ||||||||||||||||||||||||||||||||||||||||
  # read tif file as terra SpatRaster object
  # ||||||||||||||||||||||||||||||||||||||||

  # terra package by default considers the scale and offset information stored in the tiff files. Here I disable this to read the raw values as it is and later consider the scale and offset information manually. This is more safe as I found that some of the future projections do not include such information in the tiff files.
  Rstr <- terra::rast(DownPath, raw = TRUE) %>%
    stats::setNames(paste0(tools::file_path_sans_ext(InputName), ".tif")) %>%
    # crop to European boundaries
    # although it is not necessary to crop the input maps into the European boundaries, we will crop the data prior to projection. Cropping will make the values of the raster read from memory not from the file. This is a workaround to avoid wrong extreme values in the output file because of a bug in terra package (see this issue: https://github.com/rspatial/terra/issues/1356) [18.02.2023]
    terra::crop(CropExtent) %>%
    # mask by land mask
    terra::mask(LandMaskL) %>%
    # `gsp` maps contains extremely high values instead of NA; the following replace extreme values with NA
    terra::classify(cbind(420000000, Inf, NA))

  # ||||||||||||||||||||||||||||||||||||||||
  # Manually considering offset and scale
  # ||||||||||||||||||||||||||||||||||||||||

  # For `npp` layers, all tiff maps except for current climate does have a scaling factor
  # all scale and offset information were set manually
  if (VarScale != 1) Rstr <- Rstr * VarScale
  if (VarOffset != 0) Rstr <- Rstr + VarOffset

  # ||||||||||||||||||||||||||||||||||||||||
  # Projecting to reference grid EPSG 3035
  # ||||||||||||||||||||||||||||||||||||||||

  # projection method
  # Use `mode` for `kg` variables (categorical variables)
  Method <- dplyr::if_else(
    stringr::str_detect(names(Rstr), "kg[0-5]"), "mode", "average")

  Rstr <- Rstr %>%
    # project to reference grid
    terra::project(GridR, method = Method, threads = TRUE) %>%
    # mask to the reference grid
    terra::mask(GridR)

  # Ensure that the object is located in memory, not reading from temporary file
  # This may not be necessary as we save the file as .tif file not .RData
  if (magrittr::not(terra::inMemory(Rstr))) {
    terra::values(Rstr) <- terra::values(Rstr)
  }

  terra::crs(Rstr) <- "epsg:3035"

  # ||||||||||||||||||||||||||||||||||||||||
  # Write file to disk --- tiff
  # ||||||||||||||||||||||||||||||||||||||||

  OutFileTif <- dplyr::if_else(
    stringr::str_detect(OutFile, ".tif$"),
    OutFile,
    file.path(dirname(OutFile), paste0(
      tools::file_path_sans_ext(InputName), ".tif")))

  if (SaveTiff) {
    terra::writeRaster(x = Rstr, filename = OutFileTif, overwrite = TRUE)
    if (Remote && magrittr::not(KeepDownloaded)) file.remove(DownPath)
  }

  # ||||||||||||||||||||||||||||||||||||||||
  # Write file to disk --- nc
  # ||||||||||||||||||||||||||||||||||||||||

  OutFileNC <- dplyr::if_else(
    stringr::str_detect(OutFile, ".nc$"),
    OutFile,
    file.path(dirname(OutFile), paste0(tools::file_path_sans_ext(InputName), ".nc")))

  # Variable name of the output *.nc file
  VarName4NC <- c(
    VarDesc$TimePeriod, VarDesc$ClimModel, VarDesc$ClimScenario) %>%
    unique() %>%
    paste0(collapse = "__") %>%
    paste0(VarDesc$Variable, "__", .)

  # global attributes to be added to the *.nc file
  Attrs <- c(
    paste0("Var=", VarDesc$Variable),
    paste0("TimePeriod=", VarDesc$TimePeriod),
    paste0("ClimModel=", VarDesc$ClimModel),
    paste0("ClimScenario=", VarDesc$ClimScenario),
    paste0("Long_name=", VarDesc$Long_name),
    paste0("unit=", VarDesc$unit),
    paste0("explanation=", VarDesc$explanation))

  # save as *.nc file
  terra::writeCDF(
    Rstr, filename = OutFileNC, varname = VarName4NC, unit = VarDesc$unit,
    zname = VarDesc$TimePeriod, atts = Attrs, overwrite = TRUE)

  # ||||||||||||||||||||||||||||||||||||||||
  # Return map?
  # ||||||||||||||||||||||||||||||||||||||||

  if (ReturnMap) {
    return(terra::wrap(Rstr))
  } else {
    return(invisible(NULL))
  }
}
