# |---------------------------------------------------| #
# match_to_gbif.fn ----
# |---------------------------------------------------| #

#' Match taxonomy with GBIF; may return >1 match
#'
#'  Match taxonomy with GBIF; may return >1 match
#' @param taxon_name taxonomy name
#' @param taxon_id taxonomy ID
#' @param include_genus include matches at genus level; default: `FALSE`
#' @name match_to_gbif.fn
#' @importFrom rlang .data
#' @author Marina Golivets
#' @return NULL
#' @export
#' @details
#' as input, provide a vector of verbatim taxon names (preferably with authorship) and a vector of existing local identifiers for those names

match_to_gbif.fn <- function(taxon_name, taxon_id, include_genus = FALSE) {

  # perform initial matching in parallel
  no_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(no_cores)
  all_matches <- pbapply::pblapply(
    taxon_name, rgbif::name_backbone_verbose,
    kingdom = "plants", strict = TRUE, cl = cl)
  parallel::stopCluster(cl)

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
      cbind, .data,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::filter(!is.na(.data$usageKey)) %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$phylum == "Tracheophyta")

  # retrieve best matches
  best_matches <- lapply(all_matches, function(x) x$data) %>%
    mapply(
      cbind, .data,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$phylum == "Tracheophyta")

  matched <- best_matches %>%
    dplyr::filter(!(.data$matchType %in% c("NONE", "HIGHERRANK")))

  nonmatched <- best_matches %>%
    dplyr::filter(.data$matchType %in% c("NONE", "HIGHERRANK"))

  matched_alternative <- try(
    alternative_matches %>%
      dplyr::filter(.data$phylum == "Tracheophyta") %>% # use only vascular plants
      dplyr::filter(.data$confidence >= 0) %>%
      # dplyr::filter(taxon_id %in% nonmatched$taxon_id)
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
  if(nrow(accepted) > 0) {
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
    dplyr::full_join(taxon_list) %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::filter(.data$status == "SYNONYM")
  if(nrow(synonyms) > 0) {
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
    dplyr::full_join(taxon_list) %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$status == "DOUBTFUL")
  if(nrow(doubtful) > 0) {
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
    dplyr::select(-.data$has_accepted) %>%
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

Get_EASIN_Data <- function(SpKey, NSearch = 500) {
  Skip <- 0         # Skip = 0; start at the first presence grid
  ID <- 0           # iteration ID
  DT <- list()      # List object to save the data

  # Looping over data chunks
  while (TRUE) {
    ID <- ID + 1
    URL <- "https://alien.jrc.ec.europa.eu/apixg/geoxg/speciesid/{SpKey}/layertype/grid/skip/{Skip}/take/{NSearch}" %>%
      stringr::str_glue()

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
    do.call(DT, what = "rbind") %>%
      return()
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
#' @return NULL
#' @export
#' @details
#' A function to extract EASIN data for a given EASIN_ID

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

GetAcceptedName <- function(ID) {
  if(!is.na(ID)) {
    rgbif::name_usage(ID)$data$scientificName
  } else {
    NA_character_
  }
}
