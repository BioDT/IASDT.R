## |------------------------------------------------------------------------| #
# ChangeClass ----
## |------------------------------------------------------------------------| #

#' Change the class of specified columns in a GBIF data frame
#'
#' This function takes a data frame containing GBIF data and converts the data
#' types of specific columns to more appropriate types for analysis. Integer,
#' integer64, double, logical, and date classes are applied to various columns
#' as specified.
#' @param DF A data frame containing GBIF data.
#' @name ChangeClass
#' @author Ahmed El-Gabbas
#' @return A data frame with the same data as `DF` but with specified columns
#'   converted to their designated data types.
#' @note This function is currently neither exported or used by other functions.
#' @keywords internal
#' @noRd

ChangeClass <- function(DF) {

  if (is.null(DF)) {
    stop("DF cannot be NULL", call. = FALSE)
  }

  VarsInt <- c(
    "year", "month", "day", "coordinateUncertaintyInMeters",
    "acceptedNameUsageID", "taxonKey", "acceptedTaxonKey", "kingdomKey",
    "phylumKey", "classKey", "orderKey", "familyKey", "genusKey", "speciesKey")
  VarsInt64 <- c("gbifID", "catalogNumber", "recordNumber", "taxonID",
    "identificationID")
  VarsDbl <- c("decimalLatitude", "decimalLongitude")
  VarsLgl <- c("hasCoordinate", "hasGeospatialIssues", "repatriated")
  VarsDte <- c("modified", "eventDate", "dateIdentified", "lastInterpreted",
    "lastParsed", "lastCrawled")

  DF %>%
    dplyr::mutate(dplyr::across(VarsInt, as.integer)) %>%
    dplyr::mutate(dplyr::across(VarsInt64, bit64::as.integer64)) %>%
    dplyr::mutate(dplyr::across(VarsDbl, as.double))  %>%
    dplyr::mutate(dplyr::across(VarsLgl, as.logical))  %>%
    dplyr::mutate(dplyr::across(VarsDte, lubridate::as_date)) %>%
    return()
}
