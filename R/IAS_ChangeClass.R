## |------------------------------------------------------------------------| #
# ChangeClass ----
## |------------------------------------------------------------------------| #

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
