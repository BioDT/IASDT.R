## |------------------------------------------------------------------------| #
# GetAcceptedName ----
## |------------------------------------------------------------------------| #

#' Get accepted name of a taxon
#'
#' This function queries the GBIF database to retrieve the accepted scientific
#' name for a given taxa based on its GBIF ID. It uses the `rgbif` package to
#' perform the query.
#' @name GetAcceptedName
#' @param ID numeric or integer; the GBIF ID of the taxa for which the accepted
#'   name is being queried.
#' @author Ahmed El-Gabbas
#' @return A character string representing the scientific name of the accepted
#'   taxa. If the ID is valid and found in the GBIF database, the scientific
#'   name is returned. If the ID is NA, the function returns `NA_character_`.
#' @export
#' @details
#' Get accepted name of a taxa
#' @examples
#' GetAcceptedName(5372559) # https://www.gbif.org/species/5372559

GetAcceptedName <- function(ID) {

  if (is.null(ID)) {
    stop("ID cannot be NULL", call. = FALSE)
  }

  if (!is.na(ID)) {
    return(rgbif::name_usage(ID)$data$scientificName)
  } else {
    return(NA_character_)
  }
}
