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
