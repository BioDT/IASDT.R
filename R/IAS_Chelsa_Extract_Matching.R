# |---------------------------------------------------| #
# Chelsa_Extract_Matching ----
# |---------------------------------------------------| #

#' Extract time / climate models & scenarios from Chelsa URL
#'
#' A function to extract which climate model and scenario each file represent
#' @name Chelsa_Extract_Matching
#' @param String URL string
#' @param Time Time period to match
#' @param Matches Further string to match
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
