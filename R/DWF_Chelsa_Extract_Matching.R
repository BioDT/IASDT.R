## |------------------------------------------------------------------------| #
# Chelsa_Extract_Matching ----
## |------------------------------------------------------------------------| #

#' Extract time / climate models & scenarios from Chelsa URL
#'
#' This function extracts which climate model and scenario each file represents
#' based on the provided URL string, time period, and additional match criteria.
#' It is specifically designed to work with Chelsa climate data URLs.
#' @param String A character string representing the URL from which to extract
#'   information.
#' @param Time A character string representing the time period to match. The
#'   special value "1981-2010" is treated as representing current climates, and
#'   the function returns "Current" for this case.
#' @param Matches A character vector of additional strings to match within the
#'   URL. The function returns these matches if they are found in the URL.
#' @author Ahmed El-Gabbas
#' @return A character vector of matched strings if any matches are found;
#'   otherwise, returns `NA_character_`. If the time period is "1981-2010", the
#'   function returns "Current".
#' @export
#' @name Chelsa_Extract_Matching
#' @details A function to extract which climate model and scenario each file
#' represent

Chelsa_Extract_Matching <- function(String, Time, Matches) {

  # Validate inputs
  if (is.null(String) || is.null(Time) || is.null(Matches)) {
    stop("String, Time, and Matches cannot be NULL", call. = FALSE)
  }

  # assign "Current" for files represent current climates
  if (Time == "1981-2010") {
    return("Current")
  }

  # Index of matched text using sapply and which
  Which <- sapply(
    X = Matches,
    FUN = function(x) {
      grepl(x, String, ignore.case = TRUE)
    }) %>%
    which()

  # Return matched strings or NA_character_ if no match found
  if (length(Which) > 0) {
    return(Matches[Which])
  } else {
    return(NA_character_)
  }
}
