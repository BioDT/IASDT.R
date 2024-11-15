## |------------------------------------------------------------------------| #
# EASIN_Taxonomy ----
## |------------------------------------------------------------------------| #

#' Extract EASIN Taxonomy Data
#'
#' This function retrieves taxonomy data from the EASIN database for vascular
#' plants. It is also possible to get similar data for other taxonomic groups.
#' It handles pagination automatically, downloading data in chunks until all
#' available data for the specified taxa are retrieved.
#'
#' @param BaseURL character; the base URL for accessing the EASIN taxonomy API.
#'   Default is "https://easin.jrc.ec.europa.eu/apixg/catxg".
#' @param Kingdom character; the taxonomic kingdom to search within. Default is
#'   "Plantae".
#' @param Phylum character; the taxonomic phylum to filter within the specified
#'   kingdom. Default is "Tracheophyta".
#' @param NSearch integer; the number of records to attempt to retrieve per
#'   request. Default is 100.
#' @return A tibble containing the retrieved taxonomy data for the specified
#'   kingdom and phylum.
#' @author Ahmed El-Gabbas and Marina Golivets
#' @export
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [EASIN_Process] function.
#' @details This function loops through the EASIN API, retrieving data in chunks
#'   until all available data for the specified kingdom and phylum are
#'   collected. The results are returned as a tibble after filtering for the
#'   specified taxa.

EASIN_Taxonomy <- function(
    BaseURL = "https://easin.jrc.ec.europa.eu/apixg/catxg",
    Kingdom = "Plantae", Phylum = "Tracheophyta", NSearch = 100) {

  # # ..................................................................... ###

  # Checking arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("BaseURL", "Kingdom", "Phylum"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NSearch")

  # # ..................................................................... ###

  ## Download EASIN taxa -----

  # Use the API; loop until reaching the end of the list
  # Only search for plant species

  TimeStartTaxa <- lubridate::now(tzone = "CET")

  # start at the first species: skip = 0
  Skip <- 0
  # current step ID
  ID <- 0
  # a list to store the taxonomy list
  EASIN_Taxa <- list()

  while (TRUE) {
    ID <- ID + 1

    # The API has been changed in April 2023; the main change: using
    # 'apixg/catxg' instead of 'api/cat'
    URL <- stringr::str_glue(
      "{BaseURL}/kingdom/{Kingdom}/skip/{Skip}/take/{NSearch}")

    # Extract species data as tibble
    Data <- try(RCurl::getURL(URL, .mapUnicode = FALSE), silent = TRUE)
    if (inherits(Data, "try-error")) {
      break
    }

    Data <- dplyr::tibble(jsonlite::fromJSON(Data, flatten = TRUE))

    if (nrow(Data) == 0) {
      # If there is no data, break the loop
      break
    }

    EASIN_Taxa[[ID]] <- Data
    Skip <- Skip + NSearch

    # If the number of rows of the data < NSearch, break the loop
    if (nrow(Data) < NSearch) {
      break
    }
    rm(Data, URL, envir = environment())
  }

  # Merging data ----
  EASIN_Taxa <- dplyr::bind_rows(EASIN_Taxa) %>%
    # Only keep vascular plants; Although I searched only for plant taxa, there
    # are still some non-vascular plant species, even for non-plant species.
    dplyr::filter(Kingdom == !!Kingdom, Phylum == !!Phylum)


  IASDT.R::CatDiff(
    InitTime = TimeStartTaxa,
    Prefix = "Extracting EASIN taxonomy was finished in ", Level = 2)

  return(EASIN_Taxa)
}
