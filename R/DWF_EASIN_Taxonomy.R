## |------------------------------------------------------------------------| #
# EASIN_Taxonomy ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 2

EASIN_Taxonomy <- function(
    EnvFile = ".env", Kingdom = "Plantae", Phylum = "Tracheophyta",
    NSearch = 100) {

  # # ..................................................................... ###

  # Checking arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Kingdom", "Phylum", "EnvFile"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NSearch")

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "EASIN_URL", "DP_R_EASIN_taxa_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

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

  repeat {
    ID <- ID + 1

    # The API has been changed in April 2023; the main change: using
    # 'apixg/catxg' instead of 'api/cat'
    URL <- stringr::str_glue(
      "{EASIN_URL}/kingdom/{Kingdom}/skip/{Skip}/take/{NSearch}")

    # Extract species data as tibble
    Data <- try(RCurl::getURL(URL, .mapUnicode = FALSE), silent = TRUE)
    if (inherits(Data, "try-error")) {
      break
    }

    Data <- dplyr::tibble(jsonlite::fromJSON(Data, flatten = TRUE))

    # If there is no data, break the loop
    if (nrow(Data) == 0) {
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
