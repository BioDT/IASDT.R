## |------------------------------------------------------------------------| #
# EASIN_taxonomy ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name EASIN_data
#' @rdname EASIN_data
#' @order 2

EASIN_taxonomy <- function(
    env_file = ".env", kingdom = "Plantae", phylum = "Tracheophyta",
    n_search = 100) {

  # # ..................................................................... ###

  # Checking arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("kingdom", "phylum", "env_file"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_search")

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "EASIN_URL", "DP_R_EASIN_taxa_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Kingdom <- Phylum <- NULL

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
      "{EASIN_URL}/kingdom/{kingdom}/skip/{Skip}/take/{n_search}")

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
    Skip <- Skip + n_search

    # If the number of rows of the data < n_search, break the loop
    if (nrow(Data) < n_search) {
      break
    }
    rm(Data, URL, envir = environment())
  }

  # Merging data ----
  EASIN_Taxa <- dplyr::bind_rows(EASIN_Taxa) %>%
    # Only keep vascular plants; Although I searched only for plant taxa, there
    # are still some non-vascular plant species, even for non-plant species.
    dplyr::filter(Kingdom == !!kingdom, Phylum == !!phylum)


  ecokit::cat_diff(
    init_time = TimeStartTaxa,
    prefix = "Extracting EASIN taxonomy was finished in ", level = 2L)

  return(EASIN_Taxa)
}
