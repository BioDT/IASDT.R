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
    n_search = 100L) {

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::check_args(
    args_to_check = c("kingdom", "phylum"), args_type = "character")
  ecokit::check_args(args_to_check = "n_search", args_type = "numeric")

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "EASIN_URL", "DP_R_easin_taxa_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  Kingdom <- Phylum <- NULL

  # # ..................................................................... ###

  ## Download EASIN taxa -----

  # Use the API; loop until reaching the end of the list
  # Only search for plant species

  time_start_taxa <- lubridate::now(tzone = "CET")

  # start at the first species: skip = 0
  Skip <- 0
  # current step ID
  ID <- 0
  # a list to store the taxonomy list
  EASIN_taxa <- list()

  repeat {
    ID <- ID + 1

    # The API has been changed in April 2023; the main change: using
    # 'apixg/catxg' instead of 'api/cat'
    URL <- stringr::str_glue(
      "{EASIN_URL}/kingdom/{kingdom}/skip/{Skip}/take/{n_search}")

    # Extract species data as tibble
    taxa_data <- try(RCurl::getURL(URL, .mapUnicode = FALSE), silent = TRUE)
    if (inherits(taxa_data, "try-error")) {
      break
    }

    taxa_data <- dplyr::tibble(jsonlite::fromJSON(taxa_data, flatten = TRUE))

    # If there is no data, break the loop
    if (nrow(taxa_data) == 0) {
      break
    }

    EASIN_taxa[[ID]] <- taxa_data
    Skip <- Skip + n_search

    # If the number of rows of the data < n_search, break the loop
    if (nrow(taxa_data) < n_search) {
      break
    }
    rm(taxa_data, URL, envir = environment())
  }

  # Merging data ----
  EASIN_taxa <- dplyr::bind_rows(EASIN_taxa) %>%
    # Only keep vascular plants; Although I searched only for plant taxa, there
    # are still some non-vascular plant species, even for non-plant species.
    dplyr::filter(Kingdom == !!kingdom, Phylum == !!phylum)


  ecokit::cat_diff(
    init_time = time_start_taxa,
    prefix = "Extracting EASIN taxonomy was finished in ", level = 2L)

  return(EASIN_taxa)
}
