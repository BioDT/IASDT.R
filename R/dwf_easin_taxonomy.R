## |------------------------------------------------------------------------| #
# easin_taxonomy ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name easin_data
#' @rdname easin_data
#' @order 2

easin_taxonomy <- function(
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
    "easin_url", "DP_R_easin_taxa_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  Kingdom <- Phylum <- NULL

  # # ..................................................................... ###

  ## Download easin taxa -----

  # Use the API; loop until reaching the end of the list
  # Only search for plant species

  time_start_taxa <- lubridate::now(tzone = "CET")

  # start at the first species: skip = 0
  skip <- 0
  # current step ID
  id <- 0
  # a list to store the taxonomy list
  easin_taxa <- list()

  repeat {
    id <- id + 1

    # The API has been changed in April 2023; the main change: using
    # 'apixg/catxg' instead of 'api/cat'
    url <- stringr::str_glue(
      "{easin_url}/kingdom/{kingdom}/skip/{skip}/take/{n_search}")

    # Extract species data as tibble
    taxa_data <- try(RCurl::getURL(url, .mapUnicode = FALSE), silent = TRUE)
    if (inherits(taxa_data, "try-error")) {
      break
    }

    taxa_data <- dplyr::tibble(jsonlite::fromJSON(taxa_data, flatten = TRUE))

    # If there is no data, break the loop
    if (nrow(taxa_data) == 0) {
      break
    }

    easin_taxa[[id]] <- taxa_data
    skip <- skip + n_search

    # If the number of rows of the data < n_search, break the loop
    if (nrow(taxa_data) < n_search) {
      break
    }
    rm(taxa_data, url, envir = environment())
  }

  # Merging data ----
  easin_taxa <- dplyr::bind_rows(easin_taxa) %>%
    # Only keep vascular plants; Although I searched only for plant taxa, there
    # are still some non-vascular plant species, even for non-plant species.
    dplyr::filter(Kingdom == !!kingdom, Phylum == !!phylum)


  ecokit::cat_diff(
    init_time = time_start_taxa,
    prefix = "Extracting easin taxonomy was finished in ", level = 2L)

  return(easin_taxa)
}
