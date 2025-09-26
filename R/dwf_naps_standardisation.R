# # |------------------------------------------------------------------------| #
# naps_standardisation ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name naps_data
#' @rdname naps_data
#' @order 4

naps_standardisation <- function(env_file = ".env") {

  # ----------------------------------
  # Important Note:
  # ----------------------------------

  # This function should be called only once per workflow version, as it assign
  # an ID column per standardized species `ias_id` and this ID should not be
  # changed in the same workflow run. If a new standardisation file is used, new
  # `ias_id` will be assigned.

  reason_to_exclude <- taxon_name <- path_taxa_stand <- species_name <-
    species_name2 <- species_file <- path_taxa_info <- ias_id <-
    path_taxa_info_rdata <- class <- order <- family <- NULL

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_taxa_info", "DP_R_taxa_info", FALSE, FALSE,
    "path_taxa_info_rdata", "DP_R_taxa_info_rdata", FALSE, FALSE,
    "path_taxa_stand", "DP_R_taxa_stand", FALSE, TRUE)

  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)

  rm(env_vars_to_read, envir = environment())
  invisible(gc())

  # Read taxa info

  # Extract unique combination of taxon_name, class, order, and family and
  # assign a unique ID (`ias_id`) for each standardised species
  species_id <- readRDS(path_taxa_stand) %>%
    dplyr::filter(is.na(reason_to_exclude)) %>%
    dplyr::distinct(taxon_name, class, order, family) %>%
    dplyr::arrange(class, order, family, taxon_name) %>%
    dplyr::mutate(
      species_name = stringr::word(taxon_name, 1, 2),
      species_name2 = ecokit::replace_space(species_name),
      # species_file will be used for object and file names [avoid times and -
      # symbols]
      species_file = stringr::str_replace_all(species_name2, "\u00D7", "x"),
      species_file = stringr::str_replace_all(species_file, "-", ""),
      genus = stringr::word(taxon_name, 1),
      species = stringr::word(taxon_name, 2),
      ias_id = seq_len(dplyr::n())) %>%
    dplyr::select(tidyselect::all_of("ias_id"), tidyselect::everything())

  # Save species_id as a tab-separated file
  readr::write_tsv(x = species_id, file = path_taxa_info, progress = FALSE)

  selected_cols <- c("ias_id", "taxon_name")
  taxa_list <- readRDS(path_taxa_stand) %>%
    dplyr::filter(is.na(reason_to_exclude)) %>%
    dplyr::select(tidyselect::all_of(c("taxon_name", "speciesKey"))) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      species_name = stringr::word(taxon_name, 1, 2),
      species_name2 = ecokit::replace_space(species_name),
      species_file = stringr::str_replace_all(species_name2, "\u00D7", "x"),
      species_file = stringr::str_replace_all(species_file, "-", "")) %>%
    dplyr::left_join(
      y = dplyr::select(species_id, tidyselect::all_of(selected_cols)),
      by = "taxon_name") %>%
    dplyr::select(
      tidyselect::all_of(selected_cols),
      tidyselect::starts_with("species"), tidyselect::everything()) %>%
    dplyr::arrange(ias_id)

  save(taxa_list, file = path_taxa_info_rdata)

  return(invisible(NULL))

}
