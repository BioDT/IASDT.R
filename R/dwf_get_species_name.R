## |------------------------------------------------------------------------| #
# get_species_name ----
## |------------------------------------------------------------------------| #

#' Get species name or information of an `IASDT` species ID
#'
#' This function retrieves detailed information on `IASDT` species list,
#' optionally filtered by a specific `IASDT` species ID (`species_id`).
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param species_id optional IASDT species ID for which detailed information is
#'   required. If not provided, the function returns the entire species list.
#' @name get_species_name
#' @author Ahmed El-Gabbas
#' @return A data frame containing species information. If a species ID
#'   `species_id` is provided, it only returns species information for the
#'   listed species, otherwise return the full list of IAS.
#' @export

get_species_name <- function(species_id = NULL, env_file = ".env") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_pa <- ias_id <- taxa_info_file <- NULL

  # Load environment variables

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "taxa_info_file", "DP_R_taxa_info", FALSE, TRUE,
    "path_pa", "DP_R_pa", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # Reading species info

  # using encoding = "UTF-8" to keep non-ascii characters

  species_names <- utils::read.delim(
    taxa_info_file, sep = "\t", encoding = "UTF-8") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      ias_id = paste0("sp_", stringr::str_pad(ias_id, pad = "0", width = 4)))

  if (is.null(species_id)) {
    return(species_names)
  } else {

    if (is.numeric(species_id)) {
      species_id <- paste0(
        "sp_", stringr::str_pad(species_id, pad = "0", width = 4))
    }

    n_grid_cells <- fs::path(path_pa, "sp_pa_summary_df.RData")

    if (!file.exists(n_grid_cells)) {
      ecokit::stop_ctx(
        paste0(
          "`sp_pa_summary_df.RData` file does not exist in the ", path_pa,
          " folder"),
        n_grid_cells = n_grid_cells, path_pa = path_pa,
        include_backtrace = TRUE)
    }

    n_grid_cells <- ecokit::load_as(n_grid_cells) %>%
      dplyr::mutate(
        ias_id = paste0(
          "sp_", stringr::str_pad(ias_id, pad = "0", width = 4))) %>%
      dplyr::filter(ias_id == species_id) %>%
      dplyr::pull("n_cells_all")

    out <- dplyr::filter(species_names, ias_id == species_id) %>%
      dplyr::mutate(n_cells = n_grid_cells)

    return(out)
  }
}
