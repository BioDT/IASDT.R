# # |------------------------------------------------------------------------| #
# eLTER_process ----
## |------------------------------------------------------------------------| #

#' Process eLTER data for the `IASDT`
#'
#' This function processes pre-cleaned and pre-standardized Integrated European
#' Long-Term Ecosystem, critical zone and socio-ecological Research
#' ([eLTER](https://elter-ri.eu/)) data for the Invasive Alien Species Digital
#' Twin (`IASDT`).
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param start_year Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be processed. Default is `1981`, which
#'   matches the year ranges of CHELSA current climate data.
#' @return Returns `NULL` invisibly after saving the processed data.
#' @author Ahmed El-Gabbas
#' @name eLTER_process
#' @note This function processes pre-cleaned vascular plants data from eLTER
#'   sites, harmonized by Ahmed El-Gabbas. The original eLTER biodiversity data
#'   were highly heterogeneous in format and structure, requiring
#'   standardization and cleaning before use. Taxonomic standardization with the
#'   GBIF backbone was performed by Marina Golivets (Feb. 2024).
#' @export

eLTER_process <- function(env_file = ".env", start_year = 1981) {

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::check_args(args_to_check = "start_year", args_type = "numeric")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  taxa_list <- eLTER_DT <- Path_PA <- taxon_name <- speciesKey <- Year <- NULL

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "taxa_list", "DP_R_taxa_info_rdata", FALSE, TRUE,
    "Path_PA", "DP_R_PA", FALSE, FALSE,
    "eLTER_DT", "DP_R_elter_raw", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  fs::dir_create(Path_PA)
  taxa_list <- ecokit::load_as(taxa_list)

  eLTER_IAS <- readRDS(eLTER_DT) %>%
    dplyr::select(
      tidyselect::all_of(c(
        "SITE_CODE", "Lon", "Lat", "TaxaName", "Day", "Month", "Year",
        "Source", "ID", "name_to_match", "speciesKey"))) %>%
    dplyr::filter(!is.na(speciesKey)) %>%
    dplyr::left_join(taxa_list, by = "speciesKey") %>%
    dplyr::filter(!is.na(taxon_name), Year >= start_year) %>%
    sf::st_as_sf(
      coords = c("Lon", "Lat"), crs = sf::st_crs(4326), remove = FALSE) %>%
    sf::st_transform(crs = sf::st_crs(3035))

  save(eLTER_IAS, file = fs::path(Path_PA, "eLTER_IAS.RData"))

  return(invisible(NULL))
}
