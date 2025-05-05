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
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character", args_to_check = "env_file")
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "start_year")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  TaxaList <- eLTER_DT <- Path_PA <- taxon_name <- speciesKey <- Year <- NULL

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "TaxaList", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "Path_PA", "DP_R_PA", FALSE, FALSE,
    "eLTER_DT", "DP_R_eLTER_raw", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  fs::dir_create(Path_PA)
  TaxaList <- ecokit::load_as(TaxaList)

  eLTER_IAS <- readRDS(eLTER_DT) %>%
    dplyr::select(
      tidyselect::all_of(c(
        "SITE_CODE", "Lon", "Lat", "TaxaName", "Day", "Month", "Year",
        "Source", "ID", "name_to_match", "speciesKey"))) %>%
    dplyr::filter(!is.na(speciesKey)) %>%
    dplyr::left_join(TaxaList, by = "speciesKey") %>%
    dplyr::filter(!is.na(taxon_name), Year >= start_year) %>%
    sf::st_as_sf(
      coords = c("Lon", "Lat"), crs = sf::st_crs(4326), remove = FALSE) %>%
    sf::st_transform(crs = sf::st_crs(3035))

  save(eLTER_IAS, file = fs::path(Path_PA, "eLTER_IAS.RData"))

  return(invisible(NULL))
}
