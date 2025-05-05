## |------------------------------------------------------------------------| #
# get_species_name ----
## |------------------------------------------------------------------------| #

#' Get species name or information of an `IASDT` species ID
#'
#' This function retrieves detailed information on `IASDT` species list,
#' optionally filtered by a specific `IASDT` species ID (`species_ID`).
#'
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param species_ID optional IASDT species ID for which detailed information is
#'   required. If not provided, the function returns the entire species list.
#' @name get_species_name
#' @author Ahmed El-Gabbas
#' @return A data frame containing species information. If a species ID
#'   `species_ID` is provided, it only returns species information for the
#'   listed species, otherwise return the full list of IAS.
#' @export

get_species_name <- function(species_ID = NULL, env_file = ".env") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- IAS_ID <- TaxaInfoFile <- NULL

  # Load environment variables

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "TaxaInfoFile", "DP_R_Taxa_info", FALSE, TRUE,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Reading species info

  # using encoding = "UTF-8" to keep non-ascii characters

  SpNames <- utils::read.delim(TaxaInfoFile, sep = "\t", encoding = "UTF-8") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      IAS_ID = paste0("Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4)))

  if (is.null(species_ID)) {
    return(SpNames)
  } else {

    if (is.numeric(species_ID)) {
      species_ID <- paste0(
        "Sp_", stringr::str_pad(species_ID, pad = "0", width = 4))
    }

    NGridCells <- fs::path(Path_PA, "Sp_PA_Summary_DF.RData")

    if (!file.exists(NGridCells)) {
      ecokit::stop_ctx(
        paste0(
          "`Sp_PA_Summary_DF.RData` file does not exist in the ", Path_PA,
          " folder"),
        NGridCells = NGridCells, Path_PA = Path_PA)
    }

    NGridCells <- ecokit::load_as(NGridCells) %>%
      dplyr::mutate(
        IAS_ID = paste0(
          "Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4))) %>%
      dplyr::filter(IAS_ID == species_ID) %>%
      dplyr::pull("NCells_All")

    Out <- dplyr::filter(SpNames, IAS_ID == species_ID) %>%
      dplyr::mutate(NCells = NGridCells)

    return(Out)
  }
}
