## |------------------------------------------------------------------------| #
# GetSpeciesName ----

#' Get Species Name or Information Based on IASDT Species ID
#'
#' This function retrieves detailed information on IASDT Species, optionally filtered by a specific IASDT species ID. It first loads environment variables from a given file, then reads a species list from a path specified in the environment variables. If a species ID is provided, it also fetches the number of grid cells associated with that species.
#'
#' @param EnvFile EnvFile A string specifying the path to the environment variables file. Default is ".env". This file must exist in the specified path; otherwise, the function will stop with an error message.
#' @param SpID optional IASDT species ID for which detailed information is required. If not provided, the function returns the entire species list.
#' @name GetSpeciesName
#' @author Ahmed El-Gabbas
#' @return A data frame containing species information.
#' @export

GetSpeciesName <- function(EnvFile = ".env", SpID = NULL) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- NCells <- NULL

  # Load environment variables
  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
  } else {
    stop(paste0("Path for environment variables: ", EnvFile, " was not found"))
  }

  SpNamesF <- Sys.getenv("DP_R_Mod_Path_TaxaList")
  if (SpNamesF == "") {
    stop("DP_R_Mod_Path_TaxaList environment variable not set.")
  }

  SpNamesF <- file.path(SpNamesF, "Species_List_ID.txt")
  if (magrittr::not(file.exists(SpNamesF))) {
    stop("Species_List_ID.txt file does not exist")
  }

  SpNames <- utils::read.delim(SpNamesF, sep = "\t") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
      IAS_ID = paste0("Sp_", IAS_ID))

  if (is.null(SpID)) {
    return(SpNames)
  } else {
    NGridCells <- Sys.getenv("DP_R_Mod_Path_PA") %>%
      file.path("Sp_PA_Summary_DF.RData") %>%
      IASDT.R::LoadAs() %>%
      dplyr::mutate(
        IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
        IAS_ID = paste0("Sp_", IAS_ID)) %>%
      dplyr::filter(IAS_ID == SpID) %>%
      dplyr::pull(NCells)

    SpNames %>%
      dplyr::filter(IAS_ID == SpID) %>%
      dplyr::mutate(NCells = NGridCells) %>%
      return()
  }
}
