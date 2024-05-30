## |------------------------------------------------------------------------| #
# GetSpeciesName ----
## |------------------------------------------------------------------------| #

#' Get species information for specific species number
#'
#' Get species information for specific species number
#'
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param SpID Species ID
#' @name GetSpeciesName
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetSpeciesName <- function(EnvFile = ".env", SpID = NULL) {

  IAS_ID <- NCells <- NULL

  # Load environment variables
  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
  } else {
    MSG <- paste0("Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }

  SpNames <- Sys.getenv("DP_R_Mod_Path_TaxaList") %>%
    file.path("Species_List_ID.txt") %>%
    utils::read.delim(sep = "\t") %>%
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
