## |------------------------------------------------------------------------| #
# GetSpeciesName ----
## |------------------------------------------------------------------------| #

#' Get Species Name or Information Based on IASDT Species ID
#'
#' This function retrieves detailed information on IASDT Species, optionally
#' filtered by a specific IASDT species ID. It first loads environment variables
#' from a given file, then reads a species list from a path specified in the
#' environment variables. If a species ID `SpID` is provided, it only returns
#' species information for the listed species, otherwise return the full list of
#' IAS.
#' @param EnvFile EnvFile A string specifying the path to the environment
#'   variables file. Default is ".env". This file must exist in the specified
#'   path; otherwise, the function will stop with an error message.
#' @param SpID optional IASDT species ID for which detailed information is
#'   required. If not provided, the function returns the entire species list.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @name GetSpeciesName
#' @author Ahmed El-Gabbas
#' @return A data frame containing species information. If a species ID `SpID`
#'   is provided, it only returns species information for the listed species,
#'   otherwise return the full list of IAS.
#' @details The function reads the following environment variables:
#'   - **`DP_R_TaxaInfo`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_TaxaInfo_Local`** (if `FromHPC` = `FALSE`) for the location of the `Species_List_ID.txt` file containing species information.
#'    - **`DP_R_PA`** (if `FromHPC` = `TRUE`) or **`DP_R_PA_Local`** (if `FromHPC` = `FALSE`). The function reads the contents of the `Sp_PA_Summary_DF.RData` file from this path.
#' @export

GetSpeciesName <- function(EnvFile = ".env", SpID = NULL, FromHPC = TRUE) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- IAS_ID <- NCells <- TaxaInfoFile <- NULL

  # Load environment variables
  if (!file.exists(EnvFile)) {
    stop(paste0("Path for environment variables: ", EnvFile, " was not found"))
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo", FALSE, TRUE,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo_Local", FALSE, TRUE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # Reading species info
   SpNames <- utils::read.delim(TaxaInfoFile, sep = "\t") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      IAS_ID = paste0("Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4)))

  if (is.null(SpID)) {
    return(SpNames)
  } else {

    NGridCells <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")

    if (!file.exists(NGridCells)) {
      stop(
        paste0(
          "Sp_PA_Summary_DF.RData file does not exist in the ",
          Path_PA, " folder"))
    }

    NGridCells <- IASDT.R::LoadAs(NGridCells) %>%
      dplyr::mutate(
        IAS_ID = paste0(
          "Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4))) %>%
      dplyr::filter(IAS_ID == SpID) %>%
      dplyr::pull(NCells)

    dplyr::filter(SpNames, IAS_ID == SpID) %>%
      dplyr::mutate(NCells = NGridCells) %>%
      return()
  }
}
