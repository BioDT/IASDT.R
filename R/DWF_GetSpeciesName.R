## |------------------------------------------------------------------------| #
# GetSpeciesName ----
## |------------------------------------------------------------------------| #

#' Get species name or information of an `IAS-pDT` species ID
#'
#' This function retrieves detailed information on `IAS-pDT` species list,
#' optionally filtered by a specific IAS-pDT species ID (`SpeciesID`).
#'
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param SpeciesID optional IASDT species ID for which detailed information is
#'   required. If not provided, the function returns the entire species list.
#' @name GetSpeciesName
#' @author Ahmed El-Gabbas
#' @return A data frame containing species information. If a species ID
#'   `SpeciesID` is provided, it only returns species information for the listed
#'   species, otherwise return the full list of IAS.
#' @export

GetSpeciesName <- function(SpeciesID = NULL, EnvFile = ".env") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- IAS_ID <- TaxaInfoFile <- NULL

  # Load environment variables

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "TaxaInfoFile", "DP_R_Taxa_info", FALSE, TRUE,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Reading species info

  # using encoding = "UTF-8" to keep non-ascii characters

  SpNames <- utils::read.delim(TaxaInfoFile, sep = "\t", encoding = "UTF-8") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      IAS_ID = paste0("Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4)))

  if (is.null(SpeciesID)) {
    return(SpNames)
  } else {

    if (is.numeric(SpeciesID)) {
      SpeciesID <- paste0(
        "Sp_", stringr::str_pad(SpeciesID, pad = "0", width = 4))
    }

    NGridCells <- IASDT.R::Path(Path_PA, "Sp_PA_Summary_DF.RData")

    if (!file.exists(NGridCells)) {
      stop(
        "Sp_PA_Summary_DF.RData file does not exist in the ", Path_PA,
        " folder", call. = FALSE)
    }

    NGridCells <- IASDT.R::LoadAs(NGridCells) %>%
      dplyr::mutate(
        IAS_ID = paste0(
          "Sp_", stringr::str_pad(IAS_ID, pad = "0", width = 4))) %>%
      dplyr::filter(IAS_ID == SpeciesID) %>%
      dplyr::pull("NCells_All")

    Out <- dplyr::filter(SpNames, IAS_ID == SpeciesID) %>%
      dplyr::mutate(NCells = NGridCells)

    return(Out)
  }
}
