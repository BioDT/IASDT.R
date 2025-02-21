# # |------------------------------------------------------------------------| #
# elTER_Process ----
## |------------------------------------------------------------------------| #

#' eLTER data processing function
#'
#' This function processes pre-cleaned and pre-standardized eLTER data.
#' @param FromHPC Logical. Whether the processing is being done on an 
#'   High-Performance Computing (HPC) environment, to adjust file paths 
#'   accordingly. Default: `TRUE`.
#' @param EnvFile Character. Path to the environment file containing paths to 
#'   data sources. Defaults to `.env`.
#' @param StartYear Numeric. The starting year for the occurrence data. Only
#'   records from this year onward will be processed. Default is `1981`, which
#'   matches the year ranges of CHELSA current climate data.
#' @return Returns `NULL` invisibly after saving the processed data.
#' @author Ahmed El-Gabbas & Marina Golivets
#' @name elTER_Process
#' @export

elTER_Process <- function(
    FromHPC = TRUE, EnvFile = ".env", StartYear = 1981) {

  # # ..................................................................... ###

  # Checking arguments ----
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "StartYear")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  TaxaList <- eLTER_DT <- Path_PA <- taxon_name <- speciesKey <- Year <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaList", "DP_R_TaxaInfo_RData", FALSE, TRUE,
      "Path_PA", "DP_R_PA", FALSE, FALSE,
      "eLTER_DT", "DP_R_eLTER_In", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaList", "DP_R_TaxaInfo_RData_Local", FALSE, TRUE,
      "Path_PA", "DP_R_PA_Local", FALSE, FALSE,
      "eLTER_DT", "DP_R_eLTER_In_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  fs::dir_create(Path_PA)
  TaxaList <- IASDT.R::LoadAs(TaxaList)

  eLTER_IAS <- readRDS(eLTER_DT) %>%
    dplyr::select(
      tidyselect::all_of(c(
        "SITE_CODE", "Lon", "Lat", "TaxaName", "Day", "Month", "Year",
        "Source", "ID", "name_to_match", "speciesKey"))) %>%
    dplyr::filter(!is.na(speciesKey)) %>%
    dplyr::left_join(TaxaList, by = "speciesKey") %>%
    dplyr::filter(!is.na(taxon_name), Year >= StartYear) %>%
    sf::st_as_sf(
      coords = c("Lon", "Lat"), crs = sf::st_crs(4326), remove = FALSE) %>%
    sf::st_transform(crs = sf::st_crs(3035))

  save(eLTER_IAS, file = IASDT.R::Path(Path_PA, "eLTER_IAS.RData"))

  return(invisible(NULL))
}
