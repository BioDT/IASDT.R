## |------------------------------------------------------------------------| #
# Chelsa_Info ----
## |------------------------------------------------------------------------| #

#' Extracts and summarizes information from CHELSA climate model file names
#'
#' This function parses CHELSA climate model file names to extract and summarize
#' information such as the time period, climate model, climate scenario, and
#' variable represented in the file. It returns a data frame with these details
#' for each file name provided.
#' @name Chelsa_Info
#' @param FileName character vector; a vector of CHELSA file names or paths.
#'   Each element should be a string representing the file name or path of a
#'   CHELSA climate model output file.
#' @return A data frame with columns for the original file name (`FileName`),
#'   the climate variable (`Variable`), the time period (`TimePeriod`), the
#'   climate model (`ClimModel`), and the climate scenario (`ClimScenario`).
#'   Additional variable information is joined from an external source
#'   ([IASDT.R::Chelsa_Vars`).
#' @author Ahmed El-Gabbas
#' @examples
#'BaseURL <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/"
#'LFiles <- paste0(
#'  BaseURL,
#'  c(
#'    #'"2041-2070/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio14_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1.tif",
#'    "2011-2040/MRI-ESM2-0/ssp126/bio/CHELSA_bio1_2011-2040_mri-esm2-0_ssp126_V.2.1.tif",
#'    "2011-2040/IPSL-CM6A-LR/ssp370/bio/CHELSA_lgd_2011-2040_ipsl-cm6a-lr_ssp370_V.2.1.tif",
#'    "2041-2070/MRI-ESM2-0/ssp370/bio/CHELSA_bio5_2041-2070_mri-esm2-0_ssp370_V.2.1.tif",
#'    "2011-2040/MRI-ESM2-0/ssp585/bio/CHELSA_scd_2011-2040_mri-esm2-0_ssp585_V.2.1.tif"))
#'
#' Chelsa_Info(FileName = LFiles[1])
#'
#' Chelsa_Info(FileName = LFiles)
#' @export

Chelsa_Info <- function(FileName) {

  # Validate inputs
  if (is.null(FileName)) {
    stop("FileName cannot be NULL", .call = FALSE)
  }

  # Ensure FileName is a character vector
  if (!is.character(FileName)) {
    stop("FileName must be a character vector", .call = FALSE)
  }

  purrr::map_dfr(
    .x = FileName,
    .f = ~{

      file_base <- basename(.x)

      TimePeriod <- dplyr::case_when(
        stringr::str_detect(file_base, "1981-2010") ~ "Current",
        stringr::str_detect(file_base, "2011-2040") ~ "2011_2040",
        stringr::str_detect(file_base, "2041-2070") ~ "2041_2070",
        stringr::str_detect(file_base, "2071-2100") ~ "2071_2100",
        .default = NA_character_)

      ClimModel <- dplyr::case_when(
        stringr::str_detect(file_base, "1981-2010") ~ "Current",
        # National Oceanic and Atmospheric Administration, Geophysical Fluid
        # Dynamics Laboratory, Princeton, NJ 08540, USA
        stringr::str_detect(file_base, "GFDL-ESM4|gfdl-esm4") ~ "GFDL_ESM4",
        # Institut Pierre Simon Laplace, Paris 75252, France
        stringr::str_detect(file_base, "IPSL-CM6A-LR|ipsl-cm6a-lr") ~ "IPSL_CM6A",
        # Max Planck Institute for Meteorology, Hamburg 20146, Germany
        stringr::str_detect(file_base, "MPI-ESM1-2-HR|mpi-esm1-2-hr") ~ "MPI_ESM",
        # Met Office Hadley Centre, Fitzroy Road, Exeter, Devon, EX1 3PB, UK
        stringr::str_detect(file_base, "UKESM1-0-LL|ukesm1-0-ll") ~ "UKESM",
        # Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan
        stringr::str_detect(file_base, "MRI-ESM2-0|mri-esm2-0") ~ "MRI_ESM2",
        .default = NA_character_)

      ClimScenario <- dplyr::case_when(
        stringr::str_detect(file_base, "1981-2010") ~ "Current",
        # SSP1-RCP2.6 climate as simulated by the GCMs
        stringr::str_detect(file_base, "ssp126") ~ "ssp126",
        # ssp370 SSP3-RCP7 climate as simulated by the GCMs
        stringr::str_detect(file_base, "ssp370") ~ "ssp370",
        # ssp585 SSP5-RCP8.5 climate as simulated by the GCMs
        stringr::str_detect(file_base, "ssp585") ~ "ssp585",
        .default = NA_character_)

      CurrVar <-  file_base %>%
        stringr::str_remove_all(
          (c("1981-2010", "2011-2040", "2041-2070", "2071-2100",
             "GFDL-ESM4|gfdl-esm4", "IPSL-CM6A-LR|ipsl-cm6a-lr",
             "MPI-ESM1-2-HR|mpi-esm1-2-hr",
             "UKESM1-0-LL|ukesm1-0-ll", "MRI-ESM2-0|mri-esm2-0",
             "ssp126", "ssp370", "ssp585", "CHELSA") %>%
             paste0(collapse = "|") %>%
             paste0("|V.2.1.tif"))) %>%
        stringr::str_remove_all("_")

      tibble::tibble(
        FileName =  file_base, Variable = CurrVar, TimePeriod = TimePeriod,
        ClimModel = ClimModel, ClimScenario = ClimScenario)
    }) %>%
    dplyr::left_join(IASDT.R::Chelsa_Vars(), by = "Variable") %>%
    return()
}
