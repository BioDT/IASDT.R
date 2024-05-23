## |------------------------------------------------------------------------| #
# Chelsa_Info ----
## |------------------------------------------------------------------------| #

#' Extract info from Chelsa file name
#'
#' Extract info from Chelsa file name
#' @name Chelsa_Info
#' @param FileName character; URL, file path or file name
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
#' Chelsa_Info(LFiles[1])
#'
#' Chelsa_Info(LFiles)
#' @export

Chelsa_Info <- function(FileName) {

  FileName %>%
    basename() %>%
    purrr::map_dfr(
      .f = ~{
        TimePeriod <- dplyr::case_when(
          stringr::str_detect(.x, "1981-2010") ~ "Current",
          stringr::str_detect(.x, "2011-2040") ~ "2011_2040",
          stringr::str_detect(.x, "2041-2070") ~ "2041_2070",
          stringr::str_detect(.x, "2071-2100") ~ "2071_2100",
          .default = NA_character_)

        ClimModel <- dplyr::case_when(
          stringr::str_detect(.x, "1981-2010") ~ "Current",
          # National Oceanic and Atmospheric Administration, Geophysical Fluid Dynamics Laboratory, Princeton, NJ 08540, USA
          stringr::str_detect(.x, "GFDL-ESM4|gfdl-esm4") ~ "GFDL_ESM4",
          # Institut Pierre Simon Laplace, Paris 75252, France
          stringr::str_detect(.x, "IPSL-CM6A-LR|ipsl-cm6a-lr") ~ "IPSL_CM6A",
          # Max Planck Institute for Meteorology, Hamburg 20146, Germany
          stringr::str_detect(.x, "MPI-ESM1-2-HR|mpi-esm1-2-hr") ~ "MPI_ESM",
          # Met Office Hadley Centre, Fitzroy Road, Exeter, Devon, EX1 3PB, UK
          stringr::str_detect(.x, "UKESM1-0-LL|ukesm1-0-ll") ~ "UKESM",
          # Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan
          stringr::str_detect(.x, "MRI-ESM2-0|mri-esm2-0") ~ "MRI_ESM2",
          .default = NA_character_)

        ClimScenario <- dplyr::case_when(
          stringr::str_detect(.x, "1981-2010") ~ "Current",
          # SSP1-RCP2.6 climate as simulated by the GCMs
          stringr::str_detect(.x, "ssp126") ~ "ssp126",
          # ssp370 SSP3-RCP7 climate as simulated by the GCMs
          stringr::str_detect(.x, "ssp370") ~ "ssp370",
          # ssp585 SSP5-RCP8.5 climate as simulated by the GCMs
          stringr::str_detect(.x, "ssp585") ~ "ssp585",
          .default = NA_character_)

        CurrVar <-  .x %>%
          stringr::str_remove_all(
            c("1981-2010", "2011-2040", "2041-2070", "2071-2100", "GFDL-ESM4|gfdl-esm4",
              "IPSL-CM6A-LR|ipsl-cm6a-lr", "MPI-ESM1-2-HR|mpi-esm1-2-hr",
              "UKESM1-0-LL|ukesm1-0-ll", "MRI-ESM2-0|mri-esm2-0",
              "ssp126", "ssp370", "ssp585", "CHELSA") %>%
              paste0(collapse = "|", sep = "_") %>%
              paste0("|V.2.1.tif")) %>%
          stringr::str_remove_all("_")

        tibble::tibble(
          FileName =  .x, Variable = CurrVar, TimePeriod = TimePeriod,
          ClimModel = ClimModel, ClimScenario = ClimScenario)
      }) %>%
    dplyr::left_join(IASDT.R::Chelsa_Vars(), by = "Variable")
}
