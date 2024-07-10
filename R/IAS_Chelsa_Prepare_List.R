## |------------------------------------------------------------------------| #
# Chelsa_Prepare_List ----
## |------------------------------------------------------------------------| #

#' Prepare list of potential variables
#'
#' Prepare list of potential variables
#'
#' @name Chelsa_Prepare_List
#' @param Down Logical. Download Chelsa files?
#' @param DownParallel Logical. Download input files on parallel (if `Down` = `TRUE`)
#' @param DwnPath String. Download path
#' @param OutPath String. Output path
#' @param UpdateExisting Logical. Reprocess existing file?
#' @param Path_Chelsa String. Path for Chelsa analyses
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
#' @details
#' list of variables exist under current and future climates.
#'
#' 46 variables available at 46 options (current and 45 future scenarios)

Chelsa_Prepare_List <- function(
    Down = FALSE, DownParallel = TRUE, DwnPath = NULL, OutPath = NULL,
    UpdateExisting = FALSE, Path_Chelsa = "Data/Chelsa") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Variable <- DownFile <- NULL


  BaseURL <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/"
  ClimateModels <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  ClimateScenarios <- c("ssp126", "ssp370", "ssp585")

  # Variables to exclude
  Exclude <- c(
    "ai", "hurs", "clt", "sfcWind", "vpd", "rsds", "pet", "cmi", "swb",
    "pr_", "tasmax_", "tasmin_", "tas_") %>%
    stringr::str_c(collapse = "|")

  # List of Chelsa Vars we are interested in (available in current/future)
  ChelsaVars <- IASDT.R::Chelsa_Vars() %>%
    dplyr::select(-"scale", -"offset") %>%
    dplyr::rename(Var = Variable)

  ChelsaClimData <- file.path(Path_Chelsa, "DwnLinks") %>%
    # files containing download links for climatology data
    list.files(
      pattern = "DwnLinks_Climatologies_.+txt$", recursive = TRUE,
      full.names = TRUE) %>%
    dplyr::tibble(URL_File = .) %>%
    dplyr::group_by(.data$URL_File) %>%

    # Add download links
    dplyr::mutate(
      URL = purrr::map(
        .x = .data$URL_File,
        .f = ~{
          .x %>%
            readr::read_lines() %>%
            trimws() # remove trailing spaces
        }),
      URL_File = basename(.data$URL_File)) %>%
    tidyr::unnest(cols = "URL") %>%
    dplyr::ungroup() %>%

    dplyr::mutate(
      # The name of the downloaded file and folder
      Folder = stringr::str_remove_all(string = .data$URL, pattern = BaseURL),
      File = basename(.data$Folder),
      Folder = dirname(.data$Folder),

      # Extract time period
      TimePeriod = stringr::str_remove_all(
        string = .data$URL_File, pattern = "DwnLinks_Climatologies_|.txt"),

      # File extension
      Ext = tools::file_ext(.data$URL),

      # which climate model
      ClimModel = purrr::map2_chr(
        .x = .data$Folder, .y = .data$TimePeriod,
        .f = Chelsa_Extract_Matching, Matches = ClimateModels),

      # which climate scenario
      ClimScenario = purrr::map2_chr(
        .x = .data$Folder, .y = .data$TimePeriod,
        .f = Chelsa_Extract_Matching, Matches = ClimateScenarios),

      URL_File = NULL) %>%

    # Extract variable name from the file name
    dplyr::rowwise() %>%

    dplyr::mutate(
      Var = purrr::map_chr(
        .x = .data$File,
        .f = stringr::str_remove_all,
        pattern = stringr::str_glue("_r1i1p1f1_w5e5_|_norm|CHELSA_|V.2.1|_V\\.2\\.1|{TimePeriod}|.{Ext}|{ClimScenario}")),

      Var = purrr::map2_chr(
        .x = .data$Var, .y = .data$ClimModel,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue("{.y}|{tolower(.y)}"))
        }),

      Var = purrr::map2_chr(
        .x = .data$Var, .y = .data$TimePeriod,
        .f = ~{
          .x %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue('{.y}|{stringr::str_replace(.y, "-", "_")}'))
        }),
      Var = purrr::map_chr(
        .x = .data$Var, .f = stringr::str_remove_all, pattern = "__|___"),
      Var = purrr::map_chr(
        .x = .data$Var, .f = stringr::str_remove_all, pattern = "^_|_$"),

      DownFile = file.path(DwnPath, .data$File),

      OutFile = purrr::map_chr(
        .x = .data$DownFile, stringr::str_replace,
        pattern = DwnPath, replacement = OutPath),

      DownCommand = stringr::str_glue('curl -k -L "{URL}" -o "{DownFile}" --silent'),

      # Unique name for variable / time combination
      OutName = paste0(
        .data$Var, "_", .data$TimePeriod, "_",
        .data$ClimModel, "_", .data$ClimScenario),

      OutName = stringr::str_replace(
        string = .data$OutName,
        pattern = "1981-2010_Current_Current",
        replacement = "1981-2010_Current")) %>%

    dplyr::ungroup() %>%
    dplyr::filter(
      # Only tif files
      .data$Ext == "tif",

      # Exclude previously determined list of variables
      stringr::str_detect(string = .data$Var, pattern = Exclude, negate = TRUE),

      # Exclude duplicated files on the Chelsa server
      .data$Folder != "climatologies/2011-2040/UKESM1-0-LL/ssp126") %>%
    dplyr::select(-"Folder") %>%
    dplyr::left_join(ChelsaVars, by = "Var")


  if (Down) {

    if (UpdateExisting) {
      Data2Down <- ChelsaClimData
    } else {
      Data2Down <- dplyr::filter(ChelsaClimData, magrittr::not(file.exists(DownFile)))
    }

    # Download in parallel
    if (DownParallel && nrow(Data2Down) > 0) {
      Data2Down %>%
        dplyr::pull(.data$DownCommand) %>%
        furrr::future_walk(
          IASDT.R::System, RObj = FALSE,
          .options = furrr::furrr_options(seed = TRUE), .progress = FALSE)
    }

    # Download sequentially
    if (!DownParallel && nrow(Data2Down) > 0) {
      Data2Down %>%
        dplyr::pull(.data$DownCommand) %>%
        purrr::walk(IASDT.R::System, RObj = FALSE, .progress = FALSE)
    }

  }

  # Save to disk
  save(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.RData"))
  readr::write_csv(ChelsaClimData, file = file.path(Path_Chelsa, "ChelsaClimData.csv"))

  return(ChelsaClimData)
}
