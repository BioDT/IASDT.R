## |------------------------------------------------------------------------| #
# CHELSA_Process ----
## |------------------------------------------------------------------------| #

#' Process CHELSA Climate Data
#'
#' This function processes CHELSA climate data, with an option to download them
#' first. It processes each variable to the European scale and reference grid,
#' and outputs in TIFF and NetCDF (NC) formats. It also saves grouped data for
#' each of the 46 climate scenarios.
#'
#' @name CHELSA_Process
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to 8.
#' @param Download Logical, whether to download the CHELSA files. Defaults to
#'   `FALSE`.
#' @param Download_Attempts Integer. The maximum number of download attempts.
#'   Defaults to 10.
#' @param Sleep Integer. Time in seconds to wait after each download attempt.
#'   Defaults to 5.
#' @param BaseURL String, the base URL for downloading CHELSA climate data.
#' @param Download_NCores integer. Number of CPU cores to use for parallel
#'   downloading of CHELSA data. Only valid if Download = `TRUE`. Defaults to 4.
#' @param Overwrite Logical, whether to re-download files that already exist.
#'   Defaults to `FALSE`.
#' @param CompressLevel Integer; the level of compression to use when saving
#'   NetCDF files. Default is 5.
#' @param OverwriteProcessed Logical; whether to overwrite already processed
#'   files. Default is `FALSE`.
#' @inheritParams CHELSA_Prepare
#' @author Ahmed El-Gabbas
#' @export

CHELSA_Process <- function(
    FromHPC = TRUE, EnvFile = ".env", NCores = 8, Download = FALSE,
    Overwrite = FALSE, Download_Attempts = 10, Sleep = 5, OtherVars = "",
    BaseURL = paste0(
      "https://os.zhdk.cloud.switch.ch/envicloud/",
      "chelsa/chelsa_V2/GLOBAL/"),
    Download_NCores = 4, CompressLevel = 5, OverwriteProcessed = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c("EnvFile", "BaseURL")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c("FromHPC", "Download", "Overwrite", "OverwriteProcessed")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  NumericArgs <- c(
    "Download_NCores", "CompressLevel", "Sleep", "Download_Attempts", "NCores")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs, envir = environment())

  if (NCores < 1 || Download_NCores < 1) {
    stop("`NCores` must be a positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_CHELSA_In <- Path_CHELSA_Out <- Path_Out_NC <- TimePeriod <-
    ClimateModel <- ClimateScenario <- Path_Out_tif <- Name <- FilePath <-
    File_List <- Processed_Maps <- Path_Down <-
    InputOkay <- AllOkay <- Process <- Failed <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  IASDT.R::CatTime("Environment variables")

  if (!file.exists(EnvFile)) {
    stop(
      paste0("Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CHELSA_In", "DP_R_CHELSA_Input", FALSE, FALSE,
      "Path_CHELSA_Out", "DP_R_CHELSA_Output", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CHELSA_In", "DP_R_CHELSA_Input_Local", FALSE, FALSE,
      "Path_CHELSA_Out", "DP_R_CHELSA_Output_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # Ensure that the output path exists
  fs::dir_create(
    c(
      Path_CHELSA_In, Path_CHELSA_Out,
      file.path(Path_CHELSA_Out, c("Tif", "NC", "Processed"))))

  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  # Prepare CHELSA metadata / download CHELSA data -----

  IASDT.R::CatTime("Prepare CHELSA metadata / download CHELSA data")
  TimePrepare <- lubridate::now(tzone = "CET")

  # 19 Bioclimatic variables (+ OtherVars, if not empty string) * 46 CC options
  CHELSA_Data <- IASDT.R::CHELSA_Prepare(
    EnvFile = EnvFile, FromHPC = FromHPC, Download = Download,
    NCores = Download_NCores, Overwrite = Overwrite, OtherVars = OtherVars,
    BaseURL = BaseURL, Download_Attempts = Download_Attempts, Sleep = Sleep)

  IASDT.R::CatDiff(
    InitTime = TimePrepare, Level = 1,
    Prefix = "Prepare CHELSA metadata / download CHELSA data was finished in ")

  # # ..................................................................... ###

  # Check input CHELSA files -----

  if (isFALSE(Download)) {

    # If Download = TRUE, there is no need to re-check the files as these files
    # were already checked while downloading the files

    if (NCores == 1) {
      IASDT.R::CatTime("Check input CHELSA files sequentially")
      future::plan("future::sequential", gc = TRUE)
    } else {
      IASDT.R::CatTime("Check input CHELSA files on parallel")
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    CHELSA_Data_Checked <- CHELSA_Data %>%
      dplyr::select(Path_Down) %>%
      dplyr::mutate(
        InputOkay = furrr::future_map_lgl(
          .x = Path_Down,
          .f = ~ (file.exists(.x) && IASDT.R::CheckTiff(.x)),
          .options = furrr::furrr_options(seed = TRUE, packages = "IASDT.R"))
      ) %>%
      dplyr::filter(isFALSE(InputOkay))

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

    if (nrow(CHELSA_Data_Checked) > 0) {
      readr::write_lines(
        x = dplyr::pull(CHELSA_Data_Checked, "InputOkay"),
        file = file.path(Path_CHELSA_Out, "ProblematicTiffs.txt"))

      stop(
        paste0(
          "Not all input tiff files are available and valid.",
          "Check `ProblematicTiffs.txt`"),
        call. = FALSE
      )
    }

    # CHELSA files that will not be processed
    Diff <- setdiff(
      list.files(Path_CHELSA_In, pattern = "CHELSA.+.tif$", full.names = TRUE),
      CHELSA_Data$Path_Down)

    if (length(Diff) > 0) {
      message(
        paste0(
          " >> Only Bioclimatic variables and variables identified in ", "`OtherVars`, if any, will be processed (",
          nrow(CHELSA_Data), " files)\n >> ", length(Diff),
          " files will not be processed.\n",
          " >> See `NotProcessed.txt` for the list of files"))

      readr::write_lines(
        x = Diff, file = file.path(Path_CHELSA_Out, "NotProcessed.txt"))
    }

    rm(AllOkay, Diff, envir = environment())
  }

  # # ..................................................................... ###

  # Processing CHELSA data -----

  TimeProcess <- lubridate::now(tzone = "CET")

  if (NCores == 1) {
    IASDT.R::CatTime("Processing CHELSA maps sequentially")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    IASDT.R::CatTime("Processing CHELSA maps on parallel")
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # Exclude previously processed files (after checking)
  if (OverwriteProcessed) {
    CHELSA2Process <- CHELSA_Data %>%
      dplyr::select(Path_Down, Path_Out_NC, Path_Out_tif)
  } else {
    CHELSA2Process <- CHELSA_Data %>%
      dplyr::select(Path_Down, Path_Out_NC, Path_Out_tif) %>%
      dplyr::mutate(
        Process = furrr::future_map2_lgl(
          .x = Path_Out_NC, .y = Path_Out_tif,
          .f = ~ {
            NC_Okay <- file.exists(.x) && IASDT.R::CheckTiff(.x)
            Tif_Okay <- file.exists(.y) && IASDT.R::CheckTiff(.y)
            return(isFALSE(NC_Okay && Tif_Okay))
          },
          .options = furrr::furrr_options(seed = TRUE, packages = "IASDT.R"))
      ) %>%
      dplyr::filter(Process) %>%
      dplyr::select(-"Process")
  }

  # Processing CHELSA files

  if (nrow(CHELSA2Process) > 0) {
    CHELSA2Process <- CHELSA2Process %>%
      dplyr::mutate(
        Failed = furrr::future_pmap_lgl(
          .l = list(Path_Down, Path_Out_NC, Path_Out_tif),
          .f = function(Down = Path_Down, Path_Out_NC, Path_Out_tif) {
            FileMetadata <- dplyr::filter(CHELSA_Data, Path_Down == Down)

            # Set `GTIFF_SRS_SOURCE` configuration option to EPSG to use
            # official parameters (overriding the ones from GeoTIFF keys)
            # see: https://stackoverflow.com/questions/78007307
            terra::setGDALconfig("GTIFF_SRS_SOURCE", "EPSG")

            Attempt <- 0
            while (TRUE) {
              Attempt <- Attempt + 1
              Try <- try(
                IASDT.R::CHELSA_Project(
                  Metadata = FileMetadata, EnvFile = EnvFile, FromHPC = FromHPC,
                  CompressLevel = CompressLevel, ReturnMap = FALSE),
                silent = TRUE)

              if (!inherits(Try, "try-error") || Attempt >= 5) {
                break
              }

              if (inherits(Try, "try-error")) {
                # re-download the file if it fails to be processed
                system(
                  command = FileMetadata$DownCommand,
                  ignore.stdout = TRUE, ignore.stderr = TRUE)
              }
            }

            invisible(gc())

            if (inherits(Try, "try-error")) {
              return(TRUE)
            } else {
              if (IASDT.R::CheckTiff(Path_Out_tif) &&
                  IASDT.R::CheckTiff(Path_Out_NC)) {
                return(FALSE)
              } else {
                return(TRUE)
              }
            }
          },
          .options = furrr::furrr_options(
            seed = TRUE,
            packages = c("dplyr", "terra", "IASDT.R", "tibble", "ncdf4"),
            globals = c("CHELSA_Data", "EnvFile", "FromHPC", "CompressLevel"))
        )
      ) %>%
      dplyr::filter(Failed)

    if (nrow(CHELSA2Process) > 0) {
      readr::write_lines(
        x = CHELSA2Process$Path_Down,
        file = file.path(Path_CHELSA_Out, "FailedProcessing.txt"))

      stop(
        paste(
          "\n >> ", nrow(CHELSA2Process), " files failed to process.\n",
          " >> Check `FailedProcessing.txt` for more details"),
        call. = FALSE)
    }

    IASDT.R::CatTime("All tiff files were processed", Level = 1)

  } else {

    IASDT.R::CatTime("All tiff files were already processed.", Level = 1)

  }

  rm(CHELSA2Process, envir = environment())

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  IASDT.R::CatDiff(
    InitTime = TimeProcess,
    Prefix = "Processing CHELSA data was finished in ", Level = 1)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Group CHELSA data by time, climate model/scenario ----

  if (NCores == 1) {
    IASDT.R::CatTime(
      "Group CHELSA data by time and climate model/scenario sequentially")
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    IASDT.R::CatTime(
      "Group CHELSA data by time and climate model/scenario on parallel")
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # String to be matched to extract variable names
  SelectedVars <- c("bio", "Bio", OtherVars) %>%
    # Only keep non-empty strings. If `OtherVars` = "", only bioclimatic
    # variables will be processed.
    stringr::str_subset(".+") %>%
    # Combine the strings into a single string separated by "|". This matches
    # any variable starting with "bio" or "Bio" or any of the characters in
    # `OtherVars` up to the first occurrence of underscore "_".
    paste0(sep = ".*?_", collapse = "|")

  CHELSA_Processed <- CHELSA_Data %>%
    dplyr::select(TimePeriod, ClimateModel, ClimateScenario, Path_Out_tif) %>%
    dplyr::slice(gtools::mixedorder(Path_Out_tif)) %>%
    dplyr::summarise(
      File_List = list(Path_Out_tif),
      .by = c(TimePeriod, ClimateModel, ClimateScenario)) %>%

    dplyr::mutate(
      Name = paste0(
        "R_", TimePeriod, "_", ClimateModel, "_", ClimateScenario),
      Name = stringr::str_replace(
        Name, "1981-2010_Current_Current", "Current"),
      Name = stringr::str_replace_all(Name, "-", "_"),

      FilePath = file.path(
        Path_CHELSA_Out, "Processed", paste0(Name, ".RData")),

      Processed_Maps = furrr::future_pmap(
        .l = list(File_List, FilePath, Name),
        .f = function(File_List, FilePath, Name) {

          MapNames <- basename(File_List) %>%
            # Extract variable names
            stringr::str_extract(SelectedVars) %>%
            # remove the last underscore
            stringr::str_remove_all("_$")

          Map <- terra::rast(File_List) %>%
            stats::setNames(MapNames) %>%
            terra::subset(gtools::mixedsort(MapNames)) %>%
            IASDT.R::setRastCRS() %>%
            IASDT.R::setRastVals() %>%
            terra::wrap()

          # save to disk
          IASDT.R::SaveAs(InObj = Map, OutObj = Name, OutPath = FilePath)

          return(Map)

        },
        .options = furrr::furrr_options(
          seed = TRUE, packages = c("terra", "IASDT.R", "stringr"),
          globals = c("CHELSA_Data", "SelectedVars"))))

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  save(
    CHELSA_Processed,
    file = file.path(Path_CHELSA_Out, "CHELSA_Processed.RData")
  )

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  CHELSA_Processed_DT <- dplyr::select(CHELSA_Processed, -Processed_Maps)

  save(
    CHELSA_Processed_DT,
    file = file.path(Path_CHELSA_Out, "CHELSA_Processed_DT.RData"))

  readr::write_csv(
    x = CHELSA_Processed_DT,
    file = file.path(Path_CHELSA_Out, "CHELSA_Processed_DT.csv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nProcessing CHELSA data took ")

  return(invisible(NULL))
}
