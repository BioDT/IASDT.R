#' Process CHELSA Climate Data for the `IAS-pDT`
#'
#' Downloads, processes, and projects [CHELSA](https://chelsa-climate.org/)
#' climate data at the European scale for the Invasive Alien Species prototype
#' Digital Twin (`IAS-pDT`). Supports multiple climate scenarios, outputting
#' data in TIFF and NetCDF formats. Orchestrated by `CHELSA_Process()`, with
#' helper functions `CHELSA_Prepare()` and `CHELSA_Project()`.
#'
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param Download Logical. If `TRUE`, downloads CHELSA files. Default: `FALSE`.
#' @param Download_Attempts Integer. Maximum download retries. Default: `10`.
#' @param Sleep Integer. Seconds to wait between download attempts. Default:
#'   `5`.
#' @param Download_NCores Integer. Number of CPU cores to use for parallel
#'   downloading of CHELSA data. Only valid if Download = `TRUE`. Defaults to 4.
#' @param Overwrite Logical. If `TRUE`, re-downloads existing files. Default:
#'   `FALSE`.
#' @param CompressLevel Integer. NetCDF compression level (1 = least, 9 = most).
#'   Default: `5`.
#' @param OverwriteProcessed Logical. If `TRUE`, overwrites processed files.
#'   Default: `FALSE`.
#' @param OtherVars Character. Additional variables to process (e.g., `"npp"`
#'   for Net Primary Productivity alongside 19 bioclimatic variables
#'   bio1-bio19). Use `""` for bioclimatic only. See [CHELSA_Vars] for details.
#'   Default: `"npp"`.
#' @param Metadata Tibble. Single-row metadata for input files, prepared by
#'   `CHELSA_Prepare()`
#' @section Functions details:
#' - **`CHELSA_Process()`**: Main function; optionally downloads CHELSA data,
#'   processes it to the European scale and reference grid, and saves TIFF and
#'   NetCDF outputs for 46 climate scenarios.
#' - **`CHELSA_Prepare()`**: Extracts metadata from local URL text files and
#'   manages optional downloads.
#' - **`CHELSA_Project()`**: Projects data to the IAS-pDT reference grid with
#'   optional transformations.
#'
#' @note
#' - `CHELSA_Prepare()` and `CHELSA_Project()` are internal helpers, not for
#' direct use.
#' - Processes 19 bioclimatic variables (bio1–bio19) plus optional variables
#' (e.g., NPP) for 46 scenarios (1 current, 45 future).
#' - Time-intensive; depends on file size and compute resources.

## |------------------------------------------------------------------------| #
# CHELSA_Process ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 1

CHELSA_Process <- function(
    EnvFile = ".env", NCores = 8L, Download = FALSE, Overwrite = FALSE,
    Download_Attempts = 10L, Sleep = 5L, OtherVars = "npp",
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

  CharArgs <- c("EnvFile", "OtherVars")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c("Download", "Overwrite", "OverwriteProcessed")
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
      "Path to environment variables: ", EnvFile, " was not found",
      call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_CHELSA_In", "DP_R_CHELSA_raw", FALSE, FALSE,
    "Path_CHELSA_Out", "DP_R_CHELSA_processed", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # Ensure that the output path exists
  fs::dir_create(
    c(
      Path_CHELSA_In, Path_CHELSA_Out,
      IASDT.R::Path(Path_CHELSA_Out, c("Tif", "NC", "Processed"))))

  # # ..................................................................... ###

  # Prepare CHELSA metadata / download CHELSA data -----

  IASDT.R::CatTime("Prepare CHELSA metadata / download CHELSA data")
  TimePrepare <- lubridate::now(tzone = "CET")

  # 19 Bioclimatic variables (+ OtherVars, if not empty string) * 46 CC options
  CHELSA_Data <- IASDT.R::CHELSA_Prepare(
    EnvFile = EnvFile, Download = Download, NCores = Download_NCores,
    Overwrite = Overwrite, OtherVars = OtherVars,
    Download_Attempts = Download_Attempts, Sleep = Sleep)

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
          .f = IASDT.R::CheckTiff, warning = FALSE,
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
        file = IASDT.R::Path(Path_CHELSA_Out, "ProblematicTiffs.txt"))

      stop(
        "Not all input tiff files are available and valid. ",
        "Check `ProblematicTiffs.txt`", call. = FALSE
      )
    }

    # CHELSA files that will not be processed
    Diff <- setdiff(
      list.files(Path_CHELSA_In, pattern = "CHELSA.+.tif$", full.names = TRUE),
      CHELSA_Data$Path_Down)

    if (length(Diff) > 0) {
      message(
        " >> Only Bioclimatic variables and variables identified in ",
        "`OtherVars`, if any, will be processed (",
        nrow(CHELSA_Data), " files)\n >> ", length(Diff),
        " files will not be processed.\n",
        " >> See `NotProcessed.txt` for the list of files")

      readr::write_lines(
        x = Diff, file = IASDT.R::Path(Path_CHELSA_Out, "NotProcessed.txt"))
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
  IASDT.R::CatTime(
    "Exclude previously processed files (after checking)", Level = 1)

  if (OverwriteProcessed) {
    CHELSA2Process <- dplyr::select(
      .data = CHELSA_Data, Path_Down, Path_Out_NC, Path_Out_tif)
  } else {
    CHELSA2Process <- CHELSA_Data %>%
      dplyr::select(Path_Down, Path_Out_NC, Path_Out_tif) %>%
      dplyr::mutate(
        Process = furrr::future_map2_lgl(
          .x = Path_Out_NC, .y = Path_Out_tif,
          .f = ~ {
            NC_Okay <- IASDT.R::CheckTiff(.x, warning = FALSE)
            Tif_Okay <- IASDT.R::CheckTiff(.y, warning = FALSE)
            return(isFALSE(NC_Okay && Tif_Okay))
          },
          .options = furrr::furrr_options(seed = TRUE, packages = "IASDT.R"))
      ) %>%
      dplyr::filter(Process) %>%
      dplyr::select(-"Process")
  }

  # Processing CHELSA files
  IASDT.R::CatTime("Processing CHELSA files", Level = 1)

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
            repeat {
              Attempt <- Attempt + 1
              Try <- try(
                IASDT.R::CHELSA_Project(
                  Metadata = FileMetadata, EnvFile = EnvFile,
                  CompressLevel = CompressLevel),
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
              if (IASDT.R::CheckTiff(Path_Out_tif, warning = FALSE) &&
                  IASDT.R::CheckTiff(Path_Out_NC, warning = FALSE)) {
                return(FALSE)
              } else {
                return(TRUE)
              }
            }
          },
          .options = furrr::furrr_options(
            seed = TRUE,
            packages = c("dplyr", "terra", "IASDT.R", "tibble", "ncdf4"),
            globals = c("CHELSA_Data", "EnvFile", "CompressLevel"))
        )
      ) %>%
      dplyr::filter(Failed)

    if (nrow(CHELSA2Process) > 0) {
      readr::write_lines(
        x = CHELSA2Process$Path_Down,
        file = IASDT.R::Path(Path_CHELSA_Out, "FailedProcessing.txt"))

      stop(
        "\n >> ", nrow(CHELSA2Process), " files failed to process.\n",
        " >> Check `FailedProcessing.txt` for more details", call. = FALSE)
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
    paste(sep = ".*?_", collapse = "|")

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

      FilePath = IASDT.R::Path(
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
    file = IASDT.R::Path(Path_CHELSA_Out, "CHELSA_Processed.RData")
  )

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  CHELSA_Processed_DT <- dplyr::select(CHELSA_Processed, -Processed_Maps)

  save(
    CHELSA_Processed_DT,
    file = IASDT.R::Path(Path_CHELSA_Out, "CHELSA_Processed_DT.RData"))

  readr::write_csv(
    x = CHELSA_Processed_DT,
    file = IASDT.R::Path(Path_CHELSA_Out, "CHELSA_Processed_DT.csv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nProcessing CHELSA data took ")

  return(invisible(NULL))
}
