## |------------------------------------------------------------------------| #
# CHELSA_Process ----
## |------------------------------------------------------------------------| #

#' Process CHELSA Climate Data
#'
#' This function processes CHELSA climate data, with an option to download them
#' first. It processes each bioclimatic variable to the European scale and
#' reference grid, and outputs in TIFF and NetCDF (NC) formats. It also saves
#' grouped data for each of the 46 climate scenarios.
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
#' @author Ahmed El-Gabbas
#' @export

CHELSA_Process <- function(
    FromHPC = TRUE, EnvFile = ".env", NCores = 8, Download = FALSE,
    Overwrite = FALSE, Download_Attempts = 10, Sleep = 5,
    BaseURL = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/",
    Download_NCores = 4, CompressLevel = 5, OverwriteProcessed = FALSE) {

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

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs)

  if (NCores < 1 || Download_NCores < 1) {
    stop("`NCores` must be a positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_CHELSA_In <- Path_CHELSA_Out <- TimePeriod <- ClimateModel <-
    ClimScenario <- Path_Out_tif <- Processed_Name <- Processed_File <-
    File_List <- Processed_Maps <- Problematic <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  IASDT.R::CatTime("Environment variables")

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CHELSA_In", "DP_R_CHELSA_Input", TRUE, FALSE,
      "Path_CHELSA_Out", "DP_R_CHELSA_Output", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CHELSA_In", "DP_R_CHELSA_Input_Local", TRUE, FALSE,
      "Path_CHELSA_Out", "DP_R_CHELSA_Output_Local", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # Ensure that the output path exists
  fs::dir_create(
    c(Path_CHELSA_Out,
      file.path(Path_CHELSA_Out, c("Tif", "NC", "Processed"))))

  # # ..................................................................... ###

  # Prepare / download CHELSA data -----

  IASDT.R::CatTime("Prepare / download CHELSA data")

  # 874 files (19 Bioclimatic variables * 46 CC options)
  CHELSA_Data <- IASDT.R::CHELSA_Prepare(
    EnvFile = EnvFile, FromHPC = FromHPC, Download = Download,
    NCores = Download_NCores, Overwrite = Overwrite, BaseURL = BaseURL,
    Download_Attempts = Download_Attempts, Sleep = Sleep)

  # # ..................................................................... ###

  # Check CHELSA files -----

  IASDT.R::CatTime("Check CHELSA files")


  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
    timeout = 1200)

  if (NCores > 1) {
    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)
  } else {
    future::plan("sequential", gc = TRUE)
  }

  AllOkay <- future.apply::future_lapply(
    X = CHELSA_Data$Path_Down,
    FUN = function(x) {
      file.exists(x) && IASDT.R::CheckTiff(x)
    },
    future.seed = TRUE, future.scheduling = Inf,
    future.packages = c("IASDT.R", "terra", "stringr")) %>%
    unlist() %>%
    all()

  future::plan("sequential", gc = TRUE)

  if (isFALSE(AllOkay)) {
    stop("Not all tiff files are available on disk")
  }

  # CHELSA files that will not be processed
  Diff <- setdiff(
    list.files(Path_CHELSA_In, pattern = "CHELSA.+.tif$", full.names = TRUE),
    CHELSA_Data$Path_Down)
  if (length(Diff) > 0) {
    message(
      paste0(
        " >> Only Bioclimatic variables will be processed (",
        nrow(CHELSA_Data), " files)\n >> ", length(Diff),
        " files will not be processed.\n",
        " >> See `NotProcessed.txt` for the list of files"))
    readr::write_lines(
      x = Diff, file = file.path(Path_CHELSA_Out, "NotProcessed.txt"))
  }

  rm(AllOkay, Diff, EnvVars2Read)

  # # ..................................................................... ###

  # Processing CHELSA data -----

  IASDT.R::CatTime("Processing CHELSA maps")
  TimeStart <- lubridate::now(tzone = "CET")

  if (NCores > 1) {
    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)
  } else {
    future::plan("sequential", gc = TRUE)
  }


  ProcessingIssues <- future.apply::future_lapply(
    X = seq_len(nrow(CHELSA_Data)),
    FUN = function(x) {

      FileMetadata <- dplyr::slice(CHELSA_Data, x)
      FileInput <- FileMetadata$Path_Down
      Process <- FALSE

      TifOkay <- file.exists(FileMetadata$Path_Out_tif)
      if (TifOkay) {
        TifOkay <- IASDT.R::CheckTiff(FileMetadata$Path_Out_tif)
      }

      NC_Okay <- file.exists(FileMetadata$Path_Out_NC)
      if (NC_Okay) {
        NC_Okay <- IASDT.R::CheckTiff(FileMetadata$Path_Out_NC)
      }

      Processed <- NC_Okay && NC_Okay

      if (isFALSE(Processed) || OverwriteProcessed) {
        Process <- TRUE
      }

      if (Process) {
        if (!file.exists(FileInput)) {
          stop(
            paste0(
              "Input file is either does not exist or corrupted",
              FileInput),
            call. = FALSE)
        }

        Attempt <- 0
        while (TRUE) {
          Attempt <- Attempt + 1
          Try <- try(IASDT.R::CHELSA_Project(
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
          return(tibble::tibble(ID = x, Problematic = TRUE))
        } else {
          Tif_NC_Okay <- all(
            file.exists(FileMetadata$Path_Out_tif),
            IASDT.R::CheckTiff(FileMetadata$Path_Out_tif),
            file.exists(FileMetadata$Path_Out_NC),
            IASDT.R::CheckTiff(FileMetadata$Path_Out_NC))

          if (Tif_NC_Okay) {
            return(tibble::tibble(ID = x, Problematic = FALSE))
          } else {
            return(tibble::tibble(ID = x, Problematic = TRUE))
          }
        }
      } else {
        return(tibble::tibble(ID = x, Problematic = FALSE))
      }
    },
    future.seed = TRUE, future.scheduling = Inf,
    future.globals = c(
      "CHELSA_Data", "EnvFile", "FromHPC",
      "CompressLevel", "OverwriteProcessed"),
    future.packages = c("dplyr", "terra", "IASDT.R")) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(Problematic)


  if (nrow(ProcessingIssues) > 0) {
    FailedProcessing <- CHELSA_Data$Path_Down[ProcessingIssues$ID]
    readr::write_lines(
      x = FailedProcessing,
      file = file.path(Path_CHELSA_Out, "FailedProcessing.txt"))
    stop(
      paste(
        "\n >> ", length(FailedProcessing), " files failed to process.\n",
        " >> Check `FailedProcessing.txt` for more details"),
      call. = FALSE)
  }

  rm(ProcessingIssues)
  future::plan("sequential", gc = TRUE)

  IASDT.R::CatDiff(
    InitTime = TimeStart,
    Prefix = "Processing CHELSA data was finished in ", Level = 1)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Group CHELSA data by scenario

  # Export each time*CC options as separate `*.RData` file


  IASDT.R::CatTime("Group CHELSA data by scenario")

  if (NCores > 1) {
    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)
  } else {
    future::plan("sequential", gc = TRUE)
  }

  CHELSA_Summary <- CHELSA_Data %>%
    dplyr::select(TimePeriod, ClimateModel, ClimScenario, Path_Out_tif) %>%
    dplyr::slice(gtools::mixedorder(Path_Out_tif)) %>%
    dplyr::summarise(
      File_List = list(Path_Out_tif),
      .by = c(TimePeriod, ClimateModel, ClimScenario)) %>%
    dplyr::mutate(
      Processed_Name = paste0(
        "R_", TimePeriod, "_", ClimateModel, "_", ClimScenario),
      Processed_Name = stringr::str_replace(
        Processed_Name, "1981-2010_Current_Current", "Current"),
      Processed_Name = stringr::str_replace_all(Processed_Name, "-", "_"),

      Processed_File = file.path(
        Path_CHELSA_Out, "Processed", paste0(Processed_Name, ".RData")),

      Processed_Maps = furrr::future_pmap(
        .l = list(File_List, Processed_File, Processed_Name),
        .f = function(File_List, Processed_File, Processed_Name) {

          MapNames <- basename(File_List) %>%
            stringr::str_extract("bio.._|bio._") %>%
            stringr::str_remove_all("_$")

          SortedNames <- gtools::mixedsort(MapNames)

          Map <- terra::rast(File_List) %>%
            stats::setNames(MapNames) %>%
            terra::subset(SortedNames) %>%
            IASDT.R::setRastCRS() %>%
            IASDT.R::setRastVals() %>%
            terra::wrap()

          # save to disk
          IASDT.R::SaveAs(
            InObj = Map, OutObj = Processed_Name, OutPath = Processed_File)

          return(Map)
        },
        .progress = FALSE,
        .options = furrr::furrr_options(
          seed = TRUE,
          packages = c("terra", "IASDT.R", "stringr"),
          globals = "CHELSA_Data")))

  future::plan("sequential", gc = TRUE)

  save(
    CHELSA_Summary,
    file = file.path(Path_CHELSA_Out, "CHELSA_Summary.RData"))

  CHELSA_Processed_DT <- dplyr::select(CHELSA_Summary, -Processed_Maps)

  save(
    CHELSA_Processed_DT,
    file = file.path(Path_CHELSA_Out, "CHELSA_Processed_DT.RData"))

  readr::write_csv(
    x = CHELSA_Processed_DT,
    file = file.path(Path_CHELSA_Out, "CHELSA_Processed_DT.csv"))


  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    ChunkText = "Processing CHELSA data took ")

  return(invisible(NULL))
}
