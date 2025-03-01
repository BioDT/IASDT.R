## |------------------------------------------------------------------------| #
# CHELSA_Prepare ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 2

CHELSA_Prepare <- function(
    EnvFile = ".env", Download = FALSE, NCores = 8L,
    Overwrite = FALSE, Download_Attempts = 10L, Sleep = 5L, OtherVars = "npp") {

  # # ..................................................................... ###

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("EnvFile", "OtherVars"), Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("Download", "Overwrite"),
    Type = "logical")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("Sleep", "Download_Attempts", "NCores"),
    Type = "numeric")

  rm(AllArgs, envir = environment())

  if (NCores < 1) {
    stop("`NCores` must be a positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Variable <- Path_Down <- TimePeriod <- Ext <- ClimateScenario <- Variable <-
    Path_CHELSA_In <- File <- Path_Out_tif <- Path_Out_NC <- Path_CHELSA_Out <-
    Path_DwnLinks <- URL <- Folder <- URL_File <- ClimateModel <- Exclude <-
    BaseURL <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_CHELSA_Out", "DP_R_CHELSA_processed", FALSE, FALSE,
    "Path_CHELSA_In", "DP_R_CHELSA_raw", FALSE, FALSE,
    "Path_DwnLinks", "DP_R_CHELSA_links", TRUE, FALSE,
    "BaseURL", "DP_R_CHELSA_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  fs::dir_create(c(Path_CHELSA_In, Path_CHELSA_Out))

  # String to be matched in the variable names
  SelectedVars <- c("^bio|^Bio", OtherVars) %>%
    # Only keep non-empty strings. If `OtherVars` = "", only bioclimatic
    # variables will be processed.
    stringr::str_subset(".+") %>%
    # Combine the strings into a single string separated by "|". This matches
    # any variable starting with "bio" or "Bio" or any of the characters in
    # `OtherVars`.
    paste(collapse = "|")

  IASDT.R::CatTime("Prepare CHELSA metadata", Level = 1)
  CHELSA_Metadata <- list.files(
    path = Path_DwnLinks, recursive = TRUE, full.names = TRUE,
    pattern = "DwnLinks_Climatologies_.+txt$") %>%
    dplyr::tibble(URL_File = .) %>%
    # Add download links
    dplyr::mutate(
      URL = purrr::map(.x = URL_File, .f = readr::read_lines),
      URL_File = basename(URL_File)) %>%
    tidyr::unnest_longer("URL") %>%
    dplyr::mutate(
      URL = purrr::map_chr(URL, stringr::str_trim),
      Folder = purrr::map_chr(URL, stringr::str_remove_all, pattern = BaseURL),
      File = purrr::map_chr(Folder, basename),
      Folder = purrr::map_chr(Folder, dirname),

      # Extract time period
      TimePeriod = purrr::map_chr(
        URL_File, stringr::str_remove_all,
        pattern = "DwnLinks_Climatologies_|.txt"),

      # File extension
      Ext = purrr::map_chr(URL, tools::file_ext),

      ModelScenario = purrr::map2(
        .x = Folder, .y = TimePeriod,
        .f = ~ {
          # assign "Current" for files represent current climates
          if (.y == "1981-2010") {
            Out <- tibble::tibble(
              ClimateModel = "Current", ClimateScenario = "Current")

          } else {
            ClimateModels <- c(
              "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
              "MRI-ESM2-0", "UKESM1-0-LL")
            ClimateScenarios <- c("ssp126", "ssp370", "ssp585")

            ClimateModel <- stringr::str_extract(.x, ClimateModels) %>%
              na.omit() %>%
              as.character()
            ClimateScenario <- stringr::str_extract(.x, ClimateScenarios) %>%
              na.omit() %>%
              as.character()
            Out <- tibble::tibble(
              ClimateModel = ClimateModel, ClimateScenario = ClimateScenario)
          }
          return(Out)
        }
      )
    ) %>%
    tidyr::unnest_wider("ModelScenario") %>%
    dplyr::filter(
      Folder != "climatologies/2011-2040/UKESM1-0-LL/ssp126") %>%
    dplyr::mutate(
      Variable = purrr::pmap_chr(
        .l = list(File, TimePeriod, ClimateScenario, Ext, ClimateModel),
        .f = function(File, TimePeriod, ClimateScenario, Ext, ClimateModel) {
          stringr::str_remove_all(
            string = File,
            pattern = paste0(
              "_r1i1p1f1_w5e5_|_norm|CHELSA_|V.2.1|_V\\.2\\.1|", TimePeriod,
              "|.", Ext, "|", ClimateScenario)) %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue(
                "{ClimateModel}|{tolower(ClimateModel)}")) %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue(
                '{TimePeriod}|{stringr::str_replace(TimePeriod, "-", "_")}')
            ) %>%
            stringr::str_remove_all(pattern = "__|___") %>%
            stringr::str_remove_all(pattern = "^_|_$")
        }
      )
    ) %>%
    # Process only selected variables
    dplyr::filter(stringr::str_detect(Variable, SelectedVars)) %>%
    dplyr::mutate(

      Path_Down = purrr::map_chr(
        .x = File, .f = ~ IASDT.R::Path(Path_CHELSA_In, .x)),

      Path_Out_tif = purrr::map_chr(
        .x = File, .f = ~ IASDT.R::Path(Path_CHELSA_Out, "Tif", .x)),

      Path_Out_NC = purrr::map_chr(
        .x = Path_Out_tif, .f = stringr::str_replace_all,
        pattern = "Tif", replacement = "NC"),
      Path_Out_NC = purrr::map_chr(
        .x = Path_Out_NC, .f = stringr::str_replace_all,
        pattern = ".tif$", replacement = ".nc"),

      DownCommand = purrr::map2_chr(
        .x = URL, .y = Path_Down,
        .f = ~ stringr::str_glue('curl -k -L "{.x}" -o "{.y}" --silent')),

      # Unique name for variable / time combination
      OutName = purrr::pmap_chr(
        .l = list(Variable, TimePeriod, ClimateModel, ClimateScenario),
        .f = function(Variable, TimePeriod, ClimateModel, ClimateScenario) {
          paste0(
            Variable, "_", TimePeriod, "_",
            ClimateModel, "_", ClimateScenario) %>%
            stringr::str_replace(
              pattern = "1981-2010_Current_Current",
              replacement = "1981-2010_Current")
        }
      )) %>%
    dplyr::select(-"Folder") %>%
    dplyr::left_join(IASDT.R::CHELSA_Vars, by = "Variable")

  # # ..................................................................... ###

  if (Download) {

    IASDT.R::CatTime("Downloading CHELSA files", Level = 1)

    if (NCores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    # # |||||||||||||||||||||||||||||||| ###

    if (Overwrite) {

      Data2Down <- CHELSA_Metadata

    } else {

      IASDT.R::CatTime("Exclude available and valid Tiff files", Level = 1)

      Data2Down <- CHELSA_Metadata %>%
        dplyr::mutate(
          Exclude = furrr::future_map_lgl(
            .x = Path_Down,
            .f = ~ {
              if (file.exists(.x)) {
                if (isFALSE(IASDT.R::CheckTiff(.x))) {
                  fs::file_delete(.x)
                  Out <- TRUE
                } else {
                  Out <- FALSE
                }
              } else {
                Out <- TRUE
              }

              invisible(gc())
              return(Out)
            },
            .options = furrr::furrr_options(
              seed = TRUE, packages = c("IASDT.R", "fs")),
            .progress = FALSE
          )
        ) %>%
        dplyr::filter(Exclude)
    }

    # # |||||||||||||||||||||||||||||||| ###

    # Download missing files
    IASDT.R::CatTime("Download missing CHELSA files", Level = 1)

    if (nrow(Data2Down) > 0) {
      furrr::future_walk(
        .x = seq_len(nrow(Data2Down)),
        .f = ~ {
          PathOut <- Data2Down$Path_Down[.x]
          URL <- Data2Down$URL[.x]

          Try <- 0
          repeat {
            Try <- Try + 1
            Down <- paste0(
              "curl --connect-timeout 60 --max-time 1200 -o ",
              PathOut, " ", URL, " --silent") %>%
              system()

            if (file.exists(PathOut)) {
              if (IASDT.R::CheckTiff(PathOut)) {
                break
              } else {
                fs::file_delete(PathOut)
              }
            }
            if (Try > Download_Attempts) {
              break
            }
            Sys.sleep(Sleep)
          }
          Sys.sleep(Sleep)
        },
        .options = furrr::furrr_options(
          seed = TRUE, globals = c("Data2Down", "Sleep"),
          packages = c("IASDT.R", "fs")))
    }

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }

  }

  # # ..................................................................... ###

  # Save to disk -----

  IASDT.R::CatTime("Save metadata to disk", Level = 1)

  save(
    CHELSA_Metadata,
    file = IASDT.R::Path(Path_CHELSA_Out, "CHELSA_Metadata.RData"))

  readr::write_csv(
    x = CHELSA_Metadata,
    file = IASDT.R::Path(Path_CHELSA_Out, "CHELSA_Metadata.csv"))

  # # ..................................................................... ###

  return(CHELSA_Metadata)
}
