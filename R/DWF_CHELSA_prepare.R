## |------------------------------------------------------------------------| #
# CHELSA_prepare ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name CHELSA_data
#' @rdname CHELSA_data
#' @order 2

CHELSA_prepare <- function(
    env_file = ".env", download = FALSE, n_cores = 8L,
    strategy = "multisession", overwrite = FALSE,
    download_attempts = 10L, sleep = 5L, other_variables = "npp") {

  # # ..................................................................... ###

  ecokit::check_args(
    args_to_check = c("download", "overwrite"), args_type = "logical")
  ecokit::check_args(
    args_to_check = c("sleep", "download_attempts"), args_type = "numeric")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Variable <- Path_Down <- TimePeriod <- Ext <- ClimateScenario <- Variable <-
    Path_CHELSA_In <- File <- Path_Out_tif <- Path_Out_NC <- Path_CHELSA_Out <-
    Path_dwnload_links <- URL <- Folder <- URL_File <- ClimateModel <-
    Exclude <- chelsa_base_url <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_CHELSA_Out", "DP_R_chelsa_processed", FALSE, FALSE,
    "Path_CHELSA_In", "DP_R_chelsa_raw", FALSE, FALSE,
    "Path_dwnload_links", "DP_R_chelsa_links", TRUE, FALSE,
    "chelsa_base_url", "DP_R_chelsa_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("fs", "ecokit"), strategy = strategy)

  # # ..................................................................... ###

  fs::dir_create(c(Path_CHELSA_In, Path_CHELSA_Out))

  # String to be matched in the variable names
  SelectedVars <- c("^bio|^Bio", other_variables) %>%
    # Only keep non-empty strings. If `other_variables` = "", only bioclimatic
    # variables will be processed.
    stringr::str_subset(".+") %>%
    # Combine the strings into a single string separated by "|". This matches
    # any variable starting with "bio" or "Bio" or any of the characters in
    # `other_variables`.
    paste(collapse = "|")

  ecokit::cat_time("Prepare CHELSA metadata", level = 1L)

  if (!endsWith(chelsa_base_url, "/")) {
    chelsa_base_url <- paste0(chelsa_base_url, "/")
  }

  CHELSA_Metadata <- list.files(
    path = Path_dwnload_links, recursive = TRUE, full.names = TRUE,
    pattern = "dwnload_links_Climatologies_.+txt$") %>%
    dplyr::tibble(URL_File = .) %>%
    # Add download links
    dplyr::mutate(
      URL = purrr::map(.x = URL_File, .f = readr::read_lines, progress = FALSE),
      URL_File = basename(URL_File)) %>%
    tidyr::unnest_longer("URL") %>%
    dplyr::mutate(
      URL = purrr::map_chr(URL, stringr::str_trim),
      Folder = purrr::map_chr(
        URL, stringr::str_remove_all, pattern = chelsa_base_url),
      File = purrr::map_chr(Folder, basename),
      Folder = purrr::map_chr(Folder, dirname),

      # Extract time period
      TimePeriod = purrr::map_chr(
        URL_File, stringr::str_remove_all,
        pattern = "dwnload_links_Climatologies_|.txt"),

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
        })
    ) %>%
    tidyr::unnest_wider("ModelScenario") %>%
    # Files for 2011-2040/UKESM1-0-LL/ssp126 are duplicated
    dplyr::filter(
      Folder != fs::path(
        "climatologies", "2011-2040", "UKESM1-0-LL", "ssp126")) %>%
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
        })
    ) %>%
    # Process only selected variables
    dplyr::filter(stringr::str_detect(Variable, SelectedVars)) %>%
    dplyr::mutate(

      Path_Down = purrr::map_chr(
        .x = File, .f = ~ fs::path(Path_CHELSA_In, .x)),

      Path_Out_tif = purrr::map_chr(
        .x = File, .f = ~ fs::path(Path_CHELSA_Out, "Tif", .x)),

      Path_Out_NC = purrr::map_chr(
        .x = Path_Out_tif, .f = stringr::str_replace_all,
        pattern = "Tif", replacement = "NC"),

      Path_Out_NC = purrr::map_chr(
        .x = Path_Out_NC, .f = stringr::str_replace_all,
        pattern = ".tif$", replacement = ".nc"),
      DownCommand = purrr::map2_chr(
        .x = URL, .y = Path_Down,
        .f = ~ stringr::str_glue(
          'curl -k -L --connect-timeout 120 --max-time 600 --retry 5 \\
          "{.x}" -o "{.y}" --silent')),

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
        })
    ) %>%
    dplyr::select(-"Folder") %>%
    dplyr::left_join(IASDT.R::CHELSA_variables, by = "Variable")

  # # ..................................................................... ###

  if (download) {

    ecokit::cat_time("Downloading CHELSA files", level = 1L)

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # # |||||||||||||||||||||||||||||||| ###

    if (overwrite) {

      Data2Down <- CHELSA_Metadata

    } else {

      ecokit::cat_time("Exclude available and valid Tiff files", level = 1L)

      Data2Down <- CHELSA_Metadata %>%
        dplyr::mutate(
          Exclude = furrr::future_map_lgl(
            .x = Path_Down,
            .f = ~ {

              if (file.exists(.x)) {
                if (isFALSE(ecokit::check_tiff(.x, warning = FALSE))) {
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
              seed = TRUE, packages = pkg_to_export),
            .progress = FALSE)) %>%
        dplyr::filter(Exclude)
    }

    # # |||||||||||||||||||||||||||||||| ###

    # download missing files
    ecokit::cat_time("Download missing CHELSA files", level = 1L)

    if (nrow(Data2Down) > 0) {

      ecokit::cat_time(
        paste0(
          nrow(Data2Down), " of ", nrow(CHELSA_Metadata),
          " files need to be downloaded."),
        level = 2L)

      download_problems <- future.apply::future_lapply(
        X = seq_len(nrow(Data2Down)),
        FUN = function(x) {

          PathOut <- Data2Down$Path_Down[x]
          Try <- 0

          repeat {
            Try <- Try + 1
            Down <- system(Data2Down$DownCommand[x])
            Sys.sleep(5)

            if (ecokit::check_tiff(PathOut, warning = FALSE) ||
                Try > download_attempts) {
              download_problem <- FALSE
              break
            } else {
              try(fs::file_delete(PathOut), silent = TRUE)
              download_problem <- TRUE
            }
          }

          rm(Down, envir = environment())
          Sys.sleep(sleep)

          return(download_problem)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.packages = pkg_to_export,
        future.globals = c("Data2Down", "sleep", "download_attempts"))

      download_problems <- sum(unlist(download_problems))

      if (download_problems == 0) {
        ecokit::cat_time("All Tiff files were downloaded.", level = 2L)
      } else {
        ecokit::cat_time(
          paste0(download_problems, " Tiff files were not downloaded"),
          level = 2L)
      }

      rm(download_problems, envir = environment())

    } else {

      ecokit::cat_time(
        "All Tiff files were already available and valid", level = 2L)

    }

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("sequential", gc = TRUE)
    }

  } else {

    ecokit::cat_time("CHELSA files will not be downloaded", level = 1L)

  }

  # # ..................................................................... ###

  # Save to disk -----

  ecokit::cat_time("Save metadata to disk", level = 1L)

  save(
    CHELSA_Metadata,
    file = fs::path(Path_CHELSA_Out, "CHELSA_Metadata.RData"))

  readr::write_csv(
    x = CHELSA_Metadata,
    file = fs::path(Path_CHELSA_Out, "CHELSA_Metadata.csv"))

  # # ..................................................................... ###

  return(CHELSA_Metadata)
}
