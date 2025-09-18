## |------------------------------------------------------------------------| #
# chelsa_prepare ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name chelsa_data
#' @rdname chelsa_data
#' @order 2

chelsa_prepare <- function(
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
  variable <- path_down <- time_period <- extension <- climate_scenario <-
    path_chelsa_in <- path_file <- path_out_tif <- path_out_nc <-
    path_chelsa_out <- path_dwnload_links <- url <- path_dir <- url_file <-
    climate_model <- exclude <- chelsa_base_url <- NULL

  # # ..................................................................... ###

  # Environment variables -----
  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_chelsa_out", "DP_R_chelsa_processed", FALSE, FALSE,
    "path_chelsa_in", "DP_R_chelsa_raw", FALSE, FALSE,
    "path_dwnload_links", "DP_R_chelsa_links", TRUE, FALSE,
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

  fs::dir_create(c(path_chelsa_in, path_chelsa_out))

  # String to be matched in the variable names
  selected_vars <- c("^bio|^Bio", other_variables) %>%
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

  chelsa_metadata <- list.files(
    path = path_dwnload_links, recursive = TRUE, full.names = TRUE,
    pattern = "dwnload_links_climatologies_.+txt$") %>%
    dplyr::tibble(url_file = .) %>%
    # Add download links
    dplyr::mutate(
      url = purrr::map(.x = url_file, .f = readr::read_lines, progress = FALSE),
      url_file = basename(url_file)) %>%
    tidyr::unnest_longer("url") %>%
    dplyr::mutate(
      url = purrr::map_chr(url, stringr::str_trim),
      path_dir = purrr::map_chr(
        url, stringr::str_remove_all, pattern = chelsa_base_url),
      path_file = purrr::map_chr(path_dir, basename),
      path_dir = purrr::map_chr(path_dir, dirname),

      # Extract time period
      time_period = purrr::map_chr(
        url_file, stringr::str_remove_all,
        pattern = "dwnload_links_climatologies_|.txt"),

      # File extension
      extension = purrr::map_chr(url, tools::file_ext),

      model_scenario = purrr::map2(
        .x = path_dir, .y = time_period,
        .f = ~ {
          # assign "current" for files represent current climates
          if (.y == "1981-2010") {
            out <- tibble::tibble(
              climate_model = "current", climate_scenario = "current")

          } else {
            climate_models <- c(
              "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
              "MRI-ESM2-0", "UKESM1-0-LL")
            climate_scenarios <- c("ssp126", "ssp370", "ssp585")

            climate_model <- stringr::str_extract(.x, climate_models) %>%
              na.omit() %>%
              as.character()
            climate_scenario <- stringr::str_extract(.x, climate_scenarios) %>%
              na.omit() %>%
              as.character()
            out <- tibble::tibble(
              climate_model = climate_model,
              climate_scenario = climate_scenario)
          }
          return(out)
        })) %>%
    tidyr::unnest_wider("model_scenario") %>%
    # Files for 2011-2040/UKESM1-0-LL/ssp126 are duplicated
    dplyr::filter(
      path_dir != fs::path(
        "climatologies", "2011-2040", "UKESM1-0-LL", "ssp126")) %>%
    dplyr::mutate(
      variable = purrr::pmap_chr(
        .l = list(
          path_file, time_period, climate_scenario, extension, climate_model),
        .f = function(
    path_file, time_period, climate_scenario, extension, climate_model) {
          stringr::str_remove_all(
            string = path_file,
            pattern = paste0(
              "_r1i1p1f1_w5e5_|_norm|CHELSA_|V.2.1|_V\\.2\\.1|", time_period,
              "|.", extension, "|", climate_scenario)) %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue(
                "{climate_model}|{tolower(climate_model)}")) %>%
            stringr::str_remove_all(
              pattern = stringr::str_glue(
                '{time_period}|{stringr::str_replace(time_period, "-", "_")}')
            ) %>%
            stringr::str_remove_all(pattern = "__|___") %>%
            stringr::str_remove_all(pattern = "^_|_$")
        })) %>%
    # Process only selected variables
    dplyr::filter(stringr::str_detect(variable, selected_vars)) %>%
    dplyr::mutate(

      path_down = purrr::map_chr(
        .x = path_file, .f = ~ fs::path(path_chelsa_in, .x)),

      path_out_tif = purrr::map_chr(
        .x = path_file, .f = ~ fs::path(path_chelsa_out, "chelsa_tif", .x)),

      path_out_nc = purrr::map_chr(
        .x = path_out_tif, .f = stringr::str_replace_all,
        pattern = "chelsa_tif", replacement = "chelsa_netcdf"),

      path_out_nc = purrr::map_chr(
        .x = path_out_nc, .f = stringr::str_replace_all,
        pattern = ".tif$", replacement = ".nc"),
      down_command = purrr::map2_chr(
        .x = url, .y = path_down,
        .f = ~ stringr::str_glue(
          'curl -k -L --connect-timeout 240 --max-time 1200 --retry 5 \\
          "{.x}" -o "{.y}" --silent')),

      # Unique name for variable / time combination
      out_name = purrr::pmap_chr(
        .l = list(variable, time_period, climate_model, climate_scenario),
        .f = function(variable, time_period, climate_model, climate_scenario) {
          paste0(
            variable, "_", time_period, "_",
            climate_model, "_", climate_scenario) %>%
            stringr::str_replace(
              pattern = "1981-2010_current_current",
              replacement = "1981-2010_current")
        })) %>%
    dplyr::select(-"path_dir") %>%
    dplyr::left_join(IASDT.R::chelsa_variables, by = "variable")

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

      data_to_down <- chelsa_metadata

    } else {

      ecokit::cat_time("Exclude available and valid Tiff files", level = 1L)

      data_to_down <- chelsa_metadata %>%
        dplyr::mutate(
          exclude = furrr::future_map_lgl(
            .x = path_down,
            .f = ~ {

              if (file.exists(.x)) {
                if (isFALSE(ecokit::check_tiff(.x, warning = FALSE))) {
                  fs::file_delete(.x)
                  out <- TRUE
                } else {
                  out <- FALSE
                }
              } else {
                out <- TRUE
              }

              invisible(gc())
              return(out)
            },
            .options = furrr::furrr_options(
              seed = TRUE, packages = pkg_to_export),
            .progress = FALSE)) %>%
        dplyr::filter(exclude)
    }

    # # |||||||||||||||||||||||||||||||| ###

    # download missing files
    ecokit::cat_time("Download missing CHELSA files", level = 1L)

    if (nrow(data_to_down) > 0) {

      ecokit::cat_time(
        paste0(
          nrow(data_to_down), " of ", nrow(chelsa_metadata),
          " files need to be downloaded."),
        level = 2L)

      download_problems <- future.apply::future_lapply(
        X = seq_len(nrow(data_to_down)),
        FUN = function(x) {

          path_out <- data_to_down$path_down[x]
          try_n <- 0

          repeat {
            try_n <- try_n + 1
            download <- system(data_to_down$down_command[x])
            Sys.sleep(5)

            if (ecokit::check_tiff(path_out, warning = FALSE) ||
                try_n > download_attempts) {
              download_problem <- FALSE
              break
            } else {
              try(fs::file_delete(path_out), silent = TRUE)
              download_problem <- TRUE
            }
          }

          rm(download, envir = environment())
          Sys.sleep(sleep)

          download_problem
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.packages = pkg_to_export,
        future.globals = c("data_to_down", "sleep", "download_attempts"))

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
    chelsa_metadata,
    file = fs::path(path_chelsa_out, "chelsa_metadata.RData"))

  readr::write_csv(
    x = chelsa_metadata,
    file = fs::path(path_chelsa_out, "chelsa_metadata.csv"))

  # # ..................................................................... ###

  return(chelsa_metadata)
}
