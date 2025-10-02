## |------------------------------------------------------------------------| #
# Merge model chains ----
## |------------------------------------------------------------------------| #

#' Merge model chains into `Hmsc` and `coda` objects
#'
#' These functions merge posterior chains from multiple runs of `Hmsc` models
#' into unified `Hmsc` and `coda` objects, facilitating further analysis. They
#' check for missing or incomplete chains, optionally report these issues, and
#' save the processed results to disk. `mod_merge_chains` handles regular
#' models, while `mod_merge_chains_cv` is designed for cross-validated models.
#'
#' @param model_dir Character. Path to the root directory where the model was
#'   fitted. For `mod_merge_chains`, subdirectories `model_fitted` and
#'   `model_coda` are created within this path to store the merged `Hmsc` and
#'   `coda` objects, respectively. For `mod_merge_chains_cv`, merged objects are
#'   stored under `model_fitting_cv/model_fitted`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to 8L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param model_info_name Character. Name of the file (without extension) where
#'   updated model information is saved. If `NULL`, overwrites the existing
#'   `model_info.RData` file in `model_dir` directory. If specified, creates a
#'   new `.RData` file with this name in `model_dir` directory. Applies only to
#'   `mod_merge_chains`.
#' @param print_incomplete Logical. If `TRUE`, prints the names of model
#'   variants that failed to merge due to missing or incomplete chains. Defaults
#'   to `TRUE`.
#' @param from_json Logical. Whether to convert loaded models from JSON format
#'   before reading. Defaults to `FALSE`.
#' @param out_extension Character. File extension (without dot) for output files
#'   containing merged `Hmsc` and `coda` objects. Options are `qs2` (faster
#'   read/write via the `qs2` package) or `RData` (standard R format). Defaults
#'   to `qs2`.
#' @param cv_names Character vector. Names of cross-validation strategies to
#'   merge, matching those used during model setup. Defaults to `c("cv_dist",
#'   "cv_large")`. The names should be one of `cv_dist`, `cv_large`, or
#'   `cv_sac`. Applies only to `mod_merge_chains_cv`.
#' @export
#' @rdname mod_merge_chains
#' @name mod_merge_chains
#' @order 1
#' @author Ahmed El-Gabbas
#' @return Both functions return `invisible(NULL)` and save processed model
#'   information and merged objects to disk in the specified locations.
#' @details
#' - `mod_merge_chains` merges posterior chains from multiple runs of an `Hmsc`
#' model fitted without cross-validation. It checks for missing or incomplete
#' chains, aligns posteriors (using `alignPost = TRUE`, falling back to `FALSE`
#' if alignment fails), and saves a merged `Hmsc` object and a `coda` object for
#' MCMC diagnostics. It also records fitting time and memory usage from progress
#' files.
#' - `mod_merge_chains_cv` performs a similar merging process for
#' cross-validated `Hmsc` models, processing each fold of the specified
#' `cv_names` separately. It saves merged `Hmsc` objects per fold but does not
#' generate `coda` objects.

## |------------------------------------------------------------------------| #
# mod_merge_chains ----
## |------------------------------------------------------------------------| #

mod_merge_chains <- function(
    model_dir = NULL, n_cores = 8L, strategy = "multisession",
    model_info_name = NULL, print_incomplete = TRUE, from_json = FALSE,
    out_extension = "qs2") {

  .start_time <- lubridate::now(tzone = "CET")

  # Checking arguments ----

  ecokit::check_args(
    args_to_check = c("model_dir", "out_extension"), args_type = "character")
  ecokit::check_args(
    args_to_check = c("print_incomplete", "from_json"), args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_post <- m_name_fit <- path_fitted_model <- path_coda <-
    n_missing_chains <- missing_models <- model_finished <- path_mod_progress <-
    post_aligned_2 <- fitting_time <- fitting_memory <- NULL

  # # ..................................................................... ###

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!out_extension %in% c("qs2", "RData")) {
    ecokit::stop_ctx(
      "`out_extension` must be either 'qs2' or 'RData'.",
      out_extension = out_extension, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("Hmsc", "coda", "purrr", "ecokit", "dplyr", "fs"),
    strategy = strategy)

  # # ..................................................................... ###

  # Creating paths -----

  path_fitted_models <- fs::path(model_dir, "model_fitted")
  path_coda <- fs::path(model_dir, "model_coda")
  fs::dir_create(c(path_fitted_models, path_coda))

  # # ..................................................................... ###

  # Loading model info ----

  path_mod_info <- fs::path(model_dir, "model_info.RData")

  if (!file.exists(path_mod_info)) {
    ecokit::stop_ctx(
      "ModInfo file does not exist", path_mod_info = path_mod_info,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files ----

  path_model_fit <- fs::path(model_dir, "model_fitting_hpc")
  temp_files <- list.files(
    path = path_model_fit, pattern = ".rds_temp$", full.names = TRUE)

  if (length(temp_files) > 0) {
    ecokit::cat_time(
      paste0(
        "There are ", length(temp_files),
        " unsuccessful model variants to be removed"))

    temp_files_rds <- stringr::str_replace_all(temp_files, ".rds_temp$", ".rds")
    temp_files_progress <- stringr::str_replace_all(
      temp_files, "_post.rds_temp", "_progress.txt")

    purrr::walk(
      .x = c(temp_files_rds, temp_files, temp_files_progress),
      .f = ~{
        if (file.exists(.x)) {
          file.remove(.x)
        }
      })
  }

  # # ..................................................................... ###

  model_info_2 <- ecokit::load_as(path_mod_info)

  # Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(model_info_2)),
      future_max_size = 800L, strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time("Check if any posterior files is missing")
  # Check if any posterior files is missing
  model_info_2 <- model_info_2 %>%
    dplyr::mutate(
      post_missing = furrr::future_map_lgl(
        .x = path_post,
        .f = function(x) {

          purrr::map_lgl(
            .x = as.character(x),
            .f = function(y) {

              if (isFALSE(fs::file_exists(y))) {
                TRUE
              } else if (ecokit::check_data(y, warning = FALSE)) {
                FALSE
              } else {
                fs::file_delete(y)
                TRUE
              }

            }) %>%
            any()
        }),
      # delete these columns if already exist from previous function execution
      path_fitted_model = NULL, path_coda = NULL)

  # # ..................................................................... ###

  # Merge posteriors and save as Hmsc model / coda object
  ecokit::cat_time("Merge posteriors and save as Hmsc model / coda object")

  model_info_3 <- future.apply::future_lapply(
    X = seq_len(nrow(model_info_2)),
    FUN = function(x) {

      if (model_info_2$post_missing[[x]]) {
        return(
          list(
            path_fitted_model = NA_character_,
            path_coda = NA_character_, post_aligned_2 = NA))
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Merge fitted models
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      m_name_fit <- model_info_2$m_name_fit[[x]]

      path_fitted_model <- fs::path(
        path_fitted_models, paste0(m_name_fit, "_model.", out_extension))

      # Check if model exists and is valid
      model_file_okay <- ecokit::check_data(path_fitted_model, warning = FALSE)

      if (isFALSE(model_file_okay)) {

        # delete corrupted file
        if (file.exists(path_fitted_model)) {
          fs::file_delete(path_fitted_model)
        }

        # Get posteriors
        posts <- purrr::map(
          .x = as.character(model_info_2$path_post[[x]]),
          .f = IASDT.R::mod_get_posteriors, from_json = from_json)

        # Convert to Hmsc object
        # Try with `alignPost = TRUE`
        model_fit <- Hmsc::importPosteriorFromHPC(
          m = ecokit::load_as(model_info_2$path_m_init[x]),
          postList = posts, nSamples = model_info_2$m_samples[x],
          thin = model_info_2$m_thin[x],
          transient = model_info_2$m_transient[x],
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(model_fit, "try-error")) {
          model_fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = ecokit::load_as(model_info_2$path_m_init[x]),
              postList = posts, nSamples = model_info_2$m_samples[x],
              thin = model_info_2$m_thin[x],
              transient = model_info_2$m_transient[x], alignPost = FALSE),
            silent = TRUE)
          post_aligned_2 <- FALSE
        } else {
          post_aligned_2 <- TRUE
        }

        rm(posts, envir = environment())
        invisible(gc())

        if (inherits(model_fit, "try-error")) {
          return(
            list(
              path_fitted_model = NA_character_,
              path_coda = NA_character_, post_aligned_2 = NA))
        }
        ecokit::save_as(
          object = model_fit, object_name = paste0(m_name_fit, "_model"),
          out_path = path_fitted_model)

      } else {
        post_aligned_2 <- model_info_2$post_aligned[[x]]
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Convert to Coda object
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      path_coda <- fs::path(
        path_coda, paste0(m_name_fit, "_coda.", out_extension))

      if (!exists("model_fit") && !file.exists(path_fitted_model)) {
        return(
          list(
            path_fitted_model = NA_character_, path_coda = NA_character_,
            post_aligned_2 = NA))
      }

      if (!exists("model_fit")) {
        model_fit <- ecokit::load_as(path_fitted_model)
      }

      # Check if coda file exists and is valid
      if (isFALSE(ecokit::check_data(path_coda, warning = FALSE))) {

        if (file.exists(path_coda)) {
          fs::file_delete(path_coda)
        }

        mod_coda <- Hmsc::convertToCodaObject(
          model_fit, spNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE))

        ecokit::save_as(
          object = mod_coda, object_name = paste0(m_name_fit, "_coda"),
          out_path = path_coda)

        rm(mod_coda, envir = environment())
      }

      rm(model_fit, envir = environment())
      invisible(gc())

      # Return list of objects
      list(
        path_fitted_model = path_fitted_model, path_coda = path_coda,
        post_aligned_2 = post_aligned_2)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c(
      "out_extension", "model_info_2", "path_fitted_models",
      "from_json", "path_coda"))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  ecokit::cat_time("Extract information on elapsed time and memory usage")

  model_info_2 <- dplyr::mutate(model_info_2, model_posts = model_info_3) %>%
    tidyr::unnest_wider("model_posts") %>%
    dplyr::mutate(post_aligned = dplyr::coalesce(post_aligned_2)) %>%
    dplyr::select(-post_aligned_2) %>%
    dplyr::mutate(

      # Check if both merged fitted model and coda file exist
      model_finished = purrr::map2_lgl(
        .x = path_fitted_model, .y = path_coda,
        .f = ~all(file.exists(c(.x, .y)))),

      # Extract fitting time from the progress file
      fitting_time = purrr::map(
        .x = path_mod_progress,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(file) {
              if (file.exists(file)) {
                readr::read_lines(file = file, progress = FALSE) %>%
                  stringr::str_subset("Whole Gibbs sampler elapsed") %>%
                  stringr::str_remove("Whole Gibbs sampler elapsed") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(60) %>%
                  round(1)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean elapsed time
      fitting_time_mean = purrr::map2_dbl(
        .x = model_finished, .y = fitting_time,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }),

      # Maximum memory
      fitting_memory = purrr::map(
        .x = path_mod_progress,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(file) {
              if (file.exists(file)) {
                readr::read_lines(file = file, progress = FALSE) %>%
                  stringr::str_subset("Maximum resident set size") %>%
                  stringr::str_remove_all(
                    "\t|:|Maximum resident set size \\(kbytes\\)") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(1024 * 1024) %>%
                  round(2)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean memory
      fitting_memory_mean = purrr::map2_dbl(
        .x = model_finished, .y = fitting_memory,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Print to the console the name of failed models and number of missing chain
  # files

  if (print_incomplete) {
    missing_model_vars <- model_info_2 %>%
      dplyr::filter(!model_finished) %>%
      dplyr::mutate(
        n_missing_chains = purrr::map_int(
          .x = path_post, .f = ~sum(!file.exists(.x))),
        missing_models = paste0(m_name_fit, " (", n_missing_chains, " chains)")
      ) %>%
      dplyr::pull(missing_models) %>%
      gtools::mixedsort()

    if (length(missing_model_vars) > 0) {
      ecokit::cat_time("Unsuccessful models")
      purrr::walk(
        .x = missing_model_vars, .f = ecokit::cat_time,
        cat_timestamp = FALSE, level = 1L)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model_info to disk

  if (is.null(model_info_name)) {
    ecokit::save_as(
      object = model_info_2, object_name = "model_info",
      out_path = path_mod_info)
  } else {
    ecokit::save_as(
      object = model_info_2, object_name = model_info_name,
      out_path = fs::path(model_dir, paste0(model_info_name, ".RData")))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(init_time = .start_time, prefix = "Merging chains took ")

  return(invisible(NULL))
}

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_merge_chains_cv ----
# Merge chains for cross-validated models
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_merge_chains
#' @name mod_merge_chains
#' @order 2
#' @author Ahmed El-Gabbas

mod_merge_chains_cv <- function(
    model_dir = NULL, n_cores = 8L, strategy = "multisession",
    cv_names = c("cv_dist", "cv_large"), from_json = FALSE,
    out_extension = "qs2") {

  ecokit::cat_time("Merge chains for cross-validated models")
  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----

  ecokit::check_args(
    args_to_check = c("model_dir", "out_extension"), args_type = "character")
  ecokit::check_args(args_to_check = "from_json", args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  model_finished <- fitting_memory <- path_post <- path_mod_fitted <-
    fitting_time <- path_mod_progress <- cv_name <- NULL

  # # ..................................................................... ###

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!out_extension %in% c("qs2", "RData")) {
    ecokit::stop_ctx(
      "`out_extension` must be either 'qs2' or 'RData'.",
      out_extension = out_extension, include_backtrace = TRUE)
  }

  if (!all(cv_names %in% c("cv_dist", "cv_large", "cv_sac"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for cv_names argument. Valid values ",
        "are: 'cv_dist', 'cv_large', or `cv_sac`"),
      cv_names = cv_names, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "Hmsc", "purrr", "ecokit", "IASDT.R", "dplyr", "stringr", "fs"),
    strategy = strategy)

  # # ..................................................................... ###

  # Creating paths -----

  path_fitted_models <- fs::path(model_dir, "model_fitting_cv", "model_fitted")
  fs::dir_create(path_fitted_models)

  # # ..................................................................... ###

  # Loading CV model info -----

  path_cv_data <- fs::path(model_dir, "model_fitting_cv", "cv_data.RData")
  if (!file.exists(path_cv_data)) {
    ecokit::stop_ctx(
      "cv_data file does not exist", path_cv_data = path_cv_data,
      include_backtrace = TRUE)
  }
  if (isFALSE(ecokit::check_data(path_cv_data, warning = FALSE))) {
    ecokit::stop_ctx(
      "cv_data file is not a valid file", path_cv_data = path_cv_data,
      include_backtrace = TRUE)
  }

  cv_data <- ecokit::load_as(path_cv_data) %>%
    # filter only selected cross-validation strategies
    dplyr::filter(cv_name %in% stringr::str_remove(cv_names, "cv_"))

  # # ..................................................................... ###

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(cv_data)), future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  cv_data <- cv_data %>%
    dplyr::mutate(
      post_missing = furrr::future_map_lgl(
        .x = path_post,
        .f = function(x) {
          purrr::map_lgl(
            .x = as.character(x),
            .f = function(y) {

              if (isFALSE(fs::file_exists(y))) {
                TRUE
              } else if (ecokit::check_data(y, warning = FALSE)) {
                FALSE
              } else {
                fs::file_delete(y)
                TRUE
              }
            }) %>%
            any()
        },
        .options = furrr::furrr_options(seed = TRUE, packages = "fs")),

      path_mod_fitted = stringr::str_replace_all(
        path_mod_fitted, "RData$|qs2$", out_extension))

  invisible(gc())

  # # ..................................................................... ###

  # Merge posteriors and save as Hmsc object
  ecokit::cat_time("Merge posteriors and save as Hmsc object", level = 1L)

  cv_data2 <- future.apply::future_lapply(
    X = seq_len(nrow(cv_data)),
    FUN = function(x) {

      if (cv_data$post_missing[[x]]) {
        return(NA)
      }

      path_fitted <- cv_data$path_mod_fitted[[x]]

      # Check if model exists and is valid
      model_file_okay <- ecokit::check_data(path_fitted, warning = FALSE)

      if (isFALSE(model_file_okay)) {

        # delete corrupted file
        if (file.exists(path_fitted)) {
          fs::file_delete(path_fitted)
        }

        # Get posteriors
        posts <- purrr::map(
          .x = as.character(cv_data$path_post[[x]]),
          .f = IASDT.R::mod_get_posteriors, from_json = from_json)

        # Convert to Hmsc object
        model_init_rds <- ecokit::load_as(cv_data$path_mod_init_rds[x])
        model_init <- ecokit::load_as(cv_data$path_mod_init[x])
        model_fit <- Hmsc::importPosteriorFromHPC(
          m = model_init, postList = posts, nSamples = model_init_rds$samples,
          thin = model_init_rds$thin, transient = model_init_rds$transient,
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(model_fit, "try-error")) {
          model_fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = ecokit::load_as(cv_data$path_m_init[x]),
              postList = posts, nSamples = cv_data$m_samples[x],
              thin = cv_data$m_thin[x],
              transient = cv_data$m_transient[x], alignPost = FALSE),
            silent = TRUE)
          post_aligned <- FALSE
        } else {
          post_aligned <- TRUE
        }

        rm(posts, envir = environment())

        if (inherits(model_fit, "try-error")) {
          return(NA)
        }

        ecokit::save_as(
          object = model_fit,
          object_name = stringr::str_remove(
            basename(path_fitted), ".RData$|.qs2"),
          out_path = path_fitted)

      } else {
        post_aligned <- NA
      }

      invisible(gc())

      # Return list of objects
      post_aligned
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("out_extension", "cv_data", "from_json"))

  invisible(gc())

  # # ..................................................................... ###

  # Check saved Hmsc object and extract info on model fitting
  ecokit::cat_time(
    "Check saved Hmsc object and extract info on model fitting", level = 1L)

  # Check if both merged fitted model file exist
  cv_data <- cv_data %>%
    dplyr::mutate(

      post_aligned = unlist(cv_data2),

      model_finished = furrr::future_map_lgl(
        .x = path_mod_fitted, .f = ecokit::check_data, warning = FALSE,
        .options = furrr::furrr_options(seed = TRUE, packages = pkg_to_export)),

      # Extract fitting time from the progress file
      fitting_time = purrr::map(
        .x = path_mod_progress,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(file) {
              if (file.exists(file)) {
                readr::read_lines(file = file, progress = FALSE) %>%
                  stringr::str_subset("Whole Gibbs sampler elapsed") %>%
                  stringr::str_remove("Whole Gibbs sampler elapsed") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(60) %>%
                  round(1)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean elapsed time
      fitting_time_mean = purrr::map2_dbl(
        .x = model_finished, .y = fitting_time,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }),

      # Maximum memory
      fitting_memory = purrr::map(
        .x = path_mod_progress,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(file) {
              if (file.exists(file)) {
                readr::read_lines(file = file, progress = FALSE) %>%
                  stringr::str_subset("Maximum resident set size") %>%
                  stringr::str_remove_all(
                    "\t|:|Maximum resident set size \\(kbytes\\)") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(1024 * 1024) %>%
                  round(2)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean memory
      fitting_memory_mean = purrr::map2_dbl(
        .x = model_finished, .y = fitting_memory,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }))

  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # stopping the cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE)
    future::plan("sequential", gc = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model_info to disk

  ecokit::save_as(
    object = cv_data, object_name = "cv_data_fitted",
    out_path = fs::path(model_dir, "model_fitting_cv", "cv_data_fitted.RData"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(init_time = .start_time, prefix = "Merging chains took ")

  return(invisible(NULL))
}
