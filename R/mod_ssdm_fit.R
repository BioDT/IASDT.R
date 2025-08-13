
# # ========================================================================= #
# fit_sdm_models ------
# # ========================================================================= #

#' Species Distribution Modelling Workflow for Single-Species Models
#'
#' This comprehensive workflow implements single-species species distribution
#' models (sSDMs) for invasive alien plant species in Europe at the habitat
#' level. It orchestrates the entire process from data preparation to model
#' fitting, evaluation, and prediction across current and future climate
#' scenarios. The workflow employs the `sdm` R package for model fitting and
#' handles cross-validation, parallel processing, and various environmental
#' predictors.
#'
#' @param sdm_method Character. A single SDM algorithm to use for fitting
#'   models. Valid values: "glm", "glmpoly", "gam", "glmnet", "mars", "gbm",
#'   "rf", "ranger", "cart", "rpart", "maxent", "mlp", "rbf", "svm", "mda", and
#'   "fda". These correspond to selected methods supported by the `sdm` package.
#'   For details and supported options, see [sdm::getmethodNames()].
#' @param model_settings List or NULL. List of model-specific settings. If
#'   `NULL`, defaults to custom settings defined within the workflow.
#' @param model_dir Character. Path to the directory containing model data and
#'   where outputs and results will be saved. Model data are prepared using the
#'   [mod_prepare_HPC()] and [mod_prepare_data()] functions.
#' @param hab_abb Character. Abbreviation for a single SynHab habitat type.
#'   Valid values: "0", "1", "2", "3", "4a", "4b", "10", "12a", "12b". See
#'   [mod_prepare_HPC()] for more details.
#' @param cv_type Character. Cross-validation type. One of `CV_Dist` (default)
#'   or `CV_Large`. See [mod_CV_fit()] for more details.
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default
#'   is 8.
#' @param future_max_size Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See [ecokit::set_parallel()] and
#'   `future.globals.maxSize` argument of [future::future.options()] for more
#'   details.
#' @param selected_species,excluded_species Character vector or NULL. Names of
#'   species to include or exclude for modelling.
#' @param env_file Character. Path to a file with environment variable
#'   definitions for spatial datasets. Default is `".env"`.
#' @param clamp_pred Logical. Should clamping be applied to sampling efforts and
#'   river length predictors for prediction? Default is `TRUE`.
#' @param fix_efforts,fix_rivers Character or numeric (length 1). Method or
#'   fixed value for sampling effort and river length (both at log-scale) when
#'   clamping is enabled (`clamp_pred = TRUE`). Valid methods: "identity" (use
#'   observed, with no clamping), summary statistics for the sampling efforts
#'   layer ("median", "mean", "max", or "q90" (default; 90th percentile)), or a
#'   single numeric value within observed range.
#' @param climate_models Character vector or "all". Which climate change models
#'   to use for future projections. Valid values (case-sensitive): "GFDL-ESM4",
#'   "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL", or "all"
#'   (default, meaning all available models). If a subset, must be a subset of
#'   the listed valid models.
#' @param climate_scenarios Character vector or "all". Which climate change
#'   scenarios to use for future projections. Valid values: "ssp126", "ssp370",
#'   "ssp585", or "all" (default, meaning all available scenarios). If a subset,
#'   must be a subset of the listed valid scenarios.
#' @param climate_periods Character vector or "all". Time periods for
#'   prediction. Valid values are "2011-2040", "2041-2070", "2071-2100", or
#'   "all" (default), or subset of supported periods.
#' @param copy_maxent_html Logical. Whether to copy the directory containing
#'   HTML results from Maxent to the modelling directory. Default is `TRUE`.
#'
#' @return A tibble summarizing model results for each species, including:
#'   - Evaluation metrics for training and testing data (AUC, TSS, Kappa, etc.)
#'   - Variable importance scores
#'   - Response curves for each environmental variable
#'   - Prediction summaries for current and future climate scenarios
#'   - Paths to generated model files and prediction rasters
#'
#'   Additionally, the function saves various outputs to disk for future use:
#'   - Fitted model objects (as .RData files)
#'   - Extracted model information (evaluation metrics, variable importance,
#'   etc.)
#'   - Prediction rasters for each species, cross-validation fold, and climate
#'   scenario
#'   - Summary statistics across CV folds (mean, weighted mean, SD, and
#'   coefficient of variation)
#'   - Species richness maps for each climate scenario
#'
#' @details The `fit_sdm_models` function orchestrates a comprehensive workflow
#'   that handles all aspects of single-species distribution modelling for
#'   invasive alien plant species in Europe. The workflow integrates several
#'   internal components that manage different stages of the modelling process:
#'
#'   **Overall workflow:**
#'   - *Input validation*: Checks all parameters for validity
#'   - *Data preparation*: Loads and processes model data
#'   - *Parallel processing setup*: Configures computational resources
#'   - *Model fitting and prediction*: For each species and CV fold
#'   - *Results summarization*: Compiles metrics, variable importance, and
#'   predictions
#'   - *Species richness calculation*: Across all modelled species
#'
#'   **Core capabilities:**
#'   - *Data preparation*: The workflow validates and prepares necessary
#'   input data including modelling data, environmental predictors, and
#'   prediction datasets. It handles species selection, data loading, and
#'   preprocessing of spatial predictors (including clamping of sampling efforts
#'   and river length when required).
#'   - *Model parameterization*: The function provides carefully selected
#'   default settings for various SDM algorithms, ensuring consistent
#'   parameterization across models.
#'   - *Model information extraction*: After fitting, the workflow
#'   automatically extracts key information from fitted SDM objects, including
#'   evaluation metrics, variable importance, and response curves.
#'   - *Model optimization*: Technical improvements like optimizing SDM model
#'   object size by setting formula environments to the base environment address
#'   known issues in the sdm package.
#'   - *Parallel prediction*: The workflow efficiently generates predictions
#'   for each species and cross-validation fold, handling model fitting,
#'   information extraction, prediction, and file saving in parallel.
#'   - *Statistical summarization*: Summary statistics are calculated across
#'   cross-validation folds, including mean, weighted mean (by test AUC),
#'   standard deviation, and coefficient of variation of predictions.
#'
#' @author Ahmed El-Gabbas
#' @export

fit_sdm_models <- function(
    sdm_method = NULL, model_settings = NULL, model_dir = NULL,
    hab_abb = NULL, cv_type = "CV_Dist", n_cores = 8L, future_max_size = 2000L,
    selected_species = NULL, excluded_species = NULL, env_file = ".env",
    clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    climate_models = "all", climate_scenarios = "all", climate_periods = "all",
    copy_maxent_html = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  summary_data <- packages <- path_grid <- mod_method <- cv_fold <- preds <-
    pred_mean <- pred_w_mean <- species_name <- method_is_glm <- output_path <-
    evaluation_testing <- auc_test <- summary_prediction_path <- issues <-
    climate_name <- pred_mean_okay <- pred_w_mean_okay <- richness_map <-
    pred_type <- cv <- preds_summ <- NULL

  stringr::str_glue(
    "Fitting models using `{sdm_method}` method and `{cv_type}` ",
    "cross-validation type") %>%
    ecokit::info_chunk(
      cat_bold = TRUE, cat_red = TRUE, line_char_rep = 80L,
      info_lines_before = 2L)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input arguments -------

  ecokit::cat_time("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c(
      "sdm_method", "model_dir", "cv_type", "env_file", "hab_abb"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("clamp_pred", "copy_maxent_html"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = "future_max_size")
  rm(AllArgs, envir = environment())

  # |||||||||||||||||||||||||||||||||||||||||||

  # Environment file
  if (!fs::file_exists(env_file)) {
    ecokit::stop_ctx(
      "Error: Environment file is invalid or does not exist.",
      env_file = env_file, include_backtrace = TRUE)
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # cv_type
  valid_cv_types <- c("CV_Dist", "CV_Large")
  if (!cv_type %in% valid_cv_types) {
    ecokit::stop_ctx(
      "Invalid CV type", cv_type = cv_type, valid_cv_types = valid_cv_types)
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # model method / model_setting

  # rbf is not bounded; see https://github.com/babaknaimi/sdm/issues/42
  valid_sdm_methods <- c(
    "glm", "glmpoly", "gam", "glmnet", "mars", "gbm", "rf", "ranger",
    "cart", "rpart", "maxent", "mlp", "svm", "mda", "fda")
  sdm_method_valid <- any(
    is.null(sdm_method), length(sdm_method) != 1L,
    !is.character(sdm_method), !sdm_method %in% valid_sdm_methods)
  if (sdm_method_valid) {
    ecokit::stop_ctx(
      "Invalid model method", sdm_method = sdm_method,
      class_sdm_method = class(sdm_method),
      valid_sdm_methods = valid_sdm_methods)
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  if (is.null(model_settings)) {
    ecokit::cat_time("Loading default model settings")
    model_settings <- sdm_model_settings()
  }

  if (sdm_method %in% names(model_settings)) {
    model_settings <- model_settings[sdm_method]
  } else {
    ecokit::cat_time(
      "No specific model settings found, using default settings", level = 1L)
    model_settings <- list()
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # n_cores
  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      n_cores < 1L || is.na(n_cores)) {
    ecokit::stop_ctx(
      "n_cores must be a positive integer of length 1",
      n_cores = n_cores, class_n_cores = class(n_cores))
  }
  n_cores <- as.integer(n_cores)
  max_cores <- parallelly::availableCores()
  if (n_cores > max_cores) {
    warning(
      stringr::str_glue(
        "`n_cores` exceeds available cores: {n_cores}. Using all available",
        " cores: {max_cores}"),
      call. = FALSE)
    n_cores <- max_cores
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # future_max_size

  if (is.na(future_max_size) || !is.numeric(future_max_size) ||
      length(future_max_size) != 1L || future_max_size <= 0) {
    ecokit::stop_ctx(
      "future_max_size must be a positive integer of length 1",
      future_max_size = future_max_size,
      class_future_max_size = class(future_max_size))
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # model_dir

  model_dir_valid <- any(
    !inherits(model_dir, "character"),
    length(model_dir) != 1L, !nzchar(model_dir))
  if (model_dir_valid) {
    ecokit::stop_ctx(
      "model_dir must be a character string",
      model_dir = model_dir, class_model_dir = class(model_dir))
  }

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx("Model directory does not exist", model_dir = model_dir)
  }

  output_directory <- fs::path(model_dir, "sdm_models", sdm_method)
  fs::dir_create(output_directory)

  model_results_name <- paste0(sdm_method, "_results")
  model_results_path <- fs::path(
    output_directory, paste0(model_results_name, ".RData"))

  model_summary_path <- fs::path(
    output_directory, paste0(sdm_method, "_summary.qs2"))
  model_richness_path <- fs::path(
    output_directory, paste0(sdm_method, "_species_richness.qs2"))

  if (ecokit::check_data(model_summary_path, warning = FALSE) &&
      ecokit::check_data(model_richness_path, warning = FALSE)) {
    ecokit::cat_time("Models and summary data already exist")
    return(invisible(NULL))
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Environment variables ------

  ecokit::cat_time("Environment variables")
  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  path_grid_r <- fs::path(path_grid, "Grid_10_Land_Crop.RData")

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Processing models -------

  if (ecokit::check_data(model_results_path, warning = FALSE)) {

    # Loading model results -----
    ecokit::cat_time("Loading saved model results")
    model_results <- ecokit::load_as(model_results_path)

  } else {

    ## Define packages to load for each method -----
    pkg_to_load <- tibble::tribble(
      ~mod_method, ~packages,
      "glm", NA_character_,
      "glmpoly", NA_character_,
      "gam", "mgcv",
      "glmnet", "glmnet",
      "mars", "earth",
      "gbm", "gbm",
      "rf", "randomForest",
      "ranger", "ranger",
      "cart", "tree",
      "rpart", "rpart",
      "maxent", c("dismo", "rJava"),
      "mlp", "RSNNS",
      # "rbf", "RSNNS",
      "svm", "kernlab",
      "mda", "mda",
      "fda", "mda") %>%
      dplyr::filter(mod_method == sdm_method) %>%
      dplyr::pull(packages) %>%
      unlist()
    pkg_to_load <- if (all(is.na(pkg_to_load))) NULL else pkg_to_load

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Prepare input data ------

    ecokit::cat_time("Prepare input data")
    input_data <- prepare_input_data(
      model_dir = model_dir, cv_type = cv_type,
      selected_species = selected_species,
      excluded_species = excluded_species, env_file = env_file,
      hab_abb = hab_abb, clamp_pred = clamp_pred, fix_efforts = fix_efforts,
      fix_rivers = fix_rivers, climate_models = climate_models,
      climate_scenarios = climate_scenarios, climate_periods = climate_periods,
      n_cores = n_cores)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::cat_time("Loading modelling data")
    model_data <- ecokit::load_as(input_data$path_species_data)
    cv_folds <- ecokit::load_as(model_data$data_path[[1]]) %>%
      dplyr::pull(cv) %>%
      unique()
    model_data <- dplyr::mutate(
      model_data, method_is_glm = (sdm_method == "glm")) %>%
      tidyr::expand_grid(cv = cv_folds)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Fit models and making predictions ------
    ecokit::cat_time("Fit models and making predictions")

    if (n_cores == 1L) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, nrow(model_data)),
        future_max_size = future_max_size, level = 1L, cat_timestamp = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # Processing model fitting and predictions in parallel -------
    ecokit::cat_time(
      "Processing model fitting and predictions in parallel", level = 1L)

    pkgs_to_load <- c(
      "terra", "stringr", "ecokit", "tibble", "dplyr", "sf",
      "tidyr", "qs2", "fs", "sdm", "purrr", pkg_to_load)
    future_globals <- c(
      "sdm_method", "model_data", "model_settings",
      "input_data", "output_directory", "path_grid_r", "reduce_sdm_formulas",
      "fit_predict_internal", "extract_sdm_info", "copy_maxent_html")

    model_data2 <- ecokit::quietly(
      future.apply::future_lapply(
        X = seq_len(nrow(model_data)),
        FUN = function(line_id) {
          fit_predict_internal(
            line_id = line_id, sdm_method = sdm_method,
            model_data = model_data, model_settings = model_settings,
            input_data = input_data, output_directory = output_directory,
            path_grid_r = path_grid_r, copy_maxent_html = copy_maxent_html)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.packages = pkgs_to_load, future.globals = future_globals))

    ecokit::set_parallel(level = 1L, stop_cluster = TRUE, cat_timestamp = FALSE)
    future::plan("sequential", gc = TRUE)
    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Merge outputs into a single tibble -----

    ecokit::cat_time("Merge outputs into a single tibble", level = 1L)
    model_results <- dplyr::mutate(model_data, data2 = model_data2) %>%
      tidyr::unnest("data2") %>%
      dplyr::select(
        species_name, sdm_method, cv_fold = cv, tidyselect::everything(),
        -method_is_glm)

    ecokit::save_as(
      object = model_results, object_name = model_results_name,
      out_path = model_results_path)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Check model results -----
    ecokit::cat_time("Check model results", level = 1L)
    check_model_results(
      model_results = model_results, n_cores = n_cores,
      future_max_size = future_max_size)

  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Summarize evaluation, response curves, and variable importance -------

  ecokit::cat_time(
    "Summarize evaluation, response curves, and variable importance")

  if (ecokit::check_data(model_summary_path, warning = FALSE)) {

    ecokit::cat_time("Loading summary data", level = 1)
    model_summary <- ecokit::load_as(model_summary_path)

  } else {

    model_summary <- model_results %>%
      dplyr::select(species_name, sdm_method, cv_fold, output_path) %>%
      dplyr::mutate(
        summ_data = purrr::map(
          .x = output_path,
          .f = ~ {
            ecokit::load_as(.x) %>%
              magrittr::extract(
                !names(.) %in% c("fitted_model", "prediction_data")) %>%
              lapply(list) %>%
              tibble::as_tibble()
          })
      ) %>%
      tidyr::unnest("summ_data") %>%
      tidyr::nest(summary_data = -"species_name") %>%
      dplyr::mutate(
        summary1 = purrr::map(
          .x = summary_data,
          .f = ~ {

            ## Evaluation - training ------

            evaluation_training <- dplyr::bind_rows(.x$evaluation_training) %>%
              dplyr::mutate(cv_fold = as.character(cv_fold))
            evaluation_training <- evaluation_training %>%
              dplyr::select(-cv_fold) %>%
              dplyr::group_by(species_name, sdm_method) %>%
              dplyr::summarize_all(mean, na.rm = TRUE) %>%
              dplyr::mutate(cv_fold = "mean") %>%
              dplyr::ungroup() %>%
              dplyr::bind_rows(evaluation_training, .) %>%
              dplyr::arrange(cv_fold)

            # |||||||||||||||||||||||||||||||||||||||||||

            # Evaluation - testing ------

            evaluation_testing <- dplyr::bind_rows(.x$evaluation_testing) %>%
              dplyr::mutate(cv_fold = as.character(cv_fold))
            evaluation_testing <- evaluation_testing %>%
              dplyr::select(-cv_fold) %>%
              dplyr::group_by(species_name, sdm_method) %>%
              dplyr::summarize_all(mean, na.rm = TRUE) %>%
              dplyr::mutate(cv_fold = "mean") %>%
              dplyr::ungroup() %>%
              dplyr::bind_rows(evaluation_testing, .) %>%
              dplyr::arrange(cv_fold)

            # |||||||||||||||||||||||||||||||||||||||||||

            # Variable importance ------

            variable_importance <- dplyr::bind_rows(.x$variable_importance) %>%
              dplyr::mutate(cv_fold = as.character(cv_fold))
            variable_importance <- variable_importance %>%
              dplyr::select(-cv_fold) %>%
              dplyr::group_by(species_name, sdm_method, variable) %>%
              dplyr::summarize_all(mean, na.rm = TRUE) %>%
              dplyr::mutate(cv_fold = "mean") %>%
              dplyr::ungroup() %>%
              dplyr::bind_rows(variable_importance, .) %>%
              dplyr::arrange(cv_fold, variable)

            # |||||||||||||||||||||||||||||||||||||||||||

            # Response curves -------
            response_curves <- dplyr::bind_rows(.x$response_curves) %>%
              dplyr::bind_rows() %>%
              dplyr::mutate(cv_fold = as.character(cv_fold))
            response_curves <- response_curves %>%
              dplyr::select(-cv_fold) %>%
              dplyr::group_by(species_name, sdm_method, variable, x_value) %>%
              dplyr::summarize_all(mean, na.rm = TRUE) %>%
              dplyr::mutate(cv_fold = "mean") %>%
              dplyr::ungroup() %>%
              dplyr::bind_rows(response_curves, .) %>%
              dplyr::arrange(cv_fold, variable, x_value)

            # |||||||||||||||||||||||||||||||||||||||||||

            # Merge data ---------
            tibble::tibble(
              evaluation_training = list(evaluation_training),
              evaluation_testing = list(evaluation_testing),
              variable_importance = list(variable_importance),
              response_curves = list(response_curves))
          })
      ) %>%
      tidyr::unnest("summary1") %>%
      dplyr::select(-"summary_data")

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ### Summarizing predictions in parallel -----

    ecokit::cat_time("Summarizing predictions in parallel")

    dt_auc_test <- model_summary %>%
      dplyr::mutate(auc_test = purrr::map(evaluation_testing, ~.x$auc_test)) %>%
      dplyr::select(species_name, auc_test)
    model_pred_results <- model_results %>%
      dplyr::select(-tidyselect::all_of(c("data_path", "cv_fold"))) %>%
      tidyr::nest(pred_paths = output_path) %>%
      dplyr::left_join(dt_auc_test, by = "species_name")

    if (n_cores == 1L) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, nrow(model_pred_results)),
        future_max_size = future_max_size, level = 1L, cat_timestamp = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    pkgs_to_load <- c(
      "terra", "stringr", "ecokit", "tibble", "dplyr",
      "qs2", "tools", "purrr", "tidyr", "fs")

    ecokit::cat_time("Calculate summary predictions in parallel", level = 1L)

    pred_summary <- ecokit::quietly(
      future.apply::future_lapply(
        X = seq_len(nrow(model_pred_results)),
        FUN = function(line_id) {
          summarize_predictions(
            line_id, model_pred_results, output_directory)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.packages = pkgs_to_load,
        future.globals = c(
          "model_pred_results", "summarize_predictions", "output_directory")))

    ecokit::set_parallel(level = 1L, stop_cluster = TRUE, cat_timestamp = FALSE)
    future::plan("sequential", gc = TRUE)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Prepare and save summary data ----------
    ecokit::cat_time("Prepare and save summary data", level = 1L)
    model_summary <- dplyr::left_join(
      model_summary, dplyr::bind_rows(pred_summary), by = "species_name")

    ecokit::save_as(object = model_summary, out_path = model_summary_path)
    invisible(gc())

  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check for issues in summary maps ------
  ecokit::cat_time("Check for issues in summary maps", level = 1L)

  summ_issues <- model_summary %>%
    dplyr::select(species_name, summary_prediction_path) %>%
    dplyr::mutate(
      preds_summ = purrr::map(
        .x = summary_prediction_path, .f = ecokit::load_as)) %>%
    tidyr::unnest(preds_summ) %>%
    dplyr::select(species_name, climate_name, pred_mean, pred_w_mean) %>%
    dplyr::mutate(
      pred_mean_okay = purrr::map_lgl(
        pred_mean, inherits, "PackedSpatRaster"),
      pred_w_mean_okay = purrr::map_lgl(
        pred_w_mean, inherits, "PackedSpatRaster")) %>%
    dplyr::filter(!pred_mean_okay | !pred_w_mean_okay)

  if (nrow(summ_issues) > 0L) {
    summ_issues_species <- unique(summ_issues$species_name)

    ecokit::cat_sep(
      line_char_rep = 60L, sep_lines_before = 1L,
      line_char = "=", cat_bold = TRUE, cat_red = TRUE)

    paste0(
      "\n!! There are issues in summary maps (",
      length(summ_issues_species), " species) !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_bold = TRUE, cat_timestamp = FALSE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(summ_issues_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    summ_issues %>%
      dplyr::select(-pred_mean, -pred_w_mean) %>%
      tidyr::nest(issues = -climate_name) %>%
      dplyr::mutate(
        message = purrr::map_chr(
          .x = issues,
          .f = ~ {
            tidyr::pivot_longer(
              data = .x,
              cols = c("pred_mean_okay", "pred_w_mean_okay"),
              names_to = "mean_type", values_to = "mean_okay") %>%
              dplyr::filter(!mean_okay) %>%
              dplyr::mutate(
                mean_type = dplyr::if_else(
                  mean_type == "pred_mean_okay", "mean", "w_mean")) %>%
              dplyr::select(-mean_okay) %>%
              tidyr::nest(message_int = -mean_type) %>%
              dplyr::mutate(
                message_int = purrr::map2_chr(
                  .x = mean_type, .y = message_int,
                  .f = function(type, name) {
                    mean_name <- dplyr::if_else(
                      type == "w_mean", "weighted mean", "mean")
                    paste0(
                      crayon::blue(crayon::bold(mean_name)), " (",
                      paste(unique(unlist(name)), collapse = "; "), ")")
                  })
              ) %>%
              dplyr::pull(message_int) %>%
              paste(collapse = " --- ") %>%
              paste0("  >>>  ", .)
          }),
        message = purrr::map2_chr(
          .x = climate_name, .y = message,
          .f = ~ {
            paste0(crayon::bold(crayon::red(.x)), "\n", .y) %>%
              paste(collapse = "\n")
          })) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n") %>%
      ecokit::cat_time(cat_timestamp = FALSE)

    ecokit::cat_sep(
      line_char_rep = 60L, sep_lines_before = 1L, sep_lines_after = 2L,
      line_char = "=", cat_bold = TRUE, cat_red = TRUE)
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Calculate species richness at each climate option -------

  ecokit::cat_time("Calculate species richness at each climate option")

  if (nrow(summ_issues) > 0L) {
    sr_options <- vector(mode = "character")
    if (all(summ_issues$pred_mean_okay)) sr_options <- c(sr_options, "mean")
    if (all(summ_issues$pred_w_mean_okay)) sr_options <- c(sr_options, "w_mean")

    if (length(sr_options) == 1L) {
      ecokit::cat_time(
        paste0(
          "Calculating species richness only for ", sr_options, " predictions"),
        level = 1L, cat_timestamp = FALSE)

    } else if (length(sr_options) == 0L) {
      ecokit::cat_time(
        "Species richness will not be calculated",
        level = 1L, cat_timestamp = FALSE)
    }
  } else {
    sr_options <- c("mean", "w_mean")
  }


  if (length(sr_options) > 0L) {

    exclude_cols <- c(
      "pred_sd", "pred_cov", "time_period", "climate_model", "climate_scenario")

    species_richness_r <- model_summary$summary_prediction_path %>%
      purrr::map_dfr(ecokit::load_as) %>%
      dplyr::select(-tidyselect::all_of(exclude_cols)) %>%
      tidyr::pivot_longer(
        cols = c("pred_mean", "pred_w_mean"), names_to = "pred_type",
        values_to = "pred_value") %>%
      dplyr::filter(pred_type %in% paste0("pred_", sr_options)) %>%
      tidyr::nest(
        preds = "pred_value", .by = c("climate_name", "pred_type")) %>%
      dplyr::mutate(
        richness_map = purrr::map2(
          .x = preds, .y = climate_name,
          .f = ~{
            purrr::map(unlist(.x), terra::unwrap) %>%
              terra::rast() %>%
              terra::app(sum, na.rm = TRUE) %>%
              stats::setNames(paste0("richness_", .y)) %>%
              terra::wrap()
          })) %>%
      dplyr::select(-preds) %>%
      tidyr::pivot_wider(names_from = "pred_type", values_from = richness_map)

    ecokit::save_as(object = species_richness_r, out_path = model_richness_path)

  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, cat_bold = TRUE, cat_red = TRUE,
    prefix = paste0("\nProcessing ", sdm_method, " models took "))

  return(invisible(NULL))

}
