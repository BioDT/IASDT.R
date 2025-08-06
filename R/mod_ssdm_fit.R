
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
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default
#'   is 8.
#' @param model_dir Character. Path to the directory containing model data and
#'   where outputs and results will be saved. Model data are prepared using the
#'   [mod_prepare_HPC()] and [mod_prepare_data()] functions.
#' @param cv_type Character. Cross-validation type. One of `CV_Dist` (default)
#'   or `CV_Large`. See [mod_CV_fit()] for more details.
#' @param selected_species Character vector or NULL. Names of species to include
#'   for modelling.
#' @param excluded_species Character vector or NULL. Names of species to exclude
#'   from modelling.
#' @param env_file Character. Path to a file with environment variable
#'   definitions for spatial datasets. Default is `".env"`.
#' @param hab_abb Character. Abbreviation for a single SynHab habitat type.
#'   Valid values: "0", "1", "2", "3", "4a", "4b", "10", "12a", "12b". See
#'   [mod_prepare_HPC()] for more details.
#' @param clamp_pred Logical. Should clamping be applied to sampling efforts and
#'   river length predictors for prediction? Default is `TRUE`.
#' @param fix_efforts Character or numeric (length 1). Method or fixed value for
#'   sampling effort (log-scale) when clamping is enabled (`clamp_pred = TRUE`).
#'   Valid methods: "identity" (use observed, with no clamping), summary
#'   statistics for the sampling efforts layer ("median", "mean", "max", or
#'   "q90" (default; 90th percentile)), or a single numeric value within
#'   observed range. If using a numeric value, ensure it is within the range of
#'   sampling effort layer.
#' @param fix_rivers Character or numeric (length 1). Method or fixed value for
#'   river length (log-scale) when clamping is enabled (`clamp_pred = TRUE`).
#'   Valid values are similar to `fix_efforts`.
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
    sdm_method = NULL, model_settings = NULL, n_cores = 8,
    model_dir = NULL, cv_type = "CV_Dist", selected_species = NULL,
    excluded_species = NULL, env_file = ".env", hab_abb = NULL,
    clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    climate_models = "all", climate_scenarios = "all",
    climate_periods = "all") {

  .start_time <- lubridate::now(tzone = "CET")

  summary_data <- model_formula <- sdm_data <- . <- method_is_glm <-
    packages <- path_grid <- mod_method <- maps <- cv_fold <- pred_dir <- NULL

  stringr::str_glue(
    "Fitting models using `{sdm_method}` method and `{cv_type}` ",
    "cross-validation type") %>%
    ecokit::info_chunk(
      cat_bold = TRUE, cat_red = TRUE, line_char_rep = 80L,
      info_lines_before = 2)

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
    args_all = AllArgs, args_type = "logical", args_to_check = "clamp_pred")
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
    is.null(sdm_method), length(sdm_method) != 1,
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
      "No specific model settings found, using default settings", level = 1)
    model_settings <- list()
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # n_cores
  if (!is.numeric(n_cores) || length(n_cores) != 1 ||
      n_cores < 1 || is.na(n_cores)) {
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

  # model_dir

  model_dir_valid <- any(
    !inherits(model_dir, "character"),
    length(model_dir) != 1, !nzchar(model_dir))
  if (model_dir_valid) {
    ecokit::stop_ctx(
      "model_dir must be a character string",
      model_dir = model_dir, class_model_dir = class(model_dir))
  }

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx("Model directory does not exist", model_dir = model_dir)
  }

  output_directory <- fs::path(model_dir, "sdm_models", sdm_method)
  model_results_dir <- fs::path(output_directory, "model_results")
  fs::dir_create(model_results_dir)

  model_results_name <- paste0(sdm_method, "_results")
  model_results_path <- fs::path(
    output_directory, paste0(model_results_name, ".RData"))

  model_summary_name <- paste0(sdm_method, "_summary")
  model_summary_path <- fs::path(
    output_directory, paste0(model_summary_name, ".RData"))

  if (ecokit::check_data(model_summary_path, warning = FALSE)) {
    ecokit::cat_time("Loading saved model results")
    return(ecokit::load_as(model_summary_path))
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
      climate_scenarios = climate_scenarios, climate_periods = climate_periods)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::cat_time("Loading modelling data")
    model_data <- ecokit::load_as(input_data$path_species_data) %>%
      dplyr::filter(method_is_glm == (sdm_method == "glm")) %>%
      dplyr::select(-method_is_glm)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Create prediction paths -------

    ecokit::cat_time("Create prediction paths")
    paste0(
      "pred_",
      ecokit::load_as(input_data$path_prediction_data)$climate_name) %>%
      fs::path(output_directory, .) %>%
      fs::dir_create()

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Fit models and making predictions ------
    ecokit::cat_time("Fit models and making predictions")

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, nrow(model_data)), future_max_size = 800L,
        level = 1, cat_timestamp = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # Processing model fitting and predictions in parallel -------
    ecokit::cat_time(
      "Processing model fitting and predictions in parallel", level = 1)

    pkgs_to_load <- c(
      "terra", "stringr", "ecokit", "tibble", "dplyr", "sf", "tidyr",
      "fs", "sdm", "purrr", pkg_to_load)
    future_globals <- c(
      "sdm_method", "model_data", "model_settings", "model_results_dir",
      "input_data", "output_directory", "path_grid_r", "reduce_sdm_formulas",
      "fit_predict_internal", "extract_sdm_info")

    model_data2 <- withCallingHandlers(
      suppressPackageStartupMessages(
        future.apply::future_lapply(
          X = seq_len(nrow(model_data)),
          FUN = function(line_id) {
            fit_predict_internal(
              line_id = line_id, sdm_method = sdm_method,
              model_data = model_data, model_settings = model_settings,
              model_results_dir = model_results_dir,
              input_data = input_data, output_directory = output_directory,
              path_grid_r = path_grid_r)
          },
          future.scheduling = Inf, future.seed = TRUE,
          future.packages = pkgs_to_load, future.globals = future_globals)
      ),
      warning = function(w) {
        if (grepl(
          "was built under R version", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })

    ecokit::set_parallel(
      n_cores = n_cores, level = 1, stop_cluster = TRUE, cat_timestamp = FALSE)
    future::plan("sequential", gc = TRUE)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Merge outputs into a single tibble -----

    ecokit::cat_time("Merge outputs into a single tibble", level = 1)
    model_results <- dplyr::mutate(model_data, data2 = model_data2) %>%
      tidyr::unnest("data2")

    ecokit::save_as(
      object = model_results, object_name = model_results_name,
      out_path = model_results_path)
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Summarize evaluation, response curves, and variable importance -------

  ecokit::cat_time(
    "Summarize evaluation, response curves, and variable importance")

  model_summary <- model_results %>%
    dplyr::select(-model_formula, -sdm_data) %>%
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
    tidyr::unnest("summary1")

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ### Summarizing predictions in parallel -----

  ecokit::cat_time("Summarizing predictions in parallel")

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(model_summary)), future_max_size = 800L,
      level = 1, cat_timestamp = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  pkgs_to_load <- c(
    "terra", "stringr", "ecokit", "tibble", "dplyr",
    "tools", "purrr", "tidyr", "fs")

  ecokit::cat_time("Calculate summary predictions in parallel", level = 1)
  pred_summary <- withCallingHandlers(
    suppressPackageStartupMessages(
      future.apply::future_lapply(
        X = seq_len(nrow(model_summary)),
        FUN = function(line_id) {
          summarize_predictions(line_id, model_summary)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.packages = pkgs_to_load,
        future.globals = c("model_summary", "summarize_predictions"))
    ),
    warning = function(w) {
      if (grepl(
        "was built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })

  ecokit::set_parallel(
    n_cores = n_cores, level = 1, stop_cluster = TRUE, cat_timestamp = FALSE)
  future::plan("sequential", gc = TRUE)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare and save summary data ----------
  ecokit::cat_time("Prepare and save summary data", level = 1)
  model_summary <- model_summary %>%
    dplyr::mutate(prediction_summary = pred_summary) %>%
    dplyr::select(-summary_data)
  ecokit::save_as(
    object = model_summary, object_name = model_summary_name,
    out_path = model_summary_path)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Calculate species richness at each climate option

  ecokit::cat_time("Calculate species richness at each climate option")
  cols_to_remove <- c("data_okay", "tif_path", "tif_okay", "species_name")

  richness_summary <- dplyr::bind_rows(model_summary$prediction_summary) %>%
    dplyr::filter(cv_fold %in% c("mean", "weighted_mean")) %>%
    dplyr::select(-tidyselect::all_of(cols_to_remove)) %>%
    tidyr::nest(.key = "maps", .by = -"data_path") %>%
    dplyr::mutate(
      richness_map = purrr::pmap(
        .l = list(maps, pred_dir, cv_fold),
        .f = function(maps, pred_dir, cv_fold) {

          richness_name <- paste0(sdm_method, "_species_richness_", cv_fold)
          data_path <- fs::path(pred_dir, paste0(richness_name, ".RData"))
          tiff_path <- fs::path(pred_dir, paste0(richness_name, ".tif"))

          richness <- unlist(maps) %>%
            purrr::map(ecokit::load_as, unwrap_r = TRUE) %>%
            terra::rast() %>%
            terra::app(sum, na.rm = TRUE)

          terra::writeRaster(
            x = richness, overwrite = TRUE, filename = tiff_path,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
          ecokit::save_as(
            object = terra::wrap(richness), object_name = richness_name,
            out_path = data_path)
        }))

  rm(richness_summary)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, cat_bold = TRUE, cat_red = TRUE,
    prefix = paste0("\nProcessing ", sdm_method, " models took "),
    ... = "\n")

  model_summary

}
