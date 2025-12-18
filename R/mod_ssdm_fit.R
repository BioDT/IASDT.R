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
#'   models. Valid values: "glm", "glmpoly", "gam", "glmnet/glmnet2",
#'   "mars/mars2", "gbm/gbm2", "rf/rf2", "ranger/ranger2", "cart", "rpart",
#'   "maxent", "mlp", "rbf", "svm/svm2", "mda/mda2", and "fda/fda2". These
#'   correspond to selected methods supported by the `sdm` package. For details
#'   and supported options, see [sdm::getmethodNames()]. Note that some methods
#'   have custom implementations (e.g., "glmnet2", "gbm2", "mars2", "ranger2",
#'   "rf2", "svm2", "mda2", "fda2") to ensure consistent parameterisation and
#'   performance across models.
#' @param model_settings List or NULL. List of model-specific settings. If
#'   `NULL`, defaults to custom settings defined within the workflow.
#' @param model_dir Character. Path to the directory containing model data and
#'   where outputs and results will be saved. Model data are prepared using the
#'   [mod_prepare_hpc()] and [mod_prepare_data()] functions.
#' @param hab_abb Character. Abbreviation for a single SynHab habitat type.
#'   Valid values: "0", "1", "2", "3", "4a", "4b", "10", "12a", "12b". See
#'   [mod_prepare_hpc()] for more details.
#' @param cv_type Character. Cross-validation type. One of `cv_dist` (default)
#'   or `cv_large`. See [mod_cv_fit()] for more details.
#' @param n_cores,n_cores_check,n_cores_summary Integer. Number of CPU cores for
#'   parallel processing of model fitting, model checking, and summarising model
#'   outputs. Default is 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
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
    hab_abb = NULL, cv_type = "cv_dist", n_cores = 8L, n_cores_check = 8L,
    n_cores_summary = 8L, strategy = "multisession", future_max_size = 2000L,
    selected_species = NULL, excluded_species = NULL, env_file = ".env",
    clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    climate_models = "all", climate_scenarios = "all", climate_periods = "all",
    copy_maxent_html = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  terra_temp <- fs::path_temp(
    paste0("terra_temp_", format(Sys.time(), "%H%M%S")))
  fs::dir_create(terra_temp)
  withr::defer(try({
    fs::dir_delete(terra_temp)
    terra::tmpFiles(current = TRUE, remove = TRUE)
  }, silent = TRUE))

  terra::terraOptions(
    # fraction of RAM terra may use (0-0.9)
    memfrac = 0.1,
    # (GB) below which mem is assumed available
    memmin = 1L,
    # (GB) cap for terra
    memmax = 10L,
    # silence per-worker progress bars
    progress = 0L,
    todisk = TRUE,
    tempdir = terra_temp)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input arguments -------

  ecokit::check_args(
    args_to_check = c("sdm_method", "model_dir", "cv_type"),
    args_type = "character")
  ecokit::check_args(
    args_to_check = c("clamp_pred", "copy_maxent_html"), args_type = "logical")
  ecokit::check_args(args_to_check = "future_max_size", args_type = "numeric")

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  n_cores <- .validate_n_cores(n_cores)
  n_cores_check <- .validate_n_cores(n_cores_check)
  n_cores_summary <- .validate_n_cores(n_cores_summary)
  strategy <- .validate_strategy(strategy)

  if (strategy == "sequential") {
    n_cores <- n_cores_check <- n_cores_summary <- 1L
  }

  if (future_max_size <= 0) {
    ecokit::stop_ctx(
      "future_max_size must be a positive integer of length 1",
      future_max_size = future_max_size)
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # cv_type
  valid_cv_types <- c("cv_dist", "cv_large")
  if (!cv_type %in% valid_cv_types) {
    ecokit::stop_ctx(
      "Invalid CV type", cv_type = cv_type, valid_cv_types = valid_cv_types)
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  packages <- path_grid <- mod_method <- cv_fold <- preds <- pred_mean <-
    pred_w_mean <- species_name <- time_period <- climate_model <- pred <-
    climate_scenario <- climate_name <- pred_mean_okay <- pred_w_mean_okay <-
    cv <- preds_summ <- n_rep <- NULL

  stringr::str_glue(
    "Fitting models using `{sdm_method}` and `{cv_type}` ",
    "cross-validation type") %>%
    ecokit::info_chunk(cat_bold = TRUE, cat_red = TRUE, line_char_rep = 80L)

  # Ensure that all methods are registered
  methods_to_copy <- c(
    "fda2", "mda2", "glmnet2", "mars2", "ranger2", "rf2", "svm2")

  if (sdm_method %in% methods_to_copy) {
    sdm_path <- system.file("methods/sdm", package = "sdm")
    method_file <- system.file(paste0(sdm_method, ".R"), package = "IASDT.R")

    if (!nzchar(method_file)) {
      ecokit::stop_ctx(
        paste0("`", method_file, "` method file not found in IASDT.R package"))
    }
    if (!fs::file_exists(method_file)) {
      ecokit::stop_ctx(
        paste0("`", method_file, "` method file does not exist"),
        method_file = method_file)
    }

    fs::file_copy(path = method_file, new_path = sdm_path, overwrite = TRUE)
  }

  # |||||||||||||||||||||||||||||||||||||||||||

  # model method / model_setting

  # rbf is not bounded; see https://github.com/babaknaimi/sdm/issues/42
  valid_sdm_methods <- c(
    "glm", "glmpoly", "gam", "glmnet", "glmnet2", "mars", "mars2", "gbm",
    "gbm2", "rf", "rf2", "ranger", "ranger2", "cart", "rpart", "maxent", "mlp",
    "svm", "svm2", "mda", "mda2", "fda", "fda2")
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

  # model_dir

  model_dir_invalid <- any(
    !inherits(model_dir, "character"),
    length(model_dir) != 1L, !nzchar(model_dir))

  if (model_dir_invalid) {
    ecokit::stop_ctx(
      "model_dir must be a character string",
      model_dir = model_dir, class_model_dir = class(model_dir))
  }

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx("Model directory does not exist", model_dir = model_dir)
  }

  output_directory <- fs::path(model_dir, sdm_method)
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

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  path_grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!fs::file_exists(path_grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist",
      path_grid_r = path_grid_r)
  }

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
      "glmnet2", "glmnet",
      "mars", "earth",
      "mars2", "earth",
      "gbm", "gbm",
      "gbm2", "gbm",
      "gbm2", "dismo",
      "rf", "randomForest",
      "rf2", "randomForest",
      "ranger", "ranger",
      "ranger2", "ranger",
      "cart", "tree",
      "rpart", "rpart",
      "maxent", c("dismo", "rJava"),
      "mlp", "RSNNS",
      # "rbf", "RSNNS",
      "svm", "kernlab",
      "svm2", "e1071",
      "mda", "mda",
      "mda2", "mda",
      "fda", "mda",
      "fda2", "mda") %>%
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
    cv_folds <- ecokit::load_as(model_data$species_data[[1]]) %>%
      dplyr::pull(cv) %>%
      unique()
    model_data <- tidyr::expand_grid(model_data, cv = cv_folds)
    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Fit models and making predictions ------
    ecokit::cat_time("Fit models and making predictions")

    if (n_cores == 1L) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, nrow(model_data)), strategy = strategy,
        future_max_size = future_max_size, level = 1L)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # Processing model fitting and predictions in parallel -------
    ecokit::cat_time(
      "Processing model fitting and predictions in parallel", level = 1L)

    pkgs_to_load <- unique(
      c(
        "terra", "stringr", "ecokit", "tibble", "dplyr", "sf", "rlang",
        "tidyselect", "tidyr", "qs2", "fs", "sdm", "purrr", "methods", "stats",
        "magrittr", pkg_to_load))

    future_globals <- c(
      "sdm_method", "model_data", "model_settings",
      "input_data", "output_directory", "path_grid_r", "reduce_sdm_formulas",
      "fit_predict_internal", "extract_sdm_info", "copy_maxent_html")

    model_data2 <- future.apply::future_lapply(
      X = seq_len(nrow(model_data)),
      FUN = function(line_id) {
        ecokit::quietly(
          fit_predict_internal(
            line_id = line_id, sdm_method = sdm_method,
            model_data = model_data, model_settings = model_settings,
            input_data = input_data, output_directory = output_directory,
            path_grid_r = path_grid_r, copy_maxent_html = copy_maxent_html),
          "fitted probabilities numerically",
          "Using formula.+ is deprecated when")
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkgs_to_load, future.globals = future_globals) %>%
      dplyr::bind_rows()

    ecokit::set_parallel(level = 1L, stop_cluster = TRUE)
    future::plan("sequential", gc = TRUE)
    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Merge outputs into a single tibble -----

    ecokit::cat_time("Merge outputs into a single tibble", level = 1L)

    model_results <- dplyr::left_join(
      model_data, model_data2, by = c("species_name", "cv")) %>%
      dplyr::rename(cv_fold = cv) %>%
      dplyr::select(
        tidyselect::all_of(c("species_name", "sdm_method", "cv_fold")),
        tidyselect::everything())

    ecokit::save_as(
      object = model_results, object_name = model_results_name,
      out_path = model_results_path)

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Check model results -----
    ecokit::cat_time("Check model results", level = 1L)
    check_model_results(
      model_results = model_results, n_cores_check = n_cores_check,
      strategy = strategy, future_max_size = future_max_size)

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
      dplyr::group_by(species_name, sdm_method, n_rep) %>%
      dplyr::summarise(model_data = list(model_data), .groups = "drop")

    if (n_cores_summary == 1L) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores_summary, nrow(model_summary)),
        strategy = strategy, future_max_size = future_max_size, level = 1L)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    pkgs_to_load <- c(
      "magrittr", "ecokit", "tibble", "dplyr", "qs2", "tidyselect",
      "stringr", "terra", "purrr", "fs", "withr", "tidyr", "stats")

    model_summary_0 <- future.apply::future_lapply(
      X = seq_len(nrow(model_summary)),
      FUN = function(line_id) {
        ecokit::quietly(
          summarise_ssdms(
            line_id = line_id, output_directory = output_directory,
            model_summary = model_summary))
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkgs_to_load,
      future.globals = c(
        "model_summary", "output_directory", "summarise_ssdms")) %>%
      dplyr::bind_rows()

    ecokit::set_parallel(level = 1L, stop_cluster = TRUE)
    future::plan("sequential", gc = TRUE)

    model_summary <- dplyr::bind_cols(model_summary, model_summary_0) %>%
      dplyr::select(-tidyselect::all_of("model_data"))

    rm(model_summary_0, envir = environment())
    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Prepare and save summary data ----------
    ecokit::cat_time("Prepare and save summary data", level = 1L)
    ecokit::save_as(object = model_summary, out_path = model_summary_path)
    invisible(gc())

  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check for issues in summary maps ------
  ecokit::cat_time("Check for issues in summary maps", level = 1L)

  summ_issues <- model_summary %>%
    dplyr::select(tidyselect::all_of(c("species_name", "pred"))) %>%
    dplyr::mutate(
      preds_summ = purrr::map(.x = pred, .f = ecokit::load_as)) %>%
    tidyr::unnest(preds_summ) %>%
    dplyr::select(
      tidyselect::all_of(
        c("species_name", "climate_name", "pred_mean", "pred_w_mean"))) %>%
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

    model_summary_proj <- model_summary$pred
    rm(model_summary, envir = environment())
    invisible(gc())

    if (n_cores == 1L) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = min(n_cores, length(model_summary_proj)),
        strategy = strategy, future_max_size = future_max_size, level = 1L)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    pkgs_to_load <- c(
      "magrittr", "ecokit", "dplyr", "qs2", "tidyselect", "terra", "fs")

    species_richness_r <- future.apply::future_lapply(
      X = model_summary_proj,
      FUN = function(pred) {
        dplyr::select(
          ecokit::load_as(pred), -tidyselect::all_of(c("pred_sd", "pred_cov")))
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkgs_to_load)

    ecokit::set_parallel(level = 1L, stop_cluster = TRUE)
    future::plan("sequential", gc = TRUE)
    invisible(gc())

    calc_sr <- function(preds, climate_name, sdm_method) {
      map_name <- paste0("richness_", sdm_method, "_", climate_name)
      purrr::map(preds, terra::unwrap) %>%
        terra::rast() %>%
        terra::app(sum, na.rm = TRUE) %>%
        stats::setNames(map_name) %>%
        terra::toMemory() %>%
        terra::wrap()
    }

    species_richness_r <- dplyr::bind_rows(species_richness_r) %>%
      dplyr::arrange(time_period, climate_name) %>%
      dplyr::group_by(
        time_period, climate_model, climate_scenario, climate_name) %>%
      dplyr::summarize(
        pred_mean = list(pred_mean), pred_w_mean = list(pred_w_mean),
        .groups = "drop") %>%
      dplyr::mutate(
        pred_mean = purrr::map2(
          .x = pred_mean, .y = climate_name,
          .f = calc_sr, sdm_method = sdm_method))
    invisible(gc())

    species_richness_r <- dplyr::mutate(
      species_richness_r,
      pred_w_mean = purrr::map2(
        .x = pred_w_mean, .y = climate_name,
        .f = calc_sr, sdm_method = sdm_method))
    invisible(gc())

    ecokit::save_as(object = species_richness_r, out_path = model_richness_path)

    rm(species_richness_r, envir = environment())
    invisible(gc())

  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, cat_bold = TRUE, cat_red = TRUE,
    prefix = paste0("\nProcessing ", sdm_method, " models took "))

  return(invisible(NULL))

}
