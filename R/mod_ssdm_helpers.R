# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# extract_sdm_info ------
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' Extract SDM Model Information and Evaluation Metrics
#'
#' This function extracts key information, evaluation metrics, variable
#' importance, and response curves from a fitted species distribution model
#' (SDM) object of class `sdmModels`, or from a file path to a saved model. It
#' is designed as a component of the `IASDT.R` workflow for modelling the
#' distribution of invasive alien plant species in Europe. The function works
#' with models fitted using the `sdm` package and provides a tidy summary of
#' model performance for both training and independent test data, as well as
#' variable importance and response curves.
#'
#' @param model An object of class `sdmModels` (from `sdm` package) or a
#'   character path to a saved model file (produced by `sdm`). Required.
#' @param cv_fold Integer or character. Identifier for the cross-validation
#'   fold. Required, must not be `NULL`.
#'
#' @return A list with the following elements:
#' - evaluation_training`: A tibble with evaluation metrics for the training
#'   data, including AUC, TSS, Kappa, prevalence, sensitivity, specificity, and
#'   other statistics.
#' - evaluation_testing`: A tibble with evaluation metrics for the independent
#'   test data, with the same structure as `evaluation_training`.
#' - variable_importance`: A tibble summarizing the importance of each predictor
#'   variable, including correlation and AUC metrics for the test data.
#' - response_curves`: A tibble containing response curves for each predictor
#'   variable, showing the predicted response across the range of values.
#'
#'   If evaluation or variable importance data are missing, the function returns
#'   tibbles with `NA` values for the corresponding metrics. The function also
#'   handles errors gracefully and provides informative messages if extraction
#'   fails.
#'
#' @details The function is not expected to be called directly by users, but is
#'   used internally within the `IASDT.R` workflow, particularly the
#'   [fit_sdm_models()] function.
#'
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

extract_sdm_info <- function(model = NULL, cv_fold = NULL) {

  AUCtest <- corTest <- criteria <- variables <- . <- NULL

  if (is.null(cv_fold)) {
    ecokit::stop_ctx("`cv_fold` can not be empty")
  }

  if (is.null(model)) {
    ecokit::stop_ctx("`model` can not be empty")
  }

  # Validate and load model
  model_fail <- FALSE
  if (!inherits(model, "sdmModels")) {
    if (inherits(model, "character")) {
      if (ecokit::check_data(model, warning = FALSE)) {
        model <- ecokit::load_as(model)
      } else {
        model_fail <- TRUE
      }
    } else {
      model_fail <- TRUE
    }
  }

  if (model_fail) {
    ecokit::stop_ctx(
      "model must be a sdmModels object or a file path",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model, "setting")) {
    ecokit::stop_ctx(
      "model does not have a 'setting' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting, "featureFrame")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@featureFrame' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting@featureFrame, "predictors")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@featureFrame@predictors' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting, "methods")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@methods' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model, "models")) {
    ecokit::stop_ctx(
      "model does not have a 'models' slot.",
      model = model, class_model = class(model))
  }

  # Ensure predictor_names is a character vector
  predictor_names <- as.character(model@setting@featureFrame@predictors)
  if (!is.character(predictor_names)) {
    ecokit::stop_ctx(
      "'model@setting@featureFrame@predictors' must be a character vector.",
      predictors = predictor_names, class_predictors = class(predictor_names))
  }

  sdm_method <- model@setting@methods
  model_info <- tryCatch(
    sdm::getModelInfo(model),
    error = function(e) {
      ecokit::stop_ctx(
        "Failed to get model info from sdmModels object.",
        error_message = e$message)
    })

  species_name <- as.character(model_info$species)

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract model metadata
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  basic_info <- tibble::tibble(
    species_name = species_name, sdm_method = sdm_method, cv_fold)

  species_model <- tryCatch(
    model@models[[1L]][[1L]][[1L]],
    error = function(e) NULL)

  if (is.null(species_model)) {
    ecokit::stop_ctx("No fitted model found in sdmModels object.")
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation ----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # columns to keep in evaluation data
  evaluation_sort <- c(
    "TSS", "Kappa", "threshold", "prevalence", "sensitivity", "specificity")

  # Create empty evaluation tibble to handle cases where evaluation is missing
  # or not available
  empty_evaluation_train <- tibble::tibble(
    prevalence_train = NA_real_, auc_train = NA_real_,
    boyce_train = NA_real_, cor_train = NA_real_,
    cor_p_train = NA_real_, deviance_train = NA_real_,
    tss_ess_train = NA_real_, kappa_ess_train = NA_real_,
    threshold_ess_train = NA_real_, prevalence_ess_train = NA_real_,
    sensitivity_ess_train = NA_real_, specificity_ess_train = NA_real_,
    tss_mss_train = NA_real_, kappa_mss_train = NA_real_,
    threshold_mss_train = NA_real_, prevalence_mss_train = NA_real_,
    sensitivity_mss_train = NA_real_, specificity_mss_train = NA_real_) %>%
    dplyr::bind_cols(basic_info, .)

  empty_evaluation_test <- empty_evaluation_train %>%
    dplyr::rename_with(~ stringr::str_replace(., "train", "test"))

  # Ensure that `evaluation` slot exists
  if (methods::.hasSlot(species_model, "evaluation")) {

    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ## Evaluation - training ------
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if ("training" %in% names(species_model@evaluation)) {

      eval_train <- species_model@evaluation$training

      # `Statistics` slot
      if (methods::.hasSlot(eval_train, "statistics")) {
        eval_train_stats <- eval_train@statistics
        if ("cBoyce" %in% names(eval_train_stats)) {
          cont_boyce_train <- eval_train_stats$cBoyce
        } else {
          cont_boyce_train <- NA_real_
        }
        eval_train_stats <- tibble::tibble(
          prevalence_train = eval_train_stats$Prevalence,
          auc_train = eval_train_stats$AUC,
          boyce_train = cont_boyce_train,
          cor_train = eval_train_stats$COR[1L],
          cor_p_train = eval_train_stats$COR[2L],
          deviance_train = eval_train_stats$Deviance)
      } else {
        # Return empty tibble if no `statistics` slot exists
        eval_train_stats <- tibble::tibble(
          prevalence_train = NA_real_, auc_train = NA_real_,
          boyce_train = NA_real_, cor_train = NA_real_, cor_p_train = NA_real_,
          deviance_train = NA_real_)
      }

      # `threshold_based` slot
      if (methods::.hasSlot(eval_train, "threshold_based")) {
        eval_train_thr <- eval_train@threshold_based
        eval_train_thr_ess <- eval_train_thr %>%
          dplyr::filter(criteria == "sp=se") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_ess_train"))
        eval_train_thr_mss <- eval_train_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_mss_train"))
      } else {
        eval_train_thr_ess <- tibble::tibble(
          tss_ess_train = NA_real_, kappa_ess_train = NA_real_,
          threshold_ess_train = NA_real_, prevalence_ess_train = NA_real_,
          sensitivity_ess_train = NA_real_, specificity_ess_train = NA_real_)
        eval_train_thr_mss <- tibble::tibble(
          tss_mss_train = NA_real_, kappa_mss_train = NA_real_,
          threshold_mss_train = NA_real_, prevalence_mss_train = NA_real_,
          sensitivity_mss_train = NA_real_, specificity_mss_train = NA_real_)
      }

      evaluation_training <- dplyr::bind_cols(
        basic_info, eval_train_stats, eval_train_thr_ess, eval_train_thr_mss)

    } else {
      # no training evaluation data, return empty tibble
      evaluation_training <- empty_evaluation_train
    }

    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Evaluation - testing -----
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if ("test.indep" %in% names(species_model@evaluation)) {

      eval_test <- species_model@evaluation$test.indep

      # `Statistics` slot
      if (methods::.hasSlot(eval_test, "statistics")) {
        eval_test_stats <- eval_test@statistics
        if ("cBoyce" %in% names(eval_test_stats)) {
          cont_boyce_test <- eval_test_stats$cBoyce
        } else {
          cont_boyce_test <- NA_real_
        }
        eval_test_stats <- tibble::tibble(
          prevalence_test = eval_test_stats$Prevalence,
          auc_test = eval_test_stats$AUC,
          boyce_test = cont_boyce_test,
          cor_test = eval_test_stats$COR[1L],
          cor_p_test = eval_test_stats$COR[2L],
          deviance_test = eval_test_stats$Deviance)
      } else {
        # Return empty tibble if no `statistics` slot exists
        eval_test_stats <- tibble::tibble(
          prevalence_test = NA_real_, auc_test = NA_real_,
          boyce_test = NA_real_, cor_test = NA_real_, cor_p_test = NA_real_,
          deviance_test = NA_real_)
      }

      # `threshold_based` slot
      if (methods::.hasSlot(eval_test, "threshold_based")) {
        eval_test_thr <- eval_test@threshold_based
        eval_test_thr_ess <- eval_test_thr %>%
          dplyr::filter(criteria == "sp=se") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_ess_test"))
        eval_test_thr_mss <- eval_test_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_mss_test"))
      } else {
        eval_test_thr_ess <- tibble::tibble(
          tss_ess_test = NA_real_, kappa_ess_test = NA_real_,
          threshold_ess_test = NA_real_, prevalence_ess_test = NA_real_,
          sensitivity_ess_test = NA_real_, specificity_ess_test = NA_real_)
        eval_test_thr_mss <- tibble::tibble(
          tss_mss_test = NA_real_, kappa_mss_test = NA_real_,
          threshold_mss_test = NA_real_, prevalence_mss_test = NA_real_,
          sensitivity_mss_test = NA_real_, specificity_mss_test = NA_real_)
      }

      evaluation_testing <- dplyr::bind_cols(
        basic_info, eval_test_stats, eval_test_thr_ess, eval_test_thr_mss)

    } else {
      # no testing evaluation data, return empty tibble
      evaluation_testing <- empty_evaluation_test
    }

  } else {
    evaluation_training <- empty_evaluation_train
    evaluation_testing <- empty_evaluation_test
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Variable importance -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  empty_var_imp <- dplyr::bind_cols(
    basic_info, variable = predictor_names,
    cor_test = NA_real_, auc_test = NA_real_)

  if (methods::.hasSlot(species_model, "varImportance")) {

    variable_importance <- species_model@varImportance

    if ("test.indep" %in% names(variable_importance)) {
      variable_importance <- variable_importance$test.indep

      if (methods::.hasSlot(variable_importance, "varImportance")) {
        variable_importance <- variable_importance@varImportance

        if (is.null(variable_importance) || length(variable_importance) == 0L) {
          variable_importance <- empty_var_imp
        } else {
          variable_importance <- tibble::tibble(variable_importance) %>%
            dplyr::rename(
              variable = variables, cor_test = corTest, auc_test = AUCtest) %>%
            dplyr::bind_cols(basic_info, .)
        }
      } else {
        variable_importance <- empty_var_imp
      }
    } else {
      variable_importance <- empty_var_imp
    }
  } else {
    variable_importance <- empty_var_imp
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Response curves -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Extract response curves
  r_curves <- sdm::getResponseCurve(x = model, id = 1L)

  empty_r_curves <- tibble::tibble(
    variable = predictor_names, x_value = NA_real_, prediction = NA_real_) %>%
    dplyr::bind_cols(basic_info, .)

  if (methods::.hasSlot(r_curves, "variables") &&
      methods::.hasSlot(r_curves, "response")) {
    if (length(r_curves@variables) == 0L || length(r_curves@response) == 0L) {
      r_curves <- empty_r_curves
    } else {
      r_curves <- purrr::map_dfr(
        .x = r_curves@variables,
        .f = ~{
          r_curves@response[[.x]] %>%
            tibble::tibble() %>%
            stats::setNames(c("x_value", "prediction")) %>%
            dplyr::mutate(variable = .x, .before = 1L) %>%
            dplyr::bind_cols(basic_info, .)
        })
    }
  } else {
    r_curves <- empty_r_curves
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # List of outputs
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  list(
    evaluation_training = evaluation_training,
    evaluation_testing = evaluation_testing,
    variable_importance = variable_importance, response_curves = r_curves)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# prepare_input_data ------
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' Prepare Input Data for Species Distribution Modelling
#'
#' This function prepares and validates all necessary input data for single
#' species distribution modelling of invasive alien plant species in Europe
#' fitted using `sdm` R package. This includes modelling (training and testing)
#' data, environmental predictors, and prediction datasets.
#' @param model_dir Character. Path to the directory containing model data and
#'   where outputs and results will be saved. Model data are prepared using the
#'   [mod_prepare_HPC()] and [mod_prepare_data()] functions.
#' @param cv_type Character. Cross-validation type. One of `CV_Dist` (default)
#'   or `CV_Large`. See [mod_CV_fit()] for more details.
#' @param selected_species Character vector or NULL. Names of species to include
#'   in modelling; must match names in model data. Default is `NULL` (all
#'   species).
#' @param excluded_species Character vector or NULL. Names of species to exclude
#'   from modelling; must match names in model data. Default is `NULL` (none
#'   excluded).
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
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default
#'   is 8.
#'
#' @details The function performs the following steps:
#' - Validates all input arguments and required files/directories.
#' - Loads model data prepared during the jSDMs workflow and checks for
#'   required components.
#' - Handles selection and exclusion of species for modelling.
#' - Prepares modelling data, including extraction of linear and quadratic
#'   predictors.
#' - Processes and clamps spatial predictors ("EffortsLog", "RiversLog",
#'   "HabLog", "RoadRailLog") as needed.
#' - Loads and filters CHELSA climate data for current and future prediction
#'   options.
#' - Merges all predictors and prepares prediction datasets for each climate
#'   scenario and period.
#' - Saves processed species modelling data, prediction data, and prediction
#'   options to disk.
#'
#'   If processed data already exists and is valid, the function returns paths
#'   to these files without reprocessing.
#'
#'   The function is not expected to be called directly by users, but is used
#'   internally within the `IASDT.R` workflow, particularly the
#'   [fit_sdm_models()] function.
#'
#' @return (Invisibly) A named list with the following elements:
#' - `path_species_data`: Path to the saved species modelling data file.
#' - `path_prediction_data`: Path to the saved prediction data file.
#' - `path_prediction_options`: Path to the saved prediction options file.
#'
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

prepare_input_data <- function(
    model_dir = NULL, cv_type = "CV_Dist", selected_species = NULL,
    excluded_species = NULL, env_file = ".env", hab_abb = NULL,
    clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    climate_models = "all", climate_scenarios = "all",
    climate_periods = "all", n_cores = 8) {

  Name <- TimePeriod <- ClimateScenario <- ClimateModel <- pred_df <- cv <-
    CellCode <- FilePath <- path_rail <- path_roads <- path_clc <- path_bias <-
    path_rivers <- path_chelsa <- pred_data <- quadratic <- CellNum <-
    variable <- climate_name <- species_name <- valid_species <- data_path <-
    path_wetness <- path_soil <- NULL

  # # ..................................................................... ###

  # Check input arguments ----
  all_args <- ls(envir = environment())
  all_args <- purrr::map(
    all_args,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(all_args)

  ecokit::check_args(
    args_all = all_args, args_type = "character",
    args_to_check = c(
      "model_dir", "cv_type", "env_file", "hab_abb",
      "fix_efforts", "fix_rivers"))
  ecokit::check_args(
    args_all = all_args, args_type = "logical", args_to_check = "clamp_pred")
  ecokit::check_args(
    args_all = all_args, args_type = "numeric", args_to_check = "n_cores")
  rm(all_args, envir = environment())

  # |||||||||||||||||||||||||||||||||||||||||||

  ## n_cores ----
  n_cores <- .validate_n_cores(n_cores)

  # |||||||||||||||||||||||||||||||||||||||||||

  ## future climate options ------

  valid_models <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  if (!is.null(climate_models)) {
    if (climate_models == "all") {
      climate_models <- valid_models
    } else if (!all(climate_models %in% valid_models)) {
      ecokit::stop_ctx(
        "Provided climate change models do not match valid values",
        climate_models = climate_models, valid_models = valid_models)
    }
  }

  valid_scenarios <- c("ssp126", "ssp370", "ssp585")
  if (!is.null(climate_scenarios)) {
    if (climate_scenarios == "all") {
      climate_scenarios <- valid_scenarios
    } else if (!all(climate_scenarios %in% valid_scenarios)) {
      ecokit::stop_ctx(
        "Provided climate change scenarios do not match valid values",
        climate_scenarios = climate_scenarios,
        valid_scenarios = valid_scenarios)
    }
  }

  valid_periods <- c("2011-2040", "2041-2070", "2071-2100")
  if (!is.null(climate_periods)) {
    if (climate_periods == "all") {
      climate_periods <- valid_periods
    } else if (!all(climate_periods %in% valid_periods)) {
      ecokit::stop_ctx(
        "Provided future time periods do not match valid values",
        climate_periods = climate_periods, valid_periods = valid_periods)
    }
  }

  ## CV type ------
  valid_cv <- c("CV_Dist", "CV_Large")
  if (!cv_type %in% valid_cv) {
    ecokit::stop_ctx("Invalid CV type", cv_type = cv_type, valid_cv = valid_cv)
  }

  ## Model dir ------
  if (!inherits(model_dir, "character") ||
      length(model_dir) != 1L || !nzchar(model_dir)) {
    ecokit::stop_ctx(
      "model_dir must be a character string",
      model_dir = model_dir, class_model_dir = class(model_dir))
  }
  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory does not exist", model_dir = model_dir)
  }

  # path to save all models
  sdm_model_dir <- fs::path(model_dir, "sdm_models")
  species_data_dir <- fs::path(sdm_model_dir, "species_data")
  fs::dir_create(species_data_dir)

  # species modelling data
  path_species_data <- fs::path(species_data_dir, "species_modelling_data.qs2")
  # prediction datasets
  path_prediction_data <- fs::path(species_data_dir, "prediction_data.qs2")
  path_prediction_options <- fs::path(
    species_data_dir, "prediction_data_options.qs2")

  ## Check if data already exist and valid -----
  model_data_exist <- all(
    ecokit::check_data(path_species_data, warning = FALSE),
    ecokit::check_data(path_prediction_data, warning = FALSE),
    ecokit::check_data(path_prediction_options, warning = FALSE))
  if (model_data_exist) {
    outputs <- list(
      path_species_data = path_species_data,
      path_prediction_data = path_prediction_data,
      path_prediction_options = path_prediction_options)
    return(invisible(outputs))
  }

  ## env file ------

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_rail", "DP_R_Railways_processed", TRUE, FALSE,
    "path_roads", "DP_R_Roads_processed", TRUE, FALSE,
    "path_clc", "DP_R_CLC_processed", TRUE, FALSE,
    "path_bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "path_rivers", "DP_R_Rivers_processed", TRUE, FALSE,
    "path_soil", "DP_R_soil_density", TRUE, FALSE,
    "path_wetness", "DP_R_wetness_processed", TRUE, FALSE,
    "path_chelsa", "DP_R_CHELSA_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  ## Habitat types -----
  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  ## clamp_pred / fix_efforts / fix_rivers -------
  if (clamp_pred) {

    # Check if the `fix_efforts` value is valid
    if (is.null(fix_efforts)) {
      ecokit::stop_ctx(
        "`fix_efforts` can not be `NULL` when Clamping is implemented",
        fix_efforts = fix_efforts, include_backtrace = TRUE)
    }
    if (length(fix_efforts) != 1L) {
      ecokit::stop_ctx(
        "`fix_efforts` must be a single value of length 1.",
        fix_efforts = fix_efforts, include_backtrace = TRUE)
    }

    # Check if the `fix_rivers` value is valid
    if (is.null(fix_rivers)) {
      ecokit::stop_ctx(
        "`fix_rivers` can not be `NULL` when Clamping is implemented",
        fix_rivers = fix_rivers, include_backtrace = TRUE)
    }
    if (length(fix_rivers) != 1L) {
      # Check if fix_rivers is a single value of length 1
      ecokit::stop_ctx(
        "`fix_rivers` must be a single value of length 1.",
        fix_rivers = fix_rivers, include_backtrace = TRUE)
    }
  }

  # # ..................................................................... ###

  # Loading model data ----

  modelling_data <- fs::path(model_dir, "ModDT.RData")
  data_CV <- fs::path(model_dir, "CV_data.RData")
  data_subset <- fs::path(model_dir, "ModDT_subset.RData")
  data_training <- fs::path(model_dir, "ModDT_training.RData")
  data_testing <- fs::path(model_dir, "ModDT_testing.RData")

  if (!ecokit::check_data(modelling_data, warning = FALSE)) {
    ecokit::stop_ctx(
      "Model data at the full extent file does not exist or is invalid",
      modelling_data = modelling_data, include_backtrace = TRUE)
  }
  if (!ecokit::check_data(data_subset, warning = FALSE)) {
    ecokit::stop_ctx(
      "Model data at the subset extent file does not exist or is invalid",
      data_subset = data_subset, include_backtrace = TRUE)
  }
  if (!ecokit::check_data(data_CV, warning = FALSE)) {
    ecokit::stop_ctx(
      "Cross-validation file does not exist or is invalid",
      data_CV = data_CV, include_backtrace = TRUE)
  }

  data_subset <- ecokit::load_as(data_subset)
  predictor_names <- names(data_subset$DT_x)
  model_form <- data_subset$Form_x
  train_test_exist <- ecokit::check_data(data_training, warning = FALSE) &&
    ecokit::check_data(data_testing, warning = FALSE)

  if (train_test_exist) {
    train_test_data <- dplyr::bind_rows(
      ecokit::load_as(data_training), ecokit::load_as(data_testing))
    model_cell_numbers <- dplyr::pull(train_test_data, "CellNum")
    model_cv_folds <- dplyr::select(
      train_test_data, CellNum, CellCode, cv_fold = tidyselect::all_of(cv_type))
    rm(train_test_data, envir = environment())
  } else {
    model_cell_numbers <- dplyr::pull(data_subset$DT_All, "CellNum")
    model_cv_folds <- dplyr::select(
      data_subset$DT_CV,
      CellNum, CellCode, cv_fold = tidyselect::all_of(cv_type))
  }

  modelling_data <- ecokit::load_as(modelling_data) %>%
    dplyr::filter(CellNum %in% model_cell_numbers) %>%
    dplyr::left_join(model_cv_folds, by = c("CellNum", "CellCode"))

  n_cv_folds <- length(unique(modelling_data$cv_fold))
  if (n_cv_folds < 2L) {
    ecokit::stop_ctx(
      "Not enough CV folds found in model data", length_cv_folds = n_cv_folds)
  }

  # `clamp_pred` can not be TRUE when `EffortsLog` is not used as predictor
  if (clamp_pred && !("EffortsLog" %in% predictor_names)) {
    ecokit::stop_ctx(
      "`clamp_pred` can not be used when `EffortsLog` is not used as predictor",
      clamp_pred = clamp_pred, names_data = predictor_names,
      include_backtrace = TRUE)
  }

  species_names <- names(data_subset$DT_y)
  chelsa_pattern <- paste0(
    "^", IASDT.R::CHELSA_variables$Variable, collapse = "|")
  other_variables <- stringr::str_subset(
    predictor_names, chelsa_pattern, negate = TRUE)
  bio_variables <- stringr::str_subset(
    predictor_names, chelsa_pattern, negate = FALSE)

  rm(data_subset, model_cv_folds, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Selecting species -----

  if (!is.null(selected_species)) {
    if (!inherits(selected_species, "character") ||
        length(selected_species) < 1L ||
        length(selected_species) > length(species_names)) {
      ecokit::stop_ctx(
        "selected_species must be a character vector of length >= 1",
        selected_species = selected_species,
        length_selected_species = length(selected_species),
        class_selected_species = class(selected_species))
    }
    if (!all(selected_species %in% species_names)) {
      invalid_species <- setdiff(selected_species, species_names)
      ecokit::stop_ctx(
        "Some or all of the selected species not found in model data",
        selected_species = selected_species, species_names = species_names,
        invalid_species = invalid_species)
    }
    species_names <- selected_species
  }

  if (length(species_names) == 0L) {
    ecokit::stop_ctx(
      "No species left after exclusions.", excluded_species = excluded_species)
  }

  if (!is.null(excluded_species)) {
    if (!inherits(excluded_species, "character") ||
        length(excluded_species) < 1L ||
        length(excluded_species) > length(species_names)) {
      ecokit::stop_ctx(
        "excluded_species must be a character vector of length >= 1",
        excluded_species = excluded_species,
        length_excluded_species = length(excluded_species),
        class_excluded_species = class(excluded_species))
    }
    if (!all(excluded_species %in% species_names)) {
      invalid_species <- setdiff(excluded_species, species_names)
      ecokit::stop_ctx(
        "Some or all of the selected species not found in model data",
        excluded_species = excluded_species, species_names = species_names,
        invalid_species = invalid_species)
    }
    species_names <- setdiff(species_names, excluded_species)
  }

  # # ..................................................................... ###

  # Prepare species data -----

  # extract list of linear and quadratic terms from the model formula
  # Use terms() to robustly extract variable names from the formula
  term_labels <- attr(stats::terms(model_form), "term.labels")
  formula_vars <- purrr::map(
    term_labels,
    .f = ~{
      # Check for quadratic terms using poly() or I()
      if (grepl("poly\\(|I\\(", .x)) {
        # Extract variable name inside poly() or I()
        var_name <- sub(".*\\(([^,\\)]+).*", "\\1", .x)
        tibble::tibble(variable = var_name, quadratic = TRUE)
      } else {
        tibble::tibble(variable = .x, quadratic = FALSE)
      }
    }) %>%
    dplyr::bind_rows()

  # list of linear and quadratic predictors
  l_preds <- dplyr::filter(formula_vars, !quadratic) %>%
    dplyr::pull(variable)
  q_preds <- dplyr::filter(formula_vars, quadratic) %>%
    dplyr::pull(variable)

  if (length(q_preds) > 0L) {
    modelling_data <- modelling_data %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(q_preds), .fns = ~ I(.x^2L), .names = "{.col}_sq"))
  }

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, length(species_names)), show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  pkg_to_load <- c(
    "dplyr", "ecokit", "fs", "qs2", "stringr", "purrr",
    "sdm", "tibble", "tidyselect", "rlang", "tidyr")
  future_globals <- c(
    "species_data_dir", "q_preds", "l_preds", "predictor_names",
    "modelling_data", "n_cv_folds", "model_form")

  species_modelling_data2 <- ecokit::quietly(
    future.apply::future_lapply(
      X = species_names,
      FUN = function(species_name) {

        method_is_glm <- NULL

        species_data_file <- fs::path(
          species_data_dir, paste0("data_", species_name, ".qs2"))

        if (ecokit::check_data(species_data_file, warning = FALSE)) {
          return(
            tibble::tibble(
              valid_species = TRUE, species_data = species_data_file))
        }

        species_modelling_data <- tidyr::expand_grid(
          cv = seq_len(n_cv_folds), method_is_glm = c(TRUE, FALSE)) %>%
          dplyr::mutate(
            model_form = purrr::map2(
              .x = cv, .y = method_is_glm,
              .f = ~ {

                valid_species <- TRUE

                # model formula
                if (.y) {
                  if (length(q_preds) == 0L) {
                    # No quadratic terms, use linear predictors only
                    model_formula <- paste(
                      species_name, " ~ ", paste(l_preds, collapse = " + "))
                    predictor_names_local <- predictor_names
                  } else {
                    # Add quadratic terms as columns to the modelling data
                    predictor_names_local <- c(
                      predictor_names, paste0(q_preds, "_sq"))
                    # Update model formula
                    model_formula <- c(
                      l_preds, q_preds, paste0(q_preds, "_sq")) %>%
                      paste(collapse = " + ") %>%
                      paste(species_name, " ~ ", .)
                  }
                } else {
                  predictor_names_local <- predictor_names
                  str_rem_1 <- "stats::poly\\(|, degree = 2, raw = TRUE\\)"
                  model_formula <- deparse1(model_form) %>%
                    stringr::str_remove_all(str_rem_1) %>%
                    paste(species_name, .)
                }

                if (length(model_formula) > 1) {
                  model_formula <- paste(model_formula, collapse = " ")
                }

                model_formula <- stats::as.formula(
                  model_formula, env = baseenv())

                # Training and testing data
                training_data <- modelling_data %>%
                  dplyr::filter(cv_fold != .x) %>%
                  dplyr::select(
                    tidyselect::all_of(c(species_name, predictor_names_local)))
                testing_data <- modelling_data %>%
                  dplyr::filter(cv_fold == .x) %>%
                  dplyr::select(
                    tidyselect::all_of(c(species_name, predictor_names_local)))

                valid_training <- training_data %>%
                  dplyr::distinct(.data[[species_name]]) %>%
                  unlist() %>%
                  sort() %>%
                  as.integer() %>%
                  identical(c(0L, 1L))
                valid_testing <- testing_data %>%
                  dplyr::distinct(.data[[species_name]]) %>%
                  unlist() %>%
                  sort() %>%
                  as.integer() %>%
                  identical(c(0L, 1L))

                # sdm data
                if (isFALSE(valid_training) || isFALSE(valid_testing)) {
                  sdm_data <- NULL
                  valid_species <- FALSE
                } else {
                  sdm_data <- sdm::sdmData(
                    formula = model_formula,
                    train = as.data.frame(training_data),
                    test = as.data.frame(testing_data))
                }

                tibble::tibble(
                  model_formula = list(model_formula),
                  predictor_names = list(predictor_names_local),
                  sdm_data = list(sdm_data), valid_species = valid_species)
              })) %>%
          tidyr::unnest("model_form") %>%
          dplyr::mutate(species_name = species_name, .before = 1)

        ecokit::save_as(
          object = species_modelling_data, out_path = species_data_file)

        tibble::tibble(
          valid_species = all(species_modelling_data$valid_species),
          data_path = species_data_file)
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_load, future.globals = future_globals))

  ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
  future::plan("sequential", gc = TRUE)
  invisible(gc())

  species_modelling_data <- dplyr::tibble(
    species_name = species_names, species_data = species_modelling_data2) %>%
    tidyr::unnest("species_data")
  rm(species_modelling_data2, envir = environment())
  invisible(gc())


  ## Check invalid species ------
  invalid_species <- dplyr::filter(species_modelling_data, !valid_species)

  if (nrow(invalid_species) > 0) {

    excluded_species <- unique(invalid_species$species_name)
    excluded_species_n <- length(excluded_species)

    paste0(
      "\n!! There are ", excluded_species_n,
      " invalid species that will be excluded in all model types !!") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE, ... = "\n")

    invalid_species %>%
      dplyr::mutate(
        message = purrr::map2_chr(
          .x = species_name, .y = data_path,
          .f = ~ {
            ecokit::load_as(.y) %>%
              dplyr::distinct(cv, valid_species) %>%
              dplyr::filter(!valid_species) %>%
              dplyr::pull(cv) %>%
              paste(collapse = "; ") %>%
              paste0(crayon::bold(.x), " (cv: ", ., ")")
          })) %>%
      dplyr::pull(message) %>%
      paste(collapse = " + ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1, ... = "\n")

    species_modelling_data <- species_modelling_data %>%
      dplyr::filter(!species_name %in% excluded_species) %>%
      dplyr::select(-tidyselect::all_of(valid_species))

  } else {
    excluded_species <- NA_character_
  }

  ecokit::save_as(object = species_modelling_data, out_path = path_species_data)

  # # ..................................................................... ###

  # Predictor data ------

  ## CHELSA data -----

  path_chelsa <- fs::path(path_chelsa, "CHELSA_Processed_DT.RData")
  if (!fs::file_exists(path_chelsa)) {
    ecokit::stop_ctx(
      "Processed CHLESA data can not be found", path_chelsa = path_chelsa,
      include_backtrace = TRUE)
  }

  prediction_options <- ecokit::load_as(path_chelsa) %>%
    dplyr::select(-"File_List") %>%
    dplyr::filter(
      ClimateModel %in% c("Current", climate_models),
      ClimateScenario %in% c("Current", climate_scenarios),
      TimePeriod  %in% c("1981-2010", climate_periods)) %>%
    dplyr::mutate(
      Name = paste0(TimePeriod, "_", ClimateScenario, "_", ClimateModel),
      Name = stringr::str_replace(Name, "1981-2010_Current_Current", "Current"),
      Name = stringr::str_replace_all(Name, "-", "_"),
      clamp = clamp_pred)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Road and railway intensity -----

  static_predictors <- list()

  if ("RoadRailLog" %in% other_variables) {

    r_railways <- fs::path(path_rail, "Railways_Length.RData")
    if (!fs::file_exists(r_railways)) {
      ecokit::stop_ctx(
        "Railways data does not exist", r_railways = r_railways,
        include_backtrace = TRUE)
    }
    r_railways <- ecokit::load_as(r_railways, unwrap_r = TRUE) %>%
      magrittr::extract2("rail")

    r_roads <- fs::path(path_roads, "Road_Length.RData")
    if (!fs::file_exists(r_roads)) {
      ecokit::stop_ctx(
        "Roads data does not exist", r_roads = r_roads,
        include_backtrace = TRUE)
    }
    r_roads <- ecokit::load_as(r_roads, unwrap_r = TRUE) %>%
      magrittr::extract2("All")

    # Calculating the sum of road and railway intensity
    # add 1 (older versions 0.1) to get log for 0 values
    # [only for # rivers/roads/efforts, not hab/rivers]
    r_road_rail <- (r_roads + r_railways) %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("RoadRailLog")

    static_predictors <- c(static_predictors, r_road_rail)
    rm(r_road_rail, r_roads, r_railways, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Habitat information ----

  # Check if habitat information is used as predictor
  hab_predictor <- "HabLog" %in% other_variables

  if (hab_predictor) {

    r_hab <- fs::path(
      path_clc, "Summary_RData", "PercCov_SynHab_Crop.RData")
    if (!fs::file_exists(r_hab)) {
      ecokit::stop_ctx(
        "Habitat data does not exist", r_hab = r_hab, include_backtrace = TRUE)
    }

    r_hab <- ecokit::load_as(r_hab, unwrap_r = TRUE) %>%
      magrittr::extract2(paste0("SynHab_", hab_abb))

    # Models are trained and predictions are made only at grid cells with > 0 %
    # coverage. Mask layer to exclude grid cells with zero % coverage from
    # predictions.
    r_hab_mask <- terra::classify(r_hab, cbind(0L, NA), others = 1L)

    r_hab <- log10(r_hab + 0.1) %>%
      stats::setNames("HabLog")

    static_predictors <- c(static_predictors, r_hab)
    rm(r_hab, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Sampling efforts -----

  if ("EffortsLog" %in% other_variables) {

    r_efforts <- fs::path(path_bias, "Efforts_SummaryR.RData")
    if (!fs::file_exists(r_efforts)) {
      ecokit::stop_ctx(
        "Sampling efforts data does not exist", r_efforts = r_efforts,
        include_backtrace = TRUE)
    }

    # add 1 (older versions 0.1) to get log for 0 values
    # [only for rivers/roads/efforts, not hab/rivers]
    r_efforts <- ecokit::load_as(r_efforts, unwrap_r = TRUE) %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("EffortsLog")


    if (clamp_pred) {

      # Check fix_efforts value
      if (is.numeric(fix_efforts)) {

        # If `fix_efforts` is numeric value, check if it is within the range of
        # the observed efforts
        efforts_range <- terra::global(r_efforts, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_efforts, efforts_range[1L], efforts_range[2L]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_efforts` value is out of the range of observed efforts",
            fix_efforts = fix_efforts, efforts_range = round(efforts_range, 2L),
            include_backtrace = TRUE)
        }

        # Fix value
        efforts_val <- fix_efforts
        fix_efforts_lc <- "fixed"

      } else {

        # If `fix_efforts` is character, check if it is one of the valid values:
        # identity, median, mean, max, and q90
        fix_efforts_lc <- stringr::str_to_lower(fix_efforts)

        if (!(fix_efforts_lc %in%
              c("identity", "median", "mean", "max", "q90"))) {
          ecokit::stop_ctx(
            paste0(
              "`fix_efforts` has to be either NULL, single numeric ",
              "value, or one of the following: 'identity', 'median', ",
              "'mean', 'max', or 'q90'."),
            fix_efforts = fix_efforts, include_backtrace = TRUE)
        }
      }

      # Fix value
      if (fix_efforts_lc != "fixed") {
        efforts_val <- dplyr::case_when(
          # Do not fix if `fix_efforts` is "identity"
          fix_efforts_lc == "identity" ~ NA_real_,
          # Fix at 90% quantile
          fix_efforts_lc == "q90" ~ {
            terra::global(
              r_efforts,
              fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at median value
          fix_efforts_lc == "median" ~ {
            terra::global(
              r_efforts, fun = function(x) median(x, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at mean value
          fix_efforts_lc == "mean" ~ {
            terra::global(r_efforts, fun = mean, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at max value
          fix_efforts_lc == "max" ~ {
            terra::global(r_efforts, fun = max, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },
          .default = NA_real_)
      }

      if (is.na(efforts_val)) {
        r_efforts_clamp <- stats::setNames(r_efforts, "EffortsLog_clamp")
      } else {
        # Set a minimum value for efforts variable to `efforts_val`. Using
        # upper = Inf keeps efforts values > efforts_val as they are.
        r_efforts_clamp <- terra::clamp(
          x = r_efforts, lower = efforts_val, upper = Inf) %>%
          stats::setNames("EffortsLog_clamp")
      }

      static_predictors <- c(static_predictors, r_efforts, r_efforts_clamp)
      rm(r_efforts, r_efforts_clamp, envir = environment())

    } else {

      # Do not fix at single value
      static_predictors <- c(static_predictors, r_efforts)
      rm(r_efforts, envir = environment())

    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## River length ----

  if ("RiversLog" %in% other_variables) {

    r_rivers <- fs::path(path_rivers, "River_Lengths.RData")
    if (!fs::file_exists(r_rivers)) {
      ecokit::stop_ctx(
        "River length data does not exist", r_rivers = r_rivers,
        include_backtrace = TRUE)
    }

    r_rivers <- ecokit::load_as(r_rivers, unwrap_r = TRUE) %>%
      magrittr::extract2("STRAHLER_5") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("RiversLog")

    if (clamp_pred) {

      # Check fix_rivers value
      if (is.numeric(fix_rivers)) {

        # If `fix_rivers` is numeric value, check if it is within the range of
        # the observed river lengths
        rivers_range <- terra::global(r_rivers, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_rivers, rivers_range[1L], rivers_range[2L]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_rivers` value is out of the range of observed river length",
            fix_rivers = fix_rivers, rivers_range = round(rivers_range, 2L),
            include_backtrace = TRUE)
        }

        # Fix value
        rivers_value <- fix_rivers

      } else {

        # If `fix_rivers` is character, check if it is one of the valid values:
        # identity, median, mean, max, and q90
        fix_rivers_lc <- stringr::str_to_lower(fix_rivers)

        valid_river_fix_values <- c("identity", "median", "mean", "max", "q90")
        if (!(fix_rivers_lc %in% valid_river_fix_values)) {
          ecokit::stop_ctx(
            paste0(
              "`fix_rivers` has to be either NULL, single numeric ",
              "value, or one of the following: 'identity', 'median', ",
              "'mean', 'max', or 'q90'."),
            fix_rivers = fix_rivers, include_backtrace = TRUE)
        }

        # Fix value
        rivers_value <- dplyr::case_when(
          # Do not fix if `fix_rivers` is "identity"
          fix_rivers_lc == "identity" ~ NA_real_,
          # Fix at 90% quantile
          fix_rivers_lc == "q90" ~ {
            terra::global(
              r_rivers,
              fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at median value
          fix_rivers_lc == "median" ~ {
            terra::global(
              r_rivers, fun = function(x) median(x, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at mean value
          fix_rivers_lc == "mean" ~ {
            terra::global(r_rivers, fun = mean, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },
          # Fix at max value
          fix_rivers_lc == "max" ~ {
            terra::global(r_rivers, fun = max, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },
          .default = NA_real_)
      }

      if (is.na(rivers_value)) {
        r_rivers_clamp <- stats::setNames(r_rivers, "RiversLog_clamp")
      } else {
        # Set a minimum value for river length variable to `rivers_value`. Using
        # upper = Inf keeps  river length values > rivers_value as they are.
        r_rivers_clamp <- terra::clamp(
          x = r_rivers, lower = rivers_value, upper = Inf) %>%
          stats::setNames("RiversLog_clamp")
      }

      static_predictors <- c(static_predictors, r_rivers, r_rivers_clamp)
      rm(r_rivers, r_rivers_clamp, rivers_value, envir = environment())

    } else {

      # Do not fix at single value
      static_predictors <- c(static_predictors, r_rivers)
      rm(r_rivers, envir = environment())

    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Soil bulk density -----

  if ("soil" %in% other_variables) {
    r_soil <- fs::path(path_soil, "soil_density.RData")
    if (!fs::file_exists(r_soil)) {
      ecokit::stop_ctx(
        "Soil bulk density data does not exist", r_soil = r_soil,
        include_backtrace = TRUE)
    }
    r_soil <- ecokit::load_as(r_soil, unwrap_r = TRUE) %>%
      magrittr::extract2("bdod_5_15_mean") %>%
      stats::setNames("soil")

    static_predictors <- c(static_predictors, r_soil)
    rm(r_soil, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Topographic wetness index -----

  if ("wetness" %in% other_variables) {
    r_wetness <- fs::path(path_wetness, "wetness_index.RData")
    if (!fs::file_exists(r_wetness)) {
      ecokit::stop_ctx(
        "Topographic wetness index data does not exist", r_wetness = r_wetness,
        include_backtrace = TRUE)
    }
    r_wetness <- ecokit::load_as(r_wetness, unwrap_r = TRUE) %>%
      stats::setNames("wetness")

    static_predictors <- c(static_predictors, r_wetness)
    rm(r_wetness, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  static_predictors <- terra::rast(static_predictors)

  # If Habitat predictor is used and r_hab_mask exists, grid cells with zero %
  # coverage are excluded from predictions
  if (hab_predictor && exists("r_hab_mask", inherits = FALSE)) {
    static_predictors <- terra::mask(static_predictors, r_hab_mask)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  if (clamp_pred) {

    if ("EffortsLog" %in% names(static_predictors)) {
      # use clamped Effort values
      static_predictors$EffortsLog <- static_predictors$EffortsLog_clamp
      static_predictors$EffortsLog_clamp <- NULL
    }
    if (all(c("RiversLog", "RiversLog_clamp") %in% names(static_predictors))) {
      # use clamped rivers values
      static_predictors$RiversLog <- static_predictors$RiversLog_clamp
      static_predictors$RiversLog_clamp <- NULL
    }

  } else {

    # Remove clamped layers
    if ("EffortsLog_clamp" %in% names(static_predictors)) {
      static_predictors <- terra::subset(
        x = static_predictors, subset = "EffortsLog_clamp", negate = TRUE)
    }
    if ("RiversLog_clamp" %in% names(static_predictors)) {
      static_predictors <- terra::subset(
        x = static_predictors, subset = "RiversLog_clamp", negate = TRUE)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Extracting data at training and new sites ------
  prediction_data <- prediction_options %>%
    dplyr::mutate(
      pred_df = purrr::map(
        FilePath,
        ~ {
          pred_df0 <- ecokit::load_as(.x, unwrap_r = TRUE) %>%
            terra::subset(bio_variables) %>%
            c(static_predictors) %>%
            terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
            tibble::tibble()

          if (length(q_preds) > 0L) {
            pred_df0 <- pred_df0 %>%
              dplyr::mutate(
                dplyr::across(
                  tidyselect::all_of(q_preds),
                  .fns = ~ .x^2L, .names = "{.col}_sq"))
          }
          pred_df0
        })) %>%
    dplyr::select(-FilePath) %>%
    tidyr::unnest(pred_df) %>%
    dplyr::select(-tidyselect::any_of(c("clamp", "Processed_Name"))) %>%
    dplyr::rename(
      time_period = TimePeriod, climate_model = ClimateModel,
      climate_scenario = ClimateScenario, climate_name = Name) %>%
    tidyr::nest(
      pred_data = -c(
        "time_period", "climate_model",
        "climate_scenario", "climate_name")) %>%
    dplyr::mutate(climate_name = stringr::str_to_lower(climate_name))

  ecokit::save_as(object = prediction_data, out_path = path_prediction_data)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  prediction_data_options <- dplyr::select(prediction_data, -pred_data)
  ecokit::save_as(
    object = prediction_data_options, out_path = path_prediction_options)

  outputs <- list(
    path_species_data = path_species_data, excluded_species = excluded_species,
    path_prediction_data = path_prediction_data,
    path_prediction_options = path_prediction_options)

  rm(bio_variables, l_preds, q_preds, static_predictors, envir = environment())
  return(invisible(outputs))

}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# sdm_model_settings ------
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' @title Species Distribution Model (SDM) Model Settings
#'
#' @description Defines a list of default settings for various species
#'   distribution modelling algorithms. Each element in the list corresponds to
#'   a specific modelling method and contains its associated control parameters.
#'
#' @return A named list containing model-specific settings for use in SDM
#'   workflows.
#' @details The function is not expected to be called directly by users, but is
#'   used internally within the `IASDT.R` workflow, particularly the
#'   [fit_sdm_models()] function.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

sdm_model_settings <- function() {
  list(
    glm = list(control = list(maxit = 50L)),
    glmpoly = list(degree = 2L),
    gam = list(method = "REML", select = TRUE, gamma = 1.2),
    glmnet = list(maxit = 200000L),
    mars = list(pmethod = "backward"),
    gbm = list(n.trees = 3000L, interaction.depth = 2L),
    rf = list(ntree = 1000L, nodesize = 5L),
    ranger = list(
      num.trees = 1000L, importance = "impurity", min.node.size = 5L),
    cart = list(),
    rpart = list(),
    maxent = list(
      removeDuplicates = FALSE,
      args = c(
        "maximumiterations=2000",
        # "betamultiplier=1",
        "convergencethreshold=0.000005",
        "noautofeature", "hinge", "linear",
        "noproduct", "noquadratic", "nothreshold",
        "writeplotdata", "-J", "-P",
        "noremoveduplicates", "noaddsamplestobackground",
        "doclamp", "nowriteclampgrid", "outputformat=cloglog")),
    mlp = list(maxit = 250L),
    # rbf = list(maxit = 200L),
    svm = list(),
    svm2 = list(),
    mda = list(),
    fda = list())
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# reduce_sdm_formulas ------
# # ========================================================================= #

#' Reduce SDM Model Formulas to Base Environment To Reduce File Size
#'
#' This function iterates through the models stored in a `sdm` object, and
#' ensures that all formulas and their associated environments are set to the
#' base environment. This helps to reduce the fitted model's object size, see
#' [here](https://github.com/babaknaimi/sdm/issues/43).
#'
#' @param obj An `sdmModels` object containing fitted SDM models.
#' @return A copy of the input object with all model formulas and terms
#'   environments set to the base environment.
#' @author Ahmed El-Gabbas
#' @details The function is not expected to be called directly by users, but is
#'   used internally within the `IASDT.R` workflow, particularly the
#'   [fit_sdm_models()] function.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

reduce_sdm_formulas <- function(obj) {

  if (is.null(obj) || !inherits(obj, "sdmModels")) {
    ecokit::stop_ctx(
      "Input object must be of class 'sdmModels'.",
      input_object = obj, input_class = class(obj))
  }

  # Make a copy to avoid modifying the original
  model2 <- obj

  # Loop through each model using the structure of the SDM object
  for (sp_name in names(model2@models)) {
    for (method_name in names(model2@models[[sp_name]])) {
      for (run_idx in seq_along(model2@models[[sp_name]][[method_name]])) {

        # Check if there's a formula in the call
        model <- model2@models[[sp_name]][[method_name]][[run_idx]]

        if (
          !isS4(model@object) &&
          !is.null(model@object$call$formula) &&
          deparse1(model@object$call$formula) != ".f") {

          # Replace formula in the model call (if it exists and is not ".f")
          formula_text <- deparse1(model@object$call$formula)
          clean_formula <- stats::as.formula(formula_text, env = baseenv())
          model@object$call$formula <- clean_formula

          if (!is.null(model@object$formula)) {
            # Replace formula with clean version
            formula_text <- deparse1(model@object$formula)
            clean_formula <- stats::as.formula(formula_text, env = baseenv())
            model@object$formula <- clean_formula
          }

          # Fix terms environment if it exists
          if (!is.null(model@object$terms)) {
            attr(model@object$terms, ".Environment") <- baseenv()
          }
          # Update the model in the object
          model2@models[[sp_name]][[method_name]][[run_idx]] <- model
        }

      }
    }
  }

  # Force garbage collection
  invisible(gc())

  model2
}

# ............................................... ----

# # ========================================================================= #
# fit_predict_internal ------
# # ========================================================================= #

#' Fit and predict SDM models for a given species and cross-validation fold.
#'
#' This function fits species distribution models using the `sdm` package with
#' the specified method and data, saves the fitted model, extracts relevant
#' information, makes predictions on new data under current and selected future
#' climate scenarios, saves prediction rasters and data, and returns a summary
#' of results for the species and fold.
#'
#' @param line_id Integer. Index of the row in `model_data` for the species and
#'   cross-validation fold.
#' @param sdm_method Character. A single species distribution modelling
#'   algorithm to use for fitting models. Valid values: "glm", "glmpoly", "gam",
#'   "glmnet", "mars", "gbm", "rf", "ranger", "cart", "rpart", "maxent", "mlp",
#'   "rbf", "svm", "mda", and "fda". These correspond to selected methods
#'   supported by the `sdm` package. For details and supported options, see
#'   [sdm::getmethodNames()].
#' @param model_data Data frame with model formulas, SDM data, species names,
#'   and cross-validation folds.
#' @param model_settings List. Model-specific settings passed to the SDM fitting
#'   function.
#' @param input_data List. Contains paths and objects for prediction and input.
#' @param output_directory Character. Path where prediction outputs are stored.
#' @param path_grid_r Character. Path to the reference raster grid for spatial
#'   predictions.
#' @param copy_maxent_html Logical. Whether to copy the directory containing
#'   HTML results from Maxent to the modelling directory. Default is `TRUE`.
#'
#' @return A tibble summarizing the results for the species and cross-validation
#'   fold, including model paths, prediction raster paths, and status flags.
#'
#' @details
#'  - The function is not expected to be called directly by users, but is
#' used internally within the `IASDT.R` workflow, particularly the
#' [fit_sdm_models()] function.
#'  - The function checks if the fitted model and extracted data already exist;
#' if not, fits the model and extracts information.
#' - Makes predictions for each prediction dataset, saves rasters and data,
#' and checks their validity.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

fit_predict_internal <- function(
    line_id, sdm_method, model_data, model_settings,
    input_data, output_directory, path_grid_r, copy_maxent_html = TRUE) {

  pred_data <- climate_name <- method_is_glm <- cv <- NULL

  species_name <- model_data$species_name[[line_id]]
  cv_fold <- model_data$cv[[line_id]]
  model_DT <- ecokit::load_as(model_data$data_path[[line_id]]) %>%
    dplyr::filter(method_is_glm == (sdm_method == "glm"), cv == cv_fold)
  base_model_name <- paste0(sdm_method, "_", species_name, "_cv", cv_fold)

  if (nrow(model_DT) != 1) {
    ecokit::stop_ctx("Modelling data should be only one row")
  }

  output_path <- fs::path(output_directory, paste0(base_model_name, ".qs2"))
  if (ecokit::check_data(output_path, warning = FALSE)) {
    return(tibble::tibble(sdm_method = sdm_method, output_path = output_path))
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Model fit -----

  if (sdm_method == "maxent") {
    # Use a temp directory for Java preferences to avoid file lock issues
    prefs_dir <- file.path(tempdir(), paste0(".java_", Sys.getpid()))
    dir.create(prefs_dir, showWarnings = FALSE)
    Sys.setenv(     # nolint: undesirable_function_linter
      JAVA_TOOL_OPTIONS = paste0("-Djava.util.prefs.userRoot=", prefs_dir))
  }

  fitted_model <- ecokit::quietly(
    sdm::sdm(
      formula = model_DT$model_formula[[1]], data = model_DT$sdm_data[[1]],
      methods = sdm_method, modelSettings = model_settings))


  # copy model files from temp dir for maxent models
  if (sdm_method == "maxent") {

    maxent_html <- ecokit::normalize_path(
      fitted_model@models[[1L]][[1L]][[1L]]@object@html)

    if (fs::file_exists(maxent_html) && copy_maxent_html) {
      out_maxent_dir <- fs::path(
        output_directory, "maxent_html", base_model_name)
      if (fs::dir_exists(out_maxent_dir))  fs::dir_delete(out_maxent_dir)
      fs::dir_create(out_maxent_dir)
      fs::dir_copy(
        fs::path_dir(maxent_html), out_maxent_dir, overwrite = TRUE)

      # delete some not-needed files to save space
      out_maxent_dir %>%
        fs::path(c("presence", "absence", "species_explain.bat")) %>%
        fs::file_delete()

      # Overwrite the new HTML file path in the model object
      fitted_model@models[[1L]][[1L]][[1L]]@object@html <-
        fs::path(out_maxent_dir, "maxent.html")
    }

    # clean up temp dir containing maxent results
    fs::dir_delete(fs::path_dir(maxent_html))
  }

  # Reduce models objects
  fitted_model <- reduce_sdm_formulas(obj = fitted_model)

  invisible(gc())

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Extract info from fitted model object -----
  extracted_data <- ecokit::quietly(
    extract_sdm_info(model = fitted_model, cv_fold = cv_fold))

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Making predictions and extract prediction paths -----

  predictor_names <- fitted_model@setting@featureFrame@predictors
  if (is.null(predictor_names) || length(predictor_names) == 0L) {
    ecokit::stop_ctx("No predictor names found in fitted model.")
  }

  prediction_data <- ecokit::load_as(input_data$path_prediction_data) %>%
    dplyr::mutate(
      preds = purrr::map2(
        .x = pred_data,
        .y = climate_name,
        .f = ~ {
          pred2 <- unlist(
            ecokit::quietly(
              predict(
                object = fitted_model,
                newdata = dplyr::select(.x, tidyselect::all_of(predictor_names))
              )))

          prediction_r <- dplyr::mutate(.x, pred = unlist(pred2)) %>%
            dplyr::select(
              -tidyselect::all_of(c(predictor_names, "cell"))) %>%
            sf::st_as_sf(coords = c("x", "y"), crs = 3035L) %>%
            terra::rasterize(
              y = ecokit::load_as(path_grid_r, unwrap_r = TRUE),
              field = "pred", fun = "mean", na.rm = TRUE) %>%
            stats::setNames(paste0(base_model_name, "_", .y))

          prediction_okay <- prediction_r %>%
            terra::global(range, na.rm = TRUE) %>%
            unlist() %>%
            dplyr::between(-0.0000001, 1.00000001) %>%
            all()

          tibble::tibble(
            pred = list(terra::wrap(prediction_r)), pred_okay = prediction_okay)

        })) %>%
    tidyr::unnest("preds") %>%
    dplyr::select(-pred_data) %>%
    dplyr::mutate(
      species_name = species_name, sdm_method = sdm_method,
      cv_fold = cv_fold, .before = 1)

  # Merge outputs -------

  c(fitted_model = fitted_model, extracted_data,
    prediction_data = list(prediction_data)) %>%
    ecokit::save_as(out_path = output_path)

  tibble::tibble(sdm_method = sdm_method, output_path = output_path)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# summarize_predictions -------
# # ========================================================================= #

#' Summarize Model Predictions and Generate Raster Statistics
#'
#' This function processes model prediction summaries for a given line
#' (species/model combination), computes summary raster statistics (mean,
#' weighted mean, standard deviation, coefficient of variation) across
#' cross-validation folds, and ensures the existence and validity of
#' corresponding TIFF and RData files. If summary files do not exist or are
#' invalid, they are generated and saved.
#'
#' @param line_id Integer. Index (row number) in `model_pred_results` to process
#'   for summary statistics.
#' @param model_pred_results List. Model summary structure containing:
#'   - `summary_data`: List of data frames with prediction details.
#'   - `species_name`: Character vector of species names.
#'   - `evaluation_testing`: List of evaluation results, including test AUC.
#'
#' @return A data frame (tibble) combining original prediction info and summary
#'   statistics for each combination of time period, climate model, and
#'   scenario. Includes columns for:
#'   - `species_name`: Name of the species.
#'   - `time_period`, `climate_model`, `climate_scenario`, `climate_name`,
#'   `pred_dir`: Metadata columns.
#'   - `cv_fold`: Cross-validation fold or summary statistic ("mean",
#'   "weighted_mean", "sd", "cov").
#'   - `data_path`, `tif_path`: Paths to RData and TIFF files for
#'   predictions/statistics.
#'   - `data_okay`, `tif_okay`: Logical flags indicating file validity.
#'
#' @details
#' - The function is not expected to be called directly by users, but is
#' used internally within the `IASDT.R` workflow, particularly the
#' [fit_sdm_models()] function.
#' - For each group of predictions (by time period, climate model, scenario,
#' etc.), the function:
#'   - Checks validity of input TIFF files.
#'   - Computes and saves summary rasters (mean, weighted mean by AUC, standard
#' deviation, coefficient of variation) if not already present or invalid.
#'   - Returns a tidy data frame with paths and validity flags for all
#' prediction and summary files.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

summarize_predictions <- function(
    line_id, model_pred_results, output_directory) {

  preds <- output_path <- pred <- pred_summary <- time_period <-
    climate_model <- climate_scenario <- NULL

  pred_path <- fs::path(
    output_directory,
    paste0(
      unique(model_pred_results$sdm_method), "_",
      model_pred_results$species_name[[line_id]], "_summary_pred.qs2")
  )

  pred_out <- tibble::tibble(
    species_name = model_pred_results$species_name[[line_id]],
    summary_prediction_path = pred_path)

  if (ecokit::check_data(pred_path, warning = FALSE)) {
    return(pred_out)
  }

  mean_auc <- model_pred_results$auc_test[[line_id]]
  # If any value of testing AUC is NA, do not calculate weighted mean
  if (anyNA(mean_auc)) {
    calc_w_mean <- FALSE             # nolint: object_usage_linter
  } else {
    calc_w_mean <- TRUE              # nolint: object_usage_linter

    # avoid extreme cases when any of testing AUC is very small (= 0)
    n_zeros <- which(mean_auc == 0L)
    if (length(n_zeros) > 0L) {
      mean_auc[n_zeros] <- 0.001
    }
  }

  group_by_cols <- c(
    "time_period", "climate_model", "climate_scenario", "climate_name")

  pred_summ <- model_pred_results$pred_paths[[line_id]] %>%
    dplyr::mutate(
      preds = purrr::map(
        .x = output_path,
        .f = ~ {
          ecokit::load_as(.x) %>%
            magrittr::extract2("prediction_data")
        })) %>%
    tidyr::unnest("preds") %>%
    dplyr::select(-output_path) %>%
    dplyr::arrange(time_period, climate_model, climate_scenario) %>%
    dplyr::summarise(
      preds = list(pred), .by = tidyselect::all_of(group_by_cols)) %>%
    dplyr::mutate(
      pred_summary = purrr::map(
        .x = preds,
        .f = ~ {

          cv_maps <- terra::rast(purrr::map(.x, terra::unwrap))

          # |||||||||||||||||||||||||||||||||||||||||||

          # mean -------
          name_mean <- stringr::str_replace(
            names(cv_maps)[1], "cv[0-9]_", "mean_")
          pred_mean <- terra::app(cv_maps, mean, na.rm = TRUE) %>%
            stats::setNames(name_mean) %>%
            terra::wrap()

          # |||||||||||||||||||||||||||||||||||||||||||

          # weighted mean -----
          if (calc_w_mean) {
            name_w_mean <- stringr::str_replace(name_mean, "mean", "w_mean")
            pred_w_mean <- terra::weighted.mean(
              x = cv_maps, w = mean_auc, na.rm = TRUE) %>%
              stats::setNames(name_w_mean) %>%
              terra::wrap()
          } else {
            pred_w_mean <- list()
          }

          # |||||||||||||||||||||||||||||||||||||||||||

          # sd ------
          name_sd <- stringr::str_replace(name_mean, "mean", "sd")
          pred_sd <- terra::app(cv_maps, sd, na.rm = TRUE) %>%
            stats::setNames(name_sd) %>%
            terra::wrap()

          # |||||||||||||||||||||||||||||||||||||||||||

          # cov ------
          name_cov <- stringr::str_replace(name_mean, "mean", "cov")
          # Avoid division by zero
          pred_mean2 <- terra::clamp(terra::unwrap(pred_mean), lower = 1e-8)
          pred_cov <- (terra::unwrap(pred_sd) / pred_mean2) %>%
            stats::setNames(name_cov) %>%
            terra::wrap()

          # |||||||||||||||||||||||||||||||||||||||||||

          tibble::tibble(
            pred_mean = list(pred_mean), pred_w_mean = list(pred_w_mean),
            pred_sd = list(pred_sd), pred_cov = list(pred_cov))

        })
    ) %>%
    dplyr::select(-preds) %>%
    tidyr::unnest(pred_summary)

  ecokit::save_as(object = pred_summ, out_path = pred_path)

  pred_out
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# check_model_results -------
# # ========================================================================= #

#' Check Species Distribution Model Results for Issues
#'
#' This function examines a species distribution model (SDM) results object for
#' potential issues. It checks for missing or `NaN` values in:
#' - 1. Training evaluation metrics
#' - 2. Testing evaluation metrics
#' - 3. Variable importance data
#' - 4. Response curves
#' - 5. Prediction files (TIF files and data files)
#'
#' @param model_results A list object containing SDM results with the following
#'   components:
#'   - `evaluation_training`: List of tibbles with training evaluation metrics
#'   - `evaluation_testing`: List of tibbles with testing evaluation metrics
#'   - `variable_importance`: List of tibbles with variable importance values
#'   - `response_curves`: List of tibbles with response curve data
#'   - `prediction_info`: List of tibbles with prediction file information
#' @param n_cores Integer. Number of CPU cores for parallel processing.
#' @param future_max_size Numeric or character. Maximum allowed size for future
#'   package globals (see future documentation).
#'
#' @details
#' - The function outputs information about any identified issues to the
#' console, including which species are affected and what specific metrics or
#' data elements have issues.
#' - The function is not expected to be called directly by users, but is
#' used internally within the `IASDT.R` workflow, particularly the
#' [fit_sdm_models()] function.

#' @return Returns invisibly NULL. This function is used for its side effects
#'   (console output). The return value is always `invisible(NULL)`.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

check_model_results <- function(model_results, n_cores, future_max_size) {

  climate_name <- cv_fold <- x_value <- na_count <- prediction <- cor_test <-
    auc_test <- species_name <- . <- variable <- prediction_data <- times <-
    response_curves <- evaluation_training <- evaluation_testing <-
    variable_importance <- pred_okay <- output_path <- sdm_method <- NULL

  n_species <- length(unique(model_results$species_name))

  # Check if model_results is a data frame with at least one row
  if (!inherits(model_results, "data.frame") || nrow(model_results) == 0L) {
    ecokit::stop_ctx(
      "model_results must be a data frame with at least 1 row.",
      class_model_results = class(model_results))
  }

  issue_detected <- FALSE

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(model_results)),
      future_max_size = future_max_size, show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  all_model_results0 <- ecokit::quietly(
    future.apply::future_lapply(
      X = seq_len(nrow(model_results)),
      FUN = function(line_id) {

        model_results$output_path[[line_id]] %>%
          ecokit::load_as() %>%
          magrittr::extract(names(.) != "fitted_model") %>%
          lapply(list) %>%
          tibble::as_tibble()

      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = c("tibble", "sdm", "magrittr"),
      future.globals = "model_results"))

  ecokit::set_parallel(level = 1L, show_log = FALSE)
  future::plan("sequential", gc = TRUE)


  all_model_results <- model_results %>%
    dplyr::select(species_name, sdm_method, cv_fold, output_path) %>%
    dplyr::bind_cols(dplyr::bind_rows(all_model_results0))

  rm(all_model_results0)
  invisible(gc())

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation data - Training -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_eval_train <- dplyr::pull(all_model_results, evaluation_training) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(species_name) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = -c("sdm_method", "cv_fold"),
        .fns = ~ sum(is.na(.) | is.nan(.))),
      .groups = "drop") %>%
    dplyr::select(
      species_name, tidyselect::where(~is.numeric(.) && sum(.) != 0L)) %>%
    dplyr::filter(dplyr::if_any(.cols = -species_name, .fns =  ~ .x > 0L))

  if (ncol(check_eval_train) > 1L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_eval_train <- check_eval_train %>%
      tidyr::pivot_longer(
        -species_name,
        names_to = "variable", values_to = "na_count") %>%
      dplyr::filter(na_count > 0L)

    issues_eval_train_species <- unique(issues_eval_train$species_name)
    issues_eval_train_n_species <- length(issues_eval_train_species)

    paste0(
      "\n!! There are issues in training evaluation data for ",
      issues_eval_train_n_species, " / ", n_species, " species !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(issues_eval_train_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    dplyr::select(issues_eval_train, -species_name) %>%
      dplyr::mutate(
        variable = stringr::str_remove_all(
          variable, "_mss_train|_ess_train|_train"),
        variable = stringr::str_to_lower(variable)) %>%
      dplyr::pull(variable) %>%
      unique() %>%
      paste(collapse = "; ") %>%
      stringr::str_wrap(50L) %>%
      stringr::str_replace_all("\n", "\n  >>>  ") %>%
      paste0(
        crayon::bold("Affected metrics: "),
        "\n  >>>  ", ., collapse = "\n") %>%
      ecokit::cat_time(cat_timestamp = FALSE)
  }

  all_model_results <- dplyr::select(all_model_results, -evaluation_training)
  invisible(gc())

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation data - Testing -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_eval_test <- dplyr::pull(all_model_results, evaluation_testing) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(species_name) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = -c("sdm_method", "cv_fold"),
        .fns =  ~ sum(is.na(.) | is.nan(.))),
      .groups = "drop") %>%
    dplyr::select(
      species_name, tidyselect::where(~is.numeric(.) && sum(.) != 0L)) %>%
    dplyr::filter(dplyr::if_any(.cols = -species_name, .fns =  ~ .x > 0L))

  if (ncol(check_eval_test) > 1L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_eval_test <- check_eval_test %>%
      tidyr::pivot_longer(
        -species_name, names_to = "variable", values_to = "na_count") %>%
      dplyr::filter(na_count > 0L)

    issues_eval_test_species <- unique(issues_eval_test$species_name)
    issues_eval_test_n_species <- length(issues_eval_test_species)

    paste0(
      "\n!! There are issues in testing evaluation data for ",
      issues_eval_test_n_species, " / ", n_species, " species !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(issues_eval_test_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    dplyr::select(issues_eval_test, -species_name) %>%
      dplyr::mutate(
        variable = stringr::str_remove_all(
          variable, "_mss_test|_ess_test|_test"),
        variable = stringr::str_to_lower(variable)) %>%
      dplyr::pull(variable) %>%
      unique() %>%
      paste(collapse = "; ") %>%
      stringr::str_wrap(50L) %>%
      stringr::str_replace_all("\n", "\n  >>>  ") %>%
      paste0(
        crayon::bold("Affected metrics: "),
        "\n  >>>  ", ., collapse = "\n") %>%
      ecokit::cat_time(cat_timestamp = FALSE)
  }

  all_model_results <- dplyr::select(all_model_results, -evaluation_testing)
  invisible(gc())

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Variable importance -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_var_imp <- dplyr::pull(all_model_results, variable_importance) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(species_name, variable) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = c("cor_test", "auc_test"),
        .fns =  ~ sum(is.na(.) | is.nan(.))), .groups = "drop") %>%
    dplyr::filter(
      dplyr::if_any(.cols = c("cor_test", "auc_test"), .fns =  ~ .x > 0L))

  if (nrow(check_var_imp) > 0L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_var_imp_species <- unique(check_var_imp$species_name)
    issues_var_imp_n_species <- length(issues_var_imp_species)

    issues_var_imp <- check_var_imp %>%
      dplyr::select(-species_name) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = c("cor_test", "auc_test"),
          .fns =  ~ sum(.)), .groups = "drop") %>%
      dplyr::mutate(
        message = paste0(
          "cor_test: ", cor_test, "; auc_test: ", auc_test,
          " sp/cv combinations"),
        message = stringr::str_remove_all(message, "^; | $"),
        message = paste0(
          stringr::str_pad(variable, width = 12L, side = "right"),
          " --> ", message)) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n  >>>  ")

    paste0(
      "\n!! There are issues in variable importance data for ",
      issues_var_imp_n_species, " / ", n_species, " species !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(issues_var_imp_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    ecokit::cat_time(
      issues_var_imp, cat_timestamp = FALSE, level = 1L, ... = "\n")
  }

  all_model_results <- dplyr::select(all_model_results, -variable_importance)
  invisible(gc())

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Response curves -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_res_curv <- dplyr::pull(all_model_results, response_curves) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(species_name, variable) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = c("x_value", "prediction"),
        .fns = ~ sum(is.na(.) | is.nan(.))),
      .groups = "drop") %>%
    dplyr::filter(
      dplyr::if_any(.cols = c("x_value", "prediction"), .fns =  ~ .x > 0L))

  if (nrow(check_res_curv) > 0L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_res_curv_species <- unique(check_res_curv$species_name)
    issues_res_curv_n_species <- length(issues_res_curv_species)

    issues_res_curv <- check_res_curv %>%
      dplyr::select(-species_name) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = c("x_value", "prediction"),
          .fns =  ~ sum(.)), .groups = "drop") %>%
      dplyr::mutate(
        message = paste0(
          "missing x values: ", x_value, "; prediction: ", prediction),
        message = stringr::str_remove_all(message, "^; | $"),
        message = paste0(
          stringr::str_pad(variable, width = 12L, side = "right"),
          " --> ", message)) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n  >>>  ")

    paste0(
      "\n!! There are issues in response curves data for ",
      issues_res_curv_n_species, " / ", n_species, " species !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(issues_res_curv_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    ecokit::cat_time(
      issues_res_curv, cat_timestamp = FALSE, level = 1L, ... = "\n")
  }

  all_model_results <- dplyr::select(all_model_results, -response_curves)
  invisible(gc())

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Predictions -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check number of invalid tiff or data files

  check_preds <- dplyr::pull(all_model_results, prediction_data) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!pred_okay)

  if (nrow(check_preds) > 0L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_preds_species <- unique(check_preds$species_name)
    issues_preds_n_species <- length(issues_preds_species)

    paste0(
      "\n!! There are issues in prediction data for ",
      issues_preds_n_species, " / ", n_species, " species !!\n") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE)

    ecokit::cat_time(
      "Affected species: ", cat_timestamp = FALSE, cat_bold = TRUE)
    paste(issues_preds_species, collapse = "; ") %>%
      stringr::str_wrap(width = 65L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      stringr::str_replace_all(" ", " ") %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")

    ecokit::cat_time(
      "Affected climate options: ", cat_timestamp = FALSE, cat_bold = TRUE)
    check_preds %>%
      dplyr::group_by(climate_name) %>%
      dplyr::tally(name = "times") %>%
      dplyr::mutate(message = paste0(climate_name, " (", times, ")")) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1L, ... = "\n")
  }

  if (issue_detected) {
    ecokit::cat_sep(
      line_char_rep = 60L, sep_lines_before = 1L, sep_lines_after = 2L,
      line_char = "=", cat_bold = TRUE, cat_red = TRUE)
  }

  invisible(NULL)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# copy_svm2 -------
# # ========================================================================= #

#' Copy the svm2 methodInfo to the sdm package methods/sdm directory
#'
#' This function checks if the file `svm2.R` exists in the `methods/sdm`
#' subdirectory of the installed sdm package. If not, it writes a new file
#' containing the methodInfo list for registering the svm2 method (using `e1071`
#' SVM) in sdm.
#'
#' @details The methodInfo list defines the `svm2` modeling method for sdm,
#'   using the `e1071::svm` implementation.
#' @return Returns `invisible(NULL)`. The function is called for its side effect
#'   of writing a file.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

copy_svm2 <- function() {
  # Target path in the sdm package
  target_dir <- system.file("methods/sdm", package = "sdm")
  target_file <- fs::path(target_dir, "svm2.R")

  # Check if file already exists
  if (fs::file_exists(target_file)) {
    return(invisible(NULL))
  }

  # MethodInfo function text to write
  method_text <-
    'methodInfo <- list(
  name = c("svm2", "SVM2", "svm_e1071"),
  packages = "e1071",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame", v = "sdmVariables"),
  fitSettings = list(kernel = "radial", probability = TRUE),
  fitFunction = function(formula, data, v, ...) {
    x <- sdm:::.getData.sdmMatrix(
      formula, data, normalize = TRUE, frame = v@varInfo$numeric, scale = FALSE)
    y <- sdm:::.getData.sdmY(formula, data)
    e1071::svm(x = x, y = y, scale = TRUE, ...)

    # Set class weights for pa/pb (binary response)
    n0 <- sum(y == 0, na.rm = TRUE)
    n1 <- sum(y == 1, na.rm = TRUE)
    max_weight <- 20 # Upper bound for weight
    if (n0 >= n1) {
      # More absences
      class.weights <- c("0" = 1, "1" = min(n0 / n1, max_weight))
    } else {
      # More presences
      class.weights <- c("0" = min(n1 / n0, max_weight), "1" = 1)
    }

    e1071::svm(x = x, y = y, scale = TRUE, class.weights = class.weights, ...)
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula", newx = "sdmDataFrame",
    v = "sdmVariables"),
  predictSettings = list(probability = TRUE),
  predictFunction = function(object, formula, newx, v, ...) {
    newx <- sdm:::.getData.sdmMatrix(
      formula, newx, normalize = TRUE,
      frame = v@varInfo$numeric, scale = FALSE)
    predict(object, newx, ...)
  },

  #------ metadata (optional):
  title = "Support Vector Machines using e1071",
  creator = "Ahmed El-Gabbas"
  )'

  # Write to file, preserving formatting
  writeLines(method_text, target_file)
  # message("svm2.R written to sdm package methods/sdm directory.")
  return(invisible(NULL))
}
