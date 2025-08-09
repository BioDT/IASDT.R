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

  predictor_names <- model@setting@featureFrame@predictors
  # Ensure predictor_names is a character vector
  predictor_names <- as.character(predictor_names)
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
    }
  )
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
          dplyr::rename_with(~ stringr::str_c(., "_ess_train"))
        eval_train_thr_mss <- eval_train_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(., "_mss_train"))
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
          dplyr::rename_with(~ stringr::str_c(., "_ess_test"))
        eval_test_thr_mss <- eval_test_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(., "_mss_test"))
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

  Name <- TimePeriod <- ClimateScenario <- ClimateModel <- pred_df <-
    FilePath <- path_rail <- path_roads <- path_clc <- path_bias <- cv <- cvs <-
    path_rivers <- path_chelsa <- cv_fold <- . <- pred_data <- quadratic <-
    variable <- climate_name <- is_valid_option <- species_name <- NULL

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
  if (!fs::file_exists(env_file)) {
    ecokit::stop_ctx(
      "Error: Environment file is invalid or does not exist.",
      env_file = env_file, include_backtrace = TRUE)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_rail", "DP_R_Railways_processed", TRUE, FALSE,
    "path_roads", "DP_R_Roads_processed", TRUE, FALSE,
    "path_clc", "DP_R_CLC_processed", TRUE, FALSE,
    "path_bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "path_rivers", "DP_R_Rivers_processed", TRUE, FALSE,
    "path_chelsa", "DP_R_CHELSA_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  ## Habitat types -----
  hab_abb <- tolower(as.character(hab_abb))

  # Check if `hab_abb` is a single character value
  if (length(hab_abb) != 1L || !nzchar(hab_abb)) {
    ecokit::stop_ctx(
      "`hab_abb` must be a single character value",
      hab_abb = hab_abb, length_hab_abb = length(hab_abb),
      include_backtrace = TRUE)
  }

  valid_hab_abbs <- c(as.character(0L:3L), "4a", "4b", "10", "12a", "12b")
  if (!(hab_abb %in% valid_hab_abbs)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid Habitat abbreviation. Valid values are:\n >> ",
        toString(valid_hab_abbs)),
      hab_abb = hab_abb, include_backtrace = TRUE)
  }

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

  file_model_data <- fs::dir_ls(
    model_dir, regexp = ".*/ModDT_.*subset\\.RData$")
  if (length(file_model_data) != 1 || !nzchar(file_model_data)) {
    ecokit::stop_ctx(
      "Model data file not found",
      model_dir = model_dir, file_model_data = file_model_data)
  }

  if (!ecokit::check_data(file_model_data, warning = FALSE)) {
    ecokit::stop_ctx(
      "Model data at the subset extent file does not exist or invalid",
      model_dir = model_dir, file_model_data = file_model_data,
      model_data_exists = fs::file_exists(file_model_data))
  }

  model_data <- ecokit::load_as(file_model_data)

  if (!all(c("DT_x", "DT_y", "DT_CV") %in% names(model_data))) {
    ecokit::stop_ctx(
      "Loaded model data is missing required components.",
      file_model_data = file_model_data, names_model_data = names(model_data))
  }

  if (is.null(model_data$DT_x) || is.null(model_data$DT_y) ||
      is.null(model_data$DT_CV)) {
    ecokit::stop_ctx(
      "Loaded model data is missing required components.",
      file_model_data = file_model_data, names_model_data = names(model_data))
  }

  predictor_names <- names(model_data$DT_x)

  n_cv_folds <- length(unique(dplyr::pull(model_data$DT_CV, cv_type)))
  if (n_cv_folds < 2L) {
    ecokit::stop_ctx(
      "Not enough CV folds found in model data",
      length_cv_folds = n_cv_folds)
  }

  # `clamp_pred` can not be TRUE when `EffortsLog` is not used as predictor
  if (clamp_pred && !("EffortsLog" %in% predictor_names)) {
    ecokit::stop_ctx(
      "`clamp_pred` can not be used when `EffortsLog` is not used as predictor",
      clamp_pred = clamp_pred, names_data = predictor_names,
      include_backtrace = TRUE)
  }

  species_names <- names(model_data$DT_y)
  chelsa_pattern <- paste0(
    "^", IASDT.R::CHELSA_variables$Variable, collapse = "|")
  other_variables <- stringr::str_subset(
    predictor_names, chelsa_pattern, negate = TRUE)
  bio_variables <- stringr::str_subset( # nolint: object_usage_linter
    predictor_names, chelsa_pattern, negate = FALSE)

  # # ..................................................................... ###

  # Selecting species -----

  if (!is.null(selected_species)) {
    if (!inherits(selected_species, "character") ||
        length(selected_species) < 1L ||
        length(selected_species) > length(names(model_data$DT_y))) {
      ecokit::stop_ctx(
        "selected_species must be a character vector of length >= 1",
        selected_species = selected_species,
        length_selected_species = length(selected_species),
        class_selected_species = class(selected_species))
    }
    if (!all(selected_species %in% names(model_data$DT_y))) {
      invalid_species <- setdiff(selected_species, names(model_data$DT_y))
      ecokit::stop_ctx(
        "Some or all of the selected species not found in model data",
        selected_species = selected_species,
        available_species = names(model_data$DT_y),
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
        length(excluded_species) > length(names(model_data$DT_y))) {
      ecokit::stop_ctx(
        "excluded_species must be a character vector of length >= 1",
        excluded_species = excluded_species,
        length_excluded_species = length(excluded_species),
        class_excluded_species = class(excluded_species))
    }
    if (!all(excluded_species %in% names(model_data$DT_y))) {
      invalid_species <- setdiff(excluded_species, names(model_data$DT_y))
      ecokit::stop_ctx(
        "Some or all of the selected species not found in model data",
        excluded_species = excluded_species,
        available_species = names(model_data$DT_y),
        invalid_species = invalid_species)
    }
    species_names <- setdiff(species_names, excluded_species)
  }

  # # ..................................................................... ###

  # Prepare species data -----

  modelling_data <- cbind.data.frame(
    model_data$DT_y, model_data$DT_x, cv_fold = model_data$DT_CV[[cv_type]])

  # extract list of linear and quadratic terms from the model formula
  # Use terms() to robustly extract variable names from the formula
  term_labels <- attr(stats::terms(model_data$Form_x), "term.labels")
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

  species_modelling_data <- tidyr::expand_grid(
    species_name = species_names, cv = seq_len(n_cv_folds),
    method_is_glm = c(TRUE, FALSE))

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(species_modelling_data)), show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  pkg_to_load <- c(
    "dplyr", "ecokit", "fs", "qs2", "stringr", "sdm", "tibble", "tidyselect")
  future_globals <- c(
    "species_modelling_data", "species_data_dir", "q_preds", "l_preds",
    "predictor_names", "model_data", "modelling_data")

  species_modelling_data2 <- quietly(
    future.apply::future_lapply(
      X = seq_len(nrow(species_modelling_data)),
      FUN = function(line_id) {

        species_name <- species_modelling_data$species_name[[line_id]]
        cv <- species_modelling_data$cv[[line_id]]
        method_is_glm <- species_modelling_data$method_is_glm[[line_id]]

        species_data_name <- paste0(
          species_name, "_cv", cv, "_",
          dplyr::if_else(method_is_glm, "glm", "not_glm"))
        species_data_file <- fs::path(
          species_data_dir, paste0(species_data_name, ".qs2"))

        if (ecokit::check_data(species_data_file, warning = FALSE)) {
          return(
            tibble::tibble(
              data_path = species_data_file, is_valid_option = TRUE))
        }

        # model formula
        if (method_is_glm) {

          if (length(q_preds) == 0L) {
            # No quadratic terms, use linear predictors only
            model_formula <- paste(
              species_name, " ~ ", paste(l_preds, collapse = " + ")) %>%
              stats::as.formula(env = baseenv())
            predictor_names_local <- predictor_names
          } else {
            # Add quadratic terms as columns to the modelling data
            predictor_names_local <- c(
              predictor_names, paste0(q_preds, "_sq"))

            # Update model formula
            model_formula <- c(l_preds, q_preds, paste0(q_preds, "_sq")) %>%
              paste(collapse = " + ") %>%
              paste(species_name, " ~ ", .) %>%
              stats::as.formula(env = baseenv())
          }

        } else {
          predictor_names_local <- predictor_names
          str_rem_1 <- "stats::poly\\(|, degree = 2, raw = TRUE\\)"
          model_formula <- deparse1(model_data$Form_x) %>%
            stringr::str_remove_all(str_rem_1) %>%
            paste(species_name, .) %>%
            stats::as.formula(env = baseenv())
        }

        # Training and testing data
        training_data <- dplyr::filter(modelling_data, cv_fold != cv) %>%
          dplyr::select(
            tidyselect::all_of(c(species_name, predictor_names_local)))
        testing_data <- dplyr::filter(modelling_data, cv_fold == cv) %>%
          dplyr::select(
            tidyselect::all_of(c(species_name, predictor_names_local)))

        valid_training_data <- unique(training_data[, species_name]) %>%
          sort() %>%
          as.integer() %>%
          identical(c(0L, 1L))
        valid_testing_data <- unique(testing_data[, species_name]) %>%
          sort() %>%
          as.integer() %>%
          identical(c(0L, 1L))

        if (isFALSE(valid_training_data) || isFALSE(valid_testing_data)) {
          return(
            tibble::tibble(
              data_path = NA_character_, is_valid_option = FALSE))
        }

        # sdm data
        sdm_data <- sdm::sdmData(
          formula = model_formula, train = training_data, test = testing_data)

        species_fitting_data <- list(
          species_name = species_name, cv = cv,
          method_is_glm = method_is_glm, model_formula = model_formula,
          sdm_data = sdm_data)

        ecokit::save_as(
          object = species_fitting_data, out_path = species_data_file)

        tibble::tibble(data_path = species_data_file, is_valid_option = TRUE)

      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_load, future.globals = future_globals)
  )

  ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
  future::plan("sequential", gc = TRUE)
  invisible(gc())

  species_modelling_data <- dplyr::bind_cols(
    species_modelling_data, dplyr::bind_rows(species_modelling_data2))

  invalid_species <- dplyr::filter(species_modelling_data, !is_valid_option) %>%
    dplyr::distinct(species_name, cv)

  if (nrow(invalid_species) > 0) {

    excluded_species <- unique(invalid_species$species_name)
    excluded_species_n <- length(excluded_species)

    paste0(
      "\n!! There are ", excluded_species_n,
      " invalid species that will be excluded in all model types !!") %>%
      crayon::blue() %>%
      ecokit::cat_time(cat_timestamp = FALSE, cat_bold = TRUE, ... = "\n")

    tidyr::nest(invalid_species, cvs = cv) %>%
      dplyr::mutate(
        message = purrr::map2_chr(
          .x = species_name, .y = cvs,
          .f = ~ {
            collapsed_cvs <- paste(sort(unlist(.y)), collapse = " & ")
            paste0(.x, " (cv: ", collapsed_cvs, ")")
          })) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n  >>>  ") %>%
      ecokit::cat_time(cat_timestamp = FALSE, level = 1, ... = "\n")

    species_modelling_data <- species_modelling_data %>%
      dplyr::filter(!species_name %in% excluded_species) %>%
      dplyr::select(-is_valid_option)

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
    r_road_rail <- (r_roads + r_railways) %>%
      magrittr::add(0.1) %>%
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

    r_efforts <- ecokit::load_as(r_efforts, unwrap_r = TRUE) %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(0.1) %>%
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
                  .fns = ~ I(.x^2L), .names = "{.col}_sq"))
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
    rf = list(ntree = 3000L, nodesize = 5L),
    ranger = list(
      num.trees = 3000L, importance = "impurity", min.node.size = 5L),
    cart = list(),
    rpart = list(),
    maxent = list(
      removeDuplicates = FALSE,
      args = c(
        "maximumiterations=4000",
        # "betamultiplier=1",
        "convergencethreshold=0.000001",
        "noautofeature", "hinge", "linear",
        "noproduct", "noquadratic", "nothreshold",
        "writeplotdata", "-J", "-P",
        "noremoveduplicates", "noaddsamplestobackground",
        "doclamp", "nowriteclampgrid", "outputformat=cloglog")),
    mlp = list(maxit = 2000L),
    # rbf = list(maxit = 2000L),
    svm = list(),
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
#' @param model_results_dir Character. Path to the directory for saving fitted
#'   model outputs.
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
    line_id, sdm_method, model_data, model_settings, model_results_dir,
    input_data, output_directory, path_grid_r, copy_maxent_html = TRUE) {

  pred_data <- climate_name <- NULL

  species_name <- model_data$species_name[[line_id]]
  cv_fold <- model_data$cv[[line_id]]
  model_DT <- ecokit::load_as(model_data$data_path[[line_id]])
  base_model_name <- paste0(sdm_method, "_", species_name, "_cv", cv_fold)

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Model fit -----

  model_name <- paste0(base_model_name, "_model")
  model_path <- fs::path(model_results_dir, paste0(model_name, ".RData"))

  if (ecokit::check_data(model_path, warning = FALSE)) {
    fitted_model <- ecokit::load_as(model_path)
  } else {

    fitted_model <- quietly(
      sdm::sdm(
        formula = model_DT$model_formula, data = model_DT$sdm_data,
        methods = sdm_method, modelSettings = model_settings)
    )

    # copy model files from temp dir for maxent models
    if (sdm_method == "maxent" && copy_maxent_html) {
      maxent_html <- ecokit::normalize_path(
        fitted_model@models[[1L]][[1L]][[1L]]@object@html)

      if (fs::file_exists(maxent_html)) {
        out_maxent_dir <- fs::path(
          fs::path_dir(model_results_dir), "maxent_html", base_model_name)
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
    }

    # Reduce models objects
    fitted_model <- reduce_sdm_formulas(obj = fitted_model)

    ecokit::save_as(
      object = fitted_model, object_name = model_name, out_path = model_path)
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Extract model info / predict ------

  extracted_data_name <- paste0(base_model_name, "_extracted_data")
  extracted_data_path <- fs::path(
    model_results_dir, paste0(extracted_data_name, ".RData"))

  if (ecokit::check_data(extracted_data_path, warning = FALSE)) {
    extracted_data <- ecokit::load_as(extracted_data_path)
  } else {

    ## Extract info from fitted model object -----
    extracted_data <- quietly(
      extract_sdm_info(model = fitted_model, cv_fold = cv_fold)
    )

    # Making predictions and extract prediction paths -----
    predictor_names <- fitted_model@setting@featureFrame@predictors
    if (is.null(predictor_names) || length(predictor_names) == 0L) {
      ecokit::stop_ctx("No predictor names found in fitted model.")
    }

    prediction_info <- ecokit::load_as(input_data$path_prediction_data) %>%
      dplyr::mutate(
        prediction = purrr::map2(
          .x = climate_name,
          .y = pred_data,
          .f = ~{

            pred_dir <- fs::path(output_directory, paste0("pred_", .x))
            tif_path <- fs::path(pred_dir, paste0(base_model_name, ".tif"))
            tif_okay <- ecokit::check_tiff(tif_path, warning = FALSE)
            data_path <- fs::path(pred_dir, paste0(base_model_name, ".RData"))
            data_okay <- ecokit::check_data(data_path, warning = FALSE)
            pred_data_df <- as.data.frame(.y)
            gdal_options <- c("COMPRESS=DEFLATE", "TILED=YES")

            if (isFALSE(tif_okay) || isFALSE(data_okay)) {

              pred <- quietly(
                predict(
                  object = fitted_model,
                  newdata = dplyr::select(
                    pred_data_df, tidyselect::all_of(predictor_names)))
              ) %>%
                unlist() %>%
                unname()

              pred_r <- tibble::tibble(pred = pred) %>%
                dplyr::bind_cols(dplyr::select(pred_data_df, x, y)) %>%
                sf::st_as_sf(coords = c("x", "y"), crs = 3035L) %>%
                terra::rasterize(
                  y = ecokit::load_as(path_grid_r, unwrap_r = TRUE),
                  field = "pred", fun = "mean", na.rm = TRUE) %>%
                stats::setNames(base_model_name)

              terra::writeRaster(
                x = pred_r, overwrite = TRUE, filename = tif_path,
                gdal = gdal_options)
              ecokit::save_as(
                object = terra::wrap(pred_r), object_name = base_model_name,
                out_path = data_path)

              tif_okay <- ecokit::check_tiff(tif_path, warning = FALSE)
              data_okay <- ecokit::check_data(data_path, warning = FALSE)

              invisible(gc())
            }

            tibble::tibble(
              species_name = species_name, pred_dir = pred_dir,
              data_path = data_path, data_okay = data_okay,
              tif_path = tif_path, tif_okay = tif_okay) %>%
              dplyr::mutate(cv_fold = cv_fold, .before = 1L)

          })) %>%
      dplyr::select(-pred_data) %>%
      tidyr::unnest("prediction") %>%
      dplyr::select(species_name, cv_fold, tidyselect::everything())

    extracted_data <- c(extracted_data, prediction_info = list(prediction_info))

    ecokit::save_as(
      object = extracted_data, object_name = extracted_data_name,
      out_path = extracted_data_path)
  }

  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model results -----
  species_results <- lapply(extracted_data, list) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(model_path = model_path, .before = 1L)

  species_results
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
#' @param line_id Integer. Index (row number) in `model_summary` to process for
#'   summary statistics.
#' @param model_summary List. Model summary structure containing:
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

summarize_predictions <- function(line_id, model_summary) {

  climate_name <- cv_fold <- tiff_paths <- NULL

  DT <- model_summary$summary_data[[line_id]]
  species_name <- model_summary$species_name[[line_id]]
  mean_auc <- model_summary$evaluation_testing[[line_id]]$auc_test # nolint: object_usage_linter

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

  exclude_cols <- c(
    "data_path", "data_okay", "tif_okay", "cv_fold", "species_name")
  keep_cols <- c(
    "time_period", "climate_model", "climate_scenario",
    "climate_name", "pred_dir")

  preds_dt_orig <- dplyr::select(DT, "prediction_info") %>%
    tidyr::unnest("prediction_info") %>%
    dplyr::mutate(cv_fold = as.character(cv_fold))

  preds_dt <- preds_dt_orig %>%
    dplyr::select(-tidyselect::all_of(exclude_cols)) %>%
    tidyr::nest(
      .by = tidyselect::all_of(keep_cols), .key = "tiff_paths")

  pred_summ <- purrr::map(
    .x = preds_dt$tiff_paths,
    .f = ~ {

      tif_path <- NULL
      tif_paths <- unname(unlist(.x))

      # Check if all input tiff files are valid
      in_tiffs_okay <- purrr::map_lgl(
        tif_paths, ecokit::check_tiff, warning = FALSE) %>%
        all()
      if (isFALSE(in_tiffs_okay)) {
        ecokit::stop_ctx(
          "Not all input tiff files are valid", tif_paths = tif_paths)
      }

      path_mean_tif <- stringr::str_replace_all(
        tif_paths, "_cv[1-9]", "_mean") %>%
        unique()
      path_mean_data <- stringr::str_replace_all(
        path_mean_tif, ".tif$", ".RData")
      path_sd_tif <- stringr::str_replace_all(
        path_mean_tif, "mean", "sd")
      path_sd_data <- stringr::str_replace_all(
        path_mean_data, "mean", "sd")
      path_cov_tif <- stringr::str_replace_all(
        path_mean_tif, "mean", "cov")
      path_cov_data <- stringr::str_replace_all(
        path_mean_data, "mean", "cov")

      # Check if output maps already exist and valid
      tiffs_to_check <- c(path_mean_tif, path_sd_tif, path_cov_tif)
      data_to_check <- c(path_mean_data, path_sd_data, path_cov_data)
      if (calc_w_mean) {
        path_w_mean_tif <- stringr::str_replace_all(
          path_mean_tif, "_mean", "_weighted_mean")
        tiffs_to_check <- c(tiffs_to_check, path_w_mean_tif)

        path_w_mean_data <- stringr::str_replace_all(
          path_mean_data, "_mean", "_weighted_mean")
        data_to_check <- c(data_to_check, path_w_mean_data)
      } else {
        path_w_mean_tif <- path_w_mean_data <- NA_character_
      }

      tif_okay <- purrr::map_lgl(
        tiffs_to_check, ecokit::check_tiff, warning = FALSE) %>%
        all()
      data_okay <- purrr::map_lgl(
        data_to_check, ecokit::check_data, warning = FALSE) %>%
        all()

      check_data_tif <- function(DT, tiff) {
        ecokit::check_data(DT, warning = FALSE) &&
          ecokit::check_tiff(tiff, warning = FALSE)
      }

      if (isFALSE(tif_okay && data_okay)) {

        maps <- terra::rast(tif_paths)
        gdal_options <- c("COMPRESS=DEFLATE", "TILED=YES")

        # |||||||||||||||||||||||||||||||||||||||||||

        # mean -------

        if (check_data_tif(path_mean_data, path_mean_tif)) {
          # If mean data already exists, load it
          mean_pred <- ecokit::load_as(path_mean_data, unwrap_r = TRUE)
        } else {
          mean_name <- basename(path_mean_tif) %>%
            tools::file_path_sans_ext()
          mean_pred <- terra::app(maps, mean, na.rm = TRUE) %>%
            stats::setNames(mean_name)
          terra::writeRaster(
            x = mean_pred, overwrite = TRUE, filename = path_mean_tif,
            gdal = gdal_options)
          ecokit::save_as(
            object = terra::wrap(mean_pred), object_name = mean_name,
            out_path = path_mean_data)
        }

        # |||||||||||||||||||||||||||||||||||||||||||

        # weighted mean -----

        if (calc_w_mean) {

          if (check_data_tif(path_w_mean_data, path_w_mean_tif)) {
            # If weighted mean data already exists, load it
            w_mean_pred <- ecokit::load_as(path_w_mean_data, unwrap_r = TRUE)
          } else {
            w_mean_name <- basename(path_w_mean_tif) %>%
              tools::file_path_sans_ext()
            w_mean_pred <- terra::weighted.mean(
              x = maps, w = mean_auc, na.rm = TRUE) %>%
              stats::setNames(w_mean_name)
            terra::writeRaster(
              x = w_mean_pred, overwrite = TRUE, filename = path_w_mean_tif,
              gdal = gdal_options)
            ecokit::save_as(
              object = terra::wrap(w_mean_pred), object_name = w_mean_name,
              out_path = path_w_mean_data)
          }
        }

        # |||||||||||||||||||||||||||||||||||||||||||

        # sd ------

        if (check_data_tif(path_sd_data, path_sd_tif)) {
          # If standard deviation data already exists, load it
          sd_pred <- ecokit::load_as(path_sd_data, unwrap_r = TRUE)
        } else {
          sd_name <- tools::file_path_sans_ext(basename(path_sd_tif))
          sd_pred <- terra::app(maps, sd, na.rm = TRUE) %>%
            stats::setNames(sd_name)
          terra::writeRaster(
            x = sd_pred, overwrite = TRUE, filename = path_sd_tif,
            gdal = gdal_options)
          ecokit::save_as(
            object = terra::wrap(sd_pred), object_name = sd_name,
            out_path = path_sd_data)
        }

        # |||||||||||||||||||||||||||||||||||||||||||

        # cov ------

        if (check_data_tif(path_cov_data, path_cov_tif)) {
          # If coefficient of variation data already exists, load it
          cov_pred <- ecokit::load_as(path_cov_data, unwrap_r = TRUE)
        } else {
          cov_name <- basename(path_cov_tif) %>%
            tools::file_path_sans_ext()
          # Avoid division by zero
          mean_pred <- terra::clamp(mean_pred, lower = 1e-8)
          cov_pred <- stats::setNames((sd_pred / mean_pred), cov_name)
          terra::writeRaster(
            x = cov_pred, overwrite = TRUE, filename = path_cov_tif,
            gdal = gdal_options)
          ecokit::save_as(
            object = terra::wrap(cov_pred), object_name = cov_name,
            out_path = path_cov_data)
        }

      }

      # |||||||||||||||||||||||||||||||||||||||||||

      tibble::tribble(
        ~cv_fold, ~data_path, ~tif_path,
        "mean", path_mean_data, path_mean_tif,
        "weighted_mean", path_w_mean_data, path_w_mean_tif,
        "sd", path_sd_data, path_sd_tif,
        "cov", path_cov_data, path_cov_tif) %>%
        dplyr::mutate(
          data_okay = purrr::map_lgl(
            data_path, ecokit::check_data, warning = FALSE),
          tif_okay = purrr::map_lgl(
            tif_path, ecokit::check_tiff, warning = FALSE))
    })

  preds_dt <- dplyr::mutate(preds_dt, pred_summ = pred_summ) %>%
    tidyr::unnest("pred_summ") %>%
    dplyr::select(-tiff_paths) %>%
    dplyr::mutate(species_name = species_name, .before = 1L)

  dplyr::bind_rows(preds_dt_orig, preds_dt) %>%
    dplyr::arrange(climate_name, cv_fold)
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
#'
#' @details
#' - The function outputs information about any identified issues to the
#' console, including which species are affected and what specific metrics or
#' data elements have issues.
#' - The function is not expected to be called directly by users, but is
#' used internally within the `IASDT.R` workflow, particularly the
#' [fit_sdm_models()] function.

#' @return Returns invisibly NULL. This function is used for its side effects
#'   (console output).
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

check_model_results <- function(model_results) {

  data_okay <- tif_okay <- climate_name <- data_path <- cv_fold <- x_value <-
    tif_path <- climate_model <- climate_scenario <- time_period <- na_count <-
    prediction <- cor_test <- auc_test <- species_name <- . <-  pred_dir <-
    variable <- min_value <- max_value <- times <- bad_vals <- NULL

  n_species <- length(unique(model_results$species_name))

  # Check if model_results is a data frame with at least one row
  if (!inherits(model_results, "data.frame") || nrow(model_results) == 0L) {
    ecokit::stop_ctx(
      "model_results must be a data frame with at least 1 row.",
      class_model_results = class(model_results))
  }

  # Check if model_results contains the required columns
  required_cols <- c(
    "evaluation_training", "evaluation_testing", "variable_importance",
    "response_curves", "prediction_info")
  if (!all(required_cols %in% names(model_results))) {
    missing_columns <- required_cols[!required_cols %in% names(model_results)]
    ecokit::stop_ctx(
      paste0(
        "`model_results` must contain the following columns: ",
        toString(missing_columns)),
      names_model_results = names(model_results))
  }

  invalid_cols <- purrr::map_lgl(
    .x = required_cols,
    .f = ~ {
      example_data <- model_results[[.x]][[1L]]
      inherits(example_data, "tbl_df") && nrow(example_data) > 0L
    })
  invalid_columns <- required_cols[!invalid_cols]
  if (length(invalid_columns) > 0L) {
    ecokit::stop_ctx(
      paste0(
        "model_results components must be tibbles with nrow > 0: ",
        toString(invalid_columns)),
      names_model_results = names(model_results))
  }

  issue_detected <- FALSE

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation data - Training -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_eval_train <- model_results$evaluation_training %>%
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

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation data - Testing -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_eval_test <- model_results$evaluation_testing %>%
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

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Variable importance -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_var_imp <- dplyr::bind_rows(model_results$variable_importance) %>%
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
          "cor_test: ", cor_test, " sp; auc_test: ", auc_test, " sp"),
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

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Response curves -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  check_res_curv <- dplyr::bind_rows(model_results$response_curves) %>%
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

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Predictions -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check number of invalid tiff or data files

  check_preds <- dplyr::bind_rows(model_results$prediction_info) %>%
    dplyr::filter(!data_okay | !tif_okay)

  if (nrow(check_preds) > 0L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    issues_preds_species <- unique(check_preds$species_name)
    issues_preds_n_species <- length(issues_preds_species)

    issues_preds <- check_preds %>%
      dplyr::select(
        -time_period, -climate_model, -climate_scenario,
        -tif_path, -pred_dir, -data_path, -cv_fold) %>%
      dplyr::group_by(climate_name) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = c("tif_okay", "data_okay"), .fns =  ~ sum(!.)),
        .groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        message = paste0(tif_okay, " tiff files; ", data_okay, " data files"),
        message = stringr::str_remove_all(message, "^; | $"),
        message = paste0(
          stringr::str_pad(climate_name, width = 35L, side = "right"),
          " --> ", message)) %>%
      dplyr::pull(message) %>%
      paste(collapse = "\n  >>>  ")

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
    ecokit::cat_time(
      issues_preds, cat_timestamp = FALSE, level = 1L, ... = "\n")
  }

  # Check if predictions are < 0 or > 1
  pred_odd_vals <- dplyr::bind_rows(model_results$prediction_info) %>%
    dplyr::select(climate_name, species_name, tif_path) %>%
    dplyr::mutate(
      min_max = purrr::map(
        .x = tif_path,
        .f = ~{
          terra::global(terra::rast(.x), range, na.rm = TRUE) %>%
            stats::setNames(c("min_value", "max_value"))
        })) %>%
    tidyr::unnest("min_max") %>%
    dplyr::filter(min_value < -0.001 | max_value > 1.001)

  if (nrow(pred_odd_vals) > 0L) {

    if (isFALSE(issue_detected)) {
      ecokit::cat_sep(
        line_char_rep = 60L, sep_lines_before = 1L,
        line_char = "=", cat_bold = TRUE, cat_red = TRUE)
      issue_detected <- TRUE
    }

    pred_odd_vals %>%
      dplyr::group_by(climate_name, species_name) %>%
      dplyr::tally(name = "times") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(message = paste0(species_name, " (", times, ")")) %>%
      dplyr::select(climate_name, message) %>%
      tidyr::nest(bad_vals = -climate_name) %>%
      dplyr::mutate(
        bad_vals = purrr::map2_chr(
          .x = bad_vals, .y = climate_name,
          .f = ~ {
            sp_l <- paste(unlist(.x), collapse = "; ") %>%
              stringr::str_wrap(width = 40L) %>%
              stringr::str_split(pattern = "\n", simplify = TRUE) %>%
              paste0("  >>>  ", ., collapse = "\n")
            paste0(.y, ": ") %>%
              crayon::bold() %>%
              paste0("\n", sp_l)
          })) %>%
      dplyr::pull(bad_vals) %>%
      paste(collapse = "\n") %>%
      ecokit::cat_time(cat_timestamp = FALSE)
  }

  if (issue_detected) {
    ecokit::cat_sep(
      line_char_rep = 60L, sep_lines_before = 1L, sep_lines_after = 2L,
      line_char = "=", cat_bold = TRUE, cat_red = TRUE)
  }

  invisible(NULL)
}

# ............................................... ----

# # ========================================================================= #
# quietly ------
# # ========================================================================= #

#' Quietly Evaluate an Expression
#'
#' Evaluates an R expression while suppressing package startup messages and
#' selected warnings. Specifically, warnings containing "was built under R
#' version" or "Loading required namespace" are muffled, allowing for cleaner
#' output during package loading or function execution.
#'
#' @param expr An R expression to be evaluated quietly.
#'
#' @return The result of evaluating \code{expr}, with specified messages and
#' warnings suppressed.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

quietly <- function(expr) {
  withCallingHandlers(
    suppressPackageStartupMessages(expr),
    warning = function(w) {
      if (grepl(
        "was built under R version|Loading required namespace",
        conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
}
