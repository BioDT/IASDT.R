## |------------------------------------------------------------------------| #
# predict_maps ----
## |------------------------------------------------------------------------| #

#' Predict habitat suitability of `Hmsc` models
#'
#' This package provides two functions for predicting habitat suitability of
#' `Hmsc` models in the `IASDT` framework. [predict_maps] generates current and
#' future habitat suitability maps (mean, sd, cov) from a full Hmsc model fit.
#' [predict_maps_CV] predicts and evaluates cross-validated Hmsc models for
#' current climate conditions. For more details, see the respective function
#' documentation and the details section below.
#'
#' @param path_model Character. Path to fitted `Hmsc` model object.
#' @param hab_abb Character. Habitat abbreviation indicating the specific
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type. Valid
#'   values: `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`. See [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447) for details.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param clamp_pred Logical indicating whether to clamp the sampling efforts at
#'   a single value. If `TRUE` (default), the `fix_efforts` argument must be
#'   provided.
#' @param fix_efforts Numeric or character. When `clamp_pred = TRUE`, fixes the
#'   sampling efforts predictor at this value during predictions. If numeric,
#'   uses the value directly (on log<sub>10</sub> scale). If character, must be
#'   one of `identity` (i.e., do not fix), `median`, `mean`, `max`, or `q90`
#'   (90% quantile). Using `max` may reflect extreme sampling efforts from
#'   highly sampled locations, while `q90` captures high sampling areas without
#'   extremes. Required if `clamp_pred = TRUE`.
#' @param fix_rivers Numeric, character, or `NULL`. Similar to `fix_efforts`,
#'   but for the river length predictor. If `NULL`, the river length is not
#'   fixed. Default: `q90`.
#' @param pred_new_sites Logical. Whether to predict suitability at new sites.
#'   Default: `TRUE`.
#' @param LF_only Logical. Whether to predict only the latent factor. This is
#'   useful for distributing processing load between GPU and CPU. When `LF_only
#'   = TRUE`, latent factor prediction needs to be computed separately on GPU.
#'   When computations are finished on GPU, the function can later be rerun with
#'   `LF_only = FALSE` (default) to predict habitat suitability using the
#'   already-computed latent factor predictions.
#' @param CC_models Character vector. Climate models for future predictions.
#'   Available options are `c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
#'   "MRI-ESM2-0", "UKESM1-0-LL")` (default).
#' @param CC_scenario Character vector. Climate scenarios for future
#'   predictions. Available options are: `c("ssp126", "ssp370", "ssp585")`
#'   (default).
#' @param tar_predictions Logical. Whether to compress the add files into a
#'   single `*.tar` file (without compression). Default: `TRUE`.
#' @param model_dir Character. Path to the directory containing cross-validated
#'   models.
#' @param CV_name Character. Cross-validation strategy. Valid values are
#'   `CV_Dist`, `CV_Large`, or `CV_SAC`.
#' @param CV_fold Integer. The cross-validation fold number.
#' @inheritParams predict_hmsc
#' @details
#'   - **`predict_maps`**: Generates habitat suitability maps for `Hmsc` models
#' fitted on the full dataset, for both current and future climate options. It
#' produces maps for mean, standard deviation (sd), and coefficient of variation
#' (cov) of suitability for each species and overall species richness. It
#' evaluate model's explanatory power using various metrics. For future
#' predictions, it also generates anomaly maps (future - current). The function
#' supports ensemble predictions across multiple climate models and prepares
#' data for upload to the [OPeNDAP](http://opendap.biodt.eu/ias-pdt/) server for
#' use in the `IASDT` [Shiny App](https://app.biodt.eu/).
#'   - **`predict_maps_CV`**: Computes predictions for cross-validated `Hmsc`
#' models using only the testing folds. It evaluates model performance
#' (predictive power) with various metrics and plots evaluation results for
#' predictive and explanatory power. Unlike `predict_maps`, this function does
#' not perform clamping and does not generate future climate predictions.
#' @seealso predict_hmsc
#' @importFrom foreach %dopar%
#' @export
#' @name predict_maps
#' @rdname predict_maps
#' @order 1
#' @author Ahmed El-Gabbas

predict_maps <- function(
    path_model = NULL, hab_abb = NULL, env_file = ".env", n_cores = 8L,
    strategy = "multisession", clamp_pred = TRUE, fix_efforts = "q90",
    fix_rivers = "q90", pred_new_sites = TRUE, use_TF = TRUE, TF_environ = NULL,
    TF_use_single = FALSE, LF_n_cores = n_cores, LF_check = FALSE,
    LF_temp_cleanup = TRUE, LF_only = FALSE, LF_commands_only = FALSE,
    temp_dir = "TEMP_Pred", temp_cleanup = TRUE, tar_predictions = TRUE,
    spatial_model = TRUE,
    CC_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_scenario = c("ssp126", "ssp370", "ssp585")) {

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential")  n_cores <- LF_n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)
  LF_n_cores <- .validate_n_cores(LF_n_cores)

  # # ..................................................................... ###
  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  hab_name <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", hab_abb, "_")) %>%
    stringr::str_remove(paste0("^", hab_abb, "_")) %>%
    stringr::str_replace_all("_", " ")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----

  ecokit::cat_time("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("hab_abb", "env_file", "path_model", "temp_dir"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "clamp_pred", "pred_new_sites", "use_TF", "TF_use_single", "LF_check",
      "LF_temp_cleanup", "LF_only", "LF_commands_only", "temp_cleanup",
      "tar_predictions"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores"))

  rm(AllArgs, envir = environment())

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- Path_Roads <- Path_Rail <- Path_Bias <- tif_path <-
    time_period <- climate_scenario <- Path_CHELSA <- row_id <-
    Path_CLC <- ias_id <- taxon_name <- ClimateModel <- ClimateScenario <-
    Name <- File_Pred_sf <- class <- order <- family <- species_name <-
    tif_path_mean <- tif_path_cov <- tif_path_sd <- Clamp <- Path_Wetness <-
    Train <- Ensemble_File <- Ensemble_Maps <- tifs <- layer_name <-
    TimePeriod <- File_Pred_summary <- Ensemble_DT <- Dir_Ensemble <-
    File_Pred_R <- tif_path_anomaly <- Path_Rivers <- Path_Soil <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  Path_Eval <- fs::path(dirname(dirname(path_model)), "Model_Evaluation")

  Path_Prediction1 <- fs::path(
    dirname(dirname(path_model)), "Model_Prediction")
  Path_Prediction_Clamp <- fs::path(Path_Prediction1, "Clamp")
  Path_Prediction_NoClamp <- fs::path(Path_Prediction1, "NoClamp")
  fs::dir_create(c(Path_Eval, Path_Prediction_NoClamp))

  # Path for overall summary - paths for summaries of identical scenarios
  Path_Summary_RData <- fs::path(
    dplyr::if_else(
      clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary.RData")
  Path_Summary_txt <- fs::path(
    dplyr::if_else(
      clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary.txt")

  # Path for overall summary - for ShinyApp
  Path_Summary_RData_Shiny <- fs::path(
    dplyr::if_else(
      clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary_Shiny.RData")
  Path_Summary_txt_Shiny <- fs::path(
    dplyr::if_else(
      clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary_Shiny.txt")

  # Check if the prediction summary is already available on disk
  if (all(
    fs::file_exists(
      c(
        Path_Summary_RData, Path_Summary_txt,
        Path_Summary_RData_Shiny, Path_Summary_txt_Shiny)))) {
    ecokit::cat_time(
      paste0(
        "All model predictions and prediction summary are already available ",
        "on disk"))
    return(invisible(NULL))
  }


  if (clamp_pred) {

    # Check if the `fix_efforts` value is valid
    if (is.null(fix_efforts)) {
      ecokit::stop_ctx(
        "`fix_efforts` can not be `NULL` when Clamping is implemented",
        fix_efforts = fix_efforts, include_backtrace = TRUE)
    }
    if (length(fix_efforts) != 1) {
      ecokit::stop_ctx(
        "`fix_efforts` must be of length 1.",
        fix_efforts = fix_efforts, include_backtrace = TRUE)
    }

    # Check if the `fix_rivers` value is valid
    if (is.null(fix_rivers)) {
      ecokit::stop_ctx(
        "`fix_rivers` can not be `NULL` when Clamping is implemented",
        fix_rivers = fix_rivers, include_backtrace = TRUE)
    }
    if (length(fix_rivers) != 1) {
      # Check if fix_rivers is a vector or length 1
      ecokit::stop_ctx(
        "`fix_rivers` must be of length 1.", fix_rivers = fix_rivers,
        include_backtrace = TRUE)
    }

    # Create folder for clamp results only if clamp_pred == TRUE
    fs::dir_create(Path_Prediction_Clamp)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ------

  ecokit::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Rail", "DP_R_Railways_processed", TRUE, FALSE,
    "Path_Roads", "DP_R_Roads_processed", TRUE, FALSE,
    "Path_CLC", "DP_R_CLC_processed", TRUE, FALSE,
    "Path_Bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Rivers", "DP_R_Rivers_processed", TRUE, FALSE,
    "Path_Wetness", "DP_R_wetness_processed", FALSE, TRUE,
    "Path_Soil", "DP_R_soil_density", FALSE, TRUE,
    "Path_CHELSA", "DP_R_CHELSA_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading input data -----
  ecokit::cat_time("Loading input data")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Species information -----

  SpeciesInfo <- IASDT.R::get_species_name(env_file = env_file) %>%
    janitor::clean_names() %>%
    dplyr::select(ias_id, taxon_name, species_name, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Reference grid -----

  ecokit::cat_time("Reference grid", level = 1L)
  Path_GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!fs::file_exists(Path_GridR)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", Path_GridR = Path_GridR,
      include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Model object -----

  ecokit::cat_time("Model object", level = 1L)

  if (is.null(path_model) || !fs::file_exists(path_model)) {
    ecokit::stop_ctx(
      "Model path is NULL or does not exist ", path_model = path_model,
      include_backtrace = TRUE)
  }

  Model <- ecokit::load_as(path_model)

  # `clamp_pred` can not be TRUE when `EffortsLog` is not used as predictor
  if (clamp_pred && isFALSE("EffortsLog" %in% names(Model$XData))) {
    ecokit::stop_ctx(
      "`clamp_pred` can not be used when `EffortsLog` is not used as predictor",
      clamp_pred = clamp_pred, names_data = names(Model$XData),
      include_backtrace = TRUE)
  }

  other_variables <- paste0(
    "^", c(IASDT.R::CHELSA_variables$Variable, "CV"), collapse = "|") %>%
    stringr::str_subset(names(Model$XData), ., negate = TRUE)
  bio_variables <- paste0(
    "^", IASDT.R::CHELSA_variables$Variable, collapse = "|") %>%
    stringr::str_subset(names(Model$XData), ., negate = FALSE)

  # Coordinates of the model sampling units
  if (is.null(Model$rL) || is.null(Model$rL[[1]]$s)) {
    suppressWarnings({
      Model_Coords <- tibble::tibble(x = numeric(0L), y = numeric(0L)) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035L) %>%
        dplyr::mutate(Train = logical(0L))
    })
  } else {
    Model_Coords <- as.data.frame(Model$rL$sample$s) %>%
      tibble::tibble() %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      dplyr::mutate(Train = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## CHELSA data -----

  ecokit::cat_time("CHELSA data", level = 1L)

  Path_CHELSA <- fs::path(Path_CHELSA, "CHELSA_Processed_DT.RData")
  if (!fs::file_exists(Path_CHELSA)) {
    ecokit::stop_ctx(
      "Processed CHELSA data can not be found", Path_CHELSA = Path_CHELSA,
      include_backtrace = TRUE)
  }

  Prediction_Options <- ecokit::load_as(Path_CHELSA) %>%
    dplyr::select(-"File_List") %>%
    dplyr::filter(
      # Filter only selected future climate models
      ClimateModel %in% c("Current", CC_models),
      # Filter only selected future climate scenarios
      ClimateScenario %in% c("Current", CC_scenario)) %>%
    dplyr::mutate(
      Name = paste0(TimePeriod, "_", ClimateScenario, "_", ClimateModel),
      Name = stringr::str_replace(Name, "1981-2010_Current_Current", "Current"),
      Name = stringr::str_replace_all(Name, "-", "_"))

  if (clamp_pred) {
    # Add a prediction options for no clamping predictions (model evaluation)
    ClampOption <- Prediction_Options %>%
      dplyr::filter(ClimateModel == "Current") %>%
      dplyr::mutate(Clamp = FALSE)
    Prediction_Options <- dplyr::mutate(Prediction_Options, Clamp = TRUE) %>%
      dplyr::bind_rows(ClampOption, .)
  } else {
    Prediction_Options <- dplyr::mutate(Prediction_Options, Clamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Road and railway intensity -----

  StaticPredictors <- list()

  if ("RoadRailLog" %in% other_variables) {

    ecokit::cat_time("Road and railway intensity", level = 1L)

    R_Railways <- fs::path(Path_Rail, "Railways_Length.RData")
    if (!fs::file_exists(R_Railways)) {
      ecokit::stop_ctx(
        "Railways data does not exist", R_Railways = R_Railways,
        include_backtrace = TRUE)
    }
    R_Railways <- ecokit::load_as(R_Railways, unwrap_r = TRUE) %>%
      magrittr::extract2("rail")

    R_Roads <- fs::path(Path_Roads, "Road_Length.RData")
    if (!fs::file_exists(R_Roads)) {
      ecokit::stop_ctx(
        "Roads data does not exist", R_Roads = R_Roads,
        include_backtrace = TRUE)
    }
    R_Roads <- ecokit::load_as(R_Roads, unwrap_r = TRUE) %>%
      magrittr::extract2("All")

    # Calculating the sum of road and railway intensity
    # add 1 (older versions 0.1) to get log for 0 values [only for
    # rivers/roads/efforts, not hab/rivers]
    R_RoadRail <- (R_Roads + R_Railways) %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("RoadRailLog")

    StaticPredictors <- c(StaticPredictors, R_RoadRail)
    rm(R_RoadRail, R_Roads, R_Railways, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Habitat information ----

  # Check if habitat information is used as predictor
  Hab_Predictor <- "HabLog" %in% other_variables

  if (Hab_Predictor) {

    ecokit::cat_time("Habitat information", level = 1L)
    R_Hab <- fs::path(
      Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
    if (!fs::file_exists(R_Hab)) {
      ecokit::stop_ctx(
        "Habitat data does not exist", R_Hab = R_Hab, include_backtrace = TRUE)
    }

    R_Hab <- ecokit::load_as(R_Hab, unwrap_r = TRUE) %>%
      magrittr::extract2(paste0("SynHab_", hab_abb))

    # Models are trained and predictions are made only at grid cells with > 0 %
    # coverage. Mask layer to exclude grid cells with zero % coverage from
    # predictions.
    R_Hab_Mask <- terra::classify(R_Hab, cbind(0, NA), others = 1)

    R_Hab <- stats::setNames(log10(R_Hab + 0.1), "HabLog")

    StaticPredictors <- c(StaticPredictors, R_Hab)
    rm(R_Hab, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Sampling efforts -----

  if ("EffortsLog" %in% other_variables) {

    ecokit::cat_time("Sampling efforts", level = 1L)

    R_Efforts <- fs::path(Path_Bias, "Efforts_SummaryR.RData")
    if (!fs::file_exists(R_Efforts)) {
      ecokit::stop_ctx(
        "Sampling efforts data does not exist", R_Efforts = R_Efforts,
        include_backtrace = TRUE)
    }

    # add 1 (older versions 0.1) to get log for 0 values
    # [only for rivers/roads/efforts, not hab/rivers]
    R_Efforts <- ecokit::load_as(R_Efforts, unwrap_r = TRUE) %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("EffortsLog")

    if (clamp_pred) {

      ecokit::cat_time("Fixing sampling efforts values", level = 2L)

      # Check fix_efforts value
      if (is.numeric(fix_efforts)) {

        # If `fix_efforts` is numeric value, check if it is within the range of
        # the observed efforts
        EffortsRange <- terra::global(R_Efforts, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_efforts, EffortsRange[1], EffortsRange[2]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_efforts` value is out of the range of observed efforts",
            fix_efforts = fix_efforts, EffortsRange = round(EffortsRange, 2),
            include_backtrace = TRUE)
        }

        # Fix value
        EffortsVal <- fix_efforts

      } else {

        # If `fix_efforts` is character, check if it is one of the valid values:
        # identity, median, mean, max, and q90
        fix_efforts <- stringr::str_to_lower(fix_efforts)
        if (!(fix_efforts %in% c("identity", "median", "mean", "max", "q90"))) {
          ecokit::stop_ctx(
            paste0(
              "`fix_efforts` has to be either NULL, single numeric ",
              "value, or one of the following: 'identity', 'median', ",
              "'mean', 'max', or 'q90'."),
            fix_efforts = fix_efforts, include_backtrace = TRUE)
        }
      }

      # Fix value
      EffortsVal <- dplyr::case_when(

        # Do not fix if `fix_efforts` is "identity"
        fix_efforts == "identity" ~ NA_real_,

        # Fix at 90% quantile
        fix_efforts == "q90" ~ {
          terra::global(
            R_Efforts,
            fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at median value
        fix_efforts == "median" ~ {
          terra::global(
            R_Efforts, fun = function(x) median(x, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at mean value
        fix_efforts == "mean" ~ {
          terra::global(R_Efforts, fun = mean, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at max value
        fix_efforts == "max" ~ {
          terra::global(R_Efforts, fun = max, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        .default = NA_real_)


      # Print fixed value
      if (is.na(EffortsVal)) {
        if (fix_efforts == "identity") {
          ecokit::cat_time(
            "Sampling efforts predictor is not fixed",
            level = 2L, cat_timestamp = FALSE)
        }
        R_Efforts_Clamp <- stats::setNames(R_Efforts, "EffortsLog_Clamp")

      } else {

        ecokit::cat_time(
          paste0("Fixed value is ", round(EffortsVal, 2), " [log10 scale]"),
          level = 2L, cat_timestamp = FALSE)

        # Set a minimum value for efforts variable to `EffortsVal`. Using
        # upper = Inf keeps efforts values > EffortsVal as they are.
        R_Efforts_Clamp <- terra::clamp(
          x = R_Efforts, lower = EffortsVal, upper = Inf) %>%
          stats::setNames("EffortsLog_Clamp")

      }

      StaticPredictors <- c(StaticPredictors, R_Efforts, R_Efforts_Clamp)
      rm(R_Efforts, R_Efforts_Clamp, envir = environment())

    } else {

      # Do not fix at single value
      ecokit::cat_time(
        "Sampling efforts is not fixed", level = 2L, cat_timestamp = FALSE)
      StaticPredictors <- c(StaticPredictors, R_Efforts)
      rm(R_Efforts, envir = environment())

    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## River length ----

  if ("RiversLog" %in% other_variables) {

    ecokit::cat_time("River length", level = 1L)

    R_Rivers <- fs::path(Path_Rivers, "River_Lengths.RData")
    if (!fs::file_exists(R_Rivers)) {
      ecokit::stop_ctx(
        "River length data does not exist", R_Rivers = R_Rivers,
        include_backtrace = TRUE)
    }

    R_Rivers <- ecokit::load_as(R_Rivers, unwrap_r = TRUE) %>%
      magrittr::extract2("STRAHLER_5") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("RiversLog")

    if (clamp_pred) {
      ecokit::cat_time("Fixing river length values", level = 2L)

      # Check fix_rivers value
      if (is.numeric(fix_rivers)) {

        # If `fix_rivers` is numeric value, check if it is within the range of
        # the observed river lengths
        RiversRange <- terra::global(R_Rivers, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_rivers, RiversRange[1], RiversRange[2]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_rivers` value is out of the range of observed river length",
            fix_rivers = fix_rivers, RiversRange = round(RiversRange, 2),
            include_backtrace = TRUE)
        }

        # Fix value
        RiversVal <- fix_rivers

      } else {

        # If `fix_rivers` is character, check if it is one of the valid values:
        # identity, median, mean, max, and q90
        fix_rivers <- stringr::str_to_lower(fix_rivers)
        if (!(fix_rivers %in% c("identity", "median", "mean", "max", "q90"))) {
          ecokit::stop_ctx(
            paste0(
              "`fix_rivers` has to be either NULL, single numeric ",
              "value, or one of the following: 'identity', 'median', ",
              "'mean', 'max', or 'q90'."),
            fix_rivers = fix_rivers, include_backtrace = TRUE)
        }

        # Fix value
        RiversVal <- dplyr::case_when(

          # Do not fix if `fix_rivers` is "identity"
          fix_rivers == "identity" ~ NA_real_,

          # Fix at 90% quantile
          fix_rivers == "q90" ~ {
            terra::global(
              R_Rivers,
              fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at median value
          fix_rivers == "median" ~ {
            terra::global(
              R_Rivers, fun = function(x) median(x, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at mean value
          fix_rivers == "mean" ~ {
            terra::global(R_Rivers, fun = mean, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at max value
          fix_rivers == "max" ~ {
            terra::global(R_Rivers, fun = max, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          .default = NA_real_)
      }

      # Print fixed value

      if (is.na(RiversVal)) {

        if (fix_rivers == "identity") {
          ecokit::cat_time(
            "River length predictor is not fixed at a single value",
            level = 2L, cat_timestamp = FALSE)
        }
        R_Rivers_Clamp <- stats::setNames(R_Rivers, "RiversLog_Clamp")

      } else {

        ecokit::cat_time(
          paste0("Fixed value is ", round(RiversVal, 2), " [log10 scale]"),
          level = 2L, cat_timestamp = FALSE)

        # Set a minimum value for river length variable to `RiversVal`. Using
        # upper = Inf keeps  river length values > RiversVal as they are.
        R_Rivers_Clamp <- terra::clamp(
          x = R_Rivers, lower = RiversVal, upper = Inf) %>%
          stats::setNames("RiversLog_Clamp")

      }

      StaticPredictors <- c(StaticPredictors, R_Rivers, R_Rivers_Clamp)
      rm(R_Rivers, R_Rivers_Clamp, RiversVal, envir = environment())

    } else {

      # Do not fix at single value
      ecokit::cat_time(
        "River length predictor is not fixed at a single value",
        level = 2L, cat_timestamp = FALSE)

      StaticPredictors <- c(StaticPredictors, R_Rivers)
      rm(R_Rivers, envir = environment())

    }

  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Soil bulk density -----

  if ("soil" %in% other_variables) {

    ecokit::cat_time("Soil bulk density", level = 1L)

    R_Soil <- fs::path(Path_Soil, "soil_density.RData")
    if (!fs::file_exists(R_Soil)) {
      ecokit::stop_ctx(
        "Soil bulk density data does not exist", R_Soil = R_Soil,
        include_backtrace = TRUE)
    }

    R_Soil <- ecokit::load_as(R_Soil, unwrap_r = TRUE) %>%
      magrittr::extract2("bdod_5_15_mean") %>%
      stats::setNames("soil")

    StaticPredictors <- c(StaticPredictors, R_Soil)
    rm(R_Soil, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Topographic wetness index -----

  if ("wetness" %in% other_variables) {

    ecokit::cat_time("Topographic wetness index", level = 1L)

    R_wetness <- fs::path(Path_Wetness, "wetness_index.RData")
    if (!fs::file_exists(R_wetness)) {
      ecokit::stop_ctx(
        "Wetness data does not exist", R_wetness = R_wetness,
        include_backtrace = TRUE)
    }
    R_wetness <- ecokit::load_as(R_wetness, unwrap_r = TRUE) %>%
      stats::setNames("wetness")

    StaticPredictors <- c(StaticPredictors, R_wetness)
    rm(R_wetness, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  ecokit::cat_time("Merge static predictors", level = 1L)

  StaticPredictors <- terra::rast(StaticPredictors)

  # If Habitat predictor is used, grid cells with zero % coverage are
  # excluded from predictions
  if (Hab_Predictor) {
    StaticPredictors <- terra::mask(StaticPredictors, R_Hab_Mask)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict latent factor at new locations ------

  Path_Test_LF <- fs::path(Path_Prediction1, "Test_LF.qs2")

  if (!fs::file_exists(Path_Test_LF) && pred_new_sites && spatial_model) {

    ecokit::cat_time("Predict latent factor at new locations")

    ecokit::cat_time(
      "Preparing input data for predicting latent factor", level = 1L)

    Predict_DF_Test <- Prediction_Options %>%
      dplyr::filter(ClimateModel == "Current") %>%
      dplyr::pull("FilePath") %>%
      # If clamp_pred`=`TRUE`, there are two options for Current climate data
      # (with and without clamping). Two sets of predictions under current
      # climates will be produced. Predictions without clamping is used for
      # model evaluation.
      utils::head(1) %>%
      ecokit::load_as(unwrap_r = TRUE) %>%
      # Only extract predictors used in the model
      terra::subset(bio_variables) %>%
      # Combine with other static predictors
      c(StaticPredictors) %>% # nolint: consecutive_concatenation_linter
      # If Habitat predictor is used, grid cells with zero % coverage are
      # excluded from predictions [na.rm = TRUE]
      terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
      tibble::tibble() %>%
      sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
      sf::st_join(Model_Coords) %>%
      tidyr::replace_na(list(Train = FALSE)) %>%
      dplyr::filter(!Train)

    Test_XY <- sf::st_drop_geometry(Predict_DF_Test[, c("x", "y")])
    Test_X <- Predict_DF_Test %>%
      dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
      sf::st_drop_geometry()
    Gradient <- Hmsc::prepareGradient(
      hM = Model, XDataNew = as.data.frame(Test_X),
      sDataNew = list(sample = as.data.frame(Test_XY)))

    rm(Predict_DF_Test, Test_X, Test_XY, Model, envir = environment())
    invisible(gc())

    ecokit::cat_time("Predicting latent factor", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2,
      n_separators = 1L, line_char = "*")

    # Predicting latent factor only
    Preds_LF <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy,
      model_name = paste0("LF_", hab_abb, "_Test"), temp_dir = temp_dir,
      temp_cleanup = temp_cleanup, use_TF = use_TF, TF_environ = TF_environ,
      LF_out_file = Path_Test_LF, TF_use_single = TF_use_single,
      LF_only = TRUE, LF_n_cores = LF_n_cores, LF_check = LF_check,
      LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = LF_commands_only,
      evaluate = FALSE, verbose = TRUE, spatial_model = spatial_model)

    rm(Gradient, Preds_LF, envir = environment())

    ecokit::cat_time("Predicting latent factor is finished!", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2,
      n_separators = 1L, line_char = "*")

    if (LF_commands_only) {
      return(invisible(NULL))
    }

  } else {

    if (pred_new_sites & spatial_model) {
      ecokit::cat_time(
        "LF prediction is already available on disk", level = 1L)
    } else {
      ecokit::cat_time("LF prediction will NOT be made", level = 1L)
    }

    rm(Model, envir = environment())
    invisible(gc())
  }

  if (LF_only) {
    return(invisible(NULL))
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict_Internal ------

  Predict_Internal <- function(ID) {

    # ID: Index of the current prediction option

    # Whether to clamp the sampling efforts
    DoClamp <- Prediction_Options$Clamp[[ID]]

    # Name of the current option
    Option_Name <- Prediction_Options$Name[[ID]]

    # Name of the current model
    model_name <- paste0(
      Option_Name, " - ", dplyr::if_else(DoClamp, "clamping", "no clamping"))

    MSG <- paste0(
      model_name, " (", ID, "/", nrow(Prediction_Options), ")")
    ecokit::info_chunk(
      MSG, n_separators = 1L, line_char = "-", line_char_rep = 70L,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1L,
      info_lines_before = 1L)

    StaticPreds <- StaticPredictors

    if (DoClamp) {

      # Do not evaluate for options with clamping
      evaluate <- FALSE

      # Make prediction files at `Path_Prediction_Clamp`
      Path_Prediction <- Path_Prediction_Clamp

      if ("EffortsLog" %in% names(StaticPredictors)) {
        # use clamped Effort values
        StaticPreds$EffortsLog <- StaticPreds$EffortsLog_Clamp
        StaticPreds$EffortsLog_Clamp <- NULL
      }

      if ("RiversLog" %in% names(StaticPredictors)) {
        # use clamped rivers values
        StaticPreds$RiversLog <- StaticPreds$RiversLog_Clamp
        StaticPreds$RiversLog_Clamp <- NULL
      }

    } else {

      # evaluate for "Current" climates, without clamping
      evaluate <- (Option_Name == "Current")

      # Make prediction files at `Path_Prediction_NoClamp`
      Path_Prediction <- Path_Prediction_NoClamp

      # Remove clamped layers
      if ("EffortsLog_Clamp" %in% names(StaticPreds)) {
        StaticPreds <- terra::subset(
          x = StaticPreds, subset = "EffortsLog_Clamp", negate = TRUE)
      }
      if ("RiversLog_Clamp" %in% names(StaticPreds)) {
        StaticPreds <- terra::subset(
          x = StaticPreds, subset = "RiversLog_Clamp", negate = TRUE)
      }

    }

    Path_Prediction_sf <- fs::path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_sf.qs2"))
    Path_Prediction_R <- fs::path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_R.qs2"))
    Path_Prediction_summary <- fs::path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_Summary.RData"))

    # Path for saving tif files of the current option
    Path_Prediction_tif <- fs::path(Path_Prediction, Option_Name)
    fs::dir_create(Path_Prediction_tif)
    invisible(gc())

    # ______________________________________________
    # ______________________________________________

    # Making predictions if not already processed ----
    OutMissing <- fs::file_exists(
      c(Path_Prediction_R, Path_Prediction_sf, Path_Prediction_summary)) %>%
      all() %>%
      isFALSE()

    if (OutMissing) {

      # Load model
      Model <- ecokit::load_as(path_model)

      # Skip predictions if the predictions as sf object is already on disk
      if (fs::file_exists(Path_Prediction_sf)) {

        ecokit::cat_time("Loading predictions `sf` from disk", level = 1L)
        Prediction_sf <- ecokit::load_as(Path_Prediction_sf)

      } else {

        .OptionStartTime <- lubridate::now(tzone = "CET")

        # Extracting data at training and new sites ------
        ecokit::cat_time("Extracting data at training and new sites")
        Predict_DF <- Prediction_Options$FilePath[[ID]] %>%
          ecokit::load_as(unwrap_r = TRUE) %>%
          terra::subset(bio_variables) %>%
          c(StaticPreds) %>%
          terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
          tibble::tibble() %>%
          sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
          sf::st_join(Model_Coords) %>%
          tidyr::replace_na(list(Train = FALSE))

        # Training locations
        Model_Name_Train <- paste0(
          Option_Name, "_",
          dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Train")
        Predict_DF_Train <- dplyr::filter(Predict_DF, Train)
        if (nrow(Predict_DF_Train) > 0L) {
          Train_XY <- sf::st_drop_geometry(Predict_DF_Train[, c("x", "y")])
          Train_PA <- as.data.frame(Model$Y)
          Train_X <- Predict_DF_Train %>%
            dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
            sf::st_drop_geometry() %>%
            stats::model.matrix(Model$XFormula, ., xlev = NULL)
        } else {
          Train_XY <- Train_PA <- Train_X <- NULL
        }

        # Testing locations
        Model_Name_Test <- paste0(
          Option_Name, "_",
          dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Test")
        Predict_DF_Test <- dplyr::filter(Predict_DF, !Train)
        Test_XY <- Predict_DF_Test[, c("x", "y")]
        Test_X <- Predict_DF_Test %>%
          dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
          sf::st_drop_geometry()
        Gradient <- Hmsc::prepareGradient(
          hM = Model, XDataNew = as.data.frame(Test_X),
          sDataNew = list(
            sample = as.data.frame(sf::st_drop_geometry(Test_XY))))

        rm(Model, Predict_DF_Test, Predict_DF, envir = environment())
        invisible(gc())

        # ______________________________________________

        if (nrow(Predict_DF_Train) > 0L) {

          # Predicting at training sites ----
          ecokit::cat_time("Predicting at training sites")

          Path_Current_Train <- fs::path(
            Path_Prediction, paste0("Prediction_", Option_Name, "_Train.qs2"))

          if (fs::file_exists(Path_Current_Train)) {

            ecokit::cat_time("Loading predictions from disk", level = 2L)
            Preds_ModFitSites <- tibble::tibble(Pred_Path = Path_Current_Train)

          } else {

            Preds_ModFitSites <- IASDT.R::predict_hmsc(
              path_model = path_model, X = Train_X, gradient = NULL,
              expected = TRUE, n_cores = n_cores, strategy = strategy,
              model_name = Model_Name_Train, temp_dir = temp_dir,
              temp_cleanup = temp_cleanup, use_TF = use_TF,
              TF_environ = TF_environ, TF_use_single = TF_use_single,
              LF_return = TRUE, LF_n_cores = LF_n_cores, LF_check = LF_check,
              LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
              pred_directory = Path_Prediction, pred_PA = Train_PA,
              pred_XY = Train_XY, evaluate = evaluate, evaluation_name = NULL,
              evaluation_directory = Path_Eval, verbose = FALSE,
              spatial_model = spatial_model)

          }
        } else {
          Preds_ModFitSites <- tibble::tibble(Pred_Path = NULL)
        }

        # ______________________________________________

        # Predicting at new sites ----

        if (pred_new_sites) {

          ecokit::cat_time("Predicting at new sites")

          Path_Current_Test <- fs::path(
            Path_Prediction, paste0("Prediction_", Option_Name, "_Test.qs2"))

          if (fs::file_exists(Path_Current_Test)) {

            ecokit::cat_time("Loading predictions from disk", level = 2L)
            Preds_NewSites <- tibble::tibble(Pred_Path = Path_Current_Test)

          } else {

            Preds_NewSites <- IASDT.R::predict_hmsc(
              path_model = path_model, gradient = Gradient, expected = TRUE,
              n_cores = n_cores, strategy = strategy,
              model_name = Model_Name_Test, temp_dir = temp_dir,
              temp_cleanup = temp_cleanup, use_TF = use_TF,
              TF_environ = TF_environ, TF_use_single = TF_use_single,
              LF_return = TRUE, LF_inputFile = Path_Test_LF,
              LF_n_cores = LF_n_cores, LF_check = LF_check,
              LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
              verbose = FALSE, pred_directory = Path_Prediction,
              evaluate = FALSE, pred_XY = sf::st_drop_geometry(Test_XY),
              spatial_model = spatial_model)

          }
        } else {
          ecokit::cat_time(
            "Predictions at new sites will NOT be made", level = 1L)
          Preds_NewSites <- tibble::tibble(Pred_Path = NULL)
        }

        # ______________________________________________

        # Merge & save predictions - sf ------
        ecokit::cat_time("Merge & save predictions at training and new sites")
        Prediction_sf <- purrr::map_dfr(
          .x = c(Preds_ModFitSites$Pred_Path, Preds_NewSites$Pred_Path),
          .f = ~ {
            if (fs::file_exists(.x)) {
              ecokit::load_as(.x)
            } else {
              tibble::tibble()
            }
          })

        # ______________________________________________

        # Save predictions as sf object
        ecokit::save_as(object = Prediction_sf, out_path = Path_Prediction_sf)

        ecokit::cat_diff(
          init_time = .OptionStartTime, prefix = "Prediction took ", level = 1L)

      }

      # ______________________________________________
      # ______________________________________________

      ### Predictions as spatRaster / tif -----

      ecokit::cat_time("Rasterization & prepare summary data", level = 1L)

      Fields2Raster <- names(Prediction_sf) %>%
        stringr::str_subset("^Sp_|^SR_") %>%
        gtools::mixedsort()

      Grid10 <- ecokit::load_as(Path_GridR, unwrap_r = TRUE)
      Prediction_R <- terra::rasterize(
        Prediction_sf, Grid10, field = Fields2Raster)

      # Calculate prediction anomaly for future projections
      if (Option_Name != "Current") {

        # Names of current taxa (and SR) for the current model
        MeanNames <- stringr::str_subset(Fields2Raster, "_mean")

        # Loading mean predictions at current climates
        CurrentMean <- list.files(
          # use relevant folder containing the current predictions. This is
          # determined by `Path_Prediction`, which is not the same whether
          # clamping is used or not
          path = Path_Prediction, pattern = "Prediction_Current.*_R.qs2",
          full.names = TRUE) %>%
          ecokit::load_as(unwrap_r = TRUE) %>%
          terra::subset(MeanNames)

        # Calculate anomaly as difference in predicted value between future and
        # current climate (future - current)
        Preds_Anomaly <- terra::subset(Prediction_R, MeanNames) - CurrentMean
        # Assign names to anomaly maps
        AnomalyNames <- names(Preds_Anomaly) %>%
          stringr::str_replace_all("_mean", "_anomaly")
        names(Preds_Anomaly) <- AnomalyNames

        # Add anomaly maps to the list of predictions
        Fields2Raster <- c(Fields2Raster, AnomalyNames)
        Prediction_R <- c(Prediction_R, Preds_Anomaly)

        # clean up
        rm(Preds_Anomaly, CurrentMean, envir = environment())

      }

      Out_Summary <- tibble::tibble(
        hab_abb = hab_abb, hab_name = hab_name,
        layer_name = Fields2Raster,
        time_period = Prediction_Options$TimePeriod[[ID]],
        climate_model = Prediction_Options$ClimateModel[[ID]],
        climate_scenario = Prediction_Options$ClimateScenario[[ID]],
        Clamp = Prediction_Options$Clamp[[ID]],
        Path_Prediction = Path_Prediction) %>%
        dplyr::mutate(
          Stats = dplyr::case_when(
            endsWith(layer_name, "_mean") ~ "tif_path_mean",
            endsWith(layer_name, "_sd") ~ "tif_path_sd",
            endsWith(layer_name, "_cov") ~ "tif_path_cov",
            endsWith(layer_name, "_anomaly") ~ "tif_path_anomaly",
            .default = NULL),
          ias_id = stringr::str_remove(
            layer_name, "_mean$|_sd$|_cov$|_anomaly"),
          tif_path = fs::path(
            Path_Prediction_tif, paste0(layer_name, ".tif")))

      # Save as tif
      ecokit::cat_time("Save as tif", level = 1L)
      Out_Summary0 <- Out_Summary %>%
        dplyr::mutate(
          Map = purrr::map2(
            .x = layer_name, .y = tif_path,
            .f = ~ {
              terra::writeRaster(
                x = Prediction_R[[.x]], filename = .y, overwrite = TRUE,
                gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
            }))

      rm(Out_Summary0, envir = environment())

      Out_Summary <- Out_Summary %>%
        tidyr::pivot_wider(
          id_cols = c(
            "hab_abb", "hab_name", "time_period", "climate_model",
            "climate_scenario", "ias_id", "Clamp", "Path_Prediction"),
          names_from = "Stats", values_from = "tif_path") %>%
        dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
        dplyr::select(
          tidyselect::any_of(
            c(
              "hab_abb", "hab_name", "time_period", "climate_model",
              "climate_scenario", "ias_id", "taxon_name", "species_name",
              "class", "order", "family", "tif_path_mean", "tif_path_sd",
              "tif_path_cov", "tif_path_anomaly")))

      # save as spatRaster - qs2
      ecokit::cat_time("Save as spatRaster - qs2", level = 1L)
      Prediction_R <- terra::wrap(Prediction_R)
      ecokit::save_as(object = Prediction_R, out_path = Path_Prediction_R)

      # Save summary - RData
      ecokit::cat_time("Save summary - RData", level = 1L)
      ecokit::save_as(
        object = Out_Summary, object_name = paste0(Option_Name, "_Summary"),
        out_path = Path_Prediction_summary)

      # Save summary - csv
      ecokit::cat_time("Save summary - csv", level = 1L)
      utils::write.table(
        x = Out_Summary,
        file = stringr::str_replace(Path_Prediction_summary, ".RData$", ".txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = FALSE, fileEncoding = "UTF-8")
    }

    # output
    return(
      tibble::tibble(
        Clamp = DoClamp,
        Name = Option_Name,
        File_Pred_R = Path_Prediction_R, File_Pred_sf = Path_Prediction_sf,
        File_Pred_summary = Path_Prediction_summary))
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predicting ------

  ecokit::info_chunk(
    "Making spatial predictions", n_separators = 2L, level = 1L,
    line_char = "*", line_char_rep = 70L, cat_red = TRUE,
    cat_bold = TRUE, cat_timestamp = FALSE)

  Grid10 <- ecokit::load_as(Path_GridR, unwrap_r = TRUE)

  Prediction_Summary <- purrr::map_dfr(
    .x = seq_len(nrow(Prediction_Options)), .f = Predict_Internal) %>%
    dplyr::full_join(Prediction_Options, ., by = c("Name", "Clamp")) %>%
    dplyr::select(-"FilePath")

  rm(Predict_Internal, Grid10, Model_Coords, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Ensemble model predictions ------

  ecokit::info_chunk(
    "\tEnsemble model predictions", n_separators = 1L, line_char = "-",
    line_char_rep = 70L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  # Prepare input data to calculate ensemble predictions
  ecokit::cat_time(
    "Prepare input data to calculate ensemble predictions", level = 1L)

  Prediction_Ensemble <- Prediction_Summary %>%
    dplyr::filter(ClimateModel != "Current") %>%
    dplyr::select(-File_Pred_sf, -File_Pred_R, -Name, -ClimateModel) %>%
    dplyr::mutate(
      Prediction2 = purrr::map(
        .x = File_Pred_summary,
        .f = ~ {
          if (!fs::file_exists(.x)) {
            warning("File not found: ", .x, call. = FALSE)
            return(NULL)
          }

          ecokit::load_as(.x) %>%
            dplyr::select(
              -tidyselect::all_of(
                c("tif_path_sd", "tif_path_cov", "tif_path_anomaly")))

        })) %>%
    dplyr::select("Prediction2", "Clamp") %>%
    tidyr::unnest("Prediction2") %>%
    dplyr::mutate(
      climate_model = "Ensemble",
      Dir_Ensemble = fs::path(
        dirname(dirname(tif_path_mean)),
        paste0(
          stringr::str_replace(time_period, "-", "_"), "_", climate_scenario,
          "_Ensemble"))) %>%
    dplyr::group_by(dplyr::across(-tif_path_mean)) %>%
    dplyr::summarise(tifs = list(tif_path_mean), .groups = "drop") %>%
    dplyr::mutate(
      tif_path_mean = fs::path(Dir_Ensemble, paste0(ias_id, "_mean.tif")),
      tif_path_anomaly = fs::path(
        Dir_Ensemble, paste0(ias_id, "_anomaly.tif")),
      tif_path_sd = fs::path(Dir_Ensemble, paste0(ias_id, "_sd.tif")),
      tif_path_cov = fs::path(Dir_Ensemble, paste0(ias_id, "_cov.tif")))

  # --------------------------------------------------------- #

  ecokit::cat_time("Create directories for ensemble predictions", level = 1L)
  fs::dir_create(unique(Prediction_Ensemble$Dir_Ensemble))

  # --------------------------------------------------------- #

  # Loading mean predictions at current climates
  ecokit::cat_time("Loading mean predictions at current climates", level = 1L)

  CurrentMean <- list.files(
    path = dplyr::if_else(
      clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    pattern = "Prediction_Current.*_R.qs2", full.names = TRUE)

  # --------------------------------------------------------- #

  # Calculate ensemble predictions
  ecokit::cat_time("Calculate ensemble predictions", level = 1L)

  # set up parallel processing
  doParallel::registerDoParallel(cores = n_cores)
  ecokit::load_packages(package_list = "foreach")
  withr::defer(doParallel::stopImplicitCluster())

  Prediction_Ensemble_0 <- dplyr::select(
    Prediction_Ensemble,
    tifs, tif_path_mean, tif_path_anomaly, tif_path_sd, tif_path_cov, ias_id)

  calculate_ensemble <- foreach::foreach(
    row_id = seq_len(nrow(Prediction_Ensemble_0)),
    .export = c("CurrentMean", "Prediction_Ensemble_0"),
    .packages = c("terra", "purrr", "ecokit", "qs2", "magrittr", "fs")
  ) %dopar% { # nolint: object_usage_linter

    path_mean <- Prediction_Ensemble_0$tif_path_mean[[row_id]]
    path_anomaly <- Prediction_Ensemble_0$tif_path_anomaly[[row_id]]
    path_cov <- Prediction_Ensemble_0$tif_path_cov[[row_id]]
    path_sd <- Prediction_Ensemble_0$tif_path_sd[[row_id]]

    paths_all <- c(path_mean, path_anomaly, path_cov, path_sd)
    paths_exist <- purrr::map_lgl(.x = paths_all, .f = fs::file_exists)
    all_okay <- all(paths_exist)
    if (all_okay) {
      all_okay <- purrr::map_lgl(
        .x = paths_all, .f = ecokit::check_tiff, warning = FALSE) %>%
        all()
    }

    if (all_okay) {
      return(NULL)
    }

    fs::file_delete(paths_all[paths_exist])

    ias_id <- Prediction_Ensemble_0$ias_id[[row_id]]

    # load maps for future climate option
    tiffs_R <- terra::rast(Prediction_Ensemble_0$tifs[[row_id]])

    # Mean
    Ensemble_mean <- terra::app(tiffs_R, "mean", na.rm = TRUE)
    # Standard deviation
    Ensemble_sd <- terra::app(tiffs_R, "sd", na.rm = TRUE)
    rm(tiffs_R)

    # Anomaly
    CurrentMean0 <- ecokit::load_as(CurrentMean, unwrap_r = TRUE) %>%
      terra::subset(paste0(ias_id, "_mean"))
    Ensemble_anomaly <- Ensemble_mean - CurrentMean0
    rm(CurrentMean0)

    # coefficient of variation: Replace very small mean values with
    # reasonably small number to avoid overflow warning
    Ensemble_mean_0 <- terra::classify(
      x = Ensemble_mean, rcl = cbind(0, 1e-8, 1e-8))
    Ensemble_cov <- (Ensemble_sd / Ensemble_mean_0)

    gdal_o <- c("COMPRESS=DEFLATE", "TILED=YES")
    terra::writeRaster(
      x = Ensemble_mean, filename = path_mean,
      overwrite = TRUE, gdal = gdal_o)
    terra::writeRaster(
      x = Ensemble_anomaly, filename = path_anomaly,
      overwrite = TRUE, gdal = gdal_o)
    terra::writeRaster(
      x = Ensemble_cov, filename = path_cov,
      overwrite = TRUE, gdal = gdal_o)
    terra::writeRaster(
      x = Ensemble_sd, filename = path_sd,
      overwrite = TRUE, gdal = gdal_o)

    return(NULL)
  }

  rm(
    CurrentMean, Prediction_Ensemble_0, calculate_ensemble,
    envir = environment())
  invisible(gc())

  # --------------------------------------------------------- #

  # Save ensemble maps as SpatRast
  ecokit::cat_time("Save ensemble predictions as SpatRast", level = 1L)

  vars_to_select <- c(
    "ias_id", "Ensemble_File", "tif_path_mean", "tif_path_anomaly",
    "tif_path_sd", "tif_path_cov")
  Prediction_Ensemble_R <- Prediction_Ensemble %>%
    dplyr::select(-Dir_Ensemble, -tifs) %>%
    dplyr::mutate(
      Ensemble_File = fs::path(
        dplyr::if_else(
          clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
        paste0(
          "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
          climate_scenario, "_Ensemble_R.qs2"))) %>%
    dplyr::select(tidyselect::all_of(vars_to_select)) %>%
    tidyr::nest(data = -Ensemble_File)

  Prediction_Ensemble_R <- foreach::foreach(
    row_id = seq_len(nrow(Prediction_Ensemble_R)),
    .export = "Prediction_Ensemble_R",
    .packages = c("dplyr", "magrittr", "terra", "purrr", "ecokit", "fs")
  ) %dopar% { # nolint: object_usage_linter

    OutMaps <- dplyr::slice(Prediction_Ensemble_R, row_id) %>%
      dplyr::pull("data") %>%
      magrittr::extract2(1) %>%
      dplyr::mutate(
        Maps = purrr::pmap(
          .l = list(ias_id, tif_path_mean, tif_path_anomaly,
                    tif_path_sd, tif_path_cov),
          .f = function(ias_id, tif_path_mean, tif_path_anomaly,
                        tif_path_sd, tif_path_cov) {

            Mean <- terra::rast(tif_path_mean) %>%
              stats::setNames(paste0(ias_id, "_mean"))
            SD <- terra::rast(tif_path_sd) %>%
              stats::setNames(paste0(ias_id, "_sd"))
            COV <- terra::rast(tif_path_cov) %>%
              stats::setNames(paste0(ias_id, "_cov"))
            Anomaly <- terra::rast(tif_path_anomaly) %>%
              stats::setNames(paste0(ias_id, "_anomaly"))

            c(Mean, SD, COV, Anomaly) %>%
              terra::wrap()
          }
        )) %>%
      dplyr::pull("Maps") %>%
      purrr::map(terra::unwrap) %>%
      terra::rast() %>%
      terra::wrap()

    ecokit::save_as(
      object = OutMaps,
      out_path = Prediction_Ensemble_R$Ensemble_File[[row_id]])
    return(NULL)
  }

  doParallel::stopImplicitCluster()

  rm(Prediction_Ensemble_R, envir = environment())

  # --------------------------------------------------------- #

  # Save summary of ensemble predictions
  ecokit::cat_time("Save summary of ensemble predictions", level = 1L)

  Prediction_Ensemble_Summary <- Prediction_Ensemble %>%
    dplyr::select(-Ensemble_Maps) %>%
    dplyr::mutate(
      Ensemble_File = fs::path(
        dplyr::if_else(
          clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
        paste0(
          "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
          climate_scenario, "_Ensemble_Summary.RData")),
      # utils::write.table does not support `fs::path()`
      # convert paths to character and use readr::write_delim instead
      dplyr::across(
        .cols = tidyselect::all_of(
          c("Dir_Ensemble", "tifs", "tif_path_mean",
            "tif_path_anomaly", "tif_path_sd", "tif_path_cov")),
        .fns = as.character)) %>%
    tidyr::nest(Ensemble_DT = -Ensemble_File) %>%
    dplyr::mutate(
      Ensemble_Save = purrr::map2(
        .x = Ensemble_DT, .y = Ensemble_File,
        .f = ~ {
          ecokit::save_as(
            object = .x, out_path = .y,
            object_name = stringr::str_remove(basename(.y), ".RData"))
          readr::write_delim(
            x = .x, file = stringr::str_replace(.y, ".RData", ".txt"),
            delim = "\t", col_names = TRUE,
            append = FALSE, progress = FALSE)
        }),
      Ensemble_Save = NULL) %>%
    tidyr::unnest(Ensemble_DT) %>%
    dplyr::select(
      tidyselect::all_of(
        c("Ensemble_File", "time_period", "climate_model", "climate_scenario"))
    ) %>%
    dplyr::rename(File_Pred_summary = Ensemble_File) %>%
    dplyr::mutate(
      Clamp = clamp_pred,
      Name = paste0(
        stringr::str_replace(time_period, "-", "_"), "_",
        climate_scenario, "_Ensemble"),
      File_Pred_sf = NA_character_,
      File_Pred_R = stringr::str_replace(
        File_Pred_summary, "_Summary.RData", "_R.qs2")) %>%
    dplyr::distinct()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Overall summary -----

  ecokit::info_chunk(
    "\tPrepare overall summary", n_separators = 1L, line_char = "-",
    line_char_rep = 70L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  Prediction_Summary <- Prediction_Summary %>%
    dplyr::rename(
      time_period = TimePeriod, climate_model = ClimateModel,
      climate_scenario = ClimateScenario) %>%
    dplyr::bind_rows(Prediction_Ensemble_Summary) %>%
    dplyr::mutate(hab_abb = hab_abb, hab_name = hab_name, .before = 1)

  readr::write_delim(
    x = Prediction_Summary, file = Path_Summary_txt, delim = "\t",
    col_names = TRUE, append = FALSE, progress = FALSE)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Overall summary - to be uploaded to the data server; for the Shiny App -----

  if (clamp_pred) {
    Prediction_Summary_Shiny <- dplyr::filter(Prediction_Summary, Clamp)
  } else {
    Prediction_Summary_Shiny <- dplyr::filter(Prediction_Summary, !Clamp)
  }

  Prediction_Summary_Shiny <- Prediction_Summary_Shiny$File_Pred_summary %>%
    purrr::map(ecokit::load_as) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::all_of(
          c("tif_path_mean", "tif_path_sd",
            "tif_path_cov", "tif_path_anomaly")),
        .fns = ~ {
          stringr::str_remove(
            string = .x,
            pattern = paste0(dirname(Path_Summary_RData_Shiny), "/"))
        }))

  save(Prediction_Summary_Shiny, file = Path_Summary_RData_Shiny)
  utils::write.table(
    x = Prediction_Summary_Shiny, sep = "\t", row.names = FALSE,
    col.names = TRUE, file = Path_Summary_txt_Shiny, quote = FALSE,
    fileEncoding = "UTF-8")

  # Create tar file for prediction files
  if (tar_predictions) {

    ecokit::cat_time("Create tar file for prediction files", level = 1L)

    # Directory to save the tar file
    TarDir <- dirname(Path_Summary_RData_Shiny)
    # Path to the tar file
    TarFile <- fs::path(TarDir, "Predictions.tar")
    # List of directories in the prediction folder. All directories will be
    # included in the tar file
    TarFiles <- list.dirs(      # nolint: object_name_linter
      path = TarDir, full.names = FALSE, recursive = FALSE) %>%
      paste(collapse = " ") %>%
      # Add the summary files to the list
      paste(
        "Prediction_Summary_Shiny.RData", "Prediction_Summary_Shiny.txt",
        collapse = " ")

    # Command to create the tar file
    Command <- stringr::str_glue(
      "cd {fs::path_abs(TarDir)}; tar -cf {basename(TarFile)} -b 2048 \\
      {TarFiles}")

    # Create tar file
    system(Command)

    # Change the permission of the tar file
    Sys.chmod(TarFile, "755", use_umask = FALSE)
  }

  # # ................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nThe whole prediction function took ")

  return(Prediction_Summary)
}
