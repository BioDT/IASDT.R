## |------------------------------------------------------------------------| #
# predict_maps ----
## |------------------------------------------------------------------------| #

#' Predict habitat suitability of `Hmsc` models
#'
#' This package provides two functions for predicting habitat suitability of
#' `Hmsc` models in the `IASDT` framework. [predict_maps] generates current and
#' future habitat suitability maps (mean, sd, cov) from a full Hmsc model fit.
#' [predict_maps_cv] predicts and evaluates cross-validated Hmsc models for
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
#' @param future_max_size	Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details.
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
#' @param lf_only Logical. Whether to predict only the latent factor. This is
#'   useful for distributing processing load between GPU and CPU. When `lf_only
#'   = TRUE`, latent factor prediction needs to be computed separately on GPU.
#'   When computations are finished on GPU, the function can later be rerun with
#'   `lf_only = FALSE` (default) to predict habitat suitability using the
#'   already-computed latent factor predictions.
#' @param climate_models Character vector. Climate models for future
#'   predictions. Available options are `c("GFDL-ESM4", "IPSL-CM6A-LR",
#'   "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")` (default).
#' @param climate_scenario Character vector. Climate scenarios for future
#'   predictions. Available options are: `c("ssp126", "ssp370", "ssp585")`
#'   (default).
#' @param tar_predictions Logical. Whether to compress tiff files for predicted
#'   habitat suitability into a single `*.tar` file (without compression).
#'   Default: `TRUE`.
#' @param model_dir Character. Path to the directory containing cross-validated
#'   models.
#' @param cv_name Character. Cross-validation strategy. Valid values are
#'   `cv_dist`, `cv_large`, or `cv_sac`.
#' @param n_cores_pred Integer. Number of cores to use for predicting species'
#'   habitat suitability.
#' @param cv_fold Integer. The cross-validation fold number.
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
#'   - **`predict_maps_cv`**: Computes predictions for cross-validated `Hmsc`
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
    strategy = "multisession", future_max_size = 1000L, clamp_pred = TRUE,
    fix_efforts = "q90", fix_rivers = "q90", pred_new_sites = TRUE,
    use_tf = TRUE, tf_environ = NULL, tf_use_single = FALSE,
    n_cores_lf = n_cores, lf_check = FALSE, lf_temp_cleanup = TRUE,
    lf_only = FALSE, lf_commands_only = FALSE, temp_dir = "temp_pred",
    temp_cleanup = TRUE, tar_predictions = TRUE, spatial_model = TRUE,
    n_cores_pred = n_cores,
    climate_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    climate_scenario = c("ssp126", "ssp370", "ssp585")) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----

  ecokit::cat_time("Checking input arguments")

  ecokit::check_args(
    args_to_check = c("path_model", "temp_dir"), args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "clamp_pred", "pred_new_sites", "use_tf", "tf_use_single", "lf_check",
      "lf_temp_cleanup", "lf_only", "lf_commands_only", "temp_cleanup",
      "tar_predictions", "spatial_model"),
    args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential")  n_cores <- n_cores_lf <- 1L
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  hab_name <- c(
    "0_all", "1_forests", "2_open_forests", "3_scrub",
    "4a_natural_grasslands", "4b_human_maintained_grasslands",
    "10_wetland", "12a_ruderal_habitats", "12b_agricultural_habitats") %>%
    stringr::str_subset(paste0("^", hab_abb, "_")) %>%
    stringr::str_remove(paste0("^", hab_abb, "_")) %>%
    stringr::str_replace_all("_", " ")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid <- path_roads <- path_railway <- path_efforts <- tif_path <-
    time_period <- climate_scenario <- path_chelsa <- row_id <-
    path_clc <- ias_id <- taxon_name <- climate_model <- climate_scenario <-
    name <- file_pred_sf <- class <- order <- family <- species_name <-
    tif_path_mean <- tif_path_cov <- tif_path_sd <- clamp <- path_wetness <-
    train <- ensemble_file <- ensemble_maps <- tifs <- layer_name <-
    time_period <- file_pred_summary <- ensemble_data <- dir_ensemble <-
    file_pred_r <- tif_path_anomaly <- path_river <- path_soil <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  path_eval <- fs::path(dirname(dirname(path_model)), "model_evaluation")

  path_prediction_1 <- fs::path(
    dirname(dirname(path_model)), "model_prediction")
  path_prediction_clamp <- fs::path(path_prediction_1, "clamp")
  path_prediction_no_clamp <- fs::path(path_prediction_1, "no_clamp")
  fs::dir_create(c(path_eval, path_prediction_no_clamp))

  # Path for overall summary - paths for summaries of identical scenarios
  path_summary_rdata <- fs::path(
    dplyr::if_else(
      clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
    "prediction_summary.RData")
  path_summary_text <- fs::path(
    dplyr::if_else(
      clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
    "prediction_summary.txt")

  # Path for overall summary - for ShinyApp
  path_summary_rdata_shiny <- fs::path(
    dplyr::if_else(
      clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
    "prediction_summary_shiny.RData")
  path_summary_txt_shiny <- fs::path(
    dplyr::if_else(
      clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
    "prediction_summary_shiny.txt")

  # Check if the prediction summary is already available on disk
  if (all(
    fs::file_exists(
      c(
        path_summary_rdata, path_summary_text,
        path_summary_rdata_shiny, path_summary_txt_shiny)))) {
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
        "`fix_efforts` can not be `NULL` when clamping is implemented",
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
        "`fix_rivers` can not be `NULL` when clamping is implemented",
        fix_rivers = fix_rivers, include_backtrace = TRUE)
    }
    if (length(fix_rivers) != 1) {
      # Check if fix_rivers is a vector or length 1
      ecokit::stop_ctx(
        "`fix_rivers` must be of length 1.", fix_rivers = fix_rivers,
        include_backtrace = TRUE)
    }

    # Create folder for clamp results only if clamp_pred == TRUE
    fs::dir_create(path_prediction_clamp)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ------

  ecokit::cat_time("Environment variables")

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_railway", "DP_R_railway_processed", TRUE, FALSE,
    "path_roads", "DP_R_roads_processed", TRUE, FALSE,
    "path_clc", "DP_R_clc_processed", TRUE, FALSE,
    "path_efforts", "DP_R_efforts_processed", TRUE, FALSE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_river", "DP_R_rivers_processed", TRUE, FALSE,
    "path_wetness", "DP_R_wetness_processed", FALSE, TRUE,
    "path_soil", "DP_R_soil_density", FALSE, TRUE,
    "path_chelsa", "DP_R_chelsa_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading input data -----
  ecokit::cat_time("Loading input data")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Species information -----

  species_info <- IASDT.R::get_species_name(env_file = env_file) %>%
    dplyr::select(ias_id, taxon_name, species_name, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Reference grid -----

  ecokit::cat_time("Reference grid", level = 1L)
  path_grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!fs::file_exists(path_grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist",
      path_grid_r = path_grid_r, include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Model object -----

  ecokit::cat_time("Model object", level = 1L)

  if (is.null(path_model) || !fs::file_exists(path_model)) {
    ecokit::stop_ctx(
      "Model path is NULL or does not exist ", path_model = path_model,
      include_backtrace = TRUE)
  }

  model_object <- ecokit::load_as(path_model)

  # `clamp_pred` can not be TRUE when `efforts_log` is not used as predictor
  if (clamp_pred && isFALSE("efforts_log" %in% names(model_object$XData))) {
    ecokit::stop_ctx(
      paste0(
        "`clamp_pred` can not be used when `efforts_log` is",
        "not used as predictor"),
      clamp_pred = clamp_pred, names_data = names(model_object$XData),
      include_backtrace = TRUE)
  }

  other_variables <- paste0(
    "^", c(IASDT.R::chelsa_variables$variable, "cv"), collapse = "|") %>%
    stringr::str_subset(names(model_object$XData), ., negate = TRUE)
  bio_variables <- paste0(
    "^", IASDT.R::chelsa_variables$variable, collapse = "|") %>%
    stringr::str_subset(names(model_object$XData), ., negate = FALSE)

  # Coordinates of the model sampling units
  if (is.null(model_object$rL) || is.null(model_object$rL[[1]]$s)) {
    suppressWarnings({
      model_coords <- tibble::tibble(x = numeric(0L), y = numeric(0L)) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035L) %>%
        dplyr::mutate(train = logical(0L))
    })
  } else {
    model_coords <- as.data.frame(model_object$rL$sample$s) %>%
      tibble::tibble() %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      dplyr::mutate(train = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## CHELSA data -----

  ecokit::cat_time("CHELSA data", level = 1L)

  path_chelsa <- fs::path(path_chelsa, "chelsa_processed_data.RData")
  if (!fs::file_exists(path_chelsa)) {
    ecokit::stop_ctx(
      "Processed CHELSA data can not be found", path_chelsa = path_chelsa,
      include_backtrace = TRUE)
  }

  prediction_options <- ecokit::load_as(path_chelsa) %>%
    dplyr::select(-"file_list") %>%
    dplyr::filter(
      # Filter only selected future climate models
      climate_model %in% c("current", climate_models),
      # Filter only selected future climate scenarios
      climate_scenario %in% c("current", climate_scenario)) %>%
    dplyr::mutate(
      name = paste0(time_period, "_", climate_scenario, "_", climate_model),
      name = stringr::str_replace(name, "1981-2010_current_current", "current"),
      name = stringr::str_replace_all(name, "-", "_"))

  if (clamp_pred) {
    # Add a prediction options for no clamping predictions (model evaluation)
    clamp_option <- prediction_options %>%
      dplyr::filter(climate_model == "current") %>%
      dplyr::mutate(clamp = FALSE)
    prediction_options <- dplyr::mutate(prediction_options, clamp = TRUE) %>%
      dplyr::bind_rows(clamp_option, .)
  } else {
    prediction_options <- dplyr::mutate(prediction_options, clamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Road and railway intensity -----

  static_predictors <- list()

  if ("road_rail_log" %in% other_variables) {

    ecokit::cat_time("Road and railway intensity", level = 1L)

    r_railway <- fs::path(path_railway, "railway_length.RData")
    if (!fs::file_exists(r_railway)) {
      ecokit::stop_ctx(
        "Railways data does not exist", r_railway = r_railway,
        include_backtrace = TRUE)
    }
    r_railway <- ecokit::load_as(r_railway, unwrap_r = TRUE) %>%
      magrittr::extract2("rail")

    r_roads <- fs::path(path_roads, "road_length.RData")
    if (!fs::file_exists(r_roads)) {
      ecokit::stop_ctx(
        "Roads data does not exist", r_roads = r_roads,
        include_backtrace = TRUE)
    }
    r_roads <- ecokit::load_as(r_roads, unwrap_r = TRUE) %>%
      magrittr::extract2("all")

    # Calculating the sum of road and railway intensity
    # add 1 (older versions 0.1) to get log for 0 values [only for
    # rivers/roads/efforts, not hab/rivers]
    r_road_railway <- (r_roads + r_railway) %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("road_rail_log")

    static_predictors <- c(static_predictors, r_road_railway)
    rm(r_road_railway, r_roads, r_railway, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Habitat information ----

  # Check if habitat information is used as predictor
  hab_predictor <- "habitat_log" %in% other_variables

  if (hab_predictor) {

    ecokit::cat_time("Habitat information", level = 1L)
    r_hab <- fs::path(
      path_clc, "summary_rdata", "perc_cover_synhab_crop.RData")
    if (!fs::file_exists(r_hab)) {
      ecokit::stop_ctx(
        "Habitat data does not exist", r_hab = r_hab, include_backtrace = TRUE)
    }

    r_hab <- ecokit::load_as(r_hab, unwrap_r = TRUE) %>%
      magrittr::extract2(paste0("synhab_", hab_abb))

    # Models are trained and predictions are made only at grid cells with > 0 %
    # coverage. Mask layer to exclude grid cells with zero % coverage from
    # predictions.
    # r_hab_mask <- terra::classify(r_hab, cbind(0, NA), others = 1)
    # Update 2025_09: make predictions at the full study area, and if needed
    # predictions at 0% coverage could be masked

    r_hab <- stats::setNames(log10(r_hab + 0.1), "habitat_log")

    static_predictors <- c(static_predictors, r_hab)
    rm(r_hab, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Sampling efforts -----

  if ("efforts_log" %in% other_variables) {

    ecokit::cat_time("Sampling efforts", level = 1L)

    r_efforts <- fs::path(path_efforts, "efforts_summary_r.RData")
    if (!fs::file_exists(r_efforts)) {
      ecokit::stop_ctx(
        "Sampling efforts data does not exist", r_efforts = r_efforts,
        include_backtrace = TRUE)
    }

    # add 1 (older versions 0.1) to get log for 0 values
    # [only for rivers/roads/efforts, not hab/rivers]
    r_efforts <- ecokit::load_as(r_efforts, unwrap_r = TRUE) %>%
      magrittr::extract2("n_obs") %>%
      magrittr::add(1) %>%
      log10() %>%
      stats::setNames("efforts_log")

    if (clamp_pred) {

      ecokit::cat_time("Fixing sampling efforts values", level = 2L)

      # Check fix_efforts value
      if (is.numeric(fix_efforts)) {

        # If `fix_efforts` is numeric value, check if it is within the range of
        # the observed efforts
        efforts_range <- terra::global(r_efforts, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_efforts, efforts_range[1], efforts_range[2]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_efforts` value is out of the range of observed efforts",
            fix_efforts = fix_efforts, efforts_range = round(efforts_range, 2),
            include_backtrace = TRUE)
        }

        # Fix value
        efforts_values <- fix_efforts

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
      efforts_values <- dplyr::case_when(

        # Do not fix if `fix_efforts` is "identity"
        fix_efforts == "identity" ~ NA_real_,

        # Fix at 90% quantile
        fix_efforts == "q90" ~ {
          terra::global(
            r_efforts,
            fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at median value
        fix_efforts == "median" ~ {
          terra::global(
            r_efforts, fun = function(x) median(x, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at mean value
        fix_efforts == "mean" ~ {
          terra::global(r_efforts, fun = mean, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at max value
        fix_efforts == "max" ~ {
          terra::global(r_efforts, fun = max, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        .default = NA_real_)


      # Print fixed value
      if (is.na(efforts_values)) {
        if (fix_efforts == "identity") {
          ecokit::cat_time(
            "Sampling efforts predictor is not fixed",
            level = 2L, cat_timestamp = FALSE)
        }
        r_efforts_clamp <- stats::setNames(r_efforts, "efforts_log_clamp")

      } else {

        ecokit::cat_time(
          paste0("Fixed value is ", round(efforts_values, 2), " [log10 scale]"),
          level = 2L, cat_timestamp = FALSE)

        # Set a minimum value for efforts variable to `efforts_values`. Using
        # upper = Inf keeps efforts values > efforts_values as they are.
        r_efforts_clamp <- terra::clamp(
          x = r_efforts, lower = efforts_values, upper = Inf) %>%
          stats::setNames("efforts_log_clamp")

      }

      static_predictors <- c(static_predictors, r_efforts, r_efforts_clamp)
      rm(r_efforts, r_efforts_clamp, envir = environment())

    } else {

      # Do not fix at single value
      ecokit::cat_time(
        "Sampling efforts is not fixed", level = 2L, cat_timestamp = FALSE)
      static_predictors <- c(static_predictors, r_efforts)
      rm(r_efforts, envir = environment())

    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## River length ----

  if ("rivers_log" %in% other_variables) {

    ecokit::cat_time("River length", level = 1L)

    r_rivers <- fs::path(path_river, "river_lengths.RData")
    if (!fs::file_exists(r_rivers)) {
      ecokit::stop_ctx(
        "River length data does not exist", r_rivers = r_rivers,
        include_backtrace = TRUE)
    }

    r_rivers <- ecokit::load_as(r_rivers, unwrap_r = TRUE) %>%
      magrittr::extract2("STRAHLER_5") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("rivers_log")

    if (clamp_pred) {
      ecokit::cat_time("Fixing river length values", level = 2L)

      # Check fix_rivers value
      if (is.numeric(fix_rivers)) {

        # If `fix_rivers` is numeric value, check if it is within the range of
        # the observed river lengths
        rivers_range <- terra::global(r_rivers, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        invalid_value <- isFALSE(
          dplyr::between(fix_rivers, rivers_range[1], rivers_range[2]))

        if (invalid_value) {
          ecokit::stop_ctx(
            "`fix_rivers` value is out of the range of observed river length",
            fix_rivers = fix_rivers, rivers_range = round(rivers_range, 2),
            include_backtrace = TRUE)
        }

        # Fix value
        rivers_values <- fix_rivers

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
        rivers_values <- dplyr::case_when(

          # Do not fix if `fix_rivers` is "identity"
          fix_rivers == "identity" ~ NA_real_,

          # Fix at 90% quantile
          fix_rivers == "q90" ~ {
            terra::global(
              r_rivers,
              fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at median value
          fix_rivers == "median" ~ {
            terra::global(
              r_rivers, fun = function(x) median(x, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at mean value
          fix_rivers == "mean" ~ {
            terra::global(r_rivers, fun = mean, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at max value
          fix_rivers == "max" ~ {
            terra::global(r_rivers, fun = max, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          .default = NA_real_)
      }

      # Print fixed value

      if (is.na(rivers_values)) {

        if (fix_rivers == "identity") {
          ecokit::cat_time(
            "River length predictor is not fixed at a single value",
            level = 2L, cat_timestamp = FALSE)
        }
        r_rivers_clamp <- stats::setNames(r_rivers, "rivers_log_clamp")

      } else {

        ecokit::cat_time(
          paste0("Fixed value is ", round(rivers_values, 2), " [log10 scale]"),
          level = 2L, cat_timestamp = FALSE)

        # Set a minimum value for river length variable to `rivers_values`.
        # Using upper = Inf keeps  river length values > rivers_values as they
        # are.
        r_rivers_clamp <- terra::clamp(
          x = r_rivers, lower = rivers_values, upper = Inf) %>%
          stats::setNames("rivers_log_clamp")

      }

      static_predictors <- c(static_predictors, r_rivers, r_rivers_clamp)
      rm(r_rivers, r_rivers_clamp, rivers_values, envir = environment())

    } else {

      # Do not fix at single value
      ecokit::cat_time(
        "River length predictor is not fixed at a single value",
        level = 2L, cat_timestamp = FALSE)

      static_predictors <- c(static_predictors, r_rivers)
      rm(r_rivers, envir = environment())

    }

  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Soil bulk density -----

  if ("soil" %in% other_variables) {

    ecokit::cat_time("Soil bulk density", level = 1L)

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

    ecokit::cat_time("Topographic wetness index", level = 1L)

    r_wetness <- fs::path(path_wetness, "wetness_index.RData")
    if (!fs::file_exists(r_wetness)) {
      ecokit::stop_ctx(
        "Wetness data does not exist", r_wetness = r_wetness,
        include_backtrace = TRUE)
    }
    r_wetness <- ecokit::load_as(r_wetness, unwrap_r = TRUE) %>%
      stats::setNames("wetness")

    static_predictors <- c(static_predictors, r_wetness)
    rm(r_wetness, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  ecokit::cat_time("Merge static predictors", level = 1L)

  static_predictors <- terra::rast(static_predictors) %>%
    terra::mask(ecokit::load_as(path_grid_r, unwrap_r = TRUE))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict latent factor at new locations ------

  path_test_lf <- fs::path(path_prediction_1, "test_lf.qs2")

  if (!fs::file_exists(path_test_lf) && pred_new_sites && spatial_model) {

    ecokit::cat_time("Predict latent factor at new locations")

    ecokit::cat_time(
      "Preparing input data for predicting latent factor", level = 1L)

    predict_data_test <- prediction_options %>%
      dplyr::filter(climate_model == "current") %>%
      dplyr::pull("file_path") %>%
      # If clamp_pred`=`TRUE`, there are two options for current climate data
      # (with and without clamping). Two sets of predictions under current
      # climates will be produced. Predictions without clamping is used for
      # model evaluation.
      utils::head(1) %>%
      ecokit::load_as(unwrap_r = TRUE) %>%
      # Only extract predictors used in the model
      terra::subset(bio_variables) %>%
      # Combine with other static predictors
      c(static_predictors) %>% # nolint: consecutive_concatenation_linter
      # If Habitat predictor is used, grid cells with zero % coverage are
      # excluded from predictions [na.rm = TRUE]
      terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
      tibble::tibble() %>%
      sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
      sf::st_join(model_coords) %>%
      tidyr::replace_na(list(train = FALSE)) %>%
      dplyr::filter(!train)

    test_xy <- sf::st_drop_geometry(predict_data_test[, c("x", "y")])
    test_x <- predict_data_test %>%
      dplyr::select(tidyselect::all_of(names(model_object$XData))) %>%
      sf::st_drop_geometry()
    predict_gradient <- Hmsc::prepareGradient(
      hM = model_object, XDataNew = as.data.frame(test_x),
      sDataNew = list(sample = as.data.frame(test_xy)))

    rm(predict_data_test, test_x, test_xy, model_object, envir = environment())
    invisible(gc())

    ecokit::cat_time("Predicting latent factor", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2,
      n_separators = 1L, line_char = "*")

    # Predicting latent factor only
    predictions_lf <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = predict_gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy, future_max_size = future_max_size,
      model_name = paste0("lf_", hab_abb, "_test"), temp_dir = temp_dir,
      temp_cleanup = temp_cleanup, use_tf = use_tf, tf_environ = tf_environ,
      lf_out_file = path_test_lf, tf_use_single = tf_use_single,
      lf_only = TRUE, n_cores_lf = n_cores_lf, lf_check = lf_check,
      lf_temp_cleanup = lf_temp_cleanup, lf_commands_only = lf_commands_only,
      evaluate = FALSE, verbose = TRUE, spatial_model = spatial_model)

    rm(predict_gradient, predictions_lf, envir = environment())

    ecokit::cat_time("Predicting latent factor is finished!", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2,
      n_separators = 1L, line_char = "*")

    if (lf_commands_only) {
      return(invisible(NULL))
    }

  } else {

    if (pred_new_sites && spatial_model) {
      ecokit::cat_time(
        "LF prediction is already available on disk", level = 1L)
    } else {
      ecokit::cat_time("LF prediction will NOT be made", level = 1L)
    }

    rm(model_object, envir = environment())
    invisible(gc())
  }

  if (lf_only) {
    return(invisible(NULL))
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # predict_internal ------

  predict_internal <- function(id) {

    # ID: Index of the current prediction option

    # Whether to clamp the sampling efforts
    do_clamp <- prediction_options$clamp[[id]]

    # Name of the current option
    option_name <- prediction_options$name[[id]]

    # Name of the current model
    model_name <- paste0(
      option_name, " - ", dplyr::if_else(do_clamp, "clamping", "no clamping"))

    print_message <- paste0(
      model_name, " (", id, "/", nrow(prediction_options), ")")
    ecokit::info_chunk(
      print_message, n_separators = 1L, line_char = "-", line_char_rep = 70L,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1L,
      info_lines_before = 1L)

    static_preds <- static_predictors

    if (do_clamp) {

      # Do not evaluate for options with clamping
      evaluate <- FALSE

      # Make prediction files at `path_prediction_clamp`
      path_prediction <- path_prediction_clamp

      if ("efforts_log" %in% names(static_predictors)) {
        # use clamped Effort values
        static_preds$efforts_log <- static_preds$efforts_log_clamp
        static_preds$efforts_log_clamp <- NULL
      }

      if ("rivers_log" %in% names(static_predictors)) {
        # use clamped rivers values
        static_preds$rivers_log <- static_preds$rivers_log_clamp
        static_preds$rivers_log_clamp <- NULL
      }

    } else {

      # evaluate for "current" climates, without clamping
      evaluate <- (option_name == "current")

      # Make prediction files at `path_prediction_no_clamp`
      path_prediction <- path_prediction_no_clamp

      # Remove clamped layers
      if ("efforts_log_clamp" %in% names(static_preds)) {
        static_preds <- terra::subset(
          x = static_preds, subset = "efforts_log_clamp", negate = TRUE)
      }
      if ("rivers_log_clamp" %in% names(static_preds)) {
        static_preds <- terra::subset(
          x = static_preds, subset = "rivers_log_clamp", negate = TRUE)
      }

    }

    path_prediction_sf <- fs::path(
      path_prediction, paste0("prediction_", option_name, "_sf.qs2"))
    path_prediction_r <- fs::path(
      path_prediction, paste0("prediction_", option_name, "_r.qs2"))
    path_prediction_summary <- fs::path(
      path_prediction, paste0("prediction_", option_name, "_summary.RData"))

    # Path for saving tif files of the current option
    path_prediction_tif <- fs::path(path_prediction, option_name)
    fs::dir_create(path_prediction_tif)
    invisible(gc())

    # ______________________________________________
    # ______________________________________________

    # Making predictions if not already processed ----
    out_missing <- fs::file_exists(
      c(path_prediction_r, path_prediction_sf, path_prediction_summary)) %>%
      all() %>%
      isFALSE()

    if (out_missing) {

      # Load model
      model_object <- ecokit::load_as(path_model)

      # Skip predictions if the predictions as sf object is already on disk
      if (fs::file_exists(path_prediction_sf)) {

        ecokit::cat_time("Loading predictions `sf` from disk", level = 1L)
        prediction_sf <- ecokit::load_as(path_prediction_sf)

      } else {

        .option_start_time <- lubridate::now(tzone = "CET")

        # Extracting data at training and new sites ------
        ecokit::cat_time("Extracting data at training and new sites")
        predict_data <- prediction_options$file_path[[id]] %>%
          ecokit::load_as(unwrap_r = TRUE) %>%
          terra::subset(bio_variables) %>%
          c(static_preds) %>%
          terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
          tibble::tibble() %>%
          sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
          sf::st_join(model_coords) %>%
          tidyr::replace_na(list(train = FALSE))

        # Training locations
        model_name_train <- paste0(
          option_name, "_",
          dplyr::if_else(do_clamp, "clamping", "no_clamping"), "_train")
        predict_data_train <- dplyr::filter(predict_data, train)
        if (nrow(predict_data_train) > 0L) {
          train_xy <- sf::st_drop_geometry(predict_data_train[, c("x", "y")])
          train_pa <- as.data.frame(model_object$Y)
          train_x <- predict_data_train %>%
            dplyr::select(tidyselect::all_of(names(model_object$XData))) %>%
            sf::st_drop_geometry() %>%
            stats::model.matrix(model_object$XFormula, ., xlev = NULL)
        } else {
          train_xy <- train_pa <- train_x <- NULL
        }

        # Testing locations
        model_name_test <- paste0(
          option_name, "_",
          dplyr::if_else(do_clamp, "clamping", "no_clamping"), "_test")
        predict_data_test <- dplyr::filter(predict_data, !train)
        if (nrow(predict_data_test) > 0L) {
          test_xy <- predict_data_test[, c("x", "y")]
          test_x <- predict_data_test %>%
            dplyr::select(tidyselect::all_of(names(model_object$XData))) %>%
            sf::st_drop_geometry()
          predict_gradient <- Hmsc::prepareGradient(
            hM = model_object, XDataNew = as.data.frame(test_x),
            sDataNew = list(
              sample = as.data.frame(sf::st_drop_geometry(test_xy))))
        } else {
          test_xy <- test_x <- predict_gradient <- NULL
        }

        rm(model_object, predict_data, envir = environment())
        invisible(gc())

        # ______________________________________________

        if (nrow(predict_data_train) > 0L) {

          # Predicting at training sites ----
          ecokit::cat_time("Predicting at training sites")

          path_current_train <- fs::path(
            path_prediction, paste0("prediction_", option_name, "_train.qs2"))

          if (fs::file_exists(path_current_train)) {

            ecokit::cat_time("Loading predictions from disk", level = 2L)
            preds_fit_sites <- tibble::tibble(pred_path = path_current_train)

          } else {

            preds_fit_sites <- IASDT.R::predict_hmsc(
              path_model = path_model, X = train_x, gradient = NULL,
              expected = TRUE, n_cores = n_cores_pred, strategy = strategy,
              future_max_size = future_max_size, model_name = model_name_train,
              temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_tf = use_tf,
              tf_environ = tf_environ, tf_use_single = tf_use_single,
              lf_return = TRUE, n_cores_lf = n_cores_lf, lf_check = lf_check,
              lf_temp_cleanup = lf_temp_cleanup, lf_commands_only = FALSE,
              pred_directory = path_prediction, pred_pa = train_pa,
              pred_xy = train_xy, evaluate = evaluate, evaluation_name = NULL,
              evaluation_directory = path_eval, verbose = FALSE,
              spatial_model = spatial_model)

          }
        } else {
          ecokit::cat_time(
            paste0(
              "All sites are new sites; no predictions will be made at ",
              "training sites"),
            level = 1L)
          preds_fit_sites <- tibble::tibble(pred_path = NULL)
        }

        # ______________________________________________

        # Predicting at new sites ----

        if (pred_new_sites && nrow(predict_data_test) > 0L) {

          ecokit::cat_time("Predicting at new sites")

          path_current_test <- fs::path(
            path_prediction, paste0("prediction_", option_name, "_test.qs2"))

          if (fs::file_exists(path_current_test)) {

            ecokit::cat_time("Loading predictions from disk", level = 2L)
            preds_new_sites <- tibble::tibble(pred_path = path_current_test)

          } else {

            preds_new_sites <- IASDT.R::predict_hmsc(
              path_model = path_model, gradient = predict_gradient,
              expected = TRUE, n_cores = n_cores_pred, strategy = strategy,
              future_max_size = future_max_size, model_name = model_name_test,
              temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_tf = use_tf,
              tf_environ = tf_environ, tf_use_single = tf_use_single,
              lf_return = TRUE, lf_input_file = path_test_lf,
              n_cores_lf = n_cores_lf, lf_check = lf_check,
              lf_temp_cleanup = lf_temp_cleanup, lf_commands_only = FALSE,
              verbose = FALSE, pred_directory = path_prediction,
              evaluate = FALSE, pred_xy = sf::st_drop_geometry(test_xy),
              spatial_model = spatial_model)

          }
        } else {
          ecokit::cat_time(
            "Predictions at new sites will NOT be made", level = 1L)
          preds_new_sites <- tibble::tibble(pred_path = NULL)
        }

        # ______________________________________________

        # Merge & save predictions - sf ------
        ecokit::cat_time("Merge & save predictions at training and new sites")
        prediction_sf <- purrr::map_dfr(
          .x = c(preds_fit_sites$pred_path, preds_new_sites$pred_path),
          .f = ~ {
            if (fs::file_exists(.x)) {
              ecokit::load_as(.x)
            } else {
              tibble::tibble()
            }
          })

        # ______________________________________________

        # Save predictions as sf object
        ecokit::save_as(object = prediction_sf, out_path = path_prediction_sf)

        ecokit::cat_diff(
          init_time = .option_start_time,
          prefix = "Prediction took ", level = 1L)

      }

      # ______________________________________________
      # ______________________________________________

      ### Predictions as spatRaster / tif -----

      ecokit::cat_time("Rasterization & prepare summary data", level = 1L)

      fields_to_raster <- names(prediction_sf) %>%
        stringr::str_subset("^sp_|^sr_") %>%
        gtools::mixedsort()

      grid_10 <- ecokit::load_as(path_grid_r, unwrap_r = TRUE)
      prediction_r <- terra::rasterize(
        prediction_sf, grid_10, field = fields_to_raster)

      # Calculate prediction anomaly for future projections
      if (option_name != "current") {

        # Names of current taxa (and SR) for the current model
        mean_names <- stringr::str_subset(fields_to_raster, "_mean")

        # Loading mean predictions at current climates
        current_mean <- list.files(
          # use relevant folder containing the current predictions. This is
          # determined by `path_prediction`, which is not the same whether
          # clamping is used or not
          path = path_prediction, pattern = "prediction_current.*_r.qs2",
          full.names = TRUE) %>%
          ecokit::load_as(unwrap_r = TRUE) %>%
          terra::subset(mean_names)

        # Calculate anomaly as difference in predicted value between future and
        # current climate (future - current)
        preds_anomaly <- terra::subset(prediction_r, mean_names) - current_mean
        # Assign names to anomaly maps
        anomaly_names <- names(preds_anomaly) %>%
          stringr::str_replace_all("_mean", "_anomaly")
        names(preds_anomaly) <- anomaly_names

        # Add anomaly maps to the list of predictions
        fields_to_raster <- c(fields_to_raster, anomaly_names)
        prediction_r <- c(prediction_r, preds_anomaly)

        # clean up
        rm(preds_anomaly, current_mean, envir = environment())

      }

      out_summary <- tibble::tibble(
        hab_abb = hab_abb, hab_name = hab_name,
        layer_name = fields_to_raster,
        time_period = prediction_options$time_period[[id]],
        climate_model = prediction_options$climate_model[[id]],
        climate_scenario = prediction_options$climate_scenario[[id]],
        clamp = prediction_options$clamp[[id]],
        path_prediction = path_prediction) %>%
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
            path_prediction_tif, paste0(layer_name, ".tif")))

      # Save as tif
      ecokit::cat_time("Save as tif", level = 1L)
      out_summary_0 <- out_summary %>%
        dplyr::mutate(
          map = purrr::map2(
            .x = layer_name, .y = tif_path,
            .f = ~ {
              terra::writeRaster(
                x = prediction_r[[.x]], filename = .y, overwrite = TRUE,
                gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
            }))

      rm(out_summary_0, envir = environment())

      out_summary <- out_summary %>%
        tidyr::pivot_wider(
          id_cols = c(
            "hab_abb", "hab_name", "time_period", "climate_model",
            "climate_scenario", "ias_id", "clamp", "path_prediction"),
          names_from = "Stats", values_from = "tif_path") %>%
        dplyr::left_join(species_info, by = "ias_id") %>%
        dplyr::select(
          tidyselect::any_of(
            c(
              "hab_abb", "hab_name", "time_period", "climate_model",
              "climate_scenario", "ias_id", "taxon_name", "species_name",
              "class", "order", "family", "tif_path_mean", "tif_path_sd",
              "tif_path_cov", "tif_path_anomaly")))

      # save as spatRaster - qs2
      ecokit::cat_time("Save as spatRaster - qs2", level = 1L)
      prediction_r <- terra::wrap(prediction_r)
      ecokit::save_as(object = prediction_r, out_path = path_prediction_r)

      # Save summary - RData
      ecokit::cat_time("Save summary - RData", level = 1L)
      ecokit::save_as(
        object = out_summary, object_name = paste0(option_name, "_summary"),
        out_path = path_prediction_summary)

      # Save summary - csv
      ecokit::cat_time("Save summary - csv", level = 1L)
      utils::write.table(
        x = out_summary,
        file = stringr::str_replace(path_prediction_summary, ".RData$", ".txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = FALSE, fileEncoding = "UTF-8")
    }

    # output
    tibble::tibble(
      clamp = do_clamp,
      name = option_name,
      file_pred_r = path_prediction_r, file_pred_sf = path_prediction_sf,
      file_pred_summary = path_prediction_summary)
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predicting ------

  ecokit::info_chunk(
    "Making spatial predictions", n_separators = 2L, level = 1L,
    line_char = "*", line_char_rep = 70L, cat_red = TRUE,
    cat_bold = TRUE, cat_timestamp = FALSE)

  grid_10 <- ecokit::load_as(path_grid_r, unwrap_r = TRUE)

  prediction_summary <- purrr::map_dfr(
    .x = seq_len(nrow(prediction_options)), .f = predict_internal) %>%
    dplyr::full_join(prediction_options, ., by = c("name", "clamp")) %>%
    dplyr::select(-"file_path")

  rm(predict_internal, grid_10, model_coords, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  if (length(climate_models) > 1 && !is.null(climate_scenario)) {

    # ensemble model predictions ------

    ecokit::info_chunk(
      "\tensemble model predictions", n_separators = 1L, line_char = "-",
      line_char_rep = 70L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    # Prepare input data to calculate ensemble predictions
    ecokit::cat_time(
      "Prepare input data to calculate ensemble predictions", level = 1L)

    prediction_ensemble <- prediction_summary %>%
      dplyr::filter(climate_model != "current") %>%
      dplyr::select(-file_pred_sf, -file_pred_r, -name, -climate_model) %>%
      dplyr::mutate(
        prediction_2 = purrr::map(
          .x = file_pred_summary,
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
      dplyr::select("prediction_2", "clamp") %>%
      tidyr::unnest("prediction_2") %>%
      dplyr::mutate(
        climate_model = "ensemble",
        dir_ensemble = fs::path(
          dirname(dirname(tif_path_mean)),
          paste0(
            stringr::str_replace(time_period, "-", "_"), "_", climate_scenario,
            "_ensemble"))) %>%
      dplyr::group_by(dplyr::across(-tif_path_mean)) %>%
      dplyr::summarise(tifs = list(tif_path_mean), .groups = "drop") %>%
      dplyr::mutate(
        tif_path_mean = fs::path(dir_ensemble, paste0(ias_id, "_mean.tif")),
        tif_path_anomaly = fs::path(
          dir_ensemble, paste0(ias_id, "_anomaly.tif")),
        tif_path_sd = fs::path(dir_ensemble, paste0(ias_id, "_sd.tif")),
        tif_path_cov = fs::path(dir_ensemble, paste0(ias_id, "_cov.tif")))

    # --------------------------------------------------------- #

    ecokit::cat_time("Create directories for ensemble predictions", level = 1L)
    fs::dir_create(unique(prediction_ensemble$dir_ensemble))

    # --------------------------------------------------------- #

    # Loading mean predictions at current climates
    ecokit::cat_time("Loading mean predictions at current climates", level = 1L)

    current_mean <- list.files(
      path = dplyr::if_else(
        clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
      pattern = "prediction_current.*_r.qs2", full.names = TRUE)

    # --------------------------------------------------------- #

    # Calculate ensemble predictions
    ecokit::cat_time("Calculate ensemble predictions", level = 1L)

    # set up parallel processing
    doParallel::registerDoParallel(cores = n_cores)
    ecokit::load_packages(package_list = "foreach")
    withr::defer(doParallel::stopImplicitCluster())

    prediction_ensemble_0 <- dplyr::select(
      prediction_ensemble,
      tifs, tif_path_mean, tif_path_anomaly, tif_path_sd, tif_path_cov, ias_id)

    calculate_ensemble <- foreach::foreach(
      row_id = seq_len(nrow(prediction_ensemble_0)),
      .export = c("current_mean", "prediction_ensemble_0"),
      .packages = c("terra", "purrr", "ecokit", "qs2", "magrittr", "fs")
    ) %dopar% { # nolint: object_usage_linter

      path_mean <- prediction_ensemble_0$tif_path_mean[[row_id]]
      path_anomaly <- prediction_ensemble_0$tif_path_anomaly[[row_id]]
      path_cov <- prediction_ensemble_0$tif_path_cov[[row_id]]
      path_sd <- prediction_ensemble_0$tif_path_sd[[row_id]]

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

      ias_id <- prediction_ensemble_0$ias_id[[row_id]]

      # load maps for future climate option
      tiffs_r <- terra::rast(prediction_ensemble_0$tifs[[row_id]])

      # Mean
      ensemble_mean <- terra::app(tiffs_r, "mean", na.rm = TRUE)
      # Standard deviation
      ensemble_sd <- terra::app(tiffs_r, "sd", na.rm = TRUE)
      rm(tiffs_r)

      # Anomaly
      current_mean_0 <- ecokit::load_as(current_mean, unwrap_r = TRUE) %>%
        terra::subset(paste0(ias_id, "_mean"))
      ensemble_anomaly <- ensemble_mean - current_mean_0
      rm(current_mean_0)

      # coefficient of variation: Replace very small mean values with
      # reasonably small number to avoid overflow warning
      ensemble_mean_0 <- terra::classify(
        x = ensemble_mean, rcl = cbind(0, 1e-8, 1e-8))
      ensemble_cov <- (ensemble_sd / ensemble_mean_0)

      gdal_o <- c("COMPRESS=DEFLATE", "TILED=YES")
      terra::writeRaster(
        x = ensemble_mean, filename = path_mean,
        overwrite = TRUE, gdal = gdal_o)
      terra::writeRaster(
        x = ensemble_anomaly, filename = path_anomaly,
        overwrite = TRUE, gdal = gdal_o)
      terra::writeRaster(
        x = ensemble_cov, filename = path_cov,
        overwrite = TRUE, gdal = gdal_o)
      terra::writeRaster(
        x = ensemble_sd, filename = path_sd,
        overwrite = TRUE, gdal = gdal_o)

      return(NULL)
    }

    rm(
      current_mean, prediction_ensemble_0, calculate_ensemble,
      envir = environment())
    invisible(gc())

    # --------------------------------------------------------- #

    # Save ensemble maps as SpatRast
    ecokit::cat_time("Save ensemble predictions as SpatRast", level = 1L)

    vars_to_select <- c(
      "ias_id", "ensemble_file", "tif_path_mean", "tif_path_anomaly",
      "tif_path_sd", "tif_path_cov")
    prediction_ensemble_r <- prediction_ensemble %>%
      dplyr::select(-dir_ensemble, -tifs) %>%
      dplyr::mutate(
        ensemble_file = fs::path(
          dplyr::if_else(
            clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
          paste0(
            "prediction_", stringr::str_replace(time_period, "-", "_"), "_",
            climate_scenario, "_ensemble_r.qs2"))) %>%
      dplyr::select(tidyselect::all_of(vars_to_select)) %>%
      tidyr::nest(data = -ensemble_file)

    prediction_ensemble_r <- foreach::foreach(
      row_id = seq_len(nrow(prediction_ensemble_r)),
      .export = "prediction_ensemble_r",
      .packages = c("dplyr", "magrittr", "terra", "purrr", "ecokit", "fs")
    ) %dopar% { # nolint: object_usage_linter

      out_maps <- dplyr::slice(prediction_ensemble_r, row_id) %>%
        dplyr::pull("data") %>%
        magrittr::extract2(1) %>%
        dplyr::mutate(
          maps = purrr::pmap(
            .l = list(ias_id, tif_path_mean, tif_path_anomaly,
                      tif_path_sd, tif_path_cov),
            .f = function(ias_id, tif_path_mean, tif_path_anomaly,
                          tif_path_sd, tif_path_cov) {

              mean <- terra::rast(tif_path_mean) %>%
                stats::setNames(paste0(ias_id, "_mean"))
              sd <- terra::rast(tif_path_sd) %>%
                stats::setNames(paste0(ias_id, "_sd"))
              cov <- terra::rast(tif_path_cov) %>%
                stats::setNames(paste0(ias_id, "_cov"))
              anomaly <- terra::rast(tif_path_anomaly) %>%
                stats::setNames(paste0(ias_id, "_anomaly"))

              c(mean, sd, cov, anomaly) %>%
                terra::wrap()
            }
          )) %>%
        dplyr::pull("maps") %>%
        purrr::map(terra::unwrap) %>%
        terra::rast() %>%
        terra::wrap()

      ecokit::save_as(
        object = out_maps,
        out_path = prediction_ensemble_r$ensemble_file[[row_id]])
      return(NULL)
    }

    doParallel::stopImplicitCluster()

    rm(prediction_ensemble_r, envir = environment())

    # --------------------------------------------------------- #

    # Save summary of ensemble predictions
    ecokit::cat_time("Save summary of ensemble predictions", level = 1L)

    prediction_ensemble_summary <- prediction_ensemble %>%
      dplyr::select(-ensemble_maps) %>%
      dplyr::mutate(
        ensemble_file = fs::path(
          dplyr::if_else(
            clamp_pred, path_prediction_clamp, path_prediction_no_clamp),
          paste0(
            "prediction_", stringr::str_replace(time_period, "-", "_"), "_",
            climate_scenario, "_ensemble_summary.RData")),
        # utils::write.table does not support `fs::path()`
        # convert paths to character and use readr::write_delim instead
        dplyr::across(
          .cols = tidyselect::all_of(
            c("dir_ensemble", "tifs", "tif_path_mean",
              "tif_path_anomaly", "tif_path_sd", "tif_path_cov")),
          .fns = as.character)) %>%
      tidyr::nest(ensemble_data = -ensemble_file) %>%
      dplyr::mutate(
        ensemble_save = purrr::map2(
          .x = ensemble_data, .y = ensemble_file,
          .f = ~ {
            ecokit::save_as(
              object = .x, out_path = .y,
              object_name = stringr::str_remove(basename(.y), ".RData"))
            readr::write_delim(
              x = .x, file = stringr::str_replace(.y, ".RData", ".txt"),
              delim = "\t", col_names = TRUE,
              append = FALSE, progress = FALSE)
          }),
        ensemble_save = NULL) %>%
      tidyr::unnest(ensemble_data) %>%
      dplyr::select(
        tidyselect::all_of(
          c("ensemble_file", "time_period",
            "climate_model", "climate_scenario"))) %>%
      dplyr::rename(file_pred_summary = ensemble_file) %>%
      dplyr::mutate(
        clamp = clamp_pred,
        name = paste0(
          stringr::str_replace(time_period, "-", "_"), "_",
          climate_scenario, "_ensemble"),
        file_pred_sf = NA_character_,
        file_pred_r = stringr::str_replace(
          file_pred_summary, "_summary.RData", "_r.qs2")) %>%
      dplyr::distinct()
  }
  # # ..................................................................... ###
  # # ..................................................................... ###

  # Overall summary -----

  ecokit::info_chunk(
    "\tPrepare overall summary", n_separators = 1L, line_char = "-",
    line_char_rep = 70L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  prediction_summary <- dplyr::rename(
    prediction_summary,
    time_period = time_period, climate_model = climate_model,
    climate_scenario = climate_scenario)

  if (length(climate_models) > 1 && !is.null(climate_scenario)) {
    prediction_summary <- dplyr::bind_rows(
      prediction_summary, prediction_ensemble_summary)
  }

  prediction_summary <- dplyr::mutate(
    prediction_summary, hab_abb = hab_abb, hab_name = hab_name, .before = 1)

  readr::write_delim(
    x = prediction_summary, file = path_summary_text, delim = "\t",
    col_names = TRUE, append = FALSE, progress = FALSE)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Overall summary - to be uploaded to the data server; for the Shiny App -----

  if (clamp_pred) {
    prediction_summary_shiny <- dplyr::filter(prediction_summary, clamp)
  } else {
    prediction_summary_shiny <- dplyr::filter(prediction_summary, !clamp)
  }

  prediction_summary_shiny <- prediction_summary_shiny$file_pred_summary %>%
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
            pattern = paste0(dirname(path_summary_rdata_shiny), "/"))
        }))

  save(prediction_summary_shiny, file = path_summary_rdata_shiny)
  utils::write.table(
    x = prediction_summary_shiny, sep = "\t", row.names = FALSE,
    col.names = TRUE, file = path_summary_txt_shiny, quote = FALSE,
    fileEncoding = "UTF-8")

  # Create tar file for prediction files
  if (tar_predictions) {

    ecokit::cat_time("Create tar file for prediction files", level = 1L)

    # Directory to save the tar file
    tar_dir <- dirname(path_summary_rdata_shiny)
    # Path to the tar file
    tar_file <- fs::path(tar_dir, "predictions.tar")
    # List of directories in the prediction folder. All directories will be
    # included in the tar file
    tar_files <- list.dirs(      # nolint: object_name_linter
      path = tar_dir, full.names = FALSE, recursive = FALSE) %>%
      paste(collapse = " ") %>%
      # Add the summary files to the list
      paste(
        "prediction_summary_shiny.RData", "prediction_summary_shiny.txt",
        collapse = " ")

    # Command to create the tar file
    tar_command <- stringr::str_glue(
      "cd {fs::path_abs(tar_dir)}; tar -cf {basename(tar_file)} -b 2048 \\
      {tar_files}")

    # Create tar file
    system(tar_command)

    # Change the permission of the tar file
    Sys.chmod(tar_file, "755", use_umask = FALSE)
  }

  # # ................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nThe whole prediction function took ")

  return(prediction_summary)
}
