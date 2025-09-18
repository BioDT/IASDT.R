## |------------------------------------------------------------------------| #
# predict_maps_cv ----
## |------------------------------------------------------------------------| #

#' @export
#' @name predict_maps
#' @rdname predict_maps
#' @order 2
#' @author Ahmed El-Gabbas

predict_maps_cv <- function(
    model_dir = NULL, cv_name = NULL, cv_fold = NULL, n_cores = 8L,
    strategy = "multisession", env_file = ".env", use_tf = TRUE,
    tf_environ = NULL, tf_use_single = FALSE, n_cores_lf = n_cores,
    lf_check = FALSE, lf_temp_cleanup = TRUE, lf_only = FALSE,
    lf_commands_only = FALSE, temp_cleanup = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----
  ecokit::cat_time("Checking input arguments")

  ecokit::check_args(
    args_to_check = c("cv_name", "model_dir"), args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "use_tf", "tf_use_single", "lf_check", "lf_temp_cleanup", "lf_only",
      "lf_commands_only", "temp_cleanup"),
    args_type = "logical")

  cv_name <- .validate_cv_name(cv_name)
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- n_cores_lf <- 1L

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  x <- y <- path_grid <- ias_id <- taxon_name <- species_name <- sp_type <-
    type <- n_grids_pres <- n_grids_abs <- sp <- ias_id <- ncells <-
    layer_name <- tif_path <- class <- order <- family <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

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
  if (!file.exists(path_grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist",
      path_grid_r = path_grid_r, include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Reading input data + managing File and directory paths

  ecokit::cat_time("Reading input data + managing File and directory paths")

  # Loading CV data
  cv_data <- fs::path(model_dir, "model_fitting_cv", "cv_data_fitted.RData")
  if (!ecokit::check_data(cv_data, warning = FALSE)) {
    ecokit::stop_ctx(
      "CV data does not exist", cv_data = cv_data, include_backtrace = TRUE)
  }

  cv_data <- ecokit::load_as(cv_data) %>%
    dplyr::filter(.data$cv_name == cv_name, .data$cv == cv_fold)

  # model name
  hab_abb <- stringr::str_extract(model_dir, "(Hab|hab).+") %>%
    stringr::str_remove("/$") %>%
    stringr::str_remove("Hab|hab")
  # full habitat name
  hab_name <- c(
    "0_all", "1_forests", "2_open_forests", "3_scrub",
    "4a_natural_grasslands", "4b_human_maintained_grasslands",
    "10_wetland", "12a_ruderal_habitats", "12b_agricultural_habitats") %>%
    stringr::str_subset(paste0("^", hab_abb, "_")) %>%
    stringr::str_remove(paste0("^", hab_abb, "_")) %>%
    stringr::str_replace_all("_", " ")
  # model name
  model_name <- paste0(hab_abb, "_", cv_data$ModName)

  # path for model evaluation
  path_eval <- fs::path(model_dir, "model_fitting_cv", "evaluation")

  # path for model predictions
  path_prediction <- fs::path(
    model_dir, "model_fitting_cv", "model_prediction", model_name)

  temp_dir <- fs::path(
    model_dir, "model_fitting_cv", "Temp", paste0("Temp_", model_name))
  fs::dir_create(c(path_eval, path_prediction, temp_dir))

  path_preds_sf <- fs::path(
    path_prediction, paste0("prediction_", model_name, ".qs2"))
  path_eval_file <- fs::path(path_eval, paste0("eval_", model_name, ".qs2"))
  path_preds_r <- fs::path(
    path_prediction, paste0("prediction_", model_name, "_R.qs2"))
  path_preds_summary <- fs::path(
    path_prediction, paste0("prediction_", model_name, "_summary.RData"))

  # CV fitted cross-validated model
  path_model <- cv_data$path_mod_fitted
  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model data does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading data used in full model (without cross-validation)
  ecokit::cat_time("Loading data used in full model")

  model_data <- fs::path(model_dir, "model_data_subset.RData")
  if (!ecokit::check_data(model_data)) {
    ecokit::stop_ctx(
      "Model data file not found",
      model_data = model_data, include_backtrace = TRUE)
  }
  model_data <- ecokit::load_as(model_data)

  # ..................................................................... ###
  # ..................................................................... ###

  # Prepare testing data sets
  ecokit::cat_time("Prepare testing datasets")

  test_data <- dplyr::filter(model_data$data_all, .data[[cv_name]] == cv_fold)
  if (nrow(test_data) == 0) {
    ecokit::stop_ctx(
      "No data available for the current CV fold",
      cv_name = cv_name, cv_fold = cv_fold, include_backtrace = TRUE)
  }

  test_x <- test_data %>%
    dplyr::select(tidyselect::all_of(names(model_data$data_x))) %>%
    as.data.frame()
  test_xy <- as.data.frame(dplyr::select(test_data, x, y))

  model_obj <- ecokit::load_as(path_model)
  gradient <- Hmsc::prepareGradient(
    hM = model_obj, XDataNew = test_x, sDataNew = list(sample = test_xy))
  test_pa <- dplyr::select(test_data, tidyselect::starts_with("sp_")) %>%
    as.data.frame()

  n_grid_summary <- test_pa %>%
    dplyr::summarise(
      dplyr::across(
        .cols = tidyselect::everything(),
        .fns = list(n_grids_pres = sum, n_grids_abs = ~sum(.x == 0)),
        .names = "{.col}__{.fn}")) %>%
    tidyr::pivot_longer(
      col = tidyselect::everything(),
      names_to = "sp_type", values_to = "ncells") %>%
    tidyr::separate_wider_delim(
      cols = sp_type, delim = "__", names = c("ias_id", "type")) %>%
    tidyr::pivot_wider(names_from = type, values_from = ncells) %>%
    dplyr::mutate(dplyr::across(c(n_grids_pres, n_grids_abs), as.integer))

  rm(model_obj, test_x, envir = environment())
  invisible(gc())

  # ..................................................................... ###
  # ..................................................................... ###

  # Predict latent factor at new locations ------

  path_test_lf <- fs::path(
    path_prediction, paste0("lf_test_", model_name, ".qs2"))

  if (file.exists(path_test_lf)) {

    ecokit::cat_time("LF prediction is already available on disk")

  } else {

    ecokit::info_chunk("Predict latent factor at new locations")

    # Predicting latent factor only
    preds_lf <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy, model_name = model_name,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_tf = use_tf,
      tf_environ = tf_environ, lf_out_file = path_test_lf,
      tf_use_single = tf_use_single, lf_only = TRUE, n_cores_lf = n_cores_lf,
      lf_check = lf_check, lf_temp_cleanup = lf_temp_cleanup,
      lf_commands_only = lf_commands_only, evaluate = FALSE, verbose = TRUE)

    rm(preds_lf, envir = environment())

    ecokit::cat_time("Predicting latent factor is finished!", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2L, n_separators = 1L,
      line_char = "*")

    if (lf_commands_only) {
      return(invisible(NULL))
    }
  }

  if (lf_only) {
    return(invisible())
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict habitat suitability at testing cross-validation folds ------

  if (all(file.exists(path_preds_sf, path_eval_file))) {

    # Skip predictions if the predictions as sf object and evaluation data
    # already on disk
    ecokit::cat_time("Loading predictions `sf` from disk")
    prediction_sf <- ecokit::load_as(path_preds_sf)

  } else {

    ecokit::info_chunk(
      "Predict habitat suitability at testing cross-validation folds",
      line_char_rep = 75)

    prediction_sf <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy, model_name = model_name,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_tf = use_tf,
      tf_environ = tf_environ, tf_use_single = tf_use_single, lf_return = TRUE,
      lf_input_file = path_test_lf, n_cores_lf = n_cores_lf,
      lf_check = lf_check, lf_temp_cleanup = lf_temp_cleanup,
      lf_commands_only = FALSE, verbose = TRUE,
      pred_directory = path_prediction, evaluate = TRUE,
      evaluation_directory = path_eval, pred_xy = test_xy, pred_pa = test_pa)

    #  --------------------------------------------------------------------- #

    if (!file.exists(path_preds_sf)) {
      ecokit::stop_ctx(
        "Prediction file does not exist", path_preds_sf = path_preds_sf,
        include_backtrace = TRUE)

    }
    if (isFALSE(ecokit::check_data(path_preds_sf, warning = FALSE))) {
      ecokit::stop_ctx(
        "Prediction file is corrupted", path_preds_sf = path_preds_sf,
        include_backtrace = TRUE)
    }

    if (!file.exists(path_eval_file)) {
      ecokit::stop_ctx(
        "Evaluation file does not exist", path_eval_file = path_eval_file,
        include_backtrace = TRUE)
    }

    if (isFALSE(ecokit::check_data(path_eval_file, warning = FALSE))) {
      ecokit::stop_ctx(
        "Evaluation file is corrupted", path_eval_file = path_eval_file,
        include_backtrace = TRUE)
    }

    # loading evaluation data
    eval_data <- ecokit::load_as(prediction_sf$eval_path) %>%
      dplyr::select(-sp) %>%
      dplyr::rename(ias_id = ias_id)

    # loading prediction data
    prediction_sf <- ecokit::load_as(prediction_sf$pred_path)

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Rasterizing predictions and making summary ------

  if (all(file.exists(path_preds_r, path_preds_summary))) {

    out_summary <- ecokit::load_as(path_preds_summary)

  } else {

    ecokit::cat_time("\nRasterization predictions at testing sites")

    species_info <- IASDT.R::get_species_name(env_file = env_file) %>%
      janitor::clean_names() %>%
      dplyr::select(ias_id, taxon_name, species_name, class, order, family)

    fields_to_raster <- names(prediction_sf) %>%
      stringr::str_subset("^sp_|^sr_") %>%
      gtools::mixedsort()

    grid_10 <- ecokit::load_as(path_grid_r, unwrap_r = TRUE)
    prediction_r <- terra::rasterize(
      prediction_sf, grid_10, field = fields_to_raster)

    #  --------------------------------------------------------------------- #

    ecokit::cat_time("Save prediction outputs")

    # Save as tif
    ecokit::cat_time("Save predictions as tif files", level = 1L)

    selected_columns <- c(
      "hab_abb", "hab_name", "cv_name", "cv_fold", "ias_id", "taxon_name",
      "species_name", "class", "order", "family", "tif_path_mean",
      "tif_path_sd", "tif_path_cov", "n_grids_pres", "n_grids_abs")

    out_summary <- tibble::tibble(
      layer_name = fields_to_raster, hab_abb = hab_abb, hab_name = hab_name,
      cv_name = cv_name, cv_fold = cv_fold) %>%
      dplyr::mutate(
        Stats = dplyr::case_when(
          endsWith(layer_name, "_mean") ~ "tif_path_mean",
          endsWith(layer_name, "_sd") ~ "tif_path_sd",
          endsWith(layer_name, "_cov") ~ "tif_path_cov",
          .default = NULL),
        ias_id = stringr::str_remove(layer_name, "_mean$|_sd$|_cov$"),
        tif_path = fs::path(
          path_prediction, paste0(layer_name, "_", model_name, ".tif")),
        out_map = purrr::map2(
          .x = layer_name, .y = tif_path,
          .f = ~ {
            terra::writeRaster(
              x = prediction_r[[.x]], filename = .y, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
          }),
        out_map = NULL) %>%
      tidyr::pivot_wider(
        id_cols = c("hab_abb", "hab_name", "cv_name", "cv_fold", "ias_id"),
        names_from = "Stats", values_from = "tif_path") %>%
      dplyr::left_join(species_info, by = "ias_id") %>%
      dplyr::left_join(n_grid_summary, by = "ias_id") %>%
      dplyr::select(tidyselect::all_of(selected_columns)) %>%
      dplyr::left_join(eval_data, by = "ias_id") %>%
      dplyr::slice(gtools::mixedorder(ias_id))

    ecokit::cat_time("Save summary data", level = 1L)
    # save as spatRaster - qs2
    ecokit::cat_time("Save as spatRaster - qs2", level = 1L)
    prediction_r <- terra::wrap(prediction_r)
    ecokit::save_as(object = prediction_r, out_path = path_preds_r)

    # Save summary - RData
    ecokit::cat_time("Save summary - RData", level = 1L)
    ecokit::save_as(
      object = out_summary,
      object_name = paste0("prediction_", model_name, "_summary"),
      out_path = path_preds_summary)

    # Save summary - csv
    ecokit::cat_time("Save summary - csv", level = 1L)
    utils::write.table(
      x = out_summary,
      file = stringr::str_replace(path_preds_summary, ".RData$", ".txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
      fileEncoding = "UTF-8")
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nThe whole prediction function took ")

  return(path_preds_summary)
}
