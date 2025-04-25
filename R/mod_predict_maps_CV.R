## |------------------------------------------------------------------------| #
# predict_maps_CV ----
## |------------------------------------------------------------------------| #

#' Predict and evaluate cross-validated `Hmsc` models
#'
#' This function computes predicted values of cross-validated `Hmsc` models in
#' the testing cross-validation folders and evaluates the model performance
#' using different metrics. In contrast to [predict_maps], this function
#' predicts only under current climate conditions and does not use clamping.
#'
#' @param model_dir Character. Path to the root directory where the
#'   cross-validated models were fitted.
#' @param CV_name Character. Name of the cross-validation strategy used. Valid
#'   values are `CV_Dist`, `CV_Large`, or `CV_SAC`.
#' @param CV_fold Integer. The cross-validation fold number.
#' @export
#' @name predict_maps_CV
#' @author Ahmed El-Gabbas
#' @inheritParams predict_hmsc
#' @inheritParams predict_maps
#' @return A tibble containing the prediction summary and file paths for output
#'   `*.tif` files.
#' @seealso [predict_maps]
#' @export

predict_maps_CV <- function(
    model_dir = NULL, CV_name = NULL, CV_fold = NULL, n_cores = 8L,
    env_file = ".env", use_TF = TRUE, TF_environ = NULL, TF_use_single = FALSE,
    LF_n_cores = n_cores, LF_check = FALSE, LF_temp_cleanup = TRUE,
    LF_only = FALSE, LF_commands_only = FALSE, temp_cleanup = TRUE) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  x <- y <- Path_Grid <- ias_id <- taxon_name <- species_name <-
    layer_name <- tif_path <-  class <- order <- family <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----
  IASDT.R::cat_time("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("CV_name", "model_dir", "env_file"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "use_TF", "TF_use_single", "LF_check", "LF_temp_cleanup",
      "LF_temp_cleanup", "LF_only", "LF_commands_only", "temp_cleanup"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "CV_fold", "LF_n_cores"))

  rm(AllArgs, envir = environment())

  if (!CV_name %in% c("CV_Dist", "CV_Large", "CV_SAC")) {
    IASDT.R::stop_ctx(
      paste0(
        "Invalid value for CV_name argument.\nValid values are: CV_Dist, ",
        "CV_Large, or CV_SAC"),
      CV_name = CV_name)
  }

  if (!file.exists(env_file)) {
    IASDT.R::stop_ctx(
      "Environment file is invalid or does not exist.", env_file = env_file)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ------

  IASDT.R::cat_time("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Path_GridR <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    IASDT.R::stop_ctx(
      "Path for the Europe boundaries does not exist", Path_GridR = Path_GridR)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Reading input data + managing File and directory paths

  IASDT.R::cat_time("Reading input data + managing File and directory paths")

  # Loading CV data
  CV_DT <- IASDT.R::path(model_dir, "Model_Fitting_CV", "CV_DT_fitted.RData")
  if (!IASDT.R::check_data(CV_DT, warning = FALSE)) {
    IASDT.R::stop_ctx("CV data does not exist", CV_DT = CV_DT)

  }
  CV_DT <- IASDT.R::load_as(CV_DT) %>%
    dplyr::filter(.data$CV_name == CV_name, .data$CV == CV_fold)

  # model name
  habitat_abb <- stringr::str_extract(model_dir, "(Hab|hab).+") %>%
    stringr::str_remove("/$") %>%
    stringr::str_remove("Hab|hab")
  model_name <- paste0(habitat_abb, "_", CV_DT$ModName)

  # path for model evaluation
  Path_Eval <- IASDT.R::path(model_dir, "Model_Fitting_CV", "Evaluation")

  # path for model predictions
  Path_Prediction <- IASDT.R::path(
    model_dir, "Model_Fitting_CV", "Model_Prediction", model_name)

  temp_dir <- IASDT.R::path(
    model_dir, "Model_Fitting_CV", "Temp", paste0("Temp_", model_name))
  fs::dir_create(c(Path_Eval, Path_Prediction, temp_dir))

  Path_Preds_sf <- IASDT.R::path(
    Path_Prediction, paste0("Prediction_", model_name, "_sf.qs2"))
  Path_Eval_File <- IASDT.R::path(
    Path_Eval, paste0("Eval_", model_name, ".qs2"))
  Path_Preds_R <- IASDT.R::path(
    Path_Prediction, paste0("Prediction_", model_name, "_R.qs2"))
  Path_Preds_summary <- IASDT.R::path(
    Path_Prediction, paste0("Prediction_", model_name, "_Summary.RData"))

  # CV fitted cross-validated model
  path_model <- CV_DT$Path_ModFitted
  if (!file.exists(path_model)) {
    IASDT.R::stop_ctx("Model data does not exist", path_model = path_model)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading data used in full model (without cross-validation)
  IASDT.R::cat_time("Loading data used in full model")

  model_data <- list.files(
    path = model_dir, pattern = "ModDT_.+_subset.RData", full.names = TRUE)
  if (length(model_data) != 1) {
    IASDT.R::stop_ctx(
      paste0("There are ", length(model_data), " model data files"),
      model_dir = model_dir)
  }
  model_data <- IASDT.R::load_as(model_data)

  # ..................................................................... ###
  # ..................................................................... ###

  # Prepare testing data sets
  IASDT.R::cat_time("Prepare testing datasets")

  Test_DT <- dplyr::filter(model_data$DT_All, .data[[CV_name]] == CV_fold)
  if (nrow(Test_DT) == 0) {
    IASDT.R::stop_ctx(
      "No data available for the current CV fold",
      CV_name = CV_name, CV_fold = CV_fold)
  }

  Test_X <- Test_DT %>%
    dplyr::select(tidyselect::all_of(names(model_data$DT_x))) %>%
    as.data.frame()
  Test_XY <- as.data.frame(dplyr::select(Test_DT, x, y))

  Model <- IASDT.R::load_as(path_model)
  Gradient <- Hmsc::prepareGradient(
    hM = Model, XDataNew = Test_X, sDataNew = list(sample = Test_XY))
  Test_PA <- dplyr::select(Test_DT, tidyselect::starts_with("Sp_")) %>%
    as.data.frame()

  rm(Model, Test_X, envir = environment())
  invisible(gc())

  # ..................................................................... ###
  # ..................................................................... ###

  # Predict latent factor at new locations ------

  IASDT.R::cat_time("Predict latent factor at new locations")

  Path_Test_LF <- IASDT.R::path(
    Path_Prediction, paste0("LF_Test_", model_name, ".qs2"))

  if (file.exists(Path_Test_LF)) {

    IASDT.R::cat_time("LF prediction is already available on disk", level = 1)

  } else {

    IASDT.R::cat_time("Predicting latent factor", level = 1)
    IASDT.R::cat_sep(
      sep_lines_before = 1, sep_lines_after = 2,
      n_separators = 1, line_char = "*")

    # Predicting latent factor only
    Preds_LF <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, model_name = model_name,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_TF = use_TF,
      TF_environ = TF_environ, LF_out_file = Path_Test_LF,
      TF_use_single = TF_use_single, LF_only = TRUE, LF_n_cores = LF_n_cores,
      LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup,
      LF_commands_only = LF_commands_only, evaluate = FALSE, verbose = TRUE)

    IASDT.R::cat_time("Predicting latent factor is finished!", level = 1)
    IASDT.R::cat_sep(
      sep_lines_before = 1, sep_lines_after = 2,
      n_separators = 1, line_char = "*")

    if (LF_commands_only) {
      return(invisible(NULL))
    }
  }

  if (LF_only) {
    return(invisible())
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict habitat suitability at testing cross-validation folds ------



  if (all(file.exists(Path_Preds_sf, Path_Eval_File))) {

    # Skip predictions if the predictions as sf object and evaluation data
    # already on disk
    IASDT.R::cat_time("Loading predictions `sf` from disk", level = 1)
    Prediction_sf <- IASDT.R::load_as(Path_Preds_sf)

  } else {

    IASDT.R::cat_time("predicting at new sites", level = 1)
    .OptionStartTime <- lubridate::now(tzone = "CET")

    Prediction_sf <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, model_name = model_name, temp_dir = temp_dir,
      temp_cleanup = temp_cleanup, use_TF = use_TF, TF_environ = TF_environ,
      TF_use_single = TF_use_single, LF_return = TRUE,
      LF_inputFile = Path_Test_LF, LF_n_cores = LF_n_cores, LF_check = LF_check,
      LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
      verbose = FALSE, pred_directory = Path_Prediction, evaluate = TRUE,
      evaluation_directory = Path_Eval, pred_XY = Test_XY, pred_PA = Test_PA)


    if (!file.exists(Path_Preds_sf)) {
      IASDT.R::stop_ctx(
        "Prediction file does not exist", Path_Preds_sf = Path_Preds_sf)

    }
    if (isFALSE(IASDT.R::check_data(Path_Preds_sf, warning = FALSE))) {
      IASDT.R::stop_ctx(
        "Prediction file is corrupted", Path_Preds_sf = Path_Preds_sf)
    }

    if (!file.exists(Path_Eval_File)) {
      IASDT.R::stop_ctx(
        "Evaluation file does not exist", Path_Eval_File = Path_Eval_File)
    }

    if (isFALSE(IASDT.R::check_data(Path_Eval_File, warning = FALSE))) {
      IASDT.R::stop_ctx(
        "Evaluation file is corrupted", Path_Eval_File = Path_Eval_File)
    }

    # loading prediction data
    Prediction_sf <- IASDT.R::load_as(Prediction_sf$Pred_Path)

    # elapsed time
    IASDT.R::cat_diff(
      init_time = .OptionStartTime, prefix = "Prediction took ", level = 2)

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Rasterizing predictions and making summary ------

  if (all(file.exists(Path_Preds_R, Path_Preds_summary))) {

    Out_Summary <- IASDT.R::load_as(Path_Preds_summary)

  } else {

    IASDT.R::cat_time("Rasterization predictions at testing sites", level = 1)

    SpeciesInfo <- IASDT.R::get_species_name(env_file = env_file) %>%
      janitor::clean_names() %>%
      dplyr::select(ias_id, taxon_name, species_name, class, order, family)

    Fields2Raster <- names(Prediction_sf) %>%
      stringr::str_subset("^Sp_|^SR_") %>%
      gtools::mixedsort()

    Grid10 <- terra::unwrap(IASDT.R::load_as(Path_GridR))
    Prediction_R <- terra::rasterize(
      Prediction_sf, Grid10, field = Fields2Raster)

    # Save as tif
    IASDT.R::cat_time("Save predictions as tif files", level = 1)

    Out_Summary <- tibble::tibble(
      layer_name = Fields2Raster, hab_abb = habitat_abb,
      CV_name = CV_name, CV_fold = CV_fold) %>%
      dplyr::mutate(
        Stats = dplyr::case_when(
          stringr::str_detect(layer_name, "_mean$") ~ "tif_path_mean",
          stringr::str_detect(layer_name, "_sd$") ~ "tif_path_sd",
          stringr::str_detect(layer_name, "_cov$") ~ "tif_path_cov",
          .default = NULL),
        ias_id = stringr::str_remove(layer_name, "_mean$|_sd$|_cov$"),
        tif_path = IASDT.R::path(
          Path_Prediction, paste0(layer_name, "_", model_name, ".tif")),
        Outt_Map = purrr::map2(
          .x = layer_name, .y = tif_path,
          .f = ~ {
            terra::writeRaster(
              x = Prediction_R[[.x]], filename = .y, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
          }),
        Outt_Map = NULL) %>%
      tidyr::pivot_wider(
        id_cols = c("hab_abb",  "CV_name", "CV_fold", "ias_id"),
        names_from = "Stats", values_from = "tif_path") %>%
      dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
      dplyr::select(
        tidyselect::all_of(
          c("hab_abb",  "CV_name", "CV_fold", "ias_id", "taxon_name",
            "species_name", "class", "order", "family", "tif_path_mean",
            "tif_path_sd", "tif_path_cov")))

    IASDT.R::cat_time("Save summary data", level = 1)

    # save as spatRaster - qs2
    IASDT.R::cat_time("Save as spatRaster - qs2", level = 1)
    Prediction_R <- terra::wrap(Prediction_R)
    IASDT.R::save_as(object = Prediction_R, out_path = Path_Preds_R)

    # Save summary - RData
    IASDT.R::cat_time("Save summary - RData", level = 1)
    IASDT.R::save_as(
      object = Out_Summary,
      object_name = paste0("Prediction_", model_name, "_Summary"),
      out_path = Path_Preds_summary)

    # Save summary - csv
    IASDT.R::cat_time("Save summary - csv", level = 1)
    utils::write.table(
      x = Out_Summary,
      file = stringr::str_replace(Path_Preds_summary, ".RData$", ".txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
      fileEncoding = "UTF-8")
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nThe whole prediction function took ")

  return(Out_Summary)
}
