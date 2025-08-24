## |------------------------------------------------------------------------| #
# predict_maps_CV ----
## |------------------------------------------------------------------------| #

#' @export
#' @name predict_maps
#' @rdname predict_maps
#' @order 2
#' @author Ahmed El-Gabbas

predict_maps_CV <- function(
    model_dir = NULL, CV_name = NULL, CV_fold = NULL, n_cores = 8L,
    strategy = "multisession", env_file = ".env", use_TF = TRUE,
    TF_environ = NULL, TF_use_single = FALSE, n_cores_LF = n_cores,
    LF_check = FALSE, LF_temp_cleanup = TRUE, LF_only = FALSE,
    LF_commands_only = FALSE, temp_cleanup = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----
  ecokit::cat_time("Checking input arguments")

  ecokit::check_args(
    args_to_check = c("CV_name", "model_dir"), args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "use_TF", "TF_use_single", "LF_check", "LF_temp_cleanup", "LF_only",
      "LF_commands_only", "temp_cleanup"),
    args_type = "logical")

  CV_name <- .validate_cv_name(CV_name)
  n_cores <- .validate_n_cores(n_cores)
  n_cores_LF <- .validate_n_cores(n_cores_LF)
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- n_cores_LF <- 1L

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  x <- y <- Path_Grid <- ias_id <- taxon_name <- species_name <- sp_type <-
    type <- n_grids_pres <- n_grids_abs <- Sp <- IAS_ID <- ncells <-
    layer_name <- tif_path <- class <- order <- family <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ------

  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Path_GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", Path_GridR = Path_GridR,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Reading input data + managing File and directory paths

  ecokit::cat_time("Reading input data + managing File and directory paths")

  # Loading CV data
  CV_DT <- fs::path(model_dir, "Model_Fitting_CV", "CV_DT_fitted.RData")
  if (!ecokit::check_data(CV_DT, warning = FALSE)) {
    ecokit::stop_ctx(
      "CV data does not exist", CV_DT = CV_DT, include_backtrace = TRUE)
  }

  CV_DT <- ecokit::load_as(CV_DT) %>%
    dplyr::filter(.data$CV_name == CV_name, .data$CV == CV_fold)

  # model name
  hab_abb <- stringr::str_extract(model_dir, "(Hab|hab).+") %>%
    stringr::str_remove("/$") %>%
    stringr::str_remove("Hab|hab")
  # full habitat name
  hab_name <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", hab_abb, "_")) %>%
    stringr::str_remove(paste0("^", hab_abb, "_")) %>%
    stringr::str_replace_all("_", " ")
  # model name
  model_name <- paste0(hab_abb, "_", CV_DT$ModName)

  # path for model evaluation
  Path_Eval <- fs::path(model_dir, "Model_Fitting_CV", "Evaluation")

  # path for model predictions
  Path_Prediction <- fs::path(
    model_dir, "Model_Fitting_CV", "Model_Prediction", model_name)

  temp_dir <- fs::path(
    model_dir, "Model_Fitting_CV", "Temp", paste0("Temp_", model_name))
  fs::dir_create(c(Path_Eval, Path_Prediction, temp_dir))

  Path_Preds_sf <- fs::path(
    Path_Prediction, paste0("Prediction_", model_name, ".qs2"))
  Path_Eval_File <- fs::path(Path_Eval, paste0("Eval_", model_name, ".qs2"))
  Path_Preds_R <- fs::path(
    Path_Prediction, paste0("Prediction_", model_name, "_R.qs2"))
  Path_Preds_summary <- fs::path(
    Path_Prediction, paste0("Prediction_", model_name, "_Summary.RData"))

  # CV fitted cross-validated model
  path_model <- CV_DT$Path_ModFitted
  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model data does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading data used in full model (without cross-validation)
  ecokit::cat_time("Loading data used in full model")

  model_data <- fs::path(model_dir, "ModDT_subset.RData")
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

  Test_DT <- dplyr::filter(model_data$DT_All, .data[[CV_name]] == CV_fold)
  if (nrow(Test_DT) == 0) {
    ecokit::stop_ctx(
      "No data available for the current CV fold",
      CV_name = CV_name, CV_fold = CV_fold, include_backtrace = TRUE)
  }

  Test_X <- Test_DT %>%
    dplyr::select(tidyselect::all_of(names(model_data$DT_x))) %>%
    as.data.frame()
  Test_XY <- as.data.frame(dplyr::select(Test_DT, x, y))

  Model <- ecokit::load_as(path_model)
  Gradient <- Hmsc::prepareGradient(
    hM = Model, XDataNew = Test_X, sDataNew = list(sample = Test_XY))
  Test_PA <- dplyr::select(Test_DT, tidyselect::starts_with("Sp_")) %>%
    as.data.frame()

  n_grid_summary <- Test_PA %>%
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

  rm(Model, Test_X, envir = environment())
  invisible(gc())

  # ..................................................................... ###
  # ..................................................................... ###

  # Predict latent factor at new locations ------

  Path_Test_LF <- fs::path(
    Path_Prediction, paste0("LF_Test_", model_name, ".qs2"))

  if (file.exists(Path_Test_LF)) {

    ecokit::cat_time("LF prediction is already available on disk")

  } else {

    ecokit::info_chunk("Predict latent factor at new locations")

    # Predicting latent factor only
    Preds_LF <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy, model_name = model_name,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_TF = use_TF,
      TF_environ = TF_environ, LF_out_file = Path_Test_LF,
      TF_use_single = TF_use_single, LF_only = TRUE, n_cores_LF = n_cores_LF,
      LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup,
      LF_commands_only = LF_commands_only, evaluate = FALSE, verbose = TRUE)

    rm(Preds_LF, envir = environment())

    ecokit::cat_time("Predicting latent factor is finished!", level = 1L)
    ecokit::cat_sep(
      sep_lines_before = 1L, sep_lines_after = 2L, n_separators = 1L,
      line_char = "*")

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
    ecokit::cat_time("Loading predictions `sf` from disk")
    Prediction_sf <- ecokit::load_as(Path_Preds_sf)

  } else {

    ecokit::info_chunk(
      "Predict habitat suitability at testing cross-validation folds",
      line_char_rep = 75)

    Prediction_sf <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, strategy = strategy, model_name = model_name,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_TF = use_TF,
      TF_environ = TF_environ, TF_use_single = TF_use_single, LF_return = TRUE,
      LF_inputFile = Path_Test_LF, n_cores_LF = n_cores_LF, LF_check = LF_check,
      LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
      verbose = TRUE, pred_directory = Path_Prediction, evaluate = TRUE,
      evaluation_directory = Path_Eval, pred_XY = Test_XY, pred_PA = Test_PA)

    #  --------------------------------------------------------------------- #

    if (!file.exists(Path_Preds_sf)) {
      ecokit::stop_ctx(
        "Prediction file does not exist", Path_Preds_sf = Path_Preds_sf,
        include_backtrace = TRUE)

    }
    if (isFALSE(ecokit::check_data(Path_Preds_sf, warning = FALSE))) {
      ecokit::stop_ctx(
        "Prediction file is corrupted", Path_Preds_sf = Path_Preds_sf,
        include_backtrace = TRUE)
    }

    if (!file.exists(Path_Eval_File)) {
      ecokit::stop_ctx(
        "Evaluation file does not exist", Path_Eval_File = Path_Eval_File,
        include_backtrace = TRUE)
    }

    if (isFALSE(ecokit::check_data(Path_Eval_File, warning = FALSE))) {
      ecokit::stop_ctx(
        "Evaluation file is corrupted", Path_Eval_File = Path_Eval_File,
        include_backtrace = TRUE)
    }

    # loading evaluation data
    Eval_data <- ecokit::load_as(Prediction_sf$Eval_Path) %>%
      dplyr::select(-Sp) %>%
      dplyr::rename(ias_id = IAS_ID)

    # loading prediction data
    Prediction_sf <- ecokit::load_as(Prediction_sf$Pred_Path)

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Rasterizing predictions and making summary ------

  if (all(file.exists(Path_Preds_R, Path_Preds_summary))) {

    Out_Summary <- ecokit::load_as(Path_Preds_summary)

  } else {

    ecokit::cat_time("\nRasterization predictions at testing sites")

    SpeciesInfo <- IASDT.R::get_species_name(env_file = env_file) %>%
      janitor::clean_names() %>%
      dplyr::select(ias_id, taxon_name, species_name, class, order, family)

    Fields2Raster <- names(Prediction_sf) %>%
      stringr::str_subset("^Sp_|^SR_") %>%
      gtools::mixedsort()

    Grid10 <- ecokit::load_as(Path_GridR, unwrap_r = TRUE)
    Prediction_R <- terra::rasterize(
      Prediction_sf, Grid10, field = Fields2Raster)

    #  --------------------------------------------------------------------- #

    ecokit::cat_time("Save prediction outputs")

    # Save as tif
    ecokit::cat_time("Save predictions as tif files", level = 1L)

    selected_columns <- c(
      "hab_abb", "hab_name", "CV_name", "CV_fold", "ias_id", "taxon_name",
      "species_name", "class", "order", "family", "tif_path_mean",
      "tif_path_sd", "tif_path_cov", "n_grids_pres", "n_grids_abs")

    Out_Summary <- tibble::tibble(
      layer_name = Fields2Raster, hab_abb = hab_abb, hab_name = hab_name,
      CV_name = CV_name, CV_fold = CV_fold) %>%
      dplyr::mutate(
        Stats = dplyr::case_when(
          endsWith(layer_name, "_mean") ~ "tif_path_mean",
          endsWith(layer_name, "_sd") ~ "tif_path_sd",
          endsWith(layer_name, "_cov") ~ "tif_path_cov",
          .default = NULL),
        ias_id = stringr::str_remove(layer_name, "_mean$|_sd$|_cov$"),
        tif_path = fs::path(
          Path_Prediction, paste0(layer_name, "_", model_name, ".tif")),
        Out_Map = purrr::map2(
          .x = layer_name, .y = tif_path,
          .f = ~ {
            terra::writeRaster(
              x = Prediction_R[[.x]], filename = .y, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
          }),
        Out_Map = NULL) %>%
      tidyr::pivot_wider(
        id_cols = c("hab_abb", "hab_name", "CV_name", "CV_fold", "ias_id"),
        names_from = "Stats", values_from = "tif_path") %>%
      dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
      dplyr::left_join(n_grid_summary, by = "ias_id") %>%
      dplyr::select(tidyselect::all_of(selected_columns)) %>%
      dplyr::left_join(Eval_data, by = "ias_id") %>%
      dplyr::slice(gtools::mixedorder(ias_id))

    ecokit::cat_time("Save summary data", level = 1L)
    # save as spatRaster - qs2
    ecokit::cat_time("Save as spatRaster - qs2", level = 1L)
    Prediction_R <- terra::wrap(Prediction_R)
    ecokit::save_as(object = Prediction_R, out_path = Path_Preds_R)

    # Save summary - RData
    ecokit::cat_time("Save summary - RData", level = 1L)
    ecokit::save_as(
      object = Out_Summary,
      object_name = paste0("Prediction_", model_name, "_Summary"),
      out_path = Path_Preds_summary)

    # Save summary - csv
    ecokit::cat_time("Save summary - csv", level = 1L)
    utils::write.table(
      x = Out_Summary,
      file = stringr::str_replace(Path_Preds_summary, ".RData$", ".txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
      fileEncoding = "UTF-8")
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nThe whole prediction function took ")

  return(Path_Preds_summary)
}
