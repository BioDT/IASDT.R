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
    use_TF = TRUE, TF_environ = NULL, TF_use_single = FALSE,
    LF_n_cores = n_cores, LF_check = FALSE, LF_temp_cleanup = TRUE,
    LF_only = FALSE, LF_commands_only = FALSE, temp_cleanup = TRUE) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----
  #
  IASDT.R::cat_time("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("CV_name", "model_dir"))
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
    stop(
      "Invalid value for CV_name argument. Valid values ",
      "are: 'CV_Dist', 'CV_Large', or `CV_SAC`", call. = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  x <- y <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Loading CV data
  IASDT.R::cat_time("Loading CV data")
  CV_DT <- IASDT.R::path(model_dir, "Model_Fitting_CV", "CV_DT_fitted.RData")
  if (!IASDT.R::check_data(CV_DT, warning = FALSE)) {
    stop("CV data does not exist: ", CV_DT, call. = FALSE)
  }
  CV_DT <- IASDT.R::load_as(CV_DT) %>%
    dplyr::filter(.data$CV_name == CV_name, .data$CV == CV_fold)

  # Loading CV fitted cross-validated model
  IASDT.R::cat_time("Loading CV fitted cross-validated model")
  path_model <- CV_DT$Path_ModFitted
  if (!file.exists(path_model)) {
    stop("Model data does not exist: ", path_model, call. = FALSE)
  }

  # model name
  habitat_abb <- stringr::str_extract(model_dir, "(Hab|hab).+") %>%
    stringr::str_remove("/$") %>%
    stringr::str_remove("Hab|hab")
  model_name <- paste0(habitat_abb, "_", CV_DT$ModName)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading data used in full model (without cross-validation)
  IASDT.R::cat_time("Loading data used in full model")

  model_data <- list.files(
    path = model_dir, pattern = "ModDT_.+_subset.RData", full.names = TRUE)
  if (length(model_data) != 1) {
    stop(
      "There are ", length(model_data), " model data files in the path: ",
      model_dir, call. = FALSE)
  }
  model_data <- IASDT.R::load_as(model_data)

  # path for model evaluation
  Path_Eval <- IASDT.R::path(model_dir, "Model_Fitting_CV", "Model_Evaluation")

  # path for model predictions
  Path_Prediction <- IASDT.R::path(
    model_dir, "Model_Fitting_CV", "Model_Prediction")
  Path_Temp <- IASDT.R::path(
    model_dir, "Model_Fitting_CV", "Temp", paste0("Temp_", model_name))
  fs::dir_create(c(Path_Eval, Path_Prediction, Path_Temp))

  # ..................................................................... ###
  # ..................................................................... ###

  # Predict latent factor at new locations ------

  IASDT.R::cat_time("Predict latent factor at new locations")

  Path_Test_LF <- IASDT.R::path(
    Path_Prediction, paste0("LF_Test_", model_name, ".qs2"))

  if (file.exists(Path_Test_LF)) {

    IASDT.R::cat_time("LF prediction is already available on disk", level = 1)

  } else {

    IASDT.R::cat_time(
      "Preparing input data for predicting latent factor", level = 1)

    Test_DT <- dplyr::filter(model_data$DT_All, .data[[CV_name]] == CV_fold)
    if (nrow(Test_DT) == 0) {
      stop(
        "No data available for the current CV fold: ", CV_name, " - ", CV_fold,
        call. = FALSE)
    }

    Test_X <- dplyr::select(Test_DT, tidyselect::all_of(names(model_data$DT_x)))
    Test_XY <- dplyr::select(Test_DT, x, y)

    Model <- IASDT.R::load_as(path_model)
    Gradient <- Hmsc::prepareGradient(
      hM = Model, XDataNew = as.data.frame(Test_X),
      sDataNew = list(sample = as.data.frame(Test_XY)))

    rm(Model, envir = environment())
    invisible(gc())

    IASDT.R::cat_time("Predicting latent factor", level = 1)
    IASDT.R::cat_sep(
      sep_lines_before = 1, sep_lines_after = 2,
      n_separators = 1, line_char = "*")

    # Predicting latent factor only
    Preds_LF <- IASDT.R::predict_hmsc(
      path_model = path_model, gradient = Gradient, expected = TRUE,
      n_cores = n_cores, model_name = model_name,
      temp_dir = Path_Temp, temp_cleanup = temp_cleanup, use_TF = use_TF,
      TF_environ = TF_environ, LF_out_file = Path_Test_LF,
      TF_use_single = TF_use_single, LF_only = TRUE, LF_n_cores = LF_n_cores,
      LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup,
      LF_commands_only = LF_commands_only, evaluate = FALSE, verbose = TRUE)

    rm(Gradient, Preds_LF, envir = environment())

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

  # # ..................................................................... ###
  # # ..................................................................... ###

  # # Predict_Internal ------
  #
  # Predict_Internal <- function(ID) {
  #
  #   # ID: Index of the current prediction option
  #
  #   # Whether to clamp the sampling efforts
  #   DoClamp <- Prediction_Options$Clamp[[ID]]
  #
  #   # Name of the current option
  #   Option_Name <- Prediction_Options$Name[[ID]]
  #
  #   # Name of the current model
  #   model_name <- paste0(
  #     Option_Name, " - ", dplyr::if_else(DoClamp, "clamping", "no clamping"))
  #
  #   MSG <- paste0(
  #     model_name, " (", ID, "/", nrow(Prediction_Options), ")")
  #   IASDT.R::info_chunk(
  #     MSG, n_separators = 1, line_char = "-", line_char_rep = 70,
  #     cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1,
  #     info_lines_before = 1L)
  #
  #   if (DoClamp) {
  #
  #     # Do not evaluate for options with clamping
  #     evaluate <- FALSE
  #
  #     # Make prediction files at `Path_Prediction_Clamp`
  #     Path_Prediction <- Path_Prediction_Clamp
  #
  #     # use clamped Effort values
  #     StaticPreds <- terra::subset(
  #       x = StaticPredictors, subset = "EffortsLog", negate = TRUE)
  #     StaticPreds$EffortsLog <- StaticPreds$EffortsLog_Clamp
  #     StaticPreds$EffortsLog_Clamp <- NULL
  #
  #   } else {
  #
  #     # evaluate for "Current" climates, without clamping
  #     evaluate <- (Option_Name == "Current")
  #
  #     # Make prediction files at `Path_Prediction_NoClamp`
  #     Path_Prediction <- Path_Prediction_NoClamp
  #
  #     # use original effort data
  #     if ("EffortsLog_Clamp" %in% names(StaticPredictors)) {
  #       StaticPreds <- terra::subset(
  #         x = StaticPredictors, subset = "EffortsLog_Clamp", negate = TRUE)
  #     } else {
  #       StaticPreds <- StaticPredictors
  #     }
  #   }
  #
  #   Path_Prediction_sf <- IASDT.R::path(
  #     Path_Prediction, paste0("Prediction_", Option_Name, "_sf.qs2"))
  #   Path_Prediction_R <- IASDT.R::path(
  #     Path_Prediction, paste0("Prediction_", Option_Name, "_R.qs2"))
  #   Path_Prediction_summary <- IASDT.R::path(
  #     Path_Prediction, paste0("Prediction_", Option_Name, "_Summary.RData"))
  #
  #   # Path for saving tif files of the current option
  #   Path_Prediction_tif <- IASDT.R::path(Path_Prediction, Option_Name)
  #   fs::dir_create(Path_Prediction_tif)
  #   invisible(gc())
  #
  #   # ______________________________________________
  #   # ______________________________________________
  #
  #   # Making predictions if not already processed ----
  #   OutMissing <- file.exists(
  #     c(Path_Prediction_R, Path_Prediction_sf, Path_Prediction_summary)) %>%
  #     all() %>%
  #     isFALSE()
  #
  #   if (OutMissing) {
  #
  #     # Load model
  #     Model <- IASDT.R::load_as(path_model)
  #
  #     # Skip predictions if the predictions as sf object is already on disk
  #     if (file.exists(Path_Prediction_sf)) {
  #
  #       IASDT.R::cat_time("Loading predictions `sf` from disk", level = 1)
  #       Prediction_sf <- IASDT.R::load_as(Path_Prediction_sf)
  #
  #     } else {
  #
  #       .OptionStartTime <- lubridate::now(tzone = "CET")
  #
  #       # Extracting data at training and new sites ------
  #       IASDT.R::cat_time(
  #         "Extracting data at training and new sites", level = 1)
  #       Predict_DF <- Prediction_Options$FilePath[[ID]] %>%
  #         IASDT.R::load_as() %>%
  #         terra::unwrap() %>%
  #         terra::subset(bio_variables) %>%
  #         c(StaticPreds) %>%
  #         terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
  #         tibble::tibble() %>%
  #         sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
  #         sf::st_join(Model_Coords) %>%
  #         tidyr::replace_na(list(Train = FALSE))
  #
  #       # Training locations
  #       Model_Name_Train <- paste0(
  #         Option_Name, "_",
  #         dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Train")
  #       Predict_DF_Train <- dplyr::filter(Predict_DF, Train)
  #       Train_XY <- sf::st_drop_geometry(Predict_DF_Train[, c("x", "y")])
  #       Train_PA <- as.data.frame(Model$Y)
  #       Train_X <- Predict_DF_Train %>%
  #         dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
  #         sf::st_drop_geometry() %>%
  #         stats::model.matrix(Model$XFormula, ., xlev = NULL)
  #
  #       # Testing locations
  #       Model_Name_Test <- paste0(
  #         Option_Name, "_",
  #         dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Test")
  #       Predict_DF_Test <- dplyr::filter(Predict_DF, !Train)
  #       Test_XY <- Predict_DF_Test[, c("x", "y")]
  #       Test_X <- Predict_DF_Test %>%
  #         dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
  #         sf::st_drop_geometry()
  #       Gradient <- Hmsc::prepareGradient(
  #         hM = Model, XDataNew = as.data.frame(Test_X),
  #         sDataNew = list(
  #           sample = as.data.frame(sf::st_drop_geometry(Test_XY))))
  #
  #       rm(
  #         Model, Predict_DF_Test, Predict_DF_Train, Predict_DF,
  #         envir = environment())
  #       invisible(gc())
  #
  #       # ______________________________________________
  #
  #       # Predictions at training sites ----
  #       IASDT.R::cat_time("Predictions at training sites", level = 1)
  #
  #       Path_Current_Train <- IASDT.R::path(
  #         Path_Prediction, paste0("Prediction_", Option_Name, "_Train.qs2"))
  #
  #       if (file.exists(Path_Current_Train)) {
  #
  #         IASDT.R::cat_time("Loading predictions from disk", level = 2)
  #         Preds_ModFitSites <- tibble::tibble(Pred_Path = Path_Current_Train)
  #
  #       } else {
  #
  #         Preds_ModFitSites <- IASDT.R::predict_hmsc(
  #           path_model = path_model, X = Train_X, gradient = NULL,
  #           expected = TRUE, n_cores = n_cores, model_name = Model_Name_Train,
  #           temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_TF = use_TF,
  #           TF_environ = TF_environ, TF_use_single = TF_use_single,
  #           LF_return = TRUE, LF_n_cores = LF_n_cores, LF_check = LF_check,
  #           LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
  #           pred_directory = Path_Prediction, pred_PA = Train_PA,
  #           pred_XY = Train_XY, evaluate = evaluate, evaluation_name = NULL,
  #           evaluation_directory = Path_Eval, verbose = FALSE)
  #
  #       }
  #
  #       # ______________________________________________
  #
  #       # Predictions at new sites ----
  #
  #       if (pred_new_sites) {
  #
  #         IASDT.R::cat_time("Predictions at new sites", level = 1)
  #
  #         Path_Current_Test <- IASDT.R::path(
  #           Path_Prediction, paste0("Prediction_", Option_Name, "_Test.qs2"))
  #
  #         if (file.exists(Path_Current_Test)) {
  #
  #           IASDT.R::cat_time("Loading predictions from disk", level = 2)
  #           Preds_NewSites <- tibble::tibble(Pred_Path = Path_Current_Test)
  #
  #         } else {
  #
  #           Preds_NewSites <- IASDT.R::predict_hmsc(
  #             path_model = path_model, gradient = Gradient, expected = TRUE,
  #             n_cores = n_cores, model_name = Model_Name_Test,
  #             temp_dir = temp_dir, temp_cleanup = temp_cleanup, use_TF = use_TF,
  #             TF_environ = TF_environ, TF_use_single = TF_use_single,
  #             LF_return = TRUE, LF_inputFile = Path_Test_LF,
  #             LF_n_cores = LF_n_cores, LF_check = LF_check,
  #             LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
  #             verbose = FALSE, pred_directory = Path_Prediction,
  #             evaluate = FALSE, pred_XY = sf::st_drop_geometry(Test_XY))
  #
  #         }
  #
  #         # ______________________________________________
  #
  #         # Merge & save predictions - sf ------
  #         IASDT.R::cat_time(
  #           "Merge & save predictions at training and new sites", level = 1)
  #
  #         Prediction_sf <- dplyr::bind_rows(
  #           IASDT.R::load_as(Preds_ModFitSites$Pred_Path),
  #           IASDT.R::load_as(Preds_NewSites$Pred_Path))
  #
  #       } else {
  #
  #         IASDT.R::cat_time(
  #           "Predictions at new sites will NOT be made", level = 1)
  #
  #         # Get column names from predictions at training sites
  #         ColsToAdd <- IASDT.R::load_as(Preds_ModFitSites$Pred_Path) %>%
  #           names() %>%
  #           stringr::str_subset("^Sp_|^SR_|model_name")
  #
  #         # Add missing columns to predictions at new sites
  #         Preds_New_NA <- Test_XY %>%
  #           dplyr::mutate(model_name = Model_Name_Test) %>%
  #           IASDT.R::add_missing_columns(FillVal = NA_real_, ColsToAdd)
  #
  #         # combine predictions
  #         Prediction_sf <- IASDT.R::load_as(Preds_ModFitSites$Pred_Path) %>%
  #           dplyr::bind_rows(Preds_New_NA)
  #
  #         rm(Preds_New_NA, ColsToAdd, envir = environment())
  #       }
  #
  #       # ______________________________________________
  #
  #       # Save predictions as sf object
  #       IASDT.R::save_as(object = Prediction_sf, out_path = Path_Prediction_sf)
  #
  #       IASDT.R::cat_diff(
  #         init_time = .OptionStartTime, prefix = "Prediction took ", level = 1)
  #
  #     }
  #
  #     # ______________________________________________
  #     # ______________________________________________
  #
  #     ### Predictions as spatRaster / tif -----
  #
  #     IASDT.R::cat_time("Rasterization & prepare summary data", level = 1)
  #
  #     Fields2Raster <- names(Prediction_sf) %>%
  #       stringr::str_subset("^Sp_|^SR_") %>%
  #       gtools::mixedsort()
  #
  #     Grid10 <- terra::unwrap(IASDT.R::load_as(Path_GridR))
  #     Prediction_R <- terra::rasterize(
  #       Prediction_sf, Grid10, field = Fields2Raster)
  #
  #     # Calculate prediction anomaly for future projections
  #     if (Option_Name != "Current") {
  #
  #       # Names of current taxa (and SR) for the current model
  #       MeanNames <- stringr::str_subset(Fields2Raster, "_mean")
  #
  #       # Loading mean predictions at current climates
  #       CurrentMean <- list.files(
  #         # use relevant folder containing the current predictions. This is
  #         # determined by `Path_Prediction`, which is not the same whether
  #         # clamping is used or not
  #         path = Path_Prediction, pattern = "Prediction_Current.*_R.qs2",
  #         full.names = TRUE) %>%
  #         IASDT.R::load_as() %>%
  #         terra::unwrap() %>%
  #         terra::subset(MeanNames)
  #
  #       # Calculate anomaly as difference in predicted value between future and
  #       # current climate (future - current)
  #       Preds_Anomaly <- terra::subset(Prediction_R, MeanNames) - CurrentMean
  #       # Assign names to anomaly maps
  #       AnomalyNames <- names(Preds_Anomaly) %>%
  #         stringr::str_replace_all("_mean", "_anomaly")
  #       names(Preds_Anomaly) <- AnomalyNames
  #
  #       # Add anomaly maps to the list of predictions
  #       Fields2Raster <- c(Fields2Raster, AnomalyNames)
  #       Prediction_R <- c(Prediction_R, Preds_Anomaly)
  #
  #       # clean up
  #       rm(Preds_Anomaly, CurrentMean, envir = environment())
  #
  #     }
  #
  #     Out_Summary <- tibble::tibble(
  #       hab_abb = hab_abb, hab_name = Hab_Name,
  #       layer_name = Fields2Raster,
  #       time_period = Prediction_Options$TimePeriod[[ID]],
  #       climate_model = Prediction_Options$ClimateModel[[ID]],
  #       climate_scenario = Prediction_Options$ClimateScenario[[ID]],
  #       Clamp = Prediction_Options$Clamp[[ID]],
  #       Path_Prediction = Path_Prediction) %>%
  #       dplyr::mutate(
  #         Stats = dplyr::case_when(
  #           stringr::str_detect(layer_name, "_mean$") ~ "tif_path_mean",
  #           stringr::str_detect(layer_name, "_sd$") ~ "tif_path_sd",
  #           stringr::str_detect(layer_name, "_cov$") ~ "tif_path_cov",
  #           stringr::str_detect(layer_name, "_anomaly$") ~ "tif_path_anomaly",
  #           .default = NULL),
  #         ias_id = stringr::str_remove(
  #           layer_name, "_mean$|_sd$|_cov$|_anomaly"),
  #         tif_path = IASDT.R::path(
  #           Path_Prediction_tif, paste0(layer_name, ".tif")))
  #
  #     # Save as tif
  #     IASDT.R::cat_time("Save as tif", level = 1)
  #     Out_Summary0 <- Out_Summary %>%
  #       dplyr::mutate(
  #         Map = purrr::map2(
  #           .x = layer_name, .y = tif_path,
  #           .f = ~ {
  #             terra::writeRaster(
  #               x = Prediction_R[[.x]], filename = .y, overwrite = TRUE,
  #               gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  #           }))
  #
  #     rm(Out_Summary0, envir = environment())
  #
  #     Out_Summary <- Out_Summary %>%
  #       tidyr::pivot_wider(
  #         id_cols = c(
  #           "hab_abb", "hab_name", "time_period", "climate_model",
  #           "climate_scenario", "ias_id", "Clamp", "Path_Prediction"),
  #         names_from = "Stats", values_from = "tif_path") %>%
  #       dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
  #       dplyr::select(
  #         tidyselect::any_of(
  #           c(
  #             "hab_abb", "hab_name", "time_period", "climate_model",
  #             "climate_scenario", "ias_id", "taxon_name", "species_name",
  #             "class", "order", "family", "tif_path_mean", "tif_path_sd",
  #             "tif_path_cov", "tif_path_anomaly")))
  #
  #     # save as spatRaster - qs2
  #     IASDT.R::cat_time("Save as spatRaster - qs2", level = 1)
  #     Prediction_R <- terra::wrap(Prediction_R)
  #     IASDT.R::save_as(object = Prediction_R, out_path = Path_Prediction_R)
  #
  #     # Save summary - RData
  #     IASDT.R::cat_time("Save summary - RData", level = 1)
  #     IASDT.R::save_as(
  #       object = Out_Summary, object_name = paste0(Option_Name, "_Summary"),
  #       out_path = Path_Prediction_summary)
  #
  #     # Save summary - csv
  #     IASDT.R::cat_time("Save summary - csv", level = 1)
  #     utils::write.table(
  #       x = Out_Summary,
  #       file = stringr::str_replace(Path_Prediction_summary, ".RData$", ".txt"),
  #       sep = "\t", row.names = FALSE, col.names = TRUE,
  #       quote = FALSE, fileEncoding = "UTF-8")
  #   }
  #
  #   # output
  #   return(
  #     tibble::tibble(
  #       Clamp = DoClamp,
  #       Name = Option_Name,
  #       File_Pred_R = Path_Prediction_R, File_Pred_sf = Path_Prediction_sf,
  #       File_Pred_summary = Path_Prediction_summary))
  # }
  #
  # invisible(gc())
  #
  # # # ..................................................................... ###
  # # # ..................................................................... ###
  #
  # # Predicting ------
  #
  # IASDT.R::info_chunk(
  #   "Making spatial predictions", n_separators = 2, level = 1,
  #   line_char = "*", line_char_rep = 70, cat_red = TRUE,
  #   cat_bold = TRUE, cat_timestamp = FALSE)
  #
  # Grid10 <- terra::unwrap(IASDT.R::load_as(Path_GridR))
  #
  # Prediction_Summary <- purrr::map_dfr(
  #   .x = seq_len(nrow(Prediction_Options)), .f = Predict_Internal) %>%
  #   dplyr::full_join(Prediction_Options, ., by = c("Name", "Clamp")) %>%
  #   dplyr::select(-"FilePath")
  #
  # rm(Predict_Internal, Grid10, Model_Coords, envir = environment())
  # invisible(gc())
  #
  # # # ..................................................................... ###
  # # # ..................................................................... ###
  #
  # # Ensemble model predictions ------
  #
  # IASDT.R::info_chunk(
  #   "\tEnsemble model predictions", n_separators = 1, line_char = "-",
  #   line_char_rep = 70, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)
  #
  # IASDT.R::cat_time("Prepare working in parallel", level = 1)
  #
  # if (n_cores == 1) {
  #   future::plan("future::sequential", gc = TRUE)
  # } else {
  #   withr::local_options(
  #     future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
  #     future.seed = TRUE)
  #   c1 <- snow::makeSOCKcluster(min(n_cores, nrow(Prediction_Summary)))
  #   on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  #   future::plan("future::cluster", workers = c1, gc = TRUE)
  #   withr::defer(future::plan("future::sequential", gc = TRUE))
  # }
  #
  # # --------------------------------------------------------- #
  #
  # # Prepare input data to calculate ensemble predictions
  # IASDT.R::cat_time(
  #   "Prepare input data to calculate ensemble predictions", level = 1)
  #
  # Prediction_Ensemble <- Prediction_Summary %>%
  #   dplyr::filter(ClimateModel != "Current") %>%
  #   dplyr::select(-File_Pred_sf, -File_Pred_R, -Name, -ClimateModel) %>%
  #   dplyr::mutate(
  #     Prediction2 = furrr::future_map(
  #       .x = File_Pred_summary,
  #       .f = ~ {
  #         if (!file.exists(.x)) {
  #           warning("File not found: ", .x, call. = FALSE)
  #           return(NULL)
  #         }
  #
  #         IASDT.R::load_as(.x) %>%
  #           dplyr::select(
  #             -tidyselect::all_of(
  #               c("tif_path_sd", "tif_path_cov", "tif_path_anomaly")))
  #
  #       },
  #       .options = furrr::furrr_options(
  #         seed = TRUE, scheduling = 1,
  #         packages = c("tidyselect", "dplyr", "IASDT.R")))) %>%
  #   dplyr::select("Prediction2", "Clamp") %>%
  #   tidyr::unnest("Prediction2") %>%
  #   dplyr::mutate(
  #     climate_model = "Ensemble",
  #     Dir_Ensemble = IASDT.R::path(
  #       dirname(dirname(tif_path_mean)),
  #       paste0(
  #         stringr::str_replace(time_period, "-", "_"), "_", climate_scenario,
  #         "_Ensemble"))) %>%
  #   dplyr::group_by(dplyr::across(-tif_path_mean)) %>%
  #   dplyr::summarise(tifs = list(tif_path_mean), .groups = "drop") %>%
  #   dplyr::mutate(
  #     tif_path_mean = IASDT.R::path(Dir_Ensemble, paste0(ias_id, "_mean.tif")),
  #     tif_path_anomaly = IASDT.R::path(
  #       Dir_Ensemble, paste0(ias_id, "_anomaly.tif")),
  #     tif_path_sd = IASDT.R::path(Dir_Ensemble, paste0(ias_id, "_sd.tif")),
  #     tif_path_cov = IASDT.R::path(Dir_Ensemble, paste0(ias_id, "_cov.tif")))
  #
  # # --------------------------------------------------------- #
  #
  # IASDT.R::cat_time("Create directories for ensemble predictions", level = 1)
  # fs::dir_create(unique(Prediction_Ensemble$Dir_Ensemble))
  #
  # # --------------------------------------------------------- #
  #
  # # Loading mean predictions at current climates
  # IASDT.R::cat_time("Loading mean predictions at current climates", level = 1)
  #
  # CurrentMean <- list.files(
  #   path = dplyr::if_else(
  #     clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
  #   pattern = "Prediction_Current.*_R.qs2",
  #   full.names = TRUE) %>%
  #   IASDT.R::load_as() %>%
  #   terra::unwrap()
  # CurrentMean <- terra::wrap(CurrentMean["_mean"])
  #
  # # --------------------------------------------------------- #
  #
  # # Calculate ensemble predictions
  # IASDT.R::cat_time("Calculate ensemble predictions", level = 1)
  #
  # Prediction_Ensemble <- Prediction_Ensemble %>%
  #   dplyr::mutate(
  #     Ensemble_Maps = furrr::future_pmap(
  #       .l = list(
  #         tifs, tif_path_mean, tif_path_anomaly,
  #         tif_path_sd, tif_path_cov, ias_id),
  #       .f = function(tifs, tif_path_mean, tif_path_anomaly,
  #                     tif_path_sd, tif_path_cov, ias_id) {
  #
  #         # load maps for future climate option
  #         tiffs_R <- terra::rast(tifs)
  #
  #         # Mean maps
  #         Ensemble_mean <- terra::app(tiffs_R, "mean", na.rm = TRUE) %>%
  #           stats::setNames(paste0(ias_id, "_mean"))
  #         terra::writeRaster(
  #           x = Ensemble_mean, filename = tif_path_mean,
  #           overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  #
  #         # Anomaly maps: calculate prediction anomaly for ensemble future
  #         # projection
  #         CurrentMean0 <- terra::unwrap(CurrentMean) %>%
  #           terra::subset(paste0(ias_id, "_mean"))
  #         Ensemble_anomaly <- (Ensemble_mean - CurrentMean0) %>%
  #           stats::setNames(paste0(ias_id, "_anomaly"))
  #         terra::writeRaster(
  #           x = Ensemble_anomaly, filename = tif_path_anomaly,
  #           overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  #
  #         rm(CurrentMean0, envir = environment())
  #
  #         # Standard deviation
  #         Ensemble_sd <- terra::app(tiffs_R, "sd", na.rm = TRUE) %>%
  #           stats::setNames(paste0(ias_id, "_sd"))
  #         terra::writeRaster(
  #           x = Ensemble_sd, filename = tif_path_sd,
  #           overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  #
  #         # coefficient of variation
  #         Ensemble_mean0 <- Ensemble_mean
  #         # Replace very small mean values with reasonably small number to avoid
  #         # overflow warning
  #         Ensemble_mean0[Ensemble_mean0 < 1e-8] <- 1e-8
  #         Ensemble_cov <- (Ensemble_sd / Ensemble_mean0) %>%
  #           stats::setNames(paste0(ias_id, "_cov"))
  #         terra::writeRaster(
  #           x = Ensemble_cov, filename = tif_path_cov,
  #           overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  #
  #         return(NULL)
  #       },
  #       .options = furrr::furrr_options(
  #         seed = TRUE, scheduling = 1, packages = "terra",
  #         globals = "CurrentMean")),
  #     Dir_Ensemble = NULL, tifs = NULL, Ensemble_Maps = NULL)
  #
  # rm(CurrentMean, envir = environment())
  #
  # # --------------------------------------------------------- #
  #
  # # Save ensemble maps as SpatRast
  # IASDT.R::cat_time("Save ensemble predictions as SpatRast", level = 1)
  #
  # Prediction_Ensemble_R <- Prediction_Ensemble %>%
  #   dplyr::mutate(
  #     Ensemble_File = IASDT.R::path(
  #       dplyr::if_else(
  #         clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
  #       paste0(
  #         "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
  #         climate_scenario, "_Ensemble_R.qs2"))) %>%
  #   dplyr::select(tidyselect::all_of(c(
  #     "ias_id", "Ensemble_File", "tif_path_mean", "tif_path_anomaly",
  #     "tif_path_sd", "tif_path_cov"))) %>%
  #   tidyr::nest(data = -Ensemble_File)
  #
  # Prediction_Ensemble_R <- future.apply::future_lapply(
  #   X = seq_len(nrow(Prediction_Ensemble_R)),
  #   FUN = function(id) {
  #
  #     OutFile <- Prediction_Ensemble_R$Ensemble_File[[id]]
  #
  #     OutMaps <- dplyr::slice(Prediction_Ensemble_R, id) %>%
  #       dplyr::pull("data") %>%
  #       magrittr::extract2(1) %>%
  #       dplyr::mutate(
  #         Maps = purrr::pmap(
  #           .l = list(ias_id, tif_path_mean, tif_path_anomaly,
  #                     tif_path_sd, tif_path_cov),
  #           .f = function(ias_id, tif_path_mean, tif_path_anomaly,
  #                         tif_path_sd, tif_path_cov) {
  #
  #             Mean <- terra::rast(tif_path_mean) %>%
  #               stats::setNames(paste0(ias_id, "_mean"))
  #             SD <- terra::rast(tif_path_sd) %>%
  #               stats::setNames(paste0(ias_id, "_sd"))
  #             COV <- terra::rast(tif_path_cov) %>%
  #               stats::setNames(paste0(ias_id, "_cov"))
  #             Anomaly <- terra::rast(tif_path_anomaly) %>%
  #               stats::setNames(paste0(ias_id, "_anomaly"))
  #
  #             c(Mean, SD, COV, Anomaly) %>%
  #               terra::wrap()
  #           }
  #         )) %>%
  #       dplyr::pull("Maps") %>%
  #       purrr::map(terra::unwrap) %>%
  #       terra::rast() %>%
  #       terra::wrap()
  #
  #     IASDT.R::save_as(object = OutMaps, out_path = OutFile)
  #   },
  #   future.scheduling = Inf, future.seed = TRUE,
  #   future.globals = "Prediction_Ensemble_R",
  #   future.packages = c("terra", "dplyr", "purrr", "IASDT.R", "qs2"))
  #
  # snow::stopCluster(c1)
  # future::plan("future::sequential", gc = TRUE)
  #
  # rm(Prediction_Ensemble_R, envir = environment())
  #
  # # --------------------------------------------------------- #
  #
  # # Save summary of ensemble predictions
  # IASDT.R::cat_time("Save summary of ensemble predictions", level = 1)
  #
  # Prediction_Ensemble_Summary <- Prediction_Ensemble %>%
  #   dplyr::select(-Ensemble_Maps) %>%
  #   dplyr::mutate(
  #     Ensemble_File = IASDT.R::path(
  #       dplyr::if_else(
  #         clamp_pred, Path_Prediction_Clamp, Path_Prediction_NoClamp),
  #       paste0(
  #         "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
  #         climate_scenario, "_Ensemble_Summary.RData"))) %>%
  #   tidyr::nest(Ensemble_DT = -Ensemble_File) %>%
  #   dplyr::mutate(
  #     Ensemble_Save = purrr::map2(
  #       .x = Ensemble_DT, .y = Ensemble_File,
  #       .f = ~ {
  #
  #         IASDT.R::save_as(
  #           object = .x, out_path = .y,
  #           object_name = stringr::str_remove(basename(.y), ".RData"))
  #
  #         utils::write.table(
  #           x = .x, sep = "\t", row.names = FALSE, col.names = TRUE,
  #           file = stringr::str_replace(.y, ".RData", ".txt"),
  #           quote = FALSE, fileEncoding = "UTF-8")
  #
  #       }),
  #     Ensemble_Save = NULL) %>%
  #   tidyr::unnest(Ensemble_DT) %>%
  #   dplyr::select(
  #     File_Pred_summary = Ensemble_File, time_period,
  #     climate_model, climate_scenario) %>%
  #   dplyr::mutate(
  #     Clamp = clamp_pred,
  #     Name = paste0(
  #       stringr::str_replace(time_period, "-", "_"), "_",
  #       climate_scenario, "_Ensemble"),
  #     File_Pred_sf = NA_character_,
  #     File_Pred_R = stringr::str_replace(
  #       File_Pred_summary, "_Summary.RData", "_R.qs2")) %>%
  #   dplyr::distinct()
  #
  # # # ..................................................................... ###
  # # # ..................................................................... ###
  #
  # # Overall summary -----
  #
  # IASDT.R::info_chunk(
  #   "\tPrepare overall summary", n_separators = 1, line_char = "-",
  #   line_char_rep = 70, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)
  #
  # Prediction_Summary <- Prediction_Summary %>%
  #   dplyr::rename(
  #     time_period = TimePeriod, climate_model = ClimateModel,
  #     climate_scenario = ClimateScenario) %>%
  #   dplyr::bind_rows(Prediction_Ensemble_Summary) %>%
  #   dplyr::mutate(hab_abb = hab_abb, hab_name = Hab_Name, .before = 1)
  #
  # save(Prediction_Summary, file = Path_Summary_RData)
  # utils::write.table(
  #   x = Prediction_Summary, sep = "\t", row.names = FALSE, col.names = TRUE,
  #   file = Path_Summary_txt, quote = FALSE, fileEncoding = "UTF-8")
  #
  # # # ..................................................................... ###
  # # # ..................................................................... ###
  #
  # # Overall summary - to be uploaded to the data server; for the Shiny App -----
  #
  # if (clamp_pred) {
  #   Prediction_Summary_Shiny <- dplyr::filter(Prediction_Summary, Clamp)
  # } else {
  #   Prediction_Summary_Shiny <- dplyr::filter(Prediction_Summary, !Clamp)
  # }
  #
  # Prediction_Summary_Shiny <- Prediction_Summary_Shiny$File_Pred_summary %>%
  #   purrr::map(IASDT.R::load_as) %>%
  #   dplyr::bind_rows() %>%
  #   dplyr::distinct() %>%
  #   dplyr::mutate(
  #     dplyr::across(
  #       .cols = tidyselect::all_of(
  #         c("tif_path_mean", "tif_path_sd",
  #           "tif_path_cov", "tif_path_anomaly")),
  #       .fns = ~ {
  #         stringr::str_remove(
  #           string = .x,
  #           pattern = paste0(dirname(Path_Summary_RData_Shiny), "/"))
  #       }))
  #
  # save(Prediction_Summary_Shiny, file = Path_Summary_RData_Shiny)
  # utils::write.table(
  #   x = Prediction_Summary_Shiny, sep = "\t", row.names = FALSE,
  #   col.names = TRUE, file = Path_Summary_txt_Shiny, quote = FALSE,
  #   fileEncoding = "UTF-8")
  #
  # # Create tar file for prediction files
  # if (tar_predictions) {
  #
  #   IASDT.R::cat_time("Create tar file for prediction files", level = 1)
  #
  #   # Directory to save the tar file
  #   TarDir <- dirname(Path_Summary_RData_Shiny)
  #   # Path to the tar file
  #   TarFile <- IASDT.R::path(TarDir, "Predictions.tar")
  #   # List of directories in the prediction folder. All directories will be
  #   # included in the tar file
  #   TarFiles <- list.dirs(
  #     path = TarDir, full.names = FALSE, recursive = FALSE) %>%
  #     paste(collapse = " ") %>%
  #     # Add the summary files to the list
  #     paste(
  #       "Prediction_Summary_Shiny.RData", "Prediction_Summary_Shiny.txt",
  #       collapse = " ")
  #
  #   # Command to create the tar file
  #   Command <- stringr::str_glue(
  #     "cd {fs::path_abs(TarDir)}; tar -cf {basename(TarFile)} -b 2048 \\
  #     {TarFiles}")
  #
  #   # Create tar file
  #   system(Command)
  #
  #   # Change the permission of the tar file
  #   Sys.chmod(TarFile, "755", use_umask = FALSE)
  # }
  #
  # # # ................................................................... ###
  # # # ..................................................................... ###
  #
  # IASDT.R::cat_diff(
  #   init_time = .StartTime, prefix = "\nThe whole prediction function took ")
  #
  # return(Prediction_Summary)
}
