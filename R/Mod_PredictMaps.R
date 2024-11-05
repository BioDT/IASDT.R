## |------------------------------------------------------------------------| #
# Predict_Maps ----
## |------------------------------------------------------------------------| #

#' Predicts habitat suitability of `Hmsc` model at different climate scenarios
#'
#' This function generates prediction maps of `Hmsc` models under current and
#' future climates. It also predicts an ensemble predictions for different
#' climate models used. For each species and for species richness, the function
#' exports three maps representing the mean, standard deviation (sd), and
#' coefficient of variation (cov).
#'
#' @param Path_Model Character. Path to fitted `Hmsc` model object.
#' @param Hab_Abb Character. Habitat abbreviation indicating the specific
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type to
#'   prepare data for. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`,
#'   `12a`, `12b`. For more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param NCores Integer specifying the number of parallel cores for
#'   parallelization. Default: 8 cores.
#' @param ChunkSize Integer. The size of chunks used to split the input raster
#'   for parallel processing. Predictor rasters are converted into data frame
#'   (tibble) then split into small chunks for processing on parallel. Default:
#'   `500`.
#' @param RemoveChunks Logical. Whether to remove intermediate chunk files after
#'   prediction.
#' @export
#' @name Predict_Maps
#' @author Ahmed El-Gabbas
#' @inheritParams predictHmsc
#' @return Returns a tibble with the prediction summary and file paths for
#'   output tif files.
#' @export


Predict_Maps <- function(
    Path_Model = NULL, Hab_Abb = NULL, EnvFile = ".env", FromHPC = TRUE,
    NCores = 8,
    Pred_Clamp = c(TRUE, FALSE), FixEfforts = "max",
    TF_Environ = NULL, UseTF = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")
  Hab_Abb <- as.character(Hab_Abb)
  Hab_Name <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_remove(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_replace_all("_", " ")

  # # ..................................................................... ###

  # Check input arguments ----

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c("Hab_Abb", "EnvFile", "Path_Model")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("FromHPC", "RemoveChunks"), Type = "logical")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("NCores", "ChunkSize"), Type = "numeric")

  rm(AllArgs, CharArgs)

  # # ..................................................................... ###

  ValidHabAbbs <- c(as.character(0:3), "4a", "4b", "10", "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0(
        "Invalid Habitat abbreviation. Valid values are:\n >> ",
        paste0(ValidHabAbbs, collapse = ", ")),
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- Path_Roads <- Path_Rail <- Path_Bias <- FilePath <-
    time_period <- climate_model <- climate_scenario <- Path_CHELSA <-
    Path_CLC <- ias_id <- taxon_name <- ClimateModel <- filepath <-
    class <- order <- family <- species_name <- tif_path_mean <-
    tif_path_cov <- tif_path_sd <- hab_abb <- hab_name <-
    Train <- NULL

  # # ..................................................................... ###

  # Environment variables ------

  IASDT.R::CatTime("Environment variables")

  if (!file.exists(EnvFile)) {
    stop(
      paste0("Path for environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Rail", "DP_R_Railways", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts", TRUE, FALSE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_CHELSA", "DP_R_CHELSA_Output", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Rail", "DP_R_Railways_Local", TRUE, FALSE,
      "Path_Roads", "DP_R_Roads_Local", TRUE, FALSE,
      "Path_CLC", "DP_R_CLC_Local", TRUE, FALSE,
      "Path_Bias", "DP_R_Efforts_Local", TRUE, FALSE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_CHELSA", "DP_R_CHELSA_Output_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # Loading input data -----
  IASDT.R::CatTime("Loading Loading input data")

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Species information -----

  SpeciesInfo <- IASDT.R::GetSpeciesName(
    EnvFile = EnvFile, FromHPC = FromHPC) %>%
    janitor::clean_names() %>%
    dplyr::select(ias_id, taxon_name, species_name, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Reference grid -----
  IASDT.R::CatTime("Reference grid", Level = 1)

  Path_GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    stop(
      paste0("Path for the Europe boundaries does not exist: ", Path_GridR),
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Model object -----

  IASDT.R::CatTime("Model object", Level = 1)

  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop(
      paste0("Model path is NULL or does not exist: ", Path_Model),
      call. = FALSE)
  }

  Model <- IASDT.R::LoadAs(Path_Model)

  OtherVars <- stringr::str_subset(names(Model$XData), "^bio|CV", negate = TRUE)
  BioVars <- stringr::str_subset(names(Model$XData), "^bio", negate = FALSE)
  Path_Prediction <- file.path(
    dirname(dirname(Path_Model)), "Model_Prediction")
  Path_Eval <- file.path(dirname(dirname(Path_Model)), "Model_Evaluation")
  fs::dir_create(c(Path_Prediction, Path_Eval))


  if (sum(Pred_Clamp) > 0 && is.null(FixEfforts)) {
    stop("`FixEfforts` can not be NULL when Clamping is implemented")
  }

  if (!is.null(FixEfforts)) {
    Path_Prediction_Clamp <- file.path(
      dirname(dirname(Path_Model)), "Model_Prediction_Clamp")
    fs::dir_create(Path_Prediction_Clamp)
  }

  Model_Coords <- Model$rL$sample$s %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    dplyr::mutate(Train = TRUE)

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## CHELSA data -----

  IASDT.R::CatTime("CHELSA data", Level = 1)

  CHELSA_Data <- file.path(Path_CHELSA, "CHELSA_Processed_DT.RData")
  if (!file.exists(CHELSA_Data)) {
    stop(
      paste0("Processed CHLESA data can not be found at: ", CHELSA_Data),
      call. = FALSE)
  }
  CHELSA_Data <- IASDT.R::LoadAs(CHELSA_Data)

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Road and railway intensity -----

  StaticPredictors <- list()

  if ("RoadRailLog" %in% OtherVars) {

    IASDT.R::CatTime("Road and railway intensity", Level = 1)

    R_Railways <- file.path(Path_Rail, "Railways_Length.RData")
    R_Roads <- file.path(Path_Roads, "Road_Length.RData")

    if (!file.exists(R_Railways)) {
      stop(
        paste0("Railways data does not exist at: ", R_Railways), call. = FALSE)
    }

    if (!file.exists(R_Roads)) {
      stop(
        paste0("Roads data does not exist at: ", R_Roads), call. = FALSE)
    }

    R_Railways <- IASDT.R::LoadAs(R_Railways) %>%
      terra::unwrap() %>%
      magrittr::extract2("rail")

    R_Roads <- IASDT.R::LoadAs(R_Roads) %>%
      terra::unwrap() %>%
      magrittr::extract2("All")

    R_RoadRail <- (R_Roads + R_Railways) %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("RoadRailLog")

    StaticPredictors <- c(StaticPredictors, R_RoadRail)
    rm(R_RoadRail, R_Roads, R_Railways)
  }

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Habitat information ----

  if ("HabLog" %in% OtherVars) {

    IASDT.R::CatTime("Habitat information", Level = 1)

    R_Hab <- file.path(Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
    if (!file.exists(R_Hab)) {
      stop(
        paste0("Habitat data: '", R_Hab, "' does not exist"), call. = FALSE)
    }

    R_Hab <- IASDT.R::LoadAs(R_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb)) %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("HabLog")

    StaticPredictors <- c(StaticPredictors, R_Hab)

    rm(R_Hab)
  }

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Sampling efforts -----

  if ("EffortsLog" %in% OtherVars) {

    IASDT.R::CatTime("Sampling efforts", Level = 1)

    R_Efforts <- file.path(Path_Bias, "Efforts_SummaryR.RData")
    if (!file.exists(R_Efforts)) {
      stop(
        paste0("Sampling efforts data does not exist at: ", R_Efforts),
        call. = FALSE)
    }

    R_Efforts <- IASDT.R::LoadAs(R_Efforts) %>%
      terra::unwrap() %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("EffortsLog")

    # Check if the FixEfforts value is valid
    if (!is.null(FixEfforts)) {
      if (length(FixEfforts) != 1) {
        stop(
          paste0(
            "`FixEfforts` must be a single value. The current value is: ",
            paste0(FixEfforts, collapse = " & ")),
          call. = FALSE)
      }

      if (is.numeric(FixEfforts)) {
        EffortsRange <- terra::global(R_Efforts, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()
        if (!dplyr::between(FixEfforts, EffortsRange[1], EffortsRange[2])) {
          stop(
            paste0(
              "`FixEfforts` value (", FixEfforts, ") is out of the range of ",
              "the observed efforts: From ",
              paste0(round(EffortsRange, 2), collapse = " to ")),
            call. = FALSE)
        }
      } else {
        FixEfforts <- stringr::str_to_lower(FixEfforts)
        if (!(FixEfforts %in% c("median", "mean", "max"))) {
          stop(
            paste0(
              "`FixEfforts` has to be either NULL, single numeric ",
              "value, or one of the following: 'median', 'mean', or 'max'. ",
              "The current value is: ", FixEfforts),
            call. = FALSE)
        }
      }
    }

    # Value to fix efforts at
    if (is.numeric(FixEfforts)) {
      EffortsVal <- FixEfforts
    } else {
      EffortsVal <- dplyr::case_when(
        is.null(FixEfforts) ~ NA_real_,
        FixEfforts == "median" ~ {
          terra::global(
            R_Efforts,
            fun = function(x) {
              median(x, na.rm = TRUE)
            }) %>%
            unlist() %>%
            as.numeric()
        },
        FixEfforts == "mean" ~ {
          terra::global(R_Efforts, fun = mean, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },
        FixEfforts == "max" ~ {
          terra::global(R_Efforts, fun = max, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },
        .default = NA_real_
      )
    }

    # Fix at single value
    if (!is.na(EffortsVal)) {
      R_Efforts_Clamp <- terra::clamp(R_Efforts, EffortsVal, EffortsVal) %>%
        stats::setNames("EffortsLog_Clamp")
      StaticPredictors <- c(StaticPredictors, R_Efforts, R_Efforts_Clamp)
      rm(R_Efforts, R_Efforts_Clamp)
    } else {
      StaticPredictors <- c(StaticPredictors, R_Efforts)
      rm(R_Efforts)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  IASDT.R::CatTime("Merge static predictors", Level = 1)

  StaticPredictors <- terra::rast(StaticPredictors)

  # # ..................................................................... ###

  # Predict latent factor for new locations ------

  Path_Test_LF <- file.path(Path_Prediction, "Test_LF.qs")

  if (!file.exists(Path_Test_LF)) {

    IASDT.R::InfoChunk("Predict latent factor for new sites", Extra2 = 1)

    Model_Name_LF <- paste0("LF_", Hab_Abb, "_Test")

    Predict_DF_Test <- CHELSA_Data %>%
      dplyr::filter(ClimateModel == "Current") %>%
      dplyr::pull(FilePath) %>%
      IASDT.R::LoadAs() %>%
      terra::unwrap() %>%
      terra::subset(BioVars) %>%
      c(StaticPredictors) %>%
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

    rm(Model, Predict_DF_Test, Test_X, Test_XY)
    invisible(gc())

    Preds_Current_Test <- IASDT.R::predictHmsc(
      object = Path_Model,
      Gradient = Gradient,
      expected = TRUE,
      NCores = NCores,
      Model_Name = Model_Name_LF,
      Temp_Dir = "TEMP2Pred",
      UseTF = UseTF,
      TF_Environ = TF_Environ,
      LF_OutFile = Path_Test_LF,
      LF_Only = TRUE,
      Evaluate = FALSE,
      Verbose = FALSE)

    rm(Gradient, Preds_Current_Test)
  } else {
    IASDT.R::CatTime(
      "Latent factor predictions for new sites is already available on disk")
  }

  # # ..................................................................... ###

  # PredictInternal ------

  IASDT.R::CatTime("PredictInternal function")

  PredictInternal <- function(ID) {
    Name <- stringr::str_remove(Prediction_Options$Name[[ID]], "^R_")
    Model_Name <- paste0(Name, "_", Hab_Abb)
    IASDT.R::CatTime(Model_Name, Level = 1)

    FilePath <- Prediction_Options$FilePath[[ID]]
    TimePeriod <- Prediction_Options$TimePeriod[[ID]]
    ClimateModel <- Prediction_Options$ClimateModel[[ID]]
    ClimateScenario <- Prediction_Options$ClimateScenario[[ID]]
    Clamp <- Prediction_Options$Clamp[[ID]]

    if (Clamp) {
      IASDT.R::CatTime(
        "Predictions will be made with fixing sampling efforts at single value",
        Level = 2)
      Path_Pred <- Path_Prediction_Clamp
      Evaluate <- dplyr::if_else(Name == "Current", TRUE, FALSE)
      Path_Prediction_sf <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_Clamp_sf.qs"))
      Path_Prediction_R <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_Clamp_R.qs"))
      Path_Prediction_summary <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_Clamp_Summary.RData"))

      StaticPreds <- StaticPredictors %>%
        terra::subset(subset = "EffortsLog", negate = TRUE)
      StaticPreds$EffortsLog <- StaticPreds$EffortsLog_Clamp
      StaticPreds$EffortsLog_Clamp <- NULL

      Model_Name <- paste0(Model_Name, "_Clamp")

    } else {

      IASDT.R::CatTime(
        "Predictions will be made at original sampling efforts values",
        Level = 2)

      Path_Pred <- Path_Prediction
      Evaluate <- FALSE
      Path_Prediction_sf <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_sf.qs"))
      Path_Prediction_R <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_R.qs"))
      Path_Prediction_summary <- file.path(
        Path_Pred, paste0("Prediction_", Name, "_Summary.RData"))

      StaticPreds <- StaticPredictors %>%
        terra::subset(subset = "EffortsLog_Clamp", negate = TRUE)
    }

    # Path for tif files
    Path_Prediction_tif <- file.path(Path_Pred, Name)
    fs::dir_create(Path_Prediction_tif)

    # Load model
    Model <- IASDT.R::LoadAs(Path_Model)


    # Making predictions if not already processed

    OutMissing <- file.exists(
      c(Path_Prediction_R, Path_Prediction_sf, Path_Prediction_summary)) %>%
      all() %>%
      isFALSE()

    if (OutMissing) {

      # Data at training and new sites ------
      IASDT.R::CatTime("Data at training and new sites", Level = 2)
      Predict_DF <- FilePath %>%
        IASDT.R::LoadAs() %>%
        terra::unwrap() %>%
        terra::subset(BioVars) %>%
        c(StaticPreds) %>%
        terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
        tibble::tibble() %>%
        sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
        sf::st_join(Model_Coords) %>%
        tidyr::replace_na(list(Train = FALSE))

      ## training locations
      Model_Name_Train <- paste0(Model_Name, "_Train")
      Predict_DF_Train <- dplyr::filter(Predict_DF, Train)
      Train_XY <- sf::st_drop_geometry(Predict_DF_Train[, c("x", "y")])
      Train_PA <- as.data.frame(Model$Y)
      Train_X <- Predict_DF_Train %>%
        dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
        sf::st_drop_geometry() %>%
        stats::model.matrix(Model$XFormula, ., xlev = NULL)

      # Testing Locations
      Model_Name_Test <- paste0(Model_Name, "_Test")
      Predict_DF_Test <- dplyr::filter(Predict_DF, !Train)
      Test_XY <- sf::st_drop_geometry(Predict_DF_Test[, c("x", "y")])
      Test_X <- Predict_DF_Test %>%
        dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
        sf::st_drop_geometry() %>%
        stats::model.matrix(Model$XFormula, ., xlev = NULL)
      Gradient <- Hmsc::prepareGradient(
        hM = Model, XDataNew = as.data.frame(Test_X),
        sDataNew = list(sample = as.data.frame(Test_XY)))

      rm(Model, Predict_DF_Test, Predict_DF_Train, Predict_DF)
      invisible(gc())


      # Predictions at training sites ----
      IASDT.R::CatTime("Predictions at training sites", Level = 2)
      Preds_Current_Train <- IASDT.R::predictHmsc(
        object = Path_Model,
        X = Train_X,
        Gradient = NULL,
        expected = TRUE,
        NCores = NCores,
        Model_Name = Model_Name_Train,
        Temp_Dir = "TEMP2Pred",
        UseTF = UseTF,
        TF_Environ = TF_Environ,
        LF_Return = TRUE,
        Pred_Dir = Path_Pred,
        Pred_PA = Train_PA,
        Pred_XY = Train_XY,
        Evaluate = Evaluate,
        Eval_Name = NULL,
        Eval_Dir = Path_Eval,
        Verbose = FALSE) %>%
        dplyr::mutate(Scenario = Model_Name_Train)

      # Predictions at new sites ----
      IASDT.R::CatTime("Predictions at new sites", Level = 2)
      Preds_Current_Test <- IASDT.R::predictHmsc(
        object = Path_Model,
        Gradient = Gradient,
        expected = TRUE,
        NCores = NCores,
        Model_Name = Model_Name_Test,
        Temp_Dir = "TEMP2Pred",
        UseTF = UseTF,
        TF_Environ = TF_Environ,
        LF_Return = TRUE,
        LF_InputFile = Path_Test_LF,
        Pred_Dir = Path_Pred,
        Pred_XY = Test_XY,
        Evaluate = FALSE,
        Verbose = FALSE) %>%
        dplyr::mutate(Scenario = Model_Name_Test)

      # Merge & save predictions - sf ------
      IASDT.R::CatTime(
        "Merge & save predictions at training and new sites", Level = 2)
      Prediction_sf <- dplyr::bind_rows(
        IASDT.R::LoadAs(Preds_Current_Train$Pred_Path),
        IASDT.R::LoadAs(Preds_Current_Test$Pred_Path))
      qs::qsave(Prediction_sf, Path_Prediction_sf, preset = "fast")

      fs::file_delete(
        c(Preds_Current_Train$Pred_Path, Preds_Current_Test$Pred_Path))


      # predictions as spatRaster / tif -----

      # Rasterization & prepare summary table
      IASDT.R::CatTime("rasterization & prepare summary table", Level = 2)

      Prediction_R <- terra::rasterize(
        Prediction_sf, Grid10, field = Fields2Raster)

      Fields2Raster <- names(Prediction_sf) %>%
        stringr::str_subset("^Sp_|^SR_") %>%
        gtools::mixedsort()
      Out_Summary <- tibble::tibble(
        hab_abb = Hab_Abb, hab_name = Hab_Name,
        file = Fields2Raster,
        time_period = TimePeriod, climate_model = ClimateModel,
        climate_scenario = ClimateScenario) %>%
        dplyr::mutate(
          Stats = dplyr::case_when(
            stringr::str_detect(file, "_mean$") ~ "tif_path_mean",
            stringr::str_detect(file, "_sd$") ~ "tif_path_sd",
            stringr::str_detect(file, "_cov$") ~ "tif_path_cov",
            .default = NULL),
          ias_id = stringr::str_remove(file, "_mean$|_sd$|_cov$"),
          filepath = file.path(Path_Prediction_tif, paste0(file, ".tif")))

      # Save as tif
      IASDT.R::CatTime("Save as tif", Level = 2)
      Out_Summary0 <- Out_Summary %>%
        dplyr::mutate(
          Map = purrr::map2(
            .x = file, .y = filepath,
            .f = ~{
              terra::writeRaster(
                x = Prediction_R[[.x]], filename = .y, overwrite = TRUE)
            }))
      rm(Out_Summary0)

      Out_Summary <- Out_Summary %>%
        tidyr::pivot_wider(
          id_cols = c(
            "hab_abb", "hab_name", "time_period", "climate_model",
            "climate_scenario", "ias_id"),
          names_from = "Stats", values_from = "file") %>%
        dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
        dplyr::select(
          hab_abb, hab_name, time_period, climate_model, climate_scenario,
          ias_id, taxon_name, species_name, class, order, family,
          tif_path_mean, tif_path_sd, tif_path_cov)

      # save as spatRaster - qs
      IASDT.R::CatTime("Save as spatRaster - qs", Level = 2)
      Prediction_R <- terra::wrap(Prediction_R)
      qs::qsave(Prediction_R, Path_Prediction_R, preset = "fast")

      # Save summary - RData
      IASDT.R::CatTime("Save summary - RData", Level = 2)
      IASDT.R::SaveAs(
        InObj = Out_Summary, OutObj = paste0(Name, "_Summary"),
        OutPath = Path_Prediction_summary)

      # Save summary - csv
      IASDT.R::CatTime("Save summary - csv", Level = 2)
      utils::write.table(
        x = Out_Summary,
        file = stringr::str_replace(Path_Prediction_summary, ".RData$", ".txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = FALSE, fileEncoding = "UTF-8")
    }

    # output
    tibble::tibble(
      Name = Name,
      Pred_R = Path_Prediction_R,
      Pred_sf = Path_Prediction_sf,
      Pred_summary = Path_Prediction_summary) %>%
      return()
  }


  # Predicting ------
  IASDT.R::InfoChunk("Making predictions", Extra2 = 1)

  Grid10 <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))
  Prediction_Options <- tidyr::expand_grid(CHELSA_Data, Clamp = Pred_Clamp)

  Prediction_Summary <- purrr::map_dfr(
    seq_len(nrow(Prediction_Options)), PredictInternal)


  Prediction_Summary



  # ## CHECK -----
  #
  # withr::local_options(
  #   future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)
  #
  # if (NCores == 1) {
  #   future::plan("future::sequential", gc = TRUE)
  # } else {
  #   c1 <- snow::makeSOCKcluster(min(NCores, nrow(Prediction)))
  #   on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  #   future::plan("future::cluster", workers = c1, gc = TRUE)
  #   on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  # }
  #
  # Prediction <- dplyr::mutate(
  #   Prediction,
  #   Prediction2 = furrr::future_map(
  #     .x  = path_prediction,
  #     .f = ~{
  #       if (file.exists(.x)) {
  #         IASDT.R::LoadAs(.x) %>%
  #           dplyr::select(
  #             -tidyselect::all_of(c(
  #               "time_period", "climate_model", "climate_scenario",
  #               "Prediction_mean", "Prediction_sd", "Prediction_cov")))
  #       } else {
  #         warning(paste0("File not found: ", .x))
  #         return(NULL)
  #       }
  #     },
  #     .options = furrr::furrr_options(
  #       seed = TRUE, scheduling = Inf,
  #       packages = c("tidyselect", "dplyr", "IASDT.R", "terra")
  #     ))) %>%
  #   tidyr::unnest("Prediction2")
  #
  # # # ................................................................... ###
  #
  # # Ensemble of climate models ------
  #
  # IASDT.R::CatTime("\nEnsemble of climate models", ... = "\n")
  #
  # # Create ensemble paths
  # expand.grid(
  #   Time = c("2011_2040", "2041_2070", "2071_2100"),
  #   Scenario = c("ssp126", "ssp370", "ssp585")) %>%
  #   dplyr::mutate(Folder = paste0("Ensemble_", Time, "_", Scenario)) %>%
  #   dplyr::pull("Folder") %>%
  #   file.path(Path_Prediction, .) %>%
  #   fs::dir_create()
  #
  # Ensemble <- dplyr::filter(Prediction, time_period != "1981-2010") %>%
  #   dplyr::select(
  #     tidyselect::all_of(c(
  #       "time_period", "climate_scenario", "ias_id",
  #       "tif_path_mean", "hab_abb", "hab_name"))) %>%
  #   dplyr::mutate(climate_model = "Ensemble", .before = "ias_id") %>%
  #   dplyr::group_by(
  #     time_period, climate_scenario, climate_model,
  #     ias_id, hab_abb, hab_name) %>%
  #   dplyr::summarize(Tif_Path = list(tif_path_mean), .groups = "keep") %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(
  #     Data_Ensemble = furrr::future_pmap(
  #       .l = list(ias_id, time_period, climate_scenario, Tif_Path),
  #       .f = function(ias_id, time_period, climate_scenario, Tif_Path) {
  #
  #         Name_Ensemble <- paste0(
  #           "Ensemble_", time_period, "_", climate_scenario) %>%
  #           stringr::str_replace_all("-", "_")
  #
  #         Maps <- terra::rast(purrr::map(Tif_Path, terra::rast))
  #
  #         Map_mean <- terra::app(Maps, fun = "mean")
  #         Path_mean <- file.path(
  #           Path_Prediction, Name_Ensemble, paste0(ias_id, "_mean.tif"))
  #         terra::writeRaster(
  #           x = Map_mean, overwrite = TRUE, filename = Path_mean)
  #
  #         Map_sd <- terra::app(Maps, fun = "sd")
  #         Path_sd <- file.path(
  #           Path_Prediction, Name_Ensemble, paste0(ias_id, "_sd.tif"))
  #         terra::writeRaster(x = Map_sd, overwrite = TRUE, filename = Path_sd)
  #
  #         Map_cov <- Map_sd / Map_mean
  #         Path_cov <- file.path(
  #           Path_Prediction, Name_Ensemble, paste0(ias_id, "_cov.tif"))
  #         terra::writeRaster(x = Map_cov, overwrite = TRUE,
  # filename = Path_cov)
  #
  #         return(list(
  #           tif_path_mean = Path_mean, tif_path_sd = Path_sd,
  #           tif_path_cov = Path_cov))
  #       },
  #       .options = furrr::furrr_options(
  #         seed = TRUE, scheduling = Inf, globals = "Path_Prediction",
  #         packages = c("purrr", "IASDT.R", "terra")))) %>%
  #   dplyr::select(-"Tif_Path") %>%
  #   tidyr::unnest_wider("Data_Ensemble") %>%
  #   dplyr::left_join(
  #     dplyr::distinct(
  #       Prediction, time_period, climate_scenario, ias_id, taxon_name,
  #       class, order, family, species_name),
  #     by = c("time_period", "climate_scenario", "ias_id"))
  #
  # if (NCores > 1) {
  #   snow::stopCluster(c1)
  #   future::plan("future::sequential", gc = TRUE)
  # }
  #
  # Prediction_Summary <- dplyr::bind_rows(Prediction, Ensemble) %>%
  #   dplyr::select(
  #     hab_abb, hab_name, time_period, climate_model, climate_scenario, ias_id,
  #     taxon_name, species_name, class, order, family, tif_path_mean,
  #     tif_path_sd, tif_path_cov)
  #

  # # ................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPrediction took ")

  return(Prediction_Summary)
}
