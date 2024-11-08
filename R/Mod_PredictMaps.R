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
#' @param Pred_Clamp Logical indicating whether to clamp the sampling efforts at
#'   a single value. Defaults to `TRUE`. If `TRUE`, the `FixEfforts` argument
#'   must be provided.
#' @param FixEfforts Numeric or character. Value to fix the sampling efforts at.
#'   If numeric, the value will be used as is. If character, it can be one of
#'   the following: `median`, `mean`, or `max`. Default: `mean`. This argument
#'   is required when `Pred_Clamp` is `TRUE`.
#' @param Pred_NewSites Logical indicating whether to predict habitat
#'   suitability at new sites. Default: `TRUE`. This argument is only temporary
#'   to make predictions for the workflow and will be removed in the future.
#' @param CC_Models Character vector specifying the climate models to use for
#'   future predictions. Default: `c("GFDL-ESM4", "IPSL-CM6A-LR" ,
#'   "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")`. This argument is only
#'   temporary to make predictions for the workflow and will be removed in the
#'   future.
#' @param CC_Scenario Character vector specifying the climate scenarios to use
#'   for future predictions. Default: `c("ssp126", "ssp370", "ssp585")`. This
#'   argument is only temporary to make predictions for the workflow and will be
#'   removed in the future.
#' @export
#' @name Predict_Maps
#' @author Ahmed El-Gabbas
#' @inheritParams predictHmsc
#' @return Returns a tibble with the prediction summary and file paths for
#'   output tif files.
#' @export

Predict_Maps <- function(
    Path_Model = NULL, Hab_Abb = NULL, EnvFile = ".env", FromHPC = TRUE,
    NCores = 8, Pred_Clamp = TRUE, FixEfforts = "mean", Pred_NewSites = TRUE,
    TF_Environ = NULL, UseTF = TRUE,
    CC_Models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_Scenario = c("ssp126", "ssp370", "ssp585")) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  Hab_Abb <- as.character(Hab_Abb)
  Hab_Name <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_remove(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_replace_all("_", " ")

  # # ..................................................................... ###
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
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- Path_Roads <- Path_Rail <- Path_Bias <- tif_path <-
    time_period <- climate_model <- climate_scenario <- Path_CHELSA <-
    Path_CLC <- ias_id <- taxon_name <- ClimateModel <- ClimateScenario <-
    Name <- File_Pred_sf <- class <- order <- family <- species_name <-
    tif_path_mean <- tif_path_cov <- tif_path_sd <- hab_abb <- hab_name <-
    Train <- Ensemble_File <- Ensemble_Maps <- tifs <- layer_name <-
    TimePeriod <- File_Pred_summary <- Ensemble_DT <- Path_Ensemble <-
    File_Pred_R <- NULL

  # # ..................................................................... ###
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
  # # ..................................................................... ###

  # Loading input data -----
  IASDT.R::CatTime("Loading Loading input data")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Species information -----

  SpeciesInfo <- IASDT.R::GetSpeciesName(
    EnvFile = EnvFile, FromHPC = FromHPC) %>%
    janitor::clean_names() %>%
    dplyr::select(ias_id, taxon_name, species_name, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Reference grid -----

  IASDT.R::CatTime("Reference grid", Level = 1)
  Path_GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    stop(
      paste0("Path for the Europe boundaries does not exist: ", Path_GridR),
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

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

  Path_Prediction <- dplyr::if_else(
    Pred_Clamp, "Model_Prediction_Clamp", "Model_Prediction") %>%
    file.path(dirname(dirname(Path_Model)), .)
  Path_Eval <- file.path(dirname(dirname(Path_Model)), "Model_Evaluation")
  fs::dir_create(c(Path_Prediction, Path_Eval))

  if (Pred_Clamp && is.null(FixEfforts)) {
    stop("`FixEfforts` can not be NULL when Clamping is implemented")
  }

  Model_Coords <- Model$rL$sample$s %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    dplyr::mutate(Train = TRUE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## CHELSA data -----

  IASDT.R::CatTime("CHELSA data", Level = 1)

  Path_CHELSA <- file.path(Path_CHELSA, "CHELSA_Processed_DT.RData")
  if (!file.exists(Path_CHELSA)) {
    stop(
      paste0("Processed CHLESA data can not be found at: ", Path_CHELSA),
      call. = FALSE)
  }

  Prediction_Options <- IASDT.R::LoadAs(Path_CHELSA) %>%
    dplyr::select(-"File_List") %>%
    dplyr::filter(
      # filter only for the following future climate models
      ClimateModel %in% c("Current", CC_Models),
      # filter only for the following future climate scenarios
      ClimateScenario %in% c("Current", CC_Scenario)) %>%
    dplyr::mutate(Name = stringr::str_remove(Name, "^R_"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

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
            R_Efforts, fun = function(x) median(x, na.rm = TRUE)) %>%
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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  IASDT.R::CatTime("Merge static predictors", Level = 1)

  StaticPredictors <- terra::rast(StaticPredictors)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict latent factor for new locations ------

  Path_Test_LF <- file.path(Path_Prediction, "Test_LF.qs")

  if (!file.exists(Path_Test_LF) && Pred_NewSites) {

    IASDT.R::CatTime("Predict latent factor for new sites")

    Predict_DF_Test <- Prediction_Options %>%
      dplyr::filter(ClimateModel == "Current") %>%
      dplyr::pull("FilePath") %>%
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

    # Predicting latent factor
    Preds_Current_Test <- IASDT.R::predictHmsc(
      object = Path_Model, Gradient = Gradient, expected = TRUE,
      NCores = NCores, Model_Name = paste0("LF_", Hab_Abb, "_Test"),
      Temp_Dir = "TEMP2Pred", UseTF = UseTF, TF_Environ = TF_Environ,
      LF_OutFile = Path_Test_LF, LF_Only = TRUE, Evaluate = FALSE,
      Verbose = FALSE)

    rm(Gradient, Preds_Current_Test)

  } else {
    if (Pred_NewSites) {
      IASDT.R::CatTime(
        "Latent factor prediction for new sites is already available on disk")
    } else {
      IASDT.R::CatTime(
        "Latent factor prediction for new sites will NOT be made")
    }
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # PredictInternal ------

  IASDT.R::CatTime("PredictInternal function")

  PredictInternal <- function(ID) {

    Option_Name <- Prediction_Options$Name[[ID]]
    Model_Name <- paste0(Option_Name, "_", Hab_Abb)
    cat("\n")
    IASDT.R::CatTime(Model_Name, NLines = 2)

    if (Pred_Clamp) {

      Evaluate <- dplyr::if_else(Option_Name == "Current", TRUE, FALSE)
      Path_Prediction_sf <- file.path(
        Path_Prediction, paste0("Prediction_", Option_Name, "_Clamp_sf.qs"))
      Path_Prediction_R <- file.path(
        Path_Prediction, paste0("Prediction_", Option_Name, "_Clamp_R.qs"))
      Path_Prediction_summary <- file.path(
        Path_Prediction,
        paste0("Prediction_", Option_Name, "_Clamp_Summary.RData"))

      # use clamped Effort values
      StaticPreds <- terra::subset(
        x = StaticPredictors, subset = "EffortsLog", negate = TRUE)
      StaticPreds$EffortsLog <- StaticPreds$EffortsLog_Clamp
      StaticPreds$EffortsLog_Clamp <- NULL

    } else {

      Evaluate <- FALSE
      Path_Prediction_sf <- file.path(
        Path_Prediction, paste0("Prediction_", Option_Name, "_sf.qs"))
      Path_Prediction_R <- file.path(
        Path_Prediction, paste0("Prediction_", Option_Name, "_R.qs"))
      Path_Prediction_summary <- file.path(
        Path_Prediction, paste0("Prediction_", Option_Name, "_Summary.RData"))

      # use original effort data
      StaticPreds <- terra::subset(
        x = StaticPredictors, subset = "EffortsLog_Clamp", negate = TRUE)

    }

    # Path for saving tif files of the current option
    Path_Prediction_tif <- file.path(Path_Prediction, Option_Name)
    fs::dir_create(Path_Prediction_tif)

    # Making predictions if not already processed ----
    OutMissing <- file.exists(
      c(Path_Prediction_R, Path_Prediction_sf, Path_Prediction_summary)) %>%
      all() %>%
      isFALSE()

    if (OutMissing) {

      # Load model
      Model <- IASDT.R::LoadAs(Path_Model)

      if (!file.exists(Path_Prediction_sf)) {

        # Extracting data at training and new sites ------
        IASDT.R::CatTime("Extracting data at training and new sites", Level = 1)
        Predict_DF <- Prediction_Options$FilePath[[ID]] %>%
          IASDT.R::LoadAs() %>%
          terra::unwrap() %>%
          terra::subset(BioVars) %>%
          c(StaticPreds) %>%
          terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
          tibble::tibble() %>%
          sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
          sf::st_join(Model_Coords) %>%
          tidyr::replace_na(list(Train = FALSE))

        # training locations
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
        Test_XY <- Predict_DF_Test[, c("x", "y")]
        Test_X <- Predict_DF_Test %>%
          dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
          sf::st_drop_geometry()
        Gradient <- Hmsc::prepareGradient(
          hM = Model, XDataNew = as.data.frame(Test_X),
          sDataNew = list(
            sample = as.data.frame(sf::st_drop_geometry(Test_XY))))

        rm(Model, Predict_DF_Test, Predict_DF_Train, Predict_DF)
        invisible(gc())


        # Predictions at training sites ----
        IASDT.R::CatTime("Predictions at training sites", Level = 1)

        Path_Current_Train <- file.path(
          Path_Prediction, paste0("Prediction_", Model_Name_Train, ".qs"))

        if (file.exists(Path_Current_Train)) {
          IASDT.R::CatTime("Loading predictions from disk", Level = 2)
          Preds_Current_Train <- tibble::tibble(
            Pred_Path = Path_Current_Train)
        } else {
          Preds_Current_Train <- IASDT.R::predictHmsc(
            object = Path_Model, X = Train_X, Gradient = NULL, expected = TRUE,
            NCores = NCores, Model_Name = Model_Name_Train,
            Temp_Dir = "TEMP2Pred", UseTF = UseTF, TF_Environ = TF_Environ,
            LF_Return = TRUE, Pred_Dir = Path_Prediction, Pred_PA = Train_PA,
            Pred_XY = Train_XY, Evaluate = Evaluate, Eval_Name = NULL,
            Eval_Dir = Path_Eval, Verbose = FALSE)
        }

        # Predictions at new sites ----

        if (Pred_NewSites) {

          IASDT.R::CatTime("Predictions at new sites", Level = 1)

          Path_Current_Test <- file.path(
            Path_Prediction, paste0("Prediction_", Model_Name_Test, ".qs"))

          if (file.exists(Path_Current_Test)) {
            IASDT.R::CatTime("Loading predictions from disk", Level = 2)
            Preds_Current_Test <- tibble::tibble(
              Pred_Path = Path_Current_Test)

          } else {

            Preds_Current_Test <- IASDT.R::predictHmsc(
              object = Path_Model, Gradient = Gradient, expected = TRUE,
              NCores = NCores, Model_Name = Model_Name_Test,
              Temp_Dir = "TEMP2Pred", UseTF = UseTF, TF_Environ = TF_Environ,
              LF_Return = TRUE, LF_InputFile = Path_Test_LF,
              Pred_Dir = Path_Prediction,
              Pred_XY = sf::st_drop_geometry(Test_XY),
              Evaluate = FALSE, Verbose = FALSE)
          }

          # Merge & save predictions - sf ------
          IASDT.R::CatTime(
            "Merge & save predictions at training and new sites", Level = 1)
          Prediction_sf <- dplyr::bind_rows(
            IASDT.R::LoadAs(Preds_Current_Train$Pred_Path),
            IASDT.R::LoadAs(Preds_Current_Test$Pred_Path))

          try(
            fs::file_delete(
              c(Preds_Current_Train$Pred_Path, Preds_Current_Test$Pred_Path)),
            silent = TRUE)

        } else {

          IASDT.R::CatTime(
            "Predictions at new sites will NOT be made", Level = 1)

          # Get column names from predictions at training sites
          ColsToAdd <- IASDT.R::LoadAs(Preds_Current_Train$Pred_Path) %>%
            names() %>%
            stringr::str_subset("^Sp_|^SR_|Model_Name")

          # Add missing columns to predictions at new sites
          Preds_New_NA <- Test_XY %>%
            dplyr::mutate(Model_Name = Model_Name_Test) %>%
            IASDT.R::AddMissingCols(FillVal = NA_real_, ColsToAdd)

          # combine predictions
          Prediction_sf <- IASDT.R::LoadAs(Preds_Current_Train$Pred_Path) %>%
            dplyr::bind_rows(Preds_New_NA)

          try(fs::file_delete(Preds_Current_Train$Pred_Path), silent = TRUE)
          rm(Preds_New_NA, ColsToAdd)
        }

        # Save predictions as sf object
        qs::qsave(Prediction_sf, Path_Prediction_sf, preset = "fast")

        try(
          fs::file_delete(
            c(Preds_Current_Train$Pred_Path, Preds_Current_Test$Pred_Path)),
          silent = TRUE)

      } else {

        IASDT.R::CatTime("Reading predictions `sf` from disk", Level = 1)
        Prediction_sf <- IASDT.R::LoadAs(Path_Prediction_sf)

        try(
          fs::file_delete(Preds_Current_Train$Pred_Path),
          silent = TRUE)
      }


      # Predictions as spatRaster / tif -----

      IASDT.R::CatTime("Rasterization & prepare summary data", Level = 1)

      Fields2Raster <- names(Prediction_sf) %>%
        stringr::str_subset("^Sp_|^SR_") %>%
        gtools::mixedsort()

      Prediction_R <- terra::rasterize(
        Prediction_sf, Grid10, field = Fields2Raster)

      Out_Summary <- tibble::tibble(
        hab_abb = Hab_Abb, hab_name = Hab_Name,
        layer_name = Fields2Raster,
        time_period = Prediction_Options$TimePeriod[[ID]],
        climate_model = Prediction_Options$ClimateModel[[ID]],
        climate_scenario = Prediction_Options$ClimateScenario[[ID]]) %>%
        dplyr::mutate(
          Stats = dplyr::case_when(
            stringr::str_detect(layer_name, "_mean$") ~ "tif_path_mean",
            stringr::str_detect(layer_name, "_sd$") ~ "tif_path_sd",
            stringr::str_detect(layer_name, "_cov$") ~ "tif_path_cov",
            .default = NULL),
          ias_id = stringr::str_remove(layer_name, "_mean$|_sd$|_cov$"),
          tif_path = file.path(Path_Prediction_tif, paste0(layer_name, ".tif")))

      # Save as tif
      IASDT.R::CatTime("Save maps as tif", Level = 1)
      Out_Summary0 <- Out_Summary %>%
        dplyr::mutate(
          Map = purrr::map2(
            .x = layer_name, .y = tif_path,
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
          names_from = "Stats", values_from = "tif_path") %>%
        dplyr::left_join(SpeciesInfo, by = "ias_id") %>%
        dplyr::select(
          hab_abb, hab_name, time_period, climate_model, climate_scenario,
          ias_id, taxon_name, species_name, class, order, family,
          tif_path_mean, tif_path_sd, tif_path_cov)

      # save as spatRaster - qs
      IASDT.R::CatTime("Save maps as spatRaster - qs", Level = 1)
      Prediction_R <- terra::wrap(Prediction_R)
      qs::qsave(Prediction_R, Path_Prediction_R, preset = "fast")

      # Save summary - RData
      IASDT.R::CatTime("Save summary - RData", Level = 1)
      IASDT.R::SaveAs(
        InObj = Out_Summary, OutObj = paste0(Option_Name, "_Summary"),
        OutPath = Path_Prediction_summary)

      # Save summary - csv
      IASDT.R::CatTime("Save summary - csv", Level = 1)
      utils::write.table(
        x = Out_Summary,
        file = stringr::str_replace(Path_Prediction_summary, ".RData$", ".txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = FALSE, fileEncoding = "UTF-8")
    }

    # output
    tibble::tibble(
      Name = Option_Name,
      File_Pred_R = Path_Prediction_R, File_Pred_sf = Path_Prediction_sf,
      File_Pred_summary = Path_Prediction_summary) %>%
      return()
  }

  # # ..................................................................... ###

  # Predicting ------
  IASDT.R::InfoChunk("Making predictions", Extra2 = 1)

  Grid10 <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))
  Prediction_Summary <- purrr::map_dfr(
    seq_len(nrow(Prediction_Options)), PredictInternal) %>%
    dplyr::full_join(Prediction_Options, ., by = "Name") %>%
    dplyr::select(-"FilePath")

  # # ..................................................................... ###

  # Ensemble model predictions ------

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(Prediction_Summary)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }


  Prediction_Ensemble <- Prediction_Summary %>%
    dplyr::filter(ClimateModel != "Current")  %>%
    dplyr::select(-File_Pred_sf, -File_Pred_R, -Name, -ClimateModel) %>%
    dplyr::mutate(
      Prediction2 = furrr::future_map(
        .x  = File_Pred_summary,
        .f = ~{
          if (file.exists(.x)) {
            IASDT.R::LoadAs(.x) %>%
              dplyr::select(
                -tidyselect::all_of(c("tif_path_sd", "tif_path_cov")))
          } else {
            warning(paste0("File not found: ", .x))
            return(NULL)
          }
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = 1,
          packages = c("tidyselect", "dplyr", "IASDT.R")))) %>%
    dplyr::select("Prediction2") %>%
    tidyr::unnest("Prediction2") %>%
    dplyr::mutate(
      climate_model = "Ensemble",
      Path_Ensemble = file.path(
        dirname(dirname(tif_path_mean)),
        paste0(time_period, "_Ensemble_", climate_scenario))) %>%
    dplyr::group_by(dplyr::across(-tif_path_mean)) %>%
    dplyr::summarise(tifs = list(tif_path_mean), .groups = "drop") %>%
    dplyr::mutate(
      tif_path_mean = file.path(Path_Ensemble, paste0(ias_id, "_mean.tif")),
      tif_path_sd = file.path(Path_Ensemble, paste0(ias_id, "_sd.tif")),
      tif_path_cov = file.path(Path_Ensemble, paste0(ias_id, "_cov.tif")))

  fs::dir_create(unique(Prediction_Ensemble$Path_Ensemble))

  Prediction_Ensemble <- Prediction_Ensemble %>%
    dplyr::mutate(
      Ensemble_Maps = furrr::future_pmap(
        .l = list(tifs, tif_path_mean, tif_path_sd, tif_path_cov, ias_id),
        .f = function(tifs, tif_path_mean, tif_path_sd, tif_path_cov, ias_id) {

          tiffs_R <- terra::rast(tifs)

          # Mean maps
          Ensemble_mean <- terra::app(tiffs_R, "mean", na.rm = TRUE) %>%
            stats::setNames(paste0(ias_id, "_mean"))
          terra::writeRaster(
            x = Ensemble_mean, filename = tif_path_mean, overwrite = TRUE)

          # Standard deviation
          Ensemble_sd <- terra::app(tiffs_R, "sd", na.rm = TRUE) %>%
            stats::setNames(paste0(ias_id, "_sd"))
          terra::writeRaster(
            x = Ensemble_sd, filename = tif_path_sd, overwrite = TRUE)

          # coefficient of variation
          Ensemble_cov <- (Ensemble_sd / Ensemble_mean) %>%
            stats::setNames(paste0(ias_id, "_cov"))
          terra::writeRaster(
            x = Ensemble_cov, filename = tif_path_cov, overwrite = TRUE)

          c(Ensemble_mean, Ensemble_sd, Ensemble_cov) %>%
            terra::wrap() %>%
            return()
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = 1, packages = "terra")),
      Path_Ensemble = NULL, tifs = NULL)

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)


  # Save ensemble maps as SpatRast
  Prediction_Ensemble_R <- Prediction_Ensemble %>%
    dplyr::mutate(
      Ensemble_File = file.path(
        Path_Prediction,
        paste0(
          "Prediction_", time_period, "_Ensemble_",
          climate_scenario, "_R.qs"))) %>%
    dplyr::select(Ensemble_File, Ensemble_Maps) %>%
    dplyr::group_by(Ensemble_File) %>%
    tidyr::nest(Ensemble_Maps = Ensemble_Maps) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Ensemble_Maps = purrr::map(
        .x = Ensemble_Maps,
        .f = ~{
          purrr::map(.x$Ensemble_Maps, terra::unwrap) %>%
            purrr::map(c) %>%
            terra::rast() %>%
            terra::wrap()
        })) %>%
    dplyr::mutate(
      Save = purrr::map2(
        .x = Ensemble_Maps, .y = Ensemble_File,
        .f = qs::qsave, preset = "fast"))
  rm(Prediction_Ensemble_R)

  # Save summary of ensemble maps
  Prediction_Ensemble_Summary <- Prediction_Ensemble %>%
    dplyr::select(-Ensemble_Maps) %>%
    dplyr::mutate(
      Ensemble_File = file.path(
        Path_Prediction, paste0(
          "Prediction_", time_period, "_Ensemble_",
          climate_scenario, "_Summary.RData"))) %>%
    tidyr::nest(Ensemble_DT = -Ensemble_File) %>%
    dplyr::mutate(
      Ensemble_Save = purrr::map2(
        .x = Ensemble_DT, .y = Ensemble_File,
        .f = ~ {
          save(.x, file = .y)
          utils::write.table(
            x = .x, sep = "\t", row.names = FALSE, col.names = TRUE,
            file = stringr::str_replace(.y, ".RData", ".txt"),
            quote = FALSE, fileEncoding = "UTF-8")
        }),
      Ensemble_Save = NULL) %>%
    tidyr::unnest(Ensemble_DT) %>%
    dplyr::select(
      File_Pred_summary = Ensemble_File, time_period,
      climate_model, climate_scenario) %>%
    dplyr::mutate(
      Name = paste0(time_period, "_Ensemble_", climate_scenario),
      File_Pred_sf = NA_character_,
      File_Pred_R = stringr::str_replace(
        File_Pred_summary, "_Summary.RData", "_R.qs")) %>%
    dplyr::distinct()

  # Overall summary -----
  Prediction_Summary <- Prediction_Summary %>%
    dplyr::rename(
      time_period = TimePeriod, climate_model = ClimateModel,
      climate_scenario = ClimateScenario) %>%
    dplyr::bind_rows(Prediction_Ensemble_Summary) %>%
    dplyr::mutate(hab_abb = Hab_Abb, hab_name = Hab_Name, .before = 1)

  save(
    Prediction_Summary,
    file = file.path(Path_Prediction, "Prediction_Summary.RData"))
  utils::write.table(
    x = Prediction_Summary, sep = "\t", row.names = FALSE, col.names = TRUE,
    file = file.path(Path_Prediction, "Prediction_Summary.txt"),
    quote = FALSE, fileEncoding = "UTF-8")

  # # ................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPrediction took ")

  return(Prediction_Summary)
}
