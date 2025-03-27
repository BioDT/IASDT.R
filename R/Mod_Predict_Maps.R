## |------------------------------------------------------------------------| #
# Predict_Maps ----
## |------------------------------------------------------------------------| #

#' Predict habitat suitability of `Hmsc` model across different climate options
#'
#' This function generates prediction maps of `Hmsc` models for current and
#' future climate options. It also predicts an ensemble predictions for
#' different climate models. For each species and overall species richness, the
#' function exports three maps: mean, standard deviation (sd), and coefficient
#' of variation (cov). For future predictions, the function generates maps for
#' prediction anomaly (future - current). The function prepares data to be
#' uploaded to the [OPeNDAP](http://opendap.biodt.eu/ias-pdt/) data server for
#' the use of the IAS-pDT [Shiny App](https://app.biodt.eu/).
#' @param Path_Model Character. Path to fitted `Hmsc` model object.
#' @param Hab_Abb Character. Habitat abbreviation indicating the specific
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type for
#'   which data will be prepared. Valid values are `0`, `1`, `2`, `3`, `4a`,
#'   `4b`, `10`, `12a`, `12b`. For more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param Pred_Clamp Logical indicating whether to clamp the sampling efforts at
#'   a single value. If `TRUE` (default), the `Fix_Efforts` argument must be
#'   provided.
#' @param Fix_Efforts Numeric or character. If `Pred_Clamp = TRUE`, the sampling
#'   efforts predictor with values U+02264 `Fix_Efforts` is fixed at
#'   `Fix_Efforts` during predictions. If numeric, the value is directly used
#'   (log<sub>10</sub> scale). If character, it can be one of `median`, `mean`,
#'   `max`, or `q90` (90% Quantile). Using `max` can reflect extreme values
#'   caused by rare, highly sampled locations (e.g., urban centers or popular
#'   natural reserves). While using 90% quantile avoid such extreme grid cells
#'   while still capturing areas with high sampling effort. This argument is
#'   mandatory when `Pred_Clamp` is set to `TRUE`.
#' @param Fix_Rivers Numeric or character. Similar to `Fix_Efforts`, but for
#'   fixing the length of rivers. If numeric, the value is directly used
#'   (log<sub>10</sub> scale). If character, it can be one of `median`, `mean`,
#'   `max`, `q90` (90% quantile). It can be also `NULL` for not fixing the river
#'   length predictor. Defaults to `q90`.
#' @param Pred_NewSites Logical. Whether to predict habitat suitability at new
#'   sites. Default: `TRUE`. Note: This parameter is temporary and will be
#'   removed in future updates.
#' @param LF_Only Logical. Whether to predict only the latent factor. This is
#'   useful for distributing processing load between GPU and CPU. When `LF_Only
#'   = TRUE`, latent factor prediction needs to be computed separately on GPU.
#'   When computations are finished on GPU, the function can later be rerun with
#'   `LF_Only = FALSE` (default) to predict habitat suitability using the
#'   already-computed latent factor predictions.
#' @param CC_Models Character vector. Climate models for future predictions.
#'   Available options are `c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
#'   "MRI-ESM2-0", "UKESM1-0-LL")` (default).
#' @param CC_Scenario Character vector. Climate scenarios for future
#'   predictions. Available options are: `c("ssp126", "ssp370", "ssp585")`
#'   (default).
#' @param Tar Logical. Whether to compress the add files into a single `*.tar`
#'   file (without compression). Default: `TRUE`.
#' @export
#' @name Predict_Maps
#' @author Ahmed El-Gabbas
#' @inheritParams Predict_Hmsc
#' @return A tibble containing the prediction summary and file paths for output
#'   `*.tif` files.
#' @export

Predict_Maps <- function(
    Path_Model = NULL, Hab_Abb = NULL, EnvFile = ".env", NCores = 8L,
    Pred_Clamp = TRUE, Fix_Efforts = "q90", Fix_Rivers = "q90",
    Pred_NewSites = TRUE, UseTF = TRUE, TF_Environ = NULL,
    TF_use_single = FALSE, LF_NCores = NCores, LF_Check = FALSE,
    LF_Temp_Cleanup = TRUE, LF_Only = FALSE, LF_Commands_Only = FALSE,
    Temp_Dir = "TEMP_Pred", Temp_Cleanup = TRUE, Tar = TRUE,
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

  # Check if `Hab_Abb` is a single character value
  if (length(Hab_Abb) != 1) {
    stop("`Hab_Abb` must be a single character value", call. = FALSE)
  }

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

  IASDT.R::CheckArgs(
    AllArgs = AllArgs,
    Args = c("Hab_Abb", "EnvFile", "Path_Model"), Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "RemoveChunks", Type = "logical")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("NCores", "LF_NCores", "ChunkSize"),
    Type = "numeric")

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###
  # # ..................................................................... ###

  ValidHabAbbs <- c(as.character(0:3), "4a", "4b", "10", "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      "Invalid Habitat abbreviation. Valid values are:\n >> ",
      toString(ValidHabAbbs), call. = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- Path_Roads <- Path_Rail <- Path_Bias <- tif_path <-
    time_period <- climate_model <- climate_scenario <- Path_CHELSA <-
    Path_CLC <- ias_id <- taxon_name <- ClimateModel <- ClimateScenario <-
    Name <- File_Pred_sf <- class <- order <- family <- species_name <-
    tif_path_mean <- tif_path_cov <- tif_path_sd <-
    Train <- Ensemble_File <- Ensemble_Maps <- tifs <- layer_name <-
    TimePeriod <- File_Pred_summary <- Ensemble_DT <- Dir_Ensemble <-
    File_Pred_R <- tif_path_anomaly <- Path_Rivers <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  Path_Eval <- IASDT.R::Path(dirname(dirname(Path_Model)), "Model_Evaluation")

  Path_Prediction1 <- IASDT.R::Path(
    dirname(dirname(Path_Model)), "Model_Prediction")
  Path_Prediction_Clamp <- IASDT.R::Path(Path_Prediction1, "Clamp")
  Path_Prediction_NoClamp <- IASDT.R::Path(Path_Prediction1, "NoClamp")
  fs::dir_create(c(Path_Eval, Path_Prediction_NoClamp))

  # Path for overall summary - paths for summaries of identical scenarios
  Path_Summary_RData <- IASDT.R::Path(
    dplyr::if_else(
      Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary.RData")
  Path_Summary_txt <- IASDT.R::Path(
    dplyr::if_else(
      Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary.txt")

  # Path for overall summary - for ShinyApp
  Path_Summary_RData_Shiny <- IASDT.R::Path(
    dplyr::if_else(
      Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary_Shiny.RData")
  Path_Summary_txt_Shiny <- IASDT.R::Path(
    dplyr::if_else(
      Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    "Prediction_Summary_Shiny.txt")

  # Check if the prediction summary is already available on disk
  if (all(file.exists(
    Path_Summary_RData, Path_Summary_txt,
    Path_Summary_RData_Shiny, Path_Summary_txt_Shiny))) {
    IASDT.R::CatTime(
      paste0(
        "All model predictions and prediction summary are already available ",
        "on disk"))
    return(invisible(NULL))
  }


  # Check Fix_Efforts value
  if (Pred_Clamp) {

    # Fix_Efforts can not be NULL when Clamping is implemented
    if (is.null(Fix_Efforts)) {
      stop(
        "`Fix_Efforts` can not be NULL when Clamping is implemented",
        call. = FALSE)
    }

    # Check if Fix_Efforts is a vector or length 1
    if (length(Fix_Efforts) != 1) {
      stop(
        "`Fix_Efforts` must be a vector or length 1.",
        " The current value is: ",
        paste(Fix_Efforts, collapse = " & "), call. = FALSE)
    }

    # Create folder for clamp results only if Pred_Clamp == TRUE
    fs::dir_create(Path_Prediction_Clamp)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Environment variables ------

  IASDT.R::CatTime("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Rail", "DP_R_Railways_processed", TRUE, FALSE,
    "Path_Roads", "DP_R_Roads_processed", TRUE, FALSE,
    "Path_CLC", "DP_R_CLC_processed", TRUE, FALSE,
    "Path_Bias", "DP_R_Efforts_processed", TRUE, FALSE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Rivers", "DP_R_Rivers_processed", TRUE, FALSE,
    "Path_CHELSA", "DP_R_CHELSA_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Loading input data -----
  IASDT.R::CatTime("Loading input data")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Species information -----

  SpeciesInfo <- IASDT.R::GetSpeciesName(EnvFile = EnvFile) %>%
    janitor::clean_names() %>%
    dplyr::select(ias_id, taxon_name, species_name, class, order, family)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Reference grid -----

  IASDT.R::CatTime("Reference grid", Level = 1)
  Path_GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    stop(
      "Path for the Europe boundaries does not exist: ", Path_GridR,
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Model object -----

  IASDT.R::CatTime("Model object", Level = 1)

  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop("Model path is NULL or does not exist: ", Path_Model, call. = FALSE)
  }

  Model <- IASDT.R::LoadAs(Path_Model)

  # `Pred_Clamp` can not be TRUE when `EffortsLog` is not used as predictor
  if (Pred_Clamp && isFALSE("EffortsLog" %in% names(Model$XData))) {
    stop(
      "`Pred_Clamp` can not be used when `EffortsLog` is not used as predictor",
      call. = FALSE)
  }

  OtherVars <- paste0(
    "^", c(IASDT.R::CHELSA_Vars$Variable, "CV"), collapse = "|") %>%
    stringr::str_subset(names(Model$XData), ., negate = TRUE)
  BioVars <- paste0("^", IASDT.R::CHELSA_Vars$Variable, collapse = "|") %>%
    stringr::str_subset(names(Model$XData), ., negate = FALSE)

  # Coordinates of the model sampling units
  Model_Coords <- as.data.frame(Model$rL$sample$s) %>%
    tibble::tibble() %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    dplyr::mutate(Train = TRUE)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## CHELSA data -----

  IASDT.R::CatTime("CHELSA data", Level = 1)

  Path_CHELSA <- IASDT.R::Path(Path_CHELSA, "CHELSA_Processed_DT.RData")
  if (!file.exists(Path_CHELSA)) {
    stop(
      "Processed CHLESA data can not be found at: ", Path_CHELSA,
      call. = FALSE)
  }

  Prediction_Options <- IASDT.R::LoadAs(Path_CHELSA) %>%
    dplyr::select(-"File_List") %>%
    dplyr::filter(
      # Filter only selected future climate models
      ClimateModel %in% c("Current", CC_Models),
      # Filter only selected future climate scenarios
      ClimateScenario %in% c("Current", CC_Scenario)) %>%
    dplyr::mutate(
      Name = paste0(TimePeriod, "_", ClimateScenario, "_", ClimateModel),
      Name = stringr::str_replace(Name, "1981-2010_Current_Current", "Current"),
      Name = stringr::str_replace_all(Name, "-", "_"))

  if (Pred_Clamp) {
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

  if ("RoadRailLog" %in% OtherVars) {

    IASDT.R::CatTime("Road and railway intensity", Level = 1)

    R_Railways <- IASDT.R::Path(Path_Rail, "Railways_Length.RData")
    if (!file.exists(R_Railways)) {
      stop("Railways data does not exist at: ", R_Railways, call. = FALSE)
    }
    R_Railways <- IASDT.R::LoadAs(R_Railways) %>%
      terra::unwrap() %>%
      magrittr::extract2("rail")

    R_Roads <- IASDT.R::Path(Path_Roads, "Road_Length.RData")
    if (!file.exists(R_Roads)) {
      stop("Roads data does not exist at: ", R_Roads, call. = FALSE)
    }
    R_Roads <- IASDT.R::LoadAs(R_Roads) %>%
      terra::unwrap() %>%
      magrittr::extract2("All")

    # Calculating the sum of road and railway intensity
    R_RoadRail <- (R_Roads + R_Railways) %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("RoadRailLog")

    StaticPredictors <- c(StaticPredictors, R_RoadRail)
    rm(R_RoadRail, R_Roads, R_Railways, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Habitat information ----

  # Check if habitat information is used as predictor
  Hab_Predictor <- "HabLog" %in% OtherVars

  if (Hab_Predictor) {

    IASDT.R::CatTime("Habitat information", Level = 1)
    R_Hab <- IASDT.R::Path(
      Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
    if (!file.exists(R_Hab)) {
      stop("Habitat data: '", R_Hab, "' does not exist", call. = FALSE)
    }

    R_Hab <- IASDT.R::LoadAs(R_Hab) %>%
      terra::unwrap() %>%
      magrittr::extract2(paste0("SynHab_", Hab_Abb))

    # Models are trained and predictions are made only at grid cells with > 0 %
    # coverage. Mask layer to exclude grid cells with zero % coverage from
    # predictions.
    R_Hab_Mask <- terra::classify(R_Hab, cbind(0, NA), others = 1)

    R_Hab <- log10(R_Hab + 0.1) %>%
      stats::setNames("HabLog")

    StaticPredictors <- c(StaticPredictors, R_Hab)
    rm(R_Hab, envir = environment())
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Sampling efforts -----

  if ("EffortsLog" %in% OtherVars) {

    IASDT.R::CatTime("Sampling efforts", Level = 1)

    R_Efforts <- IASDT.R::Path(Path_Bias, "Efforts_SummaryR.RData")
    if (!file.exists(R_Efforts)) {
      stop(
        "Sampling efforts data does not exist at: ", R_Efforts, call. = FALSE)
    }

    R_Efforts <- IASDT.R::LoadAs(R_Efforts) %>%
      terra::unwrap() %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("EffortsLog")


    if (Pred_Clamp) {

      IASDT.R::CatTime("Fixing sampling efforts values", Level = 2)

      # Check Fix_Efforts value
      if (is.numeric(Fix_Efforts)) {

        # If `Fix_Efforts` is numeric value, check if it is within the range of
        # the observed efforts
        EffortsRange <- terra::global(R_Efforts, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        InvalidVal <- isFALSE(
          dplyr::between(Fix_Efforts, EffortsRange[1], EffortsRange[2]))

        if (InvalidVal) {
          stop(
            "`Fix_Efforts` value (", Fix_Efforts, ") is out of the range of ",
            "the observed efforts: From ",
            paste(round(EffortsRange, 2), collapse = " to "), call. = FALSE)
        }

        # Fix value
        EffortsVal <- Fix_Efforts

      } else {

        # If `Fix_Efforts` is character, check if it is one of the valid values:
        # median, mean, max, and q90
        Fix_Efforts <- stringr::str_to_lower(Fix_Efforts)
        if (!(Fix_Efforts %in% c("median", "mean", "max", "q90"))) {
          stop(
            "`Fix_Efforts` has to be either NULL, single numeric ",
            "value, or one of the following: 'median', 'mean', 'max', ",
            "or `q90`. The current value is: ", Fix_Efforts, call. = FALSE)
        }
      }

      # Fix value
      EffortsVal <- dplyr::case_when(

        # Fix at 90% quantile
        Fix_Efforts == "q90" ~ {
          terra::global(
            R_Efforts,
            fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at median value
        Fix_Efforts == "median" ~ {
          terra::global(
            R_Efforts, fun = function(x) median(x, na.rm = TRUE)) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at mean value
        Fix_Efforts == "mean" ~ {
          terra::global(R_Efforts, fun = mean, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        # Fix at max value
        Fix_Efforts == "max" ~ {
          terra::global(R_Efforts, fun = max, na.rm = TRUE) %>%
            unlist() %>%
            as.numeric()
        },

        .default = NA_real_)


      # Fix at single value
      IASDT.R::CatTime(
        paste0("Fixed value is ", round(EffortsVal, 2), " [log10 scale]"),
        Level = 2, Time = FALSE)

      # Set a minimum value for efforts variable to `EffortsVal`. Using upper =
      # Inf keeps efforts values > EffortsVal as they are.
      R_Efforts_Clamp <- terra::clamp(
        x = R_Efforts, lower = EffortsVal, upper = Inf) %>%
        stats::setNames("EffortsLog_Clamp")

      StaticPredictors <- c(StaticPredictors, R_Efforts, R_Efforts_Clamp)
      rm(R_Efforts, R_Efforts_Clamp, envir = environment())

    } else {

      StaticPredictors <- c(StaticPredictors, R_Efforts)
      rm(R_Efforts, envir = environment())

    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## River length ----

  if ("RiversLog" %in% OtherVars) {

    IASDT.R::CatTime("River length", Level = 1)

    R_Rivers <- IASDT.R::Path(Path_Rivers, "River_Lengths.RData")
    if (!file.exists(R_Rivers)) {
      stop("River length data does not exist at: ", R_Rivers, call. = FALSE)
    }

    R_Rivers <- IASDT.R::LoadAs(R_Rivers) %>%
      terra::unwrap() %>%
      magrittr::extract2("STRAHLER_5") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("RiversLog")


    if (is.null(Fix_Rivers)) {

      # Do not fix at single value
      IASDT.R::CatTime(
        "River length predictor is not fixed at a single value",
        Level = 2, Time = FALSE)

    } else {

      IASDT.R::CatTime(
        "River length predictor will be fixed at single (`Fix_Rivers`) value",
        Level = 2, Time = FALSE)


      # Check if the Fix_Rivers value is valid

      if (length(Fix_Rivers) != 1) {
        # Check if Fix_Rivers is a vector or length 1
        stop(
          "`Fix_Rivers` must be a vector or length 1. The current value is: ",
          paste(Fix_Rivers, collapse = " & "), call. = FALSE)
      }

      if (is.numeric(Fix_Rivers)) {

        # If `Fix_Rivers` is numeric value, check if it is within the range of
        # the observed river lengths
        RiversRange <- terra::global(R_Rivers, fun = range, na.rm = TRUE) %>%
          unlist() %>%
          as.vector()

        InvalidVal <- isFALSE(
          dplyr::between(Fix_Rivers, RiversRange[1], RiversRange[2]))

        if (InvalidVal) {
          stop(
            "`Fix_Rivers` value (", Fix_Rivers, ") is out of the range of ",
            "the observed river length: From ",
            paste(round(RiversRange, 2), collapse = " to "), call. = FALSE)
        }

        # Fix value
        RiversVal <- Fix_Rivers

      } else {

        # If `Fix_Rivers` is character, check if it is one of the valid values:
        # median, mean, max, and q90
        Fix_Rivers <- stringr::str_to_lower(Fix_Rivers)
        if (!(Fix_Rivers %in% c("median", "mean", "max", "q90"))) {
          stop(
            "`Fix_Rivers` has to be either NULL, single numeric ",
            "value, or one of the following: 'median', 'mean', 'max', ",
            "or 'q90'. The current value is: ", Fix_Rivers, call. = FALSE)
        }

        # Fix value
        RiversVal <- dplyr::case_when(

          # Fix at 90% quantile
          Fix_Rivers == "q90" ~ {
            terra::global(
              R_Rivers,
              fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at median value
          Fix_Rivers == "median" ~ {
            terra::global(
              R_Rivers, fun = function(x) median(x, na.rm = TRUE)) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at mean value
          Fix_Rivers == "mean" ~ {
            terra::global(R_Rivers, fun = mean, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          # Fix at max value
          Fix_Rivers == "max" ~ {
            terra::global(R_Rivers, fun = max, na.rm = TRUE) %>%
              unlist() %>%
              as.numeric()
          },

          .default = NA_real_)
      }

      IASDT.R::CatTime(
        paste0("Fixed value is ", round(RiversVal, 2), " [log10 scale]"),
        Level = 2, Time = FALSE)

      # Set a minimum value for `RiversLog` variable to `RiversVal`. Using upper
      # = Inf keeps RiversLog values > RiversVal as they are.
      R_Rivers <- terra::clamp(
        x = R_Rivers, lower = RiversVal, upper = Inf) %>%
        stats::setNames("RiversLog")

      rm(RiversVal, envir = environment())

    }

    StaticPredictors <- c(StaticPredictors, R_Rivers)
    rm(R_Rivers, envir = environment())

  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  IASDT.R::CatTime("Merge static predictors", Level = 1)

  StaticPredictors <- terra::rast(StaticPredictors)

  # If Habitat predictor is used, grid cells with zero % coverage are
  # excluded from predictions
  if (Hab_Predictor) {
    StaticPredictors <- terra::mask(StaticPredictors, R_Hab_Mask)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predict latent factor at new locations ------

  IASDT.R::CatTime("Predict latent factor at new locations")

  Path_Test_LF <- IASDT.R::Path(Path_Prediction1, "Test_LF.qs2")

  if (!file.exists(Path_Test_LF) && Pred_NewSites) {

    IASDT.R::CatTime(
      "Preparing input data for predicting latent factor", Level = 1)

    Predict_DF_Test <- Prediction_Options %>%
      dplyr::filter(ClimateModel == "Current") %>%
      dplyr::pull("FilePath") %>%
      # If Pred_Clamp`=`TRUE`, there are two options for Current climate data
      # (with and without clamping). Two sets of predictions under current
      # climates will be produced. Predictions without clamping is used for
      # model evaluation.
      utils::head(1) %>%
      IASDT.R::LoadAs() %>%
      terra::unwrap() %>%
      # Only extract predictors used in the model
      terra::subset(BioVars) %>%
      # Combine with other static predictors
      c(StaticPredictors) %>%
      # If Habitat predictor is used, grid cells with zero % coverage are
      # excluded from predictions [na.rm = TRUE]
      terra::as.data.frame(xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
      tibble::tibble() %>%
      sf::st_as_sf(remove = FALSE, coords = c("x", "y"), crs = 3035) %>%
      sf::st_join(Model_Coords) %>%
      tidyr::replace_na(list(Train = FALSE)) %>%
      dplyr::filter(Train == FALSE)

    Test_XY <- sf::st_drop_geometry(Predict_DF_Test[, c("x", "y")])
    Test_X <- Predict_DF_Test %>%
      dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
      sf::st_drop_geometry()
    Gradient <- Hmsc::prepareGradient(
      hM = Model, XDataNew = as.data.frame(Test_X),
      sDataNew = list(sample = as.data.frame(Test_XY)))

    rm(Predict_DF_Test, Test_X, Test_XY, Model, envir = environment())
    invisible(gc())

    IASDT.R::CatTime("Predicting latent factor", Level = 1)
    IASDT.R::CatSep(Extra1 = 1, Extra2 = 2, Rep = 1, Char = "*")

    # Predicting latent factor only -- no predictions are made
    Preds_LF <- IASDT.R::Predict_Hmsc(
      Path_Model = Path_Model, Gradient = Gradient, expected = TRUE,
      NCores = NCores, Model_Name = paste0("LF_", Hab_Abb, "_Test"),
      Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup, UseTF = UseTF,
      TF_Environ = TF_Environ, LF_OutFile = Path_Test_LF,
      TF_use_single = TF_use_single, LF_Only = TRUE, LF_NCores = LF_NCores,
      LF_Check = LF_Check, LF_Temp_Cleanup = LF_Temp_Cleanup,
      LF_Commands_Only = LF_Commands_Only, Evaluate = FALSE, Verbose = TRUE)

    rm(Gradient, Preds_LF, envir = environment())

    IASDT.R::CatTime("Predicting latent factor is finished!", Level = 1)
    IASDT.R::CatSep(Extra1 = 1, Extra2 = 2, Rep = 1, Char = "*")

    if (LF_Commands_Only) {
      return(invisible(NULL))
    }

  } else {
    if (Pred_NewSites) {
      IASDT.R::CatTime("LF prediction is already available on disk", Level = 1)
    } else {
      IASDT.R::CatTime("LF prediction will NOT be made", Level = 1)
    }

    rm(Model, envir = environment())
    invisible(gc())
  }


  if (LF_Only) {
    return(invisible())
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
    Model_Name <- paste0(
      Option_Name, " - ", dplyr::if_else(DoClamp, "clamping", "no clamping"))

    MSG <- paste0(
      Model_Name, " (", ID, "/", nrow(Prediction_Options), ")")
    cat("\n")
    IASDT.R::InfoChunk(
      paste0("\t", MSG), Rep = 1, Char = "-", CharReps = 70, Red = TRUE,
      Bold = TRUE, Time = FALSE)

    if (DoClamp) {

      # Do not evaluate for options with clamping
      Evaluate <- FALSE

      # Make prediction files at `Path_Prediction_Clamp`
      Path_Prediction <- Path_Prediction_Clamp

      # use clamped Effort values
      StaticPreds <- terra::subset(
        x = StaticPredictors, subset = "EffortsLog", negate = TRUE)
      StaticPreds$EffortsLog <- StaticPreds$EffortsLog_Clamp
      StaticPreds$EffortsLog_Clamp <- NULL

    } else {

      # Evaluate for "Current" climates, without clamping
      Evaluate <- (Option_Name == "Current")

      # Make prediction files at `Path_Prediction_NoClamp`
      Path_Prediction <- Path_Prediction_NoClamp

      # use original effort data
      if ("EffortsLog_Clamp" %in% names(StaticPredictors)) {
        StaticPreds <- terra::subset(
          x = StaticPredictors, subset = "EffortsLog_Clamp", negate = TRUE)
      } else {
        StaticPreds <- StaticPredictors
      }
    }

    Path_Prediction_sf <- IASDT.R::Path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_sf.qs2"))
    Path_Prediction_R <- IASDT.R::Path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_R.qs2"))
    Path_Prediction_summary <- IASDT.R::Path(
      Path_Prediction, paste0("Prediction_", Option_Name, "_Summary.RData"))

    # Path for saving tif files of the current option
    Path_Prediction_tif <- IASDT.R::Path(Path_Prediction, Option_Name)
    fs::dir_create(Path_Prediction_tif)
    invisible(gc())

    # ______________________________________________
    # ______________________________________________

    # Making predictions if not already processed ----
    OutMissing <- file.exists(
      c(Path_Prediction_R, Path_Prediction_sf, Path_Prediction_summary)) %>%
      all() %>%
      isFALSE()

    if (OutMissing) {

      # Load model
      Model <- IASDT.R::LoadAs(Path_Model)

      # Skip predictions if the predictions as sf object is already on disk
      if (file.exists(Path_Prediction_sf)) {

        IASDT.R::CatTime("Loading predictions `sf` from disk", Level = 1)
        Prediction_sf <- IASDT.R::LoadAs(Path_Prediction_sf)

      } else {

        .OptionStartTime <- lubridate::now(tzone = "CET")

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

        # Training locations
        Model_Name_Train <- paste0(
          Option_Name, "_",
          dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Train")
        Predict_DF_Train <- dplyr::filter(Predict_DF, Train)
        Train_XY <- sf::st_drop_geometry(Predict_DF_Train[, c("x", "y")])
        Train_PA <- as.data.frame(Model$Y)
        Train_X <- Predict_DF_Train %>%
          dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
          sf::st_drop_geometry() %>%
          stats::model.matrix(Model$XFormula, ., xlev = NULL)

        # Testing locations
        Model_Name_Test <- paste0(
          Option_Name, "_",
          dplyr::if_else(DoClamp, "Clamping", "NoClamping"), "_Test")
        Predict_DF_Test <- dplyr::filter(Predict_DF, isFALSE(Train))
        Test_XY <- Predict_DF_Test[, c("x", "y")]
        Test_X <- Predict_DF_Test %>%
          dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
          sf::st_drop_geometry()
        Gradient <- Hmsc::prepareGradient(
          hM = Model, XDataNew = as.data.frame(Test_X),
          sDataNew = list(
            sample = as.data.frame(sf::st_drop_geometry(Test_XY))))

        rm(
          Model, Predict_DF_Test, Predict_DF_Train, Predict_DF,
          envir = environment())
        invisible(gc())

        # ______________________________________________

        # Predictions at training sites ----
        IASDT.R::CatTime("Predictions at training sites", Level = 1)

        Path_Current_Train <- IASDT.R::Path(
          Path_Prediction, paste0("Prediction_", Option_Name, "_Train.qs2"))

        if (file.exists(Path_Current_Train)) {

          IASDT.R::CatTime("Loading predictions from disk", Level = 2)
          Preds_ModFitSites <- tibble::tibble(Pred_Path = Path_Current_Train)

        } else {

          Preds_ModFitSites <- IASDT.R::Predict_Hmsc(
            Path_Model = Path_Model, X = Train_X, Gradient = NULL,
            expected = TRUE, NCores = NCores, Model_Name = Model_Name_Train,
            Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup, UseTF = UseTF,
            TF_Environ = TF_Environ, TF_use_single = TF_use_single,
            LF_Return = TRUE, LF_NCores = LF_NCores, LF_Check = LF_Check,
            LF_Temp_Cleanup = LF_Temp_Cleanup, LF_Commands_Only = FALSE,
            Pred_Dir = Path_Prediction, Pred_PA = Train_PA, Pred_XY = Train_XY,
            Evaluate = Evaluate, Eval_Name = NULL, Eval_Dir = Path_Eval,
            Verbose = FALSE)

        }

        # ______________________________________________

        # Predictions at new sites ----

        if (Pred_NewSites) {

          IASDT.R::CatTime("Predictions at new sites", Level = 1)

          Path_Current_Test <- IASDT.R::Path(
            Path_Prediction, paste0("Prediction_", Option_Name, "_Test.qs2"))

          if (file.exists(Path_Current_Test)) {

            IASDT.R::CatTime("Loading predictions from disk", Level = 2)
            Preds_NewSites <- tibble::tibble(Pred_Path = Path_Current_Test)

          } else {

            Preds_NewSites <- IASDT.R::Predict_Hmsc(
              Path_Model = Path_Model, Gradient = Gradient, expected = TRUE,
              NCores = NCores, Model_Name = Model_Name_Test,
              Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup, UseTF = UseTF,
              TF_Environ = TF_Environ, TF_use_single = TF_use_single,
              LF_Return = TRUE, LF_InputFile = Path_Test_LF,
              LF_NCores = LF_NCores, LF_Check = LF_Check,
              LF_Temp_Cleanup = LF_Temp_Cleanup, LF_Commands_Only = FALSE,
              Verbose = FALSE, Pred_Dir = Path_Prediction, Evaluate = FALSE,
              Pred_XY = sf::st_drop_geometry(Test_XY))

          }

          # ______________________________________________

          # Merge & save predictions - sf ------
          IASDT.R::CatTime(
            "Merge & save predictions at training and new sites", Level = 1)

          Prediction_sf <- dplyr::bind_rows(
            IASDT.R::LoadAs(Preds_ModFitSites$Pred_Path),
            IASDT.R::LoadAs(Preds_NewSites$Pred_Path))

        } else {

          IASDT.R::CatTime(
            "Predictions at new sites will NOT be made", Level = 1)

          # Get column names from predictions at training sites
          ColsToAdd <- IASDT.R::LoadAs(Preds_ModFitSites$Pred_Path) %>%
            names() %>%
            stringr::str_subset("^Sp_|^SR_|Model_Name")

          # Add missing columns to predictions at new sites
          Preds_New_NA <- Test_XY %>%
            dplyr::mutate(Model_Name = Model_Name_Test) %>%
            IASDT.R::AddMissingCols(FillVal = NA_real_, ColsToAdd)

          # combine predictions
          Prediction_sf <- IASDT.R::LoadAs(Preds_ModFitSites$Pred_Path) %>%
            dplyr::bind_rows(Preds_New_NA)

          rm(Preds_New_NA, ColsToAdd, envir = environment())
        }

        # ______________________________________________

        # Save predictions as sf object
        IASDT.R::SaveAs(InObj = Prediction_sf, OutPath = Path_Prediction_sf)

        IASDT.R::CatDiff(
          InitTime = .OptionStartTime, Prefix = "Prediction took ", Level = 1)

      }

      # ______________________________________________
      # ______________________________________________

      ### Predictions as spatRaster / tif -----

      IASDT.R::CatTime("Rasterization & prepare summary data", Level = 1)

      Fields2Raster <- names(Prediction_sf) %>%
        stringr::str_subset("^Sp_|^SR_") %>%
        gtools::mixedsort()

      Grid10 <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))
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
          IASDT.R::LoadAs() %>%
          terra::unwrap() %>%
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
        hab_abb = Hab_Abb, hab_name = Hab_Name,
        layer_name = Fields2Raster,
        time_period = Prediction_Options$TimePeriod[[ID]],
        climate_model = Prediction_Options$ClimateModel[[ID]],
        climate_scenario = Prediction_Options$ClimateScenario[[ID]],
        Clamp = Prediction_Options$Clamp[[ID]],
        Path_Prediction = Path_Prediction) %>%
        dplyr::mutate(
          Stats = dplyr::case_when(
            stringr::str_detect(layer_name, "_mean$") ~ "tif_path_mean",
            stringr::str_detect(layer_name, "_sd$") ~ "tif_path_sd",
            stringr::str_detect(layer_name, "_cov$") ~ "tif_path_cov",
            stringr::str_detect(layer_name, "_anomaly$") ~ "tif_path_anomaly",
            .default = NULL),
          ias_id = stringr::str_remove(
            layer_name, "_mean$|_sd$|_cov$|_anomaly"),
          tif_path = IASDT.R::Path(
            Path_Prediction_tif, paste0(layer_name, ".tif")))

      # Save as tif
      IASDT.R::CatTime("Save as tif", Level = 1)
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
      IASDT.R::CatTime("Save as spatRaster - qs2", Level = 1)
      Prediction_R <- terra::wrap(Prediction_R)
      IASDT.R::SaveAs(InObj = Prediction_R, OutPath = Path_Prediction_R)

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
    return(
      tibble::tibble(
        Name = Option_Name,
        File_Pred_R = Path_Prediction_R, File_Pred_sf = Path_Prediction_sf,
        File_Pred_summary = Path_Prediction_summary))
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predicting ------

  IASDT.R::InfoChunk(
    paste0("\t", "Making spatial predictions"), Rep = 2, Char = "*",
    CharReps = 70, Red = TRUE, Bold = TRUE, Time = FALSE)

  Grid10 <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))

  Prediction_Summary <- purrr::map_dfr(
    .x = seq_len(nrow(Prediction_Options)), .f = Predict_Internal) %>%
    dplyr::full_join(Prediction_Options, ., by = c("Name", "Clamp")) %>%
    dplyr::select(-"FilePath")

  rm(Predict_Internal, Grid10, Model_Coords, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Ensemble model predictions ------

  IASDT.R::InfoChunk(
    "\tEnsemble model predictions", Rep = 1, Char = "-", CharReps = 70,
    Red = TRUE, Bold = TRUE, Time = FALSE)

  IASDT.R::CatTime("Prepare working on parallel", Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(Prediction_Summary)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # --------------------------------------------------------- #

  # Prepare input data to calculate ensemble predictions
  IASDT.R::CatTime(
    "Prepare input data to calculate ensemble predictions", Level = 1)

  Prediction_Ensemble <- Prediction_Summary %>%
    dplyr::filter(ClimateModel != "Current") %>%
    dplyr::select(-File_Pred_sf, -File_Pred_R, -Name, -ClimateModel) %>%
    dplyr::mutate(
      Prediction2 = furrr::future_map(
        .x = File_Pred_summary,
        .f = ~ {
          if (!file.exists(.x)) {
            warning("File not found: ", .x, call. = FALSE)
            return(NULL)
          }

          IASDT.R::LoadAs(.x) %>%
            dplyr::select(
              -tidyselect::all_of(
                c("tif_path_sd", "tif_path_cov", "tif_path_anomaly")))

        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = 1,
          packages = c("tidyselect", "dplyr", "IASDT.R")))) %>%
    dplyr::select("Prediction2") %>%
    tidyr::unnest("Prediction2") %>%
    dplyr::mutate(
      climate_model = "Ensemble",
      Dir_Ensemble = IASDT.R::Path(
        dirname(dirname(tif_path_mean)),
        paste0(
          stringr::str_replace(time_period, "-", "_"), "_", climate_scenario,
          "_Ensemble"))) %>%
    dplyr::group_by(dplyr::across(-tif_path_mean)) %>%
    dplyr::summarise(tifs = list(tif_path_mean), .groups = "drop") %>%
    dplyr::mutate(
      tif_path_mean = IASDT.R::Path(Dir_Ensemble, paste0(ias_id, "_mean.tif")),
      tif_path_anomaly = IASDT.R::Path(
        Dir_Ensemble, paste0(ias_id, "_anomaly.tif")),
      tif_path_sd = IASDT.R::Path(Dir_Ensemble, paste0(ias_id, "_sd.tif")),
      tif_path_cov = IASDT.R::Path(Dir_Ensemble, paste0(ias_id, "_cov.tif")))

  # --------------------------------------------------------- #

  IASDT.R::CatTime("Create directories for ensemble predictions", Level = 1)
  fs::dir_create(unique(Prediction_Ensemble$Dir_Ensemble))

  # --------------------------------------------------------- #

  # Loading mean predictions at current climates
  IASDT.R::CatTime("Loading mean predictions at current climates", Level = 1)

  CurrentMean <- list.files(
    path = dplyr::if_else(
      Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
    pattern = "Prediction_Current.*_R.qs2",
    full.names = TRUE) %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap()
  CurrentMean <- terra::wrap(CurrentMean["_mean"])

  # --------------------------------------------------------- #

  # Calculate ensemble predictions
  IASDT.R::CatTime("Calculate ensemble predictions", Level = 1)

  Prediction_Ensemble <- Prediction_Ensemble %>%
    dplyr::mutate(
      Ensemble_Maps = furrr::future_pmap(
        .l = list(
          tifs, tif_path_mean, tif_path_anomaly,
          tif_path_sd, tif_path_cov, ias_id),
        .f = function(tifs, tif_path_mean, tif_path_anomaly,
                      tif_path_sd, tif_path_cov, ias_id) {

          # load maps for future climate option
          tiffs_R <- terra::rast(tifs)

          # Mean maps
          Ensemble_mean <- terra::app(tiffs_R, "mean", na.rm = TRUE) %>%
            stats::setNames(paste0(ias_id, "_mean"))
          terra::writeRaster(
            x = Ensemble_mean, filename = tif_path_mean,
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

          # Anomaly maps: calculate prediction anomaly for ensemble future
          # projection
          CurrentMean0 <- terra::unwrap(CurrentMean) %>%
            terra::subset(paste0(ias_id, "_mean"))
          Ensemble_anomaly <- (Ensemble_mean - CurrentMean0) %>%
            stats::setNames(paste0(ias_id, "_anomaly"))
          terra::writeRaster(
            x = Ensemble_anomaly, filename = tif_path_anomaly,
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

          rm(CurrentMean0, envir = environment())

          # Standard deviation
          Ensemble_sd <- terra::app(tiffs_R, "sd", na.rm = TRUE) %>%
            stats::setNames(paste0(ias_id, "_sd"))
          terra::writeRaster(
            x = Ensemble_sd, filename = tif_path_sd,
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

          # coefficient of variation
          Ensemble_mean0 <- Ensemble_mean
          # Replace very small mean values with reasonably small number to avoid
          # overflow warning
          Ensemble_mean0[Ensemble_mean0 < 1e-8] <- 1e-8
          Ensemble_cov <- (Ensemble_sd / Ensemble_mean0) %>%
            stats::setNames(paste0(ias_id, "_cov"))
          terra::writeRaster(
            x = Ensemble_cov, filename = tif_path_cov,
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

          Out <- c(
            Ensemble_mean, Ensemble_sd, Ensemble_cov, Ensemble_anomaly) %>%
            terra::wrap()

          return(Out)
        },
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = 1, packages = "terra",
          globals = "CurrentMean")),
      Dir_Ensemble = NULL, tifs = NULL)

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)

  # --------------------------------------------------------- #

  # Save ensemble maps as SpatRast
  IASDT.R::CatTime("Save ensemble predictions as SpatRast", Level = 1)

  Prediction_Ensemble_R <- Prediction_Ensemble %>%
    dplyr::mutate(
      Ensemble_File = IASDT.R::Path(
        dplyr::if_else(
          Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
        paste0(
          "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
          climate_scenario, "_Ensemble_R.qs2"))) %>%
    dplyr::select(Ensemble_File, Ensemble_Maps) %>%
    dplyr::group_by(Ensemble_File) %>%
    tidyr::nest(Ensemble_Maps = Ensemble_Maps) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Ensemble_Maps = purrr::map(
        .x = Ensemble_Maps,
        .f = ~ {
          purrr::map(.x$Ensemble_Maps, terra::unwrap) %>%
            purrr::map(c) %>%
            terra::rast() %>%
            terra::wrap()
        }),
      Save = purrr::map2(
        .x = Ensemble_Maps, .y = Ensemble_File,
        .f = function(x, y) {
          IASDT.R::SaveAs(InObj = x, OutPath = y)
        }))

  rm(Prediction_Ensemble_R, envir = environment())

  # --------------------------------------------------------- #

  # Save summary of ensemble predictions
  IASDT.R::CatTime("Save summary of ensemble predictions", Level = 1)

  Prediction_Ensemble_Summary <- Prediction_Ensemble %>%
    dplyr::select(-Ensemble_Maps) %>%
    dplyr::mutate(
      Ensemble_File = IASDT.R::Path(
        dplyr::if_else(
          Pred_Clamp, Path_Prediction_Clamp, Path_Prediction_NoClamp),
        paste0(
          "Prediction_", stringr::str_replace(time_period, "-", "_"), "_",
          climate_scenario, "_Ensemble_Summary.RData"))) %>%
    tidyr::nest(Ensemble_DT = -Ensemble_File) %>%
    dplyr::mutate(
      Ensemble_Save = purrr::map2(
        .x = Ensemble_DT, .y = Ensemble_File,
        .f = ~ {

          IASDT.R::SaveAs(
            InObj = .x, OutPath = .y,
            OutObj = stringr::str_remove(basename(.y), ".RData"))

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

  IASDT.R::InfoChunk(
    "\tPrepare overall summary", Rep = 1, Char = "-", CharReps = 70,
    Red = TRUE, Bold = TRUE, Time = FALSE)

  Prediction_Summary <- Prediction_Summary %>%
    dplyr::rename(
      time_period = TimePeriod, climate_model = ClimateModel,
      climate_scenario = ClimateScenario) %>%
    dplyr::bind_rows(Prediction_Ensemble_Summary) %>%
    dplyr::mutate(hab_abb = Hab_Abb, hab_name = Hab_Name, .before = 1)

  save(Prediction_Summary, file = Path_Summary_RData)
  utils::write.table(
    x = Prediction_Summary, sep = "\t", row.names = FALSE, col.names = TRUE,
    file = Path_Summary_txt, quote = FALSE, fileEncoding = "UTF-8")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Overall summary - to be uploaded to the data server; for the Shiny App -----

  Prediction_Summary_Shiny <- Prediction_Summary %>%
    dplyr::pull(File_Pred_summary) %>%
    basename() %>%
    IASDT.R::Path(dirname(Path_Summary_RData_Shiny), .) %>% 
    purrr::map(IASDT.R::LoadAs) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct()

  save(Prediction_Summary_Shiny, file = Path_Summary_RData_Shiny)
  utils::write.table(
    x = Prediction_Summary_Shiny, sep = "\t", row.names = FALSE,
    col.names = TRUE, file = Path_Summary_txt_Shiny, quote = FALSE,
    fileEncoding = "UTF-8")

  if (Tar) {

    IASDT.R::CatTime("Create tar file for prediction files", Level = 1)

    # Directory to save the tar file
    TarDir <- dirname(Path_Summary_RData_Shiny)
    # Path to the tar file
    TarFile <- IASDT.R::Path(TarDir, "Predictions.tar")
    # List of directories in the prediction folder. All directories will be
    # included in the tar file
    TarFiles <- list.dirs(
      path = TarDir, full.names = FALSE, recursive = FALSE) %>%
      paste(collapse = " ") %>%
      # Add the summary files to the list
      paste(
        "Prediction_Summary.RData", "Prediction_Summary_Shiny.txt",
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

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nThe whole prediction function took ")

  return(Prediction_Summary)
}
