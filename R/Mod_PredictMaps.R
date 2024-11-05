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
    NCores = 8, ChunkSize = 500, RemoveChunks = TRUE, TF_Environ = NULL,
    UseTF = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")
  Hab_Abb <- as.character(Hab_Abb)

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
  Path_Grid <- Path_Roads <- Path_Rail <- Path_Bias <- Name <- FilePath <-
    time_period <- climate_model <- climate_scenario <- Path_CHELSA <-
    Path_CLC <- Tif_Path <- path_prediction <- ias_id <- taxon_name <-
    class <- order <- family <- species_name <- tif_path_mean <-
    tif_path_cov <- tif_path_sd <- hab_abb <- hab_name <- TimePeriod <-
    ClimateModel <- ClimateScenario <- Time <- Scenario <- Train <- NULL

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
  Path_Predictions <- file.path(dirname(dirname(Path_Model)), "Model_Predictions")
  Path_Eval <- file.path(dirname(dirname(Path_Model)), "Evaluation")
  fs::dir_create(c(Path_Predictions, Path_Eval))

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

    R_Bias <- file.path(Path_Bias, "Efforts_SummaryR.RData")
    if (!file.exists(R_Bias)) {
      stop(
        paste0("Sampling efforts data does not exist at: ", R_Bias),
        call. = FALSE)
    }

    R_Bias <- IASDT.R::LoadAs(R_Bias) %>%
      terra::unwrap() %>%
      magrittr::extract2("NObs") %>%
      magrittr::add(0.1) %>%
      log10() %>%
      stats::setNames("EffortsLog")

    StaticPredictors <- c(StaticPredictors, R_Bias)
    rm(R_Bias)
  }

  # # |||||||||||||||||||||||||||||||||||||||| #

  ## Merge static predictors -----
  IASDT.R::CatTime("Merge static predictors", Level = 1)

  StaticPredictors <- terra::rast(StaticPredictors)

  # # ..................................................................... ###

  # Predictions - current ------

  Name <- "Current"

  Predict_DF <- CHELSA_Data %>%
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
    tidyr::replace_na(list(Train = FALSE))

  ## training locations -----
  Model_Name_Train <- paste0(Name, "_", Hab_Abb, "_Train")
  Predict_DF_Train <- dplyr::filter(Predict_DF, Train)
  Train_XY <- sf::st_drop_geometry(Predict_DF_Train[, c("x", "y")])
  Train_PA <- as.data.frame(Model$Y)
  Train_X <- Predict_DF_Train %>%
    dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
    sf::st_drop_geometry() %>%
    stats::model.matrix(Model$XFormula, ., xlev = NULL)

  # Testing Locations ----
  Model_Name_Test <- paste0(Name, "_", Hab_Abb, "_Test")
  Predict_DF_Test <- dplyr::filter(Predict_DF, !Train)
  Path_Test_LF <- file.path(Path_Predictions, "Test_LF.qs")
  Test_XY <- sf::st_drop_geometry(Predict_DF_Test[, c("x", "y")])
  Test_X <- Predict_DF_Test %>%
    dplyr::select(tidyselect::all_of(names(Model$XData))) %>%
    sf::st_drop_geometry()
  Gradient <- Hmsc::prepareGradient(
    hM = Model, XDataNew = as.data.frame(Test_X),
    sDataNew = list(sample = as.data.frame(Test_XY)))
  Test_X <- Test_X %>%
    stats::model.matrix(Model$XFormula, ., xlev = NULL)

  rm(Model, Predict_DF_Test, Predict_DF)
  invisible(gc())


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
    Pred_Dir = Path_Predictions,
    Pred_PA = Train_PA,
    Pred_XY = Train_XY,
    Evaluate = TRUE,
    Eval_Name = "Training",
    Eval_Dir = Path_Eval,
    # Verbose = FALSE
  ) %>%
    dplyr::mutate(Scenario = Model_Name_Train)

  Preds_Current_Test <- IASDT.R::predictHmsc(
    object = Path_Model,
    # X = Test_X,
    Gradient = Gradient,
    expected = TRUE,
    NCores = NCores,
    Model_Name = Model_Name_Test,
    Temp_Dir = "TEMP2Pred",
    UseTF = UseTF,
    TF_Environ = TF_Environ,
    LF_Return = TRUE,
    LF_OutFile = Path_Test_LF,
    Pred_Dir = Path_Predictions,
    # Pred_PA = Train_PA,
    Pred_XY = Test_XY,
    Evaluate = FALSE,
    # Eval_Name = "Training",
    # Eval_Dir = Path_Eval,
    # Verbose = FALSE
  ) %>%
    dplyr::mutate(Scenario = Model_Name_Test)



# # ..................................................................... ###

# Predicting ------

IASDT.R::CatTime("\nMaking predictions", ... = "\n")

Predictions <- CHELSA_Data %>%
  dplyr::mutate(
    path_prediction = purrr::pmap_chr(
      .l = list(Name, FilePath, TimePeriod, ClimateModel, ClimateScenario),
      .f = function(Name, FilePath, TimePeriod,
                    ClimateModel, ClimateScenario) {

        Try <- 0
        while (TRUE) {
          Try <- Try + 1
          Path <- Predict_Scenario(
            Name = Name, File = FilePath,
            Hab_Abb = Hab_Abb, ChunkSize = ChunkSize, NCores = NCores,
            StaticPredictors = StaticPredictors, Path_Model = Path_Model,
            Path_Grid = Path_GridR, EnvFile = EnvFile, FromHPC = FromHPC,
            RemoveChunks = RemoveChunks, TimePeriod = TimePeriod,
            ClimateModel = ClimateModel, ClimateScenario = ClimateScenario)

          if (file.exists(Path) && IASDT.R::CheckRData(Path)) {
            break
          }

          if (Try > 5) {
            Path <- ""
            break
          }
        }
        return(Path)
      }
    )
  ) %>%
  dplyr::select(-c("File_List", "Name", "FilePath")) %>%
  dplyr::rename(
    time_period = TimePeriod, climate_model = ClimateModel,
    climate_scenario = ClimateScenario)


withr::local_options(
  future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

if (NCores == 1) {
  future::plan("future::sequential", gc = TRUE)
} else {
  c1 <- snow::makeSOCKcluster(min(NCores, nrow(Predictions)))
  on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  future::plan("future::cluster", workers = c1, gc = TRUE)
  on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
}

Predictions <- dplyr::mutate(
  Predictions,
  Predictions2 = furrr::future_map(
    .x  = path_prediction,
    .f = ~{
      if (file.exists(.x)) {
        IASDT.R::LoadAs(.x) %>%
          dplyr::select(
            -tidyselect::all_of(c(
              "time_period", "climate_model", "climate_scenario",
              "predictions_mean", "predictions_sd", "predictions_cov")))
      } else {
        warning(paste0("File not found: ", .x))
        return(NULL)
      }
    },
    .options = furrr::furrr_options(
      seed = TRUE, scheduling = Inf,
      packages = c("tidyselect", "dplyr", "IASDT.R", "terra")
    ))) %>%
  tidyr::unnest("Predictions2")

# # ................................................................... ###

# Ensemble of climate models ------

IASDT.R::CatTime("\nEnsemble of climate models", ... = "\n")

# Create ensemble paths
expand.grid(
  Time = c("2011_2040", "2041_2070", "2071_2100"),
  Scenario = c("ssp126", "ssp370", "ssp585")) %>%
  dplyr::mutate(Folder = paste0("Ensemble_", Time, "_", Scenario)) %>%
  dplyr::pull("Folder") %>%
  file.path(Path_Predictions, .) %>%
  fs::dir_create()

Ensemble <- dplyr::filter(Predictions, time_period != "1981-2010") %>%
  dplyr::select(
    tidyselect::all_of(c(
      "time_period", "climate_scenario", "ias_id",
      "tif_path_mean", "hab_abb", "hab_name"))) %>%
  dplyr::mutate(climate_model = "Ensemble", .before = "ias_id") %>%
  dplyr::group_by(
    time_period, climate_scenario, climate_model,
    ias_id, hab_abb, hab_name) %>%
  dplyr::summarize(Tif_Path = list(tif_path_mean), .groups = "keep") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Data_Ensemble = furrr::future_pmap(
      .l = list(ias_id, time_period, climate_scenario, Tif_Path),
      .f = function(ias_id, time_period, climate_scenario, Tif_Path) {

        Name_Ensemble <- paste0(
          "Ensemble_", time_period, "_", climate_scenario) %>%
          stringr::str_replace_all("-", "_")

        Maps <- terra::rast(purrr::map(Tif_Path, terra::rast))

        Map_mean <- terra::app(Maps, fun = "mean")
        Path_mean <- file.path(
          Path_Predictions, Name_Ensemble, paste0(ias_id, "_mean.tif"))
        terra::writeRaster(
          x = Map_mean, overwrite = TRUE, filename = Path_mean)

        Map_sd <- terra::app(Maps, fun = "sd")
        Path_sd <- file.path(
          Path_Predictions, Name_Ensemble, paste0(ias_id, "_sd.tif"))
        terra::writeRaster(x = Map_sd, overwrite = TRUE, filename = Path_sd)

        Map_cov <- Map_sd / Map_mean
        Path_cov <- file.path(
          Path_Predictions, Name_Ensemble, paste0(ias_id, "_cov.tif"))
        terra::writeRaster(x = Map_cov, overwrite = TRUE, filename = Path_cov)

        return(list(
          tif_path_mean = Path_mean, tif_path_sd = Path_sd,
          tif_path_cov = Path_cov))
      },
      .options = furrr::furrr_options(
        seed = TRUE, scheduling = Inf, globals = "Path_Predictions",
        packages = c("purrr", "IASDT.R", "terra")))) %>%
  dplyr::select(-"Tif_Path") %>%
  tidyr::unnest_wider("Data_Ensemble") %>%
  dplyr::left_join(
    dplyr::distinct(
      Predictions, time_period, climate_scenario, ias_id, taxon_name,
      class, order, family, species_name),
    by = c("time_period", "climate_scenario", "ias_id"))

if (NCores > 1) {
  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)
}

Predictions_Summary <- dplyr::bind_rows(Predictions, Ensemble) %>%
  dplyr::select(
    hab_abb, hab_name, time_period, climate_model, climate_scenario, ias_id,
    taxon_name, species_name, class, order, family, tif_path_mean,
    tif_path_sd, tif_path_cov)

# Save Predictions_Summary as text file
utils::write.table(
  x = Predictions_Summary,
  file = file.path(Path_Predictions, "Predictions_Summary.txt"),
  sep = "\t", row.names = FALSE, col.names = TRUE,
  quote = FALSE, fileEncoding = "UTF-8")

# Save Predictions_Summary as RData
save(
  Predictions_Summary,
  file = file.path(Path_Predictions, "Predictions_Summary.RData"))

# # ................................................................... ###

IASDT.R::CatDiff(
  InitTime = .StartTime, Prefix = "\nPrediction took ")

return(Predictions_Summary)
}
