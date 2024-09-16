
## |------------------------------------------------------------------------| #
# Predict_Scenario ----
## |------------------------------------------------------------------------| #

#' Predict from Hmsc model at a single climate option
#'
#' This function generates predictions across the full study area for a single
#' climate option, currently and in the future. The function predicts the mean,
#' standard deviation, and coefficient of variation of the level of invasion
#' and species-specific predictions. The function is intended to be used only
#' within the main prediction function: [Mod_Predict], not to be called by the
#' user.
#' @param Name Character. A combination of time / climate model / climate
#'   scenario.
#' @param File Character. Path to the `.RData` file containing raster maps for
#'   the respective climate option.
#' @param ChunkSize Integer. The size of chunks used to split the input raster
#'   for parallel processing. Predictor rasters are converted into data frame
#'   (tibble) then split into small chunks for processing on parallel.
#'   Default: `200`.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing of
#'   chunks. Defaults to 8.
#' @param StaticPredictors `SpatRaster` for static predictors (predictors
#'   other than bioclimatic predictors).
#' @param Path_Model Character. Path to the fitted `Hmsc` model object.
#' @param Path_Grid Path to the reference grid used for rasterizing the
#'   predictions.
#' @param RemoveChunks Logical. Whether to remove intermediate chunk files after
#'   prediction. Defaults to `TRUE`.
#' @param TimePeriod Character. The time frame of the climate data.
#' @param ClimateModel Character. The climate model of the climate data.
#' @param ClimateScenario Character. The climate change scenario of the climate
#'   data.
#' @return Path to the `RData` file containing prediction details.
#' @noRd

Predict_Scenario <- function(
    Name, File, Hab_Abb = NULL, ChunkSize = 200, NCores = 8,
    StaticPredictors = NULL, Path_Model = NULL, Path_Grid = NULL,
    EnvFile = ".env", FromHPC = TRUE, RemoveChunks = TRUE, TimePeriod = NULL,
    ClimateModel = NULL, ClimateScenario = NULL) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  if (is.null(Name) || is.null(File) || is.null(Path_Model) ||
      is.null(Path_Grid) || is.null(TimePeriod) || is.null(ClimateModel) ||
      is.null(ClimateScenario) || is.null(StaticPredictors)) {
    stop(
      paste0(
        "Required parameters: `Name`, `File`, `Path_Model`, ",
        "`Path_Grid`, `TimePeriod`, `ClimateModel`, and ",
        "`ClimateScenario` cannot be missing."),
      call. = FALSE)
  }

  if (!file.exists(File)) {
    stop("Input file does not exist.", call. = FALSE)
  }

  if (!file.exists(Path_Model)) {
    stop("Model file does not exist.", call. = FALSE)
  }
  if (!file.exists(Path_Grid)) {
    stop("Grid reference file does not exist.", call. = FALSE)
  }
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(NCores) || NCores <= 0) {
    stop("NCores must be a positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  Chunk <- Predictions <- LayerName <- Stats <- Chunk_File <-
    ChunkExist <- ChunkData <- Species_ID <- IAS_ID <- NULL

  # # ..................................................................... ###

  Name2 <- stringr::str_remove(Name, "^R_")
  IASDT.R::InfoChunk(
    Message = paste0("Climate option: ", Name2),
    Date = TRUE, Extra1 = 2, Extra2 = 1)

  Path_Predictions <- file.path(
    dirname(dirname(Path_Model)), "Predictions", Name2)
  fs::dir_create(Path_Predictions)

  Path_Metadata <- file.path(Path_Predictions, "Predictions_Metadata.RData")
  if (file.exists(Path_Metadata)) {
    return(Path_Metadata)
  }

  IASDT.R::CatTime("Predictions will be saved to:", Level = 2)
  IASDT.R::CatTime(Path_Predictions, Level = 3)

  # # ..................................................................... ###

  # Loading model object -----
  IASDT.R::CatTime("Loading model object", Level = 2)

  Model <- IASDT.R::LoadAs(Path_Model)
  BioVars <- stringr::str_subset(names(Model$XData), "^bio")
  rm(Model)

  # # ..................................................................... ###

  # Prepare input raster maps -----
  IASDT.R::CatTime("Prepare input raster maps", Level = 2)

  VarsR <- IASDT.R::LoadAs(File) %>%
    terra::unwrap() %>%
    terra::subset(subset = BioVars) %>%
    c(StaticPredictors)

  # # ..................................................................... ###

  # Split input maps into data.frame chunks -----
  IASDT.R::CatTime("Split input maps into data.frame chunks", Level = 2)

  Chunks <- terra::as.data.frame(x = VarsR, xy = TRUE, na.rm = TRUE) %>%
    tibble::tibble() %>%
    dplyr::mutate(
      Chunk = as.integer(ceiling(dplyr::row_number() / ChunkSize))) %>%
    tidyr::nest(ChunkData = -"Chunk") %>%
    dplyr::mutate(
      Chunk_File = purrr::map_chr(
        .x = Chunk,
        .f = ~{
          stringr::str_pad(.x, width = nchar(dplyr::n()), pad = "0") %>%
            paste0("Chunk_", ., ".RData") %>%
            file.path(Path_Predictions, .)
        }),

      ChunkExist = purrr::map_lgl(
        .x = Chunk_File, .f = ~ file.exists(.x) && IASDT.R::CheckRData(.x)))

  # # ..................................................................... ###

  # Processing chunks as simple features -----
  IASDT.R::CatTime("Processing chunks and converting to sf", Level = 2)

  # Chunks not yet processed
  IncompleteChunks <- dplyr::filter(Chunks, !ChunkExist)

  if (nrow(IncompleteChunks) > 0) {

    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)
    future::plan(
      "multisession", workers = min(NCores, nrow(IncompleteChunks)), gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)

    Predictions_sf0 <- IncompleteChunks %>%
      dplyr::mutate(
        Sf = furrr::future_map2(
          .x  = ChunkData, .y = Chunk_File,
          .f = ~{

            Model <- IASDT.R::LoadAs(Path_Model)
            ModelVars <- stringr::str_subset(
              names(Model$XData), "^CV", negate = TRUE)
            BioVars <- stringr::str_subset(names(Model$XData), "^bio")

            max_retries <- 5
            Try <- 0

            while (TRUE) {
              Try <- Try + 1
              Preds <- try(
                expr = {
                  Predict_Chunk(
                    ChunkData = .x, Model = Model, ModelVars = ModelVars)
                },
                silent = TRUE)

              if (inherits(Preds, "sf")) {
                try(save(Preds, file = .y), silent = TRUE)
                Sys.sleep(3)

                if (IASDT.R::CheckRData(.y)) {
                  invisible(gc())
                  break
                }
              }

              if (Try > max_retries) {
                warning(
                  paste0("Max retries (", max_retries, ") exceeded for Chunk ",
                         ChunkID),
                  call. = FALSE)
                invisible(gc())
                break
              }
            }

            return(invisible(NULL))
          },
          .options = furrr::furrr_options(
            seed = TRUE, scheduling = Inf,
            globals = c("Predict_Chunk", "Path_Model", "Path_Predictions"),
            packages = c(
              "Hmsc", "purrr", "dplyr", "sf", "IASDT.R", "stringr", "terra")
          )))

    future::plan("sequential")

    IncompleteChunks <- IncompleteChunks %>%
      dplyr::mutate(
        ChunkExist = purrr::map_lgl(
          .x = Chunk_File,
          .f = ~ file.exists(.x) && IASDT.R::CheckRData(.x))) %>%
      dplyr::filter(!ChunkExist)

    if (nrow(IncompleteChunks)) {
      stop(
        paste0(nrow(IncompleteChunks), " chunks failed to processed"),
        call. = FALSE)
    }
  }

  # # ..................................................................... ###

  # Merge and save chunk data into a single RData file -----
  IASDT.R::CatTime(
    "Merge and save chunk data into a single RData file", Level = 2)

  Predictions_sf <- purrr::map_dfr(Chunks$Chunk_File, IASDT.R::LoadAs)
  save(
    Predictions_sf,
    file = file.path(Path_Predictions, "Predictions_sf.RData"))

  # # ..................................................................... ###

  # Convert predictions from sf to raster -----
  IASDT.R::CatTime("Convert predictions from sf to raster", Level = 2)

  GridR <- terra::unwrap(IASDT.R::LoadAs(Path_Grid))

  Predictions_R <- terra::vect(Predictions_sf) %>%
    terra::rasterize(y = GridR, field = names(.)) %>%
    IASDT.R::setRastCRS() %>%
    IASDT.R::setRastVals()

  # # ..................................................................... ###

  # Save raster maps as tif files -----

  IASDT.R::CatTime("Save raster maps as tif files", Level = 2)

  terra::writeRaster(
    x = Predictions_R, overwrite = TRUE,
    filename = file.path(
      Path_Predictions, paste0(names(Predictions_R), ".tif")))

  # # ..................................................................... ###

  # Prepare predictions metadata -----
  IASDT.R::CatTime("Prepare predictions metadata", Level = 2)

  SpInfo <- GetSpeciesName(SpID = NULL, EnvFile = EnvFile, FromHPC = FromHPC)

  Predictions_Metadata <- terra::as.list(Predictions_R) %>%
    tibble::tibble(Predictions = .) %>%
    dplyr::mutate(
      Hab_Abb = Hab_Abb,
      TimePeriod = TimePeriod,
      ClimateModel = ClimateModel,
      ClimateScenario = ClimateScenario,
      LayerName = purrr::map_chr(.x = Predictions, .f = names),
      Predictions = purrr::map(Predictions, terra::wrap),
      Tif_Path = file.path(Path_Predictions, paste0(LayerName, ".tif")),

      Stats = purrr::map_chr(
        .x = LayerName,
        .f = ~ {
          Stats <- c("mean", "sd", "cov")
          Stats[stringr::str_detect(.x, Stats)]
        }),

      Species_ID = purrr::map2_chr(
        .x = LayerName, .y = Stats,
        .f = ~ stringr::str_remove_all(.x, paste0("_", .y)))) %>%
    dplyr::left_join(SpInfo, by = dplyr::join_by(Species_ID == IAS_ID)) %>%
    tidyr::pivot_wider(
      id_cols = c(
        "Hab_Abb", "TimePeriod", "ClimateModel", "ClimateScenario",
        "Species_ID", "taxon_name", "Class", "Order", "Family", "Species_name",
        "Species_name2", "Species_File", "Genus", "Species"),
      names_from = "Stats", values_from = c("Predictions", "Tif_Path"))

  # # ..................................................................... ###

  # Save prediction metadata without maps -----
  IASDT.R::CatTime("Save prediction metadata without maps", Level = 2)

  Predictions_Metadata_DT <- dplyr::select(
    Predictions_Metadata, -tidyselect::starts_with("Predictions"))

  save(
    Predictions_Metadata_DT,
    file = file.path(Path_Predictions, "Predictions_Metadata_DT.RData"))

  readr::write_csv(
    x = Predictions_Metadata_DT,
    file = file.path(Path_Predictions, "Predictions_Metadata_DT.csv"))

  # # ..................................................................... ###

  Predictions_R <- terra::wrap(Predictions_R)
  save(
    Predictions_R,
    file = file.path(Path_Predictions, "Predictions_R.RData"))

  # # ..................................................................... ###

  # Remove chunk files -----

  if (RemoveChunks) {
    IASDT.R::CatTime("Remove chunk files", Level = 2)

    Path_Predictions %>%
      list.files(pattern = "Chunk_.+.RData", full.names = TRUE) %>%
      fs::file_delete()
  }

  # # ..................................................................... ###

  invisible(gc())

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Prediction took ", Level = 2)

  return(Path_Metadata)
}
