
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

  NullVarsNames <- c(
    "Name", "File", "Path_Model", "Path_Grid", "TimePeriod", "ClimateModel",
    "ClimateScenario", "StaticPredictors")
  NullVars <- which(purrr::map_lgl(.x = NullVarsNames, .f = ~ is.null(get(.x))))

  if (length(NullVars) > 0) {
    NullVarsNames[NullVars]
    stop(
      paste0(
        paste0(NullVarsNames[NullVars], collapse = ", "),
        " cannot be missing."),
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

  Chunk <- predictions <- LayerName <- Stats <- Chunk_File <-
    ChunkExist <- ChunkData <- NULL

  # # ..................................................................... ###

  # Check if the results of the current climate option were already finished
  Name2 <- stringr::str_remove(Name, "^R_")
  Path_Predictions <- file.path(
    dirname(dirname(Path_Model)), "Predictions", Name2)
  fs::dir_create(Path_Predictions)

  Path_Metadata <- file.path(Path_Predictions, "Predictions_Metadata.RData")
  Path_Metadata_DT <- file.path(
    Path_Predictions, "Predictions_Metadata_DT.RData")
  Path_Predictions_R <- file.path(Path_Predictions, "Predictions_R.RData")
  Path_Predictions_sf <- file.path(Path_Predictions, "Predictions_sf.RData")

  AllExported <- purrr::map_lgl(
    .x = c(
      Path_Metadata, Path_Metadata_DT, Path_Predictions_R, Path_Predictions_sf),
    .f = ~ file.exists(.x) && IASDT.R::CheckRData(.x)) %>%
    all()

  if (AllExported) {
    return(Path_Metadata)
  }

  IASDT.R::InfoChunk(
    Message = paste0("Climate option: ", Name2),
    Date = TRUE, Extra1 = 1, Extra2 = 0)

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

    if (NCores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      c1 <- snow::makeSOCKcluster(min(NCores, nrow(IncompleteChunks)))
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
    }


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
    rm(Predictions_sf0)

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("sequential", gc = TRUE)
    }

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
  save(Predictions_sf, file = Path_Predictions_sf)

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

  IASDT.R::CatTime("Predictions metadata", Level = 2)

  IASDT.R::CatTime("Prepare predictions metadata", Level = 3)

  SpInfo <- IASDT.R::GetSpeciesName(
    SpID = NULL, EnvFile = EnvFile, FromHPC = FromHPC) %>%
    dplyr::select(
      -tidyselect::all_of(
        c("Species_name2", "Species_File", "Genus", "Species"))) %>%
    janitor::clean_names()

  hab_name <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_remove(paste0("^", as.character(Hab_Abb), "_")) %>%
    stringr::str_replace_all("_", " ")

  Predictions_Metadata <- terra::as.list(Predictions_R) %>%
    tibble::tibble(predictions = .) %>%
    dplyr::mutate(
      hab_abb = Hab_Abb,
      hab_name = hab_name,
      time_period = TimePeriod,
      climate_model = ClimateModel,
      climate_scenario = ClimateScenario,
      LayerName = purrr::map_chr(.x = predictions, .f = names),
      predictions = purrr::map(predictions, terra::wrap),
      tif_path = file.path(Path_Predictions, paste0(LayerName, ".tif")),

      Stats = purrr::map_chr(
        .x = LayerName,
        .f = ~ {
          Stats <- c("mean", "sd", "cov")
          Stats[stringr::str_detect(.x, Stats)]
        }),

      ias_id = purrr::map2_chr(
        .x = LayerName, .y = Stats,
        .f = ~ stringr::str_remove_all(.x, paste0("_", .y)))) %>%
    dplyr::left_join(SpInfo, by = "ias_id") %>%
    tidyr::pivot_wider(
      id_cols = c(
        "hab_abb", "time_period", "climate_model", "climate_scenario",
        "ias_id", "hab_name", "taxon_name", "class", "order", "family",
        "species_name"),
      names_from = "Stats", values_from = c("predictions", "tif_path"))

  IASDT.R::CatTime("Save predictions metadata", Level = 3)
  save(Predictions_Metadata, file = Path_Metadata)

  # # ..................................................................... ###

  # Save prediction metadata without maps -----
  IASDT.R::CatTime("Save prediction metadata without maps", Level = 3)

  Predictions_Metadata_DT <- dplyr::select(
    Predictions_Metadata, -tidyselect::starts_with("Predictions"))

  save(Predictions_Metadata_DT, file = Path_Metadata_DT)

  readr::write_csv(
    x = Predictions_Metadata_DT,
    file = file.path(Path_Predictions, "Predictions_Metadata_DT.csv"))

  # # ..................................................................... ###

  Predictions_R <- terra::wrap(Predictions_R)
  save(Predictions_R, file = Path_Predictions_R)

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
