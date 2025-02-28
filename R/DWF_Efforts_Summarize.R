## |------------------------------------------------------------------------| #
# Efforts_Summarize ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 4
#' @export

Efforts_Summarize <- function(
    FromHPC = TRUE, EnvFile = ".env", NCores = 6L, ChunkSize = 100000L,
    DeleteChunks = TRUE) {

  # # ..................................................................... ###

  .StartTimeProcess <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  speciesKey <- ObsN <- Path_Interim <- Taxa_Stand <- Path_Grid <-
    Path_Efforts <- NULL

  # # ..................................................................... ###

  if (!is.numeric(NCores) || length(NCores) != 1 || NCores <= 0) {
    stop("NCores must be a single positive integer.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts", TRUE, FALSE,
      "Path_Interim", "DP_R_Efforts_Interim", TRUE, FALSE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Taxa_Stand", "DP_R_TaxaStand", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Efforts", "DP_R_Efforts_Local", TRUE, FALSE,
      "Path_Interim", "DP_R_Efforts_Interim_Local", TRUE, FALSE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Taxa_Stand", "DP_R_TaxaStand_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  Path_Efforts_Cleaned <- IASDT.R::Path(Path_Interim, "CleanedData")

  ## IAS species list ----
  IAS_List <- readRDS(Taxa_Stand) %>%
    dplyr::pull("speciesKey") %>%
    unique()

  # # ..................................................................... ###

  Path_Grid_R <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  Path_Grid_SF <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf.RData")

  if (!file.exists(Path_Grid_R)) {
    stop(
      "Reference grid was not found at the specified path: ", Path_Grid_R,
      call. = FALSE)
  }

  if (!file.exists(Path_Grid_SF)) {
    stop(
      "Reference grid file was not found at the specified path: ",
      Path_Grid_SF, call. = FALSE)
  }

  Grid_SF <- IASDT.R::LoadAs(Path_Grid_SF)

  # # ..................................................................... ###

  # Prepare working on parallel -----

  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Processing data from zipped archives -----
  IASDT.R::CatTime("Processing data from zipped archives", Level = 1)

  if (DeleteChunks) {
    IASDT.R::CatTime(
      "Chunk files will be deleted after finishing processing",
      Level = 2)
  }

  # Earlier attempts with `furrr::future_map()` failed

  Path_Efforts_Request <- IASDT.R::Path(
    Path_Efforts, "Efforts_AllRequests.RData")

  if (!file.exists(Path_Efforts_Request)) {
    stop(
      "The path for the `Efforts_AllRequests` data does not exist: ",
      Path_Efforts_Request, call. = FALSE)
  }

  Efforts_AllRequests <- IASDT.R::LoadAs(Path_Efforts_Request)

  DT_Paths <- future.apply::future_lapply(
    X = seq_len(nrow(Efforts_AllRequests)),
    FUN = function(ID) {
      DownPath <- Efforts_AllRequests$DownPath[ID]
      TotalRecords <- Efforts_AllRequests$TotalRecords[ID]
      class <- Efforts_AllRequests$class[ID]
      order <- Efforts_AllRequests$order[ID]
      ClassOrder <- paste0(class, "_", order)

      # Output path to save the data
      Path_DT <- IASDT.R::Path(
        Path_Efforts_Cleaned, paste0(ClassOrder, ".RData"))

      # Should Path_DT be returned as the path of the RData file containing the
      # data or NA if there are no records in the current order or no records
      # left after processing
      ReturnNoData <- (TotalRecords == 0)

      # Check if the RData file for the current order exists and valid.
      FileOkay <- IASDT.R::CheckData(Path_DT, warning = FALSE)

      # Process current order data if the output file is not okay and the order
      # have observations
      if (isFALSE(FileOkay) && TotalRecords > 0) {
        if (file.exists(Path_DT)) {
          fs::file_delete(Path_DT)
        }

        # Check if previous chunk files for the current order exist and contain
        # the same total number of observations. If this is true, do not split
        # the data and use the chunk files directly; otherwise, split the data
        # first into small chunks. This helps to continue working on the same
        # data should previous function try failed.

        # List of chunks
        Chunks <- list.files(
          path = Path_Interim, full.names = TRUE,
          pattern = stringr::str_remove(basename(DownPath), ".zip"))

        # If there are chunk files on disk, count their total number of
        # observations
        if (length(Chunks) > 0) {
          # Total number of lines in all chunk files
          NLines <- sum(purrr::map_int(Chunks, R.utils::countLines))

          # if there are less than the total number of records, delete the chunk
          # files and recreate them
          if (NLines != TotalRecords) {
            purrr::walk(Chunks, file.remove)
            SplitChunks <- TRUE
            rm(Chunks, envir = environment())
          } else {
            SplitChunks <- FALSE
          }
        } else {
          # If there is no chunk files available, split into chunks
          SplitChunks <- TRUE
        }


        # Split data into chunks
        if (SplitChunks) {
          Chunks <- IASDT.R::Efforts_Split(
            Path_Zip = DownPath, FromHPC = TRUE, EnvFile = ".env",
            ChunkSize = ChunkSize)
        }

        # Process chunk files
        # nolint start
        AcceptedRanks <- c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")
        ColNames <- c(
          "taxonRank", "Latitude", "Longitude", "UncertainKm", "speciesKey")
        Col_Types <- readr::cols(
          UncertainKm = readr::col_double(),
          Longitude = readr::col_double(),
          Latitude = readr::col_double(),
          speciesKey = readr::col_integer(),
          taxonRank = readr::col_character(),
          .default = readr::col_double())
        # nolint end

        DT <- purrr::map_dfr(
          .x = Chunks,
          .f = ~ {
            readr::read_tsv(
              file = .x, col_names = ColNames, progress = FALSE,
              show_col_types = FALSE, col_types = Col_Types) %>%
              dplyr::mutate(UncertainKm = UncertainKm / 1000) %>%
              dplyr::filter(
                !is.na(Latitude), !is.na(Longitude),
                nzchar(speciesKey), taxonRank %in% AcceptedRanks,
                UncertainKm <= 100 | is.na(UncertainKm)) %>%
              sf::st_as_sf(
                coords = c("Longitude", "Latitude"),
                crs = 4326, remove = FALSE) %>%
              sf::st_transform(3035) %>%
              sf::st_join(Grid_SF) %>%
              dplyr::filter(magrittr::not(is.na(CellCode)))
          }
        )

        # if there are observations after the filtering, save the data to disk
        # and return the saved path, otherwise return no path
        if (nrow(DT) > 0) {
          IASDT.R::SaveAs(InObj = DT, OutObj = ClassOrder, OutPath = Path_DT)
        } else {
          ReturnNoData <- TRUE
        }

        rm(DT, envir = environment())
      }

      # delete chunk files for the current order
      if (DeleteChunks) {
        purrr::walk(Chunks, file.remove)
      }

      # Output path
      Path_DT <- dplyr::if_else(ReturnNoData, NA_character_, Path_DT)

      invisible(gc())

      return(
        tibble::tibble(
          ClassOrder = ClassOrder, Path_DT = Path_DT,
          class = class, order = order))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = c(
      "terra", "IASDT.R", "stringr", "fs", "sf", "readr", "dplyr",
      "purrr", "tibble", "R.utils"),
    future.globals = c(
      "Path_Interim", "Efforts_AllRequests", "Path_Efforts_Cleaned",
      "Grid_SF", "ChunkSize")
  ) %>%
    dplyr::bind_rows()

  # # ++++++++++++++++++++++++++++++ ###

  # only selected columns from `Efforts_AllRequests`
  RequestsCols <- c(
    "class", "order", "Request", "DownLink", "TotalRecords", "DownPath")

  # join data with requests summary
  Efforts_Summary <- Efforts_AllRequests %>%
    dplyr::select(tidyselect::all_of(RequestsCols)) %>%
    dplyr::left_join(DT_Paths, by = c("class", "order"))

  rm(DT_Paths, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Prepare summary maps per order ----
  IASDT.R::CatTime("Prepare summary maps per order", Level = 1)

  SummaryMaps <- future.apply::future_lapply(
    X = seq_len(nrow(Efforts_Summary)),
    FUN = function(ID) {
      Path_DT <- Efforts_Summary$Path_DT[ID]
      ClassOrder <- Efforts_Summary$ClassOrder[ID]

      if (is.na(Path_DT)) {
        # If there is no data for the current order
        ObsN <- ObsN_Native <- 0L
      } else {
        # Load data on current order
        DT <- IASDT.R::LoadAs(Path_DT)
        # Number of observation in the cleaned data
        ObsN <- nrow(DT)

        # Only native species (exclude IAS list)
        DT_Native <- dplyr::filter(DT, !(speciesKey %in% IAS_List))
        # Number of data for native species
        ObsN_Native <- nrow(DT_Native)
      }

      Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

      # # ++++++++++++++++++++++++++++++++++ ###

      # All species ----

      if (ObsN == 0) {
        # Create dummy maps for the number of species and records
        NObs_R <- terra::classify(Grid_R, cbind(1, 0)) %>%
          stats::setNames(paste0("NObs_", ClassOrder)) %>%
          terra::wrap()

        NSp_R <- terra::classify(Grid_R, cbind(1, 0)) %>%
          stats::setNames(paste0("NSp_", ClassOrder)) %>%
          terra::wrap()
      } else {
        # Number of observations
        NObs_R <- Efforts_SummarizeMaps(
          Data = DT, NSp = FALSE, Name = "NObs", ClassOrder = ClassOrder,
          Grid_SF = Grid_SF, Grid_R = Grid_R)

        # Number of species
        NSp_R <- Efforts_SummarizeMaps(
          Data = DT, NSp = TRUE, Name = "NSp", ClassOrder = ClassOrder,
          Grid_SF = Grid_SF, Grid_R = Grid_R)

        rm(DT, envir = environment())
      }

      # # ++++++++++++++++++++++++++++++++++ ###

      # Only native species ----

      if (ObsN_Native == 0) {
        # Create dummy maps for the number of species and records

        NObs_Native_R <- terra::classify(Grid_R, cbind(1, 0)) %>%
          stats::setNames(paste0("NObsNative_", ClassOrder)) %>%
          terra::wrap()

        NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0)) %>%
          stats::setNames(paste0("NSpNative_", ClassOrder)) %>%
          terra::wrap()
      } else {
        # Number of observations of native species
        NObs_Native_R <- Efforts_SummarizeMaps(
          Data = DT_Native, NSp = FALSE, Name = "NObs_Native",
          ClassOrder = ClassOrder, Grid_SF = Grid_SF, Grid_R = Grid_R)

        # Number of native species
        NSp_Native_R <- Efforts_SummarizeMaps(
          Data = DT_Native, NSp = TRUE, Name = "NSp_Native",
          ClassOrder = ClassOrder, Grid_SF = Grid_SF, Grid_R = Grid_R)

        rm(DT_Native, envir = environment())
      }

      # # ++++++++++++++++++++++++++++++++++ ###

      return(
        tibble::tibble(
          ClassOrder = ClassOrder,
          ObsN = ObsN, NObs_R = list(NObs_R), NSp_R = list(NSp_R),
          ObsN_Native = ObsN_Native, NObs_Native_R = list(NObs_Native_R),
          NSp_Native_R = list(NSp_Native_R)))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = c(
      "terra", "IASDT.R", "stringr", "fs", "sf", "readr", "dplyr"),
    future.globals = c("Path_Grid_R", "Grid_SF", "IAS_List")) %>%
    dplyr::bind_rows()

  # join data with requests summary
  Efforts_Summary <- dplyr::left_join(
    Efforts_Summary, SummaryMaps, by = "ClassOrder")

  rm(SummaryMaps, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }
  invisible(gc())

  # # ..................................................................... ###

  # Save summary results: `Efforts_Summary` ----
  IASDT.R::CatTime("Save summary results: `Efforts_Summary`", Level = 1)
  save(
    Efforts_Summary,
    file = IASDT.R::Path(Path_Efforts, "Efforts_Summary.RData"))

  # # ..................................................................... ###

  # Prepare summary maps - all sampling efforts ----
  IASDT.R::CatTime("Prepare summary maps - all sampling efforts", Level = 1)

  CalcNObsNSp <- function(List, Name) {
    purrr::map(.x = unlist(List), .f = terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE) %>%
      IASDT.R::setRastCRS() %>%
      IASDT.R::setRastVals() %>%
      stats::setNames(Name)
  }

  # # ..................................................................... ###

  # Exclude orders with no data
  Efforts_SummaryR <- dplyr::filter(Efforts_Summary, ObsN > 0)

  Efforts_SummaryR <- list(
    CalcNObsNSp(Efforts_SummaryR$NObs_R, "NObs"),
    CalcNObsNSp(Efforts_SummaryR$NObs_Native_R, "NObs_Native"),
    CalcNObsNSp(Efforts_SummaryR$NSp_R, "NSp"),
    CalcNObsNSp(Efforts_SummaryR$NSp_Native_R, "NSp_Native")) %>%
    terra::rast() %>%
    IASDT.R::setRastCRS() %>%
    IASDT.R::setRastVals()

  ## Save summary maps - `RData` ----
  IASDT.R::CatTime("Save summary maps", Level = 1)

  IASDT.R::CatTime("`RData`", Level = 2)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Efforts_SummaryR), OutObj = "Efforts_SummaryR",
    OutPath = IASDT.R::Path(Path_Efforts, "Efforts_SummaryR.RData"))

  ## Save summary maps - `tif` ----
  IASDT.R::CatTime("`tif`", Level = 2)
  terra::writeRaster(
    Efforts_SummaryR,
    overwrite = TRUE,
    filename = IASDT.R::Path(
      Path_Efforts, paste0("Efforts_GBIF_", names(Efforts_SummaryR), ".tif")))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeProcess,
    Prefix = "Processing Efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}



## |------------------------------------------------------------------------| #
# Efforts_SummarizeMaps ----
## |------------------------------------------------------------------------| #

#' Summarize maps for efforts data
#'
#' This function processes spatial data (as an `sf` object), summarizes it based
#' on the number of observations or distinct species, and generates a raster
#' layer.
#' @param Data An `sf` object containing spatial data, with a column named
#'   `CellCode`.
#' @param NSp Logical. Whether to generate distinct species counts (`TRUE`) or
#'   total observation counts (`FALSE`).
#' @param Name Character. Name of the count field and the prefix for the final
#'   raster layer's name.
#' @param ClassOrder Character. The class and order combination (separated by
#'   an underscore) represented in the `Data`.
#' @param Grid_SF,Grid_R Reference grid in the form of simple feature and
#'   raster.
#' @return A processed `terra` raster object representing the summarized data.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Efforts_Process] and [Efforts_Summarize]
#'   functions.
#' @author Ahmed El-Gabbas
#' @name Efforts_SummarizeMaps
#' @keywords internal
#' @noRd

Efforts_SummarizeMaps <- function(
    Data, NSp, Name, ClassOrder, Grid_SF, Grid_R) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CellCode <- speciesKey <- NULL

  # # ..................................................................... ###

  # Validate if Data is an sf object
  if (!inherits(Data, "sf")) {
    stop(
      "Input data must be a simple feature (sf) object. ",
      "Provided data is of type: ", paste(class(Data), collapse = "+"),
      call. = FALSE)
  }

  # Validate if NSp is logical
  if (!is.logical(NSp) || length(NSp) != 1) {
    stop(
      "The parameter `NSp` must be a single logical value (TRUE or FALSE). ",
      "Provided value is of type: ", paste(class(NSp), collapse = "+"),
      call. = FALSE)
  }

  # Validate the Name parameter
  if (is.null(Name)) {
    stop("The parameter `Name` can not be empty", call. = FALSE)
  }

  # # ..................................................................... ###

  # Drop geometry from Data
  Data <- sf::st_drop_geometry(Data)

  # Generate distinct species counts if NSp is TRUE
  if (NSp) {
    Data <- dplyr::distinct(Data, CellCode, speciesKey)
  }

  # Count observations or species, join with the grid, and rasterize
  Data <- Data %>%
    dplyr::count(CellCode, name = Name) %>%
    dplyr::left_join(Grid_SF, by = "CellCode") %>%
    sf::st_as_sf() %>%
    terra::rasterize(Grid_R, field = Name) %>%
    terra::classify(cbind(NA, 0)) %>%
    terra::mask(Grid_R) %>%
    IASDT.R::setRastCRS() %>%
    IASDT.R::setRastVals() %>%
    stats::setNames(paste0(Name, "_", ClassOrder)) %>%
    terra::wrap()

  return(Data)
}
