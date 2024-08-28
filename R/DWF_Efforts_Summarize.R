## |------------------------------------------------------------------------| #
# Efforts_Summarize ----
## |------------------------------------------------------------------------| #

#' Process and Summarize GBIF Data for vascular plants
#'
#' This function processes GBIF data for vascular plants by extracting,
#' summarizing, and saving the data. It also creates summary maps and saves them
#' in both `RData` and `TIFF` formats.
#'
#' @param NCores Numeric. The number of cores to use for parallel processing.
#' @param Path_Efforts Character. Path where the final processed data will be
#'   saved.
#' @param Path_Efforts_Interim Character. The directory path to save interim
#'   data.
#' @param Path_Efforts_Data Character. The directory path to save detailed
#'   effort data as `RData`.
#' @param Path_Grid Character. The directory path to load the grid data.
#' @param IAS_List A list of invasive alien species keys.
#' @param Efforts_AllRequests A data frame containing the details of the GBIF
#'   download, including paths to CSV files, zip files, order, class, and total
#'   records.
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. See [Efforts_Split] for more details.
#' @param DeleteChunks logical, indicating whether to remove file chunks after
#'   processing the data. Defaults to `TRUE`.
#' @return The function returns and saves the GBIF data summary.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Efforts_Process] function.
#' @author Ahmed El-Gabbas
#' @name Efforts_Summarize
#' @export

Efforts_Summarize <- function(
    NCores, Path_Efforts, Path_Efforts_Interim, Path_Efforts_Data, Path_Grid,
    IAS_List, Efforts_AllRequests, ChunkSize = 100000, DeleteChunks = TRUE) {

  .StartTimeProcess <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  if (!is.numeric(NCores) || length(NCores) != 1 || NCores <= 0) {
    stop("NCores must be a single positive integer.", call. = FALSE)
  }

  if (!is.character(Path_Efforts_Data) || length(Path_Efforts_Data) != 1 ||
      !dir.exists(Path_Efforts_Data)) {
    stop(
      "Path_Efforts_Data must be a single valid directory path.",
      call. = FALSE
    )
  }

  if (missing(Path_Grid) || !is.character(Path_Grid) ||
      !dir.exists(Path_Grid)) {
    stop("Path_Grid must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Efforts_Interim) || !is.character(Path_Efforts_Interim) ||
      !dir.exists(Path_Efforts_Interim)) {
    stop("Path_Efforts_Interim must be a valid directory path.", call. = FALSE)
  }

  if (missing(Path_Efforts) || !is.character(Path_Efforts) ||
      !dir.exists(Path_Efforts)) {
    stop("Path_Efforts must be a valid directory path.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  speciesKey <- CellCode <- ObsN <- year <- UncertainKm <- Latitude <-
    Longitude <- taxonRank <- ID <- DownPath <- Chunks <- TotalRecords <-
    ClassOrder <- Path_DT <- NULL

  # # ..................................................................... ###

  Path_Grid_R <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  Path_Grid_SF <- file.path(Path_Grid, "Grid_10_Land_Crop_sf.RData")

  if (!file.exists(Path_Grid_R)) {
    stop(
      paste0(
        "Reference grid was not found at the specified path: ", Path_Grid_R
      ),
      call. = FALSE
    )
  }

  if (!file.exists(Path_Grid_SF)) {
    stop(
      paste0("The path for the reference grid does not exist: ", Path_Grid_SF),
      call. = FALSE
    )
  }

  Grid_SF <- IASDT.R::LoadAs(Path_Grid_SF)

  # # ..................................................................... ###

  # Prepare working on parallel -----

  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  withr::local_options(future.globals.maxSize = 8000 * 1024^2)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  invisible(
    snow::clusterEvalQ(
      cl = c1,
      IASDT.R::LoadPackages(
        List = c("terra", "IASDT.R", "stringr", "fs", "sf", "readr", "dplyr"))))

  snow::clusterExport(
    cl = c1,
    list = c(
      "Path_Efforts", "Path_Efforts_Interim", "Efforts_AllRequests",
      "Path_Grid_R", "Path_Efforts_Data", "Grid_SF", "IAS_List", "ChunkSize"),
    envir = environment())

  # # ..................................................................... ###

  # Cleaning data from zipped archives -----
  IASDT.R::CatTime("Cleaning data from zipped archives", Level = 1)



  # Earlier attempts with `furrr::future_map()` failed

  DT_Paths <- future.apply::future_lapply(
    X = seq_len(nrow(Efforts_AllRequests)),
    FUN = function(ID) {

      DownPath <- Efforts_AllRequests$DownPath[ID]
      TotalRecords <- Efforts_AllRequests$TotalRecords[ID]
      class <- Efforts_AllRequests$class[ID]
      order <- Efforts_AllRequests$order[ID]
      ClassOrder <- paste0(class, "_", order)

      Path_DT <- file.path(Path_Efforts_Data, paste0(ClassOrder, ".RData"))
      ReturnNoData <- FALSE
      Done <- FALSE

      # Check if the RData file for the current order exists and valid
      if (file.exists(Path_DT)) {
        if (IASDT.R::CheckRData(Path_DT)) {
          Done <- TRUE
        } else {
          # If the file exists but not valid, delete the file and reprocess it
          fs::file_delete(Path_DT)
        }
      }

      if (isFALSE(Done) && TotalRecords > 0) {

        # Check if previous chunk files for the current order exist and contain
        # the total number of observations. If this is true, do not split the data
        # and use the chunk files directly; otherwise, split the data first into
        # small chunks. This helps to continue working on the same data should
        # previous function try failed.
        SplitChunks <- TRUE
        Chunks <- list.files(
          path = Path_Efforts_Interim,
          pattern = stringr::str_remove(basename(DownPath), ".zip"),
          full.names = TRUE)

        if (length(Chunks) > 0) {

          # Total number of lines in all chunk files
          NLines <- sum(purrr::map_int(Chunks, R.utils::countLines))

          # if they are less than the total records, delete them and recreate the
          # chunk files
          if (NLines < TotalRecords) {
            fs::file_delete(Chunks)
            rm(Chunks)
          }

          # if all records are in the chunk files, skip the splitting
          if (NLines == TotalRecords) {
            SplitChunks <- FALSE
          }
        } else {
          ReturnNoData <- TRUE
        }

        if (SplitChunks) {
          # Split data into chunks
          Chunks <- IASDT.R::Efforts_Split(
            Path_Zip = DownPath, Path_Output = Path_Efforts_Interim,
            ChunkSize = ChunkSize)
        }

        if (length(Chunks) == 0) {
          Done <- TRUE
        }

        # Processing order data for a maximum of 5 attempts
        attempt <- 0

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

        while (attempt <= 5 || isFALSE(Done)) {
          attempt <- attempt + 1
          DT <- try(
            expr = {
              purrr::map_dfr(
                .x = Chunks,
                .f = ~ {
                  readr::read_tsv(
                    file = .x, col_names = ColNames, progress = FALSE,
                    show_col_types = FALSE, col_types = Col_Types) %>%
                    dplyr::mutate(UncertainKm = UncertainKm / 1000) %>%
                    dplyr::filter(
                      !is.na(Latitude), !is.na(Longitude),
                      speciesKey != "", taxonRank %in% AcceptedRanks,
                      UncertainKm <= 100 | is.na(UncertainKm))
                })},
            silent = TRUE)

          if (inherits(DT, "try-error")) {
            next
          }

          if (nrow(DT) == 0) {
            Path_DT <- NA_character_
            break
          }

          DT <- try(
            expr = {
              sf::st_as_sf(
                x = DT, coords = c("Longitude", "Latitude"),
                crs = 4326, remove = FALSE) %>%
                sf::st_transform(3035) %>%
                sf::st_join(Grid_SF) %>%
                dplyr::filter(magrittr::not(is.na(CellCode)))
            },
            silent = TRUE)


          if (inherits(DT, "try-error")) {
            next
          }

          if (nrow(DT) == 0) {
            Path_DT <- NA_character_
            break
          }

          IASDT.R::SaveAs(InObj = DT, OutObj = ClassOrder, OutPath = Path_DT)

          if (file.exists(Path_DT)) {
            if (IASDT.R::CheckRData(Path_DT)) {
              break
            } else {
              next
            }
          } else {
            next
          }
        }
        # End of while loop
      }

      # delete chunk files for the current order
      if (DeleteChunks) {
        fs::file_delete(Chunks)
      }

      if (ReturnNoData) {
        return(
          tibble::tibble(
            ClassOrder = ClassOrder, Path_DT = NA_character_,
            NPoints = 0L, class = class, order = order)
        )
      } else {
        return(
          tibble::tibble(
            ClassOrder = ClassOrder, Path_DT = Path_DT,
            NPoints = nrow(Path_DT), class = class, order = order))
      }
    },
    future.scheduling = Inf, future.seed = TRUE) %>%
    dplyr::bind_rows()

  # # ++++++++++++++++++++++++++++++ ###

  # join data with requests summary
  RequestsCols <- c(
    "class", "order", "Request", "DownLink", "TotalRecords", "DownPath")

  Efforts_Summary <- Efforts_AllRequests %>%
    dplyr::select(tidyselect::all_of(RequestsCols)) %>%
    dplyr::left_join(DT_Paths, by = c("class", "order"))

  rm(DT_Paths)
  invisible(gc())

  # # ..................................................................... ###


  # Number of observations and species per order ----
  IASDT.R::CatTime("Number of observations and species per order", Level = 1)

  SummaryMaps <- future.apply::future_lapply(
    X = seq_len(nrow(Efforts_Summary)),
    FUN = function(ID) {

      Path_DT <- Efforts_Summary$Path_DT[ID]
      ClassOrder <- Efforts_Summary$ClassOrder[ID]

      if (is.na(Path_DT)) {
        ObsN <- ObsN_Native <- 0
      } else {
        DT <- IASDT.R::LoadAs(Path_DT)
        ObsN <- nrow(DT)

        DT_Native <- dplyr::filter(DT, !(speciesKey %in% IAS_List))
        ObsN_Native <- nrow(DT_Native)
      }

      Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

      # # ++++++++++++++++++++++++++++++++++ ###

      # All species ----

      if (ObsN == 0) {
        NObs_R <- NSp_R <- terra::classify(Grid_R, cbind(1, 0))
        NObs_R <- stats::setNames(NObs_R, paste0("NObs_", ClassOrder)) %>%
          terra::wrap()
        NSp_R <- stats::setNames(NSp_R, paste0("NSp_", ClassOrder)) %>%
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
      }
      rm(DT)

      # # ++++++++++++++++++++++++++++++++++ ###

      # Only native species ----
      if (ObsN_Native == 0) {
        NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
        NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
        NObs_Native_R <- stats::setNames(
          NObs_Native_R, paste0("NObsNative_", ClassOrder)) %>%
          terra::wrap()
        NSp_Native_R <- stats::setNames(NSp_Native_R, paste0("NSpNative_", ClassOrder)) %>%
          terra::wrap()
      } else {
        # Number of observations of native species
        NObs_Native_R <- Efforts_SummarizeMaps(
          Data = DT_Native, NSp = FALSE, Name = "NObs_Native",
          ClassOrder = ClassOrder,
          Grid_SF = Grid_SF, Grid_R = Grid_R)

        # Number of native species
        NSp_Native_R <- Efforts_SummarizeMaps(
          Data = DT_Native, NSp = TRUE, Name = "NSp_Native",
          ClassOrder = ClassOrder,
          Grid_SF = Grid_SF, Grid_R = Grid_R)
      }
      rm(DT_Native)

      # # ++++++++++++++++++++++++++++++++++ ###

      return(
        tibble::tibble(
          ClassOrder = ClassOrder,
          ObsN = ObsN, NObs_R = list(NObs_R), NSp_R = list(NSp_R),
          ObsN_Native = ObsN_Native, NObs_Native_R = list(NObs_Native_R),
          NSp_Native_R = list(NSp_Native_R)))
    },
    future.scheduling = Inf, future.seed = TRUE) %>%
    dplyr::bind_rows()

  # join data with requests summary
  Efforts_Summary <- dplyr::left_join(
    Efforts_Summary, SummaryMaps, by = c("ClassOrder"))

  rm(SummaryMaps)
  invisible(gc())

  # # ..................................................................... ###

  # Stopping cluster ------
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)
  invisible(gc())

  # # ..................................................................... ###

  # Save Efforts_Summary ----
  IASDT.R::CatTime("Save `Efforts_Summary`", Level = 1)
  save(Efforts_Summary, file = file.path(Path_Efforts, "Efforts_Summary.RData"))

  # # ..................................................................... ###

  # Prepare summary maps ----
  IASDT.R::CatTime("Prepare summary maps", Level = 1)
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

  # # ..................................................................... ###

  # Save summary data as RData ----
  IASDT.R::CatTime("Save as RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Efforts_SummaryR), OutObj = "Efforts_SummaryR",
    OutPath = file.path(Path_Efforts, "Efforts_SummaryR.RData"))

  # # ..................................................................... ###

  # Save summary data as tif ----
  IASDT.R::CatTime("Save as tif", Level = 1)
  terra::writeRaster(
    Efforts_SummaryR,
    overwrite = TRUE,
    filename = file.path(
      Path_Efforts, paste0("Efforts_GBIF_", names(Efforts_SummaryR), ".tif")))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTimeProcess, CatInfo = FALSE,
    Prefix = "Processing Efforts data took ", Level = 1)

  # # ..................................................................... ###

  return(invisible(NULL))
}
