## |------------------------------------------------------------------------| #
# Efforts_Process ----
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
#' @return The function returns and saves the GBIF data summary.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Sampling_Efforts] function.
#' @author Ahmed El-Gabbas
#' @name Efforts_Process
#' @export

Efforts_Process <- function(
    NCores, Path_Efforts, Path_Efforts_Interim, Path_Efforts_Data, Path_Grid,
    IAS_List, Efforts_AllRequests, ChunkSize = 100000) {
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

  # Reading data from zipped archives -----
  IASDT.R::CatTime("Reading data from zipped archives", Level = 1)

  RequestsCols <- c(
    "class", "order", "Request", "DownLink", "TotalRecords", "DownPath")

  Efforts_Summary <- Efforts_AllRequests %>%
    dplyr::select(tidyselect::all_of(RequestsCols)) %>%
    dplyr::mutate(
      ClassOrder = purrr::map2_chr(class, order, ~ paste0(.x, "_", .y)),
      Path_DT = furrr::future_pmap_chr(
        .l = list(DownPath, TotalRecords, ClassOrder),
        .f = function(DownPath, TotalRecords, ClassOrder) {
          if (TotalRecords == 0) {
            return(NA_character_)
          }

          Path_DT <- file.path(Path_Efforts_Data, paste0(ClassOrder, ".RData"))

          # Check if the RData file for the current order exists and valid
          if (file.exists(Path_DT)) {
            if (IASDT.R::CheckRData(Path_DT)) {
              # return the path if the file exists and is valid
              return(Path_DT)
            } else {
              # If the file exists but not valid, delete the file and reprocess it
              fs::file_delete(Path_DT)
            }
          }

          # Check if previous chunk files for the current order exist and
          # contain the total number of observations. If this is true, do not
          # split the data and use the chunk files directly; otherwise, split
          # the data first into small chunks. This helps to continue working on
          # the same data should previous function try failed.
          SplitChunks <- TRUE
          Chunks <- list.files(
            path = Path_Efforts_Interim,
            pattern = stringr::str_remove(basename(DownPath), ".zip"),
            full.names = TRUE)

          if (length(Chunks) > 0) {

            # Total number of lines in all chunk files
            NLines <- sum(purrr::map_int(Chunks, R.utils::countLines))

            # if they are less than the total records, delete them and recreate
            # the chunk files
            if (NLines < TotalRecords) {
              fs::file_delete(Chunks)
              rm(Chunks)
            }

            # if all records are in the chunk files, skip the splitting
            if (NLines == TotalRecords) {
              SplitChunks <- FALSE
            }
          }

          if (SplitChunks) {
            # Split data into chunks
            Chunks <- IASDT.R::Efforts_Split(
              Path_Zip = DownPath, Path_Output = Path_Efforts_Interim,
              ChunkSize = ChunkSize)
          }

          if (length(Chunks) == 0) {
            return(NA_character_)
          }

          # Initialize the attempt counter
          max_attempts <- 3
          attempt <- 1
          OrderOkay <- FALSE
          AcceptedRanks <- c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")

          while (isFALSE(OrderOkay) && (attempt <= max_attempts)) {
            tryCatch(
              {
                DT <- purrr::map_dfr(
                  .x = Chunks,
                  .f = ~ {
                    readr::read_tsv(
                      file = .x,
                      col_names = c(
                        "taxonRank", "Latitude", "Longitude",
                        "UncertainKm", "speciesKey"),
                      progress = FALSE, show_col_types = FALSE,
                      col_types = readr::cols(
                        UncertainKm = readr::col_double(),
                        Longitude = readr::col_double(),
                        Latitude = readr::col_double(),
                        speciesKey = readr::col_integer(),
                        taxonRank = readr::col_character(),
                        .default = readr::col_double())) %>%
                      dplyr::mutate(UncertainKm = UncertainKm / 1000) %>%
                      dplyr::filter(
                        !is.na(Latitude), !is.na(Longitude),
                        speciesKey != "", taxonRank %in% AcceptedRanks,
                        UncertainKm <= 100 | is.na(UncertainKm))
                  }
                )

                if (nrow(DT) == 0) {
                  return(NA_character_)
                }

                DT <- sf::st_as_sf(
                  x = DT, coords = c("Longitude", "Latitude"),
                  crs = 4326, remove = FALSE) %>%
                  # project to 3035
                  sf::st_transform(3035) %>%
                  sf::st_join(Grid_SF) %>%
                  dplyr::filter(magrittr::not(is.na(CellCode)))

                if (nrow(DT) == 0) {
                  return(NA_character_)
                }

                IASDT.R::SaveAs(
                  InObj = DT, OutObj = ClassOrder, OutPath = Path_DT)

                if (file.exists(Path_DT)) {
                  if (IASDT.R::CheckRData(Path_DT)) {
                    OrderOkay <- TRUE
                  } else {
                    OrderOkay <- FALSE
                  }
                } else {
                  OrderOkay <- FALSE
                }
              },
              error = function(e) {
                attempt <- attempt + 1
                if (attempt > max_attempts) {
                  stop("Failed after 3 attempts. Aborting.")
                }
              }
            )
          }

          fs::file_delete(Chunks)
          return(Path_DT)
        },
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
      )
    )

  invisible(gc())

  # # ..................................................................... ###

  # Number of observations and species per order ----
  IASDT.R::CatTime("Number of observations and species per order", Level = 1)

  Efforts_Summary <- Efforts_Summary %>%
    dplyr::mutate(
      NObs_NSp_R = furrr::future_map2(
        .x = Path_DT, .y = ClassOrder,
        .f = ~ {

          if (is.na(.x)) {
            ObsN <- 0
          } else {
            DT <- IASDT.R::LoadAs(.x)
            ObsN <- nrow(DT)
          }

          Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

          if (ObsN == 0) {

            NObs_R <- NSp_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_R <- setNames(NObs_R, paste0("NObs_", .y)) %>%
              terra::wrap()
            NSp_R <- setNames(NSp_R, paste0("NSp_", .y)) %>%
              terra::wrap()

          } else {

            # Number of observations
            NObs_R <- sf::st_drop_geometry(DT) %>%
              dplyr::count(CellCode, name = "NObs") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(y = Grid_R, field = "NObs") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NObs_", .y)) %>%
              terra::wrap()

            # Number of species
            NSp_R <- sf::st_drop_geometry(DT) %>%
              dplyr::distinct(CellCode, speciesKey) %>%
              dplyr::count(CellCode, name = "NSp") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(y = Grid_R, field = "NSp") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NSp_", .y)) %>%
              terra::wrap()
          }

          return(
            tibble::tibble(
              ObsN = ObsN, NObs_R = list(NObs_R), NSp_R = list(NSp_R)))

        },
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
      )
    ) %>%
    tidyr::unnest_wider("NObs_NSp_R") %>%
    tidyr::unnest(c("NObs_R", "NSp_R"))

  invisible(gc())

  # # ..................................................................... ###

  # Number of observations and species per order - native species ----
  IASDT.R::CatTime(
    "Number of observations and species per order - native species",
    Level = 1
  )

  Efforts_Summary <- Efforts_Summary %>%
    dplyr::mutate(
      Native = furrr::future_map2(
        .x = Path_DT, .y = ClassOrder,
        .f = ~ {
          if (is.na(.x)) {
            ObsN_Native <- 0
          } else {
            # Data on native species
            DT_Native <- IASDT.R::LoadAs(.x) %>%
              dplyr::filter(!(speciesKey %in% IAS_List))
            ObsN_Native <- nrow(DT_Native)
          }

          Grid_R <- terra::unwrap(IASDT.R::LoadAs(Path_Grid_R))

          if (ObsN_Native == 0) {
            NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_Native_R <- NSp_Native_R <- terra::classify(Grid_R, cbind(1, 0))
            NObs_Native_R <- setNames(
              NObs_Native_R, paste0("NObsNative_", .y)) %>%
              terra::wrap()
            NSp_Native_R <- setNames(NSp_Native_R, paste0("NSpNative_", .y)) %>%
              terra::wrap()
          } else {
            # Number of observations of native species
            NObs_Native_R <- sf::st_drop_geometry(DT_Native) %>%
              dplyr::count(CellCode, name = "NObs_Native") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(Grid_R, field = "NObs_Native") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NObsNative_", .y)) %>%
              terra::wrap()

            # # ................................... ###

            # Number of native species
            NSp_Native_R <- sf::st_drop_geometry(DT_Native) %>%
              dplyr::distinct(CellCode, speciesKey) %>%
              dplyr::count(CellCode, name = "NSp_Native") %>%
              dplyr::left_join(Grid_SF, by = "CellCode") %>%
              sf::st_as_sf() %>%
              terra::rasterize(Grid_R, field = "NSp_Native") %>%
              terra::classify(cbind(NA, 0)) %>%
              terra::mask(Grid_R) %>%
              IASDT.R::setRastCRS() %>%
              IASDT.R::setRastVals() %>%
              stats::setNames(paste0("NSpNative_", .y)) %>%
              terra::wrap()
          }

          tibble::tibble(
            ObsN_Native = ObsN_Native, NObs_Native_R = list(NObs_Native_R),
            NSp_Native_R = list(NSp_Native_R)) %>%
            return()
        },
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
      )
    ) %>%
    tidyr::unnest_wider("Native") %>%
    tidyr::unnest(c("NObs_Native_R", "NSp_Native_R"))

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
