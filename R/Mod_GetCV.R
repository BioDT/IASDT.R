## |------------------------------------------------------------------------| #
# GetCV ----
## |------------------------------------------------------------------------| #

#' Prepare spatial-block cross-validation folds for spatial analysis
#'
#' This function assign modelling input data into spatial-block cross-validation
#' folds. The function uses [blockCV::cv_spatial] to assign grid cells to
#' cross-validation folds using three strategies (see below). The function is
#' planned to be used inside the [IASDT.R::Mod_Prep4HPC] function.
#' @param DT A data frame or tibble containing the input dataset. This data
#'   frame should include two columns for `x` and `y` coordinates as long as
#'   other columns matching the names of predictors listed in `XVars` argument.
#' @param EnvFile String specifying the path to read environment variables from,
#'   with a default value of `.env`.
#' @param XVars Vector of strings specifying variables to be used in the model.
#'   This argument is mandatory and can not be empty.
#' @param NGrids For `CV_Dist` cross-validation strategy (see below), this
#'   argument determines the size of the blocks (how many grid cells in both
#'   directions).
#' @param NFolds Number of cross-validation folds. Default: 4.
#' @param NR,NC Integer, the number of rows and columns used in the `CV_Large`
#'   cross-validation strategy (see below), in which the study area is divided
#'   into large blocks given the provided `NR` and `NC` values. Both default to
#'   2 which means to split the study area into four large blocks at the median
#'   latitude and longitude.
#' @param OutPath String specifying the folder path to save the cross-validation
#'   results. Default: `NULL`.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param PlotCV Logical. Indicating whether to plot the block cross-validation
#'   folds.
#' @name GetCV
#' @author Ahmed El-Gabbas
#' @return The function returns a modified version of the input dataset `DT`
#'   with 3 additional numeric columns (integer) indicating the cross-validation
#'   folds:
#'
#'   1) `CV_SAC` in which the size of the blocks is determined by the median
#'   spatial autocorrelation range in the predictor data (estimated using
#'   [blockCV::cv_spatial_autocor]). This requires the availability of the
#'   `automap` R package.
#'
#'   2) `CV_Dist` in which the size of spatial cross-validation blocks is
#'   determined by the `NGrids` argument. The default `NGrids` value is 20,
#'   which means blocks of 20x20 grid cell each.
#'
#'   3) `CV_Large` which splits the study area into large blocks, as determined
#'   by the  `NR` and `NC` arguments. if `NR = NC` = 2 (default), four large
#'   blocks will be used, split the study area at the median coordinates..
#' @details The function reads the following environment variable:
#'    - **`DP_R_Grid`** (if `FromHPC = TRUE`) or
#'    **`DP_R_Grid_Local`** (if `FromHPC = FALSE`). The function reads
#'   the content of the `Grid_10_Land_Crop.RData` file from this path.
#'   - **`DP_R_EUBound_sf`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_EUBound_sf_Local`** (if `FromHPC` = `FALSE`): path for the
#'     `RData` file containing the country boundaries (`sf` object).
#' @export

GetCV <- function(
    DT, EnvFile = ".env", XVars, NFolds = 4, NGrids = 20,
    NR = 2, NC = 2, OutPath = NULL, FromHPC = TRUE, PlotCV = TRUE) {

  # # |||||||||||||||||||||||||||||||||||
  # # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- EU_Bound <- NULL

  if (is.null(DT) || is.null(EnvFile) || is.null(OutPath) || is.null(XVars)) {
    stop("DT, EnvFile, OutPath, and XVars can not be empty")
  }

  if (magrittr::not(file.exists(EnvFile))) {
    stop(paste0("Path for environment variables: ", EnvFile, " was not found"))
  }

  AllVars <- c("x", "y", XVars)
  AllVarsInDT <- all(AllVars %in% names(DT))
  if (magrittr::not(AllVarsInDT)) {
    MissingVars <- setdiff(AllVars, names(DT))
    stop(
      paste0(
        "Data frame DT must contain 'x' and 'y' columns and all ",
        "environmental predictors in the XVars argument.\nMissing vars are ",
        paste0(MissingVars, collapse = "; ")),
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Reference grid -----
  # # |||||||||||||||||||||||||||||||||||

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  Path_Grid <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (magrittr::not(file.exists(Path_Grid))) {
    stop("Path for reference grid does not exist")
  }
  RefGrid <- terra::unwrap(IASDT.R::LoadAs(Path_Grid))


  # # |||||||||||||||||||||||||||||||||||
  # # Coordinates as raster -----
  # # |||||||||||||||||||||||||||||||||||
  DT_R <- dplyr::select(DT, "x", "y") %>%
    as.matrix() %>%
    terra::rasterize(RefGrid) %>%
    terra::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # input data as sf -----
  # # |||||||||||||||||||||||||||||||||||
  XY_sf <- dplyr::select(DT, "x", "y") %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  # # |||||||||||||||||||||||||||||||||||
  # # data as raster stack -----
  # # |||||||||||||||||||||||||||||||||||
  DT_R <- dplyr::select(DT, dplyr::all_of(AllVars)) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(y = RefGrid, field = names(.)[-length(.)]) %>%
    raster::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # 1. CV using large blocks -----
  # # |||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("   >>>   1. CV_Large")
  CV_Large <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = NFolds,
    rows_cols = c(NR, NC), plot = FALSE, progress = FALSE, report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 2. CV based on number of grid cells -----
  # # |||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("   >>>   2. CV_Dist")
  CV_Dist <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = NFolds,
    size = NGrids * raster::res(DT_R)[1], plot = FALSE, progress = FALSE,
    report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 3. CV based on spatial autocorrelation in predictor raster files -----
  # # |||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("   >>>   3. CV_SAC")

  # Measure spatial autocorrelation in predictor raster files
  IASDT.R::CatTime("   >>>   >>>  Calculating median SAC range")
  CV_SAC_Range <- blockCV::cv_spatial_autocor(
    r = DT_R, num_sample = min(10000, nrow(DT)), plot = FALSE,
    progress = FALSE)

  # If median spatial cross-validation is very large, skip CV_SAC option. This
  # to avoid the following error if the estimated SAC range is very large that
  # disallow the assignation of grid cells into the number of CV folds in the
  # `NFolds` parameter:
  # ! 'k' is bigger than the number of spatial blocks: 1.

  # Size of the large spatial blocks that split the study area into four large
  # blocks; the same as `CV_Large` method
  MinDist <- min(
    (terra::nrow(RefGrid) * terra::res(RefGrid)[1] / 2),
    (terra::ncol(RefGrid) * terra::res(RefGrid)[1] / 2))

  if (CV_SAC_Range$range > MinDist) {
    IASDT.R::CatTime(
      paste0("   >>>   >>>  CV_SAC was NOT implemented; median SAC: ",
             CV_SAC_Range$range, "m"))
    CV_SAC <- NULL

    # Check `folds_ids` exists in each of the cross-validation strategies
    if (magrittr::not(
      "folds_ids" %in% names(CV_Dist) && "folds_ids" %in% names(CV_Large))) {
      stop(
        "Cross-validation results do not contain 'folds_ids'.",
        call. = FALSE)
    }
  } else {
    # CV based on Spatial autocorrelation
    CV_SAC <- blockCV::cv_spatial(
      x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = NFolds,
      size = CV_SAC_Range$range, plot = FALSE, progress = FALSE,
      report = FALSE)

    # Check `folds_ids` exists in each of the cross-validation strategies
    if (magrittr::not(
      "folds_ids" %in% names(CV_SAC) && "folds_ids" %in% names(CV_Dist) &&
      "folds_ids" %in% names(CV_Large))) {
      stop(
        "Cross-validation results do not contain 'folds_ids'.",
        call. = FALSE)
    }
  }

  rm(XY_sf)

  # # |||||||||||||||||||||||||||||||||||
  # # Save cross-validation results as RData -----
  # # |||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("   >>>   Save cross-validation results as RData")
  CV_data <- list(
    NGrids = NGrids, NR = NR, NC = NC, CV_SAC_Range = CV_SAC_Range,
    CV_SAC = CV_SAC, CV_Dist = CV_Dist, CV_Large = CV_Large)

  save(CV_data, file = file.path(OutPath, "CV_data.RData"))

  # # |||||||||||||||||||||||||||||||||||
  # # Plot cross-validation folds -----
  # # |||||||||||||||||||||||||||||||||||
  if (PlotCV) {

    IASDT.R::CatTime("   >>>   Plot cross-validation folds")

    DT_R <- sf::st_as_sf(DT, coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(RefGrid) %>%
      terra::trim() %>%
      terra::classify(cbind(1, 0)) %>%
      terra::as.factor() %>%
      stats::setNames("GridR")

    if (is.null(CV_SAC)) {
      PlotBox <- rbind(CV_Dist$blocks, CV_Large$blocks) %>%
        sf::st_bbox()
      CVTypes <- c("CV_Dist", "CV_Large")
    } else {
      PlotBox <- rbind(CV_SAC$blocks, CV_Dist$blocks, CV_Large$blocks) %>%
        sf::st_bbox()
      CVTypes <- c("CV_SAC", "CV_Dist", "CV_Large")
    }

    AspectRatio <- (PlotBox[3] - PlotBox[1]) / (PlotBox[4] - PlotBox[2])

    EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur_s") %>%
      magrittr::extract2("L_01") %>%
      sf::st_crop(PlotBox) %>%
      suppressWarnings()

    CV_Plots <- purrr::map(
      .x = CVTypes,
      .f = ~{

        blocks <- magrittr::extract2(CV_data, .x) %>%
          magrittr::extract2("blocks") %>%
          dplyr::mutate(folds = factor(folds))

        Plot <- ggplot2::ggplot() +
          tidyterra::geom_spatraster(data = DT_R) +
          ggplot2::geom_sf(
            data = EU_Bound, fill = "transparent", colour = "darkgrey",
            linewidth = 1.5) +
          ggplot2::geom_sf(
            data = blocks, inherit.aes = FALSE, alpha = 0.5, ,
            mapping = ggplot2::aes(fill = folds), linewidth = 0.4) +
          ggplot2::scale_fill_manual(
            values = c(
              "darkgrey", "transparent", "red", "green", "blue", "yellow"),
            na.value = "transparent") +
          ggplot2::scale_x_continuous(
            limits = PlotBox[c(1, 3)], expand = c(0, 0)) +
          ggplot2::scale_y_continuous(
            limits = PlotBox[c(2, 4)], expand = c(0, 0)) +
          ggplot2::geom_sf_text(
            data = blocks, ggplot2::aes(label = folds), size = 16,
            fontface = "bold") +
          ggplot2::ggtitle(
            paste0("<span style='font-size: 40pt;'>", .x, "</font>")) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.title = ggtext::element_markdown(), legend.position = "none")

        return(Plot)
      })

    ggplot2::ggsave(
      filename = file.path(OutPath, "CV_Blocks.pdf"),
      plot = gridExtra::marrangeGrob(
        CV_Plots, nrow = 1, ncol = 1, top = NULL),
      width = 25 * (AspectRatio), height = 25, unit = "in")
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Add cross-validation columns to the data -----
  # # |||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("   >>>   Add cross-validation columns to the data")

  if (is.null(CV_SAC)) {
    DT$CV_SAC <- NA_integer_
  } else {
    DT$CV_SAC <- magrittr::extract2(CV_SAC, "folds_ids")
  }
  DT$CV_Dist <- magrittr::extract2(CV_Dist, "folds_ids")
  DT$CV_Large <- magrittr::extract2(CV_Large, "folds_ids")

  return(DT)
}
