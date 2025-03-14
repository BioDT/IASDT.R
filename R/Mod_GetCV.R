## |------------------------------------------------------------------------| #
# Mod_GetCV ----
## |------------------------------------------------------------------------| #

#' Prepare spatial-block cross-validation folds for spatial analysis
#'
#' This function assign modelling input data into spatial-block cross-validation
#' folds using three strategies (see below) using [blockCV::cv_spatial]. The
#' function is planned to be used inside the [IASDT.R::Mod_Prep4HPC] function.
#' @param Data `data.frame`. A data frame or tibble containing the input
#'   dataset. This data frame should include two columns for `x` and `y`
#'   coordinates as long as other columns matching the names of predictors
#'   listed in `XVars` argument. This argument is mandatory and can not be
#'   empty.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param XVars Character vector. Variables to be used in the model.
#'   This argument is mandatory and can not be empty.
#' @param CV_NGrids Integer. Number of grid cells in both directions used in
#'   the `CV_Dist` cross-validation strategy (see below). Default: 20L.
#' @param CV_NFolds Integer. Number of cross-validation folds. Default: 4L.
#' @param CV_NR,CV_NC Integer. Number of rows and columns used in the `CV_Large`
#'   cross-validation strategy (see below), in which the study area is divided
#'   into large blocks given the provided `CV_NR` and `CV_NC` values. Both
#'   default to 2L which means to split the study area into four large blocks at
#'   the median latitude and longitude.
#' @param CV_SAC Logical. Whether to use the spatial autocorrelation
#'   to determine the block size. Defaults to `FALSE`,
#' @param OutPath Character. Path for directory to save the
#'   cross-validation results. This argument is mandatory and can not be empty.
#' @param CV_Plot Logical. Indicating whether to plot the block cross-validation
#'   folds.
#' @name Mod_GetCV
#' @author Ahmed El-Gabbas
#' @return The function returns a modified version of the input dataset with
#'   additional numeric columns (integer) indicating the cross-validation
#'   strategy used.
#' @note The function uses the following cross-validation strategies:
#' - `CV_Dist` in which the size of spatial cross-validation blocks is
#'   determined by the `CV_NGrids` argument. The default `CV_NGrids` value is
#'   20L, which means blocks of 20&times;20 grid cell each.
#' - `CV_Large` which splits the study area into large blocks, as determined
#'   by the  `CV_NR` and `CV_NC` arguments. if `CV_NR = CV_NC` = 2L (default),
#'   four large blocks will be used, split the study area at the median
#'   coordinates.
#' - `CV_SAC` in which the size of the blocks is determined by the median
#'   spatial autocorrelation range in the predictor data (estimated using
#'   [blockCV::cv_spatial_autocor]). This requires the availability of the
#'   `automap` R package. This strategy is currently skipped by default.

#' @export

Mod_GetCV <- function(
    Data = NULL, EnvFile = ".env", XVars = NULL, CV_NFolds = 4L,
    CV_NGrids = 20L, CV_NR = 2, CV_NC = 2L, CV_SAC = FALSE, OutPath = NULL,
    CV_Plot = TRUE) {

  # # |||||||||||||||||||||||||||||||||||
  # # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- EU_Bound <- NULL

  if (is.null(Data) || is.null(EnvFile) || is.null(OutPath) || is.null(XVars)) {
    stop("Data, EnvFile, OutPath, and XVars can not be empty", call. = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(
      "Path to environment variables: ", EnvFile, " was not found",
      call. = FALSE)
  }

  AllVars <- c("x", "y", XVars)
  AllVarsInDT <- all(AllVars %in% names(Data))
  if (isFALSE(AllVarsInDT)) {
    MissingVars <- setdiff(AllVars, names(Data))
    stop(
      "Data frame Data must contain 'x' and 'y' columns and all ",
      "environmental predictors in the XVars argument.\nMissing vars are ",
      paste(MissingVars, collapse = "; "), call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Reference grid -----
  # # |||||||||||||||||||||||||||||||||||

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Path_Grid <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_Grid)) {
    stop("Path for reference grid does not exist", call. = FALSE)
  }
  RefGrid <- terra::unwrap(IASDT.R::LoadAs(Path_Grid))

  # # |||||||||||||||||||||||||||||||||||
  # # Coordinates as raster -----
  # # |||||||||||||||||||||||||||||||||||

  DT_R <- dplyr::select(Data, "x", "y") %>%
    as.matrix() %>%
    terra::rasterize(RefGrid) %>%
    terra::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # input data as sf -----
  # # |||||||||||||||||||||||||||||||||||

  XY_sf <- dplyr::select(Data, "x", "y") %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  # # |||||||||||||||||||||||||||||||||||
  # # data as raster stack -----
  # # |||||||||||||||||||||||||||||||||||

  DT_R <- dplyr::select(Data, dplyr::all_of(AllVars)) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(y = RefGrid, field = names(.)[-length(.)]) %>%
    raster::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # 1. CV using large blocks -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("1. CV_Large", Level = 1, Time = FALSE)
  CV_Large <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_NFolds,
    rows_cols = c(CV_NR, CV_NC), plot = FALSE, progress = FALSE, report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 2. CV based on number of grid cells -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("2. CV_Dist", Level = 1, Time = FALSE)
  CV_Dist <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_NFolds,
    size = CV_NGrids * raster::res(DT_R)[1], plot = FALSE, progress = FALSE,
    report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 3. CV based on spatial autocorrelation in predictor raster files -----
  # # |||||||||||||||||||||||||||||||||||

  if (CV_SAC) {

    IASDT.R::CatTime("3. CV_SAC", Level = 1)

    # Measure spatial autocorrelation in predictor raster files
    IASDT.R::CatTime("Calculating median SAC range", Level = 2)
    CV_SAC_Range <- blockCV::cv_spatial_autocor(
      r = DT_R, num_sample = min(10000, nrow(Data)), plot = FALSE,
      progress = FALSE)

    # If median spatial cross-validation is very large, skip CV_SAC option. This
    # to avoid the following error if the estimated SAC range is very large that
    # disallow the assignation of grid cells into the number of CV folds in the
    # `CV_NFolds` parameter:
    # error in `blockCV::cv_spatial()`: 'k' is bigger than the number of spatial
    # blocks: 1.

    # Size of the large spatial blocks that split the study area into four large
    # blocks; the same as `CV_Large` method
    MinDist <- min(
      (terra::nrow(RefGrid) * terra::res(RefGrid)[1] / 2),
      (terra::ncol(RefGrid) * terra::res(RefGrid)[1] / 2))

    if (CV_SAC_Range$range > MinDist) {
      IASDT.R::CatTime(
        paste0(
          "`CV_SAC` was NOT implemented; median SAC: ",
          round(CV_SAC_Range$range / 1000, 2), "km"), Level = 2)
      CV_SAC <- NULL

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!("folds_ids" %in% names(CV_Dist) &&
            "folds_ids" %in% names(CV_Large))) {
        stop(
          "Cross-validation results do not contain 'folds_ids'.",
          call. = FALSE)
      }
    } else {
      # CV based on Spatial autocorrelation
      CV_SAC <- blockCV::cv_spatial(
        x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_NFolds,
        size = CV_SAC_Range$range, plot = FALSE, progress = FALSE,
        report = FALSE)

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!(("folds_ids" %in% names(CV_SAC)) &&
            ("folds_ids" %in% names(CV_Dist)) &&
            ("folds_ids" %in% names(CV_Large)))) {
        stop(
          "Cross-validation results do not contain 'folds_ids'.",
          call. = FALSE)
      }
    }
  } else {
    CV_SAC <- CV_SAC_Range <- NULL
  }

  rm(XY_sf, envir = environment())

  # # |||||||||||||||||||||||||||||||||||
  # # Save cross-validation results as RData -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save cross-validation results as RData", Level = 1)
  CV_data <- list(
    CV_NGrids = CV_NGrids, CV_NR = CV_NR, CV_NC = CV_NC,
    CV_SAC_Range = CV_SAC_Range, CV_SAC = CV_SAC, CV_Dist = CV_Dist,
    CV_Large = CV_Large)

  save(CV_data, file = IASDT.R::Path(OutPath, "CV_data.RData"))

  # # |||||||||||||||||||||||||||||||||||
  # # Plot cross-validation folds -----
  # # |||||||||||||||||||||||||||||||||||

  if (CV_Plot) {

    IASDT.R::CatTime("Plot cross-validation folds", Level = 1)

    DT_R <- sf::st_as_sf(Data, coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(RefGrid) %>%
      terra::classify(cbind(1, 0)) %>%
      terra::as.factor() %>%
      stats::setNames("GridR")

    if (is.null(CV_SAC)) {
      CVTypes <- c("CV_Dist", "CV_Large")
    } else {
      CVTypes <- c("CV_SAC", "CV_Dist", "CV_Large")
    }

    EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur_s") %>%
      magrittr::extract2("L_03") %>%
      suppressWarnings()

    CV_Plots <- purrr::map(
      .x = CVTypes,
      .f = ~{

        blocks <- magrittr::extract2(CV_data, .x) %>%
          magrittr::extract2("blocks") %>%
          dplyr::mutate(folds = factor(folds))

        Plot <- ggplot2::ggplot() +
          ggplot2::geom_sf(
            data = EU_Bound, fill = "gray95", colour = "darkgrey",
            linewidth = 0.5) +
          tidyterra::geom_spatraster(data = DT_R, inherit.aes = FALSE) +
          ggplot2::geom_sf(
            data = blocks, inherit.aes = FALSE, alpha = 0.35,
            mapping = ggplot2::aes(fill = folds), linewidth = 0.3) +
          ggplot2::geom_sf_text(
            data = blocks, ggplot2::aes(label = folds), size = 8,
            fontface = "bold") +
          ggplot2::coord_sf(
            xlim = c(2600000, 6550000), ylim = c(1450000, 5420000),
            clip = "off") +
          ggplot2::scale_fill_manual(
            values = c(
              "darkgrey", "transparent", "red", "green", "blue", "yellow"),
            na.value = "transparent") +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::labs(title = .x) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.margin = ggplot2::margin(0.25, 0.25, 0.125, 0.25, "cm"),
            plot.title = ggplot2::element_text(
              size = 20, face = "bold", color = "blue", hjust = 0.5,
              margin = ggplot2::margin(0.2, 0, 0.2, 0.5, "cm")),
            legend.position = "none")

        return(Plot)
      })

    grDevices::cairo_pdf(
      filename = IASDT.R::Path(OutPath, "CV_Blocks.pdf"),
      width = 12.5, height = 13, onefile = TRUE)
    invisible(purrr::map(CV_Plots, print))
    grDevices::dev.off()
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Add cross-validation columns to the data -----
  # # |||||||||||||||||||||||||||||||||||

  if (!is.null(CV_SAC)) {
    Data$CV_SAC <- magrittr::extract2(CV_SAC, "folds_ids")
  }
  Data$CV_Dist <- magrittr::extract2(CV_Dist, "folds_ids")
  Data$CV_Large <- magrittr::extract2(CV_Large, "folds_ids")

  return(Data)
}
