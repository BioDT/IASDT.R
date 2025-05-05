## |------------------------------------------------------------------------| #
# mod_CV_prepare ----
## |------------------------------------------------------------------------| #

#' Prepare spatial-block cross-validation folds for spatial analysis
#'
#' This function assign modelling input data into spatial-block cross-validation
#' folds using three strategies (see below) using [blockCV::cv_spatial]. The
#' function is planned to be used inside the [IASDT.R::mod_prepare_HPC]
#' function.
#' @param input_data `data.frame`. A data frame or tibble containing the input
#'   dataset. This data frame should include two columns for `x` and `y`
#'   coordinates as long as other columns matching the names of predictors
#'   listed in `x_vars` argument. This argument is mandatory and can not be
#'   empty.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param x_vars Character vector. Variables to be used in the model. This
#'   argument is mandatory and can not be empty.
#' @param CV_n_grids Integer. Number of grid cells in both directions used in
#'   the `CV_Dist` cross-validation strategy (see below). Default: 20L.
#' @param CV_n_folds Integer. Number of cross-validation folds. Default: 4L.
#' @param CV_n_rows,CV_n_columns Integer. Number of rows and columns used in the
#'   `CV_Large` cross-validation strategy (see below), in which the study area
#'   is divided into large blocks given the provided `CV_n_rows` and
#'   `CV_n_columns` values. Both default to 2L which means to split the study
#'   area into four large blocks at the median latitude and longitude.
#' @param CV_SAC Logical. Whether to use the spatial autocorrelation to
#'   determine the block size. Defaults to `FALSE`,
#' @param out_path Character. Path for directory to save the cross-validation
#'   results. This argument is mandatory and can not be empty.
#' @param CV_plot Logical. Indicating whether to plot the block cross-validation
#'   folds.
#' @name mod_CV_prepare
#' @author Ahmed El-Gabbas
#' @return The function returns a modified version of the input dataset with
#'   additional numeric columns (integer) indicating the cross-validation
#'   strategy used.
#' @note The function uses the following cross-validation strategies:
#' - `CV_Dist` in which the size of spatial cross-validation blocks is
#'   determined by the `CV_n_grids` argument. The default `CV_n_grids` value is
#'   20L, which means blocks of 20&times;20 grid cell each.
#' - `CV_Large` which splits the study area into large blocks, as determined
#'   by the  `CV_n_rows` and `CV_n_columns` arguments. if `CV_n_rows =
#'   CV_n_columns` = 2L (default), four large blocks will be used, split the
#'   study area at the median coordinates.
#' - `CV_SAC` in which the size of the blocks is determined by the median
#'   spatial autocorrelation range in the predictor data (estimated using
#'   [blockCV::cv_spatial_autocor]). This requires the availability of the
#'   `automap` R package. This strategy is currently skipped by default.
#' @export

mod_CV_prepare <- function(
    input_data = NULL, env_file = ".env", x_vars = NULL, CV_n_folds = 4L,
    CV_n_grids = 20L, CV_n_rows = 2, CV_n_columns = 2L, CV_SAC = FALSE,
    out_path = NULL, CV_plot = TRUE) {

  # # |||||||||||||||||||||||||||||||||||
  # # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid <- EU_Bound <- NULL

  if (is.null(input_data) || is.null(env_file) ||
      is.null(out_path) || is.null(x_vars)) {
    ecokit::stop_ctx(
      "`input_data`, `env_file`, `out_path`, and `x_vars` can not be empty",
      input_data = input_data, env_file = env_file, out_path = out_path,
      x_vars = x_vars)
  }

  if (!file.exists(env_file)) {
    ecokit::stop_ctx(
      "Path to environment variables does not exist ", env_file = env_file)
  }

  AllVars <- c("x", "y", x_vars)
  AllVarsInDT <- all(AllVars %in% names(input_data))
  if (isFALSE(AllVarsInDT)) {
    MissingVars <- setdiff(AllVars, names(input_data))
    ecokit::stop_ctx(
      paste0(
        "input_data must contain 'x' and 'y' columns and all ",
        "environmental predictors in the x_vars argument."),
      AllVarsInDT = AllVarsInDT, MissingVars = MissingVars)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Reference grid -----
  # # |||||||||||||||||||||||||||||||||||

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Path_Grid <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_Grid)) {
    ecokit::stop_ctx(
      "Path for reference grid does not exist", Path_Grid = Path_Grid)
  }
  RefGrid <- terra::unwrap(ecokit::load_as(Path_Grid))

  # # |||||||||||||||||||||||||||||||||||
  # # Coordinates as raster -----
  # # |||||||||||||||||||||||||||||||||||

  DT_R <- dplyr::select(input_data, "x", "y") %>%
    as.matrix() %>%
    terra::rasterize(RefGrid) %>%
    terra::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # input data as sf -----
  # # |||||||||||||||||||||||||||||||||||

  XY_sf <- dplyr::select(input_data, "x", "y") %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  # # |||||||||||||||||||||||||||||||||||
  # # data as raster stack -----
  # # |||||||||||||||||||||||||||||||||||

  DT_R <- dplyr::select(input_data, dplyr::all_of(AllVars)) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(y = RefGrid, field = names(.)[-length(.)]) %>%
    raster::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # 1. CV using large blocks -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("1. CV_Large", level = 1L, cat_timestamp = FALSE)
  CV_Large <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_n_folds,
    rows_cols = c(CV_n_rows, CV_n_columns),
    plot = FALSE, progress = FALSE, report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 2. CV based on number of grid cells -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("2. CV_Dist", level = 1L, cat_timestamp = FALSE)
  CV_Dist <- blockCV::cv_spatial(
    x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_n_folds,
    size = CV_n_grids * raster::res(DT_R)[1], plot = FALSE, progress = FALSE,
    report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 3. CV based on spatial autocorrelation in predictor raster files -----
  # # |||||||||||||||||||||||||||||||||||

  if (CV_SAC) {

    ecokit::cat_time("3. CV_SAC", level = 1L)

    # Measure spatial autocorrelation in predictor raster files
    ecokit::cat_time("Calculating median SAC range", level = 2L)
    CV_SAC_Range <- blockCV::cv_spatial_autocor(
      r = DT_R, num_sample = min(10000, nrow(input_data)), plot = FALSE,
      progress = FALSE)

    # If median spatial cross-validation is very large, skip CV_SAC option. This
    # to avoid the following error if the estimated SAC range is very large that
    # disallow the assignation of grid cells into the number of CV folds in the
    # `CV_n_folds` parameter:
    # error in `blockCV::cv_spatial()`: 'k' is bigger than the number of spatial
    # blocks: 1.

    # Size of the large spatial blocks that split the study area into four large
    # blocks; the same as `CV_Large` method
    min_distance <- min(
      (terra::nrow(RefGrid) * terra::res(RefGrid)[1] / 2),
      (terra::ncol(RefGrid) * terra::res(RefGrid)[1] / 2))

    if (CV_SAC_Range$range > min_distance) {
      ecokit::cat_time(
        paste0(
          "`CV_SAC` was NOT implemented; median SAC: ",
          round(CV_SAC_Range$range / 1000, 2), "km"), level = 2L)
      CV_SAC <- NULL

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!("folds_ids" %in% names(CV_Dist) &&
            "folds_ids" %in% names(CV_Large))) {
        ecokit::stop_ctx(
          "Cross-validation results do not contain 'folds_ids'.",
          names_CV_Large = names(CV_Large), names_CV_Dist = names(CV_Dist))
      }
    } else {
      # CV based on Spatial autocorrelation
      CV_SAC <- blockCV::cv_spatial(
        x = XY_sf, r = DT_R, hexagon = FALSE, iteration = 1000, k = CV_n_folds,
        size = CV_SAC_Range$range, plot = FALSE, progress = FALSE,
        report = FALSE)

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!(("folds_ids" %in% names(CV_SAC)) &&
            ("folds_ids" %in% names(CV_Dist)) &&
            ("folds_ids" %in% names(CV_Large)))) {
        ecokit::stop_ctx(
          "Cross-validation results do not contain 'folds_ids'.",
          names_CV_Large = names(CV_Large), names_CV_Dist = names(CV_Dist),
          names_CV_SAC = names(CV_SAC))
      }
    }
  } else {
    CV_SAC <- CV_SAC_Range <- NULL
  }

  rm(XY_sf, envir = environment())

  # # |||||||||||||||||||||||||||||||||||
  # # Save cross-validation results as RData -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save cross-validation results as RData", level = 1L)
  CV_data <- list(
    CV_n_grids = CV_n_grids, CV_n_rows = CV_n_rows, CV_n_columns = CV_n_columns,
    CV_SAC_Range = CV_SAC_Range, CV_SAC = CV_SAC, CV_Dist = CV_Dist,
    CV_Large = CV_Large)

  save(CV_data, file = fs::path(out_path, "CV_data.RData"))

  # # |||||||||||||||||||||||||||||||||||
  # # Plot cross-validation folds -----
  # # |||||||||||||||||||||||||||||||||||

  if (CV_plot) {

    ecokit::cat_time("Plot cross-validation folds", level = 1L)

    DT_R <- sf::st_as_sf(input_data, coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(RefGrid) %>%
      terra::classify(cbind(1, 0)) %>%
      terra::as.factor() %>%
      stats::setNames("GridR")

    if (is.null(CV_SAC)) {
      CVTypes <- c("CV_Dist", "CV_Large")
    } else {
      CVTypes <- c("CV_SAC", "CV_Dist", "CV_Large")
    }

    EU_Bound <- ecokit::load_as(EU_Bound) %>%
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
      filename = fs::path(out_path, "CV_Blocks.pdf"),
      width = 12.5, height = 13, onefile = TRUE)
    invisible(purrr::map(CV_Plots, print))
    grDevices::dev.off()
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Add cross-validation columns to the data -----
  # # |||||||||||||||||||||||||||||||||||

  if (!is.null(CV_SAC)) {
    input_data$CV_SAC <- magrittr::extract2(CV_SAC, "folds_ids")
  }
  input_data$CV_Dist <- magrittr::extract2(CV_Dist, "folds_ids")
  input_data$CV_Large <- magrittr::extract2(CV_Large, "folds_ids")

  return(input_data)
}
