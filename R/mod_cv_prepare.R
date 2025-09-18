## |------------------------------------------------------------------------| #
# mod_cv_prepare ----
## |------------------------------------------------------------------------| #

#' Prepare spatial-block cross-validation folds for spatial analysis
#'
#' This function assign modelling input data into spatial-block cross-validation
#' folds using three strategies (see below) using [blockCV::cv_spatial]. The
#' function is planned to be used inside the [IASDT.R::mod_prepare_hpc]
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
#' @param cv_n_grids Integer. Number of grid cells in both directions used in
#'   the `cv_dist` cross-validation strategy (see below). Default: 20L.
#' @param cv_n_folds Integer. Number of cross-validation folds. Default: 4L.
#' @param cv_n_rows,cv_n_columns Integer. Number of rows and columns used in the
#'   `cv_large` cross-validation strategy (see below), in which the study area
#'   is divided into large blocks given the provided `cv_n_rows` and
#'   `cv_n_columns` values. Both default to 2L which means to split the study
#'   area into four large blocks at the median latitude and longitude.
#' @param cv_sac Logical. Whether to use the spatial autocorrelation to
#'   determine the block size. Defaults to `FALSE`,
#' @param out_path Character. Path for directory to save the cross-validation
#'   results. This argument is mandatory and can not be empty.
#' @name mod_cv_prepare
#' @author Ahmed El-Gabbas
#' @return The function returns a modified version of the input dataset with
#'   additional numeric columns (integer) indicating the cross-validation
#'   strategy used.
#' @note The function uses the following cross-validation strategies:
#' - `cv_dist` in which the size of spatial cross-validation blocks is
#'   determined by the `cv_n_grids` argument. The default `cv_n_grids` value is
#'   20L, which means blocks of 20&times;20 grid cell each.
#' - `cv_large` which splits the study area into large blocks, as determined
#'   by the  `cv_n_rows` and `cv_n_columns` arguments. if `cv_n_rows =
#'   cv_n_columns` = 2L (default), four large blocks will be used, split the
#'   study area at the median coordinates.
#' - `cv_sac` in which the size of the blocks is determined by the median
#'   spatial autocorrelation range in the predictor data (estimated using
#'   [blockCV::cv_spatial_autocor]). This requires the availability of the
#'   `automap` R package. This strategy is currently skipped by default.
#' @export

mod_cv_prepare <- function(
    input_data = NULL, env_file = ".env", x_vars = NULL, cv_n_folds = 4L,
    cv_n_grids = 20L, cv_n_rows = 2L, cv_n_columns = 2L, cv_sac = FALSE,
    out_path = NULL) {

  # # |||||||||||||||||||||||||||||||||||
  # # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid <- eu_boundaries <- NULL

  if (is.null(input_data) || is.null(out_path) || is.null(x_vars)) {
    ecokit::stop_ctx(
      "`input_data`, `out_path`, and `x_vars` can not be empty",
      input_data = input_data, out_path = out_path,
      x_vars = x_vars, include_backtrace = TRUE)
  }

  all_vars <- c("x", "y", x_vars)
  all_vars_in_data <- all(all_vars %in% names(input_data))
  if (isFALSE(all_vars_in_data)) {
    missing_vars <- setdiff(all_vars, names(input_data))
    ecokit::stop_ctx(
      paste0(
        "input_data must contain 'x' and 'y' columns and all ",
        "environmental predictors in the x_vars argument."),
      all_vars_in_data = all_vars_in_data, missing_vars = missing_vars,
      include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Reference grid -----
  # # |||||||||||||||||||||||||||||||||||

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  path_grid <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(path_grid)) {
    ecokit::stop_ctx(
      "Path for reference grid does not exist", path_grid = path_grid,
      include_backtrace = TRUE)
  }
  ref_grid <- ecokit::load_as(path_grid, unwrap_r = TRUE)

  # # |||||||||||||||||||||||||||||||||||
  # # Coordinates as raster -----
  # # |||||||||||||||||||||||||||||||||||

  data_r <- dplyr::select(input_data, "x", "y") %>%
    as.matrix() %>%
    terra::rasterize(ref_grid) %>%
    terra::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # input data as sf -----
  # # |||||||||||||||||||||||||||||||||||

  xy_sf <- dplyr::select(input_data, "x", "y") %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  # # |||||||||||||||||||||||||||||||||||
  # # predictors as raster stack -----
  # # |||||||||||||||||||||||||||||||||||

  data_r <- dplyr::select(input_data, dplyr::all_of(all_vars)) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(y = ref_grid, field = names(.)[-length(.)]) %>%
    terra::trim()

  # # |||||||||||||||||||||||||||||||||||
  # # 1. CV using large blocks -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("1. cv_large", level = 1L, cat_timestamp = FALSE)
  cv_large <- blockCV::cv_spatial(
    x = xy_sf, r = data_r, hexagon = FALSE, iteration = 1000, k = cv_n_folds,
    rows_cols = c(cv_n_rows, cv_n_columns),
    plot = FALSE, progress = FALSE, report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 2. CV based on number of grid cells -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("2. cv_dist", level = 1L, cat_timestamp = FALSE)
  cv_dist <- blockCV::cv_spatial(
    x = xy_sf, r = data_r, hexagon = FALSE, iteration = 1000, k = cv_n_folds,
    size = cv_n_grids * raster::res(data_r)[1], plot = FALSE, progress = FALSE,
    report = FALSE)

  # # |||||||||||||||||||||||||||||||||||
  # # 3. CV based on spatial autocorrelation in predictor raster files -----
  # # |||||||||||||||||||||||||||||||||||

  if (cv_sac) {

    ecokit::cat_time("3. cv_sac", level = 1L)

    # Measure spatial autocorrelation in predictor raster files
    ecokit::cat_time("Calculating median SAC range", level = 2L)
    cv_sac_range <- blockCV::cv_spatial_autocor(
      r = data_r, num_sample = min(10000, nrow(input_data)), plot = FALSE,
      progress = FALSE)

    # If median spatial cross-validation is very large, skip cv_sac option. This
    # to avoid the following error if the estimated SAC range is very large that
    # disallow the assignation of grid cells into the number of CV folds in the
    # `cv_n_folds` parameter:
    # error in `blockCV::cv_spatial()`: 'k' is bigger than the number of spatial
    # blocks: 1.

    # Size of the large spatial blocks that split the study area into four large
    # blocks; the same as `cv_large` method
    min_distance <- min(
      (terra::nrow(ref_grid) * terra::res(ref_grid)[1] / 2),
      (terra::ncol(ref_grid) * terra::res(ref_grid)[1] / 2))

    if (cv_sac_range$range > min_distance) {
      ecokit::cat_time(
        paste0(
          "`cv_sac` was NOT implemented; median SAC: ",
          round(cv_sac_range$range / 1000, 2), "km"), level = 2L)
      cv_sac <- NULL

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!("folds_ids" %in% names(cv_dist) &&
            "folds_ids" %in% names(cv_large))) {
        ecokit::stop_ctx(
          "Cross-validation results do not contain 'folds_ids'.",
          names_cv_Large = names(cv_large), names_cv_Dist = names(cv_dist),
          include_backtrace = TRUE)
      }
    } else {
      # CV based on Spatial autocorrelation
      cv_sac <- blockCV::cv_spatial(
        x = xy_sf, r = data_r, hexagon = FALSE, iteration = 1000,
        k = cv_n_folds, size = cv_sac_range$range, plot = FALSE,
        progress = FALSE, report = FALSE)

      # Check `folds_ids` exists in each of the cross-validation strategies
      if (!(("folds_ids" %in% names(cv_sac)) &&
            ("folds_ids" %in% names(cv_dist)) &&
            ("folds_ids" %in% names(cv_large)))) {
        ecokit::stop_ctx(
          "Cross-validation results do not contain 'folds_ids'.",
          names_cv_Large = names(cv_large), names_cv_Dist = names(cv_dist),
          names_cv_SAC = names(cv_sac), include_backtrace = TRUE)
      }
    }
  } else {
    cv_sac <- cv_sac_range <- NULL
  }

  rm(xy_sf, envir = environment())

  # # |||||||||||||||||||||||||||||||||||
  # # Save cross-validation results as RData -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save cross-validation results as RData", level = 1L)
  cv_data <- list(
    cv_n_grids = cv_n_grids, cv_n_rows = cv_n_rows, cv_n_columns = cv_n_columns,
    cv_sac_range = cv_sac_range, cv_sac = cv_sac, cv_dist = cv_dist,
    cv_large = cv_large)

  save(cv_data, file = fs::path(out_path, "cv_data.RData"))

  # # |||||||||||||||||||||||||||||||||||
  # # Plot cross-validation folds -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Plot cross-validation folds", level = 1L)

  data_r <- sf::st_as_sf(input_data, coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(ref_grid) %>%
    terra::classify(cbind(1, 0)) %>%
    terra::as.factor() %>%
    stats::setNames("GridR")

  if (is.null(cv_sac)) {
    cv_types <- c("cv_dist", "cv_large")
  } else {
    cv_types <- c("cv_sac", "cv_dist", "cv_large")
  }

  eu_boundaries <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03") %>%
    suppressWarnings()

  cv_plots <- purrr::map(
    .x = cv_types,
    .f = ~{

      blocks <- magrittr::extract2(cv_data, .x) %>%
        magrittr::extract2("blocks") %>%
        dplyr::mutate(folds = factor(folds))

      plot <- ggplot2::ggplot(environment = emptyenv()) +
        ggplot2::geom_sf(
          data = eu_boundaries, fill = "gray95", colour = "darkgrey",
          linewidth = 0.5) +
        tidyterra::geom_spatraster(data = data_r, inherit.aes = FALSE) +
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

      return(plot)
    })

  grDevices::cairo_pdf(
    filename = fs::path(out_path, "cv_blocks.pdf"),
    width = 12.5, height = 13, onefile = TRUE)
  invisible(purrr::map(cv_plots, print))
  grDevices::dev.off()

  # # |||||||||||||||||||||||||||||||||||
  # # Add cross-validation columns to the data -----
  # # |||||||||||||||||||||||||||||||||||

  if (!is.null(cv_sac)) {
    input_data$cv_sac <- magrittr::extract2(cv_sac, "folds_ids")
  }
  input_data$cv_dist <- magrittr::extract2(cv_dist, "folds_ids")
  input_data$cv_large <- magrittr::extract2(cv_large, "folds_ids")

  return(input_data)
}
