## |------------------------------------------------------------------------| #
# cv_merge_predictions ----
## |------------------------------------------------------------------------| #

#' Merge Cross-Validation Predictions and Generate Summary Maps
#'
#' Merges prediction results from cross-validated Hmsc models, processes clamped
#' and not-clamped model outputs, and generates summary maps and anomaly
#' visualizations for species invasion predictions. It saves processed
#' prediction data and map images to disk.
#'
#' @param model_prefix Character. Prefix for model output directories and files.
#' @param n_cv_folds Integer. Number of cross-validation folds. Default is 4L.
#' @param hab_abb Character. Habitat abbreviation to process.
#' @param spatial_model Logical. Whether to use a spatial model (TRUE) or
#'   non-spatial (FALSE).
#' @param n_cores Integer. Number of cores for parallel processing. Default is
#'   8L.
#' @param clamp Character vector. Indicates which prediction types to process:
#'   "clamp", "no_clamp", or both.
#' @author Ahmed El-Gabbas
#' @export

cv_merge_predictions <- function(
    model_prefix = NULL, n_cv_folds = 4L, hab_abb = NULL,
    spatial_model = NULL, n_cores = 8L, clamp = c("clamp", "no_clamp")) {

  files_clamp <- fold <- path_prefix <- dir_clamp <- dir_no_clamp <-
    files_no_clamp <- lyr <- path_pred_mean <- ias_id <- climate_model <- NULL

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  n_cores <- .validate_n_cores(n_cores)
  n_cv_folds <- .validate_n_cores(n_cv_folds)

  ecokit::check_args(args_to_check = "spatial_model", args_type = "logical")
  ecokit::check_args(args_to_check = "model_prefix", args_type = "character")

  dir_fitting <- "datasets/processed/model_fitting"
  path_prefix <- paste0(
    model_prefix, hab_abb, "_",
    dplyr::if_else(spatial_model, "gpp", "nonspatial"))
  cols_remove <- c(
    "tif_path_cov", "tif_path_sd", "cv_fold",
    "tif_path_anomaly", "Dir_Ensemble", "tifs")
  path_summary <- fs::path(
    dir_fitting, paste0(model_prefix, hab_abb, "_cv_summary"),
    paste0("predictions_", dplyr::if_else(spatial_model, "gpp", "nonspatial")))

  pred_options <- tidyr::expand_grid(
    habitat = hab_abb, spatial_model = spatial_model,
    path_prefix = path_prefix, fold = seq_len(n_cv_folds))

  # No clamping
  if ("no_clamp" %in% clamp) {

    file_summary_no_clamp <- fs::path(path_summary, "preds_no_clamp.RData")

    if (!ecokit::check_data(file_summary_no_clamp, warning = FALSE)) {
      cols_select <- c("spatial_model", "path_prefix", "data_no_clamp")
      preds_no_clamp <- pred_options %>%
        dplyr::mutate(
          dir_no_clamp = fs::path(
            dir_fitting,
            paste0(path_prefix, "_cv", fold),
            "Model_Prediction", "NoClamp"),
          files_no_clamp = purrr::map(
            dir_no_clamp,
            list.files, pattern = "_Summary.RData", full.names = TRUE),
          data_no_clamp = purrr::map2(
            .x = files_no_clamp, .y = fold,
            .f = ~ {
              purrr::map_dfr(.x, ecokit::load_as) %>%
                ecokit::arrange_alphanum(
                  time_period, climate_model, climate_scenario) %>%
                dplyr::mutate(
                  clamp = FALSE, cv_fold = .y,
                  Clamp = NULL, .after = "hab_name") %>%
                dplyr::mutate( # nolint
                  option_name = paste0(
                    ias_id, "_", time_period, "_",
                    climate_model, "_", climate_scenario),
                  option_name = stringr::str_replace_all(option_name, "-", "_"))
            })
        ) %>%
        dplyr::select(tidyselect::all_of(cols_select)) %>%
        tidyr::unnest("data_no_clamp") %>%
        dplyr::select(-tidyselect::any_of(cols_remove)) %>%
        tidyr::nest(paths = c("tif_path_mean", "option_name"))

      preds_no_clamp_paths <- summarize_cv_maps(
        preds_no_clamp, model_prefix, n_cores, hab_abb, spatial_model)
      preds_no_clamp <- dplyr::mutate(
        preds_no_clamp, maps = preds_no_clamp_paths) %>%
        tidyr::unnest("maps")

      ecokit::save_as(
        object = preds_no_clamp, object_name = "preds_no_clamp",
        out_path = file_summary_no_clamp)

      invisible(gc())
    }
  }

  if ("clamp" %in% clamp) {

    file_summary_clamp <- fs::path(path_summary, "preds_clamp.RData")

    if (ecokit::check_data(file_summary_clamp, warning = FALSE)) {
      preds_clamp <- ecokit::load_as(file_summary_clamp)
    } else {
      cols_select <- c("spatial_model", "path_prefix", "data_clamp")
      preds_clamp <- pred_options %>%
        dplyr::mutate(
          dir_clamp = fs::path(
            dir_fitting,
            paste0(path_prefix, "_cv", fold), "Model_Prediction", "Clamp"),
          files_clamp = purrr::map(
            dir_clamp, list.files, pattern = "_Summary.RData",
            full.names = TRUE),
          data_clamp = purrr::map2(
            .x = files_clamp, .y = fold,
            .f = ~ {
              purrr::map_dfr(.x, ecokit::load_as) %>%
                ecokit::arrange_alphanum(
                  time_period, climate_model, climate_scenario) %>%
                dplyr::mutate(
                  clamp = TRUE, cv_fold = .y,
                  Clamp = NULL, .after = "hab_name") %>%
                dplyr::mutate( # nolint
                  option_name = paste0(
                    ias_id, "_", time_period, "_",
                    climate_model, "_", climate_scenario),
                  option_name = stringr::str_replace_all(option_name, "-", "_"))
            })) %>%
        dplyr::select(tidyselect::all_of(cols_select)) %>%
        tidyr::unnest("data_clamp") %>%
        dplyr::select(-tidyselect::any_of(cols_remove)) %>%
        tidyr::nest(paths = c("tif_path_mean", "option_name"))
      preds_clamp_paths <- summarize_cv_maps(
        preds_clamp, model_prefix, n_cores, hab_abb, spatial_model)
      preds_clamp <- dplyr::mutate(preds_clamp, maps = preds_clamp_paths) %>%
        tidyr::unnest("maps")

      ecokit::save_as(
        object = preds_clamp, object_name = "preds_clamp",
        out_path = file_summary_clamp)

      invisible(gc())
    }

    preds_sr <- dplyr::filter(
      preds_clamp, ias_id == "SR",
      climate_model %in% c("Current", "Ensemble")) %>%
      dplyr::pull(path_pred_mean) %>%
      terra::rast()
    preds_anomaly <- preds_sr[[2:10]] - preds_sr[[1]]

    Xlim <- c(2600000, 6500000)
    Ylim <- c(1450000, 5420000)

    sr_map <- ggplot2::ggplot(environment = emptyenv()) +
      tidyterra::geom_spatraster(data = preds_sr[[1]], maxcell = Inf)  +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", palette = "viridis::plasma", name = NULL) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim,
        oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
      ggplot2::labs(title = "Level of invasion") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0.05, 0, 0.05, "cm"),
        plot.title = ggplot2::element_text(
          size = 9, color = "grey60", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0)),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        legend.key.size = grid::unit(0.6, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.875, 0.775),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 6),
        legend.box.spacing = grid::unit(0, "pt"),
        axis.text.x = ggplot2::element_text(size = 5),
        axis.text.y = ggplot2::element_text(size = 5, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        axis.title = ggplot2::element_blank(),
        panel.spacing = grid::unit(0.2, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.05, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        panel.ontop = TRUE,
        panel.background = ggplot2::element_rect(fill = NA))

    ragg::agg_jpeg(
      filename = fs::path(path_summary, "sr_map.jpeg"),
      width = 16, height = 16, res = 600, quality = 100, units = "cm")
    print(sr_map)
    grDevices::dev.off()

    # Build custom titles dynamically
    custom_titles <- names(preds_anomaly) %>%
      stringr::str_remove_all("SR_|_Ensemble_|_mean") %>%
      stringr::str_replace_all("ssp", " --- ssp") %>%
      stringr::str_replace_all("_", "-")
    names(custom_titles) <-  names(preds_anomaly)

    year_labeller <- function(x) {
      stringr::str_replace(stringr::str_extract(x, "\\d{4}_\\d{4}"), "_", "-")
    }
    scenario_labeller <- function(x) {
      stringr::str_extract(x, "ssp\\d{3}")
    }
    min_value <- min(terra::values(preds_anomaly), na.rm = TRUE)
    max_value <- max(terra::values(preds_anomaly), na.rm = TRUE)
    lim <- max(abs(min_value), abs(max_value))

    label_plusminus <- function(x) {
      sign_char <- dplyr::case_when(x > 0 ~ "+", x < 0 ~ "-", .default =  " ")
      num_str <- as.character(abs(x))
      max_digits <- max(nchar(as.character(abs(x)))) + 1
      num_str_padded <- stringr::str_pad(
        num_str, width = max_digits, side = "left", pad = " ")
      paste0(sign_char, " ", num_str_padded)
    }

    sr_anomaly <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = preds_anomaly, maxcell = Inf) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", palette = "viridis::plasma", name = NULL,
        breaks = pretty(c(-lim, lim), n = 7),
        labels = label_plusminus(pretty(c(-lim, lim), n = 7))) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(year = year_labeller(lyr)),
        cols = ggplot2::vars(scenario = scenario_labeller(lyr))) +
      ggplot2::labs(title = "Prediction anomaly") +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim,
        oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title = ggplot2::element_text(
          size = 14, color = "grey60", face = "bold", hjust = 0,
          margin = ggplot2::margin(2, 0, 2, 0)),
        plot.title.position = "plot",
        strip.background = ggplot2::element_blank(),
        strip.text = ggtext::element_markdown(
          face = "bold", colour = "blue", size = 11,
          margin = ggplot2::margin(0, 0, 0, 0), vjust = -0.25),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.3, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.955, 0.92),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(
          size = 8, margin = ggplot2::margin(0, 0, 0, 2)),
        legend.box.spacing = grid::unit(0, "pt"),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.spacing = grid::unit(0.2, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.05, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(
          fill = "transparent", colour = NA))

    ragg::agg_jpeg(
      filename = fs::path(path_summary, "sr_anomaly.jpeg"),
      width = 16, height = 16.5, res = 600, quality = 100, units = "cm")
    print(sr_anomaly)
    grDevices::dev.off()

  }
  return(invisible(NULL))
}

## |------------------------------------------------------------------------| #
# summarize_cv_maps ----
## |------------------------------------------------------------------------| #

#' Helper function to summarize raster prediction maps in parallel
#'
#' This function processes a data frame containing raster prediction file paths,
#' computes summary statistics (mean and standard deviation) for each
#' prediction, and writes the results to disk. The computation is performed in
#' parallel if more than one core is specified.
#'
#' @param data A data frame containing prediction metadata and file paths. Must
#'   include columns `path_prefix`, `paths`, `clamp`, and `ias_id`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#' @return A list of tibbles, each containing the paths to the mean and standard
#'   deviation raster files for each prediction.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

summarize_cv_maps <- function(
    data, model_prefix, n_cores, hab_abb, spatial_model) {

  n_cores <- .validate_n_cores(n_cores)

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(n_cores = n_cores, show_log = FALSE)
    on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
  }

  dir_fitting <- "datasets/processed/model_fitting"
  dir_summary <- fs::path(
    dir_fitting, paste0(model_prefix, hab_abb, "_cv_summary"),
    paste0("predictions_", dplyr::if_else(spatial_model, "gpp", "nonspatial")))
  fs::dir_create(dir_summary)

  summary_maps <- future.apply::future_lapply(
    X = seq_len(nrow(data)),
    FUN = function(row_id) {
      sub_data <- dplyr::slice(data, row_id)
      tif_paths <- magrittr::extract2(magrittr::extract2(sub_data, "paths"), 1)
      option_name <- stringr::str_replace_all(
        unique(tif_paths$option_name), "1981_2010_Current_Current", "1981_2010")

      dir_out <- basename(dirname(tif_paths$tif_path_mean)) %>%
        unique() %>%
        stringr::str_to_lower() %>%
        stringr::str_replace_all("-", "_")
      dir_out <- fs::path(
        dir_summary,
        dplyr::if_else(sub_data$clamp, "clamp", "no_clamp"), dir_out)

      file_mean <- fs::path(dir_out, paste0(sub_data$ias_id, "_mean.tif"))
      file_sd <- fs::path(dir_out, paste0(sub_data$ias_id, "_sd.tif"))
      out_tibble <- tibble::tibble(
        path_pred_mean = file_mean, path_pred_sd = file_sd)

      files_okay <- try({
        ecokit::check_tiff(file_mean, warning = FALSE) &&
          ecokit::check_tiff(file_sd, warning = FALSE)
      },
      silent = TRUE)

      if (!inherits(files_okay, "try-error") && files_okay) {
        return(out_tibble)
      }

      fs::dir_create(dir_out)

      preds <- terra::rast(tif_paths$tif_path_mean)
      gdal_o <- c("COMPRESS=DEFLATE", "TILED=YES")

      mean_pred <- terra::app(preds, "mean") %>%
        stats::setNames(paste0(option_name, "_mean"))
      terra::writeRaster(
        x = mean_pred, file = file_mean, overwrite = TRUE, gdal = gdal_o)

      sd_pred <- terra::app(preds, "sd") %>%
        stats::setNames(paste0(option_name, "_sd"))
      terra::writeRaster(
        x = sd_pred, file = file_sd, overwrite = TRUE, gdal = gdal_o)

      out_tibble
    },
    future.packages = c(
      "terra", "tibble", "dplyr", "fs", "ecokit", "stringr", "magrittr"),
    future.globals = c("data", "dir_summary"), future.seed = TRUE)

  future::plan("sequential", gc = TRUE)

  summary_maps
}
