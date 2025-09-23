## |------------------------------------------------------------------------| #
# mod_cv_evaluate ----
## |------------------------------------------------------------------------| #

#' Cross-validation Model Evaluation and Plotting
#'
#' performs model evaluation using cross-validation results, calculates multiple
#' metrics (AUC, Tjur R<sup>2</sup>, Boyce index, RMSE), and generates summary
#' plots for explanatory and predictive power.
#'
#' @param model_dir Character. Path to the root directory of the fitted model.
#' @param cv_type Character. Cross-validation type. One of `cv_dist` (default)
#'   or `cv_large`. See [mod_cv_fit()] for more details.
#'
#' @return Invisibly returns the path to the saved evaluation data file.
#' @author Ahmed El-Gabbas
#' @export

mod_cv_evaluate <- function(model_dir = NULL, cv_type = "cv_dist") {

  metric <- AUC <- TjurR2 <- Boyce <- RMSE <- n_pres <- n_abs <- species <-
    evaluation_type <- value <- training <- testing <- NULL

  # Validate inputs -------
  cv_type <- .validate_cv_name(cv_type)

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory does not exist",
      model_dir = model_dir, include_backtrace = TRUE)
  }

  # Training data ------
  file_training <- fs::path(model_dir, "model_data_training.RData")
  if (!ecokit::check_data(file_training)) {
    ecokit::stop_ctx(
      "Training data file does not exist", file_training = file_training,
      include_backtrace = TRUE)
  }
  training_data <- ecokit::load_as(file_training)

  # Testing data ------
  file_testing <- fs::path(model_dir, "model_data_testing.RData")
  if (!ecokit::check_data(file_testing)) {
    ecokit::stop_ctx(
      "Testing data file does not exist", file_testing = file_testing,
      include_backtrace = TRUE)
  }
  testing_data <- ecokit::load_as(file_testing)

  # cv_fold -----
  if (!(cv_type %in% colnames(testing_data))) {
    ecokit::stop_ctx(
      paste0("Column '", cv_type, "' not found in testing_data"),
      cv_type = cv_type, names_testing_data = names(testing_data),
      include_backtrace = TRUE)
  }
  cv_fold <- unique(dplyr::pull(testing_data, !!cv_type))
  if (length(cv_fold) != 1) {
    ecokit::stop_ctx(
      "`cv_fold` has to be of length 1",
      cv_fold = cv_fold, cv_type = cv_type, include_backtrace = TRUE)
  }

  # Model predictions -----
  model_predictions <- fs::path(
    model_dir, "model_prediction", "no_clamp", "prediction_current_r.qs2")
  if (!ecokit::check_data(model_predictions)) {
    ecokit::stop_ctx(
      "Model prediction file does not exist",
      model_predictions = model_predictions, include_backtrace = TRUE)
  }
  model_predictions <- ecokit::load_as(model_predictions, unwrap_r = TRUE) %>%
    terra::subset(stringr::str_detect(names(.), "^sp.+_mean$")) %>%
    stats::setNames(stringr::str_remove_all(names(.), "_mean$"))

  if (!all(names(model_predictions) %in% names(training_data))) {
    missing_species <- !names(model_predictions) %in% names(training_data)
    missing_species <- names(model_predictions)[missing_species]
    ecokit::stop_ctx(
      "Not all names in the `model_predictions` object exists in training data",
      missing_species = missing_species, include_backtrace = TRUE)
  }

  # Calculate evaluation data -----

  eval_cv_data_wide <- purrr::map_dfr(
    .x =  names(model_predictions),
    .f = ~ {
      eval_train <- cv_extract_eval(
        species = .x, pred = model_predictions,
        in_data = training_data, prefix = "training_")
      eval_test <- cv_extract_eval(
        species = .x, pred = model_predictions,
        in_data = testing_data, prefix = "testing_")

      tibble::tibble(
        species = .x, cv_type = cv_type,
        cv_fold = cv_fold, eval_train, eval_test)
    })

  eval_cv_data <- eval_cv_data_wide %>%
    tidyr::pivot_longer(
      cols = -c("species", "cv_type", "cv_fold"),
      names_to = c("evaluation_type", ".value"),
      names_pattern = "^(training|testing)_(.*)$") %>%
    dplyr::arrange(species, evaluation_type)

  evaluation_path <- fs::path(
    model_dir, "model_evaluation", "eval_cv_data.RData")
  ecokit::save_as(
    object = eval_cv_data, object_name = "eval_cv_data",
    out_path = evaluation_path)

  # Plotting ------

  plotting_data <- eval_cv_data %>%
    dplyr::mutate(
      prevalence = n_pres / (n_pres + n_abs),
      evaluation_type = factor(
        evaluation_type, levels = c("training", "testing"))) %>%
    tidyr::pivot_longer(
      cols = c(AUC, TjurR2, Boyce, RMSE),
      names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric == "TjurR2" ~ "Tjur-R<sup>2</sup>",
        metric == "Boyce" ~ "Boyce index",
        .default = metric)) %>%
    # exclude NA metrics; e.g. if there is no testing evaluation data
    dplyr::filter(!is.na(value))

  y_min <- min(plotting_data$value, na.rm = TRUE)
  y_max <- max(plotting_data$value, na.rm = TRUE)

  # Create all 6 panels in correct order
  p_list <- purrr::map(
    .x = c("n_pres", "n_abs", "prevalence"),
    .f = ~{
      list(
        cv_plot_panel1(plotting_data, .x, "training", y_min, y_max),
        cv_plot_panel1(plotting_data, .x, "testing", y_min, y_max))
    }) %>%
    unlist(recursive = FALSE)
  rm(y_min, y_max, envir = environment())

  col_titles <- list(
    patchwork::wrap_elements(
      grid::textGrob(
        "Explanatory power", vjust = 0.5, y = 0.65,
        gp = grid::gpar(fontsize = 16, fontface = "bold", col = "red"))),
    patchwork::wrap_elements(
      grid::textGrob(
        "Predictive power", vjust = 0.5, y = 0.65,
        gp = grid::gpar(fontsize = 16, fontface = "bold", col = "red"))))

  final_plot_theme <- ggplot2::theme(
    legend.position = "bottom",
    legend.key = ggplot2::element_blank(),
    legend.text = ggtext::element_markdown(size = 16),
    legend.title = ggplot2::element_text(size = 12),
    legend.margin = ggplot2::margin(t = -5))
  final_plot <- patchwork::wrap_plots(
    p_list, ncol = 2, nrow = 3, byrow = TRUE, guides = "collect")
  final_plot <- (
    patchwork::wrap_plots(col_titles, ncol = 2) / final_plot) & final_plot_theme
  final_plot <- final_plot +
    patchwork::plot_layout(heights = c(0.25, 12), widths = c(1, 1))

  ragg::agg_jpeg(
    filename = fs::path(model_dir, "model_evaluation", "evaluation_plot.jpeg"),
    width = 20, height = 30, res = 600, quality = 100, units = "cm")
  print(final_plot)
  grDevices::dev.off()

  # Explanatory vs predictive power
  plot_limits <- range(
    c(
      eval_cv_data$RMSE, eval_cv_data$TjurR2,
      eval_cv_data$AUC, eval_cv_data$Boyce),
    na.rm = TRUE)

  point_cols <- c(
    AUC = "#e66101", `Tjur-R<sup>2</sup>` = "#fdb863",
    `Boyce index` = "#b2abd2", RMSE = "#5e3c99")

  training_vs_testing <- plotting_data %>%
    dplyr::select(species, cv_type, cv_fold, metric, evaluation_type, value) %>%
    tidyr::pivot_wider(
      id_cols = c(species, cv_type, cv_fold, metric),
      names_from = evaluation_type, values_from = value) %>%
    dplyr::filter(!is.na(training), !is.na(testing)) %>%
    ggplot2::ggplot(ggplot2::aes(x = training, y = testing, colour = metric)) +
    ggplot2::geom_abline(slope = 1, linetype = 2, colour = "lightgrey") +
    ggplot2::geom_point(size = 2, alpha = 0.75, shape = 16) +
    ggplot2::labs(x = "Explanatory power", y = "Predictive power") +
    ggplot2::coord_equal(
      expand = FALSE, xlim = plot_limits, ylim = plot_limits, clip = "off") +
    ggplot2::scale_color_manual(
      values = point_cols, name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, shape = 16))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 5, r = 15, b = 2, l = 5),
      axis.title = ggplot2::element_text(size = 18, face = "bold"),
      axis.text = ggplot2::element_text(size = 14),
      legend.position = "inside",
      legend.direction = "horizontal",
      legend.position.inside = c(0.22, 0.975),
      legend.key.width = grid::unit(0.5, "lines"),
      legend.key.height = grid::unit(0.5, "lines"),
      legend.spacing.x = grid::unit(0, "pt"),
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggtext::element_markdown(size = 18))

  ragg::agg_jpeg(
    filename = fs::path(
      model_dir, "model_evaluation", "training_vs_testing.jpeg"),
    width = 30.5, height = 30, res = 600, quality = 100, units = "cm")
  print(training_vs_testing)
  grDevices::dev.off()


  return(invisible(evaluation_path))

}

## |------------------------------------------------------------------------| #
# cv_extract_eval ----
## |------------------------------------------------------------------------| #

#' Internal function to cv_extract_eval
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @noRd

cv_extract_eval <- function(species, pred, in_data, prefix) {

  selected_cols <- c("x", "y", species)
  pred_species <- terra::subset(pred, species)
  pres_abs <- dplyr::select(in_data, tidyselect::all_of(selected_cols)) %>%
    dplyr::mutate(
      pred = terra::extract(
        x = pred_species, y = .[, c("x", "y")], ID = FALSE)[, 1L])
  pred <- pres_abs$pred
  pres_abs <- pres_abs[, species, drop = TRUE]
  n_pres <- sum(pres_abs == 1)
  n_abs <- sum(pres_abs == 0)

  if (n_pres == 0 || n_abs == 0) {
    # Do not calculate evaluation if there is 0 presences or absences
    RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
  } else {
    RMSE <- caret::RMSE(pres_abs, pred)
    mean_pres <- mean(pred[which(pres_abs == 1)])
    mean_abs <- mean(pred[which(pres_abs == 0)])
    TjurR2 <- mean_pres - mean_abs
    AUC <- pROC::auc(
      response = pres_abs, predictor = pred,
      levels = c(0, 1), direction = "<") %>%
      as.numeric()
    Boyce <- ecospat::ecospat.boyce(
      fit = pred, obs = pred[pres_abs == 1], PEplot = FALSE)$cor
  }

  tibble::tibble(
    n_pres = n_pres, n_abs = n_abs, RMSE = RMSE, TjurR2 = TjurR2,
    AUC = AUC, Boyce = Boyce) %>%
    dplyr::rename_with(~ paste0(prefix, .x))
}

## |------------------------------------------------------------------------| #
# cv_plot_panel1 ----
## |------------------------------------------------------------------------| #

#' Internal function to generate plot for a given x variable and evaluation type
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @noRd

cv_plot_panel1 <- function(data, x_var, eval_type, y_min, y_max) {

  evaluation_type <- value <- metric <- NULL

  if (eval_type == "training") {
    axis_title <- ggplot2::element_text(
      margin = ggplot2::margin(r = 8), vjust = 0.5)
    axis_text <- ggplot2::element_text(size = 11)
    left_margin <- 3
  }  else {
    axis_title <- axis_text <- ggplot2::element_blank()
    left_margin <- 20
  }
  x_lab <- dplyr::case_when(
    x_var == "n_pres" ~ paste0("Number of ", eval_type, " presences"),
    x_var == "n_abs" ~ paste0("Number of ", eval_type, " absences"),
    x_var == "prevalence" ~
      stringr::str_to_sentence(paste0(eval_type, " prevalence")),
    .default = x_var)

  point_cols <- c(
    AUC = "#e66101", `Tjur-R<sup>2</sup>` = "#fdb863",
    `Boyce index` = "#b2abd2", RMSE = "#5e3c99")

  plot <- data %>%
    dplyr::filter(evaluation_type == eval_type) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data[[x_var]], y = value, color = metric)) +
    ggplot2::geom_point(size = 1.25, alpha = 0.75, shape = 16) +
    ggplot2::scale_color_manual(
      values = point_cols, name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, shape = 16))) +
    ggplot2::labs(x = x_lab, y = NULL) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, y_max), expand = c(0.0125, 0)) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(
        fill = "transparent", color = NA),
      panel.grid.major = ggplot2::element_line(
        color = "grey80", linewidth = 0.35, linetype = 2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        t = 3, l = left_margin, b = 3, r = 2, unit = "pt"),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = 13, colour = "blue", face = "bold"),
      axis.title.y = axis_title,
      axis.text.y  = axis_text)

  plot <- plot +
    ggplot2::scale_x_continuous(
      labels = if (x_var %in% c("n_pres", "n_abs")) {
        scales::comma
      } else {
        ggplot2::waiver()
      },
      expand = c(0.01, 0))
  plot
}

## |------------------------------------------------------------------------| #
# mod_evaluate_cv_plot ----
## |------------------------------------------------------------------------| #

#' Plot Evaluation Results for Cross-Validated Hmsc Models
#'
#' This function evaluates cross-validation results for Hmsc models (spatial
#' and/or non-spatial), summarizes performance metrics (RMSE, TjurR2, AUC, and
#' Boyce), and generates diagnostic plots comparing model types.
#' @param model_prefix Character. Prefix for model directory name.
#' @param hab_abb Character. Habitat abbreviation.
#' @param n_cv_folds Integer. Number of cross-validation folds (default: 4L).
#' @param spatial_model Character vector. Specifies which models to evaluate:
#'   "gpp", "nonspatial", or both.
#' @author Ahmed El-Gabbas
#' @export

mod_cv_evaluate_plot <- function(
    model_prefix = NULL, hab_abb = NULL, n_cv_folds = 4,
    spatial_model = c("gpp", "nonspatial")) {

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  n_cv_folds <- .validate_n_cores(n_cv_folds)
  ecokit::check_args(args_to_check = "model_prefix", args_type = "character")

  if (!inherits(spatial_model, "character") ||
      length(spatial_model) > 2 || (length(spatial_model) == 0)) {
    ecokit::stop_ctx("`spatial_model` has to be of length 1 or 2")
  }
  if (!all(spatial_model %in% c("gpp", "nonspatial"))) {
    ecokit::stop_ctx(
      "`spatial_model` has to be either of gpp or nonspatial or both")
  }

  spatial <- value <- metric <- evaluation_type <- cv_fold <- nonspatial <-
    path_prefix <- fold <- gpp <- stat <- NULL

  dir_fit <- "datasets/processed/model_fitting"
  metric_names <- c("RMSE", "TjurR2", "AUC", "Boyce")
  dir_eval <- fs::path(dir_fit, paste0(model_prefix, hab_abb, "_cv_summary"))

  eval_summary <- tidyr::expand_grid(
    hab_abb = hab_abb, spatial = spatial_model, fold =  seq_len(n_cv_folds)) %>%
    dplyr::mutate(
      path_prefix = paste0(model_prefix, hab_abb, "_", spatial),
      eval_file = fs::path(
        dir_fit, paste0(path_prefix, "_cv", fold),
        "model_evaluation", "eval_cv_data.RData"))

  eval_files <- eval_summary$eval_file
  eval_okay <- purrr::map_lgl(eval_files, ecokit::check_data, warning = FALSE)
  if (!all(eval_okay)) {
    ecokit::stop_ctx(
      "Not all model directories exist",
      eval_files = eval_files, missing_paths = eval_files[!eval_okay],
      include_backtrace = TRUE)
  }

  eval_summary <- eval_summary %>%
    dplyr::select(-fold, -path_prefix) %>%
    tidyr::nest(eval = -c("hab_abb", "spatial")) %>%
    dplyr::mutate(
      eval = purrr::map(eval, ~ .x$eval_file),
      eval_data = purrr::map(
        .x = eval,
        .f = ~ {
          data <- purrr::map_dfr(.x, ecokit::load_as)
          summary_data <- dplyr::select(data, -cv_fold) %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::everything(),
                .fns = list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
              .by = c("species", "cv_type", "evaluation_type")) %>%
            tidyr::pivot_longer(
              cols = -c(species, cv_type, evaluation_type),
              names_to = c(".value", "cv_fold"),
              names_pattern = "(n_pres|n_abs|RMSE|TjurR2|AUC|Boyce)_(mean|sd)")

          dplyr::mutate(data, cv_fold = as.character(cv_fold)) %>%
            dplyr::bind_rows(summary_data) %>%
            dplyr::arrange(species, evaluation_type, cv_fold)
        }))

  metric_ranges <- dplyr::select(eval_summary, -tidyselect::all_of("eval")) %>%
    tidyr::unnest("eval_data") %>%
    dplyr::filter(cv_fold %in% c("mean", "sd")) %>%
    dplyr::select(
      -tidyselect::all_of(c("n_pres", "n_abs", "cv_type", "hab_abb"))) %>%
    tidyr::pivot_wider(
      names_from = cv_fold, values_from = tidyselect::all_of(metric_names)) %>%
    dplyr::select(
      -tidyselect::all_of(c("spatial", "species", "evaluation_type"))) %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::ends_with("_mean"),
        .fns = list(
          plus = ~ . + get(
            stringr::str_replace(dplyr::cur_column(), "_mean", "_sd")),
          minus = ~ . - get(
            stringr::str_replace(dplyr::cur_column(), "_mean", "_sd"))),
        .names = "{stringr::str_replace(.col, '_mean', '')}_{fn}")) %>%
    dplyr::select(
      tidyselect::ends_with("plus"), tidyselect::ends_with("minus")) %>%
    dplyr::summarise(
      dplyr::across(
        tidyselect::ends_with("_plus"),
        ~ max(.x, na.rm = TRUE), .names = "{.col}_max"),
      dplyr::across(
        tidyselect::ends_with("_minus"),
        ~ min(.x, na.rm = TRUE), .names = "{.col}_min")) %>%
    tidyr::pivot_longer(
      tidyselect::everything(), names_to = c("metric", "stat"),
      names_pattern = "(.*)_(plus|max|minus|min)", values_to = "value") %>%
    dplyr::mutate(metric = stringr::str_remove(metric, "_plus|_minus")) %>%
    tidyr::pivot_wider(names_from = stat, values_from = value)

  if ("gpp" %in% spatial_model) {
    ragg::agg_jpeg(
      filename = fs::path(dir_eval, "eval_spatial.jpeg"),
      width = 30, height = 30, res = 600, quality = 100, units = "cm")
    plot <- cv_plot_panel2(
      plot_ex = eval_summary$eval_data[[1]], metric_ranges = metric_ranges,
      plot_title = "Spatial Hmsc", metrics = metric_names)
    print(plot)
    grDevices::dev.off()
    rm(plot, envir = environment())
  }

  if ("nonspatial" %in% spatial_model) {
    ragg::agg_jpeg(
      filename = fs::path(dir_eval, "eval_nonspatial.jpeg"),
      width = 30, height = 30, res = 600, quality = 100, units = "cm")
    plot <- cv_plot_panel2(
      plot_ex = eval_summary$eval_data[[2]], metric_ranges = metric_ranges,
      plot_title = "Non-spatial Hmsc", metrics = metric_names)
    print(plot)
    grDevices::dev.off()
    rm(plot, envir = environment())
  }

  if (length(spatial_model) == 2 &&
      identical(spatial_model, c("gpp", "nonspatial"))) {
    cols_remove <- c(
      "spatial", "species", "cv_type", "cv_fold", "evaluation_type",
      metric_names)
    eval_summary2 <- dplyr::select(eval_summary, -eval) %>%
      tidyr::unnest("eval_data") %>%
      dplyr::filter(cv_fold != "sd") %>%
      dplyr::select(tidyselect::all_of(cols_remove)) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(metric_names),
        names_to = "metric", values_to = "value") %>%
      dplyr::mutate(
        metric = dplyr::case_when(
          metric == "Boyce" ~ "Continuous Boyce index",
          metric == "AUC" ~ "Area under the ROC curve",
          metric == "TjurR2" ~ "Tjur-R<sup>2</sup>",
          .default = metric))

    metric_ranges <- dplyr::group_by(eval_summary2, metric) %>%
      dplyr::summarise(
        min = min(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE), .groups = "drop")

    spatial_nonspatial_theme <- ggplot2::theme(
      plot.margin = ggplot2::margin(4, 10, 6, 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "lightgrey", linewidth = 0.15),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggtext::element_markdown(
        face = "bold", size = 20, colour = "blue"),
      plot.title.position = "plot",
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggtext::element_markdown(
        size = 16, margin = ggplot2::margin(2, 0, 2, 0),
        hjust = 0, face = "bold"),
      strip.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 24, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(
        colour = "grey90", fill = NA, linewidth = 0.6))

    spatial_nonspatial_training <- eval_summary2 %>%
      dplyr::filter(evaluation_type == "training") %>%
      tidyr::pivot_wider(names_from = spatial, values_from = "value") %>%
      dplyr::left_join(metric_ranges, by = "metric") %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(x = nonspatial, y = gpp, colour = cv_fold),
        shape = 16, size = 1.125, show.legend = FALSE, alpha = 0.8) +
      ggplot2::geom_abline(
        slope = 1, linewidth = 0.75, linetype = "dashed", colour = "grey") +
      ggplot2::geom_blank(ggplot2::aes(x = min, y = min)) +
      ggplot2::geom_blank(ggplot2::aes(x = max, y = max)) +
      ggplot2::facet_wrap(~metric, scales = "free") +
      ggplot2::coord_cartesian(expand = FALSE) +
      ggplot2::labs(
        x = "Non-spatial Hmsc", y = "Spatial Hmsc",
        title = "Spatial vs non-spatial Hmsc --- exploratory power") +
      spatial_nonspatial_theme

    ragg::agg_jpeg(
      filename = fs::path(dir_eval, "spatial_nonspatial_training.jpeg"),
      width = 30.5, height = 30, res = 600, quality = 100, units = "cm")
    print(spatial_nonspatial_training)
    grDevices::dev.off()

    spatial_nonspatial_testing <- eval_summary2 %>%
      dplyr::filter(evaluation_type == "testing") %>%
      tidyr::pivot_wider(names_from = spatial, values_from = "value") %>%
      dplyr::left_join(metric_ranges, by = "metric") %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(x = nonspatial, y = gpp, colour = cv_fold),
        shape = 16, size = 1.125, show.legend = FALSE, alpha = 0.8) +
      ggplot2::geom_abline(
        slope = 1, linewidth = 0.75, linetype = "dashed", colour = "grey") +
      ggplot2::geom_blank(ggplot2::aes(x = min, y = min)) +
      ggplot2::geom_blank(ggplot2::aes(x = max, y = max)) +
      ggplot2::facet_wrap(~metric, scales = "free") +
      ggplot2::coord_cartesian(expand = FALSE) +
      ggplot2::labs(
        x = "Non-spatial Hmsc", y = "Spatial Hmsc",
        title = "Spatial vs non-spatial Hmsc --- predictive power") +
      spatial_nonspatial_theme

    ragg::agg_jpeg(
      filename = fs::path(dir_eval, "spatial_nonspatial_testing.jpeg"),
      width = 30.5, height = 30, res = 600, quality = 100, units = "cm")
    print(spatial_nonspatial_testing)
    grDevices::dev.off()
  }
  invisible(NULL)
}

## |------------------------------------------------------------------------| #
# cv_plot_panel2 ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @keywords internal
#' @noRd

cv_plot_panel2 <- function(
    plot_ex, metric_ranges, plot_title = "",
    metrics = c("RMSE", "TjurR2", "AUC", "Boyce")) {

  training <- NULL

  plot_list <- purrr::map(
    .x = metrics,
    .f = function(m) {

      cv_fold <- value <- evaluation_type <- testing_mean <- testing_sd <-
        training_mean <- training_sd <- testing <- training_minus <- metric <-
        testing_minus <- training_plus <- testing_plus <- NULL

      axis_range <- dplyr::filter(metric_ranges, metric == m)
      axis_range <- c(axis_range$min, axis_range$max)

      m_data <- plot_ex %>%
        dplyr::select(
          tidyselect::all_of(
            c("species", "cv_type", "cv_fold", "evaluation_type", m))) %>%
        dplyr::rename(value = !!m) %>%
        tidyr::pivot_wider(
          names_from = evaluation_type, values_from = value) %>%
        # Exclude missing values; e.g. for species with no testing data in
        # spatial cross-validation
        dplyr::filter(!is.na(testing), !is.na(training))

      points_df <- dplyr::filter(m_data, !cv_fold %in% c("mean", "sd"))

      error_bars <- dplyr::filter(m_data, cv_fold %in% c("mean", "sd")) %>%
        tidyr::pivot_wider(
          names_from = cv_fold, values_from = c(training, testing)) %>%
        dplyr::mutate(
          training_plus = training_mean + training_sd,
          training_minus = training_mean - training_sd,
          testing_plus = testing_mean + testing_sd,
          testing_minus = testing_mean - testing_sd)

      if (m %in% c("AUC", "TjurR2")) {
        axis_range <- pmin(pmax(axis_range, 0), 1)
      }
      if (m %in% c("AUC", "Boyce")) {
        axis_range <- pmin(pmax(axis_range, -1), 1)
      }
      if (m == "RMSE") {
        axis_range <- pmax(axis_range, 0)
      }

      plot_subtitle <- dplyr::case_when(
        m == "Boyce" ~ "Continuous Boyce index",
        m == "AUC" ~ "Area under the ROC curve",
        m == "TjurR2" ~ "Tjur-R<sup>2</sup>",
        .default = m)

      ggplot2::ggplot() +
        ggplot2::geom_errorbar(
          data = error_bars,
          ggplot2::aes(
            x = training_mean, ymin = testing_minus, ymax = testing_plus),
          color = "red", linewidth = 0.25, alpha = 0.65) +
        ggplot2::geom_errorbarh(
          data = error_bars,
          ggplot2::aes(
            y = testing_mean, xmin = training_minus, xmax = training_plus),
          color = "red", linewidth = 0.25, alpha = 0.65) +
        ggplot2::geom_abline(
          slope = 1, linewidth = 0.5, linetype = "dashed",
          colour = "lightgrey") +
        ggplot2::geom_point(
          data = points_df, ggplot2::aes(x = training, y = testing),
          color = "blue", size = 1, shape = 16, alpha = 0.65) +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = -Inf, y = Inf, label = plot_subtitle), colour = "grey60",
          hjust = 0, vjust = 1, size = 8, fill = "white", label.color = NA) +
        ggplot2::coord_equal(
          xlim = axis_range, ylim = axis_range, expand = FALSE, clip = "on") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          plot.margin = ggplot2::margin(4, 3, 3, 3),
          axis.text = ggplot2::element_text(size = 12))
    })

  combined_patch <- patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = plot_title,
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          size = 25, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(t = 5, b = 5))))

  final_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(
      combined_patch, x = 0.02, y = 0.02, width = 0.98, height = 0.99) +
    cowplot::draw_label(
      "Explanatory power", x = 0.5, y = 0.01,
      vjust = 0, fontface = "bold", size = 20) +
    cowplot::draw_label(
      "Predictive power", x = 0.01, y = 0.5, angle = 90,
      vjust = 1, fontface = "bold", size = 20)

  return(final_plot)
}
