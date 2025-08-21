#' Cross-validation Model Evaluation and Plotting
#'
#' performs model evaluation using cross-validation results, calculates multiple
#' metrics (AUC, Tjur R<sup>2</sup>, Boyce index, RMSE), and generates summary
#' plots for explanatory and predictive power.
#'
#' @param model_dir Character. Path to the root directory of the fitted model.
#' @param cv_type Character. Cross-validation type. One of `CV_Dist` (default)
#'   or `CV_Large`. See [mod_CV_fit()] for more details.
#'
#' @return Invisibly returns the path to the saved evaluation data file.
#' @author Ahmed El-Gabbas
#' @export

mod_evaluate_cv <- function(model_dir = NULL, cv_type = "CV_Dist") {

  metric <- AUC <- TjurR2 <- Boyce <- RMSE <- n_pres <- n_abs <- species <-
    evaluation_type <- training_RMSE <- testing_RMSE <- training_AUC <-
    testing_AUC <- training_TjurR2 <- testing_TjurR2 <- training_Boyce <-
    testing_Boyce <- NULL

  # Validate inputs -------
  cv_type <- .validate_cv_name(cv_type)

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory does not exist",
      model_dir = model_dir, include_backtrace = TRUE)
  }

  # Training data ------
  file_training <- fs::path(model_dir, "ModDT_training.RData")
  if (!ecokit::check_data(file_training)) {
    ecokit::stop_ctx(
      "Training data file does not exist", file_training = file_training,
      include_backtrace = TRUE)
  }
  training_data <- ecokit::load_as(file_training)

  # Testing data ------
  file_testing <- fs::path(model_dir, "ModDT_testing.RData")
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
    model_dir, "Model_Prediction", "NoClamp", "Prediction_Current_R.qs2")
  if (!ecokit::check_data(model_predictions)) {
    ecokit::stop_ctx(
      "Model prediction file does not exist",
      model_predictions = model_predictions, include_backtrace = TRUE)
  }
  model_predictions <- ecokit::load_as(model_predictions, unwrap_r = TRUE) %>%
    terra::subset(stringr::str_detect(names(.), "^Sp.+_mean$")) %>%
    stats::setNames(stringr::str_remove_all(names(.), "_mean$"))

  if (!all(names(model_predictions) %in% names(training_data))) {
    missing_species <- !names(model_predictions) %in% names(training_data)
    missing_species <- names(model_predictions)[missing_species]
    ecokit::stop_ctx(
      "Not all names in the `model_predictions` object exists in training data",
      missing_species = missing_species, include_backtrace = TRUE)
  }

  # Calculate evaluation data -----

  evaluation_data_wide <- purrr::map_dfr(
    .x =  names(model_predictions),
    .f = ~ {
      eval_train <- extract_evaluation(
        species = .x, pred = model_predictions,
        in_data = training_data, prefix = "training_")
      eval_test <- extract_evaluation(
        species = .x, pred = model_predictions,
        in_data = testing_data, prefix = "testing_")

      tibble::tibble(
        species = .x, cv_type = cv_type,
        cv_fold = cv_fold, eval_train, eval_test)
    })

  evaluation_data <- evaluation_data_wide %>%
    tidyr::pivot_longer(
      cols = -c("species", "cv_type", "cv_fold"),
      names_to = c("evaluation_type", ".value"),
      names_pattern = "^(training|testing)_(.*)$") %>%
    dplyr::arrange(species, evaluation_type)

  evaluation_path <- fs::path(
    model_dir, "Model_Evaluation", "Eval_testing_data.RData")
  ecokit::save_as(
    object = evaluation_data, object_name = "evaluation_data",
    out_path = evaluation_path)

  # Plotting ------

  plotting_data <- evaluation_data %>%
    dplyr::mutate(
      prevalence = n_pres / (n_pres + n_abs),
      evaluation_type = factor(
        evaluation_type, levels = c("training", "testing"))) %>%
    tidyr::pivot_longer(
      cols = c(AUC, TjurR2, Boyce, RMSE),
      names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric == "TjurR2" ~ "Tjur-r<sup>2</sup>",
        metric == "Boyce" ~ "Boyce index",
        .default = metric))

  y_min <- min(plotting_data$value, na.rm = TRUE)
  y_max <- max(plotting_data$value, na.rm = TRUE)

  # Create all 6 panels in correct order
  p_list <- purrr::map(
    .x = c("n_pres", "n_abs", "prevalence"),
    .f = ~{
      list(
        plot_panel(plotting_data, .x, "training", y_min, y_max),
        plot_panel(plotting_data, .x, "testing", y_min, y_max))
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
    filename = fs::path(model_dir, "Model_Evaluation", "Evaluation_plot.jpeg"),
    width = 20, height = 30, res = 600, quality = 100, units = "cm")
  print(final_plot)
  grDevices::dev.off()

  plot_limits <- range(
    c(
      evaluation_data$RMSE, evaluation_data$TjurR2,
      evaluation_data$AUC, evaluation_data$Boyce))
  ggplot2::ggplot(evaluation_data_wide) +
    ggplot2::geom_point(
      ggplot2::aes(training_AUC, testing_AUC), colour = "#E69F00") +
    ggplot2::geom_point(
      ggplot2::aes(training_TjurR2, testing_TjurR2), colour = "#56B4E9") +
    ggplot2::geom_point(
      ggplot2::aes(training_Boyce, testing_Boyce), colour = "#009E73") +
    ggplot2::geom_point(
      ggplot2::aes(training_RMSE, testing_RMSE), colour = "#D55E00") +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(limits = plot_limits) +
    ggplot2::scale_y_continuous(limits = plot_limits) +
    ggplot2::geom_abline(slope = 1, linetype = 2) +
    ggplot2::labs(x = "Explanatory power", y = "Predictive power") +
    ggplot2::theme_bw()

  return(invisible(evaluation_path))

}


# extract_evaluation ------
# Internal function to extract_evaluation

#' @keywords internal
#' @noRd

extract_evaluation <- function(species, pred, in_data, prefix) {

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

# plot_panel -------
# Internal function to generate plot for a given x variable and evaluation type

#' @keywords internal
#' @noRd

plot_panel <- function(data, x_var, eval_type, y_min, y_max) {

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
    AUC = "#E69F00", `Tjur-r<sup>2</sup>` = "#56B4E9",
    `Boyce index` = "#009E73", RMSE = "#D55E00")

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
