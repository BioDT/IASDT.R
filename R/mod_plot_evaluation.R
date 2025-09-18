## |------------------------------------------------------------------------| #
# plot_evaluation ----
## |------------------------------------------------------------------------| #

#' Generate plots for the explanatory power of Hmsc models
#'
#' This function generates four diagnostic plots (RMSE, AUC, Continuous Boyce
#' Index, and Tjur-R²) to evaluate the performance (explanatory power without
#' cross-validation) of Hmsc models.
#'
#' @param model_dir Character. Path to the model directory containing
#'   predictions.
#' @name plot_evaluation
#' @author Ahmed El-Gabbas
#' @export

plot_evaluation <- function(model_dir) {

  ias_id <- AUC <- RMSE <- Boyce <- TjurR2 <- n_cells <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  if (is.null(model_dir) || !dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Invalid or missing `model_dir`", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  mod_eval <- fs::path(
    model_dir, "model_evaluation", "eval_current_no_clamping.qs2.qs2")
  if (!fs::file_exists(mod_eval)) {
    ecokit::stop_ctx(
      "Model evaluation file does not exist", mod_eval_File = mod_eval,
      include_backtrace = TRUE)
  }

  # Load model data
  file_model_data <- fs::path(model_dir, "model_data.RData")
  file_model_data_subset <- fs::path(model_dir, "model_data_subset.RData")

  if (fs::file_exists(file_model_data_subset)) {
    model_data <- ecokit::load_as(file_model_data_subset)
  } else {
    model_data <- ecokit::load_as(file_model_data)
  }
  # Summarize the number of cells for each ias_id
  model_n_cells <- tibble::tibble(model_data$data_y) %>%
    dplyr::summarise_all(.funs = sum) %>%
    tidyr::pivot_longer(
      cols = tidyselect::everything(),
      names_to = "ias_id", values_to = "n_cells")
  rm(model_data, file_model_data, file_model_data_subset)

  mod_eval <- ecokit::load_as(mod_eval) %>%
    dplyr::filter(ias_id != "sr") %>%
    # add n_cells
    dplyr::left_join(model_n_cells, by = "ias_id")

  # # ..................................................................... ###
  # # ..................................................................... ###

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_line(
        linewidth = 0.025, colour = "grey40", linetype = 2),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.05, colour = "grey40", linetype = 2),
      axis.title = ggtext::element_markdown(face = "bold"))

  plot_rmse <- ggplot2::ggplot(
    data = mod_eval, mapping = ggplot2::aes(x = n_cells, y = RMSE),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005)) +
    plot_theme

  plot_auc <- ggplot2::ggplot(
    data = mod_eval, mapping = ggplot2::aes(x = n_cells, y = AUC),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005), limits = c(NA, 1)) +
    plot_theme

  plot_boyce <- ggplot2::ggplot(
    data = mod_eval, mapping = ggplot2::aes(x = n_cells, y = Boyce),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(
      x = "Number of presence grid cells", y = "Continuous Boyce index") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0), limits = c(NA, 1.01)) +
    plot_theme

  plot_tjur2 <- ggplot2::ggplot(
    data = mod_eval, mapping = ggplot2::aes(x = n_cells, y = TjurR2),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(
      x = "Number of presence grid cells", y = "Tjur-R<sup>2</sup>") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(NA, 1.035)) +
    plot_theme

  plots <- patchwork::wrap_plots(
    plot_rmse, plot_auc, plot_boyce, plot_tjur2, nrow = 2, ncol = 2) +
    patchwork::plot_layout(axes = "collect")

  ragg::agg_jpeg(
    filename = fs::path(model_dir, "model_evaluation", "Eval_explanatory.jpeg"),
    width = 24, height = 20, res = 600,
    quality = 100, units = "cm")
  print(plots)
  grDevices::dev.off()

  invisible(NULL)
}
