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

  Mod_Eval <- fs::path(
    model_dir, "Model_Evaluation", "Eval_Current_NoClamping.qs2")
  if (!fs::file_exists(Mod_Eval)) {
    ecokit::stop_ctx(
      "Model evaluation file does not exist", Mod_Eval_File = Mod_Eval,
      include_backtrace = TRUE)
  }

  # Load model data
  file_model_data <- fs::path(model_dir, "ModDT.RData")
  file_model_data_subset <- fs::path(model_dir, "ModDT_subset.RData")

  if (fs::file_exists(file_model_data_subset)) {
    model_data <- ecokit::load_as(file_model_data_subset)
  } else {
    model_data <- ecokit::load_as(file_model_data)
  }
  # Summarize the number of cells for each ias_id
  model_n_cells <- tibble::tibble(model_data$DT_y) %>%
    dplyr::summarise_all(.funs = sum) %>%
    tidyr::pivot_longer(
      cols = tidyselect::everything(),
      names_to = "ias_id", values_to = "n_cells")
  rm(model_data, file_model_data, file_model_data_subset)

  Mod_Eval <- ecokit::load_as(Mod_Eval) %>%
    dplyr::filter(ias_id != "SR") %>%
    # add n_cells
    dplyr::left_join(model_n_cells, by = "ias_id")

  # # ..................................................................... ###
  # # ..................................................................... ###

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_line(
        linewidth = 0.025, colour = "grey40", linetype = 2),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.05, colour = "grey40", linetype = 2),
      axis.title = ggtext::element_markdown(face = "bold"))

  Plot_RMSE <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = n_cells, y = RMSE),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005)) +
    PlottingTheme

  Plot_AUC <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = n_cells, y = AUC),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005), limits = c(NA, 1)) +
    PlottingTheme

  Plot_Boyce <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = n_cells, y = Boyce),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(
      x = "Number of presence grid cells", y = "Continuous Boyce index") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0), limits = c(NA, 1.01)) +
    PlottingTheme

  Plot_Tjur2 <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = n_cells, y = TjurR2),
    environment = emptyenv()) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(
      x = "Number of presence grid cells", y = "Tjur-R<sup>2</sup>") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(NA, 1.035)) +
    PlottingTheme

  Plots <- patchwork::wrap_plots(
    Plot_RMSE, Plot_AUC, Plot_Boyce, Plot_Tjur2, nrow = 2, ncol = 2) +
    patchwork::plot_layout(axes = "collect")

  ragg::agg_jpeg(
    filename = fs::path(model_dir, "Model_Evaluation", "Eval_explanatory.jpeg"),
    width = 24, height = 20, res = 600,
    quality = 100, units = "cm")
  print(Plots)
  grDevices::dev.off()

  return(invisible(NULL))
}
