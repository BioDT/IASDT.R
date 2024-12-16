## |------------------------------------------------------------------------| #
# Mod_Evaluation_Plot ----
## |------------------------------------------------------------------------| #

#' Generate Evaluation Plots for Hmsc models
#'
#' This function generates four diagnostic plots (RMSE, AUC, Continuous Boyce
#' Index, and Tjur-R²) to evaluate the performance of Hmsc models.
#'
#' @param ModelDir Path to the model directory containing predictions.
#' @param EnvFile Path to the environment file (`.env`) for setting paths.
#' @param FromHPC Boolean indicating whether the environment is an HPC system.
#' @name Mod_Evaluation_Plot
#' @author Ahmed El-Gabbas
#' @export

Mod_Evaluation_Plot <- function(ModelDir, EnvFile = ".env", FromHPC = TRUE) {

  Path_PA <- NCells_Naturalized <- IAS_ID <- AUC <- NCells <- RMSE <- Boyce <-
    TjurR2 <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  # # Load species summary
  IASDT.R::CatTime("Load species summary")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  SpSummary <- file.path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    stop(paste("Species summary file not found at:", SpSummary), call. = FALSE)
  }
  SpSummary <- readr::read_csv(SpSummary, show_col_types = FALSE) %>%
    dplyr::select(tidyselect::all_of(c("IAS_ID", "NCells_Naturalized"))) %>%
    dplyr::rename(NCells = NCells_Naturalized) %>%
    dplyr::mutate(
      IAS_ID = paste0(
        "Sp_", stringr::str_pad(IAS_ID, width = 4, pad = "0")))

  # # ..................................................................... ###
  # # ..................................................................... ###

  if (is.null(ModelDir) || !dir.exists(ModelDir)) {
    stop("Invalid or missing `ModelDir`", call. = FALSE)
  }

  Mod_Eval <- file.path(
    ModelDir, "Model_Evaluation", "Eval_Current_NoClamping.qs2")

  Mod_Eval <- IASDT.R::LoadAs(Mod_Eval) %>%
    dplyr::filter(IAS_ID != "SR") %>%
    dplyr::left_join(SpSummary, by = "IAS_ID")

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
    data = Mod_Eval, mapping = ggplot2::aes(x = NCells, y = RMSE)) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005)) +
    PlottingTheme
  Plot_AUC <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = NCells, y = AUC)) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(x = "Number of presence grid cells") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.001, 0.005), limits = c(NA, 1)) +
    PlottingTheme
  Plot_Boyce <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = NCells, y = Boyce)) +
    ggplot2::geom_point(alpha = 0.5, color = "blue", size = 1, shape = 19) +
    ggplot2::labs(
      x = "Number of presence grid cells", y = "Continuous Boyce index") +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0), limits = c(NA, 1.01)) +
    PlottingTheme
  Plot_Tjur2 <- ggplot2::ggplot(
    data = Mod_Eval, mapping = ggplot2::aes(x = NCells, y = TjurR2)) +
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
    filename = file.path(
      ModelDir, "Model_Evaluation", "Eval_explanatory.jpeg"),
    width = 24, height = 20, res = 600,
    quality = 100, units = "cm")
  print(Plots)
  grDevices::dev.off()

  return(invisible(NULL))
}
