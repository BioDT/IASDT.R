## |------------------------------------------------------------------------| #
# Plot_VarPar ----
## |------------------------------------------------------------------------| #

#' Plot variance partitioning of the model
#'
#' Plot variance partitioning of the model
#'
#' @param Model an object of class `Hmsc` or a path for the model. Only needed if one of `ModelEval` and `VarPar` arguments are `NULL`.
#' @param PlotPath String. Path to save the output file.
#' @param ModelName String. Prefix to add to the title of the plot. Default: `NULL`, which means only use 'Variance partitioning' in the title.
#' @param ModelEval Result of the `Hmsc::evaluateModelFit` function. If `ModelEval = NULL` (default), `Hmsc::evaluateModelFit` will be executed on the model object to compute measures of model fit.
#' @param ModelEvalPar Integer. Number of parallel computations for computing predicted values. This is used as the `nParallel` argument of the `Hmsc::computePredictedValues` function.
#' @param VarPar Variance partitioning. An object resulted from `Hmsc::computeVariancePartitioning`.
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @name Plot_VarPar
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Plot_VarPar <- function(
    Model, PlotPath = NULL, ModelName = NULL, ModelEval = NULL,
    ModelEvalPar = 1, VarPar = NULL, EnvFile = ".env") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- variable <- value <- NULL

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading species list
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_TaxaList <- Sys.getenv("DP_R_Mod_Path_TaxaList")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }
  SpList <- file.path(Path_TaxaList, "TaxaList.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(Species = IAS_ID, Species_name) %>%
    dplyr::mutate(
      Species = stringr::str_pad(string = Species, width = 4, pad = 0),
      Species = paste0("Sp_", Species))


  if (purrr::some(list(ModelEval, VarPar), is.null)) {

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Loading model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (magrittr::not(inherits(Model, "Hmsc"))) {
      if (inherits(Model, "character")) {
        if (magrittr::not(file.exists(Model))) {
          MSG <- "The provided path for the model does not exist"
          stop(MSG)
        }
        Model <- IASDT.R::LoadAs(Model)
        if (inherits(Model, "Hmsc")) {
        } else {
          MSG <- "The loaded model is not of an Hmsc object"
          stop(MSG)
        }
      } else {
        MSG <- "The Model has to be an Hmsc model or a path to a saved model"
        stop(MSG)
      }
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Calculate variance partitioning
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (is.null(VarPar)) {
      VarPar <- Hmsc::computeVariancePartitioning(Model)
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Compute predicted Values
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (is.null(ModelEval)) {
      preds <- Hmsc::computePredictedValues(
        hM = Model, nParallel = ModelEvalPar) %>%
        suppressWarnings()

      # Evaluate model fit
      ModelEval <- Hmsc::evaluateModelFit(hM = Model, predY = preds) %>%
        suppressWarnings()

      rm(preds)
    }
    rm(Model)
    invisible(gc())
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plot theme
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Theme <- ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.position = "none",
    axis.text.y = ggplot2::element_text(face = "italic"),
    axis.title.y = ggplot2::element_text(size = 30),
    axis.title.x = ggplot2::element_text(size = 30),
    legend.text = ggplot2::element_text(size = 30),
    legend.key.spacing.x = ggplot2::unit(0.75, "cm"))

  g_legend <- function(Plot) {
    # https://stackoverflow.com/a/13650878/3652584
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(Plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }


  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Relative variance partitioning
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  vp_df <- VarPar$vals %>%
    tibble::as_tibble(rownames = "variable") %>%
    tidyr::pivot_longer(
      cols = -variable, names_to = "Species", values_to = "value")

  VarOrder <- vp_df %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(value = mean(value)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::pull(variable)

  SpOrder <- vp_df %>%
    dplyr::left_join(SpList, by = "Species") %>%
    dplyr::mutate(variable = factor(variable, VarOrder)) %>%
    dplyr::summarise(value = sum(value), .by = c(Species_name, variable)) %>%
    dplyr::arrange(variable, dplyr::desc(value)) %>%
    dplyr::pull(Species_name) %>%
    unique()

  Plot_DT <- vp_df %>%
    dplyr::left_join(SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable,
        "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder))

  Plot <- Plot_DT %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, 1.0125), expand = FALSE) +
    Theme

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Raw variance partitioning
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  for (ii in seq_along(ModelEval$TjurR2)) {
    VarPar$vals[, ii] <- ModelEval$TjurR2[ii] * VarPar$vals[, ii]
  }

  vp_df_R <- VarPar$vals %>%
    tibble::as_tibble(rownames = "variable") %>%
    tidyr::pivot_longer(
      cols = -variable, names_to = "Species", values_to = "value")

  MaxVal <- vp_df_R %>%
    dplyr::group_by(Species) %>%
    dplyr::summarize(value = sum(value)) %>%
    dplyr::pull(value) %>%
    max() %>%
    magrittr::multiply_by(100) %>%
    ceiling() %>%
    magrittr::divide_by(100)

  Plot_R_DT <- vp_df_R %>%
    dplyr::left_join(SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable, "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder))

  Plot_raw <- Plot_R_DT %>%
    ggplot2::ggplot(ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Variance explained (raw)") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, MaxVal), expand = FALSE) +
    Theme


  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Combine plots
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Legend <- Plot +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 10)))
  Legend <- g_legend(Legend)

  if (is.null(ModelName)) {
    PlotTitle <- "Variance partitioning"
  } else {
    PlotTitle <- paste0("Variance partitioning (", ModelName, ")")
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Save plot to disk
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  gridExtra::arrangeGrob(Plot, Plot_raw, nrow = 1) %>%
    gridExtra::grid.arrange(
      Legend, nrow = 2, heights = c(15, 1),
      top = grid::textGrob(
        PlotTitle, gp = grid::gpar(fontsize = 30, font = 2))) %>%
    ggplot2::ggsave(
      filename = PlotPath, height = 16, width = 24, dpi = 600)
}
