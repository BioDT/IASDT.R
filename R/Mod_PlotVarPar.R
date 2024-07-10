## |------------------------------------------------------------------------| #
# PlotVarPar ----
## |------------------------------------------------------------------------| #

#' Plot variance partitioning of the model
#'
#' Plot variance partitioning of the model
#'
#' @param Model an object of class `Hmsc` or a path for the model. Only needed if one of `ModelEval` and `VarPar` arguments are `NULL`.
#' @param Path_Plot String. Path to save the output file.
#' @param PlotTitlePrefix String. Prefix to add to the title of the plot. Default: `NULL`, which means only use 'Variance partitioning' in the title.
#' @param ModelEval Result of the `Hmsc::evaluateModelFit` function. If `ModelEval = NULL` (default), `Hmsc::evaluateModelFit` will be executed on the model object to compute measures of model fit.
#' @param NCores Integer. Number of parallel computations for computing predicted values. This is used as the `nParallel` argument of the `Hmsc::computePredictedValues` function.
#' @param VarPar Variance partitioning. An object resulted from `Hmsc::computeVariancePartitioning`. if `ModelEval = NULL` (default), `Hmsc::computeVariancePartitioning` will be executed on the model object.
#'
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param SaveVarPar Logical. If `VarPar = NULL`, should the calculated variance partitioning be saved as RData file?
#' @param SaveModelEval Logical. If `ModelEval = NULL`, should the calculated model evaluation be be saved as RData file?
#' @param ReturnGG Logical. Return the plot object. Default: `FALSE`, which does not return anything
#' @param SaveGG Logical. Should the plots be exported as RData object?
#' @name PlotVarPar
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotVarPar <- function(
    Model = NULL, Path_Plot = NULL, PlotTitlePrefix = NULL, ModelEval = NULL,
    NCores = NULL, VarPar = NULL, EnvFile = ".env",
    SaveVarPar = TRUE, SaveModelEval = TRUE, ReturnGG = FALSE, SaveGG = TRUE) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- variable <- value <- NULL

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Check input arguments ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check input arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Path_Plot", "EnvFile"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("SaveVarPar", "SaveModelEval", "ReturnGG", "SaveGG"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")

  rm(AllArgs)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading species list
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  fs::dir_create(Path_Plot)

  IASDT.R::CatTime("Loading species list")
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

    IASDT.R::CatTime("Either ModelEval or VarPar was not provided and has to be calculated")

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Loading model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (magrittr::not(inherits(Model, "Hmsc"))) {

      IASDT.R::CatTime("  >>  Loading model object from file")

      if (inherits(Model, "character")) {
        if (magrittr::not(file.exists(Model))) {
          MSG <- "The provided path for the model does not exist"
          stop(MSG)
        }
        Model <- IASDT.R::LoadAs(Model)
        invisible(gc())

        if (magrittr::not(inherits(Model, "Hmsc"))) {
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
      IASDT.R::CatTime("  >>  Calculate variance partitioning")
      VarPar <- Hmsc::computeVariancePartitioning(Model)

      if (SaveVarPar) {
        save(VarPar,
             file = file.path(Path_Plot, "VariancePartitioning_DT.RData"))
      }
      invisible(gc())
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Compute predicted Values
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (is.null(ModelEval)) {

      IASDT.R::CatTime("  >>  Compute R2")

      IASDT.R::CatTime("  >>  >>  Compute predicted Values")
      # 06.07.2024 - This uses the updated predict function, currently available on my forked version of the package github.com/elgabbas/Hmsc
      # The `Hmsc::evaluateModelFit` function expects an array object returned from `Hmsc::computePredictedValues`. The `computePredictedValues` function does not work on parallel, so I used the updated predict function on parallel then converted the out put to array

      preds <- IASDT.R::Mod_Pred2Array(
        Predict = TRUE, Model = Model, NCores = NCores)

      IASDT.R::CatTime("  >>  >>  Evaluate model fit")
      ModelEval <- Hmsc::evaluateModelFit(hM = Model, predY = preds) %>%
        suppressWarnings()
      if (SaveModelEval) {
        save(ModelEval, file = file.path(Path_Plot, "ModelEval_explanatory.RData"))
      }

      rm(preds)
    }
  }

  rm(Model)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plot theme
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Theme <- ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.position = "none",
    axis.text.y = ggplot2::element_text(face = "italic"),
    axis.title.y = ggplot2::element_text(size = 25),
    axis.title.x = ggplot2::element_text(size = 25),
    legend.text = ggplot2::element_text(size = 28),
    legend.key.spacing.x = ggplot2::unit(0.7, "cm"))

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

  IASDT.R::CatTime("Plot 1 - Relative variance partitioning")

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

  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Raw variance partitioning
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plot 2 - Raw variance partitioning")

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

  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Combine plots
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Combine plots")

  Legend <- Plot +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 10)))
  Legend <- g_legend(Legend)

  if (is.null(PlotTitlePrefix)) {
    PlotTitle <- "Variance partitioning"
  } else {
    PlotTitle <- paste0("Variance partitioning (", PlotTitlePrefix, ")")
  }

  if (SaveGG) {
    IASDT.R::CatTime("Save plots as RData")
    PlotList <- list(Plot = Plot, Plot_raw = Plot_raw)
    save(PlotList, file = file.path(Path_Plot, "VariancePartitioning_GG.RData"))
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Save plot to disk
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save plot to disk")

  PlotPath2 <- file.path(Path_Plot, "VariancePartitioning.jpeg")
  VarParPlot <- gridExtra::arrangeGrob(Plot, Plot_raw, nrow = 1) %>%
    gridExtra::grid.arrange(
      Legend, nrow = 2, heights = c(15, 1),
      top = grid::textGrob(
        PlotTitle, gp = grid::gpar(fontsize = 30, font = 2)))
  ggplot2::ggsave(
    plot = VarParPlot, filename = PlotPath2, height = 16, width = 24, dpi = 600)

  if (ReturnGG) {
    return(VarParPlot)
  } else {
    return(invisible(NULL))
  }
}
