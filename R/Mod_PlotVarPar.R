## |------------------------------------------------------------------------| #
# PlotVarPar ----
## |------------------------------------------------------------------------| #

#' Plot variance partitioning of Hmsc model using ggplot2.
#'
#' This function generates and optionally saves plots visualizing the variance
#' partitioning of a Hmsc model. It can calculate variance partitioning and
#' model evaluation if not provided and supports parallel computation for model
#' predictions.
#' @param Model an object of class `Hmsc` or a path for the model file. Only
#'   needed if one of `ModelEval` and `VarPar` arguments are `NULL`.
#' @param Path_Plot String. Path where the output file(s) will be saved.
#' @param PlotTitlePrefix String (optional). Prefix to add to the title of the
#'   plot. Default: `NULL`, which means only 'Variance partitioning' will be
#'   used in the title.
#' @param ModelEval Result of the [Hmsc::evaluateModelFit] function. If
#'   `ModelEval = NULL` (default), [Hmsc::evaluateModelFit] will be executed on
#'   the model object to compute measures of model fit.
#' @param NCores Integer. Number of parallel computations for computing
#'   predicted values. This is used as the `nParallel` argument of the
#'   [Hmsc::computePredictedValues] function.
#' @param VarPar Variance partitioning object resulted from
#'   [Hmsc::computeVariancePartitioning]. If `VarPar = NULL` (default),
#'   [Hmsc::computeVariancePartitioning] will be executed on the model object.
#' @param EnvFile String. Path to read the environment variables. Default value:
#'   `.env`
#' @param SaveVarPar Logical. If `VarPar = NULL`, should the calculated variance
#'   partitioning be saved as an RData file?  Default: `TRUE`.
#' @param SaveModelEval Logical. If `ModelEval = NULL`, should the calculated
#'   model evaluation be saved as an RData file?  Default: `TRUE`.
#' @param ReturnGG Logical. Return the ggplot object. Default: `FALSE`.
#' @param SaveGG Logical. Should the plots be exported as an RData object?
#'   Default: `TRUE`.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @name PlotVarPar
#' @author Ahmed El-Gabbas
#' @details The function reads the following environment variables:
#'   - **`DP_R_TaxaInfo_RData`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_TaxaInfo_RData_Local`** (if `FromHPC` = `FALSE`) for the location of the `TaxaList.RData` file containing species information.
#' @return If `ReturnGG` is `TRUE`, returns a ggplot object of the variance
#'   partitioning plots. Otherwise, returns `NULL`.
#' @export

PlotVarPar <- function(
    Model, Path_Plot, PlotTitlePrefix = NULL, ModelEval = NULL, NCores,
    VarPar = NULL, EnvFile = ".env", SaveVarPar = TRUE, SaveModelEval = TRUE,
    ReturnGG = FALSE, SaveGG = TRUE, FromHPC = TRUE) {

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Model) || is.null(Path_Plot) || is.null(NCores)) {
    stop("Model, Path_Plot, and NCores cannot be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- variable <- value <- TaxaInfoFile <- NULL

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

  if (magrittr::not(file.exists(EnvFile))) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo_RData", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo_RData_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  SpList <- IASDT.R::LoadAs(TaxaInfoFile) %>%
    dplyr::select(Species = IAS_ID, Species_name) %>%
    dplyr::mutate(
      Species = stringr::str_pad(string = Species, width = 4, pad = "0"),
      Species = paste0("Sp_", Species))


  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Processing ModelEval and VarPar
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (purrr::some(list(ModelEval, VarPar), is.null)) {
    if (is.null(ModelEval)) {
      IASDT.R::CatTime(
        "ModelEval was not provided and has to be calculated")
    }
    if (is.null(VarPar)) {
      IASDT.R::CatTime(
        "VarPar was not provided and has to be calculated")
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Loading model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (magrittr::not(inherits(Model, "Hmsc"))) {

      IASDT.R::CatTime("  >>  Loading model object from file")

      if (inherits(Model, "character")) {
        if (magrittr::not(file.exists(Model))) {
          stop("The provided path for the model does not exist")
        }
        Model <- IASDT.R::LoadAs(Model)

        if (magrittr::not(inherits(Model, "Hmsc"))) {
          stop("The loaded model is not of an Hmsc object")
        }

      } else {
        stop("The Model has to be an Hmsc model or a path to a saved model")
      }
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Calculate variance partitioning
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (is.null(VarPar)) {
      IASDT.R::CatTime("  >>  Calculate variance partitioning")
      VarPar <- Hmsc::computeVariancePartitioning(hM = Model)

      if (SaveVarPar) {
        save(VarPar,
             file = file.path(Path_Plot, "VariancePartitioning_DT.RData"))
      }
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Compute predicted Values
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (is.null(ModelEval)) {

      IASDT.R::CatTime("  >>  Compute R2")

      IASDT.R::CatTime("  >>  >>  Compute predicted Values")
      # 06.07.2024 - This uses the updated predict function, currently available
      # on my forked version of the package github.com/elgabbas/Hmsc The
      # Hmsc::evaluateModelFit function expects an array object returned from
      # Hmsc::computePredictedValues. The `computePredictedValues` function does
      # not work on parallel, so I used the updated predict function on parallel
      # then converted the output to array

      preds <- IASDT.R::Mod_Pred2Array(
        Predict = TRUE, Model = Model, expected = TRUE, NCores = NCores)

      IASDT.R::CatTime("  >>  >>  Evaluate model fit")
      # Suppress the warning: In cor(lbeta[[i]][k, ], lmu[[i]][k, ]) : the
      # standard deviation is zero
      invisible(ModelEval <- suppressWarnings(Hmsc::evaluateModelFit(hM = Model, predY = preds)))

      if (SaveModelEval) {
        save(
          ModelEval, file = file.path(Path_Plot, "ModelEval_explanatory.RData"))
      }

      rm(preds)
    }
    rm(Model)
  }

  if (inherits(ModelEval, "character")) {
    IASDT.R::CatTime("Loading ModelEval from RData file")
    ModelEval <- IASDT.R::LoadAs(ModelEval)
  }

  if (inherits(VarPar, "character")) {
    IASDT.R::CatTime("Loading VarPar from RData file")
    VarPar <- IASDT.R::LoadAs(VarPar)
  }

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

  # Add a common Legend for combined ggplots
  # Source: https://stackoverflow.com/a/13650878/3652584
  g_legend <- function(Plot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(Plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Relative variance partitioning
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plotting relative variance partitioning")

  vp_df <- tibble::as_tibble(VarPar$vals, rownames = "variable") %>%
    tidyr::pivot_longer(
      cols = -variable, names_to = "Species", values_to = "value")

  VarOrder <- dplyr::group_by(vp_df, variable) %>%
    dplyr::summarise(value = mean(value)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::pull(variable)

  # Order species by the mean variance partitioning per variable
  SpOrder <- dplyr::left_join(vp_df, SpList, by = "Species") %>%
    dplyr::mutate(variable = factor(variable, VarOrder)) %>%
    dplyr::summarise(value = sum(value), .by = c(Species_name, variable)) %>%
    dplyr::arrange(variable, dplyr::desc(value)) %>%
    dplyr::pull(Species_name) %>%
    unique()

  Plot_DT <- dplyr::left_join(vp_df, SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable,
        "Habitat coverage" = "HabLog",
        "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder)) %>%
    # suppress warning: ! Unknown levels in `f`: HabLog; if habitat is not used
    # as predictor
    suppressWarnings()

  Plot <- ggplot2::ggplot(
    data = Plot_DT,
    mapping = ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, 1.0125), expand = FALSE) +
    Theme



  # Use original species order

  SpOrder_Orig <- dplyr::left_join(vp_df, SpList, by = "Species") %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name) %>%
    unique()

  Plot_DT_Orig <- dplyr::left_join(vp_df, SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable,
        "Habitat coverage" = "HabLog",
        "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder_Orig)) %>%
    # suppress warning: ! Unknown levels in `f`: HabLog; if habitat is not used
    # as predictor
    suppressWarnings()

  Plot_Orig <- ggplot2::ggplot(
    data = Plot_DT_Orig,
    mapping = ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species (original order)") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, 1.0125), expand = FALSE) +
    Theme

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Raw variance partitioning
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Plotting raw variance partitioning")

  if (!is.null(ModelEval$TjurR2) &&
      length(ModelEval$TjurR2) == ncol(VarPar$vals)) {
    for (ii in seq_along(ModelEval$TjurR2)) {
      VarPar$vals[, ii] <- ModelEval$TjurR2[ii] * VarPar$vals[, ii]
    }
  } else {
    stop(paste0("Mismatch between the length of ModelEval$TjurR2 and the ",
                "number of columns in VarPar$vals"))
  }

  vp_df_R <- tibble::as_tibble(VarPar$vals, rownames = "variable") %>%
    tidyr::pivot_longer(
      cols = -variable, names_to = "Species", values_to = "value")

  MaxVal <- dplyr::group_by(vp_df_R, Species) %>%
    dplyr::summarize(value = sum(value)) %>%
    dplyr::pull(value) %>%
    max() %>%
    magrittr::multiply_by(100) %>%
    ceiling(x = .) %>%
    magrittr::divide_by(100)

  # Order species by the mean variance partitioning per variable

  Plot_R_DT <- dplyr::left_join(vp_df_R, SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable,
        "Habitat coverage" = "HabLog",
        "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder)) %>%
    # suppress warning: ! Unknown levels in `f`: HabLog; if habitat is not used
    # as predictor
    suppressWarnings()

  Plot_raw <- ggplot2::ggplot(
    data = Plot_R_DT,
    ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Variance explained (raw)") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, MaxVal), expand = FALSE) +
    Theme

  # Use original species order

  Plot_R_DT_Orig <- dplyr::left_join(vp_df_R, SpList, by = "Species") %>%
    dplyr::arrange(variable, value) %>%
    dplyr::mutate(
      variable = factor(variable, VarOrder),
      variable = forcats::fct_recode(
        variable,
        "Habitat coverage" = "HabLog",
        "Road + Rail intensity" = "RoadRailLog",
        "Sampling intensity" = "BiasLog",
        "Spatial random effect" = "Random: sample"),
      Species_name = factor(Species_name, SpOrder_Orig)) %>%
    # suppress warning: ! Unknown levels in `f`: HabLog; if habitat is not used
    # as predictor
    suppressWarnings()

  Plot_raw_Orig <- ggplot2::ggplot(
    data = Plot_R_DT_Orig,
    ggplot2::aes(x = value, y = Species_name, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species (original order)") +
    ggplot2::xlab("Variance explained (raw)") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_cartesian(xlim = c(0, MaxVal), expand = FALSE) +
    Theme

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
    PlotList <- list(
      Plot = Plot, Plot_raw = Plot_raw,
      Plot_Orig = Plot_Orig, Plot_raw_Orig = Plot_raw_Orig)
    save(PlotList, file = file.path(Path_Plot, "VariancePartitioning_GG.RData"))
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Save plot to disk
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save plot to disk")

  PlotPath2 <- file.path(Path_Plot, "VariancePartitioning.jpeg")
  VarParPlot <- gridExtra::arrangeGrob(Plot, Plot_raw, nrow = 1) %>%
    gridExtra::arrangeGrob(
      Legend, nrow = 2, heights = c(15, 1),
      top = grid::textGrob(
        PlotTitle, gp = grid::gpar(fontsize = 30, font = 2)))
  ggplot2::ggsave(
    plot = VarParPlot, filename = PlotPath2, height = 16, width = 24, dpi = 600)

  PlotPath2 <- file.path(Path_Plot, "VariancePartitioning_Orig.jpeg")
  VarParPlot <- gridExtra::arrangeGrob(Plot_Orig, Plot_raw_Orig, nrow = 1) %>%
    gridExtra::arrangeGrob(
      Legend, nrow = 2, heights = c(15, 1),
      top = grid::textGrob(
        PlotTitle, gp = grid::gpar(fontsize = 30, font = 2)))
  ggplot2::ggsave(
    plot = VarParPlot, filename = PlotPath2, height = 16, width = 24, dpi = 600)

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  if (ReturnGG) {
    return(VarParPlot)
  } else {
    return(invisible(NULL))
  }
}
