## |------------------------------------------------------------------------| #
# VarPar_Plot ----
## |------------------------------------------------------------------------| #

#' Plot variance partitioning of Hmsc model
#'
#' This function generates plots for the variance partitioning of an Hmsc model.
#' It can optionally compute the variance partitioning if not available on disk,
#' supporting parallel computation and using TensorFlow; see [VarPar_Compute].
#' It plots the relative variance partitioning sorted by the mean value per
#' predictor or by the original species order. It also plots the raw variance
#' partitioning (relative variance partitioning multiplied by the Tjur-R^2
#' value). The plots are saved as JPEG files.
#'
#' @param Path_Model Character path for the model file.
#' @param NCores Integer. Number of parallel computations for computing variance
#'   partitioning using TensorFlow. See [VarPar_Compute] for more details.
#'   Default: `1`.
#' @param EnvFile String. Path to read the environment variables. Default value:
#'   `.env`
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param Fig_width,Fig_height Numeric. Width and height of the output plot in
#'   centimeters. Default: `30` and `15`, respectively.
#' @name VarPar_Plot
#' @author Ahmed El-Gabbas
#' @details The function reads the following environment variables:
#'   - **`DP_R_TaxaInfo_RData`** (if `FromHPC` = `TRUE`) or
#'     **`DP_R_TaxaInfo_RData_Local`** (if `FromHPC` = `FALSE`) for the
#'   location of the `TaxaList.RData` file containing species information.
#' @inheritParams Predict_Hmsc
#' @export

VarPar_Plot <- function(
    Path_Model, EnvFile = ".env", FromHPC = TRUE, UseTF = TRUE,
    TF_Environ = NULL, NCores = 1, Fig_width = 30, Fig_height = 15) {

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- Variable <- VP_Value <-
    TaxaInfoFile <- Sp <- TjurR2 <- Label <- VP_Sum <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ------

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs, function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("EnvFile", "Path_Model"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical", Args = c("FromHPC", "UseTF"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Species info -----

  IASDT.R::CatTime("Loading species info")

  if (!file.exists(EnvFile)) {
    stop(
      paste0("Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
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
    dplyr::distinct() %>%
    dplyr::mutate(
      Species = stringr::str_pad(string = Species, width = 4, pad = "0"),
      Species = paste0("Sp_", Species))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Model evaluation ----

  IASDT.R::CatTime("Loading model evaluation")

  Path_Root <- dirname(dirname(Path_Model))
  Path_VarPar <- file.path(
    Path_Root, "Model_Postprocessing/Variance_Partitioning")
  Path_Eval <- file.path(Path_Root, "Model_Evaluation") %>%
    list.files("Eval_.+.qs", full.names = TRUE)

  if (length(Path_Eval) != 1) {
    stop(
      paste0(
        "The number of model evaluation files in the directory: ", Path_Root,
        " is not equal to 1"),
      call. = FALSE)
  }

  Model_Eval <- IASDT.R::LoadAs(Path_Eval) %>%
    # filter out the species that are not in the model
    dplyr::filter(stringr::str_starts(IAS_ID, "Sp_")) %>%
    dplyr::rename(Species = IAS_ID) %>%
    dplyr::select(-Sp) %>%
    dplyr::left_join(SpList, by = "Species")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Compute variance partitioning ----

  IASDT.R::CatTime("Compute/load variance partitioning")

  File_VarPar <- file.path(Path_VarPar, "VarPar.RData")

  if (!file.exists(File_VarPar)) {

    IASDT.R::CatTime("Variance partitioning will be calculated", Level = 1)

    IASDT.R::CatTime("Loading fitted model", Level = 1)
    if (!file.exists(Path_Model)) {
      stop("The provided path for the model does not exist", call. = FALSE)
    }
    Model <- IASDT.R::LoadAs(Path_Model)

    if (!inherits(Model, "Hmsc")) {
      stop("The loaded model is not of an Hmsc object", call. = FALSE)
    }

    # Keep only selected list items
    Model <- Model[c("XData", "X")]
    invisible(gc())

  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  if (!file.exists(File_VarPar)) {

    IASDT.R::CatTime("Computing variance partitioning", Level = 1)

    # names of variables used in the model
    groupnames <- names(Model$XData)

    # actual variables used in the model, including the intercept and quadratic
    # terms
    ModelVars <- dimnames(Model$X)[[2]][-1]

    # group variable to combine variable and its quadratic terms together
    group <- purrr::map(
      ModelVars, ~ which(stringr::str_detect(.x, groupnames))) %>%
      unlist() %>%
      # add intercept to the first group
      c(.[1], .)

    VarPar <- IASDT.R::VarPar_Compute(
      Path_Model = Path_Model, group = group, groupnames = groupnames,
      NCores = NCores, UseTF = UseTF, TF_Environ = TF_Environ,
      Verbose = FALSE, OutFileName = "VarPar")

    rm(ModelVars, groupnames, group, envir = environment())

  } else {

    IASDT.R::CatTime("Loading variance partitioning", Level = 2)
    VarPar <- IASDT.R::LoadAs(File_VarPar)

  }

  if (exists("Model")) {
    rm(Model, envir = environment())
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plot theme ----

  Theme <- ggplot2::theme(
    plot.title = ggtext::element_markdown(
      size = 14, hjust = 0.5, vjust = 0),
    plot.title.position = "plot",
    legend.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0.2, 0.15, 0.2, 0.05, "cm"),
    axis.text.x = ggplot2::element_text(
      face = "italic", size = 5, angle = 90, hjust = 1, vjust = 0.3),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.text = ggtext::element_markdown(
      size = 8, hjust = 0.5, lineheight = 1.15),
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.key.height = ggplot2::unit(0.4, "cm"),
    legend.key.width = ggplot2::unit(0.6, "cm"),
    legend.key.spacing.x = ggplot2::unit(0.4, "cm"),
    legend.position = "bottom",
    legend.margin = ggplot2::margin(0, 2, 0, 1, unit = "cm"),
    legend.box.spacing = ggplot2::unit(1, "pt"))

  # # ..................................................................... ###
  # # ..................................................................... ###

  IASDT.R::CatTime("Prepare plotting data")

  # Relative variance partitioning ----

  ## Plotting data ----

  VarPar_DF <- tibble::as_tibble(VarPar$vals, rownames = "Variable") %>%
    tidyr::pivot_longer(
      cols = -Variable, names_to = "Species", values_to = "VP_Value") %>%
    dplyr::left_join(SpList, by = "Species")

  # Calculate mean VarPar per variable and prepare labels for the plot
  VarPar_Mean <- VarPar_DF %>%
    dplyr::summarise(VP_Value = mean(VP_Value), .by = "Variable") %>%
    dplyr::arrange(dplyr::desc(VP_Value)) %>%
    dplyr::mutate(
      VP_Value = VP_Value * 100,
      Label = forcats::fct_recode(
        Variable,
        "Habitat coverage" = "HabLog",
        "Road+Rail intensity" = "RoadRailLog",
        "Sampling efforts" = "EffortsLog",
        "Spatial effect" = "Random: sample"),
      Label = paste0(
        "<b>", stringr::str_to_title(as.character(Label)),
        "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Order variables by the mean variance partitioning
  VarOrder <- dplyr::pull(VarPar_Mean, Variable)

  # Order species by the mean variance partitioning per Variable
  SpOrder <- VarPar_DF %>%
    dplyr::mutate(Variable = factor(Variable, VarOrder)) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c(Species_name, Variable)) %>%
    dplyr::arrange(Variable, dplyr::desc(VP_Value)) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name)

  # Order species by the mean variance partitioning per Variable
  SpOrder_TjurR2<- Model_Eval %>%
    dplyr::arrange(dplyr::desc(TjurR2)) %>%
    dplyr::pull(Species_name)


  # Plotting data for variance partitioning

  # 1. ordered by mean variance partitioning per Variable
  DT_Relative <- VarPar_DF %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder)) %>%
    dplyr::left_join(VarPar_Mean, by = "Variable")

  # 2. original species order
  SpOrder_Orig <- VarPar_DF %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name)
  DT_Relative_Orig <- DT_Relative %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Orig))

  # 3. ordered by Tjur-R2
  DT_Relative_TjurR2 <- DT_Relative %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_TjurR2))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----

  IASDT.R::CatTime("Plotting")

  # 1. ordered by mean variance partitioning

  Title_Relative <- paste0(
    "<b>Proportion of explained variance</b> ",
    '<span style="font-size:12pt">',
    '(species sorted by mean value per predictor)</span>')

  Plot_Relative <- DT_Relative %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Relative) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Relative_Ordered.jpeg"),
    width = Fig_width, height = Fig_height,
    units = "cm", quality = 100, res = 600)
  plot(Plot_Relative)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 2. ordered by original species order
  Title_Relative_Orig <- paste0(
    "<b>Proportion of explained variance</b> ",
    '<span style="font-size:12pt">',
    '(species sorted by species taxonomy)</span>')

  Plot_Relative_Orig <- DT_Relative_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Relative_Orig) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Relative_Original.jpeg"),
    width = Fig_width, height = Fig_height,
    units = "cm", quality = 100, res = 600)
  plot(Plot_Relative_Orig)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 3. ordered by Tjur-R2

  Title_Relative_TjurR2 <- paste0(
    "<b>Proportion of explained variance</b> ",
    '<span style="font-size:12pt">',
    '(species sorted by Tjur-R<sup>2</sup>)</span>')

  Plot_Relative_TjurR2 <- DT_Relative_TjurR2 %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Relative_TjurR2) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Relative_TjurR2.jpeg"),
    width = Fig_width, height = Fig_height,
    units = "cm", quality = 100, res = 600)
  plot(Plot_Relative_TjurR2)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Raw variance partitioning ----

  ## Plotting data ----

  if (!is.null(Model_Eval$TjurR2) &&
      length(Model_Eval$TjurR2) == ncol(VarPar$vals)) {

    VarPar_DF_Raw <- tibble::as_tibble(VarPar$vals, rownames = "Variable") %>%
      tidyr::pivot_longer(
        cols = -Variable, names_to = "Species", values_to = "VP_Value") %>%
      dplyr::left_join(Model_Eval, by = "Species") %>%
      dplyr::mutate(VP_Value = VP_Value * TjurR2)

  } else {

    stop(
      paste0(
        "Mismatch between the length of Model_Eval$TjurR2 and the number of ",
        " columns in VarPar$vals"),
      call. = FALSE)
  }

  # Calculate the maximum value for the x-axis
  # MaxVal <- dplyr::group_by(VarPar_DF_Raw, Species) %>%
  #   dplyr::summarize(VP_Value = sum(VP_Value)) %>%
  #   dplyr::pull(VP_Value) %>%
  #   max() %>%
  #   magrittr::multiply_by(100) %>%
  #   ceiling(x = .) %>%
  #   magrittr::divide_by(100)

  VarPar_Raw_Mean <- VarPar_DF_Raw %>%
    dplyr::summarise(VP_Value = mean(VP_Value), .by = "Variable") %>%
    dplyr::arrange(dplyr::desc(VP_Value)) %>%
    dplyr::mutate(VP_Value = VP_Value * 100) %>%
    dplyr::mutate(
      Label = forcats::fct_recode(
        Variable,
        "Habitat coverage" = "HabLog",
        "Road+Rail intensity" = "RoadRailLog",
        "Sampling efforts" = "EffortsLog",
        "Spatial effect" = "Random: sample"),
      Label = paste0(
        "<b>", stringr::str_to_title(as.character(Label)),
        "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Plotting data for relative variance partitioning - ordered by mean variance
  # partitioning per Variable
  DT_Raw <- VarPar_DF_Raw %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder)) %>%
    dplyr::left_join(VarPar_Raw_Mean, by = "Variable")

  # Plotting data for relative variance partitioning - original species order
  DT_Raw_Orig <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Orig))

  # Plotting data for relative variance partitioning - original species order
  SpOrder_TotalRaw <- DT_Raw %>%
    dplyr::summarize(
      VP_Sum = sum(VP_Value), .by = c(Species, Species_name)) %>%
    dplyr::arrange(dplyr::desc(VP_Sum)) %>%
    dplyr::mutate(Species_name = as.character(Species_name)) %>%
    dplyr::pull(Species_name)
  DT_Raw_TotalRaw <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_TotalRaw))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----

  ## 1. ordered by mean variance partitioning

  Title_Raw <- paste0(
    "<b>Proportion of raw variance</b> ",
    '<span style="font-size:12pt">',
    '(proportion of explained variance &times; Tjur-R<sup>2</sup>',
    ' --- species sorted by species taxonomy)</span>')

  Plot_Raw <- DT_Raw_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Raw) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Raw_Ordered.jpeg"),
    width = Fig_width, height = Fig_height,
    units = "cm", quality = 100, res = 600)
  plot(Plot_Raw)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## 2. ordered by _Total Variance Explained

  Title_Raw_TotalRaw <- paste0(
    "<b>Proportion of raw variance</b> ",
    '<span style="font-size:12pt">',
    '(proportion of explained variance &times; Tjur-R<sup>2</sup>',
    ' --- species sorted by total explained variance)</span>')

  Plot_Raw_TotalRaw <- DT_Raw_TotalRaw %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Raw_TotalRaw) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Raw_TotalVarExplained.jpeg"),
    width = Fig_width, height = Fig_height,
    units = "cm", quality = 100, res = 600)
  plot(Plot_Raw_TotalRaw)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Save plot to disk -----

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "Computing and plotting variance partitioning took ")

  return(invisible(NULL))
}
