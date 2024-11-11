## |------------------------------------------------------------------------| #
# VarPar_Plot ----
## |------------------------------------------------------------------------| #

#' Plot variance partitioning of Hmsc model using ggplot2.
#'
#' This function generates and optionally saves plots visualizing the variance
#' partitioning of a Hmsc model. It can calculate variance partitioning and
#' model evaluation if not provided and supports parallel computation for model
#' predictions.
#' @param Path_Model Character path for the model file.
#' @param NCores Integer. Number of parallel computations for computing variance
#'   partitioning using TensorFlow. See [VarPar_Compute] for more details.
#'   Default: `1`.
#' @param EnvFile String. Path to read the environment variables. Default value:
#'   `.env`
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
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
    TF_Environ = NULL, NCores = 1) {
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Model) || is.null(NCores)) {
    stop("`Path_Model` and `NCores` cannot be empty", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- Variable <- VP_Value <-
    TaxaInfoFile <- Sp <- TjurR2 <- Label <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs, function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("EnvFile")
  )
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")
  rm(AllArgs)

  # # ..................................................................... ###

  # Loading species list -----

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

  # # ..................................................................... ###

  # Loading species info -----

  SpList <- IASDT.R::LoadAs(TaxaInfoFile) %>%
    dplyr::select(Species = IAS_ID, Species_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Species = stringr::str_pad(string = Species, width = 4, pad = "0"),
      Species = paste0("Sp_", Species))

  # # ..................................................................... ###

  # Loading model evaluation ----

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

  Model_Eval <- Path_Eval %>%
    IASDT.R::LoadAs() %>%
    # filter out the species that are not in the model
    dplyr::filter(stringr::str_starts(IAS_ID, "Sp_")) %>%
    dplyr::rename(Species = IAS_ID) %>%
    dplyr::select(-Sp)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Calculate variance partitioning ----

  File_VarPar <- file.path(Path_VarPar, "VarPar_DT.RData")
  File_VarPar_all <- file.path(Path_VarPar, "VarPar_DT_all.RData")

  if (!file.exists(File_VarPar) || !file.exists(File_VarPar_all)) {

    if (!file.exists(Path_Model)) {
      stop("The provided path for the model does not exist", call. = FALSE)
    }
    Model <- IASDT.R::LoadAs(Path_Model)

    if (!inherits(Model, "Hmsc")) {
      stop("The loaded model is not of an Hmsc object", call. = FALSE)
    }
  }

  if (!file.exists(File_VarPar)) {
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
      Verbose = FALSE, OutFileName = "VarPar_DT")

    rm(ModelVars, groupnames, group)

    } else {

    VarPar <- IASDT.R::LoadAs(File_VarPar)

    }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  if (!file.exists(File_VarPar)) {
    groupnames_all <- dimnames(Model$X)[[2]][-1] %>%
      stringr::str_remove_all("stats::poly\\(") %>%
      stringr::str_replace_all(", degree = 2, raw = TRUE\\)", "_")
    group_all <- c(1, seq_along(groupnames_all))

    VarPar_all <- IASDT.R::VarPar_Compute(
      Path_Model = Path_Model, group = group_all, groupnames = groupnames_all,
      NCores = NCores, UseTF = UseTF, TF_Environ = TF_Environ,
      Verbose = FALSE, OutFileName = "VarPar_DT_all")

    rm(groupnames_all, group_all)

  } else {

    VarPar_all <- IASDT.R::LoadAs(File_VarPar_all)

  }

  if (exists("Model")) {
    rm(Model)
  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plot theme ----

  Theme <- ggplot2::theme(
    plot.title = ggtext::element_markdown(size = 14, hjust = 0.5),
    plot.title.position = "plot",
    legend.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0.1, 0.3, 0.2, 0.05, "cm"),
    axis.text.x = ggplot2::element_text(
      face = "italic", size = 5, angle = 90, hjust = 1, vjust = 0.3),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.text = ggtext::element_markdown(
      size = 9, hjust = 0.5, lineheight = 1.15),
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.key.height = ggplot2::unit(0.4, "cm"),
    legend.key.width = ggplot2::unit(0.6, "cm"),
    legend.key.spacing.x = ggplot2::unit(0.4, "cm"),
    legend.position = "bottom",
    legend.margin = ggplot2::margin(0, 2, 0, 1, unit = "cm"))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Grouped data and plots ------

  ## Plotting data ----

  ### Relative variance partitioning ----
  VarPar_DF <- tibble::as_tibble(VarPar$vals, rownames = "Variable") %>%
    tidyr::pivot_longer(
      cols = -Variable, names_to = "Species", values_to = "VP_Value")

  # Calculate the mean variance partitioning per variable
  VarPar_Mean <- VarPar_DF %>%
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
      Label = stringr::str_to_title(Label),
      Label = paste0(
        "<b>", as.character(Label), "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Order variables by the mean variance partitioning
  VarOrder <- dplyr::pull(VarPar_Mean, Variable)

  # Order species by the mean variance partitioning per Variable
  SpOrder <- dplyr::left_join(VarPar_DF, SpList, by = "Species") %>%
    dplyr::mutate(Variable = factor(Variable, VarOrder)) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c(Species_name, Variable)) %>%
    dplyr::arrange(Variable, dplyr::desc(VP_Value)) %>%
    dplyr::pull(Species_name) %>%
    unique()

  # Plotting data for raw variance partitioning - ordered by mean variance
  # partitioning per Variable
  DT_Relative <- dplyr::left_join(VarPar_DF, SpList, by = "Species") %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder)) %>%
    dplyr::left_join(VarPar_Mean, by = "Variable")

  # Plotting data for raw variance partitioning - original species order
  SpOrder_Orig <- dplyr::left_join(VarPar_DF, SpList, by = "Species") %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name) %>%
    unique()

  DT_Relative_Orig <- DT_Relative %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Orig))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ### Raw variance partitioning ----


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
  MaxVal <- dplyr::group_by(VarPar_DF_Raw, Species) %>%
    dplyr::summarize(VP_Value = sum(VP_Value)) %>%
    dplyr::pull(VP_Value) %>%
    max() %>%
    magrittr::multiply_by(100) %>%
    ceiling(x = .) %>%
    magrittr::divide_by(100)

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
      Label = stringr::str_to_title(Label),
      Label = paste0(
        "<b>", as.character(Label), "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Plotting data for relative variance partitioning - ordered by mean variance
  # partitioning per Variable
  DT_Raw <- dplyr::left_join(VarPar_DF_Raw, SpList, by = "Species") %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder)) %>%
    dplyr::left_join(VarPar_Raw_Mean, by = "Variable")

  # Plotting data for relative variance partitioning - original species order
  SpOrder_Raw_Orig <- dplyr::left_join(
    VarPar_DF_Raw, SpList, by = "Species") %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name) %>%
    unique()

  DT_Raw_Orig <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Raw_Orig))

  # # ..................................................................... ###

  ## Plotting ----

  ### Relative variance partitioning plot -----

  # ordered by mean variance partitioning
  Title_Relative <- c(
    "<b>Variance partitioning</b> (sorted by mean value per predictor)")

  Plot_Relative <- DT_Relative %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Relative) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_Grouped_Relative_Ordered.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Relative)
  grDevices::dev.off()

  # ordered by original species order
  Title_Relative_Orig <- paste0(
    "<b>Variance partitioning</b> (original species order)")

  Plot_Relative_Orig <- DT_Relative_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Relative_Orig) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  grDevices::jpeg(
    filename = file.path(
      Path_VarPar, "VarPar_Grouped_Relative_Original.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Relative_Orig)
  grDevices::dev.off()

  # # ..................................................................... ###

  ### Raw variance partitioning plot ----

  ## ordered by mean variance partitioning
  Title_Raw <- paste0(
    "<b>Raw variance partitioning</b> ",
    "(relative variance partitioning  &times; Tjur-R<sup>2</sup>)")

  Plot_Raw <- DT_Raw_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, MaxVal), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Raw) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(
      Path_VarPar, "VarPar_Grouped_Raw_Ordered.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Raw)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###
  # # ..................................................................... ###

  # Ungrouped data and plots -----

  ## Plotting data ----

  ### Relative variance partitioning ----
  VarPar_all_DF <- tibble::as_tibble(
    VarPar_all$vals, rownames = "Variable") %>%
    tidyr::pivot_longer(
      cols = -Variable, names_to = "Species", values_to = "VP_Value")

  # Calculate the mean variance partitioning per variable
  VarPar_all_Mean <- VarPar_all_DF %>%
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
      Label = stringr::str_to_title(Label),
      Label = paste0(
        "<b>", as.character(Label), "</b><br>(", round(VP_Value, 1), "%)"),
      Label = stringr::str_replace_all(Label, "_1", ""),
      Label = stringr::str_replace_all(Label, "_2", "^2"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Order variables by the mean variance partitioning
  VarOrder_all <- dplyr::pull(VarPar_all_Mean, Variable)

  # Order species by the mean variance partitioning per Variable
  SpOrder_all <- dplyr::left_join(VarPar_all_DF, SpList, by = "Species") %>%
    dplyr::mutate(Variable = factor(Variable, VarOrder_all)) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c(Species_name, Variable)) %>%
    dplyr::arrange(Variable, dplyr::desc(VP_Value)) %>%
    dplyr::pull(Species_name) %>%
    unique()

  # Plotting data for raw variance partitioning - ordered by mean variance
  # partitioning per Variable
  DT_Relative_all <- dplyr::left_join(
    VarPar_all_DF, SpList, by = "Species") %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_all)) %>%
    dplyr::left_join(VarPar_all_Mean, by = "Variable")

  # Plotting data for raw variance partitioning - original species order
  SpOrder_all_Orig <- VarPar_all_DF %>%
    dplyr::left_join(SpList, by = "Species") %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name) %>%
    unique()

  DT_Relative_all_Orig <- DT_Relative_all %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_all_Orig))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ### Raw variance partitioning ----

  if (!is.null(Model_Eval$TjurR2) &&
      length(Model_Eval$TjurR2) == ncol(VarPar_all$vals)) {

    VarPar_all_DF_Raw <- tibble::as_tibble(
      VarPar_all$vals, rownames = "Variable") %>%
      tidyr::pivot_longer(
        cols = -Variable, names_to = "Species", values_to = "VP_Value") %>%
      dplyr::left_join(Model_Eval, by = "Species") %>%
      dplyr::mutate(VP_Value = VP_Value * TjurR2)

  } else {

    stop(
      paste0(
        "Mismatch between the length of Model_Eval$TjurR2 and the number of ",
        " columns in VarPar_all$vals"),
      call. = FALSE)

  }


  # Calculate the maximum value for the x-axis
  MaxVal <- dplyr::group_by(VarPar_all_DF_Raw, Species) %>%
    dplyr::summarize(VP_Value = sum(VP_Value)) %>%
    dplyr::pull(VP_Value) %>%
    max() %>%
    magrittr::multiply_by(100) %>%
    ceiling(x = .) %>%
    magrittr::divide_by(100)

  VarPar_all_Raw_Mean <- VarPar_all_DF_Raw %>%
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
      Label = stringr::str_to_title(Label),
      Label = paste0(
        "<b>", as.character(Label), "</b><br>(", round(VP_Value, 1), "%)"),
      Label = stringr::str_replace_all(Label, "_1", ""),
      Label = stringr::str_replace_all(Label, "_2", "^2"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Plotting data for relative variance partitioning - ordered by mean variance
  # partitioning per Variable
  DT_Raw_all <- dplyr::left_join(VarPar_all_DF_Raw, SpList, by = "Species") %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_all)) %>%
    dplyr::left_join(VarPar_all_Raw_Mean, by = "Variable")

  # Plotting data for relative variance partitioning - original species order
  SpOrder_all_Raw_Orig <- dplyr::left_join(
    VarPar_all_DF_Raw, SpList, by = "Species") %>%
    dplyr::distinct(Species, Species_name) %>%
    dplyr::arrange(Species) %>%
    dplyr::pull(Species_name) %>%
    unique()

  DT_Raw_all_Orig <- DT_Raw_all %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_all_Raw_Orig))

  # # ..................................................................... ###

  ## Plotting ----

  ### Relative variance partitioning plot -----

  # ordered by mean variance partitioning
  Title_Relative_all <- paste0(
    "<b>Variance partitioning</b> (sorted by mean value per predictor)")

  Colours <- paletteer::paletteer_d(
    "ggthemes::Tableau_20", n = nrow(VarPar_all_Raw_Mean))

  Plot_Relative <- DT_Relative_all %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_manual(values = Colours) +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Relative_all) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_all_Relative_Ordered.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Relative)
  grDevices::dev.off()

  # ordered by original species order
  Title_Relative_all_Orig <- paste0(
    "<b>Variance partitioning</b> (original species order)")

  Plot_Relative_Orig <- DT_Relative_all_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_manual(values = Colours) +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Relative_all_Orig) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_all_Relative_Original.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Relative_Orig)
  grDevices::dev.off()

  # # ..................................................................... ###

  ### Raw variance partitioning plot ----

  ## ordered by mean variance partitioning
  Title_Raw <- paste0(
    "<b>Raw variance partitioning</b> ",
    "(relative variance partitioning  &times; Tjur-R<sup>2</sup>)")

  Plot_Raw <- DT_Raw_all_Orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_manual(values = Colours) +
    ggplot2::coord_flip(xlim = c(0, MaxVal), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    ggplot2::labs(title = Title_Raw) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::jpeg(
    filename = file.path(Path_VarPar, "VarPar_all_Raw_Ordered.jpeg"),
    width = 30, height = 20, units = "cm", quality = 100, res = 600)
  plot(Plot_Raw)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Save plot to disk -----

  IASDT.R::CatDiff(InitTime = .StartTime)
  return(invisible(NULL))
}


# setwd("D:/BioDT_IAS/")
#
# VarPar_Plot(
#     Path_Model = "Z:/datasets/processed/model_fitting/DE_SW_CV_12b/Model_Fitted/GPP25_Tree_samp1000_th100_Model.RData",
#     EnvFile = ".env", FromHPC = FALSE, UseTF = TRUE,
#     TF_Environ = "D:/r-tensorflow/", NCores = 1)
