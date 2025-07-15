## |------------------------------------------------------------------------| #
# variance_partitioning_plot ----
## |------------------------------------------------------------------------| #

#' @export
#' @name variance_partitioning
#' @rdname variance_partitioning
#' @order 2

variance_partitioning_plot <- function(
    path_model, env_file = ".env", VP_file = "VarPar", use_TF = TRUE,
    TF_environ = NULL, n_cores = 1L, width = 30, height = 15, Axis_text = 4) {

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  IAS_ID <- Species_name <- Species <- Variable <- VP_Value <-
    TaxaInfoFile <- Sp <- TjurR2 <- Label <- VP_Sum <- NULL

  # Set null device for `cairo`. This is to properly render the plots using
  # ggtext - https://github.com/wilkelab/cowplot/issues/73
  cowplot::set_null_device("cairo")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ------

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs, function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("env_file", "path_model"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "use_TF")
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_cores")
  rm(AllArgs, envir = environment())

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Species info -----

  ecokit::cat_time("Loading species info")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "TaxaInfoFile", "DP_R_Taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  SpList <- ecokit::load_as(TaxaInfoFile) %>%
    dplyr::select(Species = IAS_ID, Species_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Species = stringr::str_pad(string = Species, width = 4, pad = "0"),
      Species = paste0("Sp_", Species))

  # # ..................................................................... ###

  Path_Root <- dirname(dirname(path_model))
  Path_VarPar <- fs::path(
    Path_Root, "Model_Postprocessing", "Variance_Partitioning")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Model evaluation ----

  ecokit::cat_time("Loading model evaluation")

  Path_Eval <- fs::path(Path_Root, "Model_Evaluation") %>%
    list.files("Eval_.+.qs2", full.names = TRUE)

  if (length(Path_Eval) != 1) {
    ecokit::stop_ctx(
      paste0(
        "The number of model evaluation files in the directory: ",
        Path_Root, " is not equal to 1"),
      Path_Root = Path_Root, Path_Eval = Path_Eval,
      length_Path_Eval = length(Path_Eval), include_backtrace = TRUE)
  }

  Model_Eval <- ecokit::load_as(Path_Eval) %>%
    # filter out the species that are not in the model
    dplyr::filter(stringr::str_starts(IAS_ID, "Sp_")) %>%
    dplyr::rename(Species = IAS_ID) %>%
    dplyr::select(-Sp) %>%
    dplyr::left_join(SpList, by = "Species")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Compute or load variance partitioning ----

  ecokit::cat_time("Compute or load variance partitioning")

  if (is.null(VP_file)) {
    VP_file <- "VarPar"
  }

  File_VarPar <- fs::path(Path_VarPar, paste0(VP_file, ".RData"))

  if (file.exists(File_VarPar)) {

    ecokit::cat_time("Loading variance partitioning", level = 2L)
    VarPar <- ecokit::load_as(File_VarPar)

  } else {

    ecokit::cat_time(
      paste0(
        "Variance partitioning will be computed using ", n_cores, " cores ",
        dplyr::if_else(use_TF, "and", "without"), " `TensorFlow`."),
      level = 1L)

    VarPar <- IASDT.R::variance_partitioning_compute(
      path_model = path_model, n_cores = n_cores, use_TF = use_TF,
      TF_environ = TF_environ, verbose = TRUE, VP_file = VP_file)

  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plot theme ----

  Theme <- ggplot2::theme(
    plot.title = ggtext::element_markdown(size = 14, hjust = 0, vjust = 0),
    plot.subtitle = ggtext::element_markdown(size = 9, hjust = 0, vjust = 0),
    plot.title.position = "plot",
    legend.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0.2, 0.15, 0.2, 0.05, "cm"),
    axis.text.y = ggplot2::element_text(
      size = 7, margin = ggplot2::margin(r = 0)),
    axis.text.x = ggplot2::element_text(
      face = "italic", size = Axis_text, angle = 90, hjust = 1, vjust = 0.3,
      margin = ggplot2::margin(t = 0)),
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

  ecokit::cat_time("Prepare plotting data")

  # Relative variance partitioning ----

  ## Plotting data ----

  VarPar_DF <- tibble::as_tibble(VarPar$vals, rownames = "Variable") %>%
    tidyr::pivot_longer(
      cols = -Variable, names_to = "Species", values_to = "VP_Value") %>%
    dplyr::left_join(SpList, by = "Species")

  # Calculate mean Variance partitioning per variable and prepare labels for the
  # plot
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
        "River length" = "RiversLog",
        "Spatial effect" = "Random: sample"),
      Label = paste0(
        "<b>", stringr::str_to_title(as.character(Label)),
        "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Order variables by the mean variance partitioning
  VarOrder <- dplyr::pull(VarPar_Mean, Variable)

  # Order species by the mean variance partitioning per Variable
  SpOrder <- dplyr::mutate(VarPar_DF, Variable = factor(Variable, VarOrder)) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c(Species_name, Variable)) %>%
    dplyr::arrange(Variable, dplyr::desc(VP_Value)) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name)

  # Order species by total variance, excluding the spatial random effect
  SpOrder_NonSpatial <- VarPar_DF %>%
    dplyr::filter(!startsWith(Variable, "Random")) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c("Species", "Species_name")) %>%
    dplyr::arrange(dplyr::desc(VP_Value)) %>%
    dplyr::distinct(Species_name) %>%
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
  SpOrder_TjurR2 <- Model_Eval %>%
    dplyr::arrange(dplyr::desc(TjurR2)) %>%
    dplyr::pull(Species_name)
  DT_Relative_TjurR2 <- DT_Relative %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_TjurR2))

  # 4. ordered by total variance partitioning, excluding spatial random effect
  DT_Relative_NonSpatial <- DT_Relative %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_NonSpatial))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----

  ecokit::cat_time("Plotting")

  # 1. ordered by mean variance partitioning

  Title_Relative <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:10pt">Species are sorted by mean value \\
    per predictor</span>')

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
  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Relative_ByMean.jpeg"),
    width = width, height = height, units = "cm", quality = 100, res = 600)
  print(Plot_Relative)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 2. ordered by original species order
  Title_Relative_Orig <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:10pt">Species sorted by species \\
    taxonomy</span>')

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

  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Relative_ByTaxonomy.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(Plot_Relative_Orig)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 3. ordered by Tjur-R2

  Title_Relative_TjurR2 <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
      --- </span><span style="font-size:10pt">Species are sorted by \\
      explanatory power (Tjur-R<sup>2</sup>)</span>')

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

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Relative_ByTjurR2.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(Plot_Relative_TjurR2)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 4. ordered by total variance partitioning, excluding spatial random effect

  Title_NonSpatial <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:10pt">Species are sorted by total \\
    explained variance, excluding spatial random effect</span>')

  Plot_Relative_NonSpatial <- DT_Relative_NonSpatial %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_NonSpatial) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(
      Path_VarPar, "VarPar_Relative_ByTotalNonSpatial.jpeg"),
    width = width, height = height, units = "cm", quality = 100, res = 600)
  print(Plot_Relative_NonSpatial)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Raw variance partitioning ----

  ## Plotting data ----

  if (is.null(Model_Eval$TjurR2) ||
      length(Model_Eval$TjurR2) != ncol(VarPar$vals)) {
    ecokit::stop_ctx(
      paste0(
        "Mismatch between the length of Model_Eval$TjurR2 and the number of ",
        " columns in VarPar$vals"),
      Model_Eval_TjurR2 = Model_Eval$TjurR2,
      length_Model_Eval_TjurR2 = length(Model_Eval$TjurR2),
      ncol_VarPar_vals = ncol(VarPar$vals), include_backtrace = TRUE)
  }

  VarPar_DF_Raw <- tibble::as_tibble(VarPar$vals, rownames = "Variable") %>%
    tidyr::pivot_longer(
      cols = -Variable, names_to = "Species", values_to = "VP_Value") %>%
    dplyr::left_join(Model_Eval, by = "Species") %>%
    dplyr::mutate(VP_Value = VP_Value * TjurR2)

  VarPar_Raw_Mean <- VarPar_DF_Raw %>%
    dplyr::summarise(VP_Value = mean(VP_Value), .by = "Variable") %>%
    dplyr::arrange(dplyr::desc(VP_Value)) %>%
    dplyr::mutate(
      VP_Value = VP_Value * 100,
      Label = forcats::fct_recode(
        Variable,
        "Habitat coverage" = "HabLog",
        "Road+Rail intensity" = "RoadRailLog",
        "Sampling efforts" = "EffortsLog",
        "River length" = "RiversLog",
        "Spatial effect" = "Random: sample"),
      Label = paste0(
        "<b>", stringr::str_to_title(as.character(Label)),
        "</b><br>(", round(VP_Value, 1), "%)"),
      Label = factor(Label, .data$Label),
      VP_Value = NULL)

  # Plotting data for relative variance partitioning - ordered by mean
  # variance partitioning per Variable
  DT_Raw <- VarPar_DF_Raw %>%
    dplyr::arrange(Variable, VP_Value) %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder)) %>%
    dplyr::left_join(VarPar_Raw_Mean, by = "Variable") %>%
    dplyr::mutate(Species_name = as.character(Species_name))

  # Plotting data for relative variance partitioning - original species order
  DT_Raw_Orig <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Orig))

  # Plotting data for relative variance partitioning - original species order
  SpOrder_TotalRaw <- DT_Raw %>%
    dplyr::summarize(
      VP_Sum = sum(VP_Value), .by = c(Species, Species_name)) %>%
    dplyr::arrange(dplyr::desc(VP_Sum)) %>%
    dplyr::pull(Species_name)
  DT_Raw_TotalRaw <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_TotalRaw))

  # Order species by total variance, excluding the spatial random effect
  SpOrder_Raw_NonSpatial <- DT_Raw %>%
    dplyr::filter(!startsWith(Variable, "Random")) %>%
    dplyr::summarise(
      VP_Value = sum(VP_Value), .by = c("Species", "Species_name")) %>%
    dplyr::arrange(dplyr::desc(VP_Value)) %>%
    dplyr::distinct(Species_name) %>%
    dplyr::pull(Species_name)
  DT_Raw_NonSpatial <- DT_Raw %>%
    dplyr::mutate(Species_name = factor(Species_name, SpOrder_Raw_NonSpatial))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----

  ## 1. ordered by mean variance partitioning

  Title_Raw <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; Tjur-R<sup>2</sup>)</span> --- \\
    <span style="font-size:10pt">species are sorted by taxonomy</span>')

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

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Raw_ByTaxonomy.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")

  plot(Plot_Raw)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## 2. ordered by Total Variance Explained

  Title_Raw_TotalRaw <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; Tjur-R<sup>2</sup>)</span> --- \\
    <span style="font-size:10pt">Species are sorted by total explained \\
    variance</span>')

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

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Raw_ByMean.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(Plot_Raw_TotalRaw)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## 3. by total variance partitioning, excluding spatial random effect

  Title_Raw_NonSpatial <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance&nbsp;&times;&nbsp;Tjur-R<sup>2</sup>)\\
    </span> --- <span style="font-size:10pt">Species are sorted by \\
    total explained variance, excluding spatial random effect</span>')

  Plot_Raw_NonSpatial <- DT_Raw_NonSpatial %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = VP_Value, y = Species_name, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = Title_Raw_NonSpatial) +
    Theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(Path_VarPar, "VarPar_Raw_ByTotalNonSpatial.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(Plot_Raw_NonSpatial)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Computing and plotting variance partitioning took ")

  return(invisible(NULL))
}
