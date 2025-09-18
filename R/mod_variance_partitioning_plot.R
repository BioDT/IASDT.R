## |------------------------------------------------------------------------| #
# variance_partitioning_plot ----
## |------------------------------------------------------------------------| #

#' @export
#' @name variance_partitioning
#' @rdname variance_partitioning
#' @order 2

variance_partitioning_plot <- function(
    path_model = NULL, env_file = ".env", vp_file = "varpar", use_tf = TRUE,
    tf_environ = NULL, n_cores = 1L, width = 30, height = 15, axis_text = 4,
    spatial_model = TRUE, is_cv_model = FALSE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::check_args(
    args_to_check = c("path_model", "vp_file"), args_type = "character")
  ecokit::check_args(
    args_to_check = c("use_tf", "spatial_model", "is_cv_model"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c("width", "height", "axis_text"), args_type = "numeric")
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ias_id <- species_name <- species <- variable <- vp_value <- species <-
    taxa_info_file <- sp <- TjurR2 <- plot_label <- vp_sum <-
    evaluation_type <- NULL

  # Set null device for `cairo`. This is to properly render the plots using
  # ggtext - https://github.com/wilkelab/cowplot/issues/73
  cowplot::set_null_device("cairo")

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_time("Loading data", cat_timestamp = FALSE)

  # Species info -----

  ecokit::cat_time("Loading species info", cat_timestamp = FALSE, level = 1L)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "taxa_info_file", "DP_R_taxa_info_rdata", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  sp_list <- ecokit::load_as(taxa_info_file) %>%
    dplyr::select(species = ias_id, species_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      species = stringr::str_pad(string = species, width = 4, pad = "0"),
      species = paste0("sp_", species))

  # # ..................................................................... ###

  path_root <- ecokit::parent_dir(
    path = path_model, levels = 1L, check_dir = TRUE, warning = FALSE)
  path_varpar <- fs::path(
    path_root, "model_postprocessing", "variance_partitioning")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Model evaluation ----

  eval_dir <- fs::path(path_root, "model_evaluation")
  if (!fs::dir_exists(eval_dir)) {
    ecokit::stop_ctx(
      "Model evaluation directory is not found: ", path_root = path_root,
      eval_dir = eval_dir, include_backtrace = TRUE)
  }

  if (is_cv_model) {

    ecokit::cat_time(
      "Loading training and testing evaluation data", level = 1L,
      cat_timestamp = FALSE)
    path_eval <- fs::path(eval_dir, "eval_cv_data.RData")
    if (!ecokit::check_data(path_eval, warning = FALSE)) {
      ecokit::stop_ctx(
        "Model evaluation file is not found: ", path_root = path_root,
        path_eval = path_eval, include_backtrace = TRUE)
    }
    mod_eval <- ecokit::load_as(path_eval)

    ecokit::cat_time(
      "Loading model training data", level = 1L, cat_timestamp = FALSE)
    mod_eval_train <- mod_eval %>%
      dplyr::filter(evaluation_type == "training") %>%
      dplyr::select(
        tidyselect::all_of(
          c("species", "RMSE", "AUC", "Boyce", "TjurR2"))) %>%
      dplyr::left_join(sp_list, by = "species")

    ecokit::cat_time(
      "Loading model testing data", level = 1L, cat_timestamp = FALSE)
    mod_eval_test <- mod_eval %>%
      dplyr::filter(evaluation_type == "testing") %>%
      dplyr::select(
        tidyselect::all_of(
          c("species", "RMSE", "AUC", "Boyce", "TjurR2"))) %>%
      dplyr::left_join(sp_list, by = "species")
    rm(mod_eval, envir = environment())

  } else {

    ecokit::cat_time(
      "Loading model training evaluation data", level = 1L,
      cat_timestamp = FALSE)

    path_eval_train <- fs::path(eval_dir, "eval_current_no_clamping.qs2")
    if (!ecokit::check_data(path_eval_train, warning = FALSE)) {
      ecokit::stop_ctx(
        "Model evaluation file is not found: ",
        path_root = path_root, path_eval_train = path_eval_train,
        include_backtrace = TRUE)
    }

    mod_eval_train <- ecokit::load_as(path_eval_train) %>%
      # filter out the species that are not in the model
      dplyr::filter(stringr::str_starts(ias_id, "sp_")) %>%
      dplyr::rename(species = ias_id) %>%
      dplyr::select(-sp) %>%
      dplyr::left_join(sp_list, by = "species")

    mod_eval_test <- NULL
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Compute or load variance partitioning ----

  ecokit::cat_time(
    "Compute or load variance partitioning", cat_timestamp = FALSE)

  if (is.null(vp_file)) {
    vp_file <- "varpar"
  }

  file_varpar <- fs::path(path_varpar, paste0(vp_file, ".RData"))

  if (ecokit::check_data(file_varpar, warning = FALSE)) {

    ecokit::cat_time(
      "Loading variance partitioning data", level = 1L, cat_timestamp = FALSE)
    varpar <- ecokit::load_as(file_varpar)

  } else {

    ecokit::cat_time(
      paste0(
        "Variance partitioning will be computed using ", n_cores, " cores ",
        dplyr::if_else(use_tf, "and", "without"), " `TensorFlow`."),
      level = 1L, cat_timestamp = FALSE)

    varpar <- IASDT.R::variance_partitioning_compute(
      path_model = path_model, n_cores = n_cores, use_tf = use_tf,
      tf_environ = tf_environ, verbose = TRUE, vp_file = vp_file)

  }

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plot theme ----

  plotting_theme <- ggplot2::theme(
    plot.title = ggtext::element_markdown(size = 14L, hjust = 0L, vjust = 0L),
    plot.subtitle = ggtext::element_markdown(size = 9L, hjust = 0L, vjust = 0L),
    plot.title.position = "plot",
    legend.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0.2, 0.15, 0.2, 0.05, "cm"),
    axis.text.y = ggplot2::element_text(
      size = 7L, margin = ggplot2::margin(r = 0L)),
    axis.text.x = ggplot2::element_text(
      face = "italic", size = axis_text, angle = 90L, hjust = 1L, vjust = 0.3,
      margin = ggplot2::margin(t = 0L)),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.text = ggtext::element_markdown(
      size = 8L, hjust = 0.5, lineheight = 1.15),
    legend.key.height = ggplot2::unit(0.4, "cm"),
    legend.key.width = ggplot2::unit(0.55, "cm"),
    legend.key.spacing.x = ggplot2::unit(0.2, "cm"),
    legend.position = "bottom",
    legend.margin = ggplot2::margin(0L, 2L, 0L, 1L, unit = "cm"),
    legend.box.spacing = ggplot2::unit(1L, "pt"))

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_time("Relative variance partitioning", cat_timestamp = FALSE)

  # Relative variance partitioning ----

  ## Plotting data ----
  ecokit::cat_time("Preparing data", level = 1L, cat_timestamp = FALSE)
  varpar_data <- tidyr::pivot_longer(
    data = varpar$vals, cols = -variable,
    names_to = "species", values_to = "vp_value") %>%
    dplyr::left_join(sp_list, by = "species")

  # Calculate mean Variance partitioning per variable and prepare labels for the
  # plot
  varpar_mean <- varpar_data %>%
    dplyr::summarise(
      vp_value = mean(vp_value, na.rm = TRUE), .by = "variable") %>%
    dplyr::arrange(dplyr::desc(vp_value)) %>%
    dplyr::mutate(
      plot_label = dplyr::case_when(
        variable == "habitat_log" ~ "Habitat coverage",
        variable == "road_rail_log" ~ "Road+Rail intensity",
        variable == "efforts_log" ~ "Sampling efforts",
        variable == "rivers_log" ~ "River length",
        variable == "npp" ~ "NPP",
        variable == "wetness" ~ "Wetness index",
        variable == "soil" ~ "Soil density",
        variable == "Random: sample" ~
          dplyr::if_else(spatial_model, "Spatial effect", "Random effect"),
        .default = variable),
      plot_label = paste0(
        "<b>", plot_label, "</b><br>(", round(vp_value * 100, 1), "%)"),
      plot_label = factor(plot_label, .data$plot_label),
      vp_value = NULL)

  # Order variables by the mean variance partitioning
  var_order <- dplyr::pull(varpar_mean, variable)

  # Order species by the mean variance partitioning per variable
  sp_order <- varpar_data %>%
    dplyr::mutate(variable = factor(variable, var_order)) %>%
    dplyr::summarise(
      vp_value = sum(vp_value), .by = c(species_name, variable)) %>%
    dplyr::arrange(variable, dplyr::desc(vp_value)) %>%
    dplyr::distinct(species_name) %>%
    dplyr::pull(species_name)

  # Order species by total variance, excluding the spatial random effect
  sp_order_nonspatial <- varpar_data %>%
    dplyr::filter(!startsWith(variable, "Random")) %>%
    dplyr::summarise(
      vp_value = sum(vp_value), .by = c("species", "species_name")) %>%
    dplyr::arrange(dplyr::desc(vp_value)) %>%
    dplyr::distinct(species_name) %>%
    dplyr::pull(species_name)


  # Plotting data for variance partitioning

  # 1. ordered by mean variance partitioning per variable
  data_relative <- dplyr::arrange(varpar_data, variable, vp_value) %>%
    dplyr::mutate(species_name = factor(species_name, sp_order)) %>%
    dplyr::left_join(varpar_mean, by = "variable")

  # 2. original species order
  sp_order_orig <- dplyr::distinct(varpar_data, species, species_name) %>%
    dplyr::arrange(species) %>%
    dplyr::pull(species_name)
  data_relative_orig <- data_relative %>%
    dplyr::mutate(species_name = factor(species_name, sp_order_orig))

  # 3a. ordered by Tjur-R2 - training evaluation
  sp_order_tjurr2_train <- mod_eval_train %>%
    dplyr::arrange(dplyr::desc(TjurR2)) %>%
    dplyr::pull(species_name)
  data_relative_tjurr2_train <- data_relative %>%
    dplyr::mutate(species_name = factor(species_name, sp_order_tjurr2_train))

  # 3b. ordered by Tjur-R2 - testing evaluation
  if (is_cv_model) {
    sp_order_tjurr2_test <- mod_eval_test %>%
      dplyr::arrange(dplyr::desc(TjurR2)) %>%
      dplyr::pull(species_name)
    data_relative_tjurr2_test <- data_relative %>%
      dplyr::mutate(species_name = factor(species_name, sp_order_tjurr2_test))
  }

  # 4. ordered by total variance partitioning, excluding spatial random effect
  data_relative_nonspatial <- data_relative %>%
    dplyr::mutate(species_name = factor(species_name, sp_order_nonspatial))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----

  ecokit::cat_time(
    "Plotting relative variance partitioning",
    level = 1L, cat_timestamp = FALSE)

  # 1. ordered by mean variance partitioning

  ecokit::cat_time(
    "1. ordered by mean variance partitioning",
    level = 2L, cat_timestamp = FALSE)

  title_relative <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:9.5pt">species are sorted by \\
    value of the most-important predictor</span>')

  plot_relative <- data_relative %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = vp_value, y = species_name, fill = plot_label)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_relative) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(path_varpar, "varpar_relative_by_mean.jpeg"),
    width = width, height = height, units = "cm", quality = 100, res = 600)
  print(plot_relative)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 2. ordered by original species order

  ecokit::cat_time(
    "2. ordered by original species order", level = 2L, cat_timestamp = FALSE)

  title_relative_orig <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:9.5pt">species are sorted by species \\
    taxonomy</span>')

  plot_relative_orig <- data_relative_orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_relative_orig) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  ragg::agg_jpeg(
    filename = fs::path(path_varpar, "varpar_relative_by_taxonomy.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(plot_relative_orig)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 3a. ordered by Tjur-R2 - training

  ecokit::cat_time(
    "3a. ordered by Tjur-R2 - training", level = 2L, cat_timestamp = FALSE)

  title_relative_tjurr2_train <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b> \\
      --- </span><span style="font-size:9.5pt">species are sorted by \\
      explanatory power (training Tjur-R<sup>2</sup>)</span>')

  plot_relative_tjurr2_train <- data_relative_tjurr2_train %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_relative_tjurr2_train) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(path_varpar, "varpar_relative_by_tjurr2_training.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(plot_relative_tjurr2_train)
  grDevices::dev.off()
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 3b. ordered by Tjur-R2 - testing

  if (is_cv_model) {

    ecokit::cat_time(
      "3b. ordered by Tjur-R2 - testing", level = 2L, cat_timestamp = FALSE)

    title_relative_tjurr2_test <- stringr::str_glue(
      '<span style="font-size:12pt"><b>Proportion of explained variance</b> \\
      --- </span><span style="font-size:9.5pt">species are sorted by \\
      predictive power (testing Tjur-R<sup>2</sup>)</span>')

    plot_relative_tjurr2_test <- data_relative_tjurr2_test %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = vp_value, y = species_name, fill = plot_label),
        environment = emptyenv()) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::theme_classic() +
      ggplot2::ylab("species") +
      ggplot2::xlab("Proportion of variance explained") +
      ggplot2::scale_fill_brewer(palette = "Paired") +
      ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
      ggplot2::labs(title = title_relative_tjurr2_test) +
      plotting_theme +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = fs::path(
        path_varpar, "varpar_relative_by_tjurr2_testing.jpeg"),
      width = width, height = height, res = 600, quality = 100, units = "cm")
    plot(plot_relative_tjurr2_test)
    grDevices::dev.off()
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # 4. ordered by total variance partitioning, excluding spatial random effect

  ecokit::cat_time(
    "4. ordered by total variance partitioning, excluding random effect",
    level = 2L, cat_timestamp = FALSE)
  label0 <- dplyr::if_else(spatial_model, "spatial random", "random")

  title_nonspatial <- stringr::str_glue(
    '<span style="font-size:12pt"><b>Proportion of explained variance</b>\\
    </span> --- <span style="font-size:9.5pt">species are sorted by total \\
    explained variance, excluding {label0} effect</span>')

  plot_relative_nonspatial <- data_relative_nonspatial %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Proportion of variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_nonspatial) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  ragg::agg_jpeg(
    filename = fs::path(
      path_varpar, "varpar_relative_by_total_nonspatial.jpeg"),
    width = width, height = height, units = "cm", quality = 100, res = 600)
  print(plot_relative_nonspatial)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Raw variance partitioning - training ----

  ecokit::cat_time(
    "Raw variance partitioning - training", cat_timestamp = FALSE)

  ## Plotting data ----
  ecokit::cat_time("Preparing data", level = 1L, cat_timestamp = FALSE)

  if (is.null(mod_eval_train$TjurR2) ||
      length(mod_eval_train$TjurR2) != ncol(varpar$vals) - 1) {
    ecokit::stop_ctx(
      paste0(
        "Mismatch between the length of mod_eval_train$TjurR2 and the ",
        " number of columns in varpar$vals"),
      mod_eval_train_tjurr2 = mod_eval_train$TjurR2,
      length_mod_eval_train_tjurr2 = length(mod_eval_train$TjurR2),
      ncol_varpar_vals = ncol(varpar$vals) - 1L, include_backtrace = TRUE)
  }

  varpar_data_raw_train <- tidyr::pivot_longer(
    data = varpar$vals, cols = -variable,
    names_to = "species", values_to = "vp_value") %>%
    dplyr::left_join(mod_eval_train, by = "species") %>%
    dplyr::mutate(vp_value = vp_value * TjurR2)

  varpar_raw_mean_train <- varpar_data_raw_train %>%
    dplyr::summarise(
      vp_value = mean(vp_value, na.rm = TRUE), .by = "variable") %>%
    dplyr::arrange(dplyr::desc(vp_value)) %>%
    dplyr::mutate(
      plot_label = dplyr::case_when(
        variable == "habitat_log" ~ "Habitat coverage",
        variable == "road_rail_log" ~ "Road+Rail intensity",
        variable == "efforts_log" ~ "Sampling efforts",
        variable == "rivers_log" ~ "River length",
        variable == "npp" ~ "NPP",
        variable == "wetness" ~ "Wetness index",
        variable == "soil" ~ "Soil density",
        variable == "Random: sample" ~
          dplyr::if_else(spatial_model, "Spatial effect", "Random effect"),
        .default = variable),
      plot_label = paste0(
        "<b>", plot_label, "</b><br>(", round(vp_value * 100, 1), "%)"),
      plot_label = factor(plot_label, .data$plot_label),
      vp_value = NULL)

  # Plotting data for relative variance partitioning - ordered by mean
  # variance partitioning per variable
  data_raw_train <- varpar_data_raw_train %>%
    dplyr::arrange(variable, vp_value) %>%
    dplyr::mutate(species_name = factor(species_name, sp_order)) %>%
    dplyr::left_join(varpar_raw_mean_train, by = "variable") %>%
    dplyr::mutate(species_name = as.character(species_name))

  # Plotting data for relative variance partitioning - original species order
  data_raw_train_orig <- data_raw_train %>%
    dplyr::mutate(species_name = factor(species_name, sp_order_orig))

  # Plotting data for relative variance partitioning - original species order
  sp_order_total_raw_train <- data_raw_train %>%
    dplyr::summarize(
      vp_sum = sum(vp_value), .by = c(species, species_name)) %>%
    dplyr::arrange(dplyr::desc(vp_sum)) %>%
    dplyr::pull(species_name)
  data_raw_train_total_raw_train <- data_raw_train %>%
    dplyr::mutate(species_name = factor(species_name, sp_order_total_raw_train))

  # Order species by total variance, excluding the spatial random effect
  sp_order_raw_train_nonspatial <- data_raw_train %>%
    dplyr::filter(!startsWith(variable, "Random")) %>%
    dplyr::summarise(
      vp_value = sum(vp_value), .by = c("species", "species_name")) %>%
    dplyr::arrange(dplyr::desc(vp_value)) %>%
    dplyr::distinct(species_name) %>%
    dplyr::pull(species_name)
  data_raw_train_nonspatial <- data_raw_train %>%
    dplyr::mutate(
      species_name = factor(species_name, sp_order_raw_train_nonspatial))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## Plotting ----
  ecokit::cat_time("Plotting", level = 1L, cat_timestamp = FALSE)

  ## 1. ordered by mean variance partitioning
  ecokit::cat_time(
    "1. ordered by mean variance partitioning", level = 2L,
    cat_timestamp = FALSE)

  title_raw_train <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; explanatory Tjur-R<sup>2\\
    </sup>)</span> --- <span style="font-size:9.5pt">species are sorted by \\
    taxonomy</span>')

  plot_raw_train <- data_raw_train_orig %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_raw_train) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(path_varpar, "varpar_raw_train_by_taxonomy.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(plot_raw_train)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## 2. ordered by Total Variance Explained
  ecokit::cat_time(
    "2. ordered by total variance explained", level = 2L, cat_timestamp = FALSE)

  title_raw_train_total_raw <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; explanatory Tjur-R<sup>2\\
    </sup>)</span> --- <span style="font-size:9.5pt">species are sorted \\
    by total explained variance</span>')

  plot_raw_train_total_raw <- data_raw_train_total_raw_train %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_raw_train_total_raw) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(path_varpar, "varpar_raw_train_by_mean.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(plot_raw_train_total_raw)
  grDevices::dev.off()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## 3. by total variance partitioning, excluding spatial random effect
  ecokit::cat_time(
    "3. by total variance partitioning, excluding spatial random effect",
    level = 2L, cat_timestamp = FALSE)

  title_raw_train_nonspatial <- stringr::str_glue(
    '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance&nbsp;&times;&nbsp; explanatory \\
    Tjur-R<sup>2</sup>)</span> --- <span style="font-size:9.5pt">species \\
    are sorted by total explained variance, excluding {label0} \\
    effect</span>')

  plot_raw_train_nonspatial <- data_raw_train_nonspatial %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = vp_value, y = species_name, fill = plot_label),
      environment = emptyenv()) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("species") +
    ggplot2::xlab("Raw variance explained") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::labs(title = title_raw_train_nonspatial) +
    plotting_theme +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(
      path_varpar, "varpar_raw_train_by_total_nonspatial.jpeg"),
    width = width, height = height, res = 600, quality = 100, units = "cm")
  plot(plot_raw_train_nonspatial)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  if (is_cv_model) {

    # Raw variance partitioning - testing ----

    ## Plotting data ----
    ecokit::cat_time(
      "Raw variance partitioning - testing", cat_timestamp = FALSE)
    ecokit::cat_time("Preparing data", level = 1L, cat_timestamp = FALSE)

    if (is.null(mod_eval_test$TjurR2) ||
        length(mod_eval_test$TjurR2) != ncol(varpar$vals) - 1) {
      ecokit::stop_ctx(
        paste0(
          "Mismatch between the length of mod_eval_test$TjurR2 and the ",
          " number of columns in varpar$vals"),
        mod_eval_test_tjurr2 = mod_eval_test$TjurR2,
        length_mod_eval_test_tjurr2 = length(mod_eval_test$TjurR2),
        ncol_varpar_vals = ncol(varpar$vals) - 1L, include_backtrace = TRUE)
    }

    varpar_data_raw_test <- tidyr::pivot_longer(
      data = varpar$vals, cols = -variable,
      names_to = "species", values_to = "vp_value") %>%
      dplyr::left_join(mod_eval_test, by = "species") %>%
      dplyr::mutate(vp_value = vp_value * TjurR2)

    varpar_raw_mean_test <- varpar_data_raw_test %>%
      dplyr::summarise(
        vp_value = mean(vp_value, na.rm = TRUE), .by = "variable") %>%
      dplyr::arrange(dplyr::desc(vp_value)) %>%
      dplyr::mutate(
        plot_label = dplyr::case_when(
          variable == "habitat_log" ~ "Habitat coverage",
          variable == "road_rail_log" ~ "Road+Rail intensity",
          variable == "efforts_log" ~ "Sampling efforts",
          variable == "rivers_log" ~ "River length",
          variable == "npp" ~ "NPP",
          variable == "wetness" ~ "Wetness index",
          variable == "soil" ~ "Soil density",
          variable == "Random: sample" ~
            dplyr::if_else(spatial_model, "Spatial effect", "Random effect"),
          .default = variable),
        plot_label = paste0(
          "<b>", plot_label, "</b><br>(", round(vp_value * 100, 1), "%)"),
        plot_label = factor(plot_label, .data$plot_label),
        vp_value = NULL)

    # Plotting data for relative variance partitioning - ordered by mean
    # variance partitioning per variable
    data_raw_test <- varpar_data_raw_test %>%
      dplyr::arrange(variable, vp_value) %>%
      dplyr::mutate(species_name = factor(species_name, sp_order)) %>%
      dplyr::left_join(varpar_raw_mean_test, by = "variable") %>%
      dplyr::mutate(
        species_name = as.character(species_name),
        # replace NAs with zero for cases for which there is no testing
        # evaluation data
        vp_value = ifelse(is.na(vp_value), 0, vp_value))

    # Plotting data for relative variance partitioning - original species order
    data_raw_test_orig <- data_raw_test %>%
      dplyr::mutate(species_name = factor(species_name, sp_order_orig))

    # Plotting data for relative variance partitioning - original species order
    sp_order_total_raw_test <- data_raw_test %>%
      dplyr::summarize(
        vp_sum = sum(vp_value), .by = c(species, species_name)) %>%
      dplyr::arrange(dplyr::desc(vp_sum)) %>%
      dplyr::pull(species_name)
    data_raw_test_total_raw_test <- data_raw_test %>%
      dplyr::mutate(
        species_name = factor(species_name, sp_order_total_raw_test))

    # Order species by total variance, excluding the spatial random effect
    sp_order_raw_test_nonspatial <- data_raw_test %>%
      dplyr::filter(!startsWith(variable, "Random")) %>%
      dplyr::summarise(
        vp_value = sum(vp_value), .by = c("species", "species_name")) %>%
      dplyr::arrange(dplyr::desc(vp_value)) %>%
      dplyr::distinct(species_name) %>%
      dplyr::pull(species_name)
    data_raw_test_nonspatial <- data_raw_test %>%
      dplyr::mutate(
        species_name = factor(species_name, sp_order_raw_test_nonspatial))

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

    ## Plotting ----
    ecokit::cat_time("Plotting", level = 1L, cat_timestamp = FALSE)

    ## 1. ordered by mean variance partitioning
    ecokit::cat_time(
      "1. ordered by mean variance partitioning", level = 2L,
      cat_timestamp = FALSE)

    title_raw_test <- stringr::str_glue(
      '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; predictive Tjur-R<sup>2\\
    </sup>)</span> --- <span style="font-size:9.5pt">species are sorted by \\
    taxonomy</span>')

    plot_raw_test <- data_raw_test_orig %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = vp_value, y = species_name, fill = plot_label),
        environment = emptyenv()) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::theme_classic() +
      ggplot2::ylab("species") +
      ggplot2::xlab("Raw variance explained") +
      ggplot2::scale_fill_brewer(palette = "Paired") +
      ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
      ggplot2::labs(title = title_raw_test) +
      plotting_theme +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = fs::path(path_varpar, "varpar_raw_test_by_taxonomy.jpeg"),
      width = width, height = height, res = 600, quality = 100, units = "cm")
    plot(plot_raw_test)
    grDevices::dev.off()

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

    ## 2. ordered by Total Variance Explained
    ecokit::cat_time(
      "2. ordered by total variance explained", level = 2L,
      cat_timestamp = FALSE)

    title_raw_test_total_raw <- stringr::str_glue(
      '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance &times; predictive Tjur-R<sup>2\\
    </sup>)</span> --- <span style="font-size:9.5pt">species are sorted \\
    by total explained variance</span>')

    plot_raw_test_total_raw <- data_raw_test_total_raw_test %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = vp_value, y = species_name, fill = plot_label),
        environment = emptyenv()) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::theme_classic() +
      ggplot2::ylab("species") +
      ggplot2::xlab("Raw variance explained") +
      ggplot2::scale_fill_brewer(palette = "Paired") +
      ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
      ggplot2::labs(title = title_raw_test_total_raw) +
      plotting_theme +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = fs::path(path_varpar, "varpar_raw_test_by_mean.jpeg"),
      width = width, height = height, res = 600, quality = 100, units = "cm")
    plot(plot_raw_test_total_raw)
    grDevices::dev.off()

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

    ## 3. by total variance partitioning, excluding spatial random effect
    ecokit::cat_time(
      "3. by total variance partitioning, excluding spatial random effect",
      level = 2L, cat_timestamp = FALSE)

    title_raw_test_nonspatial <- stringr::str_glue(
      '<b>Proportion of raw variance</b> <span style="font-size:12pt">\\
    (proportion of explained variance&nbsp;&times;&nbsp; predictive \\
    Tjur-R<sup>2</sup>)</span> --- <span style="font-size:9.5pt">species \\
    are sorted by total explained variance, excluding {label0} \\
    effect</span>')

    plot_raw_test_nonspatial <- data_raw_test_nonspatial %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = vp_value, y = species_name, fill = plot_label),
        environment = emptyenv()) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::theme_classic() +
      ggplot2::ylab("species") +
      ggplot2::xlab("Raw variance explained") +
      ggplot2::scale_fill_brewer(palette = "Paired") +
      ggplot2::coord_flip(xlim = c(0, 1), expand = FALSE) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
      ggplot2::labs(title = title_raw_test_nonspatial) +
      plotting_theme +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = fs::path(
        path_varpar, "varpar_raw_test_by_total_nonspatial.jpeg"),
      width = width, height = height, res = 600, quality = 100, units = "cm")
    plot(plot_raw_test_nonspatial)
    grDevices::dev.off()

  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Computing and plotting variance partitioning took ")

  rm(label0, envir = environment())
  return(invisible(NULL))
}
