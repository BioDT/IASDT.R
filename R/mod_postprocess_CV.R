## |------------------------------------------------------------------------| #
# Model pipeline for post-processing cross-validated Hmsc models
## |------------------------------------------------------------------------| #

#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 5
#' @author Ahmed El-Gabbas
#' @export

mod_postprocess_CV_1_CPU <- function(
    model_dir = NULL, CV_names = NULL, n_cores = 8L, env_file = ".env",
    from_JSON = FALSE, use_TF = TRUE, TF_use_single = FALSE,
    TF_environ = NULL, LF_n_cores = n_cores, LF_only = TRUE,
    LF_temp_cleanup = TRUE, LF_check = FALSE, LF_runtime = "01:00:00",
    temp_cleanup = TRUE,  n_batch_files = 210L, working_directory = NULL,
    partition_name = "small-g") {

  # ****************************************************************

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CV <- CV_name <- NULL

  # ****************************************************************

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "from_JSON", "use_TF", "TF_use_single", "LF_only", "LF_temp_cleanup",
      "LF_check", "temp_cleanup"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "env_file", "partition_name", "LF_runtime"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores", "n_batch_files"))
  rm(AllArgs, envir = environment())

  if (n_batch_files <= 0) {
    ecokit::stop_ctx(
      "`n_batch_files` must be a positive integer.",
      n_batch_files = n_batch_files, include_backtrace = TRUE)
  }
  if (n_cores <= 0) {
    ecokit::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }
  if (LF_n_cores <= 0) {
    ecokit::stop_ctx(
      "`LF_n_cores` must be a positive integer.", LF_n_cores = LF_n_cores,
      include_backtrace = TRUE)
  }

  if (!file.exists(env_file)) {
    ecokit::stop_ctx(
      "Error: Environment file is invalid or does not exist.",
      env_file = env_file, include_backtrace = TRUE)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory is invalid or does not exist.", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  valid_CVs <- c("CV_Dist", "CV_Large", "CV_SAC")
  if (!all(CV_names %in% valid_CVs)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for CV_names argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_names = CV_names, include_backtrace = TRUE)
  }

  # ****************************************************************

  # # Load environment variables, for project ID
  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "ProjectID", "DP_R_LUMI_gpu", FALSE, FALSE)

  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # ****************************************************************

  # Paths to files and directories

  CV_dir <- fs::path(model_dir, "Model_Fitting_CV")

  Temp_dir <- fs::path(CV_dir, "Temp")
  # Path to store TF commands
  Path_TF <- fs::path(CV_dir, "LF_TF_commands")
  # Path to store log files
  Path_Log <- fs::path(Path_TF, "log")

  fs::dir_create(c(Path_TF, Path_Log))

  CV_DT_fitted <- fs::path(CV_dir, "CV_DT_fitted.RData")
  if (!file.exists(CV_DT_fitted)) {
    ecokit::stop_ctx(
      "CV_DT_fitted file not found.", CV_DT_fitted = CV_DT_fitted,
      include_backtrace = TRUE)
  }

  Path_LF_SLURM <- fs::path(Path_TF, "LF_SLURM.slurm")
  path_out <- fs::path(Path_Log, "%x-%A-%a.out")      # nolint: object_name_linter

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Merge chains -----

  ecokit::info_chunk(
    "Merge chains", line_char = "|", line_char_rep = 60L, cat_red = TRUE,
    cat_bold = TRUE, cat_timestamp = FALSE, info_lines_before = 2L)

  IASDT.R::mod_merge_chains_CV(
    model_dir = model_dir, n_cores = n_cores, CV_names = CV_names,
    from_JSON = FALSE, out_extension = "qs2")

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Prepare scripts for latent factor processing -----

  ecokit::info_chunk(
    "Prepare scripts for latent factor processing",
    line_char = "|", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
    cat_timestamp = FALSE, info_lines_before = 2L)

  CV_DT_fitted <- ecokit::load_as(CV_DT_fitted) %>%
    dplyr::mutate(
      LF = purrr::map2(
        .x = CV_name, .y = CV,
        .f = ~ {

          ecokit::info_chunk(
            message = paste0(.x, "_", .y), level = 1L, line_char = "+",
            line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE)

          IASDT.R::predict_maps_CV(
            model_dir = model_dir, CV_name = paste0("CV_", .x),
            CV_fold = .y, n_cores = n_cores, use_TF = use_TF,
            TF_environ = TF_environ, TF_use_single = TF_use_single,
            LF_n_cores = LF_n_cores, LF_check = LF_check,
            LF_temp_cleanup = LF_temp_cleanup, LF_only = LF_only,
            LF_commands_only = TRUE, temp_cleanup = temp_cleanup)
        }
      ))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Merge TensorFlow commands into batch files -----

  ecokit::info_chunk(
    "Merge TensorFlow commands into batch files", line_char = "|",
    line_char_rep = 60L, cat_red = TRUE,
    cat_bold = TRUE, info_lines_before = 2L)

  ecokit::cat_time(
    paste0(
      "Merge and organise TensorFlow commands for LF predictions ",
      "into a maximum of ", n_batch_files, " files"),
    level = 1L, cat_timestamp = FALSE)

  # Basic commands for TensorFlow setup
  BasicCommands <- c(
    "#!/bin/bash\n",
    "# Load TensorFlow module and configure environment",
    "ml use /appl/local/csc/modulefiles",
    "ml tensorflow\n",
    "export TF_CPP_MIN_LOG_LEVEL=3",
    "export TF_ENABLE_ONEDNN_OPTS=0\n",
    "# Verify GPU availability",
    paste0(
      'python3 -c "import tensorflow as tf; ',
      'print(\\\"Num GPUs Available:\\\", ',
      'len(tf.config.list_physical_devices(\\\"GPU\\\")))"'),
    "")

  # Change working directory if specified
  if (!is.null(working_directory)) {
    working_directory <- ecokit::normalize_path(
      working_directory, must_work = TRUE)
    BasicCommands <- c(
      BasicCommands, "# Change to working directory",
      paste0("cd ", working_directory), "")
  }

  # Find list of files matching the pattern
  LF_InFiles <- list.files(
    path = Temp_dir, pattern = "^LF_NewSites_Commands_.+.txt",
    recursive = TRUE, full.names = TRUE) %>%
    gtools::mixedsort()

  if (length(LF_InFiles) == 0) {
    ecokit::stop_ctx(
      "No command files were found in the temp directory",
      Temp_dir = Temp_dir, LF_InFiles = LF_InFiles,
      length_LF_InFiles = length(LF_InFiles), include_backtrace = TRUE)
  }

  ecokit::cat_time(
    paste0("Found ", length(LF_InFiles), " files"),
    level = 2L, cat_timestamp = FALSE)
  stringr::str_remove_all(LF_InFiles, paste0(Temp_dir, "/")) %>%
    purrr::walk(ecokit::cat_time, level = 3L, cat_timestamp = FALSE)

  # Read and merge commands from input files
  LF_commands <- purrr::map(LF_InFiles, readr::read_lines, progress = FALSE) %>%
    unlist() %>%
    gtools::mixedsort()

  ecokit::cat_time(
    paste0("Total number of commands to be executed: ", length(LF_commands)),
    level = 2L, cat_timestamp = FALSE)

  if (length(LF_commands) < n_batch_files) {
    ecokit::cat_time(
      paste0(
        "Fewer commands than the requested number of files. ",
        "Setting `n_batch_files=", n_batch_files, "`."),
      level = 2L, cat_timestamp = FALSE)
    n_batch_files <- length(LF_commands)
  }

  ecokit::cat_time(
    paste0("Splitting commands into ", n_batch_files, " files"),
    cat_timestamp = FALSE, level = 2L)
  LF_commands <- ecokit::split_vector(LF_commands, n_splits = n_batch_files)
  line_sep <- strrep("-", 60)      # nolint: object_name_linter

  purrr::walk(
    .x = seq_len(length(LF_commands)),
    .f = ~ {

      Chunk_number <- stringr::str_pad(
        .x, pad = "0", width = nchar(n_batch_files))
      File <- fs::path(Path_TF, paste0("TF_Chunk_", Chunk_number, ".txt"))
      time_now <- lubridate::now(tzone = "CET") %>%
        format(format = "%Y-%m-%d %H:%M:%S")

      readr::write_lines(x = BasicCommands, file = File, append = FALSE)
      readr::write_lines(
        x = paste0(
          "# ", length(LF_commands[[.x]]), " commands to be executed:"),
        file = File, append = TRUE)
      readr::write_lines(x = LF_commands[[.x]], file = File, append = TRUE)
      readr::write_lines(
        x = c(
          paste0("\n#", line_sep),
          paste0("# This script was created on: ", time_now),
          paste0("#", line_sep)),
        file = File, append = TRUE)

      return(invisible(NULL))
    })

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Prepare LF batch file -----
  ecokit::cat_time("Prepare LF batch file")

  time_now <- lubridate::now(tzone = "CET") %>%      # nolint: object_name_linter
    format(format = "%Y-%m-%d %H:%M:%S")

  LF_slurm_script <- stringr::str_glue(
    "#!/bin/bash\n",
    "#SBATCH --job-name=PP_LF\n",
    "#SBATCH --ntasks=1\n",
    "#SBATCH --ntasks-per-node=1\n",
    "#SBATCH --account={ProjectID}\n",
    "#SBATCH --cpus-per-task=1\n",
    "#SBATCH --gpus-per-node=1\n",
    "#SBATCH --time={LF_runtime}\n",
    "#SBATCH --partition={partition_name}\n",
    "#SBATCH --output={path_out}\n",
    "#SBATCH --error={path_out}\n",
    "#SBATCH --array=1-{n_batch_files}\n\n",

    "# Define the directory where TF_Chunk_*.txt scripts exist.",
    "Each array task will pick one of these files to execute.\n",
    'OutputDir="{Path_TF}"\n\n',

    '# Generate filename by zero-padding the "SLURM_ARRAY_TASK_ID" to ',
    'three digits (001, 002, ..., 210) and prefixing with "TF_Chunk_"\n',
    'SplitFile="${{OutputDir}}/TF_Chunk_$(printf "%03d" ',
    '"$SLURM_ARRAY_TASK_ID").txt"\n\n',

    "# Check if the target file exists\n",
    'if [[ ! -f "$SplitFile" ]]; then\n',
    '    echo "Error: File $SplitFile not found." >&2\n',
    "    exit 1\nfi\n\n",

    "# Report which file is being processed\n",
    'echo "Processing file: $SplitFile"\n\n',
    '# Execute the selected chunk script\nbash "$SplitFile"\n\n',

    "# Print a timestamped end-of-job message\n",
    'echo "End of program at $(date)"\n\n# {line_sep}\n',

    "# This script was created on: {time_now} CET\n",
    "# {line_sep}")

  ecokit::cat_time(
    paste0("Writing SLURM script to: `", Path_LF_SLURM, "`"),
    level = 2L, cat_timestamp = FALSE)
  # Write the content to a file
  readr::write_lines(LF_slurm_script, Path_LF_SLURM, append = FALSE)
  # Make the file executable
  Sys.chmod(Path_LF_SLURM, mode = "755")

  # ****************************************************************
  # ****************************************************************

  return(invisible(NULL))

}


# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_postprocess_CV_2_CPU ----
## |------------------------------------------------------------------------| #

#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 6
#' @author Ahmed El-Gabbas
#' @export

mod_postprocess_CV_2_CPU <- function(
    model_dir = NULL, CV_names = NULL, n_cores = 8L, env_file = ".env",
    use_TF = TRUE, TF_use_single = FALSE, temp_cleanup = TRUE,
    LF_temp_cleanup = TRUE, TF_environ = NULL, LF_n_cores = n_cores,
    LF_check = FALSE) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CV <- CV_name <- n_grids <- Sp <- IAS_ID <- pred_mean <- pred <- pred_sd <-
    mean_minus <- mean_plus <- summary_vals <- ias_id <- n_grids_pres_mean <-
    data <- hab_abb <- hab_name <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("use_TF", "TF_use_single", "LF_temp_cleanup", "LF_check"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "env_file"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores"))
  rm(AllArgs, envir = environment())

  if (n_cores <= 0) {
    ecokit::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }
  if (LF_n_cores <= 0) {
    ecokit::stop_ctx(
      "`LF_n_cores` must be a positive integer.", LF_n_cores = LF_n_cores,
      include_backtrace = TRUE)
  }

  valid_CVs <- c("CV_Dist", "CV_Large", "CV_SAC")
  if (!all(CV_names %in% valid_CVs)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for CV_names argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_names = CV_names, include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Load model's explanatory power ------

  Eval_explain <- list.files(
    path = fs::path(model_dir, "Model_Evaluation"),
    pattern = "Eval_.+.qs2", full.names = TRUE)

  if (length(Eval_explain) != 1L) {
    ecokit::stop_ctx(
      "There should be only one evaluation file in the Model_Evaluation folder",
      Eval_explain = Eval_explain, include_backtrace = TRUE)
  }

  Eval_explain <- ecokit::load_as(Eval_explain) %>%
    dplyr::select(-Sp) %>%
    dplyr::rename(ias_id = IAS_ID) %>%
    dplyr::filter(ias_id != "SR") %>%
    dplyr::rename_with(.cols = -ias_id, .fn = ~ paste0(.x, "_exp"))

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Predicting habitat suitability at testing cross-validation folds -------

  ecokit::info_chunk(
    message = "Predicting habitat suitability at testing sites",
    line_char = "*", line_char_rep = 70L, cat_red = TRUE, cat_bold = TRUE)

  path_CV <- fs::path(model_dir, "Model_Fitting_CV")
  path_CV_DT_fitted <- fs::path(path_CV, "CV_DT_fitted.RData")
  if (!file.exists(path_CV_DT_fitted)) {
    ecokit::stop_ctx(
      "Required file 'CV_DT_fitted.RData' not found",
      path_CV_DT_fitted = path_CV_DT_fitted, include_backtrace = TRUE)
  }

  CV_DT_fitted <- ecokit::load_as(path_CV_DT_fitted)
  if (nrow(CV_DT_fitted) < 1) {
    ecokit::stop_ctx(
      "CV_DT_fitted data should contain at least one row",
      CV_DT_fitted = CV_DT_fitted, include_backtrace = TRUE)
  }

  CV_DT_fitted <- CV_DT_fitted %>%
    dplyr::mutate(
      Path_Preds_summary = purrr::map2_chr(
        .x = CV_name, .y = CV,
        .f = ~ {

          ecokit::info_chunk(
            message = paste0(.x, "_", .y), line_char = "~", level = 1L,
            line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE)

          IASDT.R::predict_maps_CV(
            model_dir = model_dir, CV_name = paste0("CV_", .x),
            CV_fold = .y, n_cores = n_cores, env_file = env_file,
            use_TF = use_TF, TF_environ = TF_environ,
            TF_use_single = TF_use_single, LF_n_cores = LF_n_cores,
            LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup,
            LF_only = FALSE, LF_commands_only = FALSE,
            temp_cleanup = temp_cleanup)
        }
      ))

  ecokit::save_as(
    object = CV_DT_fitted, object_name = "CV_DT_fitted",
    out_path = path_CV_DT_fitted)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Merge prediction summary data ------

  ecokit::cat_time("Merge prediction summary data")

  # ID columns to be used for nesting
  cols_id <- c(
    "hab_abb", "hab_name", "CV_name", "ias_id", "taxon_name", "species_name",
    "class", "order", "family")

  # Evaluation metrics to be used for summary
  eval_metrics <- c("RMSE", "AUC", "Boyce", "TjurR2")

  # Summary columns: define columns for summary statistics (mean and sd) of
  # evaluation metrics plus presence count
  cols_summary <- outer(
    X = eval_metrics, Y = c("mean", "sd"),
    FUN = function(x, y) paste0(x, "_", y)) %>%
    as.vector() %>%
    sort() %>%
    c("n_grids_pres_mean")

  summary_all <- CV_DT_fitted$Path_Preds_summary %>%
    # load summary data for predictions
    purrr::map(ecokit::load_as) %>%
    dplyr::bind_rows() %>%
    # nest selected columns into `data` column
    tidyr::nest(.by = tidyselect::all_of(cols_id)) %>%
    # exclude species richness results
    dplyr::filter(ias_id != "SR") %>%
    dplyr::mutate(
      data = purrr::map(
        .x = data,
        .f = ~{
          .x %>%
            dplyr::arrange(CV_fold) %>%
            # convert each column to a names list
            dplyr::summarise(
              dplyr::across(
                .cols = tidyselect::everything(),
                .fns = function(y) {
                  as.vector(y) %>%
                    stats::setNames(paste0("CV_", seq_len(length(.)))) %>%
                    list()
                })) %>%
            dplyr::mutate(
              # calculate the mean and sd for each evaluation metric
              dplyr::across(
                .cols = tidyselect::all_of(eval_metrics),
                list(mean = ~mean(unlist(.)), sd = ~sd(unlist(.))),
                .names = "{.col}_{.fn}"),
              # mean number of testing presences
              n_grids_pres_mean = mean(unlist(n_grids_pres)))
        }),
      # a list column containing the summary values
      summary_vals = purrr::map(
        .x = data, .f = dplyr::select, tidyselect::all_of(cols_summary)),
      # update `data` list column to exclude summary values
      data = purrr::map(
        .x = data, .f = dplyr::select, -tidyselect::all_of(cols_summary))) %>%
    # unnest data for summary columns
    tidyr::unnest_wider(summary_vals) %>%
    # add explained power columns
    dplyr::left_join(Eval_explain, by = "ias_id")

  hab_type <- dplyr::distinct(summary_all, hab_abb, hab_name)
  if (nrow(hab_type) != 1) {
    ecokit::stop_ctx(
      "There should be only one habitat type in the summary data",
      hab_type = hab_type, include_backtrace = TRUE)
  }
  hab_type <- paste(unlist(hab_type), collapse = ": ") %>%
    stringr::str_to_lower()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting predictive power -----

  # plot different evaluation metrics vs mean number of testing presences (log10
  # scale)

  plot_cv_metric <- function(metric = NULL, data = summary_all) {

    # Ensure metric is non-null, character strings with content
    if (is.null(metric) || !is.character(metric) || !nzchar(metric))  {
      ecokit::stop_ctx(
        "`metric` have to be character of length > 1", metric = metric,
        include_backtrace = TRUE)
    }
    # Ensure data is a data frame tibble with more than 1 row
    if (!is.data.frame(data) || nrow(data) < 1) {
      ecokit::stop_ctx(
        "`data` has to be a data frame tibble with more than 1 row",
        data = data, include_backtrace = TRUE)
    }

    # metric value has to be exist in the eval_metrics vector
    if (!metric %in% eval_metrics) {
      ecokit::stop_ctx(
        paste0(
          "`metric` should be one of the following: ", toString(eval_metrics)),
        metric = metric, include_backtrace = TRUE)
    }

    # column names for mean and sd of evaluation metrics
    mean_col <- paste0(metric, "_mean")
    sd_col <- paste0(metric, "_sd")

    data2 <- data %>%
      dplyr::select(n_grids_pres_mean, dplyr::all_of(c(mean_col, sd_col))) %>%
      stats::setNames(c("n_grids", "pred_mean", "pred_sd")) %>%
      # calculate values for error bars
      dplyr::mutate(
        mean_plus = pred_mean + pred_sd,
        mean_minus = pred_mean - pred_sd)

    # if the metric is AUC/TjurR2, ensure that mean_plus <=1 and mean_minus >=0
    if (metric %in% c("AUC", "TjurR2")) {
      data2 <- dplyr::mutate(
        data2,
        mean_plus = pmin(mean_plus, 1), mean_minus = pmax(mean_minus, 0))
    }

    # if the metric is RMSE, ensure that the mean_minus >= 0
    if (metric == "RMSE") {
      data2 <- dplyr::mutate(data2, mean_minus = pmax(mean_minus, 0))
    }

    # if the metric is Boyce, ensure that mean_minus >=-1 and mean_plus <=1
    if (metric == "Boyce") {
      data2 <- dplyr::mutate(
        data2,
        mean_plus = pmin(mean_plus, 1), mean_minus = pmax(mean_minus, -1))
    }

    # expand value for x and y axes
    expand_vals <- c(0.005, 0.005)
    # Plot title
    plot_title <- dplyr::case_when(
      metric == "RMSE" ~ "Root mean square error (RMSE)",
      metric == "AUC" ~ "Area under the ROC curve (AUC)",
      metric == "Boyce" ~ "Continuous Boyce index",
      metric == "TjurR2" ~ "Tjur R-squared",
      .default = NULL)

    data2 %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = n_grids, y = pred_mean, ymin = mean_minus, ymax = mean_plus)) +
      # plot error bar then the mean point on top of it
      ggplot2::geom_errorbar(linewidth = 0.25, colour = "grey55") +
      ggplot2::geom_point(size = 0.75, colour = "blue", pch = 19) +
      # log transform the x-axis
      ggplot2::scale_x_continuous(transform = "log10", expand = expand_vals) +
      ggplot2::scale_y_continuous(expand = expand_vals) +
      ggplot2::labs(
        x = stringr::str_glue(
          "Mean number of cross-validation testing presences \\",
          "(log<sub>10</sub>scale)"),
        y = "Evaluation metric value", title = plot_title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.title = ggtext::element_markdown(
          colour = "red", face = "bold",
          margin = ggplot2::margin(-20, 0, -50, 0)),
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.major = ggplot2::element_line(
          colour = "lightgrey", linetype = 2, linewidth = 0.25))
  }

  # ****************************************************************
  # ****************************************************************

  plot_eval_metrics <- purrr::map(eval_metrics, plot_cv_metric) %>%
    # arrange plots in a single plot
    patchwork::wrap_plots(nrow = 2, ncol = 2, widths = 1, heights = 1) +
    # add a common title
    patchwork::plot_annotation(
      title = stringr::str_glue(
        "<b><span style='color:blue'>Cross-validated predictive power</span>\\
        </b><i><span style='color:grey'> --- {hab_type} </span></i>"),
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          size = 16, hjust = 0.5, margin = ggplot2::margin(1, 0, 0, 0)),
        plot.margin = ggplot2::margin(2, 0.5, -5, 0.5))) +
    patchwork::plot_layout(ncol = 2, axis_titles = "collect", axes = "collect")

  ragg::agg_jpeg(
    filename = fs::path(path_CV, "eval_metrics.jpeg"),
    width = 25, height = 23, res = 600, quality = 100, units = "cm")
  print(plot_eval_metrics)
  grDevices::dev.off()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting explanatory vs predictive power -----

  explained_vs_predictive <- function(metric, data = summary_all) {

    # Ensure metric is non-null, character strings with content
    if (is.null(metric) || !is.character(metric) || !nzchar(metric))  {
      ecokit::stop_ctx(
        "`metric` have to be character of length > 1", metric = metric,
        include_backtrace = TRUE)
    }
    # Ensure data is a data frame tibble with more than 1 row
    if (!is.data.frame(data) || nrow(data) < 1) {
      ecokit::stop_ctx(
        "`data` has to be a data frame tibble with more than 1 row",
        data = data, include_backtrace = TRUE)
    }

    # metric value has to be exist in the eval_metrics vector
    if (!metric %in% eval_metrics) {
      ecokit::stop_ctx(
        paste0(
          "`metric` should be one of the following: ", toString(eval_metrics)),
        metric = metric, include_backtrace = TRUE)
    }

    # column name for explained power
    exp_col <- paste0(metric, "_exp")
    # column names for mean and sd of predictive power
    pred_col <- paste0(metric, "_mean")
    sd_col <- paste0(metric, "_sd")

    data2 <- data %>%
      dplyr::select(dplyr::all_of(c(pred_col, exp_col, sd_col))) %>%
      stats::setNames(c("pred", "exp", "pred_sd")) %>%
      dplyr::mutate(
        mean_plus = pmin((pred + pred_sd), 1),
        mean_minus = pmax((pred - pred_sd), 0))

    # if the metric is AUC or TjurR2, ensure that mean_plus <= 1 and mean_minus
    # >= 0
    if (metric %in% c("AUC", "TjurR2")) {
      data2 <- dplyr::mutate(
        data2,
        mean_plus = pmin(mean_plus, 1), mean_minus = pmax(mean_minus, 0))
    }

    # if the metric is RMSE, ensure that the mean_minus >= 0
    if (metric == "RMSE") {
      data2 <- dplyr::mutate(data2, mean_minus = pmax(mean_minus, 0))
    }

    # if the metric is Boyce, ensure that mean_minus >=-1 and mean_plus <=1
    if (metric == "Boyce") {
      data2 <- dplyr::mutate(
        data2,
        mean_plus = pmin(mean_plus, 1), mean_minus = pmax(mean_minus, -1))
    }

    # range for x and y axes based on prediction bounds, excluding NA values
    range_xy <- c(data2$mean_plus, data2$mean_minus, data2$exp) %>%
      range(na.rm = TRUE)
    # expand value for x and y axes
    expand_vals <- c(0.0025, 0.0025)
    # Plot title
    plot_title <- dplyr::case_when(
      metric == "RMSE" ~ "Root mean square error (RMSE)",
      metric == "AUC" ~ "Area under the ROC curve (AUC)",
      metric == "Boyce" ~ "Continuous Boyce index",
      metric == "TjurR2" ~ "Tjur R-squared",
      .default = NULL)

    data2 %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = exp, y = pred, ymin = mean_minus, ymax = mean_plus)) +
      ggplot2::geom_errorbar(linewidth = 0.25, colour = "grey55") +
      ggplot2::geom_point(size = 0.75, colour = "blue", pch = 19) +
      ggplot2::geom_abline(
        intercept = 0, slope = 1, colour = "grey50", linetype = "dashed") +
      ggplot2::scale_x_continuous(expand = expand_vals, limits = range_xy) +
      ggplot2::scale_y_continuous(expand = expand_vals, limits = range_xy) +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::labs(x = NULL, y = NULL, title = plot_title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.major = ggplot2::element_line(
          colour = "lightgrey", linetype = 2, linewidth = 0.25))
  }

  # ****************************************************************
  # ****************************************************************

  plot_gp <- grid::gpar(col = "red", fontsize = 16)
  lab_y <- "<b>Predictive power</b> (spatia-block cross-validation)"
  lab_x <- "<b>Explanatory power</b> (model's internal evaluation)"

  plot_explained_vs_predictive <- eval_metrics %>%
    purrr::map(explained_vs_predictive) %>%
    # arrange plots in a single plot
    patchwork::wrap_plots(nrow = 2, ncol = 2, widths = 1, heights = 1) +
    # add a common title
    patchwork::plot_annotation(
      title = stringr::str_glue(
        "<b><span style='color:blue'>Explanatory <i>vs</i> predictive power\\
        </span></b><i><span style='color:grey'> --- {hab_type} </span></i>"),
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          size = 16, hjust = 0.5, margin = ggplot2::margin(5, 0, 2, 0)),
        plot.margin = ggplot2::margin(1, 0.5, 0, 0.5)))

  ragg::agg_jpeg(
    filename = fs::path(path_CV, "explained_vs_predictive.jpeg"),
    width = 25, height = 26.5, res = 600, quality = 100, units = "cm")
  plot_explained_vs_predictive %>%
    patchwork::patchworkGrob() %>%
    gridExtra::grid.arrange(
      left = gridtext::richtext_grob(
        lab_y, rot = 90, x = grid::unit(0.5, "cm"), gp = plot_gp),
      bottom = gridtext::richtext_grob(
        lab_x, y = grid::unit(0.4, "cm"), gp = plot_gp),
      padding = grid::unit(0.15, "cm"))
  grDevices::dev.off()

  # ****************************************************************

  ecokit::cat_diff(init_time = .start_time, prefix = "\nPost-processing took ")

  return(invisible(NULL))
}
