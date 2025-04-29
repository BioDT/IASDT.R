## |------------------------------------------------------------------------| #
# Model pipeline for post-processing cross-validated Hmsc models
## |------------------------------------------------------------------------| #

#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 5
#' @author Ahmed El-Gabbas
#' @export

mod_CV_postprocess_1_CPU <- function(
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

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "from_JSON", "use_TF", "TF_use_single", "LF_only", "LF_temp_cleanup",
      "LF_check", "temp_cleanup"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "env_file", "partition_name", "LF_runtime"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores", "n_batch_files"))
  rm(AllArgs, envir = environment())

  if (n_batch_files <= 0) {
    IASDT.R::stop_ctx(
      "`n_batch_files` must be a positive integer.",
      n_batch_files = n_batch_files)
  }
  if (n_cores <= 0) {
    IASDT.R::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores)
  }
  if (LF_n_cores <= 0) {
    IASDT.R::stop_ctx(
      "`LF_n_cores` must be a positive integer.", LF_n_cores = LF_n_cores)
  }

  if (!file.exists(env_file)) {
    IASDT.R::stop_ctx(
      "Error: Environment file is invalid or does not exist.",
      env_file = env_file)
  }

  if (!dir.exists(model_dir)) {
    IASDT.R::stop_ctx(
      "Model directory is invalid or does not exist.", model_dir = model_dir)
  }

  if (!all(CV_names %in% c("CV_Dist", "CV_Large", "CV_SAC"))) {
    IASDT.R::stop_ctx(
      paste0(
        "Invalid value for CV_names argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_names = CV_names)
  }

  # ****************************************************************

  # # Load environment variables, for project ID
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "ProjectID", "DP_R_LUMI_gpu", FALSE, FALSE)

  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # ****************************************************************

  # Paths to files and directories

  CV_dir <- IASDT.R::path(model_dir, "Model_Fitting_CV")

  Temp_dir <- IASDT.R::path(CV_dir, "Temp")
  # Path to store TF commands
  Path_TF <- IASDT.R::path(CV_dir, "LF_TensorFlow_commands")
  # Path to store log files
  Path_Log <- IASDT.R::path(Path_TF, "log")

  fs::dir_create(c(Path_TF, Path_Log))

  CV_DT_fitted <- IASDT.R::path(CV_dir, "CV_DT_fitted.RData")
  if (!file.exists(CV_DT_fitted)) {
    IASDT.R::stop_ctx(
      "CV_DT_fitted file not found.", CV_DT_fitted = CV_DT_fitted)
  }

  Path_LF_SLURM <- IASDT.R::path(Path_TF, "LF_SLURM.slurm")
  path_out <- IASDT.R::path(Path_Log, "%x-%A-%a.out")

  # ****************************************************************

  # Merge chains -----

  IASDT.R::info_chunk(
    "Merge chains", line_char = "|", line_char_rep = 60, cat_red = TRUE,
    cat_bold = TRUE, cat_timestamp = FALSE, info_lines_before = 2)

  IASDT.R::mod_merge_chains_CV(
    model_dir = model_dir, n_cores = n_cores, CV_names = CV_names,
    from_JSON = FALSE, out_extension = "qs2")

  # ****************************************************************

  # Prepare scripts for latent factor processing -----

  IASDT.R::info_chunk(
    "Prepare scripts for latent factor processing",
    line_char = "|", line_char_rep = 60, cat_red = TRUE, cat_bold = TRUE,
    cat_timestamp = FALSE, info_lines_before = 2)

  CV_DT_fitted <- IASDT.R::load_as(CV_DT_fitted) %>%
    dplyr::mutate(
      LF = purrr::map2(
        .x = CV_name, .y = CV,
        .f = ~ {

          IASDT.R::info_chunk(
            message = paste0(.x, "_", .y), level = 1, line_char = "+",
            line_char_rep = 60, cat_red = TRUE, cat_bold = TRUE)

          IASDT.R::predict_maps_CV(
            model_dir = model_dir, CV_name = paste0("CV_", .x),
            CV_fold = .y, n_cores = n_cores, use_TF = use_TF,
            TF_environ = TF_environ, TF_use_single = TF_use_single,
            LF_n_cores = LF_n_cores, LF_check = LF_check,
            LF_temp_cleanup = LF_temp_cleanup, LF_only = LF_only,
            LF_commands_only = TRUE, temp_cleanup = temp_cleanup)
        }
      ))

  # ****************************************************************

  # Merge TensorFlow commands into batch files -----

  IASDT.R::info_chunk(
    "Merge TensorFlow commands into batch files",
    line_char = "|", line_char_rep = 60, cat_red = TRUE, cat_bold = TRUE,
    info_lines_before = 2)

  IASDT.R::cat_time(
    paste0(
      "Merge and organize TensorFlow commands for LF predictions ",
      "into a maximum of ", n_batch_files, " files"),
    level = 1, cat_timestamp = FALSE)

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
    working_directory <- IASDT.R::normalize_path(
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
    IASDT.R::stop_ctx(
      "No command files were found in the temp directory",
      Temp_dir = Temp_dir, LF_InFiles = LF_InFiles,
      length_LF_InFiles = length(LF_InFiles))
  }

  IASDT.R::cat_time(
    paste0("Found ", length(LF_InFiles), " files"),
    level = 2, cat_timestamp = FALSE)
  stringr::str_remove_all(LF_InFiles, paste0(Temp_dir, "/")) %>%
    purrr::walk(IASDT.R::cat_time, level = 3, cat_timestamp = FALSE)

  # Read and merge commands from input files
  LF_commands <- purrr::map(LF_InFiles, readr::read_lines, progress = FALSE) %>%
    unlist() %>%
    gtools::mixedsort()

  IASDT.R::cat_time(
    paste0("Total number of commands to be executed: ", length(LF_commands)),
    level = 2, cat_timestamp = FALSE)

  if (length(LF_commands) < n_batch_files) {
    IASDT.R::cat_time(
      paste0(
        "Fewer commands than the requested number of files. ",
        "Setting `n_batch_files=", n_batch_files, "`."),
      level = 2, cat_timestamp = FALSE)
    n_batch_files <- length(LF_commands)
  }

  IASDT.R::cat_time(
    paste0("Splitting commands into ", n_batch_files, " files"),
    cat_timestamp = FALSE, level = 2)
  LF_commands <- IASDT.R::split_vector(LF_commands, n_splits = n_batch_files)
  line_sep <- strrep("-", 60)

  purrr::walk(
    .x = seq_len(length(LF_commands)),
    .f = ~ {

      Chunk_number <- stringr::str_pad(
        .x, pad = "0", width = nchar(n_batch_files))
      File <- IASDT.R::path(Path_TF, paste0("TF_Chunk_", Chunk_number, ".txt"))
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

  # ****************************************************************

  # Prepare LF batch file -----
  IASDT.R::cat_time("Prepare LF batch file")

  time_now <- lubridate::now(tzone = "CET") %>%
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

  IASDT.R::cat_time(
    paste0("Writing SLURM script to: `", Path_LF_SLURM, "`"),
    level = 2, cat_timestamp = FALSE)
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
# mod_CV_postprocess_2_CPU ----
## |------------------------------------------------------------------------| #

#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 6
#' @author Ahmed El-Gabbas
#' @export

mod_CV_postprocess_2_CPU <- function(
    model_dir = NULL, CV_names = NULL, n_cores = 8L, env_file = ".env",
    use_TF = TRUE, TF_use_single = FALSE, temp_cleanup = TRUE,
    LF_temp_cleanup = TRUE, TF_environ = NULL, LF_n_cores = n_cores,
    LF_check = FALSE) {

  # ****************************************************************

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  CV <- CV_name <- n_grids <- Sp <- IAS_ID <- metric <- y_label <- pred_mean <-
    pred <- pred_sd <- mean_minus <- mean_plus <- label <- summary_vals <-
    ias_id <- n_grids_pres_mean <- data <- title <- NULL

  # ****************************************************************

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("use_TF", "TF_use_single", "LF_temp_cleanup", "LF_check"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "env_file"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores"))
  rm(AllArgs, envir = environment())

  if (n_cores <= 0) {
    IASDT.R::stop_ctx(
      "`n_cores` must be a positive integer.", n_cores = n_cores)
  }
  if (LF_n_cores <= 0) {
    IASDT.R::stop_ctx(
      "`LF_n_cores` must be a positive integer.", LF_n_cores = LF_n_cores)
  }

  if (!all(CV_names %in% c("CV_Dist", "CV_Large", "CV_SAC"))) {
    IASDT.R::stop_ctx(
      paste0(
        "Invalid value for CV_names argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_names = CV_names)
  }

  # ****************************************************************

  # load model's explanatory power ------

  Eval_explain <- list.files(
    path = IASDT.R::path(model_dir, "Model_Evaluation"),
    pattern = "Eval_.+.qs2", full.names = TRUE)

  if (length(Eval_explain) != 1) {
    IASDT.R::stop_ctx(
      "There should be only one evaluation file in the Model_Evaluation folder",
      Eval_explain = Eval_explain)
  }

  Eval_explain <- IASDT.R::load_as(Eval_explain) %>%
    dplyr::select(-Sp) %>%
    dplyr::rename(ias_id = IAS_ID) %>%
    dplyr::filter(ias_id != "SR") %>%
    dplyr::rename_with(.cols = -ias_id, .fn = ~ paste0(.x, "_exp"))

  # ****************************************************************

  # Predicting habitat suitability at testing cross-validation folds -------

  IASDT.R::info_chunk(
    message = "Predicting habitat suitability at testing sites",
    line_char = "*", line_char_rep = 70, cat_red = TRUE, cat_bold = TRUE)

  path_CV_DT_fitted <- IASDT.R::path(
    model_dir, "Model_Fitting_CV", "CV_DT_fitted.RData")
  if (!file.exists(path_CV_DT_fitted)) {
    IASDT.R::stop_ctx(
      "path_CV_DT_fitted file not found.",
      path_CV_DT_fitted = path_CV_DT_fitted)
  }

  CV_DT_fitted <- IASDT.R::load_as(path_CV_DT_fitted)
  if (nrow(CV_DT_fitted) < 1) {
    IASDT.R::stop_ctx(
      "CV_DT_fitted data should contain at least one row",
      CV_DT_fitted = CV_DT_fitted)
  }

  CV_DT_fitted <- CV_DT_fitted %>%
    dplyr::mutate(
      Path_Preds_summary = purrr::map2_chr(
        .x = CV_name, .y = CV,
        .f = ~ {

          IASDT.R::info_chunk(
            message = paste0(.x, "_", .y), line_char = "~", level = 1,
            line_char_rep = 60, cat_red = TRUE, cat_bold = TRUE)

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

  IASDT.R::save_as(
    object = CV_DT_fitted, object_name = "CV_DT_fitted",
    out_path = path_CV_DT_fitted)

  invisible(gc())

  # ****************************************************************

  # Merge prediction summary data ------

  IASDT.R::cat_time("Merge prediction summary data")

  # ID columns to be used for nesting
  cols_id <- c(
    "hab_abb", "Hab_Name", "CV_name", "ias_id", "taxon_name", "species_name",
    "class", "order", "family")

  # Evaluation metrics to be used for summary
  eval_metrics <- c("RMSE", "AUC", "Boyce", "TjurR2")

  # Summary columns
  cols_summary <- outer(
    X = eval_metrics, Y = c("mean", "sd"),
    FUN = function(x, y) paste0(x, "_", y)) %>%
    as.vector() %>%
    sort() %>%
    c("n_grids_pres_mean")

  summary_all_CV <- purrr::map(
    CV_DT_fitted$Path_Preds_summary, IASDT.R::load_as) %>%
    dplyr::bind_rows() %>%
    tidyr::nest(.by = tidyselect::all_of(cols_id)) %>%
    dplyr::filter(ias_id != "SR") %>%
    dplyr::mutate(
      data = purrr::map(
        .x = data,
        .f = ~{
          .x %>%
            dplyr::arrange(CV_fold) %>%
            dplyr::summarise(
              dplyr::across(
                .cols = tidyselect::everything(),
                .fns = function(y) {
                  as.vector(y) %>%
                    stats::setNames(paste0("CV_", seq_len(length(.)))) %>%
                    list()
                })) %>%
            dplyr::mutate(
              dplyr::across(
                .cols = tidyselect::all_of(eval_metrics),
                list(mean = ~mean(unlist(.)), sd = ~sd(unlist(.))),
                .names = "{.col}_{.fn}"),
              n_grids_pres_mean = mean(unlist(n_grids_pres)))
        }),
      summary_vals = purrr::map(
        .x = data, .f = dplyr::select, tidyselect::all_of(cols_summary)),
      data = purrr::map(
        .x = data, .f = dplyr::select, -tidyselect::all_of(cols_summary))) %>%
    tidyr::unnest_wider(summary_vals) %>%
    dplyr::left_join(Eval_explain, by = "ias_id")

  # ****************************************************************

  # Plotting testing evaluation -----

  plot_cv_metric <- function(
    metric, y_label, title, data = summary_all_CV, linewidth = 0.25) {

    mean_col <- paste0(metric, "_mean")
    sd_col <- paste0(metric, "_sd")

    data %>%
      dplyr::select(n_grids_pres_mean, dplyr::all_of(c(mean_col, sd_col))) %>%
      stats::setNames(c("n_grids", "pred_mean", "pred_sd")) %>%
      dplyr::mutate(
        mean_plus = pred_mean + pred_sd,
        mean_minus = pred_mean - pred_sd) %>%
      ggplot2::ggplot(ggplot2::aes(x = n_grids, y = pred_mean)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mean_plus, ymax = mean_minus),
        linewidth = linewidth) +
      ggplot2::scale_x_continuous(transform = "log10") +
      ggplot2::labs(
        x = "Mean number of testing presences (log<sub>10</sub>scale)",
        y = y_label, title = title) +
      ggplot2::theme(axis.title = ggtext::element_markdown())
  }

  Plots <- tibble::tibble(
    metric = eval_metrics,
    y_label = c(
      "Root mean square error (RMSE)", "Area under the ROC curve (AUC)",
      "Continuous Boyce index", "Tjur R-Squared"),
    title = c("RMSE", "AUC", "continuous Boyce index", "Tjur R-Squared")) %>%
    dplyr::mutate(
      title = paste0("Cross-validated predictive power (", title, ")"),
      plot = purrr::pmap(
        .l = list(metric = metric, y_label = y_label, title = title),
        .f = plot_cv_metric))


  # Plots$plot[[1]]
  # Plots$plot[[2]]
  # Plots$plot[[3]]
  # Plots$plot[[4]]

  plot_pred_vs_exp <- function(
    metric, label, data = summary_all_CV, linewidth = 0.25) {

    pred_col <- paste0(metric, "_mean")
    exp_col <- paste0(metric, "_exp")
    sd_col <- paste0(metric, "_sd")

    data2 <- data %>%
      dplyr::select(dplyr::all_of(c(pred_col, exp_col, sd_col))) %>%
      stats::setNames(c("pred", "exp", "pred_sd")) %>%
      dplyr::mutate(
        mean_plus = pred + pred_sd,
        mean_minus = pred - pred_sd)

    range_xy <- range(c(data2$pred, data2$exp), na.rm = TRUE)

    data2 %>%
      ggplot2::ggplot(ggplot2::aes(x = exp, y = pred)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mean_minus, ymax = mean_plus),
        linewidth = linewidth) +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::coord_equal() +
      ggplot2::scale_x_continuous(limits = range_xy) +
      ggplot2::scale_y_continuous(limits = range_xy) +
      ggplot2::labs(
        x = paste0("Model's explanatory power (", label, ")"),
        y = paste0("Model's predictive power (", label, ")")) +
      ggplot2::theme_minimal()
  }

  plots2 <- tibble::tibble(
    metric = eval_metrics,
    label  = c("RMSE", "AUC", "Boyce index", "Tjur R-Squared")) %>%
    dplyr::mutate(
      plots = purrr::map2(.x = metric, .y = label, .f = plot_pred_vs_exp))

  # plots2$plots[[1]]
  # plots2$plots[[2]]
  # plots2$plots[[3]]
  # plots2$plots[[4]]

  # ****************************************************************

  IASDT.R::cat_diff(
    init_time = .StartTime, prefix = "\nPost-processing took ")

  return(invisible(NULL))
}
