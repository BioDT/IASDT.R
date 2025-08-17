## |------------------------------------------------------------------------| #
# Model pipeline for post-processing fitted Hmsc models
## |------------------------------------------------------------------------| #

#' Model pipeline for post-processing fitted Hmsc models
#'
#' These functions post-process fitted Hmsc models on both CPU and GPU. The main
#' functions in the pipeline includes `mod_postprocess_1_CPU`, `mod_prepare_TF`,
#' and `mod_postprocess_2_CPU` for full models without cross-validation, as well
#' as `mod_postprocess_CV_1_CPU` and `mod_postprocess_CV_2_CPU` for
#' cross-validated models. See details for more information.
#' @param model_dir Character. Path to the root directory of the fitted model.
#' @param job_runtime Character. Maximum allowed runtime for jobs for refitting
#'   the models (if needed) and cross validating models. Defaults to "01:00:00"
#'   for one hour. If not provided, the function throws an error.
#' @param GPP_dist Integer. Distance in *kilometres* between knots for the
#'   selected model.
#' @param use_trees Character. Whether a phylogenetic tree was used in the
#'   selected model. Accepts "Tree" (default) or "NoTree".
#' @param MCMC_thin,MCMC_n_samples Integer. Thinning value and the number of
#'   MCMC samples of the selected model.
#' @param n_cores_VP Integer. Number of cores to use for processing variance
#'   partitioning. Defaults to 10L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param CV_name `NULL` or character vector. Column name(s) in the model input
#'   data to be used to cross-validate the models (see [mod_prepare_data] and
#'   [mod_CV_prepare]). If `CV_name = NULL`, no cross-validation data
#'   preparation is done. See [mod_CV_fit] for valid options.
#' @param n_batch_files Integer. Number of output batch files to create. Must be
#'   less than or equal to the maximum job limit of the HPC environment.
#' @param working_directory Character. Optionally sets the working directory in
#'   batch scripts to this path. If `NULL`, the directory remains unchanged.
#' @param model_prefix Character. Prefix for the model name. A directory named
#'   `model_prefix_TF` is created in the `model_dir` to store the `TensorFlow`
#'   running commands. Defaults to `NULL`. This can not be `NULL`.
#' @param partition_name Character. Name of the partition to submit the SLURM
#'   jobs to. Default is `small-g`.
#' @param LF_runtime,VP_runtime Character. Time limit for latent factor
#'   prediction and variance partitioning processing jobs, respectively.
#'   Defaults are `01:00:00` and `02:00:00` respectively.
#' @param RC_n_cores Integer. The number of cores to use for response curve
#'   prediction. Defaults to `8`.
#' @param width_omega,height_omega,width_beta,height_beta Integer. The width and
#'   height of the generated heatmaps of the Omega and Beta parameters in
#'   centimetres.
#' @param spatial_model Logical. Whether the model is spatial (`TRUE`) or not
#'   (`FALSE`). Defaults to `TRUE`.
#' @param RC_prepare Logical. Whether to prepare the data for response curve
#'   prediction (using [resp_curv_prepare_data]). Defaults to `TRUE`.
#' @param RC_plot Logical. Whether to plot the response curves as JPEG files
#'   (using [resp_curv_plot_SR], [resp_curv_plot_species], and
#'   [resp_curv_plot_species_all]). Defaults to `TRUE`.
#' @param VP_prepare Logical. Whether to prepare the data for variance
#'   partitioning (using [variance_partitioning_compute]). Defaults to `TRUE`.
#' @param VP_plot Logical. Whether to plot the variance partitioning results
#'   (using [variance_partitioning_plot]). Defaults to `TRUE`.
#' @param predict_suitability Logical. Whether to predict habitat suitability
#'   across different climate options (using [predict_maps]). Defaults to
#'   `TRUE`.
#' @param plot_predictions Logical. Whether to plot species and species richness
#'   predictions as JPEG files (using [plot_prediction]). Defaults to `TRUE`.
#' @param plot_LF Logical. Whether to plot latent factors as JPEG files (using
#'   [plot_latent_factor]). Defaults to `TRUE`.
#' @param plot_internal_evaluation Logical. Whether to compute and visualise
#'   model internal evaluation (explanatory power) using [plot_evaluation].
#'   Defaults to `TRUE`.
#' @param process_VP Logical. Whether to prepares batch scripts for variance
#'   partitioning GPU computations on GPUs. Defaults to `TRUE`.
#' @param process_LF Logical. Whether to prepares batch scripts for latent
#'   factor predictions GPU computations on GPUs. Defaults to `TRUE`.
#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 1
#' @importFrom foreach %dopar%
#' @inheritParams predict_maps
#' @inheritParams mod_CV_fit
#' @inheritParams mod_merge_chains
#' @inheritParams mod_inputs
#' @inheritParams response_curves
#' @inheritParams coda_to_tibble
#' @inheritParams mod_merge_chains
#' @author Ahmed El-Gabbas
#' @export
#' @details
#'
#' **mod_postprocess_1_CPU**
#'
#' This function performs the initial post-processing step for habitat-specific
#' fitted models, automating the following tasks:
#'
#' - check unsuccessful models: [mod_SLURM_refit]
#' - merge chains and save R objects (fitted model object and coda object) to
#' `qs2` or `RData` files: [mod_merge_chains]
#' - visualise the convergence of all model variants fitted
#' [convergence_plot_all]
#' - visualise the convergence of selected model, including plotting
#' Gelman-Rubin-Brooks [plot_gelman] and [convergence_plot] for model
#' convergence diagnostics of the `rho`, `alpha`, `omega`, and `beta`
#' parameters.
#' - extract and save model summary: [mod_summary]
#' - plotting model parameters: [mod_heatmap_omega], [mod_heatmap_beta]
#' - prepare data for cross-validation and fit initial cross-validated models:
#' [mod_CV_fit]
#' - Prepare scripts for GPU processing, including:
#'    - predicting latent factors of the response curves:
#' [resp_curv_prepare_data]
#'    - predicting latent factors for new sampling units: [predict_maps]
#'    - computing variance partitioning: [variance_partitioning_compute]
#'
#' <br/>
#'
#' **mod_prepare_TF**
#'
#' After running `mod_postprocess_1_CPU` for all habitat types, this function
#' prepares batch scripts for GPU computations of all habitat types:
#'    - for <u>variance partitioning</u>, the function matches all files with
#' the pattern ` "VP_.+Command.txt"` (created by [variance_partitioning_compute]
#' and merges their contents into a single file
#' (`model_prefix_TF/VP_Commands.txt`). Then, it prepares a SLURM script for
#' variance partitioning computations (`model_prefix_TF/VP_SLURM.slurm`).
#'    - for <u>latent factor predictions</u>, the function matches all files
#' with the pattern `"^LF_NewSites_Commands_.+.txt|^LF_RC_Commands_.+txt"` and
#' split their contents into multiple scripts at the `model_prefix_TF` directory
#' for processing as a batch job. The function prepares a SLURM script for
#' latent factor predictions (`LF_SLURM.slurm`).
#'
#' This function is tailored for the LUMI HPC environment and assumes that the
#' `tensorflow` module is installed and correctly configured with all required
#' Python packages. On other HPC systems, users may need to modify the function
#' to load a Python virtual environment or install the required dependencies for
#' `TensorFlow` and related packages.
#'
#' <br/><br/>
#'
#' **mod_postprocess_2_CPU**
#'
#' This function continues running the analysis pipeline for post-processing
#' Hmsc by automating the following steps:
#' - process and visualise response curves: [response_curves]
#' - predict habitat suitability across different climate options:
#' [predict_maps]
#' - plot species & SR predictions as JPEG: [plot_prediction]
#' - plot latent factors as JPEG: [plot_latent_factor]
#' - process and visualise variance partitioning:
#' [variance_partitioning_compute] and [variance_partitioning_plot]
#' - compute and visualizing model internal evaluation (explanatory power):
#' [plot_evaluation]
#'  - initiate post-processing of fitted cross-validated models: prepare
#' commands for latent factor predictions on GPU --- **Ongoing**
#'
#' This function should be run after:
#' - completing `mod_postprocess_1_CPU` and `mod_prepare_TF` on CPU,
#' - running `VP_SLURM.slurm` and `LF_SLURM.slurm` on GPU to process response
#' curves and latent factor predictions (both scripts are generated by
#' `mod_prepare_TF`).
#' - submitting SLURM jobs for cross-validated model fitting.
#'
#' <br/>
#'
#' **mod_postprocess_CV_1_CPU**
#'
#' This function is similar to `mod_postprocess_1_CPU`, but it is specifically
#' designed for cross-validated models. It automates merging fitted
#' cross-validated model chains into `Hmsc` model objects and prepare scripts
#' for latent factor prediction on `TensorFlow` using [predict_maps_CV].
#'
#' <br/><br/>
#'
#' **mod_postprocess_CV_2_CPU**
#'
#' The function 1) processes `*.feather` files resulted from Latent Factor
#' predictions (using `TensorFlow`) and saves LF predication to disk; 2)
#' predicts species-specific mean habitat suitability at testing
#' cross-validation folds and calculates testing evaluation metrics; 3)
#' generates plots of the evaluation metrics.
#'

## |------------------------------------------------------------------------| #
# mod_postprocess_1_CPU ----
## |------------------------------------------------------------------------| #

mod_postprocess_1_CPU <- function(
    model_dir = NULL, hab_abb = NULL, n_cores = 8L, strategy = "multisession",
    env_file = ".env", path_Hmsc = NULL, memory_per_cpu = "64G",
    job_runtime = "01:00:00", from_JSON = FALSE, GPP_dist = NULL,
    use_trees = "Tree", MCMC_n_samples = 1000L, MCMC_thin = NULL,
    n_omega = 1000L, CV_name = c("CV_Dist", "CV_Large"), n_grid = 50L,
    use_TF = TRUE, TF_use_single = FALSE, LF_n_cores = n_cores,
    LF_temp_cleanup = TRUE, LF_check = FALSE, temp_cleanup = TRUE,
    TF_environ = NULL, pred_new_sites = TRUE, n_cores_VP = 10L,
    width_omega = 26, height_omega = 22.5, width_beta = 25, height_beta = 35,
    spatial_model = TRUE, clamp_pred = TRUE, fix_efforts = "q90",
    fix_rivers = "q90",
    CC_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0",
      "UKESM1-0-LL"),
    CC_scenario = c("ssp126", "ssp370", "ssp585")) {

  .start_time <- lubridate::now(tzone = "CET")

  species_name <- non_focal_variables <- NULL

  # ****************************************************************

  # Check input arguments ----
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "use_trees", "path_Hmsc"))

  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c(
      "from_JSON", "pred_new_sites", "temp_cleanup",
      "use_TF", "LF_check", "LF_temp_cleanup", "TF_use_single"))

  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c(
      "n_cores", "n_omega", "GPP_dist", "MCMC_n_samples",
      "MCMC_thin", "n_cores_VP", "LF_n_cores", "n_grid"))
  rm(AllArgs, envir = environment())

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") {
    n_cores <- LF_n_cores <- n_cores_VP <- 1L
  }
  n_cores <- .validate_n_cores(n_cores)
  LF_n_cores <- .validate_n_cores(LF_n_cores)
  n_cores_VP <- .validate_n_cores(n_cores_VP)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory was not found", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!(use_trees %in% c("Tree", "NoTree"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for Tree argument. Valid values ",
        "are: 'Tree' or 'NoTree'"),
      use_trees = use_trees, include_backtrace = TRUE)
  }

  if (!is.null(CV_name) &&
      !all(CV_name %in% c("CV_Dist", "CV_Large", "CV_SAC"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for CV_name argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_name = CV_name, include_backtrace = TRUE)
  }

  # Validate memory_per_cpu
  memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
  job_runtime <- .validate_slurm_runtime(job_runtime)

  # ****************************************************************

  ecokit::record_arguments(
    out_path = fs::path(model_dir, "Args_mod_postprocess_1_CPU.RData"))

  # ****************************************************************

  loaded_packages <- paste(
    sort(ecokit::loaded_packages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "\n  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", ecokit::os(),
      "\n  >>> Model root: ", model_dir,
      "\n  >>> n_cores: ", n_cores,
      "\n  >>> env_file: ", env_file,
      "\n  >>> path_Hmsc: ", path_Hmsc,
      "\n  >>> SLURM memory_per_cpu: ", memory_per_cpu,
      "\n  >>> SLURM Time: ", job_runtime,
      "\n  >>> n_omega: ", n_omega,
      "\n  >>> hab_abb: ", hab_abb,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"), # nolint: undesirable_function_linter
      "\n  >>> Loaded packages: \n", loaded_packages, "\n"))

  # ****************************************************************

  # Check unsuccessful models -----

  ecokit::info_chunk(
    "Check unsuccessful models", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_SLURM_refit(
    model_dir = model_dir,
    job_name = stringr::str_remove(basename(model_dir), "Mod_"),
    memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
    env_file = env_file, path_Hmsc = path_Hmsc)

  invisible(gc())

  # ****************************************************************

  ecokit::info_chunk(
    "Merge chains and saving RData or qs2 files", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_merge_chains(
    model_dir = model_dir, n_cores = n_cores, strategy = strategy,
    from_JSON = from_JSON)

  invisible(gc())

  # ****************************************************************

  ecokit::info_chunk(
    "Convergence of all model variants", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::convergence_plot_all(
    model_dir = model_dir, n_omega = n_omega, n_cores = n_cores,
    strategy = strategy, margin_type = "histogram")

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Path of selected model -----
  ecokit::info_chunk(
    "Path of selected model", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  path_model <- fs::dir_ls(
    fs::path(model_dir, "Model_Fitted"),
    regexp = paste0(
      ".+", use_trees, "_samp", MCMC_n_samples, "_th", MCMC_thin, "_Model.qs2"))
  if (length(path_model) > 1) {
    path_model <- stringr::str_subset(path_model, paste0("^GPP", GPP_dist))
  }
  if (length(path_model) == 0) {
    ecokit::stop_ctx("`path_model` is of length 0")
  }

  path_coda <- fs::dir_ls(
    fs::path(model_dir, "Model_Coda"),
    regexp = paste0(
      ".+", use_trees, "_samp", MCMC_n_samples, "_th", MCMC_thin, "_Coda.qs2"))
  if (length(path_coda) > 1) {
    path_coda <- stringr::str_subset(
      path_coda, paste0("^GPP", GPP_dist))
  }
  if (length(path_coda) == 0) {
    ecokit::stop_ctx("`path_coda` is of length 0")
  }

  cat(
    paste0("path_model:\n\t", path_model, "\nPath_Coda:\n\t", path_coda, "\n"))

  if (!all(fs::file_exists(c(path_model, path_coda)))) {
    ecokit::stop_ctx(
      "Selected model files not found",
      path_model = path_model, path_coda = path_coda, include_backtrace = TRUE)
  }

  temp_dir <- fs::path(model_dir, "TEMP_Pred")

  # ****************************************************************

  # Gelman_Plot -----
  ecokit::info_chunk(
    "Gelman_Plot", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::plot_gelman(
    path_coda = path_coda, alpha = TRUE, beta = TRUE, omega = TRUE, rho = TRUE,
    n_omega = n_omega, env_file = env_file)

  invisible(gc())

  # ****************************************************************

  # Convergence plots ----
  ecokit::info_chunk(
    "Convergence plots", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::convergence_plot(
    path_coda = path_coda, path_model = path_model, env_file = env_file,
    n_omega = n_omega, n_cores = n_cores, strategy = strategy, n_RC = c(2, 2),
    beta_n_RC = c(3, 3), margin_type = "histogram")

  invisible(gc())

  # ****************************************************************

  # Model summary ------
  ecokit::info_chunk(
    "Model summary", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_summary(path_coda = path_coda, env_file = env_file)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Omega -----
  ecokit::info_chunk(
    "Plotting model parameters - Omega", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_heatmap_omega(
    path_model = path_model, support_level = 0.95,
    width = width_omega, height = height_omega)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Beta -----
  ecokit::info_chunk(
    "Plotting model parameters - Beta", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_heatmap_beta(
    path_model = path_model, support_level = 0.95,
    width = width_beta, height = height_beta)

  invisible(gc())

  # ****************************************************************

  # Prepare input data for cross-validation -------

  if (!is.null(CV_name)) {
    ecokit::info_chunk(
      "Prepare input data for cross-validation", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::mod_CV_fit(
      path_model = path_model, CV_name = CV_name, env_file = env_file,
      job_name = paste0("CV_", hab_abb), memory_per_cpu = memory_per_cpu,
      job_runtime = job_runtime, path_Hmsc = path_Hmsc)

    invisible(gc())
  }

  # ****************************************************************

  if (spatial_model) {

    # latent factors of the response curves -------
    ecokit::info_chunk(
      "Prepare scripts for predicting latent factors of the response curves",
      line_char = "+", line_char_rep = 90L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)

    IASDT.R::resp_curv_prepare_data(
      path_model = path_model, n_grid = n_grid, n_cores = n_cores,
      strategy = strategy, use_TF = use_TF, TF_environ = TF_environ,
      TF_use_single = TF_use_single, LF_n_cores = LF_n_cores,
      LF_temp_cleanup = LF_temp_cleanup, LF_check = LF_check,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, verbose = TRUE,
      LF_commands_only = TRUE)

    invisible(gc())

    # ****************************************************************

    # latent factors for new sampling units -------

    ecokit::info_chunk(
      "Prepare scripts for predicting latent factors for new sampling units",
      line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)

    IASDT.R::predict_maps(
      path_model = path_model, hab_abb = hab_abb, env_file = env_file,
      n_cores = n_cores, strategy = strategy, clamp_pred = FALSE,
      pred_new_sites = pred_new_sites, use_TF = use_TF, TF_environ = TF_environ,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup,
      TF_use_single = TF_use_single, LF_n_cores = LF_n_cores,
      LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup, LF_only = TRUE,
      LF_commands_only = TRUE, spatial_model = spatial_model,
      fix_efforts = fix_efforts, fix_rivers = fix_rivers, CC_models = CC_models,
      CC_scenario = CC_scenario)

  } else {

    # Response curves -------

    ecokit::info_chunk(
      "Prepare response curve data", line_char = "+", line_char_rep = 90L,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1L)

    rc_dir <- fs::path(
      dirname(dirname(path_model)), "Model_Postprocessing", "RespCurv_DT")
    fs::dir_create(rc_dir)
    path_rc_data <- fs::path(rc_dir, "response_curve_data.qs2")
    path_observed_data <- fs::path(rc_dir, "observed_data.qs2")

    if (!ecokit::check_data(path_observed_data, warning = FALSE)) {
      Model <- ecokit::load_as(path_model)
      model_vars <- names(Model$XData)
      observed_data <- dplyr::bind_cols(
        x_value = Model$XData, Model$Y,
        species_richness = rowSums(Model$Y)) %>%
        tidyr::pivot_longer(
          cols = -tidyselect::all_of(model_vars), names_to = "species_name",
          values_to = "observed_value") %>%
        dplyr::select("species_name", tidyselect::everything()) %>%
        dplyr::arrange(species_name)

      ecokit::save_as(object = observed_data, out_path = path_observed_data)
      rm(Model, observed_data, envir = environment())
      invisible(gc())
    }

    if (!ecokit::check_data(path_rc_data, warning = FALSE)) {
      Model <- ecokit::load_as(path_model)
      model_vars <- names(Model$XData)
      rm(Model, envir = environment())

      response_curve_data <- tidyr::expand_grid(
        variable = model_vars, non_focal_variables = c(1, 2))

      # Prepare working in parallel
      if (n_cores == 1) {
        future::plan("sequential", gc = TRUE)
      } else {
        ecokit::set_parallel(
          n_cores = min(n_cores, nrow(response_curve_data)), level = 2L,
          strategy = strategy, cat_timestamp = FALSE, future_max_size = 1500L)
        withr::defer(future::plan("sequential", gc = TRUE))
      }

      response_curve_data2 <- future.apply::future_lapply(
        X = seq_len(nrow(response_curve_data)),
        FUN = function(x) {

          Model <- ecokit::load_as(path_model)
          variable <- response_curve_data$variable[x]
          nfv <- response_curve_data$non_focal_variables[x]

          Gradient <- Hmsc::constructGradient(
            hM = Model, focalVariable = variable,
            non.focalVariables = nfv, ngrid  = 50)
          pred <- stats::predict(
            object = Model, XData = Gradient$XDataNew,
            studyDesign = Gradient$studyDesignNew,
            ranLevels = Gradient$rLNew, expected = TRUE) %>%
            abind::abind(along = 3)
          x_vals <- Gradient$XDataNew[, variable]
          rm(Model, Gradient, envir = environment())

          pred_mean <- apply(pred, c(1, 2), mean)
          rich_array <- apply(pred, 3, rowSums)
          rich_mean <- rowMeans(rich_array)
          pred_mean <- cbind(pred_mean, species_richness = rich_mean)

          pred_sd <- apply(pred, c(1, 2), stats::sd)
          rich_sd <- apply(rich_array, 1, stats::sd)
          pred_sd <- cbind(pred_sd, species_richness = rich_sd)

          pred_mean_plus  <- pred_mean + pred_sd
          pred_mean_minus <- pred_mean - pred_sd

          pred_mean <- tibble::as_tibble(pred_mean) %>%
            dplyr::mutate(x_value = x_vals, .before = 1) %>%
            tidyr::pivot_longer(
              cols = -"x_value", names_to = "species_name",
              values_to = "mean")
          pred_mean_plus <- tibble::as_tibble(pred_mean_plus) %>%
            dplyr::mutate(x_value = x_vals, .before = 1) %>%
            tidyr::pivot_longer(
              cols = -"x_value", names_to = "species_name",
              values_to = "mean_plus_sd")
          pred_mean_minus <- tibble::as_tibble(pred_mean_minus) %>%
            dplyr::mutate(x_value = x_vals, .before = 1) %>%
            tidyr::pivot_longer(
              cols = -"x_value", names_to = "species_name",
              values_to = "mean_minus_sd")

          col_order <- c(
            "species_name", "variable", "non_focal_variables", "x_value")
          rc_data <- dplyr::left_join(
            pred_mean, pred_mean_plus, by = c("x_value", "species_name")) %>%
            dplyr::left_join(
              pred_mean_minus, by = c("x_value", "species_name")) %>%
            dplyr::mutate(
              variable = variable, non_focal_variables = as.integer(nfv)) %>%
            dplyr::select(
              tidyselect::all_of(col_order), tidyselect::everything()) %>%
            dplyr::arrange(species_name, variable, non_focal_variables)

          rc_data

        },
        future.seed = TRUE,
        future.globals = c("response_curve_data", "path_model"),
        future.packages = c(
          "Hmsc", "dplyr", "magrittr", "Hmsc", "abind", "tibble",
          "tidyr", "ecokit", "tidyselect"))

      if (n_cores > 1) {
        ecokit::set_parallel(
          stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
        future::plan("sequential", gc = TRUE)
      }

      response_curve_data <- dplyr::bind_rows(response_curve_data2)
      ecokit::save_as(object = response_curve_data, out_path = path_rc_data)
      rm(response_curve_data, envir = environment())
      invisible(gc())
    }

    # Predictions -------

    ecokit::info_chunk(
      "Making predictions", line_char = "+", line_char_rep = 60L,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1L)

    IASDT.R::predict_maps(
      path_model = path_model, hab_abb = hab_abb, env_file = env_file,
      n_cores = n_cores, strategy = strategy, clamp_pred = clamp_pred,
      pred_new_sites = pred_new_sites, use_TF = use_TF, TF_environ = TF_environ,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup,
      TF_use_single = TF_use_single, LF_n_cores = LF_n_cores,
      LF_check = LF_check, LF_temp_cleanup = LF_temp_cleanup, LF_only = FALSE,
      LF_commands_only = FALSE, spatial_model = TRUE)
  }

  # ****************************************************************

  # Prepare scripts for computing variance partitioning -------
  ecokit::info_chunk(
    "Prepare scripts for computing variance partitioning", level = 1L,
    line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
    cat_timestamp = FALSE)

  IASDT.R::variance_partitioning_compute(
    path_model = path_model, n_cores = n_cores_VP, use_TF = use_TF,
    TF_environ = TF_environ, TF_use_single = TF_use_single,
    temp_cleanup = temp_cleanup, chunk_size = 50L, verbose = TRUE,
    VP_file = "VarPar", VP_commands_only = TRUE)

  # ****************************************************************

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nPost-processing using CPU took ")

  return(invisible(NULL))
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_prepare_TF ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 2
#' @author Ahmed El-Gabbas

mod_prepare_TF <- function(
    process_VP = TRUE, process_LF = TRUE,
    n_batch_files = 210L, env_file = ".env", working_directory = NULL,
    partition_name = "small-g", LF_runtime = "01:00:00", model_prefix = NULL,
    VP_runtime = "02:00:00") {

  # ****************************************************************

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjectID <- path_model <- NULL

  # ****************************************************************

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c(
      "LF_runtime", "VP_runtime", "partition_name", "env_file", "model_prefix"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_batch_files")
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("process_VP", "process_LF"))
  rm(AllArgs, envir = environment())

  if (n_batch_files <= 0) {
    ecokit::stop_ctx(
      "`n_batch_files` must be a positive integer.",
      n_batch_files = n_batch_files, include_backtrace = TRUE)
  }

  # # Load environment variables, for project ID

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "ProjectID", "DP_R_LUMI_gpu", FALSE, FALSE,
    "path_model", "DP_R_Model_path", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # ****************************************************************

  # Path to store TF commands
  Path_TF <- fs::path(path_model, paste0(model_prefix, "_TF"))
  # Path to store log files
  Path_Log <- ecokit::normalize_path(fs::path(Path_TF, "log"))
  fs::dir_create(c(Path_TF, Path_Log))

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  # VARIANCE PARTITIONING ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  if (process_VP) {

    # Prepare post-processing data for calculating variance partitioning
    ecokit::cat_time(
      "Prepare post-processing data for calculating variance partitioning")

    # Create paths for VP SLURM script and commands
    Path_VP_SLURM <- fs::path(Path_TF, "VP_SLURM.slurm")
    Path_VP_Commands <- fs::path(Path_TF, "VP_Commands.txt")

    # Merge and organise `TensorFlow` commands for computing variance
    # partitioning
    ecokit::cat_time(
      "Merge and organise `TensorFlow` for computing variance partitioning",
      level = 1L, cat_timestamp = FALSE)

    # Find list of files matching the pattern
    VP_InFiles <- fs::dir_ls(
      path = path_model, recurse = TRUE, type = "file",
      regexp = paste0(model_prefix, ".+/TEMP_VP/VP_.+Command.txt")) %>%
      purrr::map(readr::read_lines, progress = FALSE) %>%
      unlist() %>%
      gtools::mixedsort() %>%
      unname()

    n_VP_InFiles <- length(VP_InFiles)

    # Save all VP commands to single file for batch processing
    readr::write_lines(x = VP_InFiles, file = Path_VP_Commands, append = FALSE)

    # ****************************************************************

    # Prepare batch files ----
    ecokit::cat_time("Prepare batch files", level = 1L, cat_timestamp = FALSE)

    VP_Commands <- c(
      "#!/bin/bash",
      "#SBATCH --job-name=VP_TF",
      "#SBATCH --ntasks=1",
      paste0("#SBATCH --account=", ProjectID),
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --gpus-per-node=1",
      paste0("#SBATCH --time=", VP_runtime),
      paste0("#SBATCH --partition=", partition_name),
      "#SBATCH --mem=64G",
      paste0("#SBATCH --output=", fs::path(Path_Log, "%x-%A-%a.out")),
      paste0("#SBATCH --error=", fs::path(Path_Log, "%x-%A-%a.out")),
      paste0("#SBATCH --array=1-", length(VP_InFiles)),
      "\n# File containing commands to be executed",
      paste0("File=", Path_VP_Commands),
      "\n# Load TensorFlow module and configure environment",
      "ml use /appl/local/csc/modulefiles",
      "ml tensorflow",
      "export TF_CPP_MIN_LOG_LEVEL=3",
      "export TF_ENABLE_ONEDNN_OPTS=0\n",
      "# Verify GPU availability",
      paste0(
        'python3 -c "import tensorflow as tf; ',
        'print(\\\"Num GPUs Available:\\\", ',
        'len(tf.config.list_physical_devices(\\\"GPU\\\")))"'),
      "\n# Run array job",
      "head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash",
      "\necho End of program at `date`\n",
      paste0("# ", strrep("-", 50)),
      paste0(
        "# This script was created on: ",
        format(lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"),
        " CET"),
      paste0("# ", strrep("-", 50)))

    ecokit::cat_time(
      paste0("Writing SLURM script to: `", Path_VP_SLURM, "`"),
      level = 1L, cat_timestamp = FALSE)

    readr::write_lines(x = VP_Commands, file = Path_VP_SLURM, append = FALSE)
    # Make the file executable
    Sys.chmod(Path_VP_SLURM, mode = "755")

    ecokit::cat_time(
      paste0(
        "To submit variance partitioning SLURM script:\n\tsbatch ",
        Path_VP_SLURM),
      cat_timestamp = FALSE)

  } else {

    n_VP_InFiles <- 0L

  }

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  # LF PREDICTIONS ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  if (process_LF) {

    # Prepare post-processing data for LF predictions
    ecokit::cat_time("\nPrepare post-processing data for LF predictions")

    # Ensure that the total number of simultaneous jobs (LF + VP) =
    # n_batch_files; so both can be run on the same time.
    n_batch_files <- n_batch_files - n_VP_InFiles

    # Merge and organise `TensorFlow` commands for LF predictions ----
    ecokit::cat_time(
      paste0(
        "Merge and organise `TensorFlow` commands for LF predictions ",
        "into a maximum of ", n_batch_files, " files"),
      level = 1L, cat_timestamp = FALSE)

    # Basic commands for `TensorFlow` setup
    BasicCommands <- c(
      "#!/bin/bash\n",
      "# Load TensorFlow module and configure environment",
      "ml use /appl/local/csc/modulefiles",
      "ml tensorflow",
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
    # Regex pattern to match input files
    LF_Pattern <- "(LF_NewSites_Commands_.+.txt|LF_RC_Commands_.+txt)"
    LF_InFiles <- fs::dir_ls(
      path = path_model, recurse = TRUE, type = "file",
      regexp = paste0(model_prefix, ".+/TEMP_Pred/", LF_Pattern))

    if (length(LF_InFiles) == 0) {
      ecokit::stop_ctx(
        "No files found matching the pattern",
        LF_Pattern = LF_Pattern, path_model = path_model,
        include_backtrace = TRUE)
    }

    ecokit::cat_time(
      paste0(
        "Found ", length(LF_InFiles), " files matching the pattern `",
        LF_Pattern, "`"),
      level = 1L, cat_timestamp = FALSE)
    purrr::walk(LF_InFiles, ecokit::cat_time, level = 2L, cat_timestamp = FALSE)

    # Read and merge commands from input files
    LF_commands <- purrr::map(
      .x = LF_InFiles, .f = readr::read_lines, progress = FALSE) %>%
      unlist() %>%
      gtools::mixedsort() %>%
      unname()

    ecokit::cat_time(
      paste0(
        "Total number of commands to be executed: ", length(LF_commands)),
      level = 1L, cat_timestamp = FALSE)

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
      cat_timestamp = FALSE, level = 1L)
    LF_commands <- ecokit::split_vector(LF_commands, n_splits = n_batch_files)

    purrr::walk(
      .x = seq_len(length(LF_commands)),
      .f = ~ {

        File <- fs::path(
          Path_TF,
          paste0(
            "TF_Chunk_",
            stringr::str_pad(
              .x, width = nchar(n_batch_files), pad = "0"), ".txt"))

        readr::write_lines(x = BasicCommands, file = File, append = FALSE)
        readr::write_lines(
          x = paste0(
            "# ", length(LF_commands[[.x]]), " commands to be executed:"),
          file = File, append = TRUE)
        readr::write_lines(x = LF_commands[[.x]], file = File, append = TRUE)
        readr::write_lines(
          x = c(
            paste0("\n#", strrep("_", 60)),
            paste0(
              "# This script was created on: ",
              format(lubridate::now(tzone = "CET"), "%Y-%m-%d %H:%M:%S")),
            paste0("#", strrep("_", 60))),
          file = File, append = TRUE)

        return(invisible(NULL))
      })

    # ****************************************************************

    # Prepare LF batch files ----
    ecokit::cat_time("Prepare batch files", level = 1L, cat_timestamp = FALSE)

    LF_slurm_script <- c(
      "#!/bin/bash",
      "#SBATCH --job-name=PP_LF",
      "#SBATCH --ntasks=1",
      paste0("#SBATCH --account=", ProjectID),
      "#SBATCH --mem=64G",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --gpus-per-node=1",
      paste0("#SBATCH --time=", LF_runtime),
      paste0("#SBATCH --partition=", partition_name),
      paste0("#SBATCH --output=", fs::path(Path_Log, "%x-%A-%a.out")),
      paste0("#SBATCH --error=", fs::path(Path_Log, "%x-%A-%a.out")),
      paste0("#SBATCH --array=1-", n_batch_files),
      "",
      "# Define directories",
      paste0('OutputDir="', Path_TF, '"'),
      "",
      "# Find all the split files and sort them explicitly",
      paste0(
        'SplitFiles=($(find "$OutputDir" -type f ',
        '-name "TF_Chunk_*.txt" | sort -V))'),
      "",
      "# Check if files were found",
      "if [ ${#SplitFiles[@]} -eq 0 ]; then",
      '    echo "Error: No files matching TF_Chunk_*.txt found in $OutputDir"',
      "    exit 1",
      "fi",
      "",
      paste0("# Ensure no more than `, n_batch_files, ` files are processed"),
      paste0("MaxFiles=", n_batch_files),
      "if [ ${#SplitFiles[@]} -gt $MaxFiles ]; then",
      '    SplitFiles=("${SplitFiles[@]:0:$MaxFiles}")',
      paste0(
        '    echo "More than $MaxFiles files found, ',
        'limiting to the first $MaxFiles files."'),
      "fi",
      "",
      "# Get the index of the current task based on SLURM_ARRAY_TASK_ID",
      "TaskIndex=$((SLURM_ARRAY_TASK_ID - 1))",
      "",
      "# Validate TaskIndex",
      "if [ $TaskIndex -ge ${#SplitFiles[@]} ] || [ $TaskIndex -lt 0 ]; then",
      paste0(
        '    echo "Error: TaskIndex $TaskIndex is out of range. ',
        'Valid range: 0 to $((${#SplitFiles[@]} - 1))"'),
      "    exit 1",
      "fi",
      "",
      "# Get the specific split file to process based on the job array task ID",
      'SplitFile="${SplitFiles[$TaskIndex]}"',
      "",
      "# Verify the selected split file",
      'if [ -z "$SplitFile" ] || [ ! -f "$SplitFile" ]; then',
      '    echo "Error: File $SplitFile does not exist or is invalid."',
      "    exit 1",
      "fi",
      "",
      "# Processing file",
      'echo "Processing file: $SplitFile"',
      "",
      "# Run the selected split file",
      'bash "$SplitFile"',
      "\necho End of program at `date`\n",
      paste0("# ", strrep("-", 50)),
      paste0(
        "# This script was created on: ",
        format(
          lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"), " CET"),
      paste0("# ", strrep("-", 50)))

    Path_LF_SLURM <- fs::path(Path_TF, "LF_SLURM.slurm")
    ecokit::cat_time(
      paste0("Writing SLURM script to: `", Path_LF_SLURM, "`"),
      level = 1L, cat_timestamp = FALSE)
    # Write the content to a file
    readr::write_lines(LF_slurm_script, Path_LF_SLURM, append = FALSE)
    # Make the file executable
    Sys.chmod(Path_LF_SLURM, mode = "755")

    ecokit::cat_time(
      paste0(
        "\nTo submit LF prediction SLURM script:\n\tsbatch ", Path_LF_SLURM),
      cat_timestamp = FALSE)

  }

  # ****************************************************************
  # ****************************************************************

  return(invisible(NULL))
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_postprocess_2_CPU ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 3
#' @author Ahmed El-Gabbas

mod_postprocess_2_CPU <- function(
    model_dir = NULL, hab_abb = NULL, n_cores = 8L,
    strategy = "multisession", env_file = ".env", GPP_dist = NULL,
    use_trees = "Tree", MCMC_n_samples = 1000L, MCMC_thin = NULL, use_TF = TRUE,
    TF_environ = NULL, TF_use_single = FALSE, LF_n_cores = n_cores,
    LF_check = FALSE, LF_temp_cleanup = TRUE, temp_cleanup = TRUE, n_grid = 50L,
    CC_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_scenario = c("ssp126", "ssp370", "ssp585"),
    RC_n_cores = 8L, clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    pred_new_sites = TRUE, tar_predictions = TRUE,
    RC_prepare = TRUE, RC_plot = TRUE, VP_prepare = TRUE, VP_plot = TRUE,
    predict_suitability = TRUE, plot_predictions = TRUE, plot_LF = TRUE,
    plot_internal_evaluation = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  # ****************************************************************

  # Check input arguments ----

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- LF_n_cores <- RC_n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)
  LF_n_cores <- .validate_n_cores(LF_n_cores)
  RC_n_cores <- .validate_n_cores(RC_n_cores)

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c(
      "hab_abb", "env_file", "model_dir", "use_trees", "strategy"))

  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("use_TF", "clamp_pred", "pred_new_sites"))

  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c(
      "n_cores", "GPP_dist", "MCMC_n_samples", "MCMC_thin", "n_grid",
      "LF_n_cores", "RC_n_cores"))
  rm(AllArgs, envir = environment())

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory is invalid or does not exist.", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!(use_trees %in% c("Tree", "NoTree"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for Tree argument. Valid values ",
        "are: 'Tree' or 'NoTree'"),
      use_trees = use_trees, include_backtrace = TRUE)
  }

  if (clamp_pred && is.null(fix_efforts)) {
    ecokit::stop_ctx(
      "`fix_efforts` can not be NULL when Clamping is implemented",
      clamp_pred = clamp_pred, fix_efforts = fix_efforts,
      include_backtrace = TRUE)
  }

  ValidModels <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
    "MRI-ESM2-0", "UKESM1-0-LL")

  if (!all(CC_models %in% ValidModels)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate models. Valid values are:\n >> ",
        toString(ValidModels)),
      CC_models = CC_models, include_backtrace = TRUE)
  }

  if (!all(CC_scenario %in% c("ssp126", "ssp370", "ssp585"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate scenarios. Valid values are:\n >> ",
        toString(c("ssp126", "ssp370", "ssp585"))),
      CC_scenario = CC_scenario, include_backtrace = TRUE)
  }

  # ****************************************************************

  ecokit::record_arguments(
    out_path = fs::path(model_dir, "Args_Mod_Postprocess_2_CPU.RData"))

  # ****************************************************************

  loaded_packages <- paste(
    sort(ecokit::loaded_packages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "\n  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", ecokit::os(),
      "\n  >>> Model root: ", model_dir,
      "\n  >>> n_cores: ", n_cores,
      "\n  >>> RC_n_cores: ", RC_n_cores,
      "\n  >>> LF_n_cores: ", LF_n_cores,
      "\n  >>> env_file: ", env_file,
      "\n  >>> hab_abb: ", hab_abb,
      "\n  >>> use_TF: ", use_TF,
      "\n  >>> Python environment: ", TF_environ,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"), # nolint: undesirable_function_linter
      "\n  >>> Loaded packages: \n", loaded_packages, "\n"))

  # ****************************************************************

  temp_dir <- fs::path(model_dir, "TEMP_Pred")

  ModelData <- list.files(
    path = model_dir, full.names = TRUE, pattern = "^ModDT_.*subset.RData")

  if (length(ModelData) != 1) {
    ecokit::stop_ctx(
      "Expected one model data file", length_model_data = length(ModelData),
      include_backtrace = TRUE)
  }

  # ****************************************************************

  # Path of selected model -----
  ecokit::info_chunk(
    "Path of selected model", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  path_model <- fs::dir_ls(
    fs::path(model_dir, "Model_Fitted"),
    regexp = paste0(
      ".+", use_trees, "_samp", MCMC_n_samples, "_th", MCMC_thin, "_Model.qs2"))
  if (length(path_model) > 1) {
    path_model <- stringr::str_subset(path_model, paste0("^GPP", GPP_dist))
  }
  if (length(path_model) == 0) {
    ecokit::stop_ctx("`path_model` is of length 0")
  }

  path_coda <- fs::dir_ls(
    fs::path(model_dir, "Model_Coda"),
    regexp = paste0(
      ".+", use_trees, "_samp", MCMC_n_samples, "_th", MCMC_thin, "_Coda.qs2"))
  if (length(path_coda) > 1) {
    path_coda <- stringr::str_subset(
      path_coda, paste0("^GPP", GPP_dist))
  }
  if (length(path_coda) == 0) {
    ecokit::stop_ctx("`path_coda` is of length 0")
  }

  cat(
    paste0("path_model:\n\t", path_model, "\nPath_Coda:\n\t", path_coda, "\n"))

  if (!all(fs::file_exists(c(path_model, path_coda)))) {
    ecokit::stop_ctx(
      "Selected model files not found",
      path_model = path_model, path_coda = path_coda, include_backtrace = TRUE)
  }

  # ****************************************************************

  # Prepare response curve data -----

  if (RC_prepare) {
    ecokit::info_chunk(
      "Prepare response curve data", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE,
      cat_bold = TRUE, cat_timestamp = FALSE)

    IASDT.R::resp_curv_prepare_data(
      path_model = path_model, n_grid = n_grid, n_cores = RC_n_cores,
      strategy = strategy, use_TF = use_TF, TF_environ = TF_environ,
      TF_use_single = TF_use_single, LF_n_cores = LF_n_cores,
      LF_temp_cleanup = LF_temp_cleanup, LF_check = LF_check,
      temp_dir = temp_dir, temp_cleanup = temp_cleanup, verbose = TRUE,
      LF_commands_only = FALSE, return_data = FALSE,
      probabilities = c(0.025, 0.5, 0.975))

    invisible(gc())
  }

  # ****************************************************************

  # Plotting response curves ----

  if (RC_plot) {

    ## Species richness -----
    ecokit::info_chunk(
      "Plotting response curves - species richness", level = 1L,
      line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::resp_curv_plot_SR(
      model_dir = model_dir, verbose = TRUE, n_cores = RC_n_cores,
      strategy = strategy)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Individual species -----
    ecokit::info_chunk(
      "Plotting response curves - species", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::resp_curv_plot_species(
      model_dir = model_dir, n_cores = RC_n_cores, env_file = env_file)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # All species together -----
    ecokit::info_chunk(
      "Plotting response curves - all species together", level = 1L,
      line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::resp_curv_plot_species_all(
      model_dir = model_dir, n_cores = RC_n_cores)

    invisible(gc())
  }

  # ****************************************************************

  # Predicting habitat suitability across different climate options -------

  if (predict_suitability) {

    ecokit::info_chunk(
      "Predicting habitat suitability across different climate options",
      line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)

    IASDT.R::predict_maps(
      path_model = path_model, hab_abb = hab_abb, env_file = env_file,
      n_cores = n_cores, strategy = strategy, clamp_pred = clamp_pred,
      fix_efforts = fix_efforts, fix_rivers = fix_rivers,
      pred_new_sites = pred_new_sites, use_TF = use_TF,
      TF_environ = TF_environ, TF_use_single = TF_use_single,
      LF_n_cores = LF_n_cores, LF_check = LF_check,
      LF_temp_cleanup = LF_temp_cleanup, LF_only = FALSE,
      LF_commands_only = FALSE, temp_dir = temp_dir,
      temp_cleanup = temp_cleanup, tar_predictions = tar_predictions,
      CC_models = CC_models, CC_scenario = CC_scenario)

    invisible(gc())
  }

  # ****************************************************************

  # Compute variance partitioning ------

  if (VP_prepare) {
    ecokit::info_chunk(
      "Compute variance partitioning", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::variance_partitioning_compute(
      path_model = path_model, n_cores = n_cores, use_TF = use_TF,
      TF_environ = TF_environ, TF_use_single = TF_use_single,
      temp_cleanup = temp_cleanup, chunk_size = 50L, verbose = TRUE,
      VP_file = "VarPar", VP_commands_only = FALSE)
  }

  # ****************************************************************

  # Plot Variance partitioning ------

  if (VP_plot) {
    ecokit::info_chunk(
      "Plot Variance partitioning", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::variance_partitioning_plot(
      path_model = path_model, env_file = env_file, VP_file = "VarPar",
      use_TF = use_TF, TF_environ = TF_environ, n_cores = n_cores, width = 30,
      height = 15)
  }

  # ****************************************************************

  # Plot species & SR predictions as JPEG ------

  if (plot_predictions) {
    ecokit::info_chunk(
      "Plot species & SR predictions as JPEG", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::plot_prediction(
      model_dir = model_dir, env_file = env_file, n_cores = n_cores)
  }

  # ****************************************************************

  # Plot latent factors as JPEG ------

  if (plot_LF) {
    ecokit::info_chunk(
      "Plot latent factors as JPEG", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::plot_latent_factor(path_model = path_model, env_file = env_file)
  }

  # ****************************************************************

  # Plot explanatory Power ------

  if (plot_internal_evaluation) {
    ecokit::info_chunk(
      "Plot explanatory Power", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::plot_evaluation(model_dir = model_dir)
  }

  # ****************************************************************

  ecokit::cat_diff(init_time = .start_time, prefix = "\nPost-processing took ")

  return(invisible(NULL))
}
