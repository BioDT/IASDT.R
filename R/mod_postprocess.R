## |------------------------------------------------------------------------| #
# Model pipeline for post-processing fitted Hmsc models
## |------------------------------------------------------------------------| #

#' Model pipeline for post-processing fitted Hmsc models
#'
#' These functions post-process fitted Hmsc models on both CPU and GPU. The main
#' functions in the pipeline includes `mod_postprocess_1_cpu`, `mod_prepare_tf`,
#' and `mod_postprocess_2_cpu` for full models without cross-validation, as well
#' as `mod_postprocess_cv_1_cpu` and `mod_postprocess_cv_2_cpu` for
#' cross-validated models. See details for more information.
#' @param model_dir Character. Path to the root directory of the fitted model.
#' @param job_runtime Character. Maximum allowed runtime for jobs for refitting
#'   the models (if needed) and cross validating models. Defaults to "01:00:00"
#'   for one hour. If not provided, the function throws an error.
#' @param gpp_dist Integer. Distance in *kilometres* between knots for the
#'   selected model.
#' @param use_trees Character. Whether a phylogenetic tree was used in the
#'   selected model. Accepts "tree" (default) or "no_tree".
#' @param mcmc_thin,mcmc_n_samples Integer. Thinning value and the number of
#'   MCMC samples of the selected model.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param future_max_size	Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details.
#' @param n_cores,n_cores_pred,n_cores_lf,n_cores_vp,n_cores_rc Integer. Number
#'   of cores to use for parallel processing. They are used for different
#'   processing steps: `n_cores` for merging chains and plotting convergence
#'   convergence diagnostics; `n_cores_pred` for predicting species' habitat
#'   suitability; `n_cores_lf` for predicting latent factors; `n_cores_vp` for
#'   processing variance partitioning; and `n_cores_rc` for response curve
#'   prediction. All default to `8L`. If `strategy = "sequential"`, all of these
#'   arguments are set to `1L`.
#' @param cv_name `NULL` or character vector. Column name(s) in the model input
#'   data to be used to cross-validate the models (see [mod_prepare_data] and
#'   [mod_cv_prepare]). If `cv_name = NULL`, no cross-validation data
#'   preparation is done. See [mod_cv_fit] for valid options.
#' @param n_batch_files Integer. Number of output batch files to create. Must be
#'   less than or equal to the maximum job limit of the HPC environment.
#' @param working_directory Character. Optionally sets the working directory in
#'   batch scripts to this path. If `NULL`, the directory remains unchanged.
#' @param model_prefix Character. Prefix for the model name. A directory named
#'   `model_prefix_TF` is created in the `model_dir` to store the `TensorFlow`
#'   running commands. Defaults to `NULL`. This can not be `NULL`.
#' @param partition_name Character. Name of the partition to submit the SLURM
#'   jobs to. Default is `small-g`.
#' @param lf_runtime,vp_runtime Character. Time limit for latent factor
#'   prediction and variance partitioning processing jobs, respectively.
#'   Defaults are `01:00:00` and `02:00:00` respectively.
#' @param temp_cleanup,lf_temp_cleanup Logical. Whether to delete temporary
#'   files after finishing predicting latent factor or species distribution.
#'   Default: `TRUE`.
#' @param width_omega,height_omega,width_beta,height_beta Integer. The width and
#'   height of the generated heatmaps of the Omega and Beta parameters in
#'   centimetres.
#' @param spatial_model,plot_lf Logical. Whether the model is spatial (`TRUE`)
#'   or not (`FALSE`) and whether to plot latent factors for spatial models as
#'   JPEG files (using [plot_latent_factor]). Defaults to `TRUE`.
#' @param rc_prepare,rc_plot Logical. Whether to prepare the data for response
#'   curve prediction (using [rc_prepare_data]) and plot the response curves as
#'   JPEG files. (using [rc_plot_sr], [rc_plot_species], and
#'   [rc_plot_species_all]). Defaults to `TRUE`.
#' @param vp_prepare,vp_plot Logical. Whether to prepare the data for variance
#'   partitioning (using [variance_partitioning_compute]) and plot its results
#'   (using [variance_partitioning_plot]). Defaults to `TRUE`.
#' @param predict_suitability,tar_predictions,plot_predictions Logical. Whether
#'   to predict habitat suitability across different climate options (using
#'   [predict_maps]), compress the resulted files into a single `*.tar` file
#'   (without compression), or to plot species and species richness predictions
#'   as JPEG files (using [plot_prediction]). Defaults to `TRUE`.
#' @param plot_internal_evaluation Logical. Whether to compute and visualise
#'   model internal evaluation (explanatory power) using [plot_evaluation].
#'   Defaults to `TRUE`.
#' @param process_vp,process_lf Logical. Whether to prepares batch scripts for
#'   variance partitioning computations and latent factor predictions on GPUs.
#'   Defaults to `TRUE`.
#' @param is_cv_model Logical. Whether the model is a cross-validated model
#'   (`TRUE`) or fitted with the full dataset (`FALSE`; default). If `TRUE`, the
#'   explanatory and predictive power of the model will be computed.
#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 1
#' @importFrom foreach %dopar%
#' @inheritParams predict_maps
#' @inheritParams mod_cv_fit
#' @inheritParams mod_merge_chains
#' @inheritParams mod_inputs
#' @inheritParams response_curves
#' @inheritParams coda_to_tibble
#' @inheritParams mod_merge_chains
#' @author Ahmed El-Gabbas
#' @export
#' @details
#'
#' **mod_postprocess_1_cpu**
#'
#' This function performs the initial post-processing step for habitat-specific
#' fitted models, automating the following tasks:
#'
#' - check unsuccessful models: [mod_slurm_refit]
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
#' [mod_cv_fit]
#' - Prepare scripts for GPU processing, including:
#'    - predicting latent factors of the response curves:
#' [rc_prepare_data]
#'    - predicting latent factors for new sampling units: [predict_maps]
#'    - computing variance partitioning: [variance_partitioning_compute]
#'
#' <br/>
#'
#' **mod_prepare_tf**
#'
#' After running `mod_postprocess_1_cpu` for all habitat types, this function
#' prepares batch scripts for GPU computations of all habitat types:
#'    - for <u>variance partitioning</u>, the function matches all files with
#' the pattern ` "vp_.+command.txt"` (created by [variance_partitioning_compute]
#' and merges their contents into a single file
#' (`model_prefix_TF/vp_commands.txt`). Then, it prepares a SLURM script for
#' variance partitioning computations (`model_prefix_TF/vp_slurm.slurm`).
#'    - for <u>latent factor predictions</u>, the function matches all files
#' with the pattern `"^lf_new_sites_commands_.+.txt|^lf_rc_commands_.+txt"` and
#' split their contents into multiple scripts at the `model_prefix_TF` directory
#' for processing as a batch job. The function prepares a SLURM script for
#' latent factor predictions (`lf_slurm.slurm`).
#'
#' This function is tailored for the LUMI HPC environment and assumes that the
#' `tensorflow` module is installed and correctly configured with all required
#' Python packages. On other HPC systems, users may need to modify the function
#' to load a Python virtual environment or install the required dependencies for
#' `TensorFlow` and related packages.
#'
#' <br/><br/>
#'
#' **mod_postprocess_2_cpu**
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
#' - completing `mod_postprocess_1_cpu` and `mod_prepare_tf` on CPU,
#' - running `vp_slurm.slurm` and `lf_slurm.slurm` on GPU to process response
#' curves and latent factor predictions (both scripts are generated by
#' `mod_prepare_tf`).
#' - submitting SLURM jobs for cross-validated model fitting.
#'
#' <br/>
#'
#' **mod_postprocess_cv_1_cpu**
#'
#' This function is similar to `mod_postprocess_1_cpu`, but it is specifically
#' designed for cross-validated models. It automates merging fitted
#' cross-validated model chains into `Hmsc` model objects and prepare scripts
#' for latent factor prediction on `TensorFlow` using [predict_maps_cv].
#'
#' <br/><br/>
#'
#' **mod_postprocess_cv_2_cpu**
#'
#' The function 1) processes `*.feather` files resulted from Latent Factor
#' predictions (using `TensorFlow`) and saves LF predication to disk; 2)
#' predicts species-specific mean habitat suitability at testing
#' cross-validation folds and calculates testing evaluation metrics; 3)
#' generates plots of the evaluation metrics.
#'

## |------------------------------------------------------------------------| #
# mod_postprocess_1_cpu ----
## |------------------------------------------------------------------------| #

mod_postprocess_1_cpu <- function(
    model_dir = NULL, hab_abb = NULL, strategy = "multisession",
    future_max_size = 1500L, n_cores = 8L, n_cores_pred = n_cores,
    n_cores_lf = n_cores, n_cores_vp = n_cores, env_file = ".env",
    path_hmsc = NULL, memory_per_cpu = "64G", job_runtime = "01:00:00",
    from_json = FALSE, gpp_dist = NULL, use_trees = "tree",
    mcmc_n_samples = 1000L, mcmc_thin = NULL, n_omega = 1000L,
    cv_name = c("cv_dist", "cv_large"), n_grid = 50L, use_tf = TRUE,
    tf_use_single = FALSE,  lf_temp_cleanup = TRUE, lf_check = FALSE,
    temp_cleanup = TRUE, tf_environ = NULL, pred_new_sites = TRUE,
    width_omega = 26, height_omega = 22.5, width_beta = 25, height_beta = 35,
    spatial_model = TRUE, tar_predictions = TRUE, plot_predictions = TRUE,
    is_cv_model = FALSE, clamp_pred = TRUE, fix_efforts = "q90",
    fix_rivers = "q90",
    climate_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    climate_scenario = c("ssp126", "ssp370", "ssp585")) {

  .start_time <- lubridate::now(tzone = "CET")

  # ****************************************************************

  # Check input arguments ----
  ecokit::check_args(
    args_to_check = c("model_dir", "use_trees", "path_hmsc"),
    args_type = "character")

  ecokit::check_args(
    args_to_check = c(
      "from_json", "use_tf", "tf_use_single", "lf_temp_cleanup", "lf_check",
      "temp_cleanup", "pred_new_sites", "spatial_model", "tar_predictions",
      "plot_predictions", "clamp_pred", "is_cv_model"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c(
      "future_max_size", "n_omega", "mcmc_n_samples", "mcmc_thin", "n_grid",
      "height_omega", "width_omega", "width_beta", "height_beta"),
    args_type = "numeric")

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") {
    n_cores <- n_cores_lf <- n_cores_vp <- n_cores_pred <- 1L
  }
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)
  n_cores_vp <- .validate_n_cores(n_cores_vp)
  n_cores_pred <- .validate_n_cores(n_cores_pred)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory was not found", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!(use_trees %in% c("tree", "no_tree"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for tree argument. Valid values ",
        "are: 'tree' or 'no_tree'"),
      use_trees = use_trees, include_backtrace = TRUE)
  }

  if (!is.null(cv_name) &&
      !all(cv_name %in% c("cv_dist", "cv_large", "cv_sac"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for cv_name argument. Valid values ",
        "are: 'cv_dist', 'cv_large', or `cv_sac`"),
      cv_name = cv_name, include_backtrace = TRUE)
  }

  # Validate memory_per_cpu
  memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
  job_runtime <- .validate_slurm_runtime(job_runtime)

  # ****************************************************************

  species_name <- non_focal_variables <- NULL

  ecokit::record_arguments(
    out_path = fs::path(model_dir, "args_mod_postprocess_1_cpu.RData"))

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
      "\n  >>> path_hmsc: ", path_hmsc,
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

  IASDT.R::mod_slurm_refit(
    model_dir = model_dir,
    job_name = stringr::str_remove(basename(model_dir), "mod_"),
    memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
    env_file = env_file, path_hmsc = path_hmsc)

  invisible(gc())

  # ****************************************************************

  ecokit::info_chunk(
    "Merge chains and saving RData or qs2 files", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_merge_chains(
    model_dir = model_dir, n_cores = n_cores, strategy = strategy,
    future_max_size = future_max_size, from_json = from_json)

  invisible(gc())

  # ****************************************************************

  ecokit::info_chunk(
    "Convergence of all model variants", level = 1L, line_char = "+",
    line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::convergence_plot_all(
    model_dir = model_dir, n_omega = n_omega, margin_type = "histogram",
    spatial_model = spatial_model, n_cores = n_cores,
    strategy = strategy, future_max_size = future_max_size)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Path of selected model -----
  ecokit::info_chunk(
    "Path of selected model", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  path_model <- fs::dir_ls(
    fs::path(model_dir, "model_fitted"),
    regexp = paste0(
      ".+", use_trees, "_samp", mcmc_n_samples, "_th", mcmc_thin, "_model.qs2"))
  if (length(path_model) > 1) {
    path_model <- stringr::str_subset(path_model, paste0("^gpp", gpp_dist))
  }
  if (length(path_model) == 0) {
    ecokit::stop_ctx("`path_model` is of length 0")
  }

  path_coda <- fs::dir_ls(
    fs::path(model_dir, "model_coda"),
    regexp = paste0(
      ".+", use_trees, "_samp", mcmc_n_samples, "_th", mcmc_thin, "_coda.qs2"))
  if (length(path_coda) > 1) {
    path_coda <- stringr::str_subset(
      path_coda, paste0("^gpp", gpp_dist))
  }
  if (length(path_coda) == 0) {
    ecokit::stop_ctx("`path_coda` is of length 0")
  }

  cat(
    paste0("path_model:\n\t", path_model, "\npath_coda:\n\t", path_coda, "\n"))

  if (!all(fs::file_exists(c(path_model, path_coda)))) {
    ecokit::stop_ctx(
      "Selected model files not found",
      path_model = path_model, path_coda = path_coda, include_backtrace = TRUE)
  }

  temp_dir_pred <- fs::path(model_dir, "temp_pred")
  temp_dir_vp <- fs::path(model_dir, "temp_vp")
  fs::dir_create(c(temp_dir_pred, temp_dir_vp))

  # ****************************************************************

  # gelman_plot -----
  ecokit::info_chunk(
    "gelman_plot", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::plot_gelman(
    path_coda = path_coda, alpha = spatial_model,
    beta = TRUE, omega = TRUE, rho = TRUE,
    n_omega = n_omega, env_file = env_file)

  invisible(gc())

  # ****************************************************************

  # Convergence plots ----
  ecokit::info_chunk(
    "Convergence plots", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::convergence_plot(
    path_coda = path_coda, env_file = env_file, n_omega = n_omega,
    n_cores = n_cores, strategy = strategy,
    n_rc = list(alpha = c(2L, 2L), omega = c(2L, 2L), beta = c(3L, 3L)),
    margin_type = "histogram", spatial_model = spatial_model,
    future_max_size = future_max_size)

  invisible(gc())

  # ****************************************************************

  # Model summary ------
  ecokit::info_chunk(
    "Model summary", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  IASDT.R::mod_summary(
    path_coda = path_coda, env_file = env_file, spatial_model = spatial_model)

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

  if (!is.null(cv_name)) {
    ecokit::info_chunk(
      "Prepare input data for cross-validation", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::mod_cv_fit(
      path_model = path_model, cv_name = cv_name, env_file = env_file,
      job_name = paste0("cv_", hab_abb), memory_per_cpu = memory_per_cpu,
      job_runtime = job_runtime, path_hmsc = path_hmsc)

    invisible(gc())
  }

  # ****************************************************************

  ## Response curves -------

  if (spatial_model) {

    ecokit::info_chunk(
      "Prepare scripts for predicting latent factors of the response curves",
      line_char = "+", line_char_rep = 90L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)

    IASDT.R::rc_prepare_data(
      path_model = path_model, n_grid = n_grid, n_cores = n_cores,
      strategy = strategy, future_max_size = future_max_size, use_tf = use_tf,
      tf_environ = tf_environ, tf_use_single = tf_use_single,
      n_cores_lf = n_cores_lf, lf_temp_cleanup = lf_temp_cleanup,
      lf_check = lf_check, temp_dir = temp_dir_pred,
      temp_cleanup = temp_cleanup, verbose = TRUE, lf_commands_only = TRUE)

    invisible(gc())

  } else {

    ecokit::info_chunk(
      "Prepare response curve data", line_char = "+", line_char_rep = 90L,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, level = 1L)

    rc_dir <- fs::path(
      dirname(dirname(path_model)),
      "model_postprocessing", "response_curves_data")
    fs::dir_create(rc_dir)
    path_rc_data <- fs::path(rc_dir, "response_curve_data.qs2")
    path_observed_data <- fs::path(rc_dir, "observed_data.qs2")

    # Observed data
    if (!ecokit::check_data(path_observed_data, warning = FALSE)) {
      model_obj <- ecokit::load_as(path_model)
      model_vars <- names(model_obj$XData)
      observed_data <- dplyr::bind_cols(
        x_value = model_obj$XData, model_obj$Y,
        species_richness = rowSums(model_obj$Y)) %>%
        tidyr::pivot_longer(
          cols = -tidyselect::all_of(model_vars), names_to = "species_name",
          values_to = "observed_value") %>%
        dplyr::select("species_name", tidyselect::everything()) %>%
        dplyr::arrange(species_name)

      ecokit::save_as(object = observed_data, out_path = path_observed_data)
      rm(model_obj, observed_data, envir = environment())
      invisible(gc())
    }

    if (!ecokit::check_data(path_rc_data, warning = FALSE)) {
      model_obj <- ecokit::load_as(path_model)
      model_vars <- names(model_obj$XData)
      rm(model_obj, envir = environment())

      response_curve_data <- tidyr::expand_grid(
        variable = model_vars, non_focal_variables = c(1, 2))

      # Prepare working in parallel
      if (n_cores == 1) {
        future::plan("sequential", gc = TRUE)
      } else {
        ecokit::set_parallel(
          n_cores = min(n_cores, nrow(response_curve_data)), level = 2L,
          strategy = strategy, future_max_size = future_max_size)
        withr::defer(future::plan("sequential", gc = TRUE))
      }

      response_curve_data2 <- future.apply::future_lapply(
        X = seq_len(nrow(response_curve_data)),
        FUN = function(x) {

          model_obj <- ecokit::load_as(path_model)
          variable <- response_curve_data$variable[x]
          nfv <- response_curve_data$non_focal_variables[x]

          Gradient <- Hmsc::constructGradient(
            hM = model_obj, focalVariable = variable,
            non.focalVariables = nfv, ngrid  = 50)
          pred <- stats::predict(
            object = model_obj, XData = Gradient$XDataNew,
            studyDesign = Gradient$studyDesignNew,
            ranLevels = Gradient$rLNew, expected = TRUE) %>%
            abind::abind(along = 3)
          x_vals <- Gradient$XDataNew[, variable]
          rm(model_obj, Gradient, envir = environment())
          invisible(gc())

          pred_mean <- apply(pred, c(1, 2), mean)
          rich_array <- apply(pred, 3, rowSums)
          rich_mean <- rowMeans(rich_array)
          pred_mean <- cbind(pred_mean, species_richness = rich_mean)

          pred_sd <- apply(pred, c(1, 2), stats::sd)
          rich_sd <- apply(rich_array, 1, stats::sd)
          pred_sd <- cbind(pred_sd, species_richness = rich_sd)

          pred_mean_plus  <- pred_mean + pred_sd
          pred_mean_minus <- pred_mean - pred_sd
          rm(rich_array, rich_mean, pred_sd, rich_sd, pred)
          invisible(gc())

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
          "Hmsc", "dplyr", "magrittr", "abind", "tibble", "qs2",
          "tidyr", "ecokit", "tidyselect"))

      if (n_cores > 1) {
        ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
        future::plan("sequential", gc = TRUE)
      }

      response_curve_data <- dplyr::bind_rows(response_curve_data2)
      ecokit::save_as(object = response_curve_data, out_path = path_rc_data)
      rm(response_curve_data, envir = environment())
      invisible(gc())
    }

  }

  # ****************************************************************

  ## Prepare scripts for computing variance partitioning -------
  ecokit::info_chunk(
    "Prepare scripts for computing variance partitioning",
    line_char = "+", line_char_rep = 90L, cat_red = TRUE, cat_bold = TRUE,
    cat_timestamp = FALSE, level = 1L)
  IASDT.R::variance_partitioning_compute(
    path_model = path_model, n_cores = n_cores_vp, use_tf = use_tf,
    tf_environ = tf_environ, tf_use_single = tf_use_single,
    temp_cleanup = temp_cleanup, chunk_size = 50L, verbose = TRUE,
    vp_file = "varpar", vp_commands_only = TRUE, temp_dir = temp_dir_vp)

  # ****************************************************************

  # Predicting latent factors or habitat suitability -------

  if (spatial_model) {
    ecokit::info_chunk(
      "Prepare scripts for predicting latent factors for new sampling units",
      line_char = "+", line_char_rep = 90L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)
  } else {
    ecokit::info_chunk(
      "Predicting habitat suitability across different climate options",
      line_char = "+", line_char_rep = 90L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE, level = 1L)
  }

  IASDT.R::predict_maps(
    path_model = path_model, hab_abb = hab_abb, env_file = env_file,
    n_cores = n_cores, n_cores_pred = n_cores_pred, strategy = strategy,
    future_max_size = future_max_size,
    # Do not clamp predictions for predicting latent factors
    clamp_pred = dplyr::if_else(spatial_model, FALSE, clamp_pred),
    fix_efforts = fix_efforts, fix_rivers = fix_rivers,
    pred_new_sites = pred_new_sites, use_tf = use_tf, tf_environ = tf_environ,
    temp_dir = temp_dir_pred, temp_cleanup = temp_cleanup,
    tf_use_single = tf_use_single, n_cores_lf = n_cores_lf,
    lf_check = lf_check, lf_temp_cleanup = lf_temp_cleanup,
    lf_only = spatial_model, lf_commands_only = spatial_model,
    climate_models = climate_models, climate_scenario = climate_scenario,
    spatial_model = spatial_model, tar_predictions = tar_predictions)

  if (isFALSE(spatial_model)) {

    # Plotting prediction maps as JPEG files -------
    if (plot_predictions) {
      ecokit::info_chunk(
        "Plot species & SR predictions as JPEG", level = 1L, line_char = "+",
        line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)
      IASDT.R::plot_prediction(
        model_dir = model_dir, env_file = env_file, n_cores = n_cores,
        is_cv_model = is_cv_model)
    }

    # ************************************

    # Evaluate cross-validated models -----
    if (is_cv_model) {
      ecokit::info_chunk(
        "Evaluating explanatory and predictive power of the model",
        level = 1L, line_char = "+", line_char_rep = 60L, cat_red = TRUE,
        cat_bold = TRUE, cat_timestamp = FALSE)

      IASDT.R::mod_cv_evaluate(model_dir = model_dir, cv_type = "cv_dist")
    }

    # ************************************

    # Tar prediction directories and delete TIFF files -----
    ecokit::cat_time(
      "Create tar file for each prediction option and delete TIFF files",
      level = 1L)

    path_prediction <- fs::path(model_dir, "model_prediction")
    path_prediction_clamp <- fs::path(path_prediction, "clamp")
    path_prediction_no_clamp <- fs::path(path_prediction, "no_clamp")
    l_dirs_preds <- c(
      list.dirs(path_prediction_clamp, recursive = FALSE),
      list.dirs(path_prediction_no_clamp, recursive = FALSE))

    if (length(l_dirs_preds) > 0L) {
      ecokit::set_parallel(
        n_cores = n_cores, show_log = FALSE, strategy = strategy,
        future_max_size = future_max_size)
      tar_pred_tiffs0 <- future.apply::future_lapply(
        X = l_dirs_preds, FUN = tar_pred_tiffs, delete_tiff = TRUE,
        overwrite_tar = TRUE, future.scheduling = Inf, future.seed = TRUE,
        future.packages = c("stringr", "fs"),
        future.globals = c("l_dirs_preds", "tar_pred_tiffs"))
      rm(tar_pred_tiffs0, l_dirs_preds, envir = environment())
      ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
    }

  }

  # ****************************************************************

  ecokit::cat_diff(
    init_time = .start_time, prefix = "\nPost-processing using CPU took ")

  return(invisible(NULL))
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_prepare_tf ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_postprocessing
#' @name mod_postprocessingS
#' @order 2
#' @author Ahmed El-Gabbas

mod_prepare_tf <- function(
    process_vp = TRUE, process_lf = TRUE,
    n_batch_files = 210L, env_file = ".env", working_directory = NULL,
    partition_name = "small-g", lf_runtime = "01:00:00", model_prefix = NULL,
    vp_runtime = "02:00:00") {

  # ****************************************************************

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  project_id <- path_model <- NULL

  # ****************************************************************

  # Check input arguments ----
  ecokit::check_args(
    args_to_check = c(
      "lf_runtime", "vp_runtime", "partition_name", "model_prefix"),
    args_type = "character")
  ecokit::check_args(args_to_check = "n_batch_files", args_type = "numeric")
  ecokit::check_args(
    args_to_check = c("process_vp", "process_lf"), args_type = "logical")

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

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "project_id", "DP_R_lumi_gpu", FALSE, FALSE,
    "path_model", "DP_R_model_root_path", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # ****************************************************************

  # Path to store TF commands
  path_tf <- fs::path(path_model, paste0(model_prefix, "_TF"))
  # Path to store log files
  path_log <- ecokit::normalize_path(fs::path(path_tf, "log"))
  fs::dir_create(c(path_tf, path_log))

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  # VARIANCE PARTITIONING ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  if (process_vp) {

    # Prepare post-processing data for calculating variance partitioning
    ecokit::cat_time(
      "Prepare post-processing data for calculating variance partitioning")

    # Create paths for VP SLURM script and commands
    path_vp_slurm <- fs::path(path_tf, "vp_slurm.slurm")
    path_vp_commands <- fs::path(path_tf, "vp_commands.txt")

    # Merge and organise `TensorFlow` commands for computing variance
    # partitioning
    ecokit::cat_time(
      "Merge and organise `TensorFlow` for computing variance partitioning",
      level = 1L, cat_timestamp = FALSE)

    # Find list of files matching the pattern

    # fs::dir_ls gives `Segmentation fault (core dumped)` error when too many
    # files exist in the path_model directory
    # vp_in_files <- fs::dir_ls(
    #   path = path_model, recurse = TRUE, type = "file",
    #   regexp = paste0(model_prefix, ".+/temp_vp/vp_.+command.txt")) %>%
    #   purrr::map(readr::read_lines, progress = FALSE) %>%
    #   unlist() %>%
    #   gtools::mixedsort() %>%
    #   unique() %>%
    #   unname()

    vp_in_files <- list.dirs(path = path_model, recursive = FALSE) %>%
      stringr::str_subset(pattern = model_prefix) %>%
      purrr::map(
        .f = list.files, recursive = TRUE, full.names = TRUE,
        pattern = "vp_.+command.txt") %>%
      unlist() %>%
      gtools::mixedsort() %>%
      unique() %>%
      unname()

    n_vp_in_files <- length(vp_in_files)

    # Save all VP commands to single file for batch processing
    readr::write_lines(x = vp_in_files, file = path_vp_commands, append = FALSE)

    # ****************************************************************

    # Prepare batch files ----
    ecokit::cat_time("Prepare batch files", level = 1L, cat_timestamp = FALSE)

    vp_commands <- c(
      "#!/bin/bash",
      "#SBATCH --job-name=vp_tf",
      "#SBATCH --ntasks=1",
      paste0("#SBATCH --account=", project_id),
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --gpus-per-node=1",
      paste0("#SBATCH --time=", vp_runtime),
      paste0("#SBATCH --partition=", partition_name),
      "#SBATCH --mem=64G",
      paste0("#SBATCH --output=", fs::path(path_log, "%x-%A-%a.out")),
      paste0("#SBATCH --error=", fs::path(path_log, "%x-%A-%a.out")),
      paste0("#SBATCH --array=1-", length(vp_in_files)),
      "\n# file containing commands to be executed",
      paste0("file=", path_vp_commands),
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
      "head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | bash",
      "\necho End of program at `date`\n",
      paste0("# ", strrep("-", 50)),
      paste0(
        "# This script was created on: ",
        format(lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"),
        " CET"),
      paste0("# ", strrep("-", 50)))

    ecokit::cat_time(
      paste0("Writing SLURM script to: `", path_vp_slurm, "`"),
      level = 1L, cat_timestamp = FALSE)

    readr::write_lines(x = vp_commands, file = path_vp_slurm, append = FALSE)
    # Make the file executable
    Sys.chmod(path_vp_slurm, mode = "755")

    ecokit::cat_time(
      paste0(
        "To submit variance partitioning SLURM script:\n\tsbatch ",
        path_vp_slurm),
      cat_timestamp = FALSE)

  } else {

    n_vp_in_files <- 0L

  }

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  ## LF PREDICTIONS ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  if (process_lf) {

    # Prepare post-processing data for LF predictions
    ecokit::cat_time("\nPrepare post-processing data for LF predictions")

    # Ensure that the total number of simultaneous jobs (LF + VP) =
    # n_batch_files; so both can be run on the same time.
    n_batch_files <- n_batch_files - n_vp_in_files

    ## Merge and organise `TensorFlow` commands for LF predictions ----
    ecokit::cat_time(
      paste0(
        "Merge and organise `TensorFlow` commands for LF predictions ",
        "into a maximum of ", n_batch_files, " files"),
      level = 1L, cat_timestamp = FALSE)

    # Basic commands for `TensorFlow` setup
    basic_commands <- c(
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
      basic_commands <- c(
        basic_commands, "# Change to working directory",
        paste0("cd ", working_directory), "")
    }

    # Find list of files matching the pattern
    # Regex pattern to match input files
    lf_pattern <- "(lf_new_sites_commands_.+.txt|lf_rc_commands_.+txt)"

    # fs::dir_ls gives `Segmentation fault (core dumped)` error when too many
    # files exist in the path_model directory
    # lf_in_files <- fs::dir_ls(
    #   path = path_model, recurse = TRUE, type = "file",
    #   regexp = paste0(model_prefix, ".+/temp_pred/", lf_pattern))

    lf_in_files <- list.dirs(path = path_model, recursive = FALSE) %>%
      stringr::str_subset(pattern = model_prefix) %>%
      purrr::map(
        .f = list.files, recursive = TRUE, full.names = TRUE,
        pattern = lf_pattern) %>%
      unlist() %>%
      gtools::mixedsort() %>%
      unique() %>%
      unname()

    if (length(lf_in_files) == 0) {
      ecokit::stop_ctx(
        "No files found matching the pattern",
        lf_pattern = lf_pattern, path_model = path_model,
        include_backtrace = TRUE)
    }

    ecokit::cat_time(
      paste0(
        "Found ", length(lf_in_files), " files matching the pattern `",
        lf_pattern, "`"),
      level = 1L, cat_timestamp = FALSE)
    purrr::walk(
      lf_in_files, ecokit::cat_time, level = 2L, cat_timestamp = FALSE)

    # Read and merge commands from input files
    lf_commands <- purrr::map(
      .x = lf_in_files, .f = readr::read_lines, progress = FALSE) %>%
      unlist() %>%
      gtools::mixedsort() %>%
      unname()

    ecokit::cat_time(
      paste0(
        "Total number of commands to be executed: ", length(lf_commands)),
      level = 1L, cat_timestamp = FALSE)

    if (length(lf_commands) < n_batch_files) {
      ecokit::cat_time(
        paste0(
          "Fewer commands than the requested number of files. ",
          "Setting `n_batch_files=", n_batch_files, "`."),
        level = 2L, cat_timestamp = FALSE)
      n_batch_files <- length(lf_commands)
    }

    ecokit::cat_time(
      paste0("Splitting commands into ", n_batch_files, " files"),
      cat_timestamp = FALSE, level = 1L)
    lf_commands <- ecokit::split_vector(lf_commands, n_splits = n_batch_files)

    purrr::walk(
      .x = seq_len(length(lf_commands)),
      .f = ~ {

        file <- fs::path(
          path_tf,
          paste0(
            "tf_chunk_",
            stringr::str_pad(
              .x, width = nchar(n_batch_files), pad = "0"), ".txt"))

        readr::write_lines(x = basic_commands, file = file, append = FALSE)
        readr::write_lines(
          x = paste0(
            "# ", length(lf_commands[[.x]]), " commands to be executed:"),
          file = file, append = TRUE)
        readr::write_lines(x = lf_commands[[.x]], file = file, append = TRUE)
        readr::write_lines(
          x = c(
            paste0("\n#", strrep("_", 60)),
            paste0(
              "# This script was created on: ",
              format(lubridate::now(tzone = "CET"), "%Y-%m-%d %H:%M:%S")),
            paste0("#", strrep("_", 60))),
          file = file, append = TRUE)
        Sys.chmod(file, mode = "755")

        return(invisible(NULL))
      })

    # ****************************************************************

    ## Prepare LF batch files ----
    ecokit::cat_time("Prepare batch files", level = 1L, cat_timestamp = FALSE)

    lf_slurm_script <- c(
      "#!/bin/bash",
      "#SBATCH --job-name=PP_lf",
      "#SBATCH --ntasks=1",
      paste0("#SBATCH --account=", project_id),
      "#SBATCH --mem=64G",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --gpus-per-node=1",
      paste0("#SBATCH --time=", lf_runtime),
      paste0("#SBATCH --partition=", partition_name),
      paste0("#SBATCH --output=", fs::path(path_log, "%x-%A-%a.out")),
      paste0("#SBATCH --error=", fs::path(path_log, "%x-%A-%a.out")),
      paste0("#SBATCH --array=1-", n_batch_files),
      "",
      "# Define directories",
      paste0('output_directory="', path_tf, '"'),
      "",
      "# Find all the split files and sort them explicitly",
      paste0(
        'split_files=($(find "$output_directory" -type f ',
        '-name "tf_chunk_*.txt" | sort -V))'),
      "",
      "# Check if files were found",
      "if [ ${#split_files[@]} -eq 0 ]; then",
      '    echo "Error: No files matching tf_chunk_*.txt found in $output_directory"',
      "    exit 1",
      "fi",
      "",
      paste0("# Ensure no more than `, n_batch_files, ` files are processed"),
      paste0("MaxFiles=", n_batch_files),
      "if [ ${#split_files[@]} -gt $MaxFiles ]; then",
      '    split_files=("${split_files[@]:0:$MaxFiles}")',
      paste0(
        '    echo "More than $MaxFiles files found, ',
        'limiting to the first $MaxFiles files."'),
      "fi",
      "",
      "# Get the index of the current task based on SLURM_ARRAY_TASK_ID",
      "TaskIndex=$((SLURM_ARRAY_TASK_ID - 1))",
      "",
      "# Validate TaskIndex",
      "if [ $TaskIndex -ge ${#split_files[@]} ] || [ $TaskIndex -lt 0 ]; then",
      paste0(
        '    echo "Error: TaskIndex $TaskIndex is out of range. ',
        'Valid range: 0 to $((${#split_files[@]} - 1))"'),
      "    exit 1",
      "fi",
      "",
      "# Get the specific split file to process based on the job array task ID",
      'split_file="${split_files[$TaskIndex]}"',
      "",
      "# Verify the selected split file",
      'if [ -z "$split_file" ] || [ ! -f "$split_file" ]; then',
      '    echo "Error: file $split_file does not exist or is invalid."',
      "    exit 1",
      "fi",
      "",
      "# Processing file",
      'echo "Processing file: $split_file"',
      "",
      "# Run the selected split file",
      'bash "$split_file"',
      "\necho End of program at `date`\n",
      paste0("# ", strrep("-", 50)),
      paste0(
        "# This script was created on: ",
        format(
          lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"), " CET"),
      paste0("# ", strrep("-", 50)))

    path_lf_slurm <- fs::path(path_tf, "lf_slurm.slurm")
    ecokit::cat_time(
      paste0("Writing SLURM script to: `", path_lf_slurm, "`"),
      level = 1L, cat_timestamp = FALSE)
    # Write the content to a file
    readr::write_lines(lf_slurm_script, path_lf_slurm, append = FALSE)
    # Make the file executable
    Sys.chmod(path_lf_slurm, mode = "755")

    ecokit::cat_time(
      paste0(
        "\nTo submit LF prediction SLURM script:\n\tsbatch ", path_lf_slurm),
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
# mod_postprocess_2_cpu ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_postprocessing
#' @name mod_postprocessing
#' @order 3
#' @author Ahmed El-Gabbas

mod_postprocess_2_cpu <- function(
    model_dir = NULL, hab_abb = NULL, strategy = "multisession",
    future_max_size = 1500L, n_cores = 8L, n_cores_pred = n_cores,
    n_cores_lf = n_cores, n_cores_rc = n_cores, n_cores_vp = n_cores,
    env_file = ".env", gpp_dist = NULL, use_trees = "tree",
    mcmc_n_samples = 1000L, mcmc_thin = NULL, use_tf = TRUE, tf_environ = NULL,
    tf_use_single = FALSE, lf_check = FALSE, lf_temp_cleanup = TRUE,
    temp_cleanup = TRUE, n_grid = 50L,
    climate_models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    climate_scenario = c("ssp126", "ssp370", "ssp585"),
    clamp_pred = TRUE, fix_efforts = "q90", fix_rivers = "q90",
    pred_new_sites = TRUE, tar_predictions = TRUE,
    rc_prepare = TRUE, rc_plot = TRUE, vp_prepare = TRUE, vp_plot = TRUE,
    predict_suitability = TRUE, plot_predictions = TRUE, plot_lf = TRUE,
    plot_internal_evaluation = TRUE, spatial_model = TRUE,
    is_cv_model = FALSE) {

  .start_time <- lubridate::now(tzone = "CET")

  # ****************************************************************

  # Check input arguments ----

  ecokit::check_args(
    args_to_check = c("model_dir", "use_trees"), args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "use_tf", "clamp_pred", "temp_cleanup", "rc_prepare", "tar_predictions",
      "rc_plot", "vp_prepare", "vp_plot", "predict_suitability", "lf_check",
      "plot_predictions", "plot_lf", "plot_internal_evaluation",
      "spatial_model", "is_cv_model",  "lf_temp_cleanup", "pred_new_sites"),
    args_type = "logical")
  ecokit::check_args(
    args_to_check = c(
      "gpp_dist", "mcmc_n_samples", "mcmc_thin", "n_grid", "future_max_size"),
    args_type = "numeric")

  hab_abb <- .validate_hab_abb(as.character(hab_abb))

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") {
    n_cores <- n_cores_lf <- n_cores_rc <- n_cores_pred <- 1L
  }
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)
  n_cores_rc <- .validate_n_cores(n_cores_rc)
  n_cores_pred <- .validate_n_cores(n_cores_pred)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory is invalid or does not exist.", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (!(use_trees %in% c("tree", "no_tree"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for tree argument. Valid values ",
        "are: 'tree' or 'no_tree'"),
      use_trees = use_trees, include_backtrace = TRUE)
  }

  if (clamp_pred && is.null(fix_efforts)) {
    ecokit::stop_ctx(
      "`fix_efforts` can not be NULL when clamping is implemented",
      clamp_pred = clamp_pred, fix_efforts = fix_efforts,
      include_backtrace = TRUE)
  }

  valid_models <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
    "MRI-ESM2-0", "UKESM1-0-LL")

  if (!all(climate_models %in% valid_models)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate models. Valid values are:\n >> ",
        toString(valid_models)),
      climate_models = climate_models, include_backtrace = TRUE)
  }

  if (!all(climate_scenario %in% c("ssp126", "ssp370", "ssp585"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate scenarios. Valid values are:\n >> ",
        toString(c("ssp126", "ssp370", "ssp585"))),
      climate_scenario = climate_scenario, include_backtrace = TRUE)
  }

  # ****************************************************************

  ecokit::record_arguments(
    out_path = fs::path(model_dir, "args_mod_postprocess_2_cpu.RData"))

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
      "\n  >>> n_cores_rc: ", n_cores_rc,
      "\n  >>> n_cores_lf: ", n_cores_lf,
      "\n  >>> env_file: ", env_file,
      "\n  >>> hab_abb: ", hab_abb,
      "\n  >>> use_tf: ", use_tf,
      "\n  >>> Python environment: ", tf_environ,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"), # nolint: undesirable_function_linter
      "\n  >>> Loaded packages: \n", loaded_packages, "\n"))

  # ****************************************************************

  temp_dir_pred <- fs::path(model_dir, "temp_pred")
  temp_dir_vp <- fs::path(model_dir, "temp_vp")
  fs::dir_create(c(temp_dir_pred, temp_dir_vp))

  model_data <- fs::path(model_dir, "model_data_subset.RData")
  if (!ecokit::check_data(model_data)) {
    ecokit::stop_ctx(
      "Model data file not found",
      model_data = model_data, include_backtrace = TRUE)
  }
  model_data <- ecokit::load_as(model_data)

  # ****************************************************************

  # Path of selected model -----
  ecokit::info_chunk(
    "Path of selected model", level = 1L, line_char = "+", line_char_rep = 60L,
    cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

  path_model <- fs::dir_ls(
    fs::path(model_dir, "model_fitted"),
    regexp = paste0(
      ".+", use_trees, "_samp", mcmc_n_samples, "_th", mcmc_thin, "_model.qs2"))
  if (length(path_model) > 1) {
    path_model <- stringr::str_subset(path_model, paste0("^gpp", gpp_dist))
  }
  if (length(path_model) == 0) {
    ecokit::stop_ctx("`path_model` is of length 0")
  }

  path_coda <- fs::dir_ls(
    fs::path(model_dir, "model_coda"),
    regexp = paste0(
      ".+", use_trees, "_samp", mcmc_n_samples, "_th", mcmc_thin, "_coda.qs2"))
  if (length(path_coda) > 1) {
    path_coda <- stringr::str_subset(
      path_coda, paste0("^gpp", gpp_dist))
  }
  if (length(path_coda) == 0) {
    ecokit::stop_ctx("`path_coda` is of length 0")
  }

  cat(
    paste0("path_model:\n\t", path_model, "\npath_coda:\n\t", path_coda, "\n"))

  if (!all(fs::file_exists(c(path_model, path_coda)))) {
    ecokit::stop_ctx(
      "Selected model files not found",
      path_model = path_model, path_coda = path_coda, include_backtrace = TRUE)
  }

  # ****************************************************************

  if (spatial_model) {

    # Prepare response curve data -----

    if (rc_prepare) {
      ecokit::info_chunk(
        "Prepare response curve data", level = 1L, line_char = "+",
        line_char_rep = 60L, cat_red = TRUE,
        cat_bold = TRUE, cat_timestamp = FALSE)

      IASDT.R::rc_prepare_data(
        path_model = path_model, n_grid = n_grid, n_cores = n_cores_rc,
        strategy = strategy, future_max_size = future_max_size,
        use_tf = use_tf, tf_environ = tf_environ,
        tf_use_single = tf_use_single, n_cores_lf = n_cores_lf,
        lf_temp_cleanup = lf_temp_cleanup, lf_check = lf_check,
        temp_dir = temp_dir_pred, temp_cleanup = temp_cleanup, verbose = TRUE,
        lf_commands_only = FALSE, return_data = FALSE,
        probabilities = c(0.025, 0.5, 0.975))

      invisible(gc())
    }

    # ****************************************************************

    # Plotting response curves ----

    if (rc_plot) {

      ## Species richness -----
      ecokit::info_chunk(
        "Plotting response curves - species richness", level = 1L,
        line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)

      IASDT.R::rc_plot_sr(
        model_dir = model_dir, verbose = TRUE, n_cores = n_cores_rc,
        strategy = strategy, future_max_size = future_max_size)

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      # Individual species -----
      ecokit::info_chunk(
        "Plotting response curves - species", level = 1L, line_char = "+",
        line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)

      IASDT.R::rc_plot_species(
        model_dir = model_dir, n_cores = n_cores_rc, env_file = env_file)

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      # All species together -----
      ecokit::info_chunk(
        "Plotting response curves - all species together", level = 1L,
        line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)

      IASDT.R::rc_plot_species_all(
        model_dir = model_dir, n_cores = n_cores_rc)

      invisible(gc())
    }

    # ****************************************************************

    # Predicting habitat suitability across different climate options -------

    if (predict_suitability) {

      ecokit::info_chunk(
        "Predicting habitat suitability across different climate options",
        line_char = "+", line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE, level = 1L)

      # ************************************

      IASDT.R::predict_maps(
        path_model = path_model, hab_abb = hab_abb, env_file = env_file,
        n_cores = n_cores, n_cores_pred = n_cores_pred, strategy = strategy,
        future_max_size = future_max_size, clamp_pred = clamp_pred,
        fix_efforts = fix_efforts, fix_rivers = fix_rivers,
        pred_new_sites = pred_new_sites, use_tf = use_tf,
        tf_environ = tf_environ, tf_use_single = tf_use_single,
        n_cores_lf = n_cores_lf, lf_check = lf_check,
        lf_temp_cleanup = lf_temp_cleanup, lf_only = FALSE,
        lf_commands_only = FALSE, temp_dir = temp_dir_pred,
        temp_cleanup = temp_cleanup, tar_predictions = tar_predictions,
        climate_models = climate_models, climate_scenario = climate_scenario,
        spatial_model = spatial_model)

      invisible(gc())

      # ************************************

      # Plot species & SR predictions as JPEG ------

      if (plot_predictions) {
        ecokit::info_chunk(
          "Plot species & SR predictions as JPEG", level = 1L, line_char = "+",
          line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
          cat_timestamp = FALSE)

        IASDT.R::plot_prediction(
          model_dir = model_dir, env_file = env_file, n_cores = n_cores,
          is_cv_model = is_cv_model)
      }

      # Evaluate cross-validated models -----
      if (is_cv_model) {
        ecokit::info_chunk(
          "Evaluating explanatory and predictive power of the model",
          level = 1L, line_char = "+", line_char_rep = 60L, cat_red = TRUE,
          cat_bold = TRUE, cat_timestamp = FALSE)
        IASDT.R::mod_cv_evaluate(model_dir = model_dir, cv_type = "cv_dist")
      }
      invisible(gc())
      # ************************************

      # Tar prediction directories and delete TIFF files -----

      ecokit::cat_time(
        "Create tar file for each prediction option and delete TIFF files",
        level = 1L)

      path_prediction <- fs::path(model_dir, "model_prediction")
      path_prediction_clamp <- fs::path(path_prediction, "clamp")
      path_prediction_no_clamp <- fs::path(path_prediction, "no_clamp")
      l_dirs_preds <- c(
        list.dirs(path_prediction_clamp, recursive = FALSE),
        list.dirs(path_prediction_no_clamp, recursive = FALSE))

      if (length(l_dirs_preds) > 0L) {
        ecokit::set_parallel(
          n_cores = n_cores, show_log = FALSE, strategy = strategy,
          future_max_size = future_max_size)
        tar_pred_tiffs0 <- future.apply::future_lapply(
          X = l_dirs_preds, FUN = tar_pred_tiffs, delete_tiff = TRUE,
          overwrite_tar = TRUE, future.scheduling = Inf, future.seed = TRUE,
          future.packages = c("stringr", "fs"),
          future.globals = c("l_dirs_preds", "tar_pred_tiffs"))
        rm(tar_pred_tiffs0, l_dirs_preds, envir = environment())
        ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
      }
    }

    # ****************************************************************

    ## Plot latent factors as JPEG ------

    if (plot_lf) {
      ecokit::info_chunk(
        "Plot latent factors as JPEG", level = 1L, line_char = "+",
        line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)

      IASDT.R::plot_latent_factor(path_model = path_model, env_file = env_file)
    }

    # ****************************************************************

    ## Plot explanatory Power ------

    if (isFALSE(is_cv_model) && plot_internal_evaluation) {
      ecokit::info_chunk(
        "Plot explanatory Power", level = 1L, line_char = "+",
        line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
        cat_timestamp = FALSE)

      IASDT.R::plot_evaluation(model_dir = model_dir)
    }

  }

  # ****************************************************************

  # Compute variance partitioning ------

  if (vp_prepare) {
    ecokit::info_chunk(
      "Compute variance partitioning", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::variance_partitioning_compute(
      path_model = path_model, n_cores = n_cores_vp, use_tf = use_tf,
      tf_environ = tf_environ, tf_use_single = tf_use_single,
      temp_cleanup = temp_cleanup, chunk_size = 50L, verbose = TRUE,
      vp_file = "varpar", vp_commands_only = FALSE, temp_dir = temp_dir_vp)
  }

  # ****************************************************************

  # Plot Variance partitioning ------

  if (vp_plot) {
    ecokit::info_chunk(
      "Plot Variance partitioning", level = 1L, line_char = "+",
      line_char_rep = 60L, cat_red = TRUE, cat_bold = TRUE,
      cat_timestamp = FALSE)

    IASDT.R::variance_partitioning_plot(
      path_model = path_model, env_file = env_file, vp_file = "varpar",
      use_tf = use_tf, tf_environ = tf_environ, n_cores = n_cores_vp,
      width = 30, height = 15, spatial_model = spatial_model,
      is_cv_model = is_cv_model, temp_dir = temp_dir_vp)

  }

  # ****************************************************************

  ecokit::cat_diff(init_time = .start_time, prefix = "\nPost-processing took ")

  return(invisible(NULL))
}





#' Create a tar archive of .tif files in a directory and optionally delete
#' originals
#'
#' This function creates a tar archive containing all `.tif` files in the
#' specified directory. The archive is named after the last part of the
#' directory path. Optionally, the original `.tif` files can be deleted after
#' archiving. The function also allows control over whether to overwrite an
#' existing tar file.
#'
#' @param pred_dir Character. Path to the directory containing `.tif` files to
#'   be archived.
#' @param delete_tiff Logical. If `TRUE` (default), deletes the original `.tif`
#'   files after archiving.
#' @param overwrite_tar Logical. If `TRUE` (default), overwrites the existing
#'   tar file if it exists.
#'
#' @return Invisibly returns `NULL`. Used for its side effects of creating a tar
#'   archive and optionally deleting files.
#'
#' @details
#' - The tar archive is created in the same directory as `pred_dir` and is
#'   named after the last part of the directory.
#' - If no `.tif` files are found, the function exits with a message.
#' - File permissions for the tar archive are set to `755`.
#'
#' @examples
#' \dontrun{
#' tar_pred_tiffs(
#'   pred_dir = "path/to/predictions",
#'   delete_tiff = TRUE, overwrite_tar = FALSE)
#' }
#'
#' @noRd
#' @author Ahmed El-Gabbas

tar_pred_tiffs <- function(
    pred_dir = NULL, delete_tiff = TRUE, overwrite_tar = TRUE) {

  # get last directory of the path
  last_part <- stringr::str_extract(
    stringr::str_replace(pred_dir, "/+$", ""), "[^/]+$")
  pred_tar <- fs::path(pred_dir, paste0(last_part, ".tar"))

  if (fs::file_exists(pred_tar) && !overwrite_tar) {
    message(
      "The tar file already exists and `overwrite_tar` is `FALSE`: \n",
      "  >>>  ", pred_tar)
    return(invisible(NULL))
  }

  files <- list.files(path = pred_dir, pattern = "\\.tif$")
  if (length(files) == 0) {
    message("No .tif files found in the specified directory:\n", pred_dir)
    return(invisible(NULL))
  }

  # Command to create the tar file
  tar_command <- stringr::str_glue(
    'cd {fs::path_abs(pred_dir)}; tar -cf {basename(pred_tar)} \\
    -b 2048 {paste(files, collapse = " ")}')

  # Create tar file
  system(tar_command)

  # Change the permission of the tar file
  Sys.chmod(pred_tar, "755", use_umask = FALSE)

  # Delete all tiff files
  if (delete_tiff) {
    try(fs::file_delete(fs::path(pred_dir, files)), silent = TRUE)
  }
  invisible(NULL)
}
