## |------------------------------------------------------------------------| #
# rc_prepare_data ----
## |------------------------------------------------------------------------| #

#' Prepare and plot response curve data for Hmsc models
#'
#' The `rc_*()` functions process and visualise response curves for Hmsc
#' models. They support parallel computation and optionally return processed
#' data. There are four functions in this group:
#' - `rc_prepare_data()`: Prepares response curve data for analysis
#' - `rc_plot_species()`: Generates response curve plots for
#' individual species
#' - `rc_plot_species_all()`: Generates response curves for all
#' species together in a single plot
#' - `rc_plot_sr()`: Plots response curves for species richness.
#' @param path_model Character. Path to the file containing the fitted Hmsc
#'   model.
#' @param n_grid Integer. Number of points along the gradient for continuous
#'   focal variables. Higher values result in smoother curves. Default: 50. See
#'   [Hmsc::constructGradient] for details.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to 8L for all functions, except for `rc_plot_species`, in
#'   which it defaults to 20L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param return_data Logical. If `TRUE`, the function returns processed data as
#'   an R object. Default: `FALSE`.
#' @param probabilities Numeric vector. Quantiles to calculate in response curve
#'   predictions. Default: `c(0.025, 0.5, 0.975)`. See [stats::quantile] for
#'   details.
#' @param model_dir Character. Path to the root directory containing fitted
#'   models. The function reads data from the `response_curves_data`
#'   subdirectory, which is created by `rc_prepare_data`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param plotting_alpha Numeric. Opacity level for response curve lines (0 =
#'   fully transparent, 1 = fully opaque). Default: 0.3.
#' @export
#' @inheritParams predict_hmsc
#' @rdname response_curves
#' @name response_curves
#' @importFrom foreach %dopar%
#' @order 1
#' @author Ahmed El-Gabbas

rc_prepare_data <- function(
    path_model = NULL, n_grid = 50L, n_cores = 8L, strategy = "multisession",
    return_data = FALSE, probabilities = c(0.025, 0.5, 0.975), use_tf = TRUE,
    tf_environ = NULL, tf_use_single = FALSE, n_cores_lf = n_cores,
    lf_check = FALSE, lf_temp_cleanup = TRUE, lf_commands_only = FALSE,
    temp_dir = "temp_pred", temp_cleanup = TRUE, verbose = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::cat_time("Check input arguments", verbose = verbose)

  ecokit::check_args(
    args_to_check = c("path_model", "temp_dir"), args_type = "character")
  ecokit::check_args(
    args_to_check = c("n_grid", "probabilities"),
    args_type = "numeric", arg_length = c(1L, 3L))
  ecokit::check_args(
    args_to_check = c(
      "lf_check", "use_tf", "verbose", "lf_temp_cleanup",
      "lf_commands_only", "temp_cleanup"),
    args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- n_cores_lf <- 1L
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  rc_data <- variable <- rc_data_name <- sample_id <- species <- sr <- mm <-
    nfv <- path_rc_data_orig <- var_name <- path_rc_data_prob <-
    path_rc_data_sr <- coords <- NULL

  if (any(probabilities > 1) || any(probabilities < 0)) {
    ecokit::stop_ctx(
      "`probabilities` must be between 0 and 1", probabilities = probabilities,
      include_backtrace = TRUE)
  }
  probabilities <- sort(probabilities)

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "purrr", "tidyr", "abind", "Hmsc", "parallel", "magrittr"),
    strategy = strategy)

  # # ..................................................................... ###

  # Loading model object ------

  ecokit::cat_time("Loading model object", verbose = verbose)
  if (file.exists(path_model)) {
    model_obj <- ecokit::load_as(path_model)
    if (!inherits(model_obj, "Hmsc")) {
      ecokit::stop_ctx(
        "Model object is not of class 'hmsc'",
        class_model = class(model_obj), include_backtrace = TRUE)
    }
  } else {
    ecokit::stop_ctx(
      "The model file does not exist or is not a `.RData` or `.qs2` file.",
      path_model = path_model, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # prep_rc_data -------

  prep_rc_data <- function(id, file_lf) {

    variable <- rc_data$var_name[[id]]
    rc_data_name <- rc_data$rc_data_name[[id]]
    coords <- rc_data$coords[[id]]
    nfv <- rc_data$nfv[[id]]

    # Path for original prediction values
    path_rc_data_orig <- rc_data$path_rc_data_orig[[id]]
    # Path for plotting data: probability of occurrence
    path_rc_data_prob <- rc_data$path_rc_data_prob[[id]]
    path_rc_data_prob_samples <- stringr::str_replace(
      path_rc_data_prob, ".qs2$", "_samples.qs2")

    # Path for plotting data: species richness
    path_rc_data_sr <- rc_data$path_rc_data_sr[[id]]
    path_rc_data_sr_samples <- stringr::str_replace(
      path_rc_data_sr, ".qs2$", "_samples.qs2")

    output_tibble <- tibble::tibble(
      variable = variable, nfv = nfv, coords = coords,
      path_rc_orig = path_rc_data_orig,
      path_rc_prob = path_rc_data_prob,
      path_rc_data_prob_samples = path_rc_data_prob_samples,
      rc_path_sr = path_rc_data_sr,
      path_rc_data_sr_samples = path_rc_data_sr_samples)

    out_files_exists <- c(
      path_rc_data_orig, path_rc_data_prob, path_rc_data_sr) %>%
      file.exists() %>%
      all()

    if (isFALSE(out_files_exists)) {

      if (file.exists(path_rc_data_orig)) {

        rc_data_orig <- ecokit::load_as(path_rc_data_orig)
        gradient <- rc_data_orig$gradient
        x_vals <- gradient$XDataNew[, variable]
        preds <- rc_data_orig$preds
        pred_sr <- rc_data_orig$pred_sr

      } else {

        model_obj <- ecokit::load_as(path_model)

        # constructGradient
        gradient <- Hmsc::constructGradient(
          hM = model_obj, focalVariable = variable, non.focalVariables = nfv,
          ngrid = n_grid, coordinates = list(sample = coords))

        # Values of the current predictor
        x_vals <- gradient$XDataNew[, variable]

        rm(model_obj, envir = environment())

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        # Predicting probability of occurrence
        preds <- IASDT.R::predict_hmsc(
          path_model = path_model, gradient = gradient, expected = TRUE,
          n_cores = 1, strategy = strategy, model_name = paste0("rc_", coords),
          prediction_type = coords, use_tf = use_tf, tf_environ = tf_environ,
          lf_input_file = file_lf, n_cores_lf = 1, lf_check = lf_check,
          lf_temp_cleanup = lf_temp_cleanup, lf_commands_only = FALSE,
          tf_use_single = tf_use_single, temp_dir = temp_dir,
          temp_cleanup = temp_cleanup, verbose = FALSE)

        # species richness
        pred_sr <- abind::abind(lapply(preds, rowSums), along = 2)

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # Save gradient and original prediction values
        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        rc_data_orig <- list(
          variable = variable, nfv = rc_data$nfv[[id]],
          gradient = gradient, preds = preds, pred_sr = pred_sr)

        ecokit::save_as(
          object = rc_data_orig, object_name = paste0(rc_data_name, "_orig"),
          out_path = path_rc_data_orig)
      }

      rm(rc_data_orig, envir = environment())
      invisible(gc())

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: probability of occurrence
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      model_obj <- ecokit::load_as(path_model)

      rc_data_prob <- purrr::map_dfr(
        .x = seq_len(length(preds)),
        .f = function(sample) {
          tibble::as_tibble(preds[[sample]]) %>%
            dplyr::mutate(x_vals = x_vals, sample_id = sample)
        }) %>%
        tidyr::pivot_longer(
          cols = c(-x_vals, -sample_id),
          names_to = "species", values_to = "pred") %>%
        dplyr::arrange(species, x_vals, sample_id) %>%
        tidyr::nest(samples_data = -species) %>%
        dplyr::mutate(
          plot_data_quant = purrr::map(
            .x = samples_data,
            .f = ~ {
              dplyr::reframe(
                .x, pred = stats::quantile(pred, probabilities),
                quantile = probabilities, .by = x_vals)
            }),

          # Values at observed presence and absences
          observed_pa = purrr::map(
            .x = species,
            .f = ~ tibble::tibble(
              x_vals = model_obj$XData[, variable], pred = model_obj$Y[, .x])),

          # Positive trend probability
          positive_trend_prob = purrr::map_dbl(
            .x = samples_data,
            .f = ~ {
              dplyr::group_by(.x, sample_id) %>%
                dplyr::reframe(mm = dplyr::last(pred) > dplyr::first(pred)) %>%
                dplyr::pull(mm) %>%
                mean()
            }),
          variable = variable, nfv = rc_data$nfv[[id]], .before = 1)

      # Save data
      rc_data_prob_samples <- rc_data_prob
      ecokit::save_as(
        object = rc_data_prob_samples,
        object_name = paste0(rc_data_name, "_prob_samples"),
        out_path = path_rc_data_prob_samples)

      rc_data_prob <- dplyr::select(rc_data_prob, -samples_data)
      ecokit::save_as(
        object = rc_data_prob, object_name = paste0(rc_data_name, "_prob"),
        out_path = path_rc_data_prob)


      rm(preds, rc_data_prob, rc_data_prob_samples, envir = environment())

      # CHECK
      # rm(samples_data, envir = environment())


      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: species richness
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      # predicted species richness
      samples_data <- purrr::map_dfr(
        .x = seq_len(ncol(pred_sr)),
        .f = ~ tibble::tibble(
          x_vals = x_vals, sample_id = .x, sr = pred_sr[, .x])) %>%
        dplyr::arrange(x_vals, sample_id)

      # Quantiles of species richness
      rc_data_sr_quant <- dplyr::reframe(
        samples_data,
        sr = stats::quantile(sr, probabilities),
        quantile = probabilities, .by = x_vals)

      # Trend of the species richness
      sr_positive_trend_prob <- samples_data %>%
        dplyr::group_by(sample_id) %>%
        dplyr::reframe(mm = dplyr::last(sr) > dplyr::first(sr)) %>%
        dplyr::pull(mm) %>%
        mean()

      # Values at observed species richness
      observed_sr <- tibble::tibble(
        x_vals = model_obj$XData[, variable],
        pred = rowSums(model_obj$Y, na.rm = TRUE))

      # Save species richness data
      rc_data_sr <- list(
        variable = variable, nfv = rc_data$nfv[[id]],
        rc_data_sr_quant = rc_data_sr_quant, observed_sr = observed_sr,
        sr_positive_trend_prob = sr_positive_trend_prob)

      ecokit::save_as(
        object = rc_data_sr, object_name = paste0(rc_data_name, "_sr"),
        out_path = path_rc_data_sr)

      rc_data_sr_samples <- list(
        variable = variable, nfv = rc_data$nfv[[id]],
        rc_data_sr = rc_data_sr,
        rc_data_sr_quant = rc_data_sr_quant, observed_sr = observed_sr,
        sr_positive_trend_prob = sr_positive_trend_prob)

      ecokit::save_as(
        object = rc_data_sr_samples,
        object_name = paste0(rc_data_name, "_sr_samples"),
        out_path = path_rc_data_sr_samples)


      rm(
        rc_data_sr, rc_data_sr_quant, observed_sr,
        rc_data_sr_samples, sr_positive_trend_prob, envir = environment())

      # CHECK
      # rm(pred_sr, envir = environment())
      invisible(gc())
    }

    invisible(gc())

    output_tibble
  }

  # # ..................................................................... ###

  # Prepare response curve data -------

  path_rc <- fs::path(dirname(dirname(path_model)), "model_postprocessing")
  path_rc_data <- fs::path(path_rc, "response_curves_data")
  fs::dir_create(path_rc_data)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract names of the variables
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Extract names of the variables", verbose = verbose)
  model_vars <- stringr::str_split(
    as.character(model_obj$XFormula)[2], "\\+", simplify = TRUE) %>%
    stringr::str_trim()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Prediction variants
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # `coords`: Value of the `coordinates` argument of the `constructGradient`
  # function. `coordinates = "c"` for mean of coordinates (default);
  # `coordinates = "i"` for infinite coordinates without effect of spatial
  # dependence.
  #
  # nfv: Value of the `non.focalVariables` argument of `constructGradient`.
  # non.focalVariables = 1 sets the values of the non-focal variable to the most
  # likely value (defined as expected value for covariates, mode for factors).
  # non.focalVariables = 2 sets the values of the non-focal variable to most
  # likely value, given the value of focal variable, based on a linear
  # relationship. non.focalVariables = 3 fixes to the value given

  rc_data <- tidyr::expand_grid(
    variable = model_vars, coords = c("c", "i"), nfv = c(1, 2)) %>%
    dplyr::mutate(
      var_name = purrr::map_chr(
        .x = variable, .f = stringr::str_remove_all,
        pattern = "stats::poly\\(|, degree = 2, raw = TRUE\\)"),
      rc_data_name = paste0("rc_", var_name, "_coord_", coords, "_nfv", nfv),
      path_rc_data_orig = fs::path(
        path_rc_data, paste0(rc_data_name, "_orig.qs2")),
      path_rc_data_prob = fs::path(
        path_rc_data, paste0(rc_data_name, "_prob.qs2")),
      path_rc_data_sr = fs::path(path_rc_data, paste0(rc_data_name, "_sr.qs2")),
      FileExists = purrr::pmap_lgl(
        .l = list(path_rc_data_orig, path_rc_data_prob, path_rc_data_sr),
        .f = function(path_rc_data_orig, path_rc_data_prob, path_rc_data_sr) {
          c(path_rc_data_orig, path_rc_data_prob, path_rc_data_sr) %>%
            file.exists() %>%
            all()
        }))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Checking file existence
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  missing_rows <- sum(!rc_data$FileExists)
  file_lf <- fs::path(path_rc_data, "ResCurv_lf.qs2")

  if (missing_rows == 0) {

    ecokit::cat_time(
      "All response curve data files were already available on disk",
      level = 1L, verbose = verbose, cat_timestamp = FALSE)
    rc_data <- purrr::map_dfr(
      .x = seq_len(nrow(rc_data)), .f = prep_rc_data, file_lf = file_lf)

  } else {

    if (any(rc_data$FileExists)) {
      ecokit::cat_time(
        paste0(
          "Some response curve data files (", missing_rows, " of ",
          length(rc_data$FileExists), ") were missing"),
        level = 1L, verbose = verbose)
    } else {
      ecokit::cat_time(
        paste0(
          "All response curve data (", missing_rows, ") need to be prepared"),
        level = 1L, verbose = verbose, cat_timestamp = FALSE)
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Get LF prediction for the model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (isFALSE(ecokit::check_data(file_lf, warning = FALSE))) {

      ecokit::info_chunk(
        message = "Get LF prediction at mean coordinates", cat_date = FALSE,
        cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE,
        verbose = verbose)

      ecokit::cat_time("Create gradient", verbose = verbose)
      gradient_c <- Hmsc::constructGradient(
        hM = model_obj, focalVariable = rc_data$variable[1],
        non.focalVariables = 1, ngrid = 20, coordinates = list(sample = "c"))

      # The `model_obj` object is distributed twice to cores when available on
      # the function environment. Here, I delete the model_obj object and it
      # will be loaded later after when using `predict_hmsc` function.
      rm(model_obj, envir = environment())
      invisible(gc())

      ecokit::cat_time("Predicting LF", verbose = verbose)
      model_lf <- IASDT.R::predict_hmsc(
        path_model = path_model, gradient = gradient_c, expected = TRUE,
        n_cores = n_cores, strategy = strategy, temp_dir = temp_dir,
        temp_cleanup = temp_cleanup, model_name = "rc_c",
        prediction_type = "c", use_tf = use_tf, tf_environ = tf_environ,
        lf_out_file = file_lf, n_cores_lf = n_cores_lf, lf_check = lf_check,
        lf_return = FALSE, lf_only = TRUE, lf_temp_cleanup = lf_temp_cleanup,
        lf_commands_only = lf_commands_only, tf_use_single = tf_use_single,
        verbose = verbose, pred_directory = temp_dir)

      if (lf_commands_only) {
        return(invisible(NULL))
      }

      rm(model_lf, gradient_c, envir = environment())
      invisible(gc())

    } else {
      ecokit::cat_time(
        paste0(
          "LF prediction will be loaded from available file: \n   >>>  ",
          file_lf),
        verbose = verbose)
    }


    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare working in parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::info_chunk(
      message = "Prepare response curve data", cat_date = FALSE,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE, verbose = verbose)

    n_cores <- max(min(n_cores, missing_rows), 1)

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, future_max_size = 800L, strategy = strategy,
        cat_timestamp = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare response curve data in parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::cat_time(
      "Prepare response curve data in parallel", verbose = verbose)

    rc_data <- future.apply::future_lapply(
      X = seq_len(nrow(rc_data)),
      FUN = prep_rc_data, file_lf = file_lf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c(
        "rc_data", "path_model", "prep_rc_data", "n_grid", "probabilities",
        "file_lf", "use_tf", "tf_environ", "temp_dir", "lf_check",
        "lf_commands_only")) %>%
      dplyr::bind_rows()

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, cat_timestamp = FALSE)
      future::plan("sequential", gc = TRUE)
    }

    invisible(gc())
  }

  # # ..................................................................... ###

  ecokit::cat_time("Saving data to desk", verbose = verbose)
  save(rc_data, file = fs::path(path_rc_data, "rc_data.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Preparing response curve data took ",
    verbose = verbose)

  if (return_data) {
    return(rc_data)
  } else {
    return(invisible(NULL))
  }
}
