## |------------------------------------------------------------------------| #
# predict_latent_factor ----
## |------------------------------------------------------------------------| #

#' Draws samples from the conditional predictive distribution of latent factors
#'
#' This function is optimized for speed using parallel processing and optionally
#' `TensorFlow` for matrix operations. This function is adapted from
#' [Hmsc::predictLatentFactor] with equivalent results to the original function
#' when `predictMean = TRUE`.
#' @param units_pred a factor vector with random level units for which
#'   predictions are to be made
#' @param units_model a factor vector with random level units that are
#'   conditioned on
#' @param post_eta Character. Path of `post_eta`; a list containing samples of
#'   random factors at conditioned units
#' @param post_alpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param lf_rl a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param n_cores_lf Integer. Number of cores to use for parallel processing of
#'   latent factor prediction. Defaults to 8L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param temp_dir Character. Path for temporary storage of intermediate files.
#' @param lf_temp_cleanup Logical. Whether to delete temporary files in the
#'   `temp_dir` directory after finishing the LF predictions.
#' @param model_name Character. Prefix for temporary file names. Defaults to
#'   `NULL`, in which case no prefix is used.
#' @param use_tf Logical. Whether to use `TensorFlow` for calculations. Defaults
#'   to `TRUE`.
#' @param tf_environ Character. Path to the Python environment. This argument is
#'   required if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.
#' @param lf_commands_only Logical. If `TRUE`, returns the command to run the
#'   Python script. Default is `FALSE`.
#' @param tf_use_single Logical. Whether to use single precision for the
#'   `TensorFlow` calculations. Defaults to `FALSE`.
#' @param lf_out_file Character. Path to save the outputs. If `NULL` (default),
#'   the predicted latent factors are not saved to a file. This should end with
#'   either `*.qs2` or `*.RData`.
#' @param lf_return Logical. Whether the output should be returned. Defaults to
#'   `FALSE`. If `lf_out_file` is `NULL`, this parameter cannot be set to
#'   `FALSE` because the function needs to return the result if it is not saved
#'   to a file.
#' @param lf_check Logical. If `TRUE`, the function checks if the output files
#'   are already created and valid. If `FALSE`, the function will only check if
#'   the files exist without checking their integrity. Default is `FALSE`.
#' @param verbose Logical. If `TRUE`, logs detailed information during
#'   execution. Default is `TRUE`.
#' @param solve_max_attempts Integer. Maximum number of attempts to run solve
#'   and crossprod internal function [run_crossprod_solve]. Default is 5L.
#' @param solve_chunk_size Integer. Chunk size for `solve_and_multiply` Python
#'   function. Default is 50L.
#' @export
#' @seealso [Hmsc::predictLatentFactor]
#' @name predict_latent_factor
#' @details The function is expected to be faster than the original function in
#'   the `Hmsc` package, especially when using `TensorFlow` for calculations and
#'   when working in parallel.
#'
#'   The main difference is that this function:
#' - allow for parallel processing (`n_cores_lf` argument);
#' - when `TensorFlow` is used (`use_tf = TRUE`), matrix
#'   calculations are much faster, particularly when used on GPU. The following
#'   Python modules are needed: `numpy`, `tensorflow`, `rdata`, `xarray`, and
#'   `pandas`. To use `TensorFlow` under Windows, the argument `tf_environ`
#'   should be set to the path of a Python environment with `TensorFlow`
#'   installed;
#' - if `use_tf` is set to `FALSE`, the function uses `R` (supported by
#'   relatively faster `CPP` functions) in the calculations;
#' - `d11` and `d12` matrices are processed only once and saved to disk and
#'   called when needed.

predict_latent_factor <- function(
    units_pred, units_model, post_eta, post_alpha, lf_rl, n_cores_lf = 8L,
    strategy = "multisession", temp_dir = "temp_pred",
    lf_temp_cleanup = TRUE, model_name = NULL, use_tf = TRUE, tf_environ = NULL,
    tf_use_single = FALSE, lf_out_file = NULL, lf_return = FALSE,
    lf_check = FALSE, lf_commands_only = FALSE, solve_max_attempts = 5L,
    solve_chunk_size = 50L, verbose = TRUE) {

  # # ..................................................................... ###

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- n_cores_lf <- 1L
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)

  .start_time <- lubridate::now(tzone = "CET")

  ecokit::cat_time(
    "Starting `predict_latent_factor` function", level = 1L, verbose = verbose)

  # # ..................................................................... ###

  # Check inputs

  if (is.null(lf_out_file) && isFALSE(lf_return)) {
    ecokit::stop_ctx(
      "`lf_return` must be `TRUE` when `lf_out_file` is NULL.",
      lf_return = lf_return, lf_out_file = lf_out_file,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble", "arrow", "ecokit",
      "Matrix", "Hmsc", "qs2", "fs", "purrr", "IASDT.R", "magrittr"),
    strategy = strategy)

  # # ..................................................................... ###

  # Load post_eta if it is a file path

  if (inherits(post_eta, "character")) {
    ecokit::cat_time("Load post_eta", level = 1L, verbose = verbose)
    if (!file.exists(post_eta)) {
      ecokit::stop_ctx(
        "The specified path for `post_eta` does not exist. ",
        post_eta = post_eta, include_backtrace = TRUE)
    }
    post_eta <- ecokit::load_as(post_eta)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  sample_id <- lf <- lf_id <- sample_ids <- alpha_id <- file_eta_pred <-
    n_samples <- chunk_id <- file_post_eta <- file_eta_pred_tf <- eta_data <-
    path_samp_lf <- NULL

  # # ..................................................................... ###

  # indices of units_pred in units_model
  ind_old <- (units_pred %in% units_model)
  # indices of new units_pred
  ind_new <- !(ind_old)

  # In the original Hmsc::predictLatentFactor function, the function is used
  # irrespective if the provided locations are new or not, then values at new
  # sites were replaced from model posterior. This can make the computations
  # complex when predicting at too many locations at large scale. In the
  # modified functions here, LF predictions are made only for new sites and
  # loaded at locations for model fitting. This distinction facilitate the
  # processing.

  all_training <- sum(ind_new) == 0
  all_new <- sum(ind_old) == 0

  # Either all_training or all_new should be TRUE
  if (sum(all_training, all_new) != 1) {
    ecokit::stop_ctx(
      "The input sites should be either all training sites or all new sites.",
      sum = sum(all_training, all_new), include_backtrace = TRUE)
  }

  if (all_training) {
    # If all input sites are for training sites, use LF info from the model
    # directly
    ecokit::cat_time(
      "All input sites are training sites", level = 1L, verbose = verbose)

    post_eta_pred <- purrr::map(
      .x = post_eta,
      .f = ~ {
        out <- .x[match(units_pred[ind_old], units_model), ]
        rownames(out) <- units_model
        return(out)
      })

  } else {

    ecokit::cat_time(
      "All input sites are new sites", level = 1L, verbose = verbose)

    # Check `TensorFlow` settings

    if (use_tf) {

      python_script <- system.file("crossprod_solve.py", package = "IASDT.R")

      # Check if python_script exists
      if (!file.exists(python_script)) {
        ecokit::stop_ctx(
          "Necessary Python script does not exist",
          python_script = python_script, include_backtrace = TRUE)
      }

      # Suppress `TensorFlow` warnings and disable optimizations
      Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0") # nolint: undesirable_function_linter

      ecokit::cat_time(
        "Computations will be made using `TensorFlow`", level = 1L,
        verbose = verbose)
    } else {
      ecokit::cat_time(
        "Computations will be made using R/CPP", level = 1L, verbose = verbose)
    }

    # # .................................................................... ###

    # Adjust model_name prefix

    if (is.null(model_name)) {
      model_name <- ""
    } else {
      model_name <- paste0(model_name, "_")
    }

    # ensure only one underscore at the end of model_name
    model_name <- stringr::str_replace_all(model_name, "__$", "_")

    # Create a temporary directory to store intermediate results. This directory
    # will be used to save s1/s2 or d11/d12, and intermediate post_eta files,
    # reducing memory usage.
    temp_dir_lf <- fs::path(temp_dir, "lf_prediction")
    fs::dir_create(c(temp_dir_lf, temp_dir))

    # # .................................................................... ###

    # Calculate d11 and d12 only once

    ecokit::cat_time(
      "Calculate/save necessary matrices", level = 1L, verbose = verbose)

    alphapw <- lf_rl$alphapw      # nolint: object_name_linter

    if (use_tf) {

      # Save s1 and s2 for coordinates at training and testing sites as feather
      # files, if not already exist on disk
      path_s1 <- fs::path(temp_dir, paste0(model_name, "s1.feather"))
      path_s2 <- fs::path(temp_dir, paste0(model_name, "s2.feather"))

      s1_s2_okay <- ecokit::check_data(path_s1, warning = FALSE) &&
        ecokit::check_data(path_s2, warning = FALSE)

      if (s1_s2_okay) {
        ecokit::cat_time(
          "s1 and s2 matrices were already saved",
          level = 2L, verbose = verbose, cat_timestamp = FALSE)
      } else {

        ecokit::cat_time(
          "Saving s1 and s2 matrices", level = 2L, verbose = verbose)

        # s1
        s1 <- as.data.frame(lf_rl$s[units_model, , drop = FALSE])
        ecokit::save_as(object = s1, out_path = path_s1)

        # s2
        s2 <- as.data.frame(lf_rl$s[units_pred[ind_new], , drop = FALSE])
        ecokit::save_as(object = s2, out_path = path_s2)

        rm(s1, s2, envir = environment())
      }

      rm(lf_rl, envir = environment())
      path_d11 <- path_d12 <- NULL

    } else {

      # Save d11 and d12 as feather files, if not already exist on disk
      path_d11 <- fs::path(temp_dir, paste0(model_name, "d11.qs2"))
      path_d12 <- fs::path(temp_dir, paste0(model_name, "d12.qs2"))

      if (file.exists(path_d11) && file.exists(path_d12)) {

        ecokit::cat_time(
          "d11 and d12 distance matrices are already saved",
          level = 2L, verbose = verbose)

      } else {

        s1 <- lf_rl$s[units_model, , drop = FALSE]
        s2 <- lf_rl$s[units_pred[ind_new], , drop = FALSE]

        # d11
        d11 <- Rfast::Dist(s1)
        ecokit::save_as(object = d11, out_path = path_d11)

        # d12
        d12 <- Rfast::dista(s1, s2)
        ecokit::save_as(object = d12, out_path = path_d12)

        # Clean up
        rm(lf_rl, s1, s2, d11, d12, envir = environment())

      }

      path_s1 <- path_s2 <- NULL

    }

    invisible(gc())

    # # .................................................................... ###

    # Convert post_alpha to tibble

    ecokit::cat_time(
      "Splitting and saving `post_alpha` to small chunks",
      level = 1L, verbose = verbose)

    # Unique combination of LF / alphapw / sample IDs
    lf_data <- do.call(rbind, post_alpha) %>%
      as.data.frame() %>%
      tibble::tibble() %>%
      stats::setNames(paste0("lf_", seq_len(ncol(.)))) %>%
      # ID column represents the original row number
      dplyr::mutate(sample_ids = dplyr::row_number()) %>%
      tidyr::nest(sample_ids = sample_ids) %>%
      dplyr::mutate(
        sample_ids = purrr::map(sample_ids, ~ as.vector(unlist(.x)))) %>%
      tidyr::pivot_longer(
        cols = -sample_ids, values_to = "alpha_id", names_to = "lf") %>%
      dplyr::mutate(
        lf_id = as.integer(stringr::str_remove(lf, "lf_")),
        denom = purrr::map_dbl(alpha_id, ~ alphapw[.x, 1])) %>%
      tidyr::nest(sample_ids = -c("denom", "lf_id", "alpha_id", "lf")) %>%
      dplyr::mutate(
        sample_ids = purrr::map(sample_ids, ~unname(sort(unlist(.x)))),
        # number of samples to be processed in each file
        n_samples = purrr::map_int(sample_ids, length)) %>%
      dplyr::arrange(dplyr::desc(n_samples)) %>%
      dplyr::mutate(
        chunk_id = dplyr::row_number(),
        file_post_eta = purrr::map_chr(
          .x = chunk_id,
          .f = ~ {
            chunk_id0 <- stringr::str_pad(
              .x, width = nchar(dplyr::n()), pad = 0)
            fs::path(
              temp_dir,
              paste0(model_name, "post_eta_ch", chunk_id0, ".feather"))
          }),
        file_eta_pred_tf = stringr::str_replace_all(
          file_post_eta, "_post_eta_ch", "_eta_pred_ch"),
        file_eta_pred = stringr::str_replace_all(
          file_eta_pred_tf, ".feather", ".qs2"),
        Export = purrr::pmap(
          .l = list(sample_ids, lf_id, file_post_eta),
          .f = function(sample_ids, lf_id, file_post_eta) {

            # do not export file if already exists
            if (!file.exists(file_post_eta)) {
              out <- post_eta[sample_ids] %>%
                purrr::map(~ .x[, lf_id, drop = FALSE]) %>%
                simplify2array() %>%
                as.data.frame() %>%
                stats::setNames(paste0("Eta_", sample_ids))

              ecokit::save_as(object = out, out_path = file_post_eta)
            }

            NULL

          }),
        Export = NULL)

    rm(post_eta, post_alpha, envir = environment())
    invisible(gc())

    # # .................................................................... ###

    # Internal functions to predict latent factors

    eta_preds_f <- function(row_num, units_pred, lf_check = FALSE) {

      # do not use scientific notation
      withr::local_options(scipen = 99)

      # Current denominator
      denom <- lf_data$denom[[row_num]]
      # ID for latent factor
      lf_id <- lf_data$lf_id[[row_num]]
      # ID for posterior sample
      sample_id <- lf_data$sample_ids[[row_num]]

      # File path for current data
      file_post_eta <- lf_data$file_post_eta[[row_num]]
      file_eta_pred <- lf_data$file_eta_pred[[row_num]]
      file_eta_pred_tf <- lf_data$file_eta_pred_tf[[row_num]]

      calc_pred_lf <- !file.exists(file_eta_pred)

      if (lf_check && isFALSE(calc_pred_lf)) {
        eta_n_cols <- ncol(ecokit::load_as(file_eta_pred))
        calc_pred_lf <- (eta_n_cols != lf_data$n_samples[[row_num]])
      }

      if (calc_pred_lf) {

        # If the denominator is positive, perform calculations; otherwise, set
        # `eta_indNew` to zero.

        if (denom > 0) {

          if (use_tf) {

            # Use `TensorFlow`

            # Suppress `TensorFlow` warnings and disable optimizations
            Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0") # nolint: undesirable_function_linter

            if (file.exists(file_eta_pred_tf)) {
              eta_ind_new_0 <- file_eta_pred_tf
            } else {
              eta_ind_new_0 <- run_crossprod_solve(
                tf_environ = tf_environ, s1 = path_s1, s2 = path_s2,
                denom = denom, post_eta = file_post_eta,
                path_out = file_eta_pred_tf,
                tf_use_single = tf_use_single,
                lf_commands_only = lf_commands_only,
                solve_max_attempts = solve_max_attempts,
                solve_chunk_size = solve_chunk_size)
            }

            if (lf_commands_only) {
              return(eta_ind_new_0)
            }

            rm(eta_ind_new_0, envir = environment())

            eta_pred <- tibble::tibble(
              sample_id = sample_id,
              lf = lf_id,
              path_samp_lf = fs::path(
                temp_dir_lf,
                paste0(
                  model_name, "Samp_",
                  stringr::str_pad(sample_id, width = 4, pad = "0"),
                  "_lf", lf, ".qs2")),
              eta_pred = as.list(ecokit::load_as(file_eta_pred_tf))) %>%
              tidyr::nest(eta_data = -path_samp_lf) %>%
              dplyr::mutate(
                eta_data = purrr::map(
                  .x = eta_data,
                  .f = ~ {
                    tidyr::unnest_longer(.x, "eta_pred") %>%
                      dplyr::mutate(lf = NULL, units_pred = units_pred) %>%
                      stats::setNames(
                        c("sample_id", paste0("lf_", lf_id), "units_pred"))
                  }),
                file_eta_pred = file_eta_pred,
                chunk_id = lf_data$chunk_id[[row_num]],
                sample_id = sample_id,
                Save = purrr::map2(
                  .x = eta_data, .y = path_samp_lf,
                  .f = ~qs2::qs_save(.x, .y, nthreads = 5)),
                lf = lf_id, Save = NULL, eta_data = NULL)

            ecokit::save_as(object = eta_pred, out_path = file_eta_pred)

          } else {

            # Use R / CPP

            # Reading post_eta from file
            post_eta0 <- ecokit::load_as(file_post_eta)

            # Read d11 and d12
            d11 <- ecokit::load_as(path_d11)
            d12 <- ecokit::load_as(path_d12)

            k11 <- IASDT.R::exp_neg_div(d11, denom)
            k12 <- IASDT.R::exp_neg_div(d12, denom)

            eta_pred <- purrr::map_chr(
              .x = seq_along(sample_id),
              .f = function(id) {

                path_samp_lf <- fs::path(
                  temp_dir_lf,
                  paste0(
                    model_name, "Samp_",
                    stringr::str_pad(id, width = 4, pad = "0"),
                    "_lf", lf_id, ".qs2"))

                  data <- as.matrix(post_eta0[, id]) %>%
                  IASDT.R::solve2vect(k11, .) %>%
                  as.vector() %>%
                  Matrix::crossprod(k12, .) %>%
                  as.vector() %>%
                  tibble::tibble(
                    sample_id = sample_id[id],
                    eta_pred = ., units_pred = units_pred) %>%
                  stats::setNames(
                    c("sample_id", paste0("lf_", lf_id), "units_pred"))

                qs2::qs_save(object = data, file = path_samp_lf, nthreads = 5)

                path_samp_lf

              }) %>%
              tibble::tibble(
                path_samp_lf = ., file_eta_pred = file_eta_pred,
                chunk_id = lf_data$chunk_id[[row_num]], sample_id = sample_id,
                lf = lf_id)

            ecokit::save_as(object = eta_pred, out_path = file_eta_pred)
          }

        } else {

          # When denom is zero, set `eta_indNew` to zero

          if (isFALSE(lf_commands_only)) {   # nolint: unneeded_nesting_linter

            eta_pred <- tibble::tibble(

              path_samp_lf = fs::path(
                temp_dir_lf,
                paste0(
                  model_name, "Samp_",
                  stringr::str_pad(sample_id, width = 4L, pad = "0"),
                  "_lf", lf_id, ".qs2")),

              file_eta_pred = file_eta_pred,
              chunk_id = lf_data$chunk_id[[row_num]],
              sample_id = sample_id) %>%
              dplyr::mutate(
                Save = purrr::map2(
                  .x = path_samp_lf, .y = sample_id,
                  .f = ~ {
                    data <- tibble::tibble(
                      sample_id = .y, lf = 0,
                      units_pred = units_pred) %>%
                      stats::setNames(
                        c("sample_id", paste0("lf_", lf_id), "units_pred"))
                    qs2::qs_save(object = data, file = .x, nthreads = 5)
                    return(NULL)
                  }),
                Save = NULL, lf = lf_id)

            ecokit::save_as(object = eta_pred, out_path = file_eta_pred)
          }
        }

      } else {
        eta_pred <- ecokit::load_as(file_eta_pred)
      }

      if (isFALSE(lf_commands_only)) {
        return(eta_pred)
      }
    }

    invisible(gc())

    # # .................................................................... ###
    # # .................................................................... ###

    # Predict latent factors

    if (all(file.exists(lf_data$file_eta_pred))) {
      ecokit::cat_time(
        "All lf prediction files were already created",
        level = 1L, verbose = verbose)
    } else {
      if (n_cores_lf == 1 || lf_commands_only) {

        if (lf_commands_only) {
          ecokit::cat_time(
            "Prepare commands for predicting latent factors",
            level = 1L, verbose = verbose)
        } else {
          # Sequential processing
          ecokit::cat_time(
            "Predicting Latent Factor sequentially",
            level = 1L, verbose = verbose)
        }

        # Making predictions sequentially
        eta_preds <- purrr::map(
          .x = seq_len(nrow(lf_data)),
          .f = function(x) {
            result <- try(
              expr = eta_preds_f(
                row_num = x, units_pred = units_pred, lf_check = lf_check),
              silent = FALSE)

            if (inherits(result, "try-error")) {
              NULL
            } else {
              result
            }

          })

        if (lf_commands_only) {

          command_file_prefix <- dplyr::if_else(
            startsWith(model_name, "rc_c_"),
            "lf_rc_commands_", "lf_new_sites_commands_")

          # Function to save commands to files
          save_commands_to_file <- function(commands, max_lines = 210) {
            # Determine how many files we need
            num_files <- ceiling(length(commands) / max_lines)

            # Loop through and save chunks of commands to separate files
            for (i in 1:num_files) {
              # Calculate the start and end line indices for this chunk
              start_line <- ((i - 1) * max_lines) + 1
              end_line <- min(i * max_lines, length(commands))

              # Get the chunk of commands
              chunk <- commands[start_line:end_line]

              # Define the filename
              file_name <- fs::path(
                temp_dir, paste0(command_file_prefix, i, ".txt"))

              # Write the chunk to a file with Linux line endings
              writeLines(chunk, file_name, useBytes = TRUE)
            }
          }

          # Call the function to save commands to file
          save_commands_to_file(unlist(eta_preds))

          return(NULL)
        }

      } else {

        # Parallel processing
        ecokit::cat_time(
          "Predicting Latent Factor in parallel", level = 1L, verbose = verbose)

        ecokit::set_parallel(
          n_cores = min(n_cores_lf, nrow(lf_data)), level = 2L,
          future_max_size = 800L, strategy = strategy, cat_timestamp = FALSE)
        withr::defer(future::plan("sequential", gc = TRUE))

        ecokit::cat_time(
          "Making predictions in parallel", level = 2L, verbose = verbose)
        eta_preds <- future.apply::future_lapply(
          X = seq_len(nrow(lf_data)),
          FUN = function(x) {
            result <- try(
              eta_preds_f(
                row_num = x, units_pred = units_pred, lf_check = lf_check),
              silent = FALSE)

            invisible(gc())
            if (inherits(result, "try-error")) {
              NULL
            } else {
              result
            }
          },
          future.seed = TRUE, future.packages = pkg_to_export,
          future.globals = c(
            "lf_data", "path_d11", "path_d12", "path_s1", "path_s2", "ind_new",
            "units_pred", "ind_old", "units_model", "tf_environ", "use_tf",
            "tf_use_single", "eta_preds_f", "lf_check", "run_crossprod_solve",
            "lf_commands_only", "solve_max_attempts", "solve_chunk_size",
            "temp_dir_lf"))

        # Stop the cluster
        ecokit::set_parallel(
          stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
      }

      # Check if all files are created
      ecokit::cat_time(
        "Check if all files are created", level = 1L, verbose = verbose)
      all_eta_files <- lf_data$file_eta_pred
      all_eta_files_exist <- all(file.exists(all_eta_files))

      if (!all_eta_files_exist) {
        failed_files <- all_eta_files[!file.exists(all_eta_files)]
        ecokit::stop_ctx(
          paste0(length(failed_files), " files are missing"),
          failed_files = basename(failed_files), include_backtrace = TRUE)
      }
      ecokit::cat_time(
        "All files were created", level = 2L, verbose = verbose,
        cat_timestamp = FALSE)

    }

    invisible(gc())

    # # .................................................................... ###

    # Merge results
    ecokit::cat_time("Merge results in parallel", level = 1L, verbose = verbose)

    post_eta_pred_samp <- eta_preds %>%
      dplyr::bind_rows() %>%
      dplyr::select(-lf, -chunk_id, -file_eta_pred) %>%
      dplyr::arrange(sample_id) %>%
      dplyr::mutate(
        path_sample = fs::path(
          temp_dir_lf,
          paste0(
            model_name, "Samp_",
            stringr::str_pad(sample_id, width = 4, pad = "0"), ".qs2"))) %>%
      tidyr::nest(data = -c("sample_id", "path_sample"))


    ecokit::set_parallel(
      n_cores = n_cores_lf, level = 2L, future_max_size = 800L,
      strategy = strategy, cat_timestamp = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))

    ecokit::cat_time(
      "Process results for MCMC samples in parallel",
      level = 2L, verbose = verbose)

    post_eta_pred <- future.apply::future_lapply(
      X = seq_len(nrow(post_eta_pred_samp)),
      FUN = function(x) {

        path_sample <- post_eta_pred_samp$path_sample[[x]]
        path_lf <- post_eta_pred_samp$data[[x]]$path_samp_lf

        sample_data_0 <- lapply(sort(path_lf), qs2::qs_read) %>%
          purrr::reduce(
            .f = dplyr::left_join, by = c("sample_id", "units_pred")) %>%
          dplyr::arrange(units_pred) %>%
          dplyr::select(sort(tidyselect::peek_vars())) %>%
          dplyr::select(-sample_id) %>%
          magrittr::set_rownames(NULL) %>%
          tibble::column_to_rownames("units_pred") %>%
          unname() %>%
          as.matrix()

        ecokit::save_as(sample_data_0, out_path = path_sample)
        try(fs::file_delete(path_lf), silent = TRUE)

        sample_data_0
      },
      future.seed = TRUE, future.packages = pkg_to_export,
      future.globals = c("post_eta_pred_samp", "lf_return"))

    # Stop the cluster
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)

    invisible(gc())

  }

  # # ..................................................................... ###

  # Save post_eta_pred

  if (!is.null(lf_out_file)) {
    ecokit::cat_time(
      "Saving post_eta_pred to disk", level = 1L, verbose = verbose)
    ecokit::cat_time(
      lf_out_file, cat_timestamp = FALSE, level = 2L, verbose = verbose)
    fs::dir_create(fs::path_dir(lf_out_file))
    ecokit::save_as(
      object = post_eta_pred, out_path = lf_out_file,
      n_threads = 5L, compress_level = 6L)

    ecokit::save_as(
      object = post_eta_pred_samp,
      out_path = stringr::str_replace(lf_out_file, ".qs2", "_Samp.qs2"))
  }

  # # ..................................................................... ###

  # Clean up temporary files after finishing calculations
  if (lf_temp_cleanup) {

    ecokit::cat_time(
      "Cleaning up temporary files", level = 1L, verbose = verbose)

    try(
      expr = {

        file_pattern <- paste0(
          "^", model_name,
          "(post_eta|r[0-9]|eta_pred|s1|s2|post).+(feather|qs2|log)$")
        file_paths <- list.files(
          path = ecokit::normalize_path(temp_dir),
          pattern = file_pattern, full.names = TRUE)
        if (length(file_paths) > 0) {
          try(fs::file_delete(file_paths), silent = TRUE)
        }

        file_paths_2 <- list.files(
          path = ecokit::normalize_path(temp_dir),
          pattern = "(lf_.+_test|rc_c)_Samp_.+.qs2",
          full.names = TRUE, recursive = TRUE)
        if (length(file_paths_2) > 0) {
          try(fs::file_delete(file_paths_2), silent = TRUE)
        }

        # delete temp files for cross-validated models
        file_paths_3 <- list.files(
          path = ecokit::normalize_path(temp_dir_lf),
          pattern = paste0("^", model_name, "Samp_.+.qs2"),
          full.names = TRUE, recursive = TRUE)

        if (length(file_paths_3) > 0) {
          try(fs::file_delete(file_paths_3), silent = TRUE)
        }

      },
      silent = TRUE)

  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "predict_latent_factor was finished in ",
    level = 1L, verbose = verbose)

  # # ..................................................................... ###

  if (lf_return) {
    return(post_eta_pred)
  } else {
    return(lf_out_file)
  }

}

# # ========================================================================== #
# # ========================================================================== #


## |------------------------------------------------------------------------| #
# run_crossprod_solve ----
## |------------------------------------------------------------------------| #

#' run_crossprod_solve
#'
#' Internal function to executes a Python script that performs matrix
#' computations using `TensorFlow` with provided inputs. Retries up to three
#' times if the output file validation fails.
#'
#' @param s1 Character. Path to the input file containing s1 coordinates.
#' @param s2 Character Path to the input file containing s2 coordinates.
#' @param post_eta Character. Path to the file containing the `post_eta` matrix
#'   data.
#' @param path_out Character. Path to rds file where the output results will be
#'   saved.
#' @param denom Numeric. The denominator value used in the computation.
#' @param chunk_size Numeric (Optional). Size of chunks to process at a time.
#'   Default is 1000.
#' @param threshold_mb Numeric (Optional). Memory threshold (in MB) to manage
#'   processing. Default is 2000.
#' @return Returns the `path_out` if successful. Returns `NULL` if all attempts
#'   fail.
#' @details
#' - The function checks for the existence of required input files and the
#' Python executable in the specified virtual environment.
#' - Executes the Python script using `system2`.
#' - Verifies the output file validity using `ecokit::check_data`. Retries up
#' to 3 times if the output is invalid.
#' - Generates detailed logs if `verbose` is set to `TRUE`.
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @inheritParams predict_latent_factor

run_crossprod_solve <- function(
    tf_environ, s1, s2, post_eta, path_out, denom,
    chunk_size = 1000L, threshold_mb = 2000L, tf_use_single = TRUE,
    verbose = TRUE, solve_chunk_size = 50L, solve_max_attempts = 5L,
    lf_commands_only = FALSE) {

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")     # nolint: undesirable_function_linter

  # do not use scientific notation
  withr::local_options(scipen = 99)

  script_path <- system.file("crossprod_solve.py", package = "IASDT.R")
  if (!file.exists(script_path)) {
    ecokit::stop_ctx(
      "Necessary Python script `crossprod_solve.py` does not exist",
      script_path = script_path, include_backtrace = TRUE)
  }

  # Ensure the paths are valid
  paths <- list(script_path, s1, s2, post_eta)
  names(paths) <- c("Python Script", "s1", "s2", "post_eta")
  purrr::walk(
    .x = names(paths),
    .f = function(p) {
      if (!file.exists(paths[[p]])) {
        ecokit::stop_ctx(
          paste0(p, " does not exist"), path = paths[[p]],
          include_backtrace = TRUE)
      }
    })

  # Determine the Python executable path

  # On Windows, the TF calculations has to be done through a valid virtual
  # environment; the path to the virtual environment must be specified in
  # `tf_environ`. On LUMI, this is not needed as the compatible Python
  # installation is loaded automatically when loading `tensorflow` module. When
  # using another HPC system, the function needs to be adapted accordingly.

  if (.Platform$OS.type == "windows") {

    if (is.null(tf_environ)) {
      ecokit::stop_ctx(
        "When running on Windows, `tf_environ` must be specified",
        tf_environ = tf_environ, include_backtrace = TRUE)
    }
    if (!dir.exists(tf_environ)) {
      ecokit::stop_ctx(
        "The specified `tf_environ` directory does not exist",
        tf_environ = ecokit::normalize_path(tf_environ),
        include_backtrace = TRUE)
    }

    python_executable <- ecokit::normalize_path(
      fs::path(tf_environ, "Scripts", "python.exe"), must_work = TRUE)

    if (!file.exists(python_executable)) {
      ecokit::stop_ctx(
        "Python executable not found in the virtual environment.",
        python_executable = python_executable, include_backtrace = TRUE)
    }

  } else {
    # Use `python3`` directly - on LUMI, compatible Python installation is
    # loaded automatically when loading `tensorflow`
    python_executable <- "/usr/bin/time -v python3" # nolint: absolute_paths_linter
  }

  # Construct the command to run the Python script
  lf_args <- c(
    python_executable,
    script_path,
    "--s1", ecokit::normalize_path(s1, must_work = TRUE),
    "--s2", ecokit::normalize_path(s2, must_work = TRUE),
    "--post_eta", ecokit::normalize_path(post_eta, must_work = TRUE),
    "--path_out", ecokit::normalize_path(path_out),
    "--denom", as.character(denom),
    "--chunk_size", as.character(chunk_size),
    "--threshold_mb", as.character(threshold_mb),
    "--solve_chunk_size", as.character(solve_chunk_size))

  # Add boolean flags conditionally
  if (tf_use_single) {
    lf_args <- c(lf_args, "--use_single")
  }

  if (verbose) {
    lf_args <- c(lf_args, "--verbose")
  }

  if (.Platform$OS.type != "windows") {
    path_log <- stringr::str_replace(
      fs::path(getwd(), path_out), ".feather", ".log")
    # Redirect results of time to log file
    lf_args <- c(lf_args, paste0(" >> ", path_log, " 2>&1"))
  }

  lf_args <- paste(lf_args, collapse = " ")

  if (lf_commands_only) {
    return(lf_args)
  } else {

    path_log <- stringr::str_replace(path_out, ".feather", ".log")
    f <- file(path_log, open = "a")
    on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
    cat(
      "Running command:\n", paste(lf_args, "\n\n"),
      sep = "\n", file = f, append = TRUE)

    # Initialise retry logic
    attempt <- 1
    success <- FALSE

    while (attempt <= solve_max_attempts && !success) {

      cat(
        paste0("Attempt ", attempt, " of ", solve_max_attempts, "\n\n"),
        file = f, append = TRUE)

      # Run the command and capture stdout/stderr to a log file
      result <- system(lf_args, intern = TRUE)

      # Check for errors
      if (!inherits(result, "error") || length(result) != 0 ||
          result == "Done") {
        # Check the integrity of the file
        success <- ecokit::check_data(path_out, warning = FALSE)
      }
      attempt <- attempt + 1
    }

    # If all attempts fail, return NULL
    if (success) {
      # close connection to the file
      close(f)
      return(path_out)
    } else {
      if (verbose) {
        cat("All attempts failed. Returning NULL.\n",
            sep = "\n", file = f, append = TRUE)
        # close connection to the file
        close(f)
      }
      return(NULL)
    }
  }
}

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# plot_latent_factor ----
## |------------------------------------------------------------------------| #

#' Plot spatial variation in site loadings of HMSC models
#'
#' Generate and save spatial variation in site loadings of HMSC models' latent
#' factors as a JPEG file.
#'
#' @param path_model Character. Path to fitted `Hmsc` model object.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @export
#' @author Ahmed El-Gabbas
#' @name plot_latent_factor

plot_latent_factor <- function(path_model = NULL, env_file = ".env") {

  # # ..................................................................... ###

  path_grid <- NULL

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  grid_10 <- fs::path(path_grid, "grid_10_land_crop.RData") %>%
    ecokit::load_as(unwrap_r = TRUE)

  # # ..................................................................... ###

  # Check if the model file exists
  if (is.null(path_model) || !file.exists(path_model)) {
    ecokit::stop_ctx(
      "Selected model files not found", path_model = path_model,
      include_backtrace = TRUE)
  }

  model_obj <- ecokit::load_as(path_model)
  model_coords <- model_obj$ranLevels$sample$s
  post_eta <- Hmsc::getPostEstimate(model_obj, parName = "Eta")
  n_lf <- ncol(post_eta$mean)
  eta_mean <- as.data.frame(post_eta$mean) %>%
    stats::setNames(paste("lf", seq_len(n_lf), sep = "_")) %>%
    cbind.data.frame(model_coords, .) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)
  eta_mean_r <- terra::rasterize(
    eta_mean, grid_10, field = names(eta_mean)[-ncol(eta_mean)],
    fun = "mean") %>%
    stats::setNames(stringr::str_remove(names(.), "_mean")) %>%
    stats::setNames(stringr::str_replace(names(.), "lf_", "Latent factor "))
  rm(model_obj, model_coords, eta_mean, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  x_lim <- c(2600000, 6700000)
  y_lim <- c(1450000, 5420000)

  n_layers <- terra::nlyr(eta_mean_r)
  if (n_layers > 6) {
    ncols <- 4
  } else if (n_layers > 4) {
    ncols <- 3
  } else {
    ncols <- 2
  }

  lf_plot <- ggplot2::ggplot(environment = emptyenv()) +
    tidyterra::geom_spatraster(data = eta_mean_r, maxcell = Inf) +
    ggplot2::facet_wrap(~lyr, ncol = ncols, nrow = 2) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks(), name = NULL) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0.05, "cm"),
      plot.title = ggplot2::element_text(
        size = 12, color = "blue", face = "bold", hjust = 0.5,
        margin = ggplot2::margin(0, 0, 0, 0)),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "transparent", color = "transparent"),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.position = "inside",
      legend.position.inside = c(0.94, 0.9),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 8),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.title = ggplot2::element_text(
        color = "blue", size = 7, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.text.y = ggplot2::element_text(size = 7, hjust = 0.5, angle = 90),
      axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.04, "cm"),
      panel.spacing = grid::unit(0.3, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.1, colour = "grey40", linetype = 2),
      panel.border = ggplot2::element_blank(),
      panel.ontop = TRUE, panel.background = ggplot2::element_rect(fill = NA))

  # Adapt height proportionally based on the number of layers
  plot_height <- 20L
  if (ncols == 4) plot_width <- 40L
  if (ncols == 3) plot_width <- 30L
  if (ncols == 2) plot_width <- 20L

  ragg::agg_jpeg(
    filename = fs::path(
      dirname(dirname(path_model)), "model_prediction", "lf_plot.jpeg"),
    width = plot_width, height = plot_height, res = 600L,
    quality = 100L, units = "cm")
  print(lf_plot)
  grDevices::dev.off()

  return(invisible(NULL))

}
