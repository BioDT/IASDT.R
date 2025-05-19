## |------------------------------------------------------------------------| #
# predict_latent_factor ----
## |------------------------------------------------------------------------| #

#' Draws samples from the conditional predictive distribution of latent factors
#'
#' This function is optimized for speed using parallel processing and optionally
#' TensorFlow for matrix operations. This function is adapted from
#' [Hmsc::predictLatentFactor] with equivalent results to the original function
#' when `predictMean = TRUE`.
#' @param units_pred a factor vector with random level units for which
#'   predictions are to be made
#' @param units_model a factor vector with random level units that are
#'   conditioned on
#' @param postEta Character. Path of `postEta`; a list containing samples of
#'   random factors at conditioned units
#' @param post_alpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param LF_rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param LF_n_cores Integer. Number of cores to use for parallel processing of
#'   latent factor prediction. Defaults to 8L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param temp_dir Character. Path for temporary storage of intermediate files.
#' @param LF_temp_cleanup Logical. Whether to delete temporary files in the
#'   `temp_dir` directory after finishing the LF predictions.
#' @param model_name Character. Prefix for temporary file names. Defaults to
#'   `NULL`, in which case no prefix is used.
#' @param use_TF Logical. Whether to use TensorFlow for calculations. Defaults
#'   to `TRUE`.
#' @param TF_environ Character. Path to the Python environment. This argument is
#'   required if `use_TF` is `TRUE` under Windows. Defaults to `NULL`.
#' @param LF_commands_only Logical. If `TRUE`, returns the command to run the
#'   Python script. Default is `FALSE`.
#' @param TF_use_single Logical. Whether to use single precision for the
#'   TensorFlow calculations. Defaults to `FALSE`.
#' @param LF_out_file Character. Path to save the outputs. If `NULL` (default),
#'   the predicted latent factors are not saved to a file. This should end with
#'   either `*.qs2` or `*.RData`.
#' @param LF_return Logical. Whether the output should be returned. Defaults to
#'   `FALSE`. If `LF_out_file` is `NULL`, this parameter cannot be set to
#'   `FALSE` because the function needs to return the result if it is not saved
#'   to a file.
#' @param LF_check Logical. If `TRUE`, the function checks if the output files
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
#' - allow for parallel processing (`LF_n_cores` argument);
#' - when `TensorFlow` is used (`use_TF = TRUE`), matrix
#'   calculations are much faster, particularly when used on GPU. The following
#'   Python modules are needed: `numpy`, `tensorflow`, `rdata`, `xarray`, and
#'   `pandas`. To use `TensorFlow` under Windows, the argument `TF_environ`
#'   should be set to the path of a Python environment with `TensorFlow`
#'   installed;
#' - if `use_TF` is set to `FALSE`, the function uses `R` (supported by
#'   relatively faster `CPP` functions) in the calculations;
#' - `D11` and `D12` matrices are processed only once and saved to disk and
#'   called when needed.

predict_latent_factor <- function(
    units_pred, units_model, postEta, post_alpha, LF_rL, LF_n_cores = 8L,
    strategy = "multisession", temp_dir = "TEMP_Pred",
    LF_temp_cleanup = TRUE, model_name = NULL, use_TF = TRUE, TF_environ = NULL,
    TF_use_single = FALSE, LF_out_file = NULL, LF_return = FALSE,
    LF_check = FALSE, LF_commands_only = FALSE, solve_max_attempts = 5L,
    solve_chunk_size = 50L, verbose = TRUE) {

  # # ..................................................................... ###

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  .start_time <- lubridate::now(tzone = "CET")

  ecokit::cat_time(
    "Starting `predict_latent_factor` function", level = 1L, verbose = verbose)

  # # ..................................................................... ###

  # Check inputs

  if (is.null(LF_out_file) && isFALSE(LF_return)) {
    ecokit::stop_ctx(
      "`LF_return` must be `TRUE` when `LF_out_file` is NULL.",
      LF_return = LF_return, LF_out_file = LF_out_file,
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

  # Load postEta if it is a file path

  if (inherits(postEta, "character")) {
    ecokit::cat_time("Load postEta", level = 1L, verbose = verbose)
    if (!file.exists(postEta)) {
      ecokit::stop_ctx(
        "The specified path for `postEta` does not exist. ", postEta = postEta,
        include_backtrace = TRUE)
    }
    postEta <- ecokit::load_as(postEta)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SampleID <- LF <- LF_ID <- Sample_IDs <- Alpha_ID <- File_etaPred <-
    NSamples <- ChunkID <- File_postEta <- File_etaPred_TF <- eta_DT <-
    Path_Samp_LF <- NULL

  # # ..................................................................... ###

  # indices of units_pred in units_model
  indOld <- (units_pred %in% units_model)
  # indices of new units_pred
  indNew <- !(indOld)

  # In the original Hmsc::predictLatentFactor function, the function is used
  # irrespective if the provided locations are new or not, then values at new
  # sites were replaced from model posterior. This can make the computations
  # complex when predicting at too many locations at large scale. In the
  # modified functions here, LF predictions are made only for new sites and
  # loaded at locations for model fitting. This distinction facilitate the
  # processing.

  AllTraining <- sum(indNew) == 0
  AllNew <- sum(indOld) == 0

  # Either AllTraining or AllNew should be TRUE
  if (sum(AllTraining, AllNew) != 1) {
    ecokit::stop_ctx(
      "The input sites should be either all training sites or all new sites.",
      sum = sum(AllTraining, AllNew), include_backtrace = TRUE)
  }

  if (AllTraining) {
    # If all input sites are for training sites, use LF info from the model
    # directly
    ecokit::cat_time(
      "All input sites are training sites", level = 1L, verbose = verbose)

    postEtaPred <- purrr::map(
      .x = postEta,
      .f = ~ {
        Out <- .x[match(units_pred[indOld], units_model), ]
        rownames(Out) <- units_model
        return(Out)
      })

  } else {

    ecokit::cat_time(
      "All input sites are new sites", level = 1L, verbose = verbose)

    # Check TensorFlow settings

    if (use_TF) {

      PythonScript <- system.file("crossprod_solve.py", package = "IASDT.R")

      # Check if PythonScript exists
      if (!file.exists(PythonScript)) {
        ecokit::stop_ctx(
          "Necessary Python script does not exist",
          PythonScript = PythonScript, include_backtrace = TRUE)
      }

      # Suppress TensorFlow warnings and disable optimizations
      Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

      ecokit::cat_time(
        "Computations will be made using TensorFlow", level = 1L,
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
    # will be used to save s1/s2 or D11/D12, and intermediate postEta files,
    # reducing memory usage.
    Temp_Dir_LF <- fs::path(temp_dir, "LF_Prediction")
    fs::dir_create(c(Temp_Dir_LF, temp_dir))

    # # .................................................................... ###

    # Calculate D11 and D12 only once

    ecokit::cat_time(
      "Calculate/save necessary matrices", level = 1L, verbose = verbose)

    alphapw <- LF_rL$alphapw      # nolint: object_name_linter

    if (use_TF) {

      # Save s1 and s2 for coordinates at training and testing sites as feather
      # files, if not already exist on disk
      Path_s1 <- fs::path(temp_dir, paste0(model_name, "s1.feather"))
      Path_s2 <- fs::path(temp_dir, paste0(model_name, "s2.feather"))

      s1_s2_Okay <- ecokit::check_data(Path_s1, warning = FALSE) &&
        ecokit::check_data(Path_s2, warning = FALSE)

      if (s1_s2_Okay) {
        ecokit::cat_time(
          "s1 and s2 matrices were already saved",
          level = 2L, verbose = verbose)
      } else {

        ecokit::cat_time(
          "Saving s1 and s2 matrices", level = 2L, verbose = verbose)

        # s1
        s1 <- as.data.frame(LF_rL$s[units_model, , drop = FALSE])
        ecokit::save_as(object = s1, out_path = Path_s1)

        # s2
        s2 <- as.data.frame(LF_rL$s[units_pred[indNew], , drop = FALSE])
        ecokit::save_as(object = s2, out_path = Path_s2)

        rm(s1, s2, envir = environment())
      }

      rm(LF_rL, envir = environment())
      Path_D11 <- Path_D12 <- NULL

    } else {

      # Save D11 and D12 as feather files, if not already exist on disk
      Path_D11 <- fs::path(temp_dir, paste0(model_name, "D11.qs2"))
      Path_D12 <- fs::path(temp_dir, paste0(model_name, "D12.qs2"))

      if (file.exists(Path_D11) && file.exists(Path_D12)) {

        ecokit::cat_time(
          "D11 and D12 distance matrices are already saved",
          level = 2L, verbose = verbose)

      } else {

        s1 <- LF_rL$s[units_model, , drop = FALSE]
        s2 <- LF_rL$s[units_pred[indNew], , drop = FALSE]

        # D11
        D11 <- Rfast::Dist(s1)
        ecokit::save_as(object = D11, out_path = Path_D11)

        # D12
        D12 <- Rfast::dista(s1, s2)
        ecokit::save_as(object = D12, out_path = Path_D12)

        # Clean up
        rm(LF_rL, s1, s2, D11, D12, envir = environment())

      }

      Path_s1 <- Path_s2 <- NULL

    }

    invisible(gc())

    # # .................................................................... ###

    # Convert post_alpha to tibble

    ecokit::cat_time(
      "Splitting and saving `post_alpha` to small chunks",
      level = 1L, verbose = verbose)

    # Unique combination of LF / alphapw / sample IDs
    LF_Data <- do.call(rbind, post_alpha) %>%
      as.data.frame() %>%
      tibble::tibble() %>%
      stats::setNames(paste0("LF_", seq_len(ncol(.)))) %>%
      # ID column represents the original row number
      dplyr::mutate(Sample_IDs = dplyr::row_number()) %>%
      tidyr::nest(Sample_IDs = Sample_IDs) %>%
      dplyr::mutate(
        Sample_IDs = purrr::map(Sample_IDs, ~ as.vector(unlist(.x)))) %>%
      tidyr::pivot_longer(
        cols = -Sample_IDs, values_to = "Alpha_ID", names_to = "LF") %>%
      dplyr::mutate(
        LF_ID = as.integer(stringr::str_remove(LF, "LF_")),
        Denom = purrr::map_dbl(Alpha_ID, ~ alphapw[.x, 1])) %>%
      tidyr::nest(Sample_IDs = -c("Denom", "LF_ID", "Alpha_ID", "LF")) %>%
      dplyr::mutate(
        Sample_IDs = purrr::map(Sample_IDs, ~unname(sort(unlist(.x)))),
        # number of samples to be processed in each file
        NSamples = purrr::map_int(Sample_IDs, length)) %>%
      dplyr::arrange(dplyr::desc(NSamples)) %>%
      dplyr::mutate(
        ChunkID = dplyr::row_number(),
        File_postEta = purrr::map_chr(
          .x = ChunkID,
          .f = ~ {
            ChunkID0 <- stringr::str_pad(
              .x, width = nchar(dplyr::n()), pad = 0)
            fs::path(
              temp_dir,
              paste0(model_name, "postEta_ch", ChunkID0, ".feather"))
          }),
        File_etaPred_TF = stringr::str_replace_all(
          File_postEta, "_postEta_ch", "_etaPred_ch"),
        File_etaPred = stringr::str_replace_all(
          File_etaPred_TF, ".feather", ".qs2"),
        Export = purrr::pmap(
          .l = list(Sample_IDs, LF_ID, File_postEta),
          .f = function(Sample_IDs, LF_ID, File_postEta) {

            # do not export file if already exists
            if (!file.exists(File_postEta)) {
              Out <- postEta[Sample_IDs] %>%
                purrr::map(~ .x[, LF_ID, drop = FALSE]) %>%
                simplify2array() %>%
                as.data.frame() %>%
                stats::setNames(paste0("Eta_", Sample_IDs))

              ecokit::save_as(object = Out, out_path = File_postEta)
            }

            return(NULL)

          }),
        Export = NULL)

    rm(postEta, post_alpha, envir = environment())
    invisible(gc())

    # # .................................................................... ###

    # Internal functions to predict latent factors

    etaPreds_F <- function(RowNum, units_pred, LF_check = FALSE) {

      # do not use scientific notation
      withr::local_options(scipen = 99)

      # Current denominator
      Denom <- LF_Data$Denom[[RowNum]]
      # ID for latent factor
      LF_ID <- LF_Data$LF_ID[[RowNum]]
      # ID for posterior sample
      SampleID <- LF_Data$Sample_IDs[[RowNum]]

      # File path for current data
      File_postEta <- LF_Data$File_postEta[[RowNum]]
      File_etaPred <- LF_Data$File_etaPred[[RowNum]]
      File_etaPred_TF <- LF_Data$File_etaPred_TF[[RowNum]]

      CalcPredLF <- !file.exists(File_etaPred)

      if (LF_check && isFALSE(CalcPredLF)) {
        Eta_NCols <- ncol(ecokit::load_as(File_etaPred))
        CalcPredLF <- (Eta_NCols != LF_Data$NSamples[[RowNum]])
      }

      if (CalcPredLF) {

        # If the denominator is positive, perform calculations; otherwise, set
        # `eta_indNew` to zero.

        if (Denom > 0) {

          if (use_TF) {

            # Use TensorFlow

            # Suppress TensorFlow warnings and disable optimizations
            Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

            if (file.exists(File_etaPred_TF)) {
              eta_indNew0 <- File_etaPred_TF
            } else {
              eta_indNew0 <- run_crossprod_solve(
                TF_environ = TF_environ, s1 = Path_s1, s2 = Path_s2,
                denom = Denom, postEta = File_postEta,
                path_out = File_etaPred_TF,
                TF_use_single = TF_use_single,
                LF_commands_only = LF_commands_only,
                solve_max_attempts = solve_max_attempts,
                solve_chunk_size = solve_chunk_size)
            }

            if (LF_commands_only) {
              return(eta_indNew0)
            }

            rm(eta_indNew0, envir = environment())

            etaPred <- tibble::tibble(
              SampleID = SampleID,
              LF = LF_ID,
              Path_Samp_LF = fs::path(
                Temp_Dir_LF,
                paste0(
                  model_name, "Samp_",
                  stringr::str_pad(SampleID, width = 4, pad = "0"),
                  "_LF", LF, ".qs2")),
              etaPred = as.list(ecokit::load_as(File_etaPred_TF))) %>%
              tidyr::nest(eta_DT = -Path_Samp_LF) %>%
              dplyr::mutate(
                eta_DT = purrr::map(
                  .x = eta_DT,
                  .f = ~ {
                    tidyr::unnest_longer(.x, "etaPred") %>%
                      dplyr::mutate(LF = NULL, units_pred = units_pred) %>%
                      stats::setNames(
                        c("SampleID", paste0("LF_", LF_ID), "units_pred"))
                  }),
                File_etaPred = File_etaPred,
                ChunkID = LF_Data$ChunkID[[RowNum]],
                SampleID = SampleID,
                Save = purrr::map2(
                  .x = eta_DT, .y = Path_Samp_LF,
                  .f = ~qs2::qs_save(.x, .y, nthreads = 5)),
                LF = LF_ID, Save = NULL, eta_DT = NULL)

            ecokit::save_as(object = etaPred, out_path = File_etaPred)

          } else {

            # Use R / CPP

            # Reading postEta from file
            postEta0 <- ecokit::load_as(File_postEta)

            # Read D11 and D12
            D11 <- ecokit::load_as(Path_D11)
            D12 <- ecokit::load_as(Path_D12)

            K11 <- IASDT.R::exp_neg_div(D11, Denom)
            K12 <- IASDT.R::exp_neg_div(D12, Denom)

            etaPred <- purrr::map_chr(
              .x = seq_along(SampleID),
              .f = function(ID) {

                Path_Samp_LF <- fs::path(
                  Temp_Dir_LF,
                  paste0(
                    model_name, "Samp_",
                    stringr::str_pad(ID, width = 4, pad = "0"),
                    "_LF", LF_ID, ".qs2"))

                DT <- as.matrix(postEta0[, ID]) %>%
                  IASDT.R::solve2vect(K11, .) %>%
                  as.vector() %>%
                  Matrix::crossprod(K12, .) %>%
                  as.vector() %>%
                  tibble::tibble(
                    SampleID = SampleID[ID],
                    etaPred = ., units_pred = units_pred) %>%
                  stats::setNames(
                    c("SampleID", paste0("LF_", LF_ID), "units_pred"))

                qs2::qs_save(object = DT, file = Path_Samp_LF, nthreads = 5)

                return(Path_Samp_LF)
              }) %>%
              tibble::tibble(
                Path_Samp_LF = ., File_etaPred = File_etaPred,
                ChunkID = LF_Data$ChunkID[[RowNum]], SampleID = SampleID,
                LF = LF_ID)

            ecokit::save_as(object = etaPred, out_path = File_etaPred)
          }

        } else {

          # When Denom is zero, set `eta_indNew` to zero

          if (isFALSE(LF_commands_only)) {   # nolint: unneeded_nesting_linter

            etaPred <- tibble::tibble(

              Path_Samp_LF = fs::path(
                Temp_Dir_LF,
                paste0(
                  model_name, "Samp_",
                  stringr::str_pad(SampleID, width = 4L, pad = "0"),
                  "_LF", LF_ID, ".qs2")),

              File_etaPred = File_etaPred,
              ChunkID = LF_Data$ChunkID[[RowNum]],
              SampleID = SampleID) %>%
              dplyr::mutate(
                Save = purrr::map2(
                  .x = Path_Samp_LF, .y = SampleID,
                  .f = ~ {
                    DT <- tibble::tibble(
                      SampleID = .y, LF = 0,
                      units_pred = units_pred) %>%
                      stats::setNames(
                        c("SampleID", paste0("LF_", LF_ID), "units_pred"))
                    qs2::qs_save(object = DT, file = .x, nthreads = 5)
                    return(NULL)
                  }),
                Save = NULL, LF = LF_ID)

            ecokit::save_as(object = etaPred, out_path = File_etaPred)
          }
        }

      } else {
        etaPred <- ecokit::load_as(File_etaPred)
      }

      if (isFALSE(LF_commands_only)) {
        return(etaPred)
      }
    }

    invisible(gc())

    # # .................................................................... ###
    # # .................................................................... ###

    # Predict latent factors

    if (all(file.exists(LF_Data$File_etaPred))) {
      ecokit::cat_time(
        "All LF prediction files were already created",
        level = 1L, verbose = verbose)
    } else {
      if (LF_n_cores == 1 || LF_commands_only) {

        if (LF_commands_only) {
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
        etaPreds <- purrr::map(
          .x = seq_len(nrow(LF_Data)),
          .f = function(x) {
            result <- try(
              expr = etaPreds_F(
                RowNum = x, units_pred = units_pred, LF_check = LF_check),
              silent = FALSE)

            if (inherits(result, "try-error")) {
              return(NULL)
            } else {
              return(result)
            }

          })

        if (LF_commands_only) {

          CommandFilePrefix <- dplyr::if_else(
            startsWith(model_name, "RC_c_"),
            "LF_RC_Commands_", "LF_NewSites_Commands_")

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
                temp_dir, paste0(CommandFilePrefix, i, ".txt"))

              # Write the chunk to a file with Linux line endings
              writeLines(chunk, file_name, useBytes = TRUE)
            }
          }

          # Call the function to save commands to file
          save_commands_to_file(unlist(etaPreds))

          return(NULL)
        }

      } else {

        # Parallel processing
        ecokit::cat_time(
          "Predicting Latent Factor in parallel", level = 1L, verbose = verbose)

        ecokit::set_parallel(
          n_cores = min(LF_n_cores, nrow(LF_Data)), level = 2L,
          future_max_size = 800L, strategy = strategy)
        withr::defer(future::plan("sequential", gc = TRUE))

        ecokit::cat_time(
          "Making predictions in parallel", level = 2L, verbose = verbose)
        etaPreds <- future.apply::future_lapply(
          X = seq_len(nrow(LF_Data)),
          FUN = function(x) {
            result <- try(
              etaPreds_F(
                RowNum = x, units_pred = units_pred, LF_check = LF_check),
              silent = FALSE)

            invisible(gc())
            if (inherits(result, "try-error")) {
              return(NULL)
            } else {
              return(result)
            }
          },
          future.seed = TRUE, future.packages = pkg_to_export,
          future.globals = c(
            "LF_Data", "Path_D11", "Path_D12", "Path_s1", "Path_s2", "indNew",
            "units_pred", "indOld", "units_model", "TF_environ", "use_TF",
            "TF_use_single", "etaPreds_F", "LF_check", "run_crossprod_solve",
            "LF_commands_only", "solve_max_attempts", "solve_chunk_size",
            "Temp_Dir_LF"))

        # Stop the cluster
        ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
      }

      # Check if all files are created
      ecokit::cat_time(
        "Check if all files are created", level = 1L, verbose = verbose)
      AllEtaFiles <- LF_Data$File_etaPred
      AllEtaFilesExist <- all(file.exists(AllEtaFiles))

      if (!AllEtaFilesExist) {
        FailedFiles <- AllEtaFiles[!file.exists(AllEtaFiles)]
        ecokit::stop_ctx(
          paste0(length(FailedFiles), " files are missing"),
          FailedFiles = basename(FailedFiles), include_backtrace = TRUE)
      }
      ecokit::cat_time("All files were created", level = 2L, verbose = verbose)

    }

    invisible(gc())

    # # .................................................................... ###

    # Merge results
    ecokit::cat_time("Merge results in parallel", level = 1L, verbose = verbose)

    postEtaPred_Samp <- etaPreds %>%
      dplyr::bind_rows() %>%
      dplyr::select(-LF, -ChunkID, -File_etaPred) %>%
      dplyr::arrange(SampleID) %>%
      dplyr::mutate(
        Path_Sample = fs::path(
          Temp_Dir_LF,
          paste0(
            model_name, "Samp_",
            stringr::str_pad(SampleID, width = 4, pad = "0"), ".qs2"))) %>%
      tidyr::nest(data = -c("SampleID", "Path_Sample"))


    ecokit::set_parallel(
      n_cores = LF_n_cores, level = 2L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))

    ecokit::cat_time(
      "Process results for MCMC samples in parallel",
      level = 2L, verbose = verbose)

    postEtaPred <- future.apply::future_lapply(
      X = seq_len(nrow(postEtaPred_Samp)),
      FUN = function(x) {

        Path_Sample <- postEtaPred_Samp$Path_Sample[[x]]
        Path_LF <- postEtaPred_Samp$data[[x]]$Path_Samp_LF

        SampleDT0 <- lapply(sort(Path_LF), qs2::qs_read) %>%
          purrr::reduce(
            .f = dplyr::left_join, by = c("SampleID", "units_pred")) %>%
          dplyr::arrange(units_pred) %>%
          dplyr::select(sort(tidyselect::peek_vars())) %>%
          dplyr::select(-SampleID) %>%
          magrittr::set_rownames(NULL) %>%
          tibble::column_to_rownames("units_pred") %>%
          unname() %>%
          as.matrix()

        ecokit::save_as(SampleDT0, out_path = Path_Sample)
        try(fs::file_delete(Path_LF), silent = TRUE)
        return(SampleDT0)
      },
      future.seed = TRUE, future.packages = pkg_to_export,
      future.globals = c("postEtaPred_Samp", "LF_return"))

    # Stop the cluster
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)

    invisible(gc())

  }

  # # ..................................................................... ###

  # Save postEtaPred

  if (!is.null(LF_out_file)) {
    ecokit::cat_time(
      "Saving postEtaPred to disk", level = 1L, verbose = verbose)
    ecokit::cat_time(
      LF_out_file, cat_timestamp = FALSE, level = 2L, verbose = verbose)
    fs::dir_create(fs::path_dir(LF_out_file))
    ecokit::save_as(
      object = postEtaPred, out_path = LF_out_file,
      n_threads = 5L, compress_level = 6L)

    LF_OutFile_Samp <- stringr::str_replace(LF_out_file, ".qs2", "_Samp.qs2")
    ecokit::save_as(object = postEtaPred_Samp, out_path = LF_OutFile_Samp)
  }

  # # ..................................................................... ###

  # Clean up temporary files after finishing calculations
  if (LF_temp_cleanup) {

    ecokit::cat_time(
      "Cleaning up temporary files", level = 1L, verbose = verbose)

    try(
      expr = {

        Pattern <- paste0(
          "^", model_name,
          "(postEta|r[0-9]|etaPred|s1|s2|post).+(feather|qs2|log)$")
        file_paths <- list.files(
          path = ecokit::normalize_path(temp_dir),
          pattern = Pattern, full.names = TRUE)
        if (length(file_paths) > 0) {
          try(fs::file_delete(file_paths), silent = TRUE)
        }

        file_paths2 <- list.files(
          path = ecokit::normalize_path(temp_dir),
          pattern = "(LF_.+_Test|RC_c)_Samp_.+.qs2",
          full.names = TRUE, recursive = TRUE)
        if (length(file_paths2) > 0) {
          try(fs::file_delete(file_paths2), silent = TRUE)
        }

        # delete temp files for cross-validated models
        file_paths3 <- list.files(
          path = ecokit::normalize_path(Temp_Dir_LF),
          pattern = paste0("^", model_name, "Samp_.+.qs2"),
          full.names = TRUE, recursive = TRUE)

        if (length(file_paths3) > 0) {
          try(fs::file_delete(file_paths3), silent = TRUE)
        }

      },
      silent = TRUE)

  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "predict_latent_factor was finished in ",
    level = 1L, verbose = verbose)

  # # ..................................................................... ###

  if (LF_return) {
    return(postEtaPred)
  } else {
    return(LF_out_file)
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
#' computations using TensorFlow with provided inputs. Retries up to three times
#' if the output file validation fails.
#'
#' @param s1 Character. Path to the input file containing s1 coordinates.
#' @param s2 Character Path to the input file containing s2 coordinates.
#' @param postEta Character. Path to the file containing the `postEta` matrix
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
    TF_environ, s1, s2, postEta, path_out, denom,
    chunk_size = 1000L, threshold_mb = 2000L, TF_use_single = TRUE,
    verbose = TRUE, solve_chunk_size = 50L, solve_max_attempts = 5L,
    LF_commands_only = FALSE) {

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")

  # do not use scientific notation
  withr::local_options(scipen = 99)

  script_path <- system.file("crossprod_solve.py", package = "IASDT.R")
  if (!file.exists(script_path)) {
    ecokit::stop_ctx(
      "Necessary Python script `crossprod_solve.py` does not exist",
      script_path = script_path, include_backtrace = TRUE)
  }

  # Ensure the paths are valid
  paths <- list(script_path, s1, s2, postEta)
  names(paths) <- c("Python Script", "s1", "s2", "postEta")
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
  # `TF_environ`. On LUMI, this is not needed as the compatible Python
  # installation is loaded automatically when loading tensorflow module. When
  # using another HPC system, the function needs to be adapted accordingly.

  if (.Platform$OS.type == "windows") {

    if (is.null(TF_environ)) {
      ecokit::stop_ctx(
        "When running on Windows, `TF_environ` must be specified",
        TF_environ = TF_environ, include_backtrace = TRUE)
    }
    if (!dir.exists(TF_environ)) {
      ecokit::stop_ctx(
        "The specified `TF_environ` directory does not exist",
        TF_environ = ecokit::normalize_path(TF_environ),
        include_backtrace = TRUE)
    }

    python_executable <- ecokit::normalize_path(
      fs::path(TF_environ, "Scripts", "python.exe"), must_work = TRUE)

    if (!file.exists(python_executable)) {
      ecokit::stop_ctx(
        "Python executable not found in the virtual environment.",
        python_executable = python_executable, include_backtrace = TRUE)
    }

  } else {
    # Use `python3`` directly - on LUMI, compatible Python installation is
    # loaded automatically when loading tensorflow
    python_executable <- "/usr/bin/time -v python3"
  }

  # Construct the command to run the Python script
  LF_Args <- c(
    python_executable,
    script_path,
    "--s1", ecokit::normalize_path(s1, must_work = TRUE),
    "--s2", ecokit::normalize_path(s2, must_work = TRUE),
    "--post_eta", ecokit::normalize_path(postEta, must_work = TRUE),
    "--path_out", ecokit::normalize_path(path_out),
    "--denom", as.character(denom),
    "--chunk_size", as.character(chunk_size),
    "--threshold_mb", as.character(threshold_mb),
    "--solve_chunk_size", as.character(solve_chunk_size))

  # Add boolean flags conditionally
  if (TF_use_single) {
    LF_Args <- c(LF_Args, "--use_single")
  }

  if (verbose) {
    LF_Args <- c(LF_Args, "--verbose")
  }

  if (.Platform$OS.type != "windows") {
    path_log <- stringr::str_replace(
      fs::path(getwd(), path_out), ".feather", ".log")
    # Redirect results of time to log file
    LF_Args <- c(LF_Args, paste0(" >> ", path_log, " 2>&1"))
  }

  LF_Args <- paste(LF_Args, collapse = " ")

  if (LF_commands_only) {
    return(LF_Args)
  } else {

    path_log <- stringr::str_replace(path_out, ".feather", ".log")
    f <- file(path_log, open = "a")
    on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
    cat(
      "Running command:\n", paste(LF_Args, "\n\n"),
      sep = "\n", file = f, append = TRUE)

    # Initialise retry logic
    attempt <- 1
    success <- FALSE

    while (attempt <= solve_max_attempts && !success) {

      cat(
        paste0("Attempt ", attempt, " of ", solve_max_attempts, "\n\n"),
        file = f, append = TRUE)

      # Run the command and capture stdout/stderr to a log file
      result <- system(LF_Args, intern = TRUE)

      # Check for errors
      if (!inherits(result, "error") || length(result) != 0 ||
          result == "Done") {
        # Check the integrity of the file
        success <- ecokit::check_data(path_out)
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
#' @param path_model Path to the model file.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param plot_width,plot_height Numeric. The width and height of the
#'   output plot in cm. Default is 20&times;21 cm.
#' @export
#' @author Ahmed El-Gabbas
#' @name plot_latent_factor

plot_latent_factor <- function(
    path_model = NULL, env_file = ".env", plot_width = 20, plot_height = 21) {

  # # ..................................................................... ###

  Path_Grid <- NULL

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Grid10 <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData") %>%
    ecokit::load_as() %>%
    terra::unwrap()

  # # ..................................................................... ###

  # Check if the model file exists
  if (is.null(path_model) || !file.exists(path_model)) {
    ecokit::stop_ctx(
      "Selected model files not found", path_model = path_model,
      include_backtrace = TRUE)
  }

  Model <- ecokit::load_as(path_model)
  Model_Coords <- Model$ranLevels$sample$s
  postEta <- Hmsc::getPostEstimate(Model, parName = "Eta")
  N_LF <- ncol(postEta$mean)
  Eta_Mean <- as.data.frame(postEta$mean) %>%
    stats::setNames(paste("LF", seq_len(N_LF), sep = "_")) %>%
    cbind.data.frame(Model_Coords, .) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)
  Eta_Mean_R <- terra::rasterize(
    Eta_Mean, Grid10, field = names(Eta_Mean)[-ncol(Eta_Mean)],
    fun = "mean") %>%
    stats::setNames(stringr::str_remove(names(.), "_mean")) %>%
    stats::setNames(stringr::str_replace(names(.), "LF_", "Latent factor "))
  rm(Model, Model_Coords, Eta_Mean, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  Xlim <- c(2600000, 6700000)
  Ylim <- c(1450000, 5420000)

  LF_Plot <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = Eta_Mean_R, maxcell = Inf) +
    ggplot2::facet_wrap(~lyr, ncol = 2) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = ecokit::integer_breaks(), name = NULL) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
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

  ragg::agg_jpeg(
    filename = fs::path(
      dirname(dirname(path_model)), "Model_Prediction", "LF_Plot.jpeg"),
    width = plot_width, height = plot_height, res = 600,
    quality = 100, units = "cm")
  print(LF_Plot)
  grDevices::dev.off()

  return(invisible(NULL))

}
