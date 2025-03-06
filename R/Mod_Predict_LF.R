## |------------------------------------------------------------------------| #
# Mod_Predict_LF ----
## |------------------------------------------------------------------------| #

#' Draws samples from the conditional predictive distribution of latent factors
#'
#' This function is optimized for speed using parallel processing and optionally
#' TensorFlow for matrix operations. This function is adapted from
#' [Hmsc::predictLatentFactor] with equivalent results to the original function
#' when `predictMean = TRUE`.
#' @param Units_pred a factor vector with random level units for which
#'   predictions are to be made
#' @param Units_model a factor vector with random level units that are
#'   conditioned on
#' @param postEta Character. Path of `postEta`; a list containing samples of
#'   random factors at conditioned units
#' @param postAlpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param LF_rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param LF_NCores Integer. Number of cores to use for parallel processing of
#'   latent factor prediction. Defaults to 8L.
#' @param Temp_Dir Character. Path for temporary storage of intermediate files.
#' @param LF_Temp_Cleanup Logical. Whether to delete temporary files in the
#'   `Temp_Dir` directory after finishing the LF predictions.
#' @param Model_Name Character. Prefix for temporary file names. Defaults to
#'   `NULL`, in which case no prefix is used.
#' @param UseTF Logical. Whether to use TensorFlow for calculations. Defaults to
#'   `TRUE`.
#' @param TF_Environ Character. Path to the Python environment. This argument is
#'   required if `UseTF` is `TRUE`.
#' @param LF_Commands_Only Logical. If `TRUE`, returns the command to run the
#'   Python script. Default is `FALSE`.
#' @param TF_use_single Logical. Whether to use single precision for the
#'   TensorFlow calculations. Defaults to `FALSE`.
#' @param LF_OutFile Character. Path to save the outputs. If `NULL` (default),
#'   the predicted latent factors are not saved to a file. This should end with
#'   either `*.qs2` or `*.RData`.
#' @param LF_Return Logical. Whether the output should be returned. Defaults to
#'   `FALSE`. If `LF_OutFile` is `NULL`, this parameter cannot be set to `FALSE`
#'   because the function needs to return the result if it is not saved to a
#'   file.
#' @param LF_Check Logical. If `TRUE`, the function checks if the output files
#'   are already created and valid. If `FALSE`, the function will only check if
#'   the files exist without checking their integrity. Default is `FALSE`.
#' @param Verbose Logical. If `TRUE`, logs detailed information during
#'   execution. Default is `TRUE`.
#' @param solve_max_attempts Integer. Maximum number of attempts to run solve
#'   and crossprod internal function [run_crossprod_solve]. Default is 5L.
#' @param solve_chunk_size Integer. Chunk size for `solve_and_multiply` Python
#'   function. Default is 50L.
#' @export
#' @seealso [Hmsc::predictLatentFactor]
#' @name Mod_Predict_LF
#' @details The function is expected to be faster than the original function in
#'   the `Hmsc` package, especially when using `TensorFlow` for calculations and
#'   when working on parallel.
#'
#'   The main difference is that this function:
#' - allow for parallel processing (`LF_NCores` argument);
#' - when TensorFlow is used (`UseTF = TRUE`), matrix
#'   calculations are much faster, particularly when used on GPU. The following
#'   Python modules are needed: `numpy`, `tensorflow`, `rdata`, `xarray`, and
#'   `pandas`. To use `TensorFlow`, the argument `TF_Environ` should be set to
#'   the path of a Python environment with TensorFlow installed;
#' - if `UseTF` is set to `FALSE`, the function uses `R` (supported by
#'   relatively faster `CPP` functions) in the calculations;
#' - `D11` and `D12` matrices are processed only once and saved to disk and
#'   called when needed.

Mod_Predict_LF <- function(
    Units_pred, Units_model, postEta, postAlpha, LF_rL, LF_NCores = 8L,
    Temp_Dir = "TEMP_Pred", LF_Temp_Cleanup = TRUE, Model_Name = NULL,
    UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE, LF_OutFile = NULL,
    LF_Return = FALSE, LF_Check = FALSE, LF_Commands_Only = FALSE,
    solve_max_attempts = 5L, solve_chunk_size = 50L, Verbose = TRUE) {

  # # ..................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  IASDT.R::CatTime("Starting `Mod_Predict_LF` function", Level = 1)

  # # ..................................................................... ###

  # Check inputs

  if (is.null(LF_OutFile) && isFALSE(LF_Return)) {
    stop(
      "`LF_Return` must be TRUE when `LF_OutFile` is NULL.",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Load postEta if it is a file path

  if (inherits(postEta, "character")) {
    IASDT.R::CatTime("Load postEta", Level = 1)
    if (!file.exists(postEta)) {
      stop("The specified path for `postEta` does not exist. ", call. = FALSE)
    }
    postEta <- IASDT.R::LoadAs(postEta)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SampleID <- LF <- LF_ID <- Sample_IDs <- Alpha_ID <- File_etaPred <-
    NSamples <- ChunkID <- File_postEta <- File_etaPred_TF <- eta_DT <-
    Path_Samp_LF <- NULL

  # # ..................................................................... ###

  # indices of Units_pred in Units_model
  indOld <- (Units_pred %in% Units_model)
  # indices of new Units_pred
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
    stop(
      "The input sites should be either all training sites or all new sites.",
      call. = FALSE)
  }

  if (AllTraining) {
    # If all input sites are for training sites, use LF info from the model
    # directly
    IASDT.R::CatTime("All input sites are training sites", Level = 1)

    postEtaPred <- purrr::map(
      .x = postEta,
      .f = ~ {
        Out <- .x[match(Units_pred[indOld], Units_model), ]
        rownames(Out) <- Units_model
        return(Out)
      })

  } else {

    IASDT.R::CatTime("All input sites are new sites", Level = 1)

    # Check TensorFlow settings

    if (UseTF) {

      PythonScript <- system.file("crossprod_solve.py", package = "IASDT.R")

      # Check if PythonScript exists
      if (!file.exists(PythonScript)) {
        stop("Necessary Python script does not exist", call. = FALSE)
      }

      # Suppress TensorFlow warnings and disable optimizations
      Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

      IASDT.R::CatTime("Computations will be made using TensorFlow", Level = 1)
    } else {
      IASDT.R::CatTime("Computations will be made using R/CPP", Level = 1)
    }

    # # .................................................................... ###

    # Adjust Model_Name prefix

    if (is.null(Model_Name)) {
      Model_Name <- ""
    } else {
      Model_Name <- paste0(Model_Name, "_")
    }

    # Create a temporary directory to store intermediate results. This directory
    # will be used to save s1/s2 or D11/D12, and intermediate postEta files,
    # reducing memory usage.
    Temp_Dir_LF <- IASDT.R::Path(Temp_Dir, "LF_Prediction")
    fs::dir_create(c(Temp_Dir_LF, Temp_Dir))

    # # .................................................................... ###

    # Calculate D11 and D12 only once

    IASDT.R::CatTime("Calculate/save necessary matrices", Level = 1)

    alphapw <- LF_rL$alphapw

    if (UseTF) {

      # Save s1 and s2 for coordinates at training and testing sites as feather
      # files, if not already exist on disk
      Path_s1 <- IASDT.R::Path(Temp_Dir, paste0(Model_Name, "s1.feather"))
      Path_s2 <- IASDT.R::Path(Temp_Dir, paste0(Model_Name, "s2.feather"))

      s1_s2_Okay <- IASDT.R::CheckData(Path_s1, warning = FALSE) &&
        IASDT.R::CheckData(Path_s2, warning = FALSE)

      if (s1_s2_Okay) {
        IASDT.R::CatTime("s1 and s2 matrices were already saved", Level = 2)
      } else {

        IASDT.R::CatTime("Saving s1 and s2 matrices", Level = 2)

        # s1
        s1 <- as.data.frame(LF_rL$s[Units_model, , drop = FALSE])
        IASDT.R::SaveAs(InObj = s1, OutPath = Path_s1)

        # s2
        s2 <- as.data.frame(LF_rL$s[Units_pred[indNew], , drop = FALSE])
        IASDT.R::SaveAs(InObj = s2, OutPath = Path_s2)

        rm(s1, s2, envir = environment())
      }

      rm(LF_rL, envir = environment())
      Path_D11 <- Path_D12 <- NULL

    } else {

      # Save D11 and D12 as feather files, if not already exist on disk
      Path_D11 <- IASDT.R::Path(Temp_Dir, paste0(Model_Name, "D11.qs2"))
      Path_D12 <- IASDT.R::Path(Temp_Dir, paste0(Model_Name, "D12.qs2"))

      if (file.exists(Path_D11) && file.exists(Path_D12)) {

        IASDT.R::CatTime(
          "D11 and D12 distance matrices are already saved", Level = 2)

      } else {

        s1 <- LF_rL$s[Units_model, , drop = FALSE]
        s2 <- LF_rL$s[Units_pred[indNew], , drop = FALSE]

        # D11
        D11 <- Rfast::Dist(s1)
        IASDT.R::SaveAs(InObj = D11, OutPath = Path_D11)

        # D12
        D12 <- Rfast::dista(s1, s2)
        IASDT.R::SaveAs(InObj = D12, OutPath = Path_D12)

        # Clean up
        rm(LF_rL, s1, s2, D11, D12, envir = environment())

      }

      Path_s1 <- Path_s2 <- NULL

    }

    invisible(gc())

    # # .................................................................... ###

    # Convert postAlpha to tibble

    IASDT.R::CatTime(
      "Splitting and saving `postAlpha` to small chunks", Level = 1)

    # Unique combination of LF / alphapw / sample IDs
    LF_Data <- do.call(rbind, postAlpha) %>%
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
            IASDT.R::Path(
              Temp_Dir,
              paste0(Model_Name, "postEta_ch", ChunkID0, ".feather"))
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

              IASDT.R::SaveAs(InObj = Out, OutPath = File_postEta)
            }

            return(NULL)

          }),
        Export = NULL)

    rm(postEta, postAlpha, envir = environment())
    invisible(gc())

    # # .................................................................... ###

    # Internal functions to predict latent factors

    etaPreds_F <- function(RowNum, LF_Check = FALSE, Units_pred) {

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

      if (LF_Check && isFALSE(CalcPredLF)) {
        Eta_NCols <- ncol(IASDT.R::LoadAs(File_etaPred))
        CalcPredLF <- (Eta_NCols != LF_Data$NSamples[[RowNum]])
      }

      if (CalcPredLF) {

        # If the denominator is positive, perform calculations; otherwise, set
        # `eta_indNew` to zero.

        if (Denom > 0) {

          if (UseTF) {

            # Use TensorFlow

            # Suppress TensorFlow warnings and disable optimizations
            Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

            if (file.exists(File_etaPred_TF)) {
              eta_indNew0 <- File_etaPred_TF
            } else {
              eta_indNew0 <- run_crossprod_solve(
                TF_Environ = TF_Environ, s1 = Path_s1, s2 = Path_s2,
                denom = Denom, postEta = File_postEta,
                path_out = File_etaPred_TF,
                TF_use_single = TF_use_single,
                LF_Commands_Only = LF_Commands_Only,
                solve_max_attempts = solve_max_attempts,
                solve_chunk_size = solve_chunk_size)
            }

            if (LF_Commands_Only) {
              return(eta_indNew0)
            }

            rm(eta_indNew0, envir = environment())

            etaPred <- tibble::tibble(
              SampleID = SampleID,
              LF = LF_ID,
              Path_Samp_LF = IASDT.R::Path(
                Temp_Dir_LF,
                paste0(
                  Model_Name, "Samp_",
                  stringr::str_pad(SampleID, width = 4, pad = "0"),
                  "_LF", LF, ".qs2")),
              etaPred = as.list(IASDT.R::LoadAs(File_etaPred_TF))) %>%
              tidyr::nest(eta_DT = -Path_Samp_LF) %>%
              dplyr::mutate(
                eta_DT = purrr::map(
                  .x = eta_DT,
                  .f = ~ {
                    tidyr::unnest_longer(.x, "etaPred") %>%
                      dplyr::mutate(LF = NULL, Units_pred = Units_pred) %>%
                      stats::setNames(
                        c("SampleID", paste0("LF_", LF_ID), "Units_pred"))
                  }),
                File_etaPred = File_etaPred,
                ChunkID = LF_Data$ChunkID[[RowNum]],
                SampleID = SampleID,
                Save = purrr::map2(
                  .x = eta_DT, .y = Path_Samp_LF,
                  .f = ~qs2::qs_save(.x, .y, nthreads = 5)),
                LF = LF_ID, Save = NULL, eta_DT = NULL)

            IASDT.R::SaveAs(InObj = etaPred, OutPath = File_etaPred)

          } else {

            # Use R / CPP

            # Reading postEta from file
            postEta0 <- IASDT.R::LoadAs(File_postEta)

            # Read D11 and D12
            D11 <- IASDT.R::LoadAs(Path_D11)
            D12 <- IASDT.R::LoadAs(Path_D12)

            K11 <- IASDT.R::exp_neg_div(D11, Denom)
            K12 <- IASDT.R::exp_neg_div(D12, Denom)

            etaPred <- purrr::map_chr(
              .x = seq_along(SampleID),
              .f = function(ID) {

                Path_Samp_LF <- IASDT.R::Path(
                  Temp_Dir_LF,
                  paste0(
                    Model_Name, "Samp_",
                    stringr::str_pad(ID, width = 4, pad = "0"),
                    "_LF", LF_ID, ".qs2"))

                DT <- as.matrix(postEta0[, ID]) %>%
                  IASDT.R::Solve2vect(K11, .) %>%
                  as.vector() %>%
                  Matrix::crossprod(K12, .) %>%
                  as.vector() %>%
                  tibble::tibble(
                    SampleID = SampleID[ID],
                    etaPred = ., Units_pred = Units_pred) %>%
                  stats::setNames(
                    c("SampleID", paste0("LF_", LF_ID), "Units_pred"))

                qs2::qs_save(object = DT, file = Path_Samp_LF, nthreads = 5)
                return(Path_Samp_LF)
              }) %>%
              tibble::tibble(
                Path_Samp_LF = ., File_etaPred = File_etaPred,
                ChunkID = LF_Data$ChunkID[[RowNum]], SampleID = SampleID,
                LF = LF_ID)

            IASDT.R::SaveAs(InObj = etaPred, OutPath = File_etaPred)
          }

        } else {

          # When Denom is zero, set `eta_indNew` to zero

          if (isFALSE(LF_Commands_Only)) {

            etaPred <- tibble::tibble(

              Path_Samp_LF = IASDT.R::Path(
                Temp_Dir_LF,
                paste0(
                  Model_Name, "Samp_",
                  stringr::str_pad(SampleID, width = 4, pad = "0"),
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
                      Units_pred = Units_pred) %>%
                      stats::setNames(
                        c("SampleID", paste0("LF_", LF_ID), "Units_pred"))
                    qs2::qs_save(object = DT, file = .x, nthreads = 5)
                    return(NULL)
                  }),
                Save = NULL, LF = LF_ID)

            IASDT.R::SaveAs(InObj = etaPred, OutPath = File_etaPred)
          }
        }

      } else {
        etaPred <- IASDT.R::LoadAs(File_etaPred)
      }

      if (isFALSE(LF_Commands_Only)) {
        return(etaPred)
      }
    }

    invisible(gc())

    # # .................................................................... ###
    # # .................................................................... ###

    # Predict latent factors

    if (all(file.exists(LF_Data$File_etaPred))) {
      IASDT.R::CatTime(
        "All LF prediction files were already created", Level = 1)
    } else {
      if (LF_NCores == 1 || LF_Commands_Only) {

        if (LF_Commands_Only) {
          IASDT.R::CatTime(
            "Prepare commands for predicting latent factors", Level = 1)
        } else {
          # Sequential processing
          IASDT.R::CatTime("Predicting Latent Factor sequentially", Level = 1)
        }

        # Making predictions sequentially
        etaPreds <- purrr::map(
          .x = seq_len(nrow(LF_Data)),
          .f = function(x) {

            result <- try(
              expr = etaPreds_F(
                RowNum = x, LF_Check = LF_Check, Units_pred = Units_pred),
              silent = FALSE)

            if (inherits(result, "try-error")) {
              return(NULL)
            } else {
              return(result)
            }

          })

        if (LF_Commands_Only) {

          CommandFilePrefix <- dplyr::if_else(
            stringr::str_detect(Model_Name, "^RC_c_"),
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
              file_name <- IASDT.R::Path(
                Temp_Dir, paste0(CommandFilePrefix, i, ".txt"))
              
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
        IASDT.R::CatTime(
          paste0(
            "Predicting Latent Factor in parallel using ",
            min(LF_NCores, nrow(LF_Data)), " cores"),
          Level = 1)

        IASDT.R::CatTime("Prepare for parallel processing", Level = 2)
        withr::local_options(
          future.globals.maxSize = 8000 * 1024^2, cluster.timeout = 10 * 60,
          future.gc = TRUE, future.seed = TRUE)
        c1 <- parallel::makeCluster(min(LF_NCores, nrow(LF_Data)))
        on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

        IASDT.R::CatTime("Export objects to cores", Level = 2)
        parallel::clusterExport(
          cl = c1,
          varlist = c(
            "LF_Data", "Path_D11", "Path_D12", "Path_s1", "Path_s2", "indNew",
            "Units_pred", "indOld", "Units_model", "TF_Environ", "UseTF",
            "TF_use_single", "etaPreds_F", "LF_Check", "run_crossprod_solve",
            "LF_Commands_Only", "solve_max_attempts", "solve_chunk_size",
            "Temp_Dir_LF"),
          envir = environment())

        # Load necessary libraries and load environment if using TensorFlow
        IASDT.R::CatTime("Load necessary libraries", Level = 2)
        invisible(parallel::clusterEvalQ(
          cl = c1,
          expr = {
            sapply(
              c(
                "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble", "arrow",
                "Matrix", "Hmsc", "qs2", "fs", "purrr", "IASDT.R"),
              library, character.only = TRUE)
            invisible(gc())
          }))

        # Making predictions on parallel
        IASDT.R::CatTime("Making predictions on parallel", Level = 2)

        etaPreds <- parallel::clusterApplyLB(
          cl = c1,
          x = seq_len(nrow(LF_Data)),
          fun = function(x) {
            result <- try(
              etaPreds_F(x, LF_Check = LF_Check, Units_pred = Units_pred),
              silent = FALSE)
            invisible(gc())
            if (inherits(result, "try-error")) {
              return(NULL)
            } else {
              return(result)
            }
          })

        # Stop the cluster
        IASDT.R::CatTime("Stop the cluster", Level = 2)
        parallel::stopCluster(c1)
        invisible(gc())
      }

      # Check if all files are created
      IASDT.R::CatTime("Check if all files are created", Level = 1)
      AllEtaFiles <- LF_Data$File_etaPred
      AllEtaFilesExist <- all(file.exists(AllEtaFiles))

      if (!AllEtaFilesExist) {
        FailedFiles <- AllEtaFiles[!file.exists(AllEtaFiles)]
        stop(
          length(FailedFiles), " files are missing: \n",
          paste0("  >>  ", basename(FailedFiles), collapse = "\n"),
          call. = FALSE)
      }
      IASDT.R::CatTime("All files were created", Level = 2)

    }

    # # .................................................................... ###

    # Merge results
    IASDT.R::CatTime("Merge results", Level = 1)

    postEtaPred_Samp <- etaPreds %>%
      dplyr::bind_rows() %>%
      dplyr::select(-LF, -ChunkID, -File_etaPred) %>%
      dplyr::arrange(SampleID) %>%
      dplyr::mutate(
        Path_Sample = IASDT.R::Path(
          Temp_Dir_LF,
          paste0(
            Model_Name, "Samp_",
            stringr::str_pad(SampleID, width = 4, pad = "0"), ".qs2"))) %>%
      tidyr::nest(data = -c("SampleID", "Path_Sample"))

    IASDT.R::CatTime("Prepare for parallel processing", Level = 2)
    c1 <- parallel::makeCluster(LF_NCores)
    on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

    parallel::clusterExport(
      cl = c1, varlist = c("postEtaPred_Samp", "LF_Return"),
      envir = environment())
    invisible(parallel::clusterEvalQ(
      cl = c1,
      expr = {
        sapply(
          c("dplyr", "magrittr", "tibble", "IASDT.R", "purrr"),
          library, character.only = TRUE)
        invisible(gc())
      }))

    # Merge results on parallel
    IASDT.R::CatTime("Process results for MCMC samples on parallel", Level = 2)

    postEtaPred <- parallel::parLapply(
      cl = c1,
      X = seq_len(nrow(postEtaPred_Samp)),
      fun = function(x) {

        Path_Sample <- postEtaPred_Samp$Path_Sample[[x]]
        Path_LF <- postEtaPred_Samp$data[[x]]$Path_Samp_LF

        SampleDT0 <- lapply(sort(Path_LF), qs2::qs_read) %>%
          purrr::reduce(
            .f = dplyr::left_join, by = c("SampleID", "Units_pred")) %>%
          dplyr::arrange(Units_pred) %>%
          dplyr::select(sort(tidyselect::peek_vars())) %>%
          dplyr::select(-SampleID) %>%
          magrittr::set_rownames(NULL) %>%
          tibble::column_to_rownames("Units_pred") %>%
          unname() %>%
          as.matrix()

        IASDT.R::SaveAs(SampleDT0, OutPath = Path_Sample)
        try(fs::file_delete(Path_LF), silent = TRUE)
        return(SampleDT0)
      })

    # Stop the cluster
    IASDT.R::CatTime("Stop the cluster", Level = 2)
    parallel::stopCluster(c1)
    invisible(gc())

  }

  # # ..................................................................... ###

  # Save postEtaPred

  if (!is.null(LF_OutFile)) {
    IASDT.R::CatTime("Saving postEtaPred to: ", Level = 1)
    IASDT.R::CatTime(paste0("`", LF_OutFile, "`"), Time = FALSE, Level = 2)
    fs::dir_create(fs::path_dir(LF_OutFile))
    IASDT.R::SaveAs(
      InObj = postEtaPred, OutPath = LF_OutFile,
      nthreads = 10, compress_level = 6)

    LF_OutFile_Samp <- stringr::str_replace(LF_OutFile, ".qs2", "_Samp.qs2")
    IASDT.R::SaveAs(InObj = postEtaPred_Samp, OutPath = LF_OutFile_Samp)
  }

  # # ..................................................................... ###

  # Clean up temporary files after finishing calculations
  if (LF_Temp_Cleanup) {

    IASDT.R::CatTime("Cleaning up temporary files", Level = 1)

    try(
      expr = {
        Pattern <- paste0(
          "^", Model_Name,
          "(postEta|r[0-9]|etaPred|s1|s2|post).+(feather|qs2|log)$")
        
        file_paths <- list.files(
          path = IASDT.R::NormalizePath(Temp_Dir),
          pattern = Pattern, full.names = TRUE)
        
        fs::file_delete(file_paths)
      },
      silent = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Mod_Predict_LF was finished in ",
    Level = 1)

  # # ..................................................................... ###

  if (LF_Return) {
    return(postEtaPred)
  } else {
    return(LF_OutFile)
  }

}


# # ========================================================================== #
# # ========================================================================== #


## |------------------------------------------------------------------------| #
# run_crossprod_solve ----
## |------------------------------------------------------------------------| #

#' run_crossprod_solve

#' Run Crossprod Solve
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
#' - Verifies the output file validity using `IASDT.R::CheckData`. Retries up
#' to 3 times if the output is invalid.
#' - Generates detailed logs if `verbose` is set to `TRUE`.
#' @author Ahmed El-Gabbas
#' @keywords internal
#' @inheritParams Mod_Predict_LF

run_crossprod_solve <- function(
    TF_Environ, s1, s2, postEta, path_out, denom,
    chunk_size = 1000L, threshold_mb = 2000L, TF_use_single = TRUE,
    verbose = TRUE,
    solve_chunk_size = 50L, solve_max_attempts = 5L, LF_Commands_Only = FALSE) {

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")

  # do not use scientific notation
  withr::local_options(scipen = 99)

  script_path <- system.file("crossprod_solve.py", package = "IASDT.R")
  if (!file.exists(script_path)) {
    stop(
      "Necessary Python script `crossprod_solve.py` does not exist",
      call. = FALSE)
  }

  # Ensure the paths are valid
  paths <- list(script_path, s1, s2, postEta)
  names(paths) <- c("Python Script", "s1", "s2", "postEta")
  purrr::walk(
    .x = names(paths),
    .f = function(p) {
      if (!file.exists(paths[[p]])) {
        stop(p, " does not exist: ", paths[[p]], call. = FALSE)
      }
    })

  # Determine the Python executable path

  # On Windows, the TF calculations has to be done through a valid virtual
  # environment; the path to the virtual environment must be specified in
  # `TF_Environ`. On LUMI, this is not needed as the compatible Python
  # installation is loaded automatically when loading tensorflow module. When
  # using another HPC system, the function needs to be adapted accordingly.

  if (.Platform$OS.type == "windows") {
    if (isFALSE(LF_Commands_Only)) {
      if (is.null(TF_Environ)) {
        stop(
          "When running on Windows, `TF_Environ` must be specified",
          call. = FALSE)
      }
      if (!dir.exists(TF_Environ)) {
        stop(
          "The specified `TF_Environ` directory ",
          IASDT.R::NormalizePath(TF_Environ), " does not exist", call. = FALSE)
      }
    }

    python_executable <- IASDT.R::NormalizePath(
      IASDT.R::Path(TF_Environ, "Scripts", "python.exe"),
      MustWork = TRUE)

    if (!file.exists(python_executable)) {
      stop(
        "Python executable not found in the virtual environment.",
        call. = FALSE)
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
    "--s1", IASDT.R::NormalizePath(s1, MustWork = TRUE),
    "--s2", IASDT.R::NormalizePath(s2, MustWork = TRUE),
    "--post_eta", IASDT.R::NormalizePath(postEta, MustWork = TRUE),
    "--path_out", IASDT.R::NormalizePath(path_out),
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
      IASDT.R::Path(getwd(), path_out), ".feather", ".log")
    # Redirect results of time to log file
    LF_Args <- c(LF_Args, paste0(" >> ", path_log, " 2>&1"))
  }

  LF_Args <- paste(LF_Args, collapse = " ")

  if (LF_Commands_Only) {
    return(LF_Args)
  } else {

    path_log <- stringr::str_replace(path_out, ".feather", ".log")
    f <- file(path_log, open = "a")
    on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
    cat(
      "Running command:\n", paste(LF_Args, "\n\n"),
      sep = "\n", file = f, append = TRUE)

    # Initialize retry logic
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
        success <- IASDT.R::CheckData(path_out)
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
# Mod_Plot_LF ----
## |------------------------------------------------------------------------| #

#' Plot spatial variation in site loadings of HMSC models
#'
#' Generate and save spatial variation in site loadings of HMSC models' latent
#' factors as a JPEG file.
#'
#' @param Path_Model Path to the model file.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Plot_Width,Plot_Height Numeric. The width and height of the
#'   output plot in cm. Default is 20&times;21 cm.
#' @export
#' @author Ahmed El-Gabbas
#' @name Mod_Plot_LF

Mod_Plot_LF <- function(
    Path_Model = NULL, EnvFile = ".env", Plot_Width = 20, Plot_Height = 21) {

  # # ..................................................................... ###

  Path_Grid <- NULL

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  Grid10 <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData") %>%
    IASDT.R::LoadAs() %>%
    terra::unwrap()

  # # ..................................................................... ###

  # Check if the model file exists
  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop("Selected model files not found", call. = FALSE)
  }

  Model <- IASDT.R::LoadAs(Path_Model)
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
      breaks = IASDT.R::integer_breaks(), name = NULL) +
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
    filename = IASDT.R::Path(
      dirname(dirname(Path_Model)), "Model_Prediction", "LF_Plot.jpeg"),
    width = Plot_Width, height = Plot_Height, res = 600,
    quality = 100, units = "cm")
  print(LF_Plot)
  grDevices::dev.off()

  return(invisible(NULL))

}
