## |------------------------------------------------------------------------| #
# Predict_LF ----
## |------------------------------------------------------------------------| #

#' Predict_LF
#'
#' Draws samples from the conditional predictive distribution of latent factors.
#' This function is optimized for speed using parallel processing and optionally
#' TensorFlow for matrix operations. This function is adapted from
#' [Hmsc::predictLatentFactor] with equivalent results to the original function
#' when `predictMean = TRUE`.
#'
#' @param unitsPred a factor vector with random level units for which
#'   predictions are to be made
#' @param modelunits a factor vector with random level units that are
#'   conditioned on
#' @param postEta Character string specifying the path for postEta; a list
#'   containing samples of random factors at conditioned units
#' @param postAlpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param LF_rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param LF_NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param Temp_Dir Character string specifying the path for temporary storage of
#'   intermediate files.
#' @param LF_Temp_Cleanup Logical indicating whether to delete temporary files
#'   in the `Temp_Dir` after finishing the LF predictions.
#' @param Model_Name Character string used as a prefix for temporary file names.
#'   Defaults to NULL, in which case no prefix is used.
#' @param UseTF Logical indicating whether to use TensorFlow for calculations.
#'   Defaults to TRUE.
#' @param TF_Environ Character string specifying the path to the Python
#'   environment. Defaults to NULL. This argument is required if `UseTF` is
#'   TRUE.
#' @param TF_use_single Logical indicating whether to use single precision for
#'   the TF calculations. Defaults to `FALSE`.
#' @param LF_OutFile Character string specifying the path to save the outputs.
#'   If `NULL` (default), the predicted latent factors are not saved to a file.
#'   This should end with either `*.qs2` or `*.RData`.
#' @param LF_Return Logical. Indicates if the output should be returned.
#'   Defaults to `FALSE`. If `LF_OutFile` is `NULL`, this parameter cannot be
#'   set to `FALSE` because the function needs to return the result if it is not
#'   saved to a file.
#' @param LF_Check Logical. If TRUE, the function checks if the output files are
#'   already created and valid. If FALSE, the function will only check if the
#'   files exist without checking their integrity. Default is `FALSE`.
#' @param Verbose Logical. If TRUE, detailed output is printed. Default is
#'   `FALSE`.
#' @param solve_max_attempts numeric (Optional). Maximum number of attempts to
#'   run solve and crossprod functions. Default is 5.
#' @export
#' @author This script was adapted from the [Hmsc::predictLatentFactor] function
#'   in the `Hmsc` package.
#' @seealso [Hmsc::predictLatentFactor]
#' @name Predict_LF
#' @inheritParams run_crossprod_solve
#' @details The function is expected to be faster than the original function in
#'   the `Hmsc` package, especially when using TensorFlow for calculations and
#'   when working on parallel.
#'
#'   The main difference is that this function:
#' - allow for parallel processing (`LF_NCores` argument);
#' - it is possible to use TensorFlow (`UseTF` argument) to make matrix
#'   calculations faster, particularly when used on GPU. The following modules
#'   are needed: `numpy`, `os`, `tensorflow`, `rdata`, `xarray`, and `pandas`.
#'   To use `TensorFlow`, the argument `TF_Environ` should be set to the path of
#'   a Python environment with TensorFlow installed;
#' - if `UseTF` is set to `FALSE`, the function uses R / CPP code in the
#'   calculations;
#' - calculates `D11` and `D12` matrices only once and save them to disk and
#'   call them when needed.

Predict_LF <- function(
    unitsPred, modelunits, postEta, postAlpha, LF_rL, LF_NCores = 8L,
    Temp_Dir = "TEMP2Pred", LF_Temp_Cleanup = TRUE, Model_Name = NULL,
    UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE, LF_OutFile = NULL,
    LF_Return = FALSE, LF_Check = FALSE, LF_Commands_Only = FALSE,
    solve_max_attempts = 5L, solve_chunk_size = 50L, Verbose = TRUE) {

  # # ..................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  IASDT.R::CatTime("Starting `Predict_LF` function", Level = 1)

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

  # indices of unitsPred in modelunits
  indOld <- (unitsPred %in% modelunits)
  # indices of new unitsPred
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
      call. = FALSE
    )
  }

  if (AllTraining) {
    # If all input sites are for training sites, use LF info from the model
    # directly
    IASDT.R::CatTime("All input sites are training sites", Level = 1)

    postEtaPred <- purrr::map(
      .x = postEta,
      .f = ~ {
        Out <- .x[match(unitsPred[indOld], modelunits), ]
        rownames(Out) <- modelunits
        return(Out)
      })

  } else {

    IASDT.R::CatTime("All input sites are new sites", Level = 1)

    # Check TensorFlow settings

    if (UseTF) {

      # Check if TF_Environ directory exists
      if (is.null(TF_Environ) || !dir.exists(TF_Environ)) {
        stop(
          paste0(
            "When `UseTF` is TRUE, `TF_Environ` must be specified and should ",
            "point to an existing directory with a Python environment"),
          call. = FALSE)
      }

      PythonScript <- system.file("crossprod_solve.py", package = "IASDT.R")

      # Check if PythonScript exists
      if (!file.exists(PythonScript)) {
        stop("Necessary python script does not exist", call. = FALSE)
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

    if (stringr::str_detect(Model_Name, "^RC_c_") && LF_Commands_Only) {
        CommandFilePrefix <- "LF_RC_Commands_"
    }
    if (nchar(Model_Name) > 0 && !stringr::str_detect(Model_Name, "^Rc_c") && 
      LF_Commands_Only) {
      CommandFilePrefix <- "LF_NewSites_Commands_"
    }   

    # Create a temporary directory to store intermediate results. This directory
    # will be used to save s1/s2 or D11/D12, and intermediate postEta files,
    # reducing memory usage.
    Temp_Dir_LF <- file.path(Temp_Dir, "LF_Prediction")
    fs::dir_create(c(Temp_Dir_LF, Temp_Dir))

    # # .................................................................... ###

    # Calculate D11 and D12 only once

    IASDT.R::CatTime("Calculate/save necessary matrices", Level = 1)

    alphapw <- LF_rL$alphapw

    if (UseTF) {

      # Save s1 and s2 for coordinates at training and testing sites as feather
      # files, if not already exist on disk
      Path_s1 <- file.path(Temp_Dir, paste0(Model_Name, "s1.feather"))
      Path_s2 <- file.path(Temp_Dir, paste0(Model_Name, "s2.feather"))

      s1_s2_Okay <- IASDT.R::CheckData(Path_s1, warning = FALSE) &&
        IASDT.R::CheckData(Path_s2, warning = FALSE)

      if (s1_s2_Okay) {
        IASDT.R::CatTime("s1 and s2 matrices were already saved", Level = 2)
      } else {

        IASDT.R::CatTime("Saving s1 and s2 matrices", Level = 2)

        # s1
        s1 <- as.data.frame(LF_rL$s[modelunits, , drop = FALSE])
        IASDT.R::SaveAs(InObj = s1, OutPath = Path_s1)

        # s2
        s2 <- as.data.frame(LF_rL$s[unitsPred[indNew], , drop = FALSE])
        IASDT.R::SaveAs(InObj = s2, OutPath = Path_s2)

        rm(s1, s2, envir = environment())
      }

      rm(LF_rL, envir = environment())
      Path_D11 <- Path_D12 <- NULL

    } else {

      # Save D11 and D12 as feather files, if not already exist on disk
      Path_D11 <- file.path(Temp_Dir, paste0(Model_Name, "D11.qs2"))
      Path_D12 <- file.path(Temp_Dir, paste0(Model_Name, "D12.qs2"))

      if (file.exists(Path_D11) && file.exists(Path_D12)) {

        IASDT.R::CatTime(
          "D11 and D12 distance matrices are already saved", Level = 2)

      } else {

        s1 <- LF_rL$s[modelunits, , drop = FALSE]
        s2 <- LF_rL$s[unitsPred[indNew], , drop = FALSE]

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
            file.path(
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

    etaPreds_F <- function(RowNum, LF_Check = FALSE, unitsPred) {

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
        CalcPredLF <- dplyr::if_else(
          Eta_NCols == LF_Data$NSamples[[RowNum]], FALSE, TRUE)
      }

      if (CalcPredLF) {

        # If the denominator is positive, perform calculations; otherwise, set
        # `eta_indNew` to zero.

        if (Denom > 0) {

          if (UseTF) {

            # Use TensorFlow

            # Suppress TensorFlow warnings and disable optimizations
            Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

            eta_indNew0 <- run_crossprod_solve(
              virtual_env_path = TF_Environ, s1 = Path_s1, s2 = Path_s2,
              denom = Denom, postEta = File_postEta, path_out = File_etaPred_TF,
              use_single = TF_use_single, LF_Commands_Only = LF_Commands_Only,
              solve_max_attempts = solve_max_attempts,
              solve_chunk_size = solve_chunk_size)

            if (LF_Commands_Only) {
              return(eta_indNew0)
            } else {
              rm(eta_indNew0, envir = environment())

              etaPred <- tibble::tibble(
                SampleID = SampleID,
                LF = LF_ID,
                Path_Samp_LF = file.path(
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
                        dplyr::mutate(LF = NULL, unitsPred = unitsPred) %>%
                        stats::setNames(
                          c("SampleID", paste0("LF_", LF_ID), "unitsPred"))
                    }),
                  File_etaPred = File_etaPred,
                  ChunkID = LF_Data$ChunkID[[RowNum]],
                  SampleID = SampleID,
                  Save = purrr::map2(
                    .x = eta_DT, .y = Path_Samp_LF,
                    .f = ~qs2::qs_save(.x, .y, nthreads = 5)),
                  LF = LF_ID, Save = NULL, eta_DT = NULL)

              IASDT.R::SaveAs(InObj = etaPred, OutPath = File_etaPred)
            }

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

                Path_Samp_LF <- file.path(
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
                    etaPred = ., unitsPred = unitsPred) %>%
                  stats::setNames(
                    c("SampleID", paste0("LF_", LF_ID), "unitsPred"))

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

              Path_Samp_LF = file.path(
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
                      unitsPred = unitsPred) %>%
                      stats::setNames(
                        c("SampleID", paste0("LF_", LF_ID), "unitsPred"))
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

    if (!all(file.exists(LF_Data$File_etaPred))) {

      if (LF_NCores == 1 || LF_Commands_Only) {

        if (LF_Commands_Only) {
          IASDT.R::CatTime("Prepare commands for predicting latent factors", 
          Level = 1)
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
                RowNum = x, LF_Check = LF_Check, unitsPred = unitsPred),
              silent = FALSE)

            if (inherits(result, "try-error")) {
              return(NULL)
            } else {
              return(result)
            }

          })

        if (LF_Commands_Only) {

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
              file_name <- file.path(
                dirname(Temp_Dir), paste0(CommandFilePrefix, i, ".txt"))
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
            "unitsPred", "indOld", "modelunits", "TF_Environ", "UseTF",
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
              etaPreds_F(x, LF_Check = LF_Check, unitsPred = unitsPred),
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

      if (AllEtaFilesExist) {
        IASDT.R::CatTime("All files were created", Level = 2)
      } else {
        FailedFiles <- AllEtaFiles[!file.exists(AllEtaFiles)]
        stop(
          paste0(
            length(FailedFiles), " files are missing: \n",
            paste0("  >>  ", basename(FailedFiles), collapse = "\n")),
          call. = FALSE)
      }

    } else {
      IASDT.R::CatTime(
        "All LF prediction files were already created", Level = 1)
    }

    # # .................................................................... ###

    # Merge results
    IASDT.R::CatTime("Merge results", Level = 1)

    postEtaPred_Samp <- etaPreds %>%
      dplyr::bind_rows() %>%
      dplyr::select(-LF, -ChunkID, -File_etaPred) %>%
      dplyr::arrange(SampleID) %>%
      dplyr::mutate(
        Path_Sample = file.path(
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
            .f = dplyr::left_join, by = c("SampleID", "unitsPred")) %>%
          dplyr::arrange(unitsPred) %>%
          dplyr::select(sort(tidyselect::peek_vars())) %>%
          dplyr::select(-SampleID) %>%
          magrittr::set_rownames(NULL) %>%
          tibble::column_to_rownames("unitsPred") %>%
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
          "(postEta|r[0-9]|etaPred|s1|s2|post).+(feather|qs2|log)")
        file_paths <- list.files(
          path = normalizePath(Temp_Dir, winslash = "/"),
          pattern = Pattern, full.names = TRUE)
        fs::file_delete(file_paths)
      },
      silent = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Predict_LF was finished in ", Level = 1)

  # # ..................................................................... ###

  if (LF_Return) {
    return(postEtaPred)
  } else {
    return(LF_OutFile)
  }

}

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
#' @param virtual_env_path character. Path to the virtual environment containing
#'   Python and required dependencies.
#' @param s1 character. Path to the input file containing s1 coordinates.
#' @param s2 character Path to the input file containing s2 coordinates.
#' @param postEta character. Path to the file containing the `postEta` matrix
#'   data.
#' @param path_out character. Path to rds file where the output results will be
#'   saved.
#' @param denom numeric. A denominator value used in the computation.
#' @param chunk_size numeric (Optional). Size of chunks to process at a time.
#'   Default is 1000.
#' @param threshold_mb numeric (Optional). Memory threshold (in MB) to manage
#'   processing. Default is 2000.
#' @param use_single logical. Indicates whether to use single precision.
#'   Defaults to `TRUE`.
#' @param verbose logical (Optional). If `TRUE`, logs detailed information
#'   during execution. Default is `TRUE`.
#' @param solve_max_attempts numeric (Optional). Maximum number of attempts to
#'   run solve and crossprod functions. Default is 5.
#' @param solve_chunk_size numeric. Chunk size for solve_and_multiply python
#'   function. Default is 50.
#' @param LF_Commands_Only logical. If `TRUE`, returns the command to run the
#'   Python script. Default is `FALSE`.
#' @return Returns the `path_out` if successful. Returns `NULL` if all attempts
#'   fail.
#' @details
#' - The function checks for the existence of required input files and the
#' Python executable in the specified virtual environment.
#' - Executes the Python script using `system2`.
#' - Verifies the output file validity using `IASDT.R::CheckData`. Retries up
#' to 3 times if the output is invalid.
#' - Generates detailed logs if `verbose` is set to `TRUE`.
#'
#' @keywords internal

run_crossprod_solve <- function(
    virtual_env_path, s1, s2, postEta, path_out, denom,
    chunk_size = 1000L, threshold_mb = 2000L, use_single = TRUE, verbose = TRUE,
    solve_chunk_size = 50L, solve_max_attempts = 5L, LF_Commands_Only = FALSE) {

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")

  # do not use scientific notation
  withr::local_options(scipen = 99)

  script_path <- system.file("crossprod_solve.py", package = "IASDT.R")
  if (!file.exists(script_path)) {
    stop(
      "Necessary python script `crossprod_solve.py` does not exist",
      call. = FALSE)
  }

  # Ensure the paths are valid
  paths <- list(virtual_env_path, script_path, s1, s2, postEta)
  names(paths) <- c(
    "Virtual Environment", "Python Script", "s1", "s2", "postEta")
  purrr::walk(
    .x = names(paths),
    .f = function(p) {
      if (!file.exists(paths[[p]])) {
        stop(paste0(p, " does not exist: ", paths[[p]]))
      }
    })

  # Determine the Python executable path
  if (.Platform$OS.type == "windows") {
    python_executable <- normalizePath(
      file.path(virtual_env_path, "Scripts", "python.exe"),
      winslash = "/")
    if (!file.exists(python_executable)) {
      stop(
        "Python executable not found in the virtual environment.",
        call. = FALSE)
    }
  } else {
    #python_executable <- "python3"
    #python_executable <- normalizePath(
    #file.path(virtual_env_path, "bin", "python3"), winslash = "/")
    # system("module use /appl/local/csc/modulefiles; module load tensorflow")
    # Sys.setenv(PATH = file.path(virtual_env_path, "bin"))
    # Sys.setenv(PYTHONPATH = file.path(
    #   virtual_env_path, "lib/python3.10/site-packages"))
    # Sys.setenv(VIRTUAL_ENV = virtual_env_path)
    # python_executable <- "/pfs/lustrep3/appl/local/csc/soft/ai/bin/python3"
    python_executable <- "python3"

  }

  # Construct the command to run the Python script
  LF_Args <- c(
    python_executable,
    script_path,
    "--s1", normalizePath(s1, winslash = "/"),
    "--s2", normalizePath(s2, winslash = "/"),
    "--post_eta", normalizePath(postEta, winslash = "/"),
    "--path_out", file.path(getwd(), path_out),
    "--denom", as.character(denom),
    "--chunk_size", as.character(chunk_size),
    "--threshold_mb", as.character(threshold_mb),
    "--solve_chunk_size", as.character(solve_chunk_size))

  # Add boolean flags conditionally
  if (use_single) {
    LF_Args <- c(LF_Args, "--use_single")
  }

  if (verbose) {
    LF_Args <- c(LF_Args, "--verbose")
  }

  LF_Args <- paste0(LF_Args, collapse = " ")

  if (LF_Commands_Only) {
    if (.Platform$OS.type != "windows") {
      LF_Args <- paste0("/usr/bin/time -v ", LF_Args)
    }
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
    if (!success) {
      if (verbose) {
        cat("All attempts failed. Returning NULL.\n",
            sep = "\n", file = f, append = TRUE)
        # close connection to the file
        close(f)
      }
      return(NULL)
    } else {
      # close connection to the file
      close(f)
      return(path_out)
    }
  }
}
