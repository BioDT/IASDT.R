## |------------------------------------------------------------------------| #
# Predict_LF ----
## |------------------------------------------------------------------------| #

#' Predict_LF
#'
#' Draws samples from the conditional predictive distribution of latent factors.
#' #' This function is optimized for speed using parallel processing and
#' optionally TensorFlow for matrix operations. This function is adapted from
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
#' @param rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param LF_NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param Temp_Dir Character string specifying the path for temporary storage of
#'   intermediate files.
#' @param Temp_Cleanup Logical indicating whether to delete temporary files in
#'   the `Temp_Dir` after finishing the calculations. Defaults to `FALSE`; i.e.,
#'   the temporary files will be kept, enabling continue working on unfinished
#'   calculations, if the function failed, e.g. due to OOM.
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
#'   This should end with either `qs` or `RData`.
#' @param LF_Return Logical. Indicates if the output should be returned.
#'   Defaults to `TRUE`. If `LF_OutFile` is `NULL`, this parameter cannot be set
#'   to `FALSE` because the function needs to return the result if it is not
#'   saved to a file.
#' @param LF_Check Logical. If TRUE, the function checks if the output
#'   files are already created and valid. If FALSE, the function will only check
#'   if the files exist without checking their integrity. Default is `FALSE`.
#' @param Verbose Logical. If TRUE, detailed output is printed. Default is
#'   `FALSE`.
#' @export
#' @author This script was adapted from the [Hmsc::predictLatentFactor] function
#'   in the `Hmsc` package.
#' @seealso [Hmsc::predictLatentFactor]
#' @name Predict_LF
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
    unitsPred, modelunits, postEta, postAlpha, rL, LF_NCores = 8,
    Temp_Dir = "TEMP2Pred", Temp_Cleanup = FALSE, Model_Name = NULL,
    UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE, LF_OutFile = NULL,
    LF_Return = TRUE, LF_Check = FALSE, Verbose = TRUE) {

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
    postEta_File <- postEta
    IASDT.R::CatTime("Load postEta", Level = 1)
    if (!file.exists(postEta)) {
      stop("The specified path for `postEta` does not exist. ", call. = FALSE)
    }
    postEta <- IASDT.R::LoadAs(postEta, nthreads = 5)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SampleID <- Unit_ID <- LF <- LF_ID <- etaPred <- Sample_IDs <- File <-
    Alpha_ID <- NULL

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

  AllTraining <- sum(indNew) == 0 && sum(indOld) == length(modelunits)
  AllNew <- sum(indNew) == length(modelunits) && sum(indOld) == 0

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

      # Extension for temporary files
      TempExt <- "rds"

    } else {

      IASDT.R::CatTime("Computations will be made using R/CPP", Level = 1)

      # Extension for temporary files
      TempExt <- "qs"
    }

    # # .................................................................... ###

    # Adjust Model_Name prefix

    if (is.null(Model_Name)) {
      Model_Name <- ""
    } else {
      Model_Name <- paste0(Model_Name, "_")
    }

    # Create a temporary directory to store intermediate results. This directory
    # will be used to save D11, D12, and intermediate postEta files, reducing
    # memory usage.
    fs::dir_create(Temp_Dir)

    # # .................................................................... ###

    # Calculate D11 and D12 only once

    IASDT.R::CatTime("Calculate/save necessary matrices", Level = 1)

    alphapw <- rL$alphapw

    if (UseTF) {

      # Save s1 and s2 as rds files, if not already exist on disk
      Path_s1 <- file.path(Temp_Dir, paste0(Model_Name, "s1.rds"))
      Path_s2 <- file.path(Temp_Dir, paste0(Model_Name, "s2.rds"))

      if (file.exists(Path_s1) && file.exists(Path_s2)) {
        IASDT.R::CatTime("s1 and s2 matrices are already saved", Level = 2)
      } else {
        IASDT.R::CatTime("Saving s1 and s2 matrices", Level = 2)
        s1 <- rL$s[modelunits, , drop = FALSE]
        s2 <- rL$s[unitsPred[indNew], , drop = FALSE]
        saveRDS(s1, file = Path_s1)
        saveRDS(s2, file = Path_s2)
      }
      rm(rL, envir = environment())

      Path_D11 <- Path_D12 <- NULL

    } else {

      # Save D11 and D12 as qs/rds files, if not already exist on disk
      Path_D11 <- file.path(Temp_Dir, paste0(Model_Name, "D11.", TempExt))
      Path_D12 <- file.path(Temp_Dir, paste0(Model_Name, "D12.", TempExt))

      if (file.exists(Path_D11) && file.exists(Path_D12)) {

        IASDT.R::CatTime(
          "D11 and D12 distance matrices are already saved", Level = 2)

      } else {

        s1 <- rL$s[modelunits, , drop = FALSE]
        s2 <- rL$s[unitsPred[indNew], , drop = FALSE]
        D11 <- Rfast::Dist(s1)
        D12 <- Rfast::dista(s1, s2)

        if (UseTF) {
          saveRDS(D11, file = Path_D11)
          saveRDS(D12, file = Path_D12)
        } else {
          qs::qsave(D11, file = Path_D11, preset = "fast")
          qs::qsave(D12, file = Path_D12, preset = "fast")
        }

        # Clean up
        rm(rL, s1, s2, D11, D12, envir = environment())
      }
      Path_s1 <- Path_s2 <- NULL
    }


    # Clean up temporary files after finishing calculations
    if (Temp_Cleanup) {
      on.exit(
        try(
          expr = {
            list.files(
              Temp_Dir, pattern = paste0("^", Model_Name, "_postEta"),
              full.names = TRUE) %>%
              c(Path_s1, Path_s2) %>%
              fs::file_delete()
          },
          silent = TRUE),
        add = TRUE)
    }

    invisible(gc())

    # # .................................................................... ###

    # Convert postAlpha to tibble

    IASDT.R::CatTime(
      "Splitting and saving `postAlpha` to small chunks", Level = 1)

    postAlpha_tibble <- do.call(rbind, postAlpha) %>%
      as.data.frame() %>%
      tibble::tibble() %>%
      stats::setNames(paste0("LF_", seq_len(ncol(.))))

    # Unique values alpha and their respective sample IDs
    postAlpha_unique <- postAlpha_tibble %>%
      # ID column represents the original row number
      dplyr::mutate(Sample_IDs = dplyr::row_number()) %>%
      tidyr::nest(Sample_IDs = Sample_IDs) %>%
      dplyr::mutate(
        Sample_IDs = purrr::map(Sample_IDs, ~ as.vector(sort(unlist(.x)))))

    # Prepare data for parallel processing
    Unique_Alpha <- postAlpha_unique %>%
      # This may help to distribute heavy jobs first on parallel
      dplyr::arrange(dplyr::desc(sapply(Sample_IDs, length))) %>%
      dplyr::select(-Sample_IDs) %>%
      tidyr::pivot_longer(
        cols = names(.), values_to = "Alpha_ID", names_to = "LF") %>%
      dplyr::mutate(
        LF_ID = as.integer(stringr::str_remove(LF, "LF_")),
        .before = "Alpha_ID") %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        Denom = purrr::map_dbl(Alpha_ID, ~ alphapw[.x, 1]),
        SampleID = purrr::map2(
          .x = LF, .y = Alpha_ID,
          .f = ~ {
            postAlpha_unique %>%
              dplyr::filter(.[[.x]] == .y) %>%
              dplyr::pull(Sample_IDs) %>%
              unlist() %>%
              as.vector() %>%
              sort()
          }),
        SampleID = purrr::map(SampleID, unlist),
        File = file.path(
          Temp_Dir,
          paste0(Model_Name, "postEta_ch", dplyr::row_number(), ".", TempExt)),
        File_etaPred = stringr::str_replace_all(
          File, "_postEta_ch", "_etaPred_ch"),
        Export = purrr::pmap(
          .l = list(SampleID, LF_ID, File),
          .f = function(SampleID, LF_ID, File) {

            # do not export file if already exists and is valid
            # if (isFALSE(IASDT.R::CheckData(File, warning = FALSE))) {
            if (!file.exists(File)) {
              Out <- postEta[SampleID] %>%
                purrr::map(~ .x[, LF_ID, drop = FALSE]) %>%
                simplify2array()
              if (UseTF) {
                saveRDS(Out, file = File)
              } else {
                qs::qsave(Out, file = File, preset = "fast")
              }
            }
          }),
        Export = NULL)

    rm(postEta, postAlpha, envir = environment())
    invisible(gc())

    # # .................................................................... ###

    # Internal functions to predict latent factors

    etaPreds_F <- function(RowNum, LF_Check = FALSE) {
      # Current denominator
      Denom <- Unique_Alpha$Denom[[RowNum]]
      # ID for latent factor
      LF_ID <- Unique_Alpha$LF_ID[[RowNum]]
      # ID for posterior sample
      SampleID <- Unique_Alpha$SampleID[[RowNum]]
      # File path for current alpha
      File <- Unique_Alpha$File[[RowNum]]
      File_etaPred <- Unique_Alpha$File_etaPred[[RowNum]]

      # If the denominator is positive, perform calculations; otherwise, set
      # `eta_indNew` to zero.
      if (LF_Check) {
        CalcPredLF <- isFALSE(IASDT.R::CheckData(File_etaPred, warning = FALSE))
      } else {
        CalcPredLF <- !file.exists(File_etaPred)
      }

      if (Denom > 0) {

        if (UseTF) {

          # Suppress TensorFlow warnings and disable optimizations
          Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

          # Use TensorFlow
          if (CalcPredLF) {

            PythonScript <- system.file(
              "crossprod_solve.py", package = "IASDT.R")
            
            eta_indNew0 <- run_crossprod_solve(
              virtual_env_path = TF_Environ, script_path = PythonScript,
              s1 = Path_s1, s2 = Path_s2, denom = Denom, postEta = File,
              path_out = File_etaPred, use_single = TF_use_single)
          }

          etaPred <- IASDT.R::LoadAs(File_etaPred) %>%
            purrr::map( ~ tibble::tibble(Post = unlist(.x))) %>%
            tibble::tibble(DT = .) %>%
            dplyr::mutate(SampleID = SampleID) %>%
            tidyr::unnest(DT)

        } else {

          # Use R / CPP

          if (CalcPredLF) {
            
            # Reading postEta from file
            postEta0 <- IASDT.R::LoadAs(File, nthreads = 5)

            # Read D11 and D12
            D11 <- IASDT.R::LoadAs(Path_D11, nthreads = 5)
            D12 <- IASDT.R::LoadAs(Path_D12, nthreads = 5)

            K11 <- IASDT.R::exp_neg_div(D11, Denom)
            K12 <- IASDT.R::exp_neg_div(D12, Denom)


            # NEED TO REVISE OUTPUT COLUMNS

            etaPred <- purrr::map_dfr(
              .x = seq_along(SampleID),
              .f = function(ID) {
                eta <- postEta0[, , ID]
                etaPred <- IASDT.R::Solve2vect(K11, eta) %>%
                  as.vector() %>%
                  Matrix::crossprod(K12, .) %>%
                  as.vector() %>%
                  tibble::tibble(etaPred = ., SampleID = SampleID)
              }) %>%
              dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
              dplyr::arrange(SampleID, Unit_ID, etaPred)

            qs::qsave(etaPred, file = File_etaPred, preset = "fast")

          } else {
            
            etaPred <- IASDT.R::LoadAs(File_etaPred)

          }
        }

      } else {
        # When Denom is zero, set `eta_indNew` to zero

        # NEED TO REVISE OUTPUT COLUMNS

        etaPred <- tibble::tibble(etaPred = 0, SampleID = SampleID)
      }

      return(etaPred)
    }

    invisible(gc())

    # # .................................................................... ###
    # # .................................................................... ###

    # Predict latent factors

    if (LF_NCores == 1) {

      # Sequential processing
      IASDT.R::CatTime("Predicting Latent Factor sequentially", Level = 1)

      # Making predictions sequentially
      etaPreds <- purrr::map(
        .x = seq_len(nrow(Unique_Alpha)),
        .f = function(x) {
          # maximum number of attempts
          max_tries <- 5
          attempt <- 1
          # Initialize result
          result <- NULL

          while (attempt <= max_tries) {
            result <- tryCatch({
                # Use purrr::possibly around etaPreds_F call to handle errors
                etaPreds_F(x, LF_Check = LF_Check)
              },
              error = function(e) {
                # Return NULL on error to retry
                NULL
              })

            # Exit loop if result is not NULL (successful)
            if (!is.null(result)) {
              break
            }

            attempt <- attempt + 1
          }
          # Return the result after max_tries or successful execution
          return(result)
        })

    } else {

      # Parallel processing
      IASDT.R::CatTime(
        paste0(
          "Predicting Latent Factor in parallel using ", 
          LF_NCores, " cores"),
        Level = 1)

      IASDT.R::CatTime("Prepare for parallel processing", Level = 2)
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, cluster.timeout = 10 * 60,
        future.gc = TRUE, future.seed = TRUE)
      c1 <- parallel::makeCluster(LF_NCores)
      on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

      IASDT.R::CatTime("Export objects to cores", Level = 2)
      parallel::clusterExport(
        cl = c1,
        varlist = c(
          "Unique_Alpha", "Path_D11", "Path_D12", "Path_s1", "Path_s2",
          "indNew", "unitsPred", "postEta_File", "indOld", "modelunits",
          "TF_Environ", "UseTF", "TF_use_single", "etaPreds_F",
          "LF_Check", "run_crossprod_solve"),
        envir = environment())

      # Load necessary libraries and load environment if using TensorFlow
      IASDT.R::CatTime("Load necessary libraries", Level = 2)
      invisible(parallel::clusterEvalQ(
        cl = c1,
        expr = {
          sapply(
            c(
              "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble",
              "Matrix", "Hmsc", "qs", "fs", "purrr", "IASDT.R"),
            library, character.only = TRUE)
          invisible(gc())
        }))

      # Making predictions on parallel
      IASDT.R::CatTime("Making predictions on parallel", Level = 2)
      etaPreds <- parallel::clusterApplyLB(
        cl = c1, x = seq_len(nrow(Unique_Alpha)),
        fun = function(x) {
          # maximum number of attempts
          max_tries <- 5
          attempt <- 1
          # Initialize result
          result <- NULL

          while (attempt <= max_tries) {
            result <- tryCatch({
                # Use purrr::possibly around etaPreds_F call to handle errors
                etaPreds_F(x, LF_Check = LF_Check)
              },
              error = function(e) {
                # Return NULL on error to retry
                NULL
              })

            # Exit loop if result is not NULL (successful)
            if (!is.null(result)) {
              break
            }

            attempt <- attempt + 1
          }

          # Return the result after max_tries or successful execution
          return(result)
        })

      # Stop the cluster
      IASDT.R::CatTime("Stop the cluster", Level = 2)
      parallel::stopCluster(c1)
      invisible(gc())
    }

    # Check if all files are created
    IASDT.R::CatTime("Check if all files are created", Level = 1)
    AllEtaFiles <- Unique_Alpha$File_etaPred
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

    # # .................................................................... ###

    # Merge results
    IASDT.R::CatTime("Merge results", Level = 1)

    postEtaPred <- Unique_Alpha %>%
      dplyr::select(LF, LF_ID) %>%
      dplyr::mutate(etaPred = etaPreds) %>%
      tidyr::unnest("etaPred") %>%
      tidyr::pivot_wider(
        id_cols = c(SampleID, Unit_ID),
        names_from = LF, values_from = etaPred) %>%
      dplyr::arrange(SampleID, Unit_ID) %>%
      dplyr::group_split(SampleID) %>%
      purrr::map(
        .f = ~ {
          Mat <- dplyr::select(.x, tidyselect::starts_with("LF_")) %>%
            as.matrix()
          rownames(Mat) <- unitsPred
          colnames(Mat) <- NULL
          return(Mat)
        })
  }

  # # ..................................................................... ###

  # Save postEtaPred

  if (!is.null(LF_OutFile)) {
    IASDT.R::CatTime(
      paste0("Save postEtaPred to `", LF_OutFile, "`"), Level = 1)
    fs::dir_create(fs::path_dir(LF_OutFile))
    switch(
      fs::path_ext(LF_OutFile),
      "qs" = qs::qsave(postEtaPred, file = LF_OutFile, preset = "fast"),
      "RData" = save(postEtaPred, file = LF_OutFile),
      stop("Unsupported file extension in `LF_OutFile`.", call. = FALSE))
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
#' @param script_path character. Path to the Python script to execute.
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
#' @return Returns the `path_out` if successful. Returns `NULL` if all attempts
#'   fail.
#'
#' @details
#' - The function checks for the existence of required input files and the
#' Python executable in the specified virtual environment.
#' - Executes the Python script using `system2`.
#' - Verifies the output file validity using `IASDT.R::CheckData`. Retries up
#' to 3 times if the output is invalid.
#' - Generates detailed logs if `verbose` is set to `TRUE`.
#'
#' @noRD

run_crossprod_solve <- function(
    virtual_env_path, script_path, s1, s2, postEta, path_out,
    denom, chunk_size = 1000, threshold_mb = 2000,
    use_single = TRUE, verbose = TRUE) {

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")

  if (is.null(script_path)) {
    script_path <- system.file("crossprod_solve.py", package = "IASDT.R")
  }


  # Ensure the paths are valid
  paths <- list(virtual_env_path, script_path, s1, s2, postEta)
  names(paths) <- c(
    "Virtual Environment", "Python Script", "s1", "s2", "postEta")
  lapply(names(paths), function(p) {
    if (!file.exists(paths[[p]])) {
      stop(paste0(p, " does not exist: ", paths[[p]]))
    }
  })

  # Determine the Python executable path
  python_executable <- if (.Platform$OS.type == "windows") {
    file.path(virtual_env_path, "Scripts", "python.exe")
  } else {
    file.path(virtual_env_path, "bin", "python")
  }
  if (!file.exists(python_executable)) {
    stop(
      "Python executable not found in the virtual environment.",
      call. = FALSE)
  }

  # Construct the command to run the Python script
  args <- c(
    script_path,
    "--s1", normalizePath(s1),
    "--s2", normalizePath(s2),
    "--postEta", normalizePath(postEta),
    "--path_out", normalizePath(path_out),
    "--denom", as.character(denom),
    "--chunk_size", as.character(chunk_size),
    "--threshold_mb", as.character(threshold_mb),
    "--save")

  # Add boolean flags conditionally
  if (use_single) {
    args <- c(args, "--use_single")
  }
  if (verbose) {
    args <- c(args, "--verbose")
  }

  # Initialize retry logic
  max_attempts <- 3
  attempt <- 1
  success <- FALSE

  while (attempt <= max_attempts && !success) {
    # Run the command and capture stdout/stderr
    result <- system2(
      command = python_executable, args = args, 
      stdout = TRUE, stderr = TRUE)

    # Check for errors
    if (!inherits(result, "error") || !length(result) == 0) {
      # Check if file is valid using IASDT.R::CheckData
      FileOkay <- IASDT.R::CheckData(path_out)
      if (FileOkay) {
        success <- TRUE
      }
    }

    attempt <- attempt + 1
  }

  if (verbose) {
    path_log <- stringr::str_replace(path_out, ".rds", ".log")
    f <- file(path_log, open = "wb")
    on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
    cat(
      "Running command:\n",
      paste(python_executable, paste(args, collapse = " "), "\n\n"),
      sep = "\n", file = f)
    cat(result[-length(result)], sep = "\n", file = f)
    # close connection to the file
    close(f)
  }

  # If all attempts fail, return NULL
  if (!success) {
    if (verbose) {
      cat("All attempts failed. Returning NULL.\n")
    }
    return(NULL)
  } else {
    return(path_out)
  }
}
