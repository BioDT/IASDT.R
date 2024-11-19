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
    LF_Return = TRUE, Verbose = TRUE) {

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
    crossprod_solve <- Alpha_ID <- warmup <- NULL

  # # ..................................................................... ###

  # indices of unitsPred in modelunits
  indOld <- (unitsPred %in% modelunits)
  # indices of new unitsPred
  indNew <- !(indOld)

  if (sum(indNew) == 0 && sum(indOld) == length(modelunits)) {

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

    IASDT.R::CatTime("At least some of input sites are new sites", Level = 1)

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

    IASDT.R::CatTime("Calculate/save D11 and D12 distance matrices", Level = 1)

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
        rm(rL, D11, D12, envir = environment())
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


    if (Temp_Cleanup) {
      if (Model_Name != "") {
        on.exit(
          try(
            expr = {
              list.files(
                Temp_Dir,
                pattern = paste0("^", Model_Name, "_postEta"),
                full.names = TRUE) %>%
                c(Path_s1, Path_s2) %>%
                fs::file_delete()
            },
            silent = TRUE),
          add = TRUE)
      } else {
        on.exit(
          try(fs::file_delete(c(Path_s1, Path_s2)), silent = TRUE),
          add = TRUE)
      }
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

    etaPreds_F <- function(RowNum) {
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

      if (Denom > 0) {

        if (UseTF) {
          # Use TensorFlow
          if (isFALSE(IASDT.R::CheckData(File_etaPred, warning = FALSE))) {
            eta_indNew <- crossprod_solve(
              s1 = Path_s1, s2 = Path_s2, denom = Denom,
              postEta = File, use_single = TF_use_single, save = TRUE,
              file_path = File_etaPred, verbose = FALSE)
            # saveRDS(eta_indNew, file = File_etaPred)
          } else {
            eta_indNew <- readRDS(File_etaPred)
          }

          eta_indNew <- purrr::map(
            .x = seq_along(eta_indNew),
            .f = ~ {
              tibble::tibble(etaPred = as.vector(eta_indNew[[.x]])) %>%
                dplyr::mutate(
                  Unit_ID = unitsPred[indNew], SampleID = SampleID[.x])
            }) %>%
            dplyr::bind_rows()

          postEta <- IASDT.R::LoadAs(postEta_File, nthreads = 5)

          eta_indOld <- postEta[SampleID] %>%
            purrr::map(~ .x[match(unitsPred[indOld], modelunits), LF_ID])
          eta_indOld <- purrr::map(
            .x = seq_along(eta_indOld),
            .f = ~ {
              eta_indOld[[.x]] %>%
                tibble::tibble(etaPred = .) %>%
                dplyr::mutate(
                  Unit_ID = unitsPred[indOld], SampleID = SampleID[.x])
            }) %>%
            dplyr::bind_rows()

          etaPred <- dplyr::bind_rows(eta_indOld, eta_indNew) %>%
            dplyr::select(c("SampleID", "etaPred", "Unit_ID")) %>%
            dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
            dplyr::arrange(SampleID, Unit_ID, etaPred)

        } else {

          # Use R / CPP

          # Reading postEta from file
          postEta0 <- IASDT.R::LoadAs(File, nthreads = 5)

          # Read D11 and D12
          D11 <- IASDT.R::LoadAs(Path_D11, nthreads = 5)
          D12 <- IASDT.R::LoadAs(Path_D12, nthreads = 5)

          K11 <- IASDT.R::exp_neg_div(D11, Denom)
          K12 <- IASDT.R::exp_neg_div(D12, Denom)

          etaPred <- purrr::map_dfr(
            .x = seq_along(SampleID),
            .f = function(ID) {
              eta <- postEta0[, , ID]
              eta_indNew <- IASDT.R::Solve2vect(K11, eta) %>%
                as.vector() %>%
                Matrix::crossprod(K12, .) %>%
                as.vector()

              etaPred <- rep(NA, length(unitsPred))
              etaPred[indOld] <- eta[match(unitsPred[indOld], modelunits)]
              etaPred[indNew] <- eta_indNew

              tibble::tibble(
                SampleID = SampleID[ID], etaPred = etaPred,
                Unit_ID = unitsPred)
            }) %>%
            dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
            dplyr::arrange(SampleID, Unit_ID, etaPred)
        }

      } else {

        # Handle cases where Denom is zero by setting `eta_indNew` to zero

        if (isFALSE(IASDT.R::CheckData(File_etaPred, warning = FALSE))) {

          postEta0 <- IASDT.R::LoadAs(File, nthreads = 5)

          etaPred <- purrr::map_dfr(
            .x = seq_len(length(SampleID)),
            .f = function(ID) {
              eta <- postEta0[, , ID]
              etaPred <- rep(NA, length(unitsPred))
              etaPred[indOld] <- eta[match(unitsPred[indOld], modelunits)]
              etaPred[indNew] <- 0
              tibble::tibble(
                SampleID = SampleID[ID], etaPred = etaPred,
                Unit_ID = unitsPred)
            }) %>%
            dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
            dplyr::arrange(SampleID, Unit_ID, etaPred)

          saveRDS(etaPred, file = File_etaPred)

        } else {
          etaPred <- readRDS(File_etaPred)
        }
      }

      return(etaPred)
    }

    invisible(gc())

    # # .................................................................... ###
    # # .................................................................... ###

    # Predict latent factors

    if (LF_NCores == 1) {

      if (UseTF) {
        # Suppress TensorFlow warnings and disable optimizations
        Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

        # Activate the python environment
        reticulate::use_virtualenv(TF_Environ, required = TRUE)

        # Source the script file containing the crossprod_solve function
        PythonScript <- system.file("crossprod_solve.py", package = "IASDT.R")
        reticulate::source_python(PythonScript)

        # A lightweight function to initialize necessary modules.
        warmup()
      }

      # Sequential processing
      IASDT.R::CatTime("Predicting Latent Factor sequentially", Level = 1)

      # Making predictions sequentially
      etaPreds <- purrr::map(
        .x = seq_len(nrow(Unique_Alpha)),
        .f = ~ {
          print(.x)
          purrr::possibly(etaPreds_F(.x))
        })

    } else {

      # Parallel processing
      IASDT.R::CatTime(
        paste0(
          "Predicting Latent Factor in parallel using ", LF_NCores, " cores"),
        Level = 1)

      IASDT.R::CatTime("Prepare for parallel processing", Level = 2)
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, cluster.timeout = 10 * 60,
        future.gc = TRUE, future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(LF_NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)

      IASDT.R::CatTime("Export objects to cores", Level = 2)
      snow::clusterExport(
        cl = c1,
        list = c(
          "Unique_Alpha", "Path_D11", "Path_D12", "Path_s1", "Path_s2",
          "indNew", "unitsPred", "postEta_File", "indOld", "modelunits",
          "TF_Environ", "UseTF", "TF_use_single", "etaPreds_F"),
        envir = environment())

      # Load necessary libraries and load environment if using TensorFlow
      IASDT.R::CatTime(
        "Load necessary libraries and load environment", Level = 2)
      invisible(snow::clusterEvalQ(
        cl = c1,
        expr = {
          sapply(
            c(
              "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble",
              "Matrix", "Hmsc", "qs", "fs", "purrr", "IASDT.R", "reticulate"),
            library, character.only = TRUE)

          if (UseTF) {

            # Suppress TensorFlow warnings and disable optimizations
            Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

            # Activate the python environment
            reticulate::use_virtualenv(TF_Environ, required = TRUE)

            # Source the script file containing the crossprod_solve function
            PythonScript <- system.file(
              "crossprod_solve.py", package = "IASDT.R")
            reticulate::source_python(PythonScript)

            # A lightweight function to initialize necessary modules.
            warmup()
          }
        }))

      # Making predictions on parallel
      IASDT.R::CatTime("Making predictions on parallel", Level = 2)
      etaPreds <- snow::clusterApplyLB(
        cl = c1, x = seq_len(nrow(Unique_Alpha)),
        fun = purrr::possibly(
          function(x) {
            max_tries <- 5
            attempt <- 1
            while (attempt <= max_tries) {
              result <- tryCatch({
                return(etaPreds_F(x))
                # If successful, return the result
              }, error = function(e) {
                # Return NULL on error to retry
                NULL
              })
            }

            if (!is.null(result)) {
              # If successful, exit the loop
              return(result)
            }
          }
        ))

      # Stop the cluster
      IASDT.R::CatTime("Stop the cluster", Level = 2)
      snow::stopCluster(c1)
      invisible(gc())
    }

    # Check if all files are created
    IASDT.R::CatTime("Check if all files are created", Level = 1)
    AllEtaFiles <- Unique_Alpha$File_etaPred
    AllEtaFilesExist <- all(file.exists(AllEtaFiles))
    if (isFALSE(AllEtaFilesExist)) {
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
