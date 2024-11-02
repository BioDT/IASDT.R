## |------------------------------------------------------------------------| #
# predictLF ----
## |------------------------------------------------------------------------| #

#' predictLF
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
#' @param postEta a list containing samples of random factors at conditioned
#'   units
#' @param postAlpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param TempDir Character string specifying the path for temporary storage of
#'   intermediate files.
#' @param ModelName Character string used as a prefix for temporary file names.
#'   Defaults to NULL, in which case no prefix is used.
#' @param UseTF Logical indicating whether to use TensorFlow for calculations.
#'   Defaults to TRUE.
#' @param PythonScript Character string specifying the path to the Python script
#'   file containing the `crossprod_solve` function. Defaults to the path of the
#'   file included in the package.
#' @param EnvPath Character string specifying the path to the Python
#'   environment. Defaults to NULL. This argument is required if `UseTF` is
#'   TRUE.
#' @param use_single Logical indicating whether to use single precision for the
#'   TF calculations. Defaults to `FALSE`.
#' @param Path_postEtaPred Character string specifying the path to save the
#'   outputs. If `NULL` (default), the predicted latent factors are not saved to
#'   a file.
#' @param Return_postEtaPred Logical. Indicates if the output should be
#'   returned. Defaults to `TRUE`. If `Path_postEtaPred` is `NULL`, this
#'   parameter cannot be set to `FALSE` because the function needs to return the
#'   result if it is not saved to a file.
#' @export
#' @author This script was adapted from the [Hmsc::predictLatentFactor] function
#'   in the `Hmsc` package.
#' @seealso [Hmsc::predictLatentFactor]
#' @name predictLF
#' @inheritParams LoadAs
#' @details The function is expected to be faster than the original function in
#'   the `Hmsc` package, especially when using TensorFlow for calculations and
#'   when working on parallel.
#'
#'   The main difference is that this function:
#' - allow for parallel processing (`NCores` argument);
#' - it is possible to use TensorFlow (`UseTF` argument) to make matrix
#'   calculations faster, particularly when used on GPU. The following modules
#'   are needed: `numpy`, `os`, `tensorflow`, `rdata`, `xarray`, and `pandas`.
#'   To use `TensorFlow`, the argument `EnvPath` should be set to the path of a
#'   Python environment with TensorFlow installed and `PythonScript` can
#'   optionally be set to the path of the Python script file containing the
#'   `crossprod_solve` function; otherwise, it will be loaded from the package;
#' - if `UseTF` is set to `FALSE`, the function uses R / CPP code in the
#'   calculations;
#' - the input `postEta` can be either a list or a file path for it, which can
#'   save memory, particularly when working on parallel. The function splits the
#'   `postEta` into smaller chunks, each for a combination of alpha value and
#'   latent factor and loop over them on parallel.
#' - calculates `D11` and `D12` matrices only once and save them to disk and
#'   call them when needed.

predictLF <- function(
    unitsPred, modelunits, postEta, postAlpha, rL, NCores = 8,
    TempDir = "TEMP2Pred", ModelName = NULL, UseTF = TRUE,
    PythonScript = NULL, EnvPath = NULL, use_single = FALSE,
    Path_postEtaPred = NULL, Return_postEtaPred = TRUE, nthreads = 5) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  IASDT.R::CatTime("Starting `PredictLF` function", Level = 1)

  # # ..................................................................... ###

  # Check inputs

  if (is.null(Path_postEtaPred) && isFALSE(Return_postEtaPred)) {
    stop(
      "`Return_postEtaPred` must be TRUE when `Path_postEtaPred` is NULL.",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Load postEta if it is a file path

  if (inherits(postEta, "character")) {
    IASDT.R::CatTime("Load postEta", Level = 1)
    if (!file.exists(postEta)) {
      stop(
        paste0(
          "The specified path for `postEta` does not exist. ",
          "Please verify the file path."),
        call. = FALSE)
    }
    postEta <- IASDT.R::LoadAs(postEta, nthreads = nthreads)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SampleID <- Unit_ID <- LF <- LF_ID <- etaPred <- Sample_IDs <- File <-
    crossprod_solve <- Alpha_ID <- NULL

  # # ..................................................................... ###

  # Check TensorFlow settings

  if (UseTF) {
    # Check if EnvPath directory exists
    if (is.null(EnvPath) || !dir.exists(EnvPath)) {
      stop(
        paste0(
          "When `UseTF` is TRUE, `EnvPath` must be specified and should point ",
          "to an existing directory with a Python environment"),
        call. = FALSE)
    }

    if (is.null(PythonScript)) {
      PythonScript <- system.file("crossprod_solve.py", package = "IASDT.R")
    }

    # Check if PythonScript exists
    if (!file.exists(PythonScript)) {
      stop(
        "Specified `PythonScript` does not exist at the provided path.",
        call. = FALSE)
    }

    # Suppress TensorFlow warnings and disable optimizations
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

    IASDT.R::CatTime("Computations will be made using TensorFlow", Level = 1)
  } else {
    IASDT.R::CatTime("Computations will be made using R/CPP", Level = 1)
  }

  # Extension for temporary files
  TempExt <- dplyr::if_else(UseTF, "rds", "qs")

  # # ..................................................................... ###

  # Adjust ModelName prefix

  if (is.null(ModelName)) {
    ModelName <- ""
  } else {
    ModelName <- paste0(ModelName, "_")
  }


  # Create a temporary directory to store intermediate results. This directory
  # will be used to save D11, D12, and intermediate postEta files, reducing
  # memory usage.
  fs::dir_create(TempDir)

  # indices of unitsPred in modelunits
  indOld <- (unitsPred %in% modelunits)
  # indices of new unitsPred
  indNew <- !(indOld)

  # # ..................................................................... ###

  # Calculate D11 and D12 only once

  IASDT.R::CatTime("Calculate/save D11 and D12 distance matrices", Level = 1)

  alphapw <- rL$alphapw
  s1 <- rL$s[modelunits, , drop = FALSE]
  s2 <- rL$s[unitsPred[indNew], , drop = FALSE]
  D11 <- Rfast::Dist(s1)
  D12 <- Rfast::dista(s1, s2)

  # Save D11 and D12 as qs/rds files
  Path_D11 <- file.path(TempDir, paste0(ModelName, "D11.", TempExt))
  Path_D12 <- file.path(TempDir, paste0(ModelName, "D12.", TempExt))
  if (UseTF) {
    saveRDS(D11, file = Path_D11)
    saveRDS(D12, file = Path_D12)
  } else {
    qs::qsave(D11, file = Path_D11, preset = "fast")
    qs::qsave(D12, file = Path_D12, preset = "fast")
  }

  # Clean up
  rm(rL, s1, s2, D11, D12)
  invisible(gc())

  if (!is.null(ModelName)) {
    on.exit(
      try({
        list.files(
          TempDir,
          pattern = paste0("^", ModelName, "_postEta"),
          full.names = TRUE) %>%
          c(Path_D11, Path_D12) %>%
          fs::file_delete()
      }, silent = TRUE), add = TRUE)
  } else {
    on.exit(
      try(fs::file_delete(c(Path_D11, Path_D12)), silent = TRUE),
      add = TRUE)
  }

  # # ..................................................................... ###

  # Convert postAlpha to tibble
  IASDT.R::CatTime(
    paste0("Prepare data for parallel processing using ", NCores, " cores"),
    Level = 1)
  postAlpha_tibble <- do.call(rbind, postAlpha) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    stats::setNames(paste0("LF_", seq_len(ncol(.))))
  rm(postAlpha)

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
        TempDir,
        paste0(ModelName, "postEta_ch", dplyr::row_number(), ".", TempExt)),
      Export = purrr::pmap(
        .l = list(SampleID, LF_ID, File),
        .f = function(SampleID, LF_ID, File) {
          Out <- postEta[SampleID] %>%
            purrr::map(~ .x[, LF_ID, drop = FALSE]) %>%
            simplify2array()
          if (UseTF) {
            saveRDS(Out, file = File)
          } else {
            qs::qsave(Out, file = File, preset = "fast")
          }
        }),
      Export = NULL)


  # # ..................................................................... ###

  # calculate etaPred for each unique alpha

  calc_eta_pred <- function(RowNum) {
    # Current denominator
    Denom <- Unique_Alpha$Denom[[RowNum]]
    # ID for latent factor
    LF_ID <- Unique_Alpha$LF_ID[[RowNum]]
    # ID for posterior sample
    SampleID <- Unique_Alpha$SampleID[[RowNum]]
    # File path for current alpha
    File <- Unique_Alpha$File[[RowNum]]

    # If the denominator is positive, perform calculations; otherwise, set
    # `eta_indNew` to zero.

    if (Denom > 0) {
      if (UseTF) {
        # Use TensorFlow

        # Activate the python environment
        reticulate::use_virtualenv(EnvPath, required = TRUE)
        # Source the script file containing the crossprod_solve function
        reticulate::source_python(PythonScript)

        eta_indNew <- crossprod_solve(
          Dist1 = Path_D11, Dist2 = Path_D12, Denom = Denom,
          List = File, use_single = use_single)
        eta_indNew <- purrr::map(
          .x = seq_along(eta_indNew),
          .f = ~ {
            tibble::tibble(etaPred = as.vector(eta_indNew[[.x]])) %>%
              dplyr::mutate(
                Unit_ID = unitsPred[indNew], SampleID = SampleID[.x])
          }) %>%
          dplyr::bind_rows()

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
              SampleID = SampleID[ID], etaPred = etaPred, Unit_ID = unitsPred)
          }) %>%
          dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
          dplyr::arrange(SampleID, Unit_ID, etaPred)
      }
    } else {

      # Handle cases where Denom is zero by setting `eta_indNew` to zero

      postEta0 <- IASDT.R::LoadAs(File, nthreads = nthreads)

      etaPred <- purrr::map_dfr(
        .x = seq_len(length(SampleID)),
        .f = function(ID) {
          eta <- postEta0[, , ID]
          etaPred <- rep(NA, length(unitsPred))
          etaPred[indOld] <- eta[match(unitsPred[indOld], modelunits)]
          etaPred[indNew] <- 0
          tibble::tibble(
            SampleID = SampleID[ID], etaPred = etaPred, Unit_ID = unitsPred)
        }) %>%
        dplyr::mutate(Unit_ID = factor(Unit_ID, levels = unitsPred)) %>%
        dplyr::arrange(SampleID, Unit_ID, etaPred)
    }

    # clean up
    fs::file_delete(File)

    return(etaPred)
  }

  # # ..................................................................... ###

  #
  IASDT.R::CatTime(paste0("Predicting Latent Factor in parallel"), Level = 1)

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2,
    future.gc = TRUE, future.seed = TRUE)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  future::plan("future::cluster", workers = c1, gc = TRUE)
  on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

  # Calculate etaPred
  etaPreds <- future.apply::future_lapply(
    X = seq_len(nrow(Unique_Alpha)),
    FUN = calc_eta_pred, future.seed = TRUE, future.chunk.size = 1,
    future.globals = c(
      "Unique_Alpha", "Path_D11", "Path_D12", "indNew", "unitsPred",
      "indOld", "modelunits", "EnvPath", "PythonScript", "UseTF", "use_single"),
    future.packages = c(
      "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble",
      "Matrix", "Hmsc", "qs", "fs", "purrr"))

  IASDT.R::CatTime("Merge results", Level = 1)
  # Merge results
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

  # # ..................................................................... ###

  # Save postEtaPred

  if (!is.null(Path_postEtaPred)) {
    IASDT.R::CatTime(
      paste0("Save postEtaPred to `", Path_postEtaPred, "`"), Level = 1)
    fs::dir_create(fs::path_dir(Path_postEtaPred))
    switch(
      fs::path_ext(Path_postEtaPred),
      "qs" = qs::qsave(postEtaPred, file = Path_postEtaPred, preset = "fast"),
      "RData" = save(postEtaPred, file = Path_postEtaPred),
      stop("Unsupported file extension in `Path_postEtaPred`.", call. = FALSE))
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "PredictLF was finished in ", Level = 1)

  # # ..................................................................... ###

  if (Return_postEtaPred) {
    return(postEtaPred)
  } else {
    return(Path_postEtaPred)
  }
}
