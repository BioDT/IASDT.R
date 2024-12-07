## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

## |------------------------------------------------------------------------| #
# Predict_Hmsc ----
## |------------------------------------------------------------------------| #

#' Predict_Hmsc
#'
#' Calculates predicted values from a fitted \code{Hmsc} model. This function
#' modifies the `Hmsc:::predict.Hmsc` function.
#'
#' @param Path_Model character string specifying a file name where the model
#'   object is saved.
#' @param XData a dataframe specifying the unpreprocessed covariates for the
#'   predictions to be made. Works only if the `XFormula` argument was specified
#'   in the [Hmsc::Hmsc] model constructor call. Requirements are similar to
#'   those in the `Hmsc` model constructor.
#' @param X a matrix specifying the covariates for the predictions to be made.
#'   Only one of `XData` and `X` arguments may be provided.
#' @param XRRRData a dataframe of covariates for reduced-rank regression.
#' @param XRRR a matrix of covariates for reduced-rank regression.
#' @param Gradient an object returned by [Hmsc::constructGradient]. Providing
#'   `Gradient` is an alternative for providing `XData`, `studyDesign` and
#'   `ranLevels`. Cannot be used together with `Yc`.
#' @param Yc a matrix of the outcomes that are assumed to be known for
#'   conditional predictions. Cannot be used together with `Gradient`.
#' @param mcmcStep the number of extra mcmc steps used for updating the random
#'   effects
#' @param expected boolean flag indicating whether to return the location
#'   parameter of the observation models or sample the values from those.
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param Temp_Cleanup logical, indicating whether to clean up temporary files.
#'   Defaults to `TRUE`.
#' @param Pred_Dir a character string specifying the directory where the
#'   predictions will be saved. Defaults to `NULL`, which saves model
#'   predictions to "Model_Prediction" folder of the current working directory.
#' @param Evaluate a logical flag indicating whether to evaluate the model
#'   predictions. Defaults to `FALSE`.
#' @param Eval_Name a character string specifying the name of the evaluation
#'   results. If `NULL`, the default name is used (`Eval_Model_Name.RData`).
#' @param Eval_Dir a character string specifying the directory where the
#'   evaluation results will be saved. Defaults to `Evaluation`.
#' @param RC a character string specifying the type of predictions to be made.
#'   If `NULL` (default), predictions are made for the latent factors. If `c`,
#'   predictions are made for response curves at mean coordinates. If `i`,
#'   predictions are made for response curves at infinite coordinates.
#' @param Pred_PA a matrix of presence-absence data for evaluation. If `NULL`
#'   (default), the presence-absence data from the model object is used. This
#'   argument is used only when `Evaluate` is `TRUE`.
#' @param Pred_XY a matrix of coordinates to be added to predicted values. If
#'   `NULL` (default), the coordinates from the model object is used.
#' @param LF_InputFile a character string specifying the file name where the
#'   latent factor predictions are saved. If `NULL` (default), latent factor
#'   predictions will be made. If specified, the latent factor predictions are
#'   read from the file. This allows to predicting the latent factors for new
#'   sites only once.
#' @param LF_Only a logical flag indicating whether to return the latent factor
#'   predictions only. Defaults to `FALSE`. This helps in predicting to new
#'   sites, allowing to predicting the latent factors only once, then the output
#'   can be loaded in other predictions when needed.
#' @param Verbose Logical. If TRUE, detailed output is printed. Default is
#'   `FALSE`.
#' @inheritParams Predict_LF
#' @name Predict_Hmsc
#' @export

Predict_Hmsc <- function(
    Path_Model,
    Loff = NULL, XData = NULL, X = NULL, XRRRData = NULL, XRRR = NULL,
    Gradient = NULL, Yc = NULL, mcmcStep = 1L, expected = TRUE, NCores = 8L,
    Model_Name = "Train", Temp_Dir = "TEMP_Pred", Temp_Cleanup = TRUE,
    RC = NULL, UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE,
    LF_OutFile = NULL, LF_Return = FALSE, LF_InputFile = NULL, LF_Only = FALSE,
    LF_NCores = NCores, LF_Check = FALSE, LF_Temp_Cleanup = TRUE,
    LF_Commands_Only = FALSE,
    Pred_Dir = NULL, Pred_PA = NULL, Pred_XY = NULL, Evaluate = FALSE,
    Eval_Name = NULL, Eval_Dir = "Evaluation", Verbose = TRUE) {

  # # ..................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  # To avoid non-standard evaluation
  Pred_XY <- Pred_XY
  Pred_PA <- Pred_PA

  if (LF_Only && is.null(LF_OutFile)) {
    stop(
      "`LF_OutFile` must be specified when `LF_Only` is `TRUE`",
      call. = FALSE)
  }

  if (!is.null(LF_OutFile) && !is.null(LF_InputFile)) {
    stop(
      "only one of `LF_OutFile` and `LF_InputFile` arguments can be specified",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chunk <- Sp <- IAS_ID <- Path_pred <- Sp_data <- data <- geometry <- x <-
    y <- SR_mean <- SR_sd <- SR_cov <- NULL

  # # ..................................................................... ###

  if (!is.null(RC) && !RC %in% c("c", "i")) {
    stop("`RC` must be either NULL or one of 'c' or 'i'", call. = FALSE)
  }

  if (is.null(Pred_Dir)) {
    Pred_Dir <- file.path(dirname(dirname(Path_Model)), "Model_Prediction")
  }

  if (is.null(RC) || RC == "c") {
    IASDT.R::CatTime("Creating/checking output paths")
    fs::dir_create(c(Temp_Dir, Pred_Dir))
    if (Evaluate) {
      fs::dir_create(Eval_Dir)
    }
  }

  # # ..................................................................... ###

  # Load model if it is a character
  if (inherits(Path_Model, "character")) {
    IASDT.R::CatTime("Load model object")
    Model <- IASDT.R::LoadAs(Path_Model)
  }

  # # ..................................................................... ###

  # Combines a list of single or several MCMC chains into a single chain
  IASDT.R::CatTime("Combine list of posteriors")
  post <- Hmsc::poolMcmcChains(Model$postList)
  studyDesign <- Model$studyDesign
  ranLevels <- Model$ranLevels

  # # ..................................................................... ###

  IASDT.R::CatTime("Check input parameters")

  if (!is.null(Gradient)) {
    if (!is.null(Yc)) {
      stop(
        paste0(
          "predict with arguments 'Yc' and 'Gradient' jointly is not ",
          "implemented (yet)"),
        call. = FALSE)
    }
    XData <- Gradient$XDataNew
    studyDesign <- Gradient$studyDesignNew
    ranLevels <- Gradient$rLNew
  }

  if (!is.null(XData) && !is.null(X)) {
    stop(
      "only one of XData and X arguments can be specified",
      call. = FALSE)
  }
  if (!is.null(XRRRData) && !is.null(XRRR)) {
    stop(
      "only one of XRRRData and XRRR arguments can be specified",
      call. = FALSE)
  }

  if (!is.null(XData)) {
    switch(
      class(XData)[1L],
      list = {
        if (any(unlist(lapply(XData, is.na)))) {
          stop("NA values are not allowed in 'XData'", call. = FALSE)
        }
        xlev <- lapply(Reduce(rbind, Model$XData), levels)
        xlev <- xlev[unlist(lapply(Reduce(rbind, Model$XData), is.factor))]
        X <- lapply(
          XData, function(a) {
            stats::model.matrix(Model$XFormula, a, xlev = xlev)
          })
      },
      data.frame = {
        if (any(is.na(XData))) {
          stop("NA values are not allowed in 'XData'", call. = FALSE)
        }
        xlev <- lapply(Model$XData, levels)
        xlev <- xlev[unlist(lapply(Model$XData, is.factor))]
        X <- stats::model.matrix(Model$XFormula, XData, xlev = xlev)
      })

  } else {
    if (is.null(X)) {
      X <- Model$X
    }
  }

  if (!is.null(XRRRData)) {
    xlev <- lapply(Model$XRRRData, levels)
    xlev <- xlev[unlist(lapply(Model$XRRRData, is.factor))]
    XRRR <- stats::model.matrix(Model$XRRRFormula, XRRRData, xlev = xlev)
  } else {
    if (is.null(Model$ncRRR)) {
      Model$ncRRR <- 0
    }
    if (is.null(XRRR) && Model$ncRRR > 0) {
      XRRR <- Model$XRRR
    }
  }

  switch(
    class(X)[1L],
    list = {
      nyNew <- nrow(X[[1]])
    },
    matrix = {
      nyNew <- nrow(X)
    })


  if (!is.null(Yc)) {
    if (ncol(Yc) != Model$ns) {
      stop("number of columns in Yc must be equal to ns", call. = FALSE)
    }
    if (nrow(Yc) != nyNew) {
      stop("number of rows in Yc and X must be equal", call. = FALSE)
    }
  }

  if (!is.null(Loff)) {
    if (ncol(Loff) != Model$ns) {
      stop("number of columns in Loff must be equal to ns", call. = FALSE)
    }
    if (nrow(Loff) != nyNew) {
      stop("number of rows in Loff and X must be equal", call. = FALSE)
    }
  }

  if (!all(Model$rLNames %in% colnames(studyDesign))) {
    stop(
      "dfPiNew does not contain all the necessary named columns",
      call. = FALSE)
  }

  if (!all(Model$rLNames %in% names(ranLevels))) {
    stop("rL does not contain all the necessary named levels", call. = FALSE)
  }

  if (!is.null(studyDesign)) {
    dfPiNew <- studyDesign[, Model$rLNames, drop = FALSE]
  } else {
    dfPiNew <- matrix(NA, nyNew, 0)
  }
  rL <- ranLevels[Model$rLNames]


  if (!is.null(Yc)) {
    ## object can have pre-computed data parameters, but not
    ## necessarily. These are needed only in updateEta(), but get it
    ## here anyway...
    if (is.null(Model$rLPar)) {
      rLPar <- Hmsc::computeDataParameters(Model)$rLPar
    } else {
      rLPar <- Model$rLPar
    }
  } else {
    rLPar <- NULL
  }

  predN <- length(post)

  # # ..................................................................... ###

  # free some memory
  IASDT.R::CatTime("Free some memory")
  Model$postList <- Model$YScaled <- Model$X <- Model$XScaled <- NULL

  Mod_nr <- Model$nr
  Mod_dfPi <- Model$dfPi

  # Save smaller version of the model object for later use
  Model_File_small <- file.path(Temp_Dir, "Model_small.qs2")
  if (!file.exists(Model_File_small)) {
    IASDT.R::CatTime(
      "Save smaller version of the model object to disk",
      Level = 1)
    IASDT.R::SaveAs(InObj = Model, OutPath = Model_File_small)
  }
  rm(Model, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  IASDT.R::CatTime("Predict Latent Factor")

  predPostEta <- vector("list", Mod_nr)
  PiNew <- matrix(NA, nrow(dfPiNew), Mod_nr)

  # Whether to use `Predict_LF` or read its results from file
  if (!is.null(Mod_nr)) {
    if (is.null(LF_InputFile) || length(LF_InputFile) != Mod_nr) {
      CalcLF <- TRUE
    } else {
      CalcLF <- FALSE
    }
  }

  # Do not use `Predict_LF` when predicting values for response curves when
  # using coordinates = "i" in constructGradient
  if (!is.null(RC)) {
    if (RC == "i") {
      for (r in seq_len(Mod_nr)) {
        if (r == 1) {
          IASDT.R::CatTime(
            "LF prediction for response curve with infinite coordinates",
            Level = 1)
        }
        nLF <- length(post[[1]]$Alpha[[1]])
        predPostEta[[r]] <- replicate(
          n = predN,
          expr = structure(
            rep(0, nLF), dim = c(1L, nLF),
            dimnames = list("new_unit", NULL)),
          simplify = FALSE)
      }
      CalcLF <- FALSE
    }
  }

  # Calculate latent factors
  if (CalcLF) {

    for (r in seq_len(Mod_nr)) {

      postAlpha <- lapply(post, function(c) c$Alpha[[r]])

      # Save postEta to file and load it from Predict_LF. This helps to avoid
      # the unnecessary copying of the postEta object to all cores
      postEta_file <- file.path(
        Temp_Dir, paste0(Model_Name, "_r", r, "_postEta.qs2"))

      if (isFALSE(IASDT.R::CheckData(postEta_file, warning = FALSE))) {
        IASDT.R::CatTime("Save postEta to file", Level = 1)
        postEta <- lapply(post, function(c) c$Eta[[r]])
        IASDT.R::SaveAs(InObj = postEta, OutPath = postEta_file)
        rm(postEta, envir = environment())
        invisible(gc())
      }

      # Save post to file and load it later
      if (r == Mod_nr) {
        post_file <- file.path(Temp_Dir, paste0(Model_Name, "_post.qs2"))

        if (isFALSE(IASDT.R::CheckData(post_file, warning = FALSE))) {
          # free some memory
          post <- lapply(post, function(x) {
            x$Eta <- x$Psi <- x$V <- x$Delta <- x$Gamma <- x$rho <- NULL
            x
          })
          IASDT.R::CatTime("Save post to file", Level = 1)
          IASDT.R::SaveAs(InObj = post, OutPath = post_file)
        }
        rm(post, envir = environment())
        invisible(gc())
      }

      if (r == 1) {
        IASDT.R::CatTime("LF prediction using `Predict_LF`", Level = 1)
      }

      predPostEta[[r]] <- IASDT.R::Predict_LF(
        unitsPred = levels(dfPiNew[, r]),
        modelunits = levels(Mod_dfPi[, r]),
        postEta = postEta_file, postAlpha = postAlpha, LF_rL = rL[[r]],
        LF_NCores = LF_NCores, LF_Temp_Cleanup = LF_Temp_Cleanup,
        LF_OutFile = LF_OutFile, LF_Return = LF_Return, LF_Check = LF_Check,
        LF_Commands_Only = LF_Commands_Only,
        Temp_Dir = Temp_Dir, Model_Name = Model_Name,
        UseTF = UseTF, TF_Environ = TF_Environ, TF_use_single = TF_use_single)

      rm(postEta_file, envir = environment())

      if (LF_Commands_Only) {
        return(invisible(NULL))
      }

      rowNames <- rownames(predPostEta[[r]][[1]])
      PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
    }

  } else {

    if (is.null(RC) || RC == "c") {
      IASDT.R::CatTime(
        paste0("Loading LF prediction from disk: `", LF_InputFile, "`"),
        Level = 1)

      for (r in seq_len(Mod_nr)) {
        predPostEta[[r]] <- IASDT.R::LoadAs(LF_InputFile[[r]])
        rowNames <- rownames(predPostEta[[r]][[1]])
        PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
      }
    }
  }

  if (Mod_nr > 0) {
    ppEta <- simplify2array(predPostEta)
  } else {
    ppEta <- matrix(list(), predN, 0)
  }

  # Only predicting latent factor for new sites or for response curves at
  # median site
  if (LF_Only) {
    return(LF_OutFile)
  }

  # free some memory
  rm(predPostEta, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  IASDT.R::CatTime("Predicting")

  if (!exists("post")) {
    IASDT.R::CatTime("Loading post from disk", Level = 1)
    post <- IASDT.R::LoadAs(post_file)
  }

  # Read model object from disk
  Model <- IASDT.R::LoadAs(Model_File_small)

  # prediction data for response curves
  if (!is.null(RC)) {
    IASDT.R::CatTime(
      "Predicting data for response curve (sequentially)",
      Level = 1)

    preds <- lapply(
      seq_len(predN),
      function(pN, ...) {
        get1prediction(
          Model, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
          ppEta[pN, ], PiNew, dfPiNew, nyNew, expected, mcmcStep)
      })

    IASDT.R::CatDiff(
      InitTime = .StartTime, Prefix = "Prediction was finished in ", Level = 1)
    return(preds)
  }

  # Save ppEta / post as small chunks
  IASDT.R::CatTime("Save ppEta / post as small chunks", Level = 1)
  chunk_size <- 25
  ChunkIDs <- ceiling(seq_along(post) / chunk_size)
  Chunks <- purrr::map_chr(
    .x = seq_len(max(ChunkIDs)),
    .f = ~ {
      IDs <- which(ChunkIDs == .x)
      Ch <- list(ppEta = ppEta[IDs], post = post[IDs])
      ChunkFile <- file.path(
        Temp_Dir, paste0(Model_Name, "_preds_ch", .x, ".qs2"))
      IASDT.R::SaveAs(InObj = Ch, OutPath = ChunkFile)
      return(ChunkFile)
    })

  rm(ChunkIDs, post, ppEta, envir = environment())
  invisible(gc())

  IASDT.R::CatTime(
    paste0(
      "Preparing working on parallel using ",
      min(NCores, length(Chunks)), " cores"),
    Level = 1)

  seeds <- sample.int(.Machine$integer.max, predN)

  c1 <- parallel::makePSOCKcluster(min(NCores, length(Chunks)))
  on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

  parallel::clusterExport(
    cl = c1,
    varlist = c(
      "Model", "X", "XRRR", "Yc", "Loff", "rL", "rLPar", "PiNew",
      "dfPiNew", "nyNew", "expected", "mcmcStep", "seeds", "chunk_size",
      "Chunks", "Temp_Dir", "Model_Name", "get1prediction"),
    envir = environment())

  invisible(parallel::clusterEvalQ(
    cl = c1,
    expr = {
      sapply(
        c(
          "dplyr", "Rcpp", "RcppArmadillo", "Matrix", "float", "qs2", "Hmsc",
          "purrr", "tibble", "Hmsc", "Rfast", "caret", "pROC", "ecospat", "sf"),
        library, character.only = TRUE)
    }))

  IASDT.R::CatTime("Making predictions on parallel", Level = 1)

  pred <- parallel::parLapplyLB(
    cl = c1,
    X = seq_len(length(Chunks)),
    fun = function(Chunk) {

      ChunkFile <- Chunks[Chunk]
      Ch <- IASDT.R::LoadAs(ChunkFile)
      ppEta <- Ch$ppEta
      post <- Ch$post
      rm(Ch, envir = environment())
      Seed <- (Chunk - 1) * chunk_size
      Seed <- seq(Seed + 1, Seed + chunk_size)
      Seed <- seeds[Seed]

      PredChunk <- purrr::map(
        .x = seq_len(chunk_size),
        .f = function(pN) {
          get1prediction(
            object = Model, X = X, XRRR = XRRR, Yc = Yc,
            Loff = Loff, rL = rL, rLPar = rLPar, sam = post[[pN]],
            predPostEta = ppEta[pN], PiNew = PiNew,
            dfPiNew = dfPiNew, nyNew = nyNew, expected = expected,
            mcmcStep = mcmcStep, seed = Seed[pN])
        })

      # Species richness
      ChunkSR <- simplify2array(lapply(X = PredChunk, FUN = rowSums)) %>%
        float::fl()
      dimnames(ChunkSR) <- NULL
      ChunkSR_File <- file.path(
        Temp_Dir, paste0("Pred_", Model_Name, "_ch", Chunk, "_SR.qs2"))
      IASDT.R::SaveAs(InObj = ChunkSR, OutPath = ChunkSR_File)
      
      rm(ChunkSR, envir = environment())

      save_qs2 <- function(object, file, max_retries = 5) {
            attempt <- 1
            while (attempt <= max_retries) {
              tryCatch(
                {
                  # Attempt to save the file
                  qs2::qs_save(object = object, file = file, nthreads = 5L)
                  
                  # Attempt to read the file back to verify it
                  read_back <- qs2::qs_read(file)
                  
                  # If the read is successful, return TRUE
                  return(TRUE)
                },
                error = function(e) {
                  if (attempt == max_retries) {
                    stop("Failed to save and verify file after ", max_retries, " attempts: ", file, "\nError: ", e$message)
                  }
                }
              )
              attempt <- attempt + 1
            }
            return(FALSE)
          }

      # Species predictions
      ChunkSp <- purrr::map_dfr(
        .x = seq_len(length(Model$spNames)),
        .f = function(Sp) {

          SpD <- purrr::map(PredChunk, ~ .x[, Sp], ncol = 1) %>%
            simplify2array() %>%
            float::fl()
          dimnames(SpD) <- NULL

          ChunkSp_File <- file.path(
            Temp_Dir,
            paste0("Pred_", Model_Name, "_ch", Chunk, "_taxon", Sp, ".qs2"))

          # Saving to qs2 failed for many files, breaking the loop and 
          # function. Here is a workaround to ensure the file is saved and 
          # verified
          # IASDT.R::SaveAs(InObj = SpD, OutPath = ChunkSp_File)
          # qs2::qs_save(object = SpD, file = ChunkSp_File, nthreads = 5L)
          save_qs2(object = SpD, file = ChunkSp_File)

          cbind.data.frame(
            Chunk = Chunk, Sp = Sp, IAS_ID = Model$spNames[Sp],
            ChunkSp_File = ChunkSp_File) %>%
            return()
        }) %>%
        dplyr::bind_rows(
          tibble::tibble(
            Chunk = Chunk, Sp = 0, IAS_ID = "SR", ChunkSp_File = ChunkSR_File),
          .)

      rm(PredChunk, envir = environment())

      invisible(gc())
      return(ChunkSp)
    })

  pred <- tibble::tibble(dplyr::bind_rows(pred))

  invisible(parallel::clusterEvalQ(cl = c1, expr = invisible(gc())))

  # # ..................................................................... ###

  IASDT.R::CatTime("Summarizing prediction outputs / Evaluation", Level = 1)

  Eval_DT <- dplyr::select(pred, -Chunk) %>%
    dplyr::group_nest(Sp, IAS_ID) %>%
    dplyr::mutate(data = purrr::map(data, unlist))

  parallel::clusterExport(
    cl = c1,
    varlist = c(
      "Eval_DT", "Evaluate", "Pred_Dir", "Model_Name", "Pred_PA", "Pred_XY"),
    envir = environment())

  Eval_DT <- parallel::parLapplyLB(
    cl = c1,
    X = seq_len(nrow(Eval_DT)),
    fun = function(ID) {

      Sp <- Eval_DT$Sp[[ID]]
      if (Sp == 0) {
        Sp2 <- "SR"
      } else {
        Sp2 <- paste0("taxon", Sp)
      }
      IAS_ID <- Eval_DT$IAS_ID[[ID]]
      data <- as.vector(Eval_DT$data[[ID]])

      # SpDT <- purrr::map(data, qs2::qs_read) %>%
      #   do.call(cbind, .) %>%
      #   as.double()

      SpDT <- purrr::map(data, function(x) {
        attempt <- 1
        max_retries <- 5
        while (attempt <= max_retries) {
          result <- tryCatch(
            qs2::qs_read(x),
            error = function(e) {
              if (attempt == max_retries) {
                stop(
                  paste0(
                    "Failed to read ", x, " after ", max_retries,
                    " attempts: ", x, "\nError: ", e$message), call. = FALSE)
              }
              return(NULL)
            })
          if (!is.null(result)) {
            return(result)
          }
          attempt <- attempt + 1
        }
      }) %>%
        do.call(cbind, .) %>%
        as.double()

      # Mean prediction
      SpDT_Mean <- Rfast::rowmeans(SpDT)

      # standard deviation of prediction
      SpDT_SD <- Rfast::rowVars(SpDT, std = TRUE)

      # Coefficient of variation
      SpDT_Mean0 <- SpDT_Mean
      # Replace very small mean values with reasonably small number to avoid
      # overflow warning
      SpDT_Mean0[SpDT_Mean0 < 1e-8] <- 1e-8
      SpDT_Cov <- SpDT_SD / SpDT_Mean0

      rm(SpDT, envir = environment())
      invisible(gc())

      if (is.null(Pred_XY)) {
        Pred_XY <- Model$rL$sample$s
      }

      PredSummary <- tibble::tibble(
        as.data.frame(Pred_XY),
        Mean = SpDT_Mean, SD = SpDT_SD, Cov = SpDT_Cov) %>%
        stats::setNames(
          c(
            "x", "y", paste0(IAS_ID, "_mean"),
            paste0(IAS_ID, "_sd"), paste0(IAS_ID, "_cov"))) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

      PredSummaryFile <- file.path(
        Pred_Dir, paste0("Pred_", Model_Name, "_", Sp2, ".qs2"))

      IASDT.R::SaveAs(InObj = PredSummary, OutPath = PredSummaryFile)

      if (Evaluate && Sp2 != "SR") {
        if (is.null(Pred_PA)) {
          PresAbs <- Model$Y[, Sp]
        } else {
          PresAbs <- Pred_PA[, Sp]
        }

        if (length(unique(PresAbs)) == 2) {
          # Calculate evaluation metrics if there are two both presence and
          # absence info in (testing) data
          RMSE <- caret::RMSE(PresAbs, SpDT_Mean)
          MeanPres <- mean(SpDT_Mean[which(PresAbs == 1)])
          MeanAbs <- mean(SpDT_Mean[which(PresAbs == 0)])
          TjurR2 <- MeanPres - MeanAbs
          AUC <- pROC::auc(
            response = PresAbs, predictor = SpDT_Mean,
            levels = c(0, 1), direction = "<") %>%
            as.numeric()
          Boyce <- ecospat::ecospat.boyce(
            fit = SpDT_Mean, obs = SpDT_Mean[PresAbs == 1],
            PEplot = FALSE)$cor
        } else {
          RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
        }
      } else {
        RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
      }

      tibble::tibble(
        Sp = Sp, IAS_ID = IAS_ID, Path_pred = PredSummaryFile,
        RMSE = RMSE, AUC = AUC, Boyce = Boyce, TjurR2 = TjurR2) %>%
        return()

    })

  parallel::stopCluster(c1)
  invisible(gc())

  # # ..................................................................... ###

  # Save predictions for all species in a single file
  IASDT.R::CatTime("Save predictions for all species in a single file")

  Eval_DT <- dplyr::bind_rows(Eval_DT)

  Predictions <- dplyr::select(Eval_DT, Path_pred, Sp, IAS_ID) %>%
    dplyr::mutate(
      Sp_data = purrr::map(
        .x = Path_pred,
        .f = ~ {
          IASDT.R::LoadAs(.x) %>%
            tidyr::pivot_longer(
              cols = tidyselect::starts_with(c("Sp_", "SR_")),
              names_to = "Species", values_to = "Prediction")
        })) %>%
    dplyr::pull(Sp_data) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_wider(names_from = "Species", values_from = "Prediction") %>%
    dplyr::relocate(gtools::mixedsort(names(.))) %>%
    dplyr::select(
      x, y, geometry, SR_mean, SR_sd, SR_cov, tidyselect::everything()) %>%
    dplyr::mutate(Model_Name = Model_Name, .before = "SR_mean")

  Pred_File <- file.path(
    Pred_Dir,
    paste0(
      "Prediction_",
      stringr::str_remove(Model_Name, "_Clamping|_NoClamping"), ".qs2"))

  IASDT.R::SaveAs(InObj = Predictions, OutPath = Pred_File)
  
  if (Temp_Cleanup) {
    try(fs::file_delete(Eval_DT$Path_pred), silent = TRUE)
  }
  IASDT.R::CatTime(
    paste0("Predictions were saved to `", Pred_File, "`"), Level = 1)

  if (Evaluate) {
    if (is.null(Eval_Name)) {
      Eval_Path <- file.path(
        Eval_Dir,
        paste0(
          "Eval_", stringr::str_remove(Model_Name, "_Train|_Current"), ".qs2"))
    } else {
      Eval_Path <- file.path(
        Eval_Dir,
        paste0(
          "Eval_", stringr::str_remove(Model_Name, "_Train|_Current"),
          "_", Eval_Name, ".qs2"))
    }

    Eval_DT <- dplyr::select(Eval_DT, -Path_pred)
    IASDT.R::SaveAs(InObj = Eval_DT, OutPath = Eval_Path)

    IASDT.R::CatTime(
      paste0(
        "Evaluation results were saved to `",
        file.path(Eval_Dir, "Eval_DT.qs2"), "`"),
      Level = 1)
  } else {
    Eval_Path <- NULL
  }

  # # ..................................................................... ###

  if (is.null(LF_OutFile)) {
    LF_Path <- LF_InputFile
  } else {
    LF_Path <- LF_OutFile
  }

  # # ..................................................................... ###

  # Clean up
  if (Temp_Cleanup) {

    IASDT.R::CatTime("Cleaning up temporary files", Level = 1)

    try(
      {
        Pattern <- paste0("(Pred_){0,}", Model_Name, ".+qs2")
        file_paths <- list.files(
          path = normalizePath(Temp_Dir, winslash = "/"),
          pattern = Pattern, full.names = TRUE)
        fs::file_delete(file_paths)
      },
      silent = TRUE)

    try(fs::file_delete(Model_File_small), silent = TRUE)

  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Prediction was finished in ")

  tibble::tibble(
    Pred_Path = Pred_File,  Eval_Path = Eval_Path,
    LF_Path = LF_Path) %>%
    return()
}

## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

## |------------------------------------------------------------------------| #
# get1prediction ----
## |------------------------------------------------------------------------| #

#' get1prediction
#'
## internal function to get one prediction
##
##  Needs following variables or arguments that must be passed:
##  PiNew X XRRR Yc dfPiNew expected mcmcStep nyNew object pN post
##  predPostEta rL rLPar
#' @noRd

get1prediction <- function(
    object, X, XRRR, Yc, Loff, rL, rLPar, sam, predPostEta, PiNew, dfPiNew,
    nyNew, expected, mcmcStep, seed = NULL) {

  updateZ <- updateEta <- NULL


  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (object$ncRRR > 0) {
    XB <- XRRR %*% t(sam$wRRR)
  }

  switch(
    class(X)[1L],
    matrix = {
      X1 <- X
      if (object$ncRRR > 0) {
        X1 <- cbind(X1, XB)
      }
      LFix <- X1 %*% sam$Beta
    },
    list = {
      LFix <- matrix(NA, nyNew, object$ns)
      for (j in 1:object$ns) {
        X1 <- X[[j]]
        if (object$ncRRR > 0) {
          X1 <- cbind(X1, XB)
        }
        LFix[, j] <- X1 %*% sam$Beta[, j]
      }
    }
  )

  Mod_nr <- object$nr
  LRan <- vector("list", Mod_nr)
  Eta <- vector("list", Mod_nr)

  for (r in seq_len(Mod_nr)) {
    Eta[[r]] <- predPostEta[[r]]

    if (rL[[r]]$xDim == 0) {
      LRan[[r]] <- Eta[[r]][as.character(dfPiNew[, r]), ] %*% sam$Lambda[[r]]
    } else {
      LRan[[r]] <- matrix(0, object$ny, object$ns)
      for (k in 1:rL[[r]]$xDim) {
        LRan[[r]] <- LRan[[r]] +
          (Eta[[r]][as.character(dfPiNew[, r]), ] *
             rL[[r]]$x[as.character(dfPiNew[, r]), k]) %*%
          sam$Lambda[[r]][, , k]
      }
    }
  }

  L <- Reduce("+", c(list(LFix), LRan))
  if (!is.null(Loff)) L <- L + Loff

  ## predict can be slow with Yc and especially with high mcmcStep
  if (!is.null(Yc) && any(!is.na(Yc))) {
    Z <- L
    Z <- updateZ(
      Y = Yc, Z = Z, Beta = sam$Beta, iSigma = 1 / sam$sigma, Eta = Eta,
      Lambda = sam$Lambda, Loff = Loff, X = X, Pi = PiNew, dfPi = dfPiNew,
      distr = object$distr, rL = rL
    )
    ## species CV from computePredictedValues runs this innermost
    ## loop nfolds * nfolds.sp * predN * mcmcStep times
    for (sN in seq_len(mcmcStep)) {
      Eta <- updateEta(
        Y = Yc, Z = Z, Beta = sam$Beta, iSigma = 1 / sam$sigma,
        Eta = Eta, Lambda = sam$Lambda, Alpha = sam$Alpha,
        rLPar = rLPar, Loff = Loff, X = X, Pi = PiNew, dfPi = dfPiNew, rL = rL
      )
      Z <- updateZ(
        Y = Yc, Z = Z, Beta = sam$Beta, iSigma = 1 / sam$sigma, Eta = Eta,
        Lambda = sam$Lambda, Loff = Loff, X = X, Pi = PiNew, dfPi = dfPiNew,
        distr = object$distr, rL = rL
      )
    }
    for (r in seq_len(Mod_nr)) {
      if (rL[[r]]$xDim == 0) {
        LRan[[r]] <- Eta[[r]][as.character(dfPiNew[, r]), ] %*%
          sam$Lambda[[r]]
      } else {
        LRan[[r]] <- matrix(0, object$ny, object$ns)
        for (k in 1:rL[[r]]$xDim) {
          LRan[[r]] <- LRan[[r]] +
            (Eta[[r]][as.character(dfPiNew[, r]), ] *
               rL[[r]]$x[as.character(dfPiNew[, r]), k]) %*%
            sam$Lambda[[r]][, , k]
        }
      }
    }
    L <- Reduce("+", c(list(LFix), LRan))
  }

  if (!expected) {
    Z <- L + matrix(sqrt(sam$sigma), nrow(L), object$ns, byrow = TRUE) *
      matrix(stats::rnorm(nrow(L) * object$ns), nrow(L), object$ns)
  } else {
    Z <- L
  }

  for (j in 1:object$ns) {
    if (object$distr[j, "family"] == 2) { # probit
      if (expected) {
        Z[, j] <- IASDT.R::fast_pnorm(Z[, j])
      } else {
        Z[, j] <- as.numeric(Z[, j] > 0)
      }
    }
    if (object$distr[j, "family"] == 3) { # poisson
      if (expected) {
        Z[, j] <- exp(Z[, j] + sam$sigma[j] / 2)
      } else {
        Z[, j] <- stats::rpois(nrow(Z), exp(Z[, j]))
      }
    }
  }
  colnames(Z) <- object$spNames

  for (i in 1:object$ns) {
    m <- object$YScalePar[1, i]
    s <- object$YScalePar[2, i]
    if (m != 0 || s != 1) {
      Z[, i] <- Z[, i] * s + m
    }
  }

  Z
}
