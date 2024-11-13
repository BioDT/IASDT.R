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
#' @param object a fitted `Hmsc` model object or a character string specifying a
#'   file name where the model object is saved.
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
#' @param Pred_Dir a character string specifying the directory where the
#'   predictions will be saved.
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
#'   latent factor predictions are saved. If `NULL`, the predictions are saved
#'   to a temporary file. This argument is used only when `RC` is `NULL`.
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
    object,
    Loff = NULL, XData = NULL, X = NULL, XRRRData = NULL, XRRR = NULL,
    Gradient = NULL, Yc = NULL, mcmcStep = 1, expected = TRUE, NCores = 8,
    Model_Name = "Train", Temp_Dir = "TEMP2Pred", RC = NULL,
    UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE, LF_OutFile = NULL,
    LF_Return = TRUE, LF_InputFile = NULL, LF_Only = FALSE,
    Pred_Dir = "Model_Predictions", Pred_PA = NULL, Pred_XY = NULL,
    Evaluate = FALSE, Eval_Name = NULL, Eval_Dir = "Evaluation",
    Verbose = TRUE) {

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

  if (is.null(RC) || RC == "c") {
    IASDT.R::CatTime("Creating/checking output paths")
    fs::dir_create(c(Temp_Dir, Pred_Dir))
    if (Evaluate) {
      fs::dir_create(Eval_Dir)
    }
  }

  # # ..................................................................... ###

  # Load model if it is a character
  if (inherits(object, "character")) {
    IASDT.R::CatTime("Load model object")
    object <- IASDT.R::LoadAs(object)
  }

  # # ..................................................................... ###

  # Combines a list of single or several MCMC chains into a single chain
  IASDT.R::CatTime("Combines list of posteriors")
  post <- Hmsc::poolMcmcChains(object$postList)
  studyDesign <- object$studyDesign
  ranLevels <- object$ranLevels

  # # ..................................................................... ###

  # Clean up
  on.exit({
    try(
      fs::file_delete(
        list.files(
          Temp_Dir,
          pattern = paste0("^", Model_Name, "_.+"), full.names = TRUE)),
      silent = TRUE
    )},
    add = TRUE)

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
        xlev <- lapply(Reduce(rbind, object$XData), levels)
        xlev <- xlev[unlist(lapply(Reduce(rbind, object$XData), is.factor))]
        X <- lapply(
          XData, function(a) {
            stats::model.matrix(object$XFormula, a, xlev = xlev)
          })
      },
      data.frame = {
        if (any(is.na(XData))) {
          stop("NA values are not allowed in 'XData'", call. = FALSE)
        }
        xlev <- lapply(object$XData, levels)
        xlev <- xlev[unlist(lapply(object$XData, is.factor))]
        X <- stats::model.matrix(object$XFormula, XData, xlev = xlev)
      })
  } else {
    if (is.null(X)) {
      X <- object$X
    }
  }

  if (!is.null(XRRRData)) {
    xlev <- lapply(object$XRRRData, levels)
    xlev <- xlev[unlist(lapply(object$XRRRData, is.factor))]
    XRRR <- stats::model.matrix(object$XRRRFormula, XRRRData, xlev = xlev)
  } else {
    if (is.null(object$ncRRR)) {
      object$ncRRR <- 0
    }
    if (is.null(XRRR) && object$ncRRR > 0) {
      XRRR <- object$XRRR
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
    if (ncol(Yc) != object$ns) {
      stop("number of columns in Yc must be equal to ns", call. = FALSE)
    }
    if (nrow(Yc) != nyNew) {
      stop("number of rows in Yc and X must be equal", call. = FALSE)
    }
  }

  if (!is.null(Loff)) {
    if (ncol(Loff) != object$ns) {
      stop("number of columns in Loff must be equal to ns", call. = FALSE)
    }
    if (nrow(Loff) != nyNew) {
      stop("number of rows in Loff and X must be equal", call. = FALSE)
    }
  }

  if (!all(object$rLNames %in% colnames(studyDesign))) {
    stop(
      "dfPiNew does not contain all the necessary named columns", call. = FALSE)
  }

  if (!all(object$rLNames %in% names(ranLevels))) {
    stop("rL does not contain all the necessary named levels", call. = FALSE)
  }

  if (!is.null(studyDesign)) {
    dfPiNew <- studyDesign[, object$rLNames, drop = FALSE]
  } else {
    dfPiNew <- matrix(NA, nyNew, 0)
  }
  rL <- ranLevels[object$rLNames]


  if (!is.null(Yc)) {
    ## object can have pre-computed data parameters, but not
    ## necessarily. These are needed only in updateEta(), but get it
    ## here anyway...
    if (is.null(object$rLPar)) {
      rLPar <- Hmsc::computeDataParameters(object)$rLPar
    } else {
      rLPar <- object$rLPar
    }
  } else {
    rLPar <- NULL
  }

  predN <- length(post)

  # # ..................................................................... ###

  # free some memory
  IASDT.R::CatTime("Free some memory")
  object$postList <- object$YScaled <- object$X <- object$XScaled <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime("Predict Latent Factor")

  predPostEta <- vector("list", object$nr)
  PiNew <- matrix(NA, nrow(dfPiNew), object$nr)

  # Whether to use predictLatentFactor or read its results from file
  if (!is.null(object$nr)) {
    if (is.null(LF_InputFile) || length(LF_InputFile) != object$nr) {
      CalcLF <- TRUE
    } else {
      CalcLF <- FALSE
    }
  }


  # Do not use `predictLatentFactor` when predicting values for response
  # curves when using coordinates = "i" in constructGradient
  if (!is.null(RC)) {
    if (RC == "i") {
      for (r in seq_len(object$nr)) {
        if (r == 1) {
          IASDT.R::CatTime(
            "LF prediction for response curve with infinite coordinates",
            Level = 1)
        }
        nLF <- length(post[[1]]$Alpha[[1]])
        predPostEta[[r]] <- replicate(
          n = predN,
          expr = structure(
            rep(0, nLF), dim = c(1L, nLF), dimnames = list("new_unit", NULL)),
          simplify = FALSE)
      }
      CalcLF <- FALSE
    }
  }


  if (CalcLF) {

    for (r in seq_len(object$nr)) {
      postEta <- lapply(post, function(c) c$Eta[[r]])
      postAlpha <- lapply(post, function(c) c$Alpha[[r]])

      if (r == object$nr) {
        # free some memory
        post <- lapply(post, function(x) {
          x$Eta <- x$Psi <- x$V <- x$Delta <- x$Gamma <- x$rho <- NULL
          x
        })

        # Save post to file and load it later
        IASDT.R::CatTime("Save post to file and load it later", Level = 1)
        post_file <- file.path(Temp_Dir, paste0(Model_Name, "_post.qs"))
        qs::qsave(post, file = post_file, preset = "fast")
        rm(post)
        invisible(gc())
      }

      # Save postEta to file and load it from predictLatentFactor. This helps
      # to avoid the unnecessary copying of the postEta object to all cores
      IASDT.R::CatTime("Save postEta to file", Level = 1)
      postEta_file <- file.path(Temp_Dir, paste0(Model_Name, "_postEta.qs"))
      qs::qsave(postEta, file = postEta_file, preset = "fast")
      rm(postEta)
      invisible(gc())

      if (r == 1) {
        IASDT.R::CatTime("LF prediction using `Predict_LF`", Level = 1)
      }

      predPostEta[[r]] <- IASDT.R::Predict_LF(
        unitsPred = levels(dfPiNew[, r]),
        modelunits = levels(object$dfPi[, r]),
        postEta = postEta_file, postAlpha = postAlpha, rL = rL[[r]],
        NCores = NCores, Temp_Dir = Temp_Dir, Model_Name = Model_Name,
        UseTF = UseTF, TF_Environ = TF_Environ, TF_use_single = TF_use_single,
        LF_OutFile = LF_OutFile, LF_Return = LF_Return)

      try(fs::file_delete(postEta_file), silent = TRUE)
      rm(postEta_file)

      rowNames <- rownames(predPostEta[[r]][[1]])
      PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
    }

  } else {

    if (is.null(RC) || RC == "c") {
      IASDT.R::CatTime(
        paste0("Loading LF prediction from disk: `", LF_InputFile, "`"),
        Level = 1)
      for (r in seq_len(object$nr)) {
        predPostEta[[r]] <- IASDT.R::LoadAs(LF_InputFile[[r]], nthreads = 5)
        rowNames <- rownames(predPostEta[[r]][[1]])
        PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
      }
    }
  }

  if (object$nr > 0) {
    ppEta <- simplify2array(predPostEta)
  } else {
    ppEta <- matrix(list(), predN, 0)
  }

  # free some memory
  try(rm(predPostEta), silent = TRUE)

  invisible(gc())

  if (LF_Only) {
    return(LF_OutFile)
  }

  # # ..................................................................... ###

  IASDT.R::CatTime("Predicting")

  if (!exists("post")) {
    IASDT.R::CatTime("Loading post from disk", Level = 1)
    post <- qs::qread(post_file, nthreads = 5)
  }

  if (!is.null(RC)) {
    IASDT.R::CatTime(
      "Predicting data for response curve (sequentially)", Level = 1)
    preds <- lapply(
      seq_len(predN),
      function(pN, ...) {
        get1prediction(
          object, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
          ppEta[pN, ], PiNew, dfPiNew, nyNew, expected, mcmcStep)
      })
    IASDT.R::CatDiff(
      InitTime = .StartTime, Prefix = "Prediction was finished in ", Level = 1)
    return(preds)

  } else {

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
          Temp_Dir, paste0(Model_Name, "_preds_ch", .x, ".qs"))
        qs::qsave(Ch, file = ChunkFile, preset = "fast")
        return(ChunkFile)
      })

    rm(ChunkIDs, post, ppEta)
    invisible(gc())

    IASDT.R::CatTime(
      paste0("Preparing working on parallel using ", NCores, " cores"),
      Level = 1)

    seeds <- sample.int(.Machine$integer.max, predN)

    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    snow::clusterExport(
      cl = c1,
      list = c(
        "object", "X", "XRRR", "Yc", "Loff", "rL", "rLPar", "PiNew",
        "dfPiNew", "nyNew", "expected", "mcmcStep", "seeds", "chunk_size",
        "Chunks", "Temp_Dir", "Model_Name",
        "get1prediction"),
      envir = environment())

    invisible(snow::clusterEvalQ(
      cl = c1,
      expr = {
        sapply(
          c(
            "dplyr", "Rcpp", "RcppArmadillo", "Matrix", "float",
            "qs", "Hmsc", "purrr"),
          library, character.only = TRUE)
      }))

    IASDT.R::CatTime("Making predictions on parallel", Level = 1)

    pred <- snow::parLapply(
      cl = c1,
      x = seq_len(length(Chunks)),
      fun = function(Chunk) {

        ChunkFile <- Chunks[Chunk]
        Ch <- qs::qread(ChunkFile, nthreads = 5)
        ppEta <- Ch$ppEta
        post <- Ch$post
        rm(Ch)
        Seed <- (Chunk - 1) * chunk_size
        Seed <- seq(Seed + 1, Seed + chunk_size)
        Seed <- seeds[Seed]

        PredChunk <- purrr::map(
          .x = seq_len(chunk_size),
          .f = function(pN) {
            get1prediction(
              object = object, X = X, XRRR = XRRR, Yc = Yc,
              Loff = Loff, rL = rL, rLPar = rLPar, sam = post[[pN]],
              predPostEta = ppEta[pN], PiNew = PiNew,
              dfPiNew = dfPiNew, nyNew = nyNew, expected = expected,
              mcmcStep = mcmcStep, seed = Seed[pN])
          })

        # Species richness
        ChunkSR <-  simplify2array(lapply(X = PredChunk, FUN = rowSums)) %>%
          float::fl()
        dimnames(ChunkSR) <- NULL
        ChunkSR_File <- file.path(
          Temp_Dir, paste0("Pred_", Model_Name, "_ch", Chunk, "_SR.qs"))
        qs::qsave(ChunkSR, file = ChunkSR_File, preset = "fast")
        rm(ChunkSR)

        # Species predictions
        ChunkSp <- purrr::map_dfr(
          .x = seq_len(length(object$spNames)),
          .f = function(Sp) {
            SpD <- purrr::map(PredChunk, ~ .x[, Sp], ncol = 1) %>%
              simplify2array() %>%
              float::fl()
            dimnames(SpD) <- NULL

            ChunkSp_File <- file.path(
              Temp_Dir,
              paste0("Pred_", Model_Name, "_ch", Chunk, "_taxon", Sp, ".qs"))

            qs::qsave(SpD, file = ChunkSp_File, preset = "fast")

            cbind.data.frame(
              Chunk = Chunk, Sp = Sp, IAS_ID = object$spNames[Sp],
              ChunkSp_File = ChunkSp_File) %>%
              return()
          }) %>%
          dplyr::bind_rows(
            tibble::tibble(
              Chunk = Chunk, Sp = 0, IAS_ID = "SR",
              ChunkSp_File = ChunkSR_File),
            .)

        try(fs::file_delete(ChunkFile))
        rm(PredChunk)

        return(ChunkSp)
      })

    snow::stopCluster(c1)
    invisible(gc())

    pred <- tibble::tibble(dplyr::bind_rows(pred))

  }

  # # ..................................................................... ###

  IASDT.R::CatTime("Summarizing prediction outputs / Evaluation", Level = 1)

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  future::plan("future::cluster", workers = c1, gc = TRUE)
  on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

  Eval_DT <- dplyr::select(pred, -Chunk) %>%
    dplyr::group_nest(Sp, IAS_ID) %>%
    dplyr::mutate(data = purrr::map(data, unlist))

  Eval_DT <- future.apply::future_lapply(
    X = seq_len(nrow(Eval_DT)),
    FUN = function(ID) {

      Sp <- Eval_DT$Sp[[ID]]
      if (Sp == 0) {
        Sp2 <- "SR"
      } else {
        Sp2 <- paste0("taxon", Sp)
      }
      IAS_ID <- Eval_DT$IAS_ID[[ID]]
      data <- as.vector(Eval_DT$data[[ID]])

      SpDT <- purrr::map(data, qs::qread) %>%
        do.call(cbind, .) %>%
        as.double()

      # Mean prediction
      SpDT_Mean <- Rfast::rowmeans(SpDT)

      # standard deviaion of prediction
      SpDT_SD <- Rfast::rowVars(SpDT, std = TRUE)

      # Coefficient of variation
      SpDT_Mean0 <- SpDT_Mean
      # Replace very small mean values with reasonably small number to avoid
      # overflow warning
      SpDT_Mean0[SpDT_Mean0 < 1e-8] <- 1e-8
      SpDT_Cov <- SpDT_SD / SpDT_Mean0

      rm(SpDT)

      if (is.null(Pred_XY)) {
        Pred_XY <- object$rL$sample$s
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
        Pred_Dir, paste0("Pred_", Model_Name, "_", Sp2, ".qs"))

      qs::qsave(PredSummary, file = PredSummaryFile, preset = "fast")
      try(fs::file_delete(data), silent = TRUE)

      if (Evaluate && Sp2 != "SR") {
        if (is.null(Pred_PA)) {
          PresAbs <- object$Y[, Sp]
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

      invisible(gc())

      tibble::tibble(
        Sp = Sp, IAS_ID = IAS_ID, Path_pred = PredSummaryFile,
        RMSE = RMSE, AUC = AUC, Boyce = Boyce, TjurR2 = TjurR2) %>%
        return()
    },
    future.seed = TRUE, future.chunk.size = 1,
    future.globals = c(
      "Eval_DT", "Evaluate", "object", "Pred_Dir", "Model_Name",
      "Pred_PA", "Pred_XY"),
    future.packages = c(
      "dplyr", "Matrix", "purrr", "tibble", "Hmsc", "float", "qs",
      "Rfast", "caret", "pROC", "ecospat", "sf"))

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)

  # # ..................................................................... ###

  # Save predictions for all species in a single file
  IASDT.R::CatTime("Save predictions for all species in a single file")

  Eval_DT <- dplyr::bind_rows(Eval_DT)

  Predictions <- dplyr::select(Eval_DT, Path_pred, Sp, IAS_ID) %>%
    dplyr::mutate(
      Sp_data = purrr::map(
        .x = Path_pred,
        .f = ~ {
          qs::qread(.x) %>%
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

  Pred_File <- file.path(Pred_Dir, paste0("Prediction_", Model_Name, ".qs"))
  qs::qsave(Predictions, file = Pred_File, preset = "fast")
  try(fs::file_delete(Eval_DT$Path_pred), silent = TRUE)
  IASDT.R::CatTime(
    paste0("Predictions were saved to `", Pred_File, "`"), Level = 1)

  if (Evaluate) {
    if (is.null(Eval_Name)) {
      Eval_Path <- file.path(Eval_Dir, paste0("Eval_", Model_Name, ".qs"))
    } else {
      Eval_Path <- file.path(
        Eval_Dir, paste0("Eval_", Model_Name, "_", Eval_Name, ".qs"))
    }

    Eval_DT <- dplyr::select(Eval_DT, -Path_pred)
    qs::qsave(Eval_DT, file = Eval_Path, preset = "fast")
    IASDT.R::CatTime(
      paste0(
        "Evaluation results were saved to `",
        file.path(Eval_Dir, "Eval_DT.qs"), "`"),
      Level = 1)
  } else {
    Eval_Path <- NULL
  }

  if (exists("post_file")) {
    try(fs::file_delete(post_file), silent = TRUE)
  }

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Prediction was finished in ")

  # # ..................................................................... ###

  if (is.null(LF_OutFile)) {
    LF_Path <- LF_InputFile
  } else {
    LF_Path <- LF_OutFile
  }

  tibble::tibble(
    Pred_Path = Pred_File, Eval_Path = Eval_Path, LF_Path = LF_Path) %>%
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

  LRan <- vector("list", object$nr)
  Eta <- vector("list", object$nr)

  for (r in seq_len(object$nr)) {
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
    for (r in seq_len(object$nr)) {
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
