## |------------------------------------------------------------------------| #
# predict_hmsc ----
## |------------------------------------------------------------------------| #

#' Calculates predicted values from a fitted Hmsc model
#'
#' This function modifies the `Hmsc:::predict.Hmsc` function.
#'
#' @param path_model Character. Path to the saved model object.
#' @param x_data `data.frame`. The unpreprocessed covariates for the predictions
#'   to be made. Works only if the `XFormula` argument was specified in the
#'   [Hmsc::Hmsc] model constructor call. Requirements are similar to those in
#'   the `Hmsc` model constructor.
#' @param X `matrix`. Covariates for the predictions to be made. Only one of
#'   `x_data` and `X` arguments may be provided.
#' @param XRRRData `data.frame`. Covariates for reduced-rank regression.
#' @param XRRR `matrix`. Covariates for reduced-rank regression.
#' @param gradient An object returned by [Hmsc::constructGradient]. Providing
#'   `gradient` is an alternative for providing `x_data`, `studyDesign` and
#'   `ranLevels`. Cannot be used together with `Yc`.
#' @param Yc `matrix`. Outcomes that are assumed to be known for conditional
#'   predictions. Cannot be used together with `gradient`.
#' @param mcmcStep Integer. Number of extra `mcmc` steps used for updating the
#'   random effects. Defaults to 1L.
#' @param expected Logical. Whether to return the location parameter of the
#'   observation models or sample the values from those. Defaults to `TRUE`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param temp_cleanup Logical. Whether to clean up temporary files. Defaults to
#'   `TRUE`.
#' @param pred_directory Character. Directory path indicating where the
#'   predictions will be saved. Defaults to `NULL`, which saves model
#'   predictions to "`Model_Prediction`" folder of the current working
#'   directory.
#' @param evaluate Logical. Whether to evaluate the model predictions. Defaults
#'   to `FALSE`.
#' @param evaluation_name Character. Name of the evaluation results. If `NULL`,
#'   the default name is used (`Eval_[model_name].qs2`).
#' @param evaluation_directory Character. Directory where the evaluation results
#'   will be saved. Defaults to `Evaluation`.
#' @param prediction_type Character. Type of predictions to be made. If `NULL`
#'   (default), predictions are made for the latent factors. If `c`, predictions
#'   are made for response curves at mean coordinates. If `i`, predictions are
#'   made for response curves at infinite coordinates.
#' @param pred_PA `matrix`. Presence-absence data for evaluation. If `NULL`
#'   (default), the presence-absence data from the model object is used. This
#'   argument is used only when `evaluate` is `TRUE`.
#' @param pred_XY `matrix`. Coordinates to be added to predicted values. If
#'   `NULL` (default), the coordinates from the model object is used.
#' @param LF_inputFile Character. File name where the latent factor predictions
#'   are saved. If `NULL` (default), latent factor predictions will be computed.
#'   If specified, latent factor predictions are read from this path. This
#'   allows to predicting the latent factors for new sites only once.
#' @param LF_only Logical. Whether to return the latent factor predictions only.
#'   Defaults to `FALSE`. This helps in predicting to new sites, allowing to
#'   predicting the latent factors only once, then the output can be loaded in
#'   other predictions when needed.
#' @param verbose Logical. Whether to print a message upon successful saving of
#'   files. Defaults to `FALSE`.
#' @param Loff See [Hmsc::predict.Hmsc] for more details.
#' @inheritParams predict_latent_factor
#' @name predict_hmsc
#' @export

predict_hmsc <- function(
    path_model, Loff = NULL, x_data = NULL, X = NULL, XRRRData = NULL,
    XRRR = NULL, gradient = NULL, Yc = NULL, mcmcStep = 1L, expected = TRUE,
    n_cores = 8L, model_name = "Train", temp_dir = "TEMP_Pred",
    temp_cleanup = TRUE, prediction_type = NULL, use_TF = TRUE,
    TF_environ = NULL, TF_use_single = FALSE, LF_out_file = NULL,
    LF_return = FALSE, LF_inputFile = NULL, LF_only = FALSE,
    LF_n_cores = n_cores, LF_check = FALSE, LF_temp_cleanup = TRUE,
    LF_commands_only = FALSE, pred_directory = NULL, pred_PA = NULL,
    pred_XY = NULL, evaluate = FALSE, evaluation_name = NULL,
    evaluation_directory = "Evaluation", verbose = TRUE) {

  # # ..................................................................... ###

  if (isFALSE(verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .start_time <- lubridate::now(tzone = "CET")

  # To avoid non-standard evaluation
  pred_XY <- pred_XY
  pred_PA <- pred_PA

  if (LF_only && is.null(LF_out_file)) {
    IASDT.R::stop_ctx(
      "`LF_out_file` must be specified when `LF_only` is `TRUE`",
      LF_out_file = LF_out_file, LF_only = LF_only)
  }

  if (!is.null(LF_out_file) && !is.null(LF_inputFile)) {
    IASDT.R::stop_ctx(
      "only one of `LF_out_file` and `LF_inputFile` arguments can be specified",
      LF_out_file = LF_out_file, LF_inputFile = LF_inputFile)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chunk <- Sp <- IAS_ID <- Path_pred <- Sp_data <- data <- geometry <- x <-
    y <- SR_mean <- SR_sd <- SR_cov <- NULL

  # # ..................................................................... ###

  if (!is.null(prediction_type) && !prediction_type %in% c("c", "i")) {
    IASDT.R::stop_ctx(
      "`prediction_type` must be either NULL or one of 'c' or 'i'",
      prediction_type = prediction_type)
  }

  if (is.null(pred_directory)) {
    pred_directory <- fs::path(dirname(dirname(path_model)), "Model_Prediction")
  }

  if (is.null(prediction_type) || prediction_type == "c") {
    IASDT.R::cat_time("Creating/checking output paths")
    fs::dir_create(c(temp_dir, pred_directory))
    if (evaluate) {
      fs::dir_create(evaluation_directory)
    }
  }

  # # ..................................................................... ###

  # Load model if it is a character
  if (inherits(path_model, "character")) {
    IASDT.R::cat_time("Load model object")
    Model <- IASDT.R::load_as(path_model)
  }

  # # ..................................................................... ###

  # Combines a list of single or several MCMC chains into a single chain
  IASDT.R::cat_time("Combine list of posteriors")
  post <- Hmsc::poolMcmcChains(Model$postList)
  studyDesign <- Model$studyDesign
  ranLevels <- Model$ranLevels

  # # ..................................................................... ###

  IASDT.R::cat_time("Check input parameters")

  if (!is.null(gradient)) {
    if (!is.null(Yc)) {
      IASDT.R::stop_ctx(
        paste0(
          "predict with arguments 'Yc' and 'gradient' jointly is not ",
          "implemented (yet)"),
        gradient = gradient, Yc = Yc)
    }
    x_data <- gradient$XDataNew
    studyDesign <- gradient$studyDesignNew
    ranLevels <- gradient$rLNew
  }

  if (!is.null(x_data) && !is.null(X)) {
    IASDT.R::stop_ctx(
      "only one of x_data and X arguments can be specified",
      x_data = x_data, X = X)
  }

  if (!is.null(XRRRData) && !is.null(XRRR)) {
    IASDT.R::stop_ctx(
      "only one of XRRRData and XRRR arguments can be specified",
      XRRRData = XRRRData, XRRR = XRRR)
  }

  if (!is.null(x_data)) {
    class_x <- class(x_data)[1L]

    if (class_x == "list") {
      if (any(unlist(lapply(x_data, is.na)))) {
        IASDT.R::stop_ctx(
          "NA values are not allowed in 'x_data'", class_x = class_x)
      }
      xlev <- lapply(Reduce(rbind, Model$XData), levels)
      xlev <- xlev[unlist(lapply(Reduce(rbind, Model$XData), is.factor))]
      X <- lapply(x_data, function(a) {
        stats::model.matrix(Model$XFormula, a, xlev = xlev)
      })

    } else if (class_x == "data.frame") {
      if (anyNA(x_data)) {
        IASDT.R::stop_ctx(
          "NA values are not allowed in 'x_data'", class_x = class_x)
      }
      xlev <- lapply(Model$XData, levels)
      xlev <- xlev[unlist(lapply(Model$XData, is.factor))]
      X <- stats::model.matrix(Model$XFormula, x_data, xlev = xlev)
    }

  } else if (is.null(X)) {
    X <- Model$X
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
      IASDT.R::stop_ctx(
        "number of columns in Yc must be equal to ns",
        ncol_Yc = ncol(Yc), model_ns = Model$ns)
    }
    if (nrow(Yc) != nyNew) {
      IASDT.R::stop_ctx(
        "number of rows in Yc and X must be equal",
        nrow_Yc = nrow(Yc), nyNew = nyNew)
    }
  }

  if (!is.null(Loff)) {
    if (ncol(Loff) != Model$ns) {
      IASDT.R::stop_ctx(
        "number of columns in Loff must be equal to ns",
        ncol_Loff = ncol(Loff), model_ns = Model$ns)
    }
    if (nrow(Loff) != nyNew) {
      IASDT.R::stop_ctx(
        "number of rows in Loff and X must be equal",
        nrow_Loff = nrow(Loff), nyNew = nyNew)
    }
  }

  if (!all(Model$rLNames %in% colnames(studyDesign))) {
    IASDT.R::stop_ctx(
      "dfPiNew does not contain all the necessary named columns",
      model_rLNames = Model$rLNames, names_study_design = colnames(studyDesign))
  }

  if (!all(Model$rLNames %in% names(ranLevels))) {
    IASDT.R::stop_ctx(
      "rL does not contain all the necessary named levels",
      model_rLNames = Model$rLNames, names_ranLevels = names(ranLevels))
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
  IASDT.R::cat_time("Free some memory")
  Model$postList <- Model$YScaled <- Model$X <- Model$XScaled <- NULL

  Mod_nr <- Model$nr
  Mod_dfPi <- Model$dfPi

  # Save smaller version of the model object for later use
  Model_File_small <- fs::path(temp_dir, paste0(model_name, "Model_small.qs2"))

  if (!file.exists(Model_File_small)) {
    IASDT.R::cat_time(
      "Save smaller version of the model object to disk",
      level = 1L)
    IASDT.R::save_as(object = Model, out_path = Model_File_small)
  }
  rm(Model, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  IASDT.R::cat_time("Predict Latent Factor")

  predPostEta <- vector("list", Mod_nr)
  PiNew <- matrix(NA, nrow(dfPiNew), Mod_nr)

  # Whether to use `predict_latent_factor` or read its results from file
  if (!is.null(Mod_nr)) {
    if (is.null(LF_inputFile) || length(LF_inputFile) != Mod_nr) {
      CalcLF <- TRUE
    } else {
      CalcLF <- FALSE
    }
  }

  # Do not use `predict_latent_factor` when predicting values for response
  # curves when using coordinates = "i" in constructGradient
  if (!is.null(prediction_type)) {
    if (prediction_type == "i") {
      for (r in seq_len(Mod_nr)) {
        if (r == 1) {
          IASDT.R::cat_time(
            "LF prediction for response curve with infinite coordinates",
            level = 1L)
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

      post_alpha <- lapply(post, function(c) c$Alpha[[r]])

      # Save postEta to file and load it from predict_latent_factor. This helps
      # to avoid the unnecessary copying of the postEta object to all cores
      postEta_file <- fs::path(
        temp_dir, paste0(model_name, "_r", r, "_postEta.qs2"))

      if (isFALSE(IASDT.R::check_data(postEta_file, warning = FALSE))) {
        IASDT.R::cat_time("Save postEta to file", level = 1L)
        postEta <- lapply(post, function(c) c$Eta[[r]])
        IASDT.R::save_as(object = postEta, out_path = postEta_file)
        rm(postEta, envir = environment())
        invisible(gc())
      }

      # Save post to file and load it later
      if (r == Mod_nr) {
        post_file <- fs::path(temp_dir, paste0(model_name, "_post.qs2"))

        if (isFALSE(IASDT.R::check_data(post_file, warning = FALSE))) {
          # free some memory
          post <- lapply(post, function(x) {
            x$Eta <- x$Psi <- x$V <- x$Delta <- x$Gamma <- x$rho <- NULL
            x
          })
          IASDT.R::cat_time("Save post to file", level = 1L)
          IASDT.R::save_as(object = post, out_path = post_file)
        }
        rm(post, envir = environment())
        invisible(gc())
      }

      if (r == 1) {
        IASDT.R::cat_time(
          "LF prediction using `predict_latent_factor`", level = 1L)
      }

      predPostEta[[r]] <- IASDT.R::predict_latent_factor(
        units_pred = levels(dfPiNew[, r]),
        units_model = levels(Mod_dfPi[, r]),
        postEta = postEta_file, post_alpha = post_alpha, LF_rL = rL[[r]],
        LF_n_cores = LF_n_cores, LF_temp_cleanup = LF_temp_cleanup,
        LF_out_file = LF_out_file, LF_return = LF_return, LF_check = LF_check,
        LF_commands_only = LF_commands_only, temp_dir = temp_dir,
        model_name = model_name, use_TF = use_TF, TF_environ = TF_environ,
        TF_use_single = TF_use_single)

      rm(postEta_file, envir = environment())

      if (LF_commands_only) {
        return(invisible(NULL))
      }

      rowNames <- rownames(predPostEta[[r]][[1]])
      PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
    }

  } else {

    if (is.null(prediction_type) || prediction_type == "c") {

      IASDT.R::cat_time("Loading LF prediction from disk", level = 1L)
      IASDT.R::cat_time(LF_inputFile, level = 2L, cat_timestamp = FALSE)

      for (r in seq_len(Mod_nr)) {
        predPostEta[[r]] <- IASDT.R::load_as(LF_inputFile[[r]])
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
  if (LF_only) {
    return(LF_out_file)
  }

  # free some memory
  rm(predPostEta, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  IASDT.R::cat_time("Predicting")

  if (!exists("post")) {
    IASDT.R::cat_time("Loading post from disk", level = 1L)
    post <- IASDT.R::load_as(post_file)
  }

  # Read model object from disk
  Model <- IASDT.R::load_as(Model_File_small)

  # prediction data for response curves
  if (!is.null(prediction_type)) {
    IASDT.R::cat_time(
      "Predicting data for response curve (sequentially)",
      level = 1L)

    preds <- lapply(
      seq_len(predN),
      function(pN, ...) {
        get1prediction(
          Model, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
          ppEta[pN, ], PiNew, dfPiNew, nyNew, expected, mcmcStep)
      })

    IASDT.R::cat_diff(
      init_time = .start_time, prefix = "Prediction was finished in ",
      level = 1L)
    return(preds)
  }

  # Save ppEta / post as small chunks
  IASDT.R::cat_time("Save ppEta / post as small chunks", level = 1L)
  chunk_size <- 25
  ChunkIDs <- ceiling(seq_along(post) / chunk_size)
  Chunks <- purrr::map_chr(
    .x = seq_len(max(ChunkIDs)),
    .f = ~ {
      IDs <- which(ChunkIDs == .x)
      Ch <- list(ppEta = ppEta[IDs], post = post[IDs])
      chunk_file <- fs::path(
        temp_dir, paste0(model_name, "_preds_ch", .x, ".qs2"))
      IASDT.R::save_as(object = Ch, out_path = chunk_file)
      return(chunk_file)
    })

  rm(ChunkIDs, post, ppEta, envir = environment())
  invisible(gc())

  seeds <- sample.int(.Machine$integer.max, predN)


  IASDT.R::set_parallel(n_cores = min(n_cores, length(Chunks)), level = 1L)
  withr::defer(future::plan("future::sequential", gc = TRUE))

  IASDT.R::cat_time("Making predictions in parallel", level = 1L)
  pred <- future.apply::future_lapply(
    X = seq_len(length(Chunks)),
    FUN = function(Chunk) {

      chunk_file <- Chunks[Chunk]
      Ch <- IASDT.R::load_as(chunk_file)
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
      ChunkSR_File <- fs::path(
        temp_dir, paste0("Pred_", model_name, "_ch", Chunk, "_SR.qs2"))
      IASDT.R::save_as(object = ChunkSR, out_path = ChunkSR_File)

      rm(ChunkSR, envir = environment())

      # Species predictions
      ChunkSp <- purrr::map_dfr(
        .x = seq_len(length(Model$spNames)),
        .f = function(Sp) {

          SpD <- purrr::map(PredChunk, ~ .x[, Sp], ncol = 1) %>%
            simplify2array() %>%
            float::fl()
          dimnames(SpD) <- NULL

          ChunkSp_File <- fs::path(
            temp_dir,
            paste0("Pred_", model_name, "_ch", Chunk, "_taxon", Sp, ".qs2"))

          IASDT.R::save_as(object = SpD, out_path = ChunkSp_File)

          return(
            cbind.data.frame(
              Chunk = Chunk, Sp = Sp, IAS_ID = Model$spNames[Sp],
              ChunkSp_File = ChunkSp_File))
        }) %>%
        dplyr::bind_rows(
          tibble::tibble(
            Chunk = Chunk, Sp = 0, IAS_ID = "SR", ChunkSp_File = ChunkSR_File),
          .)

      rm(PredChunk, envir = environment())

      invisible(gc())
      return(ChunkSp)
    },
    future.seed = TRUE,
    future.globals = c(
      "Model", "X", "XRRR", "Yc", "Loff", "rL", "rLPar", "PiNew",
      "dfPiNew", "nyNew", "expected", "mcmcStep", "seeds", "chunk_size",
      "Chunks", "temp_dir", "model_name", "get1prediction"),
    future.packages = c(
      "dplyr", "Rcpp", "RcppArmadillo", "Matrix", "float", "qs2", "Hmsc",
      "purrr", "tibble", "Hmsc", "Rfast", "caret", "pROC", "ecospat", "sf"))


  pred <- tibble::tibble(dplyr::bind_rows(pred))


  # # ..................................................................... ###

  IASDT.R::cat_time("Summarizing prediction outputs / Evaluation", level = 1L)

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

      SpDT <- purrr::map(data, IASDT.R::load_as) %>%
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

      if (is.null(pred_XY)) {
        pred_XY <- Model$rL$sample$s
      }

      PredSummary <- tibble::tibble(
        as.data.frame(pred_XY),
        Mean = SpDT_Mean, SD = SpDT_SD, Cov = SpDT_Cov) %>%
        stats::setNames(
          c(
            "x", "y", paste0(IAS_ID, "_mean"),
            paste0(IAS_ID, "_sd"), paste0(IAS_ID, "_cov"))) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

      PredSummaryFile <- fs::path(
        pred_directory, paste0("Pred_", model_name, "_", Sp2, ".qs2"))

      IASDT.R::save_as(object = PredSummary, out_path = PredSummaryFile)

      if (evaluate && Sp2 != "SR") {
        if (is.null(pred_PA)) {
          PresAbs <- Model$Y[, Sp]
        } else {
          PresAbs <- pred_PA[, Sp]
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

      return(
        tibble::tibble(
          Sp = Sp, IAS_ID = IAS_ID, Path_pred = PredSummaryFile,
          RMSE = RMSE, AUC = AUC, Boyce = Boyce, TjurR2 = TjurR2))

    },
    future.seed = TRUE,
    future.globals = c(
      "Eval_DT", "evaluate", "pred_directory", "model_name",
      "pred_PA", "pred_XY"),
    future.packages = c(
      "dplyr", "Rcpp", "RcppArmadillo", "Matrix", "float", "qs2", "Hmsc",
      "purrr", "tibble", "Hmsc", "Rfast", "caret", "pROC", "ecospat", "sf"))

  IASDT.R::set_parallel(stop_cluster = TRUE, level = 1L)

  invisible(gc())

  # # ..................................................................... ###

  # Save predictions for all species in a single file
  IASDT.R::cat_time("Save predictions for all species in a single file")

  Eval_DT <- dplyr::bind_rows(Eval_DT)

  Predictions <- dplyr::select(Eval_DT, Path_pred, Sp, IAS_ID) %>%
    dplyr::mutate(
      Sp_data = purrr::map(
        .x = Path_pred,
        .f = ~ {
          IASDT.R::load_as(.x) %>%
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
    dplyr::mutate(model_name = model_name, .before = "SR_mean")

  Pred_File <- fs::path(
    pred_directory,
    paste0(
      "Prediction_",
      stringr::str_remove(model_name, "_Clamping|_NoClamping"), ".qs2"))

  IASDT.R::save_as(object = Predictions, out_path = Pred_File)

  if (temp_cleanup) {
    try(fs::file_delete(Eval_DT$Path_pred), silent = TRUE)
  }

  IASDT.R::cat_time("Predictions were saved", level = 1L)
  IASDT.R::cat_time(Pred_File, level = 2L, cat_timestamp = FALSE)

  if (evaluate) {
    if (is.null(evaluation_name)) {
      Eval_Path <- fs::path(
        evaluation_directory,
        paste0(
          "Eval_", stringr::str_remove(model_name, "_Train|_Current"), ".qs2"))
    } else {
      Eval_Path <- fs::path(
        evaluation_directory,
        paste0(
          "Eval_", stringr::str_remove(model_name, "_Train|_Current"),
          "_", evaluation_name, ".qs2"))
    }

    Eval_DT <- dplyr::select(Eval_DT, -Path_pred)
    IASDT.R::save_as(object = Eval_DT, out_path = Eval_Path)

    IASDT.R::cat_time("Evaluation results were saved", level = 1L)
    IASDT.R::cat_time(
      fs::path(evaluation_directory, "Eval_DT.qs2"),
      level = 2L, cat_timestamp = FALSE)

  } else {
    Eval_Path <- NULL
  }

  # # ..................................................................... ###

  if (is.null(LF_out_file)) {
    LF_Path <- LF_inputFile
  } else {
    LF_Path <- LF_out_file
  }

  # # ..................................................................... ###

  # Clean up
  if (temp_cleanup) {

    IASDT.R::cat_time("Cleaning up temporary files", level = 1L)

    try(
      {
        Pattern <- paste0("(Pred_){0,}", model_name, ".+qs2")
        file_paths <- list.files(
          path = IASDT.R::normalize_path(temp_dir),
          pattern = Pattern, full.names = TRUE)
        try(fs::file_delete(file_paths), silent = TRUE)
      },
      silent = TRUE)

    try(fs::file_delete(Model_File_small), silent = TRUE)

  }

  # # ..................................................................... ###

  IASDT.R::cat_diff(
    init_time = .start_time, prefix = "Prediction was finished in ")

  return(
    tibble::tibble(
      Pred_Path = Pred_File,  Eval_Path = Eval_Path, LF_Path = LF_Path))
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
  if (!is.null(Yc) && !all(is.na(Yc))) {
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

  if (expected) {
    Z <- L
  } else {
    Z <- L + matrix(sqrt(sam$sigma), nrow(L), object$ns, byrow = TRUE) *
      matrix(stats::rnorm(nrow(L) * object$ns), nrow(L), object$ns)
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
