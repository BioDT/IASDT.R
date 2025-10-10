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
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param future_max_size	Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details.
#' @param temp_cleanup Logical. Whether to clean up temporary files. Defaults to
#'   `TRUE`.
#' @param pred_directory Character. Directory path indicating where the
#'   predictions will be saved. Defaults to `NULL`, which saves model
#'   predictions to "`model_prediction`" folder of the current working
#'   directory.
#' @param evaluate Logical. Whether to evaluate the model predictions. Defaults
#'   to `FALSE`.
#' @param evaluation_name Character. Name of the evaluation results. If `NULL`,
#'   the default name is used (`eval_[model_name].qs2`).
#' @param evaluation_directory Character. Directory where the evaluation results
#'   will be saved. Defaults to `evaluation`.
#' @param prediction_type Character. Type of predictions to be made. If `NULL`
#'   (default), predictions are made for the latent factors. If `c`, predictions
#'   are made for response curves at mean coordinates. If `i`, predictions are
#'   made for response curves at infinite coordinates.
#' @param pred_pa `matrix`. Presence-absence data for evaluation. If `NULL`
#'   (default), the presence-absence data from the model object is used. This
#'   argument is used only when `evaluate` is `TRUE`.
#' @param pred_xy `matrix`. Coordinates to be added to predicted values. If
#'   `NULL` (default), the coordinates from the model object is used.
#' @param lf_input_file Character. File name where the latent factor predictions
#'   are saved. If `NULL` (default), latent factor predictions will be computed.
#'   If specified, latent factor predictions are read from this path. This
#'   allows to predicting the latent factors for new sites only once.
#' @param lf_only Logical. Whether to return the latent factor predictions only.
#'   Defaults to `FALSE`. This helps in predicting to new sites, allowing to
#'   predicting the latent factors only once, then the output can be loaded in
#'   other predictions when needed.
#' @param verbose Logical. Whether to print a message upon successful saving of
#'   files. Defaults to `FALSE`.
#' @param spatial_model Logical. Whether the fitted model is a spatial model.
#'   Defaults to `TRUE`.
#' @param Loff See [Hmsc::predict.Hmsc] for more details.
#' @inheritParams predict_latent_factor
#' @name predict_hmsc
#' @export

predict_hmsc <- function(
    path_model, Loff = NULL, x_data = NULL, X = NULL, XRRRData = NULL,
    XRRR = NULL, gradient = NULL, Yc = NULL, mcmcStep = 1L, expected = TRUE,
    n_cores = 8L, strategy = "multisession", future_max_size = 1000L,
    model_name = "train", temp_dir = "temp_pred", temp_cleanup = TRUE,
    prediction_type = NULL, use_tf = TRUE, tf_environ = NULL,
    tf_use_single = FALSE, lf_out_file = NULL, lf_return = FALSE,
    lf_input_file = NULL, lf_only = FALSE, n_cores_lf = n_cores,
    lf_check = FALSE, lf_temp_cleanup = TRUE, lf_commands_only = FALSE,
    pred_directory = NULL, pred_pa = NULL, pred_xy = NULL, evaluate = FALSE,
    evaluation_name = NULL, evaluation_directory = "evaluation", verbose = TRUE,
    spatial_model = TRUE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- n_cores_lf <- 1L
  n_cores <- .validate_n_cores(n_cores)
  n_cores_lf <- .validate_n_cores(n_cores_lf)

  # To avoid non-standard evaluation
  pred_xy <- pred_xy
  pred_pa <- pred_pa

  if (lf_only && is.null(lf_out_file)) {
    ecokit::stop_ctx(
      "`lf_out_file` must be specified when `lf_only` is `TRUE`",
      lf_out_file = lf_out_file, lf_only = lf_only, include_backtrace = TRUE)
  }

  if (!is.null(lf_out_file) && !is.null(lf_input_file)) {
    ecokit::stop_ctx(
      paste0(
        "only one of `lf_out_file` and `lf_input_file` arguments",
        "can be specified"),
      lf_out_file = lf_out_file, lf_input_file = lf_input_file,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  chunk <- sp <- ias_id <- path_pred <- sp_data <- data <- geometry <- x <-
    y <- sr_mean <- sr_sd <- sr_cov <- NULL

  # # ..................................................................... ###

  if (!is.null(prediction_type) && !prediction_type %in% c("c", "i")) {
    ecokit::stop_ctx(
      "`prediction_type` must be either NULL or one of 'c' or 'i'",
      prediction_type = prediction_type, include_backtrace = TRUE)
  }

  if (is.null(pred_directory)) {
    pred_directory <- fs::path(dirname(dirname(path_model)), "model_prediction")
  }

  if (is.null(prediction_type) || prediction_type == "c") {
    ecokit::cat_time("Creating/checking output paths", verbose = verbose)
    fs::dir_create(c(temp_dir, pred_directory))
    if (evaluate) {
      fs::dir_create(evaluation_directory)
    }
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "Rcpp", "RcppArmadillo", "Matrix", "float", "qs2", "Hmsc",
      "purrr", "tibble", "Rfast", "caret", "pROC", "ecospat", "sf",
      "ecokit", "magrittr"),
    strategy = strategy)

  # # ..................................................................... ###

  # Load model if it is a character
  if (inherits(path_model, "character")) {
    ecokit::cat_time("Load model object", verbose = verbose)
    model_obj <- ecokit::load_as(path_model)
  }

  # # ..................................................................... ###

  # Combines a list of single or several MCMC chains into a single chain
  ecokit::cat_time("Combine list of posteriors", verbose = verbose)
  post <- Hmsc::poolMcmcChains(model_obj$postList)
  studyDesign <- model_obj$studyDesign
  ranLevels <- model_obj$ranLevels

  # # ..................................................................... ###

  ecokit::cat_time("Check input parameters", verbose = verbose)

  if (!is.null(gradient)) {
    if (!is.null(Yc)) {
      ecokit::stop_ctx(
        paste0(
          "predict with arguments 'Yc' and 'gradient' jointly is not ",
          "implemented (yet)"),
        gradient = gradient, Yc = Yc, include_backtrace = TRUE)
    }
    x_data <- gradient$XDataNew
    studyDesign <- gradient$studyDesignNew
    ranLevels <- gradient$rLNew
  }

  if (!is.null(x_data) && !is.null(X)) {
    ecokit::stop_ctx(
      "only one of x_data and X arguments can be specified",
      x_data = x_data, X = X, include_backtrace = TRUE)
  }

  if (!is.null(XRRRData) && !is.null(XRRR)) {
    ecokit::stop_ctx(
      "only one of XRRRData and XRRR arguments can be specified",
      XRRRData = XRRRData, XRRR = XRRR, include_backtrace = TRUE)
  }

  if (!is.null(x_data)) {
    class_x <- class(x_data)[1L]

    if (class_x == "list") {
      if (any(unlist(lapply(x_data, is.na)))) {
        ecokit::stop_ctx(
          "NA values are not allowed in 'x_data'", class_x = class_x,
          include_backtrace = TRUE)
      }
      xlev <- lapply(Reduce(rbind, model_obj$XData), levels)
      xlev <- xlev[unlist(lapply(Reduce(rbind, model_obj$XData), is.factor))]
      X <- lapply(x_data, function(a) {
        stats::model.matrix(model_obj$XFormula, a, xlev = xlev)
      })

    } else if (class_x == "data.frame") {
      if (anyNA(x_data)) {
        ecokit::stop_ctx(
          "NA values are not allowed in 'x_data'", class_x = class_x,
          include_backtrace = TRUE)
      }
      xlev <- lapply(model_obj$XData, levels)
      xlev <- xlev[unlist(lapply(model_obj$XData, is.factor))]
      X <- stats::model.matrix(model_obj$XFormula, x_data, xlev = xlev)
    }

  } else if (is.null(X)) {
    X <- model_obj$X
  }

  if (!is.null(XRRRData)) {
    xlev <- lapply(model_obj$XRRRData, levels)
    xlev <- xlev[unlist(lapply(model_obj$XRRRData, is.factor))]
    XRRR <- stats::model.matrix(model_obj$XRRRFormula, XRRRData, xlev = xlev)
  } else {
    if (is.null(model_obj$ncRRR)) {
      model_obj$ncRRR <- 0
    }
    if (is.null(XRRR) && model_obj$ncRRR > 0) {
      XRRR <- model_obj$XRRR
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
    if (ncol(Yc) != model_obj$ns) {
      ecokit::stop_ctx(
        "number of columns in Yc must be equal to ns",
        ncol_Yc = ncol(Yc), model_ns = model_obj$ns, include_backtrace = TRUE)
    }
    if (nrow(Yc) != nyNew) {
      ecokit::stop_ctx(
        "number of rows in Yc and X must be equal",
        nrow_Yc = nrow(Yc), nyNew = nyNew, include_backtrace = TRUE)
    }
  }

  if (!is.null(Loff)) {
    if (ncol(Loff) != model_obj$ns) {
      ecokit::stop_ctx(
        "number of columns in Loff must be equal to ns",
        ncol_Loff = ncol(Loff), model_ns = model_obj$ns,
        include_backtrace = TRUE)
    }
    if (nrow(Loff) != nyNew) {
      ecokit::stop_ctx(
        "number of rows in Loff and X must be equal",
        nrow_Loff = nrow(Loff), nyNew = nyNew, include_backtrace = TRUE)
    }
  }

  if (!all(model_obj$rLNames %in% colnames(studyDesign))) {
    ecokit::stop_ctx(
      "dfPiNew does not contain all the necessary named columns",
      model_rLNames = model_obj$rLNames,
      names_study_design = colnames(studyDesign), include_backtrace = TRUE)
  }

  if (!all(model_obj$rLNames %in% names(ranLevels))) {
    ecokit::stop_ctx(
      "rL does not contain all the necessary named levels",
      model_rLNames = model_obj$rLNames, names_ranLevels = names(ranLevels),
      include_backtrace = TRUE)
  }

  if (!is.null(studyDesign)) {
    dfPiNew <- studyDesign[, model_obj$rLNames, drop = FALSE]
  } else {
    dfPiNew <- matrix(NA, nyNew, 0)
  }
  rL <- ranLevels[model_obj$rLNames]


  if (!is.null(Yc)) {
    ## object can have pre-computed data parameters, but not
    ## necessarily. These are needed only in updateEta(), but get it
    ## here anyway...
    if (is.null(model_obj$rLPar)) {
      rLPar <- Hmsc::computeDataParameters(model_obj)$rLPar
    } else {
      rLPar <- model_obj$rLPar
    }
  } else {
    rLPar <- NULL
  }

  predN <- length(post)

  # # ..................................................................... ###

  # free some memory
  ecokit::cat_time("Free some memory", verbose = verbose)
  model_obj <- IASDT.R::trim_hmsc(
    model_obj, c("postList", "YScaled", "X", "XScaled"))

  mod_nr <- model_obj$nr
  mod_dfPi <- model_obj$dfPi

  # Save smaller version of the model object for later use
  model_file_small <- fs::path(temp_dir, paste0(model_name, "model_small.qs2"))

  if (!file.exists(model_file_small)) {
    ecokit::cat_time(
      "Save smaller version of the model object to disk",
      level = 1L, verbose = verbose)
    ecokit::save_as(object = model_obj, out_path = model_file_small)
  }
  rm(model_obj, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Predict Latent Factor", verbose = verbose)

  predPostEta <- vector("list", mod_nr)
  PiNew <- matrix(NA, nrow(dfPiNew), mod_nr)

  # Whether to use `predict_latent_factor` or read its results from file
  if (!is.null(mod_nr)) {
    if (is.null(lf_input_file) || length(lf_input_file) != mod_nr) {
      CalcLF <- TRUE
    } else {
      CalcLF <- FALSE
    }
  }

  # Do not use `predict_latent_factor` when predicting values for response
  # curves when using coordinates = "i" in constructGradient or for making
  # predictions for non-spatial models

  rc_predict_i <- !is.null(prediction_type) && prediction_type == "i"

  if (rc_predict_i || isFALSE(spatial_model)) {

    nLF <- length(post[[1]]$Alpha[[1]])

    for (r in seq_len(mod_nr)) {
      if (rc_predict_i) {
        if (r == 1) {
          ecokit::cat_time(
            "LF prediction for response curve with infinite coordinates",
            level = 1L, verbose = verbose)
        }
        predPostEta[[r]] <- replicate(
          n = predN,
          expr = matrix(
            0, nrow = 1L, ncol = nLF, dimnames = list("new_unit", NULL)),
          simplify = FALSE)

      } else {
        if (r == 1) {
          ecokit::cat_time(
            "No LF prediction is needed for non-spatial models",
            level = 1L, verbose = verbose)
        }
        predPostEta[[r]] <- replicate(
          n = predN,
          expr = matrix(
            0, nrow = nrow(dfPiNew), ncol = nLF,
            dimnames = list(dfPiNew$sample, NULL)),
          simplify = FALSE)
      }
    }
    CalcLF <- FALSE
  }


  # Calculate latent factors
  if (CalcLF) {

    for (r in seq_len(mod_nr)) {

      post_alpha <- lapply(post, function(c) c$Alpha[[r]])

      # Save postEta to file and load it from predict_latent_factor. This helps
      # to avoid the unnecessary copying of the postEta object to all cores
      postEta_file <- fs::path(
        temp_dir, paste0(model_name, "_r", r, "_postEta.qs2"))

      if (isFALSE(ecokit::check_data(postEta_file, warning = FALSE))) {
        ecokit::cat_time("Save postEta to file", level = 1L, verbose = verbose)
        postEta <- lapply(post, function(c) c$Eta[[r]])
        ecokit::save_as(object = postEta, out_path = postEta_file)
        rm(postEta, envir = environment())
        invisible(gc())
      }

      # Save post to file and load it later
      if (r == mod_nr) {
        post_file <- fs::path(temp_dir, paste0(model_name, "_post.qs2"))

        if (isFALSE(ecokit::check_data(post_file, warning = FALSE))) {
          # free some memory
          post <- lapply(post, function(x) {
            x$Eta <- x$Psi <- x$V <- x$Delta <- x$Gamma <- x$rho <- NULL
            x
          })
          ecokit::cat_time("Save post to file", level = 1L, verbose = verbose)
          ecokit::save_as(object = post, out_path = post_file)
        }
        rm(post, envir = environment())
        invisible(gc())
      }

      if (r == 1) {
        ecokit::cat_time(
          "LF prediction using `predict_latent_factor`",
          level = 1L, verbose = verbose)
      }

      predPostEta[[r]] <- IASDT.R::predict_latent_factor(
        units_pred = levels(dfPiNew[, r]),
        units_model = levels(mod_dfPi[, r]),
        post_eta = postEta_file, post_alpha = post_alpha, lf_rl = rL[[r]],
        n_cores_lf = n_cores_lf, strategy = strategy,
        future_max_size = future_max_size, lf_temp_cleanup = lf_temp_cleanup,
        lf_out_file = lf_out_file, lf_return = lf_return, lf_check = lf_check,
        lf_commands_only = lf_commands_only, temp_dir = temp_dir,
        model_name = model_name, use_tf = use_tf, tf_environ = tf_environ,
        tf_use_single = tf_use_single)

      rm(postEta_file, envir = environment())

      if (lf_commands_only) {
        return(invisible(NULL))
      }

      rowNames <- rownames(predPostEta[[r]][[1]])
      PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
    }

  } else {

    if (spatial_model &&
        (is.null(prediction_type) || prediction_type == "c")) {
      ecokit::cat_time(
        "Loading LF prediction from disk", level = 1L, verbose = verbose)
      ecokit::cat_time(
        lf_input_file, level = 2L, cat_timestamp = FALSE, verbose = verbose)

      for (r in seq_len(mod_nr)) {
        predPostEta[[r]] <- ecokit::load_as(lf_input_file[[r]])
        rowNames <- rownames(predPostEta[[r]][[1]])
        PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
      }
    }

    if (isFALSE(spatial_model)) {
      for (r in seq_len(mod_nr)) {
        rowNames <- rownames(predPostEta[[r]][[1]])
        PiNew[, r] <- fastmatch::fmatch(dfPiNew[, r], rowNames)
      }
    }
  }

  if (mod_nr > 0) {
    ppEta <- simplify2array(predPostEta)
  } else {
    ppEta <- matrix(list(), predN, 0)
  }

  # Only predicting latent factor for new sites or for response curves at
  # median site
  if (lf_only) {
    return(lf_out_file)
  }

  # free some memory
  rm(predPostEta, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Predicting", verbose = verbose)

  if (!exists("post")) {
    ecokit::cat_time("Loading post from disk", level = 1L, verbose = verbose)
    post <- ecokit::load_as(post_file)
  }

  # Read model object from disk
  model_obj <- ecokit::load_as(model_file_small)

  # prediction data for response curves
  if (!is.null(prediction_type)) {
    ecokit::cat_time(
      "Predicting data for response curve (sequentially)",
      level = 1L, verbose = verbose)

    preds <- lapply(
      seq_len(predN),
      function(pN, ...) {
        get1prediction(
          model_obj, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
          ppEta[pN, ], PiNew, dfPiNew, nyNew, expected, mcmcStep)
      })

    ecokit::cat_diff(
      init_time = .start_time, prefix = "Prediction was finished in ",
      level = 1L, verbose = verbose)
    return(preds)
  }

  # Save ppEta / post as small chunks
  ecokit::cat_time(
    "Save ppEta / post as small chunks", level = 1L, verbose = verbose)
  chunk_size <- 25
  chunk_ids <- ceiling(seq_along(post) / chunk_size)
  chunks <- purrr::map_chr(
    .x = seq_len(max(chunk_ids)),
    .f = ~ {
      IDs <- which(chunk_ids == .x)
      Ch <- list(ppEta = ppEta[IDs], post = post[IDs])
      chunk_file <- fs::path(
        temp_dir, paste0(model_name, "_preds_ch", .x, ".qs2"))
      ecokit::save_as(object = Ch, out_path = chunk_file)
      return(chunk_file)
    })

  rm(chunk_ids, post, ppEta, envir = environment())
  invisible(gc())

  seeds <- sample.int(.Machine$integer.max, predN)

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, length(chunks)), show_log = FALSE,
      future_max_size = future_max_size, strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time("Making predictions in parallel", level = 1L)
  pred <- future.apply::future_lapply(
    X = seq_len(length(chunks)),
    FUN = function(chunk) {

      chunk_file <- chunks[chunk]
      Ch <- ecokit::load_as(chunk_file)
      ppEta <- Ch$ppEta
      post <- Ch$post
      rm(Ch, envir = environment())
      Seed <- (chunk - 1) * chunk_size
      Seed <- seq(Seed + 1, Seed + chunk_size)
      Seed <- seeds[Seed]

      pred_chunk <- purrr::map(
        .x = seq_len(chunk_size),
        .f = function(pN) {
          get1prediction(
            object = model_obj, X = X, XRRR = XRRR, Yc = Yc,
            Loff = Loff, rL = rL, rLPar = rLPar, sam = post[[pN]],
            predPostEta = ppEta[pN], PiNew = PiNew,
            dfPiNew = dfPiNew, nyNew = nyNew, expected = expected,
            mcmcStep = mcmcStep, seed = Seed[pN])
        })

      # species richness
      chunk_sr <- simplify2array(lapply(X = pred_chunk, FUN = rowSums)) %>%
        float::fl()
      dimnames(chunk_sr) <- NULL
      chunk_sr_file <- fs::path(
        temp_dir, paste0("pred_", model_name, "_ch", chunk, "_sr.qs2"))
      ecokit::save_as(object = chunk_sr, out_path = chunk_sr_file)

      rm(chunk_sr, envir = environment())

      # species predictions
      chunk_sp <- purrr::map_dfr(
        .x = seq_len(length(model_obj$spNames)),
        .f = function(sp) {

          sp_data <- purrr::map(pred_chunk, ~ .x[, sp], ncol = 1) %>%
            simplify2array() %>%
            float::fl()
          dimnames(sp_data) <- NULL

          chunk_sp_file <- fs::path(
            temp_dir,
            paste0("pred_", model_name, "_ch", chunk, "_taxon", sp, ".qs2"))

          ecokit::save_as(object = sp_data, out_path = chunk_sp_file)

          cbind.data.frame(
            chunk = chunk, sp = sp, ias_id = model_obj$spNames[sp],
            chunk_sp_file = chunk_sp_file)
        }) %>%
        dplyr::bind_rows(
          tibble::tibble(
            chunk = chunk, sp = 0, ias_id = "sr",
            chunk_sp_file = chunk_sr_file),
          .)

      rm(pred_chunk, envir = environment())

      invisible(gc())
      chunk_sp
    },
    future.seed = TRUE,
    future.globals = c(
      "model_obj", "X", "XRRR", "Yc", "Loff", "rL", "rLPar", "PiNew",
      "dfPiNew", "nyNew", "expected", "mcmcStep", "seeds", "chunk_size",
      "chunks", "temp_dir", "model_name", "get1prediction"),
    future.packages = pkg_to_export)

  pred <- tibble::tibble(dplyr::bind_rows(pred))

  # # ..................................................................... ###

  ecokit::cat_time(
    "Summarizing prediction outputs / Evaluation",
    level = 1L, verbose = verbose)

  evaluation_data <- dplyr::select(pred, -chunk) %>%
    dplyr::group_nest(sp, ias_id) %>%
    dplyr::mutate(data = purrr::map(data, unlist))

  evaluation_data <- future.apply::future_lapply(
    X = seq_len(nrow(evaluation_data)),
    FUN = function(ID) {

      sp <- evaluation_data$sp[[ID]]
      if (sp == 0) {
        sp_2 <- "sr"
      } else {
        sp_2 <- paste0("taxon", sp)
      }
      ias_id <- evaluation_data$ias_id[[ID]]
      data <- as.vector(evaluation_data$data[[ID]])

      species_data <- purrr::map(data, ecokit::load_as) %>%
        do.call(cbind, .) %>%
        as.double()

      # Mean prediction
      species_data_mean <- Rfast::rowmeans(species_data)

      # standard deviation of prediction
      species_data_sd <- Rfast::rowVars(species_data, std = TRUE)

      # Coefficient of variation
      species_data_mean0 <- species_data_mean
      # Replace very small mean values with reasonably small number to avoid
      # overflow warning
      species_data_mean0[species_data_mean0 < 1e-8] <- 1e-8
      species_data_cov <- species_data_sd / species_data_mean0

      rm(species_data, envir = environment())
      invisible(gc())

      if (is.null(pred_xy)) {
        pred_xy <- model_obj$rL$sample$s
      }

      pred_summary <- tibble::tibble(
        as.data.frame(pred_xy),
        Mean = species_data_mean, SD = species_data_sd,
        Cov = species_data_cov) %>%
        stats::setNames(
          c(
            "x", "y", paste0(ias_id, "_mean"),
            paste0(ias_id, "_sd"), paste0(ias_id, "_cov"))) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

      pred_summary_file <- fs::path(
        pred_directory, paste0("pred_", model_name, "_", sp_2, ".qs2"))

      ecokit::save_as(object = pred_summary, out_path = pred_summary_file)

      if (evaluate && sp_2 != "sr") {
        if (is.null(pred_pa)) {
          PresAbs <- model_obj$Y[, sp]
        } else {
          PresAbs <- pred_pa[, sp]
        }

        if (length(unique(PresAbs)) == 2) {
          # Calculate evaluation metrics if there are two both presence and
          # absence info in (testing) data
          RMSE <- caret::RMSE(PresAbs, species_data_mean)
          mean_at_pres <- mean(species_data_mean[which(PresAbs == 1)])
          mean_at_abs <- mean(species_data_mean[which(PresAbs == 0)])
          TjurR2 <- mean_at_pres - mean_at_abs
          AUC <- pROC::auc(
            response = PresAbs, predictor = species_data_mean,
            levels = c(0, 1), direction = "<") %>%
            as.numeric()
          Boyce <- ecospat::ecospat.boyce(
            fit = species_data_mean, obs = species_data_mean[PresAbs == 1],
            PEplot = FALSE)$cor
        } else {
          RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
        }
      } else {
        RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
      }

      tibble::tibble(
        sp = sp, ias_id = ias_id, path_pred = pred_summary_file,
        RMSE = RMSE, AUC = AUC, Boyce = Boyce, TjurR2 = TjurR2)

    },
    future.seed = TRUE, future.packages = pkg_to_export,
    future.globals = c(
      "evaluation_data", "evaluate", "pred_directory", "model_name",
      "pred_pa", "pred_xy"))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
    future::plan("sequential", gc = TRUE)
  }

  invisible(gc())

  # # ..................................................................... ###

  # Save predictions for all species in a single file
  ecokit::cat_time(
    "Save predictions for all species in a single file", verbose = verbose)

  evaluation_data <- dplyr::bind_rows(evaluation_data)

  predictions <- dplyr::select(evaluation_data, path_pred, sp, ias_id) %>%
    dplyr::mutate(
      sp_data = purrr::map(
        .x = path_pred,
        .f = ~ {
          ecokit::load_as(.x) %>%
            tidyr::pivot_longer(
              cols = tidyselect::starts_with(c("sp_", "sr")),
              names_to = "species", values_to = "prediction")
        })) %>%
    dplyr::pull(sp_data) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_wider(names_from = "species", values_from = "prediction") %>%
    dplyr::relocate(gtools::mixedsort(names(.))) %>%
    dplyr::select(
      x, y, geometry, sr_mean, sr_sd, sr_cov, tidyselect::everything()) %>%
    dplyr::mutate(model_name = model_name, .before = "sr_mean")

  pred_file <- fs::path(
    pred_directory,
    paste0(
      "prediction_",
      stringr::str_remove(model_name, "_clamping|_no_clamping"), ".qs2"))

  ecokit::save_as(object = predictions, out_path = pred_file)

  if (temp_cleanup) {
    try(fs::file_delete(evaluation_data$path_pred), silent = TRUE)
  }

  ecokit::cat_time("Predictions were saved", level = 1L, verbose = verbose)
  ecokit::cat_time(
    pred_file, level = 2L, cat_timestamp = FALSE, verbose = verbose)

  if (evaluate) {
    if (is.null(evaluation_name)) {
      eval_path <- fs::path(
        evaluation_directory,
        paste0(
          "eval_", stringr::str_remove(model_name, "_train|_current"), ".qs2"))
    } else {
      eval_path <- fs::path(
        evaluation_directory,
        paste0(
          "eval_", stringr::str_remove(model_name, "_train|_current"),
          "_", evaluation_name, ".qs2"))
    }

    evaluation_data <- dplyr::select(evaluation_data, -path_pred)
    ecokit::save_as(object = evaluation_data, out_path = eval_path)

    ecokit::cat_time(
      "Evaluation results were saved", level = 1L, verbose = verbose)
    ecokit::cat_time(
      fs::path(evaluation_directory, "eval_data.qs2"),
      level = 2L, cat_timestamp = FALSE, verbose = verbose)

  } else {
    eval_path <- NULL
  }

  # # ..................................................................... ###

  if (is.null(lf_out_file)) {
    lf_path <- lf_input_file
  } else {
    lf_path <- lf_out_file
  }

  # # ..................................................................... ###

  # Clean up
  if (temp_cleanup) {

    ecokit::cat_time(
      "Cleaning up temporary files", level = 1L, verbose = verbose)

    try(
      {
        file_pattern <- paste0("(pred_){0,}", model_name, ".+qs2")
        file_paths <- list.files(
          path = ecokit::normalize_path(temp_dir),
          pattern = file_pattern, full.names = TRUE)
        try(fs::file_delete(file_paths), silent = TRUE)
      },
      silent = TRUE)

    try(fs::file_delete(model_file_small), silent = TRUE)

  }

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Prediction was finished in ", verbose = verbose)

  return(
    tibble::tibble(
      pred_path = pred_file, eval_path = eval_path, lf_path = lf_path))
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

  mod_nr <- object$nr
  LRan <- vector("list", mod_nr)
  Eta <- vector("list", mod_nr)

  for (r in seq_len(mod_nr)) {
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
    for (r in seq_len(mod_nr)) {
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
