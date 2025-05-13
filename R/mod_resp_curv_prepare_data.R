## |------------------------------------------------------------------------| #
# resp_curv_prepare_data ----
## |------------------------------------------------------------------------| #

#' Prepare and plot response curve data for Hmsc models
#'
#' The `RespCurv_*()` functions process and visualise response curves for Hmsc
#' models. They support parallel computation and optionally return processed
#' data. There are four functions in this group:
#' - `resp_curv_prepare_data()`: Prepares response curve data for analysis
#' - `resp_curv_plot_species()`: Generates response curve plots for
#' individual species
#' - `resp_curv_plot_species_all()`: Generates response curves for all
#' species together in a single plot
#' - `resp_curv_plot_SR()`: Plots response curves for species richness.
#' @param path_model Character. Path to the file containing the fitted Hmsc
#'   model.
#' @param n_grid Integer. Number of points along the gradient for continuous
#'   focal variables. Higher values result in smoother curves. Default: 50. See
#'   [Hmsc::constructGradient] for details.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to 8L for all functions, except for `resp_curv_plot_species`,
#'   in which it defaults to 20L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "future::sequential", "future::multisession",
#'   "future::multicore", and "future::cluster". Defaults to
#'   `"future::multicore"` (`"future::multisession"` on Windows). See
#'   [future::plan()] and [ecokit::set_parallel()] for details.
#' @param return_data Logical. If `TRUE`, the function returns processed data as
#'   an R object. Default: `FALSE`.
#' @param probabilities Numeric vector. Quantiles to calculate in response curve
#'   predictions. Default: `c(0.025, 0.5, 0.975)`. See [stats::quantile] for
#'   details.
#' @param model_dir Character. Path to the root directory containing fitted
#'   models. The function reads data from the `RespCurv_DT` subdirectory, which
#'   is created by `resp_curv_prepare_data`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param plotting_alpha Numeric. Opacity level for response curve lines (0 =
#'   fully transparent, 1 = fully opaque). Default: 0.3.
#' @export
#' @inheritParams predict_hmsc
#' @rdname response_curves
#' @name response_curves
#' @order 1
#' @author Ahmed El-Gabbas

resp_curv_prepare_data <- function(
    path_model = NULL, n_grid = 50L, n_cores = 8L,
    strategy = "future::multicore", return_data = FALSE,
    probabilities = c(0.025, 0.5, 0.975), use_TF = TRUE, TF_environ = NULL,
    TF_use_single = FALSE, LF_n_cores = n_cores, LF_check = FALSE,
    LF_temp_cleanup = TRUE, LF_commands_only = FALSE, temp_dir = "TEMP_Pred",
    temp_cleanup = TRUE, verbose = TRUE) {

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }
  if (!is.numeric(LF_n_cores) || length(LF_n_cores) != 1 || LF_n_cores <= 0) {
    ecokit::stop_ctx(
      "LF_n_cores must be a single positive integer.", LF_n_cores = LF_n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "future::sequential") {
    n_cores <- LF_n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c(
    "future::sequential", "future::multisession", "future::multicore",
    "future::cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  # # ..................................................................... ###

  if (isFALSE(verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .start_time <- lubridate::now(tzone = "CET")

  if (is.null(path_model)) {
    ecokit::stop_ctx(
      "`path_model` cannot be NULL", path_model = path_model,
      include_backtrace = TRUE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ResCurvDT <- Variable <- RC_DT_Name <- SampleID <- Species <- SR <- MM <-
    NFV <- RC_DT_Path_Orig <- VarName <- RC_DT_Path_Prob <-
    RC_DT_Path_SR <- Coords <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  ecokit::cat_time("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(.x = AllArgs, .f = get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("path_model", "temp_dir"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "LF_n_cores", "n_grid", "probabilities"))
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "use_TF")
  rm(AllArgs, envir = environment())

  if (!is.numeric(n_cores) || n_cores < 1) {
    ecokit::stop_ctx(
      "`n_cores` must be greater than 0", n_cores = n_cores,
      include_backtrace = TRUE)
  }
  if (!is.numeric(LF_n_cores) || LF_n_cores < 1) {
    ecokit::stop_ctx(
      "`LF_n_cores` must be greater than 0", LF_n_cores = LF_n_cores,
      include_backtrace = TRUE)
  }
  if (any(probabilities > 1) || any(probabilities < 0)) {
    ecokit::stop_ctx(
      "`probabilities` must be between 0 and 1", probabilities = probabilities,
      include_backtrace = TRUE)
  }

  probabilities <- sort(probabilities)

  # # ..................................................................... ###

  # Loading model object ------

  ecokit::cat_time("Loading model object")
  if (file.exists(path_model)) {
    Model <- ecokit::load_as(path_model)
    if (!inherits(Model, "Hmsc")) {
      ecokit::stop_ctx(
        "Model object is not of class 'hmsc'", class_model = class(Model),
        include_backtrace = TRUE)
    }
  } else {
    ecokit::stop_ctx(
      "The model file does not exist or is not a `.RData` or `.qs2` file.",
      path_model = path_model, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # PrepRCData -------

  PrepRCData <- function(ID, File_LF) {

    Variable <- ResCurvDT$VarName[[ID]]
    RC_DT_Name <- ResCurvDT$RC_DT_Name[[ID]]
    Coords <- ResCurvDT$Coords[[ID]]
    NFV <- ResCurvDT$NFV[[ID]]

    # Path for original prediction values
    RC_DT_Path_Orig <- ResCurvDT$RC_DT_Path_Orig[[ID]]
    # Path for plotting data: probability of occurrence
    RC_DT_Path_Prob <- ResCurvDT$RC_DT_Path_Prob[[ID]]
    RC_DT_Path_Prob_Samples <- stringr::str_replace(
      RC_DT_Path_Prob, ".qs2$", "_Samples.qs2")

    # Path for plotting data: Species richness
    RC_DT_Path_SR <- ResCurvDT$RC_DT_Path_SR[[ID]]
    RC_DT_Path_SR_Samples <- stringr::str_replace(
      RC_DT_Path_SR, ".qs2$", "_Samples.qs2")

    OutputTibble <- tibble::tibble(
      Variable = Variable, NFV = NFV, Coords = Coords,
      RC_Path_Orig = RC_DT_Path_Orig,
      RC_Path_Prob = RC_DT_Path_Prob,
      RC_DT_Path_Prob_Samples = RC_DT_Path_Prob_Samples,
      RC_Path_SR = RC_DT_Path_SR,
      RC_DT_Path_SR_Samples = RC_DT_Path_SR_Samples)

    OutFilesExists <- c(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) %>%
      file.exists() %>%
      all()

    if (isFALSE(OutFilesExists)) {

      if (file.exists(RC_DT_Path_Orig)) {

        RC_Data_Orig <- ecokit::load_as(RC_DT_Path_Orig)
        gradient <- RC_Data_Orig$gradient
        XVals <- gradient$XDataNew[, Variable]
        Preds <- RC_Data_Orig$Preds
        Pred_SR <- RC_Data_Orig$Pred_SR

      } else {

        Model <- ecokit::load_as(path_model)

        # constructGradient
        gradient <- Hmsc::constructGradient(
          hM = Model, focalVariable = Variable, non.focalVariables = NFV,
          ngrid = n_grid, coordinates = list(sample = Coords))

        # Values of the current predictor
        XVals <- gradient$XDataNew[, Variable]

        rm(Model, envir = environment())

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        # Predicting probability of occurrence
        Preds <- IASDT.R::predict_hmsc(
          path_model = path_model, gradient = gradient, expected = TRUE,
          n_cores = 1, strategy = strategy, model_name = paste0("RC_", Coords),
          prediction_type = Coords, use_TF = use_TF, TF_environ = TF_environ,
          LF_inputFile = File_LF, LF_n_cores = 1, LF_check = LF_check,
          LF_temp_cleanup = LF_temp_cleanup, LF_commands_only = FALSE,
          TF_use_single = TF_use_single, temp_dir = temp_dir,
          temp_cleanup = temp_cleanup, verbose = FALSE)

        # Species richness
        Pred_SR <- abind::abind(lapply(Preds, rowSums), along = 2)

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # Save gradient and original prediction values
        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        RC_Data_Orig <- list(
          Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
          gradient = gradient, Preds = Preds, Pred_SR = Pred_SR)

        ecokit::save_as(
          object = RC_Data_Orig, object_name = paste0(RC_DT_Name, "_Orig"),
          out_path = RC_DT_Path_Orig)
      }

      rm(RC_Data_Orig, envir = environment())
      invisible(gc())

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: probability of occurrence
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Model <- ecokit::load_as(path_model)

      RC_Data_Prob <- purrr::map_dfr(
        .x = seq_len(length(Preds)),
        .f = function(Sample) {
          tibble::as_tibble(Preds[[Sample]]) %>%
            dplyr::mutate(XVals = XVals, SampleID = Sample)
        }) %>%
        tidyr::pivot_longer(
          cols = c(-XVals, -SampleID),
          names_to = "Species", values_to = "Pred") %>%
        dplyr::arrange(Species, XVals, SampleID) %>%
        tidyr::nest(SamplesData = -Species) %>%
        dplyr::mutate(
          PlotData_Quant = purrr::map(
            .x = SamplesData,
            .f = ~ {
              dplyr::reframe(
                .x, Pred = stats::quantile(Pred, probabilities),
                Quantile = probabilities, .by = XVals)
            }),

          # Values at observed presence and absences
          Observed_PA = purrr::map(
            .x = Species,
            .f = ~ tibble::tibble(
              XVals = Model$XData[, Variable], Pred = Model$Y[, .x])),

          # Positive trend probability
          PositiveTrendProb = purrr::map_dbl(
            .x = SamplesData,
            .f = ~ {
              dplyr::group_by(.x, SampleID) %>%
                dplyr::reframe(MM = dplyr::last(Pred) > dplyr::first(Pred)) %>%
                dplyr::pull(MM) %>%
                mean()
            }),
          Variable = Variable, NFV = ResCurvDT$NFV[[ID]], .before = 1)

      # Save data
      RC_Data_Prob_Samples <- RC_Data_Prob
      ecokit::save_as(
        object = RC_Data_Prob_Samples,
        object_name = paste0(RC_DT_Name, "_Prob_Samples"),
        out_path = RC_DT_Path_Prob_Samples)

      RC_Data_Prob <- dplyr::select(RC_Data_Prob, -SamplesData)
      ecokit::save_as(
        object = RC_Data_Prob, object_name = paste0(RC_DT_Name, "_Prob"),
        out_path = RC_DT_Path_Prob)


      rm(Preds, RC_Data_Prob, RC_Data_Prob_Samples, envir = environment())

      # CHECK
      # rm(SamplesData, envir = environment())


      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: Species richness
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      # predicted species richness
      SamplesData <- purrr::map_dfr(
        .x = seq_len(ncol(Pred_SR)),
        .f = ~ tibble::tibble(
          XVals = XVals, SampleID = .x, SR = Pred_SR[, .x])) %>%
        dplyr::arrange(XVals, SampleID)

      # Quantiles of species richness
      RC_Data_SR_Quant <- dplyr::reframe(
        SamplesData,
        SR = stats::quantile(SR, probabilities),
        Quantile = probabilities, .by = XVals)

      # Trend of the species richness
      SR_PositiveTrendProb <- SamplesData %>%
        dplyr::group_by(SampleID) %>%
        dplyr::reframe(MM = dplyr::last(SR) > dplyr::first(SR)) %>%
        dplyr::pull(MM) %>%
        mean()

      # Values at observed species richness
      Observed_SR <- tibble::tibble(
        XVals = Model$XData[, Variable], Pred = rowSums(Model$Y, na.rm = TRUE))

      # Save species richness data
      RC_Data_SR <- list(
        Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
        RC_Data_SR_Quant = RC_Data_SR_Quant, Observed_SR = Observed_SR,
        SR_PositiveTrendProb = SR_PositiveTrendProb)

      ecokit::save_as(
        object = RC_Data_SR, object_name = paste0(RC_DT_Name, "_SR"),
        out_path = RC_DT_Path_SR)

      RC_Data_SR_Samples <- list(
        Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
        RC_Data_SR = RC_Data_SR,
        RC_Data_SR_Quant = RC_Data_SR_Quant, Observed_SR = Observed_SR,
        SR_PositiveTrendProb = SR_PositiveTrendProb)

      ecokit::save_as(
        object = RC_Data_SR_Samples,
        object_name = paste0(RC_DT_Name, "_SR_Samples"),
        out_path = RC_DT_Path_SR_Samples)


      rm(
        RC_Data_SR, RC_Data_SR_Quant, Observed_SR,
        RC_Data_SR_Samples, SR_PositiveTrendProb, envir = environment())

      # CHECK
      # rm(Pred_SR, envir = environment())
      invisible(gc())
    }

    invisible(gc())

    return(OutputTibble)
  }

  # # ..................................................................... ###

  # Prepare response curve data -------

  Path_RC <- fs::path(dirname(dirname(path_model)), "Model_Postprocessing")
  Path_RC_DT <- fs::path(Path_RC, "RespCurv_DT")
  fs::dir_create(Path_RC_DT)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract names of the variables
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Extract names of the variables")
  ModelVars <- stringr::str_split(
    as.character(Model$XFormula)[2], "\\+", simplify = TRUE) %>%
    stringr::str_trim()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Prediction variants
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # `Coords`: Value of the `coordinates` argument of the `constructGradient`
  # function. `coordinates = "c"` for mean of coordinates (default);
  # `coordinates = "i"` for infinite coordinates without effect of spatial
  # dependence.
  #
  # NFV: Value of the `non.focalVariables` argument of `constructGradient`.
  # non.focalVariables = 1 sets the values of the non-focal variable to the most
  # likely value (defined as expected value for covariates, mode for factors).
  # non.focalVariables = 2 sets the values of the non-focal variable to most
  # likely value, given the value of focal variable, based on a linear
  # relationship. non.focalVariables = 3 fixes to the value given

  ResCurvDT <- tidyr::expand_grid(
    Variable = ModelVars, Coords = c("c", "i"), NFV = c(1, 2)) %>%
    dplyr::mutate(
      VarName = purrr::map_chr(
        .x = Variable, .f = stringr::str_remove_all,
        pattern = "stats::poly\\(|, degree = 2, raw = TRUE\\)"),
      RC_DT_Name = paste0("RC_", VarName, "_coord_", Coords, "_NFV", NFV),
      RC_DT_Path_Orig = fs::path(Path_RC_DT, paste0(RC_DT_Name, "_Orig.qs2")),
      RC_DT_Path_Prob = fs::path(Path_RC_DT, paste0(RC_DT_Name, "_Prob.qs2")),
      RC_DT_Path_SR = fs::path(Path_RC_DT, paste0(RC_DT_Name, "_SR.qs2")),
      FileExists = purrr::pmap_lgl(
        .l = list(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR),
        .f = function(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) {
          c(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) %>%
            file.exists() %>%
            all()
        }))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Checking file existence
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  MissingRows <- sum(!ResCurvDT$FileExists)
  File_LF <- fs::path(Path_RC_DT, "ResCurv_LF.qs2")

  if (MissingRows == 0) {

    ecokit::cat_time(
      "All response curve data files were already available on disk",
      level = 1L)
    ResCurvDT <- purrr::map_dfr(
      .x = seq_len(nrow(ResCurvDT)), .f = PrepRCData, File_LF = File_LF)

  } else {

    if (any(ResCurvDT$FileExists)) {
      ecokit::cat_time(
        paste0(
          "Some response curve data files (", MissingRows, " of ",
          length(ResCurvDT$FileExists), ") were missing"),
        level = 1L)
    } else {
      ecokit::cat_time(
        paste0(
          "All response curve data (", MissingRows, ") need to be prepared"),
        level = 1L)
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Get LF prediction for the model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if (isFALSE(ecokit::check_data(File_LF, warning = FALSE))) {

      ecokit::info_chunk(
        message = "Get LF prediction at mean coordinates", cat_date = FALSE,
        cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

      ecokit::cat_time("Create gradient")
      Gradient_c <- Hmsc::constructGradient(
        hM = Model, focalVariable = ResCurvDT$Variable[1],
        non.focalVariables = 1, ngrid = 20, coordinates = list(sample = "c"))

      # The `Model` object is distributed twice to cores when available on the
      # function environment. Here, I delete the Model object and it will be
      # loaded later after when using `predict_hmsc` function.
      rm(Model, envir = environment())
      invisible(gc())

      ecokit::cat_time("Predicting LF")
      Model_LF <- IASDT.R::predict_hmsc(
        path_model = path_model, gradient = Gradient_c, expected = TRUE,
        n_cores = n_cores, strategy = strategy, temp_dir = temp_dir,
        temp_cleanup = temp_cleanup, model_name = "RC_c",
        prediction_type = "c", use_TF = use_TF, TF_environ = TF_environ,
        LF_out_file = File_LF, LF_n_cores = LF_n_cores, LF_check = LF_check,
        LF_return = FALSE, LF_only = TRUE, LF_temp_cleanup = LF_temp_cleanup,
        LF_commands_only = LF_commands_only, TF_use_single = TF_use_single,
        verbose = verbose, pred_directory = temp_dir)

      if (LF_commands_only) {
        return(invisible(NULL))
      }

      rm(Model_LF, Gradient_c, envir = environment())
      invisible(gc())

    } else {
      ecokit::cat_time(
        paste0(
          "LF prediction will be loaded from available file: \n   >>>  ",
          File_LF))
    }


    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare working in parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::info_chunk(
      message = "Prepare response curve data", cat_date = FALSE,
      cat_red = TRUE, cat_bold = TRUE, cat_timestamp = FALSE)

    n_cores <- max(min(n_cores, MissingRows), 1)

    if (n_cores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("future::sequential", gc = TRUE))
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare response curve data in parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ecokit::cat_time("Prepare response curve data in parallel")

    if (strategy == "future::multicore") {
      pkg_to_export <- NULL
    } else {
      pkg_to_export <- c("dplyr", "purrr", "tidyr", "abind", "Hmsc", "parallel")
    }

    ResCurvDT <- future.apply::future_lapply(
      X = seq_len(nrow(ResCurvDT)),
      FUN = PrepRCData, File_LF = File_LF, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c(
        "ResCurvDT", "path_model", "PrepRCData", "n_grid", "probabilities",
        "File_LF", "use_TF", "TF_environ", "temp_dir", "LF_check",
        "LF_commands_only")) %>%
      dplyr::bind_rows()

    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
      future::plan("future::sequential", gc = TRUE)
    }

    invisible(gc())
  }

  # # ..................................................................... ###

  ecokit::cat_time("Saving data to desk")
  save(ResCurvDT, file = fs::path(Path_RC_DT, "ResCurvDT.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Preparing response curve data took ")

  if (return_data) {
    return(ResCurvDT)
  } else {
    return(invisible(NULL))
  }
}
