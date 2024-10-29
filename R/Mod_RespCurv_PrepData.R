## |------------------------------------------------------------------------| #
# RespCurv_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare and Process Response Curve Data for Hmsc Models
#'
#' This function prepares and processes data for generating response curves for
#' Hmsc models. It supports parallel processing and can return the processed
#' data.
#' @param Path_Model String specifying the path to the .RData file containing
#'   the model to be used.
#' @param ngrid Integer specifying the number of points along the gradient for
#'   continuous focal variables. Defaults to 50. See [Hmsc::constructGradient]
#'   for more details.
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param ReturnData Logical indicating whether the processed response curve
#'   data should be returned as an R object. Defaults to `FALSE`.
#' @param predictEtaMean Logical; whether to predict the mean value of the
#'   latent variable. Defaults to `TRUE`. See `Hmsc:::predict.Hmsc` for more
#'   details.
#' @param Probabilities quantiles to be calculated. Defaults to `c(0.025, 0.5,
#'   0.975)`. See [stats::quantile] for more details.
#' @seealso RespCurv_PlotSp RespCurv_PlotSR
#' @return Depending on the value of `ReturnData`, either returns response curve
#'   data or `NULL` invisibly.
#' @export
#' @name RespCurv_PrepData

RespCurv_PrepData <- function(
    Path_Model = NULL, ngrid = 50, NCores = 8, ReturnData = FALSE,
    Probabilities = c(0.025, 0.5, 0.975), predictEtaMean = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be NULL", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ResCurvDT <- Variable <- RC_DT_Name <- SampleID <- Species <- SR <- MM <-
    PlotData <- NFV <- RC_DT_Path_Orig <- VarName <- RC_DT_Path_Prob <-
    RC_DT_Path_SR <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(.x = AllArgs, .f = ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "Path_Model")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "ngrid", "Probabilities"))
  rm(AllArgs)

  if (!is.numeric(NCores) || NCores < 1) {
    stop("NCores must be greater than 0", call. = FALSE)
  }
  if (any(Probabilities > 1) || any(Probabilities < 0)) {
    stop("Probabilities must be between 0 and 1", call. = FALSE)
  }

  Probabilities <- sort(Probabilities)

  # # ..................................................................... ###

  # Loading model object ------

  IASDT.R::CatTime("Loading model object")
  if (file.exists(Path_Model) &&
      stringr::str_detect(Path_Model, ".+.RData$")) {
    Model <- IASDT.R::LoadAs(Path_Model)
    if (!inherits(Model, "Hmsc")) {
      stop("Model object is not of class 'hmsc'", call. = FALSE)
    }
  } else {
    stop("Path of the model object was not found", call. = FALSE)
  }

  # # ..................................................................... ###

  # PrepRCData_Int -------

  PrepRCData_Int <- function(ID) {

    Variable <- ResCurvDT$VarName[[ID]]
    RC_DT_Name <- ResCurvDT$RC_DT_Name[[ID]]

    # Path for original prediction values
    RC_DT_Path_Orig <- ResCurvDT$RC_DT_Path_Orig[[ID]]
    # Path for plotting data: probability of occurrence
    RC_DT_Path_Prob <- ResCurvDT$RC_DT_Path_Prob[[ID]]
    # Path for plotting data: Species richness
    RC_DT_Path_SR <- ResCurvDT$RC_DT_Path_SR[[ID]]

    OutputTibble <- tibble::tibble(
      Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
      RC_Path_Orig = RC_DT_Path_Orig,
      RC_Path_Prob = RC_DT_Path_Prob, RC_Path_SR = RC_DT_Path_SR)

    OutFilesExists <- c(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) %>%
      file.exists() %>%
      all()

    if (isFALSE(OutFilesExists)) {
      if (file.exists(RC_DT_Path_Orig)) {
        RC_Data_Orig <- IASDT.R::LoadAs(RC_DT_Path_Orig)
        Gradient <- RC_Data_Orig$Gradient
        XVals <- Gradient$XDataNew[, Variable]
        Preds <- RC_Data_Orig$Preds
        Pred_SR <- RC_Data_Orig$Pred_SR
      } else {
        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        # constructGradient
        Gradient <- Hmsc::constructGradient(
          hM = Model, focalVariable = Variable,
          non.focalVariables = ResCurvDT$NFV[[ID]],
          ngrid = ngrid, coordinates = list(sample = "i"))

        # Values of the current predictor
        XVals <- Gradient$XDataNew[, Variable]

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        # Predicting probability of occurrence
        Preds <- stats::predict(
        # Preds <- predict.Hmsc(
          object = Model, Gradient = Gradient, nParallel = 1, expected = TRUE,
          predictEtaMean = predictEtaMean, RC = TRUE)

        # Species richness
        Pred_SR <- abind::abind(lapply(Preds, rowSums), along = 2)

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # Save gradient and original prediction values
        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        RC_Data_Orig <- list(
          Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
          Gradient = Gradient, Preds = Preds, Pred_SR = Pred_SR)
        IASDT.R::SaveAs(
          InObj = RC_Data_Orig, OutObj = paste0(RC_DT_Name, "_Orig"),
          OutPath = RC_DT_Path_Orig)
      }

      rm(RC_Data_Orig)

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: probability of occurrence
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
        tidyr::nest(PlotData = -Species) %>%
        dplyr::mutate(
          PlotData_Quant = purrr::map(
            .x = PlotData,
            .f = ~ {
              dplyr::reframe(
                .x,
                Pred = stats::quantile(Pred, Probabilities),
                Quantile = Probabilities, .by = XVals)
            }),

          # Values at observed presence and absences
          Observed_PA = purrr::map(
            .x = Species,
            .f = ~ tibble::tibble(
              XVals = Model$XData[, Variable], Pred = Model$Y[, .x])),

          # Positive trend probability
          PositiveTrendProb = purrr::map_dbl(
            .x = PlotData,
            .f = ~ {
              dplyr::group_by(.x, SampleID) %>%
                dplyr::reframe(MM = dplyr::last(Pred) > dplyr::first(Pred)) %>%
                dplyr::pull(MM) %>%
                mean()
            })) %>%
        dplyr::mutate(
          Variable = Variable, NFV = ResCurvDT$NFV[[ID]], .before = 1)

      # Save data
      IASDT.R::SaveAs(
        InObj = RC_Data_Prob, OutObj = paste0(RC_DT_Name, "_Prob"),
        OutPath = RC_DT_Path_Prob)

      rm(Preds, RC_Data_Prob)

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: Species richness
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      # predicted species richness
      RC_Data_SR <- purrr::map_dfr(
        .x = seq_len(ncol(Pred_SR)),
        .f = ~ tibble::tibble(
          XVals = XVals, SampleID = .x, SR = Pred_SR[, .x])) %>%
        dplyr::arrange(XVals, SampleID)

      # Quantiles of species richness
      RC_Data_SR_Quant <- dplyr::reframe(
        RC_Data_SR, SR = stats::quantile(SR, Probabilities),
        Quantile = Probabilities, .by = XVals)

      # Trend of the species richness
      SR_PositiveTrendProb <- RC_Data_SR %>%
        dplyr::group_by(SampleID) %>%
        dplyr::reframe(MM = dplyr::last(SR) > dplyr::first(SR)) %>%
        dplyr::pull(MM) %>%
        mean()

      # Values at observed species richness
      Observed_SR <- tibble::tibble(
        XVals = Model$XData[, Variable], Pred = rowSums(Model$Y, na.rm = TRUE))

      # Save species richness data
      RC_Data_SR <- list(
        Variable = Variable, NFV = ResCurvDT$NFV[[ID]], RC_Data_SR = RC_Data_SR,
        RC_Data_SR_Quant = RC_Data_SR_Quant, Observed_SR = Observed_SR,
        SR_PositiveTrendProb = SR_PositiveTrendProb)

      IASDT.R::SaveAs(
        InObj = RC_Data_SR, OutObj = paste0(RC_DT_Name, "_SR"),
        OutPath = RC_DT_Path_SR)

      rm(RC_Data_SR, RC_Data_SR_Quant, Observed_SR, SR_PositiveTrendProb)
    }
    return(OutputTibble)
  }

  # # ..................................................................... ###

  # Prepare response curve data -------

  IASDT.R::CatTime("Prepare response curve data")

  Path_ResCurve <- dirname(dirname(Path_Model)) %>%
    file.path("Model_Postprocessing")
  Path_RespCurvDT <- file.path(Path_ResCurve, "RespCurv_DT")
  fs::dir_create(Path_RespCurvDT)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract names of the variables
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  IASDT.R::CatTime("Extract names of the variables", Level = 1)
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
  # For `coordinates = "c"`, the `predictLatentFactor` will be used, which will
  # take too much time to run for large-scale studies. Here, I Only use
  # `coordinates = "i"` which is much faster after the explicit disabling the
  # use of predictLatentFactor in the updated `predict` function.
  #
  # NFV: Value of the `non.focalVariables` argument of `constructGradient`.
  # non.focalVariables = 1 sets the values of the non-focal variable to the most
  # likely value (defined as expected value for covariates, mode for factors).
  # non.focalVariables = 2 sets the values of the non-focal variable to most
  # likely value, given the value of focal variable, based on a linear
  # relationship. non.focalVariables = 3 fixes to the value given

  ResCurvDT <- tidyr::expand_grid(
    Variable = ModelVars, Coords = "i", NFV = c(1, 2)) %>%
    dplyr::mutate(
      VarName = purrr::map_chr(
        .x = Variable,
        .f = ~ {
          stringr::str_remove_all(
            .x ,
            "stats::poly\\(|, degree = 2, raw = TRUE\\)")
        }
      ),
    RC_DT_Name = paste0("RC_", VarName, "_NFV", NFV),
    RC_DT_Path_Orig = file.path(
      Path_RespCurvDT, paste0(RC_DT_Name, "_Orig.RData")),
    RC_DT_Path_Prob = file.path(
      Path_RespCurvDT, paste0(RC_DT_Name, "_Prob.RData")),
    RC_DT_Path_SR = file.path(
      Path_RespCurvDT, paste0(RC_DT_Name, "_SR.RData")),
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

  if (MissingRows == 0) {
    IASDT.R::CatTime(
      "All response curve data files were already available on disk",
      Level = 1)
    ResCurvDT <- purrr::map_dfr(
      .x = seq_len(nrow(ResCurvDT)), .f = PrepRCData_Int)
  } else {
    if (all(!ResCurvDT$FileExists)) {
      IASDT.R::CatTime(
        paste0(
          "All response curve data (", MissingRows, ") need to be prepared"),
        Level = 1)
    } else {
      IASDT.R::CatTime(
        paste0(
          "Some response curve data files (", MissingRows, " of ",
          length(ResCurvDT$FileExists), ") were missing"),
        Level = 1)
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare working on parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    NCores <- max(min(NCores, MissingRows), 1)

    IASDT.R::CatTime(
      paste0("Prepare working on parallel, using ", NCores, " cores"),
      Level = 2)

    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

    if (NCores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare response curve data on parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    IASDT.R::CatTime("Prepare response curve data on parallel", Level = 2)

    ResCurvDT <- future.apply::future_lapply(
      X = seq_len(nrow(ResCurvDT)),
      FUN = PrepRCData_Int,
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = c(
        "dplyr", "purrr", "tidyr", "abind", "Hmsc", "parallel"),
      future.globals = c(
        "ResCurvDT", "Model", "PrepRCData_Int", "ngrid", "Probabilities")) %>%
      dplyr::bind_rows()

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  IASDT.R::CatTime("Saving data to desk")
  save(ResCurvDT, file = file.path(Path_RespCurvDT, "ResCurvDT.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  if (ReturnData) {
    return(ResCurvDT)
  } else {
    return(invisible(NULL))
  }
}
