## |------------------------------------------------------------------------| #
# RespCurv_PrepData ----
## |------------------------------------------------------------------------| #

#' Prepare and Process Response Curve Data for Hmsc Models
#'
#' This function prepares and processes data for generating response curves for
#' Hmsc models. It supports parallel processing and can return the processed
#' data.
#' @param Path_Model Character. Path file containing the model.
#' @param ngrid Integer specifying the number of points along the gradient for
#'   continuous focal variables. Defaults to 50. See [Hmsc::constructGradient]
#'   for more details.
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param ReturnData Logical. If `TRUE`, returns processed response curve data
#'   as an R object. Default is `FALSE`.
#' @param Probabilities quantiles to be calculated. Defaults to `c(0.025, 0.5,
#'   0.975)`. See [stats::quantile] for more details.
#' @seealso RespCurv_PlotSp RespCurv_PlotSR
#' @return Depending on the value of `ReturnData`, either returns response curve
#'   data or `NULL` invisibly.
#' @export
#' @inheritParams Predict_Hmsc
#' @name RespCurv_PrepData

RespCurv_PrepData <- function(
    Path_Model = NULL, ngrid = 50, NCores = 8, ReturnData = FALSE,
    Probabilities = c(0.025, 0.5, 0.975), UseTF = TRUE, TF_Environ = NULL,
    Temp_Dir = "TEMP2Pred", Verbose = TRUE) {

  # # ..................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be NULL", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ResCurvDT <- Variable <- RC_DT_Name <- SampleID <- Species <- SR <- MM <-
    NFV <- RC_DT_Path_Orig <- VarName <- RC_DT_Path_Prob <-
    RC_DT_Path_SR <- Coords <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(.x = AllArgs, .f = ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Path_Model", "TF_Environ"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "ngrid", "Probabilities"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "UseTF")
  rm(AllArgs, envir = environment())

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
  if (file.exists(Path_Model)) {
    Model <- IASDT.R::LoadAs(Path_Model)
    if (!inherits(Model, "Hmsc")) {
      stop("Model object is not of class 'hmsc'", call. = FALSE)
    }
  } else {
    stop(
      "The model file does not exist or is not a .RData file.", call. = FALSE)
  }

  # # ..................................................................... ###

  # PrepRCData -------

  PrepRCData <- function(ID) {

    Variable <- ResCurvDT$VarName[[ID]]
    RC_DT_Name <- ResCurvDT$RC_DT_Name[[ID]]
    Coords <- ResCurvDT$Coords[[ID]]
    NFV <- ResCurvDT$NFV[[ID]]

    # Path for original prediction values
    RC_DT_Path_Orig <- ResCurvDT$RC_DT_Path_Orig[[ID]]
    # Path for plotting data: probability of occurrence
    RC_DT_Path_Prob <- ResCurvDT$RC_DT_Path_Prob[[ID]]
    RC_DT_Path_Prob_Samples <- stringr::str_replace(
      RC_DT_Path_Prob, ".RData$", "_Samples.RData")

    # Path for plotting data: Species richness
    RC_DT_Path_SR <- ResCurvDT$RC_DT_Path_SR[[ID]]
    RC_DT_Path_SR_Samples <- stringr::str_replace(
      RC_DT_Path_SR, ".RData$", "_Samples.RData")

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

        RC_Data_Orig <- IASDT.R::LoadAs(RC_DT_Path_Orig)
        Gradient <- RC_Data_Orig$Gradient
        XVals <- Gradient$XDataNew[, Variable]
        Preds <- RC_Data_Orig$Preds
        Pred_SR <- RC_Data_Orig$Pred_SR

      } else {

        Model <- IASDT.R::LoadAs(Path_Model)

        # constructGradient
        Gradient <- Hmsc::constructGradient(
          hM = Model, focalVariable = Variable, non.focalVariables = NFV,
          ngrid = ngrid, coordinates = list(sample = Coords))

        # Values of the current predictor
        XVals <- Gradient$XDataNew[, Variable]

        rm(Model, envir = environment())

        # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        # Predicting probability of occurrence
        Preds <- IASDT.R::Predict_Hmsc(
          Path_Model = Path_Model, Gradient = Gradient, expected = TRUE, 
          NCores = 1,
          Model_Name = paste0("RC_", Coords), RC = Coords, UseTF = UseTF,
          TF_Environ = TF_Environ, LF_InputFile = File_LF, Verbose = FALSE)

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

      rm(RC_Data_Orig, envir = environment())

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Prepare plotting data: probability of occurrence
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Model <- IASDT.R::LoadAs(Path_Model)

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
                .x, Pred = stats::quantile(Pred, Probabilities),
                Quantile = Probabilities, .by = XVals)
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
            })) %>%
        dplyr::mutate(
          Variable = Variable, NFV = ResCurvDT$NFV[[ID]], .before = 1)

      # Save data
      RC_Data_Prob_Samples <- RC_Data_Prob
      IASDT.R::SaveAs(
        InObj = RC_Data_Prob_Samples,
        OutObj = paste0(RC_DT_Name, "_Prob_Samples"),
        OutPath = RC_DT_Path_Prob_Samples)

      RC_Data_Prob <- dplyr::select(RC_Data_Prob, -SamplesData)
      IASDT.R::SaveAs(
        InObj = RC_Data_Prob, OutObj = paste0(RC_DT_Name, "_Prob"),
        OutPath = RC_DT_Path_Prob)

      rm(Preds, RC_Data_Prob, RC_Data_Prob_Samples, envir = environment())

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
        SR = stats::quantile(SR, Probabilities),
        Quantile = Probabilities, .by = XVals)

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

      IASDT.R::SaveAs(
        InObj = RC_Data_SR, OutObj = paste0(RC_DT_Name, "_SR"),
        OutPath = RC_DT_Path_SR)

      RC_Data_SR_Samples <- list(
        Variable = Variable, NFV = ResCurvDT$NFV[[ID]],
        RC_Data_SR = RC_Data_SR,
        RC_Data_SR_Quant = RC_Data_SR_Quant, Observed_SR = Observed_SR,
        SR_PositiveTrendProb = SR_PositiveTrendProb)

      IASDT.R::SaveAs(
        InObj = RC_Data_SR_Samples, OutObj = paste0(RC_DT_Name, "_SR_Samples"),
        OutPath = RC_DT_Path_SR_Samples)


      rm(
        RC_Data_SR, RC_Data_SR_Quant, Observed_SR,
        RC_Data_SR_Samples, SR_PositiveTrendProb, envir = environment())
    }

    invisible(gc())

    return(OutputTibble)
  }

  # # ..................................................................... ###

  # Prepare response curve data -------

  Path_RC <- dirname(dirname(Path_Model)) %>%
    file.path("Model_Postprocessing")
  Path_RC_DT <- file.path(Path_RC, "RespCurv_DT")
  fs::dir_create(Path_RC_DT)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract names of the variables
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Extract names of the variables")
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
      RC_DT_Path_Orig = file.path(
        Path_RC_DT, paste0(RC_DT_Name, "_Orig.RData")),
      RC_DT_Path_Prob = file.path(
        Path_RC_DT, paste0(RC_DT_Name, "_Prob.RData")),
      RC_DT_Path_SR = file.path(Path_RC_DT, paste0(RC_DT_Name, "_SR.RData")),
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
    ResCurvDT <- purrr::map_dfr(.x = seq_len(nrow(ResCurvDT)), .f = PrepRCData)
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
    # Get LF prediction for the model
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    File_LF <- file.path(Path_RC_DT, "ResCurv_LF.RData")

    if (isFALSE(IASDT.R::CheckData(File_LF, warning = FALSE))) {

      IASDT.R::InfoChunk(
        Message = "Get LF prediction at mean coordinates",
        Date = FALSE, Extra2 = 1)

      IASDT.R::CatTime("Create gradient")
      Gradient_c <- Hmsc::constructGradient(
        hM = Model, focalVariable = ResCurvDT$Variable[1],
        non.focalVariables = 1, ngrid = 20, coordinates = list(sample = "c"))

      # The `Model` object is distributed twice to cores when available on the
      # function environment. Here, I delete the Model object and it will be
      # loaded later after when using `Predict_Hmsc` function.
      rm(Model, envir = environment())
      invisible(gc())

      Model_LF <- IASDT.R::Predict_Hmsc(
        Path_Model = Path_Model, Gradient = Gradient_c, expected = TRUE,
        NCores = NCores, Temp_Dir = Temp_Dir, Model_Name = "RC_c", RC = "c",
        UseTF = UseTF, TF_Environ = TF_Environ, LF_OutFile = File_LF,
        Verbose = Verbose)

      invisible(gc())
      rm(Model_LF, Gradient_c, envir = environment())

    } else {
      IASDT.R::CatTime(
        paste0(
          "LF prediction will be loaded from available file: \n   >>>  ",File_LF))
    }


    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Prepare working on parallel
    # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    IASDT.R::InfoChunk(
      Message = "Prepare response curve data", Date = FALSE, Extra2 = 1)

    NCores <- max(min(NCores, MissingRows), 1)

    IASDT.R::CatTime(
      paste0("Prepare working on parallel, using ", NCores, " cores"))

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

    IASDT.R::CatTime("Prepare response curve data on parallel")

    ResCurvDT <- future.apply::future_lapply(
      X = seq_len(nrow(ResCurvDT)),
      FUN = PrepRCData, future.seed = TRUE,
      future.packages = c(
        "dplyr", "purrr", "tidyr", "abind", "Hmsc", "parallel"),
      future.globals = c(
        "ResCurvDT", "Path_Model", "PrepRCData", "ngrid", "Probabilities",
        "File_LF", "UseTF", "TF_Environ")) %>%
      dplyr::bind_rows()

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }
    invisible(gc())
  }

  # # ..................................................................... ###

  IASDT.R::CatTime("Saving data to desk", ... = "\n")
  save(ResCurvDT, file = file.path(Path_RC_DT, "ResCurvDT.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Preparing response curve data took ")

  if (ReturnData) {
    return(ResCurvDT)
  } else {
    return(invisible(NULL))
  }
}
