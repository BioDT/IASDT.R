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
#'   continuous focal variables. Defaults to 50.
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 15.
#' @param ReturnData Logical indicating whether the processed response curve
#'   data should be returned as an R object. Defaults to `FALSE`.
#' @seealso RespCurv_PlotSp RespCurv_PlotSR
#' @return Depending on the value of `ReturnData`, either returns response curve
#'   data or `NULL` invisibly.
#' @export
#' @name RespCurv_PrepData

RespCurv_PrepData <- function(
    Path_Model = NULL, ngrid = 50, NCores = 15, ReturnData = FALSE) {

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be NULL", call. = FALSE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ResCurvDT <- Variable <- RC_DT_Name <- SampleID <- Species <-
    SR <- MM <- PlotData <- NFV <- Coords <- RC_DT_Path_Orig <-
    RC_DT_Path_Prob <- RC_DT_Path_SR <- NULL

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Check input arguments ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(.x = AllArgs, .f = ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "Path_Model")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("NCores", "ngrid"))
  rm(AllArgs)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Loading model object ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading model object")
  if (file.exists(Path_Model) &&
      stringr::str_detect(Path_Model, ".+.RData$")) {
    Model <- IASDT.R::LoadAs(Path_Model)
  } else {
    stop("Path of the model object was not found", call. = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # PrepRCData_Int -------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  PrepRCData_Int <- function(ID) {

    Variable <- ResCurvDT$Variable[[ID]]
    Coords <- ResCurvDT$Coords[[ID]]
    NFV <- ResCurvDT$NFV[[ID]]
    RC_DT_Name <- ResCurvDT$RC_DT_Name[[ID]]

    # original prediction values
    RC_DT_Path_Orig <- ResCurvDT$RC_DT_Path_Orig[[ID]]
    # plotting data: probability of occurrence
    RC_DT_Path_Prob <- ResCurvDT$RC_DT_Path_Prob[[ID]]
    # plotting data: Species richness
    RC_DT_Path_SR <- ResCurvDT$RC_DT_Path_SR[[ID]]

    OutputTibble <- tibble::tibble(
      Variable = Variable, Coords = Coords, NFV = NFV,
      RC_Path_Orig = RC_DT_Path_Orig, RC_Path_Prob = RC_DT_Path_Prob,
      RC_Path_SR = RC_DT_Path_SR)

    OutFilesExists <- c(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) %>%
      file.exists() %>%
      all()

    if (isFALSE(OutFilesExists)) {

      if (file.exists(RC_DT_Path_Orig)) {
        RC_Data_Orig <- IASDT.R::LoadAs(RC_DT_Path_Orig)
        Gradient <- RC_Data_Orig$Gradient
        XVals <- Gradient$XDataNew[, Variable]
        Pred_Expected <- RC_Data_Orig$Pred_Expected
        Pred_SR <- RC_Data_Orig$Pred_SR

      } else {

        # +++++++++++++++++++++++++++++++++
        # constructGradient
        # +++++++++++++++++++++++++++++++++
        Gradient <- Hmsc::constructGradient(
          hM = Model, focalVariable = Variable, non.focalVariables = NFV,
          ngrid = ngrid, coordinates = list(sample = Coords))

        XVals <- Gradient$XDataNew[, Variable]

        # +++++++++++++++++++++++++++++++++
        # Predicting
        # +++++++++++++++++++++++++++++++++

        # probability of occurrence
        Pred_Expected <- stats::predict(
          object = Model, Gradient = Gradient, expected = TRUE)

        # Species richness
        Pred_SR <- stats::predict(
          object = Model, Gradient = Gradient, expected = FALSE)

        # +++++++++++++++++++++++++++++++++
        # Save gradient and original prediction values
        # +++++++++++++++++++++++++++++++++

        RC_Data_Orig <- list(
          Variable = Variable, Coords = Coords, NFV = NFV,
          Gradient = Gradient, Pred_Expected = Pred_Expected, Pred_SR = Pred_SR)
        IASDT.R::SaveAs(
          InObj = RC_Data_Orig, OutObj = paste0(RC_DT_Name,  "_Orig"),
          OutPath = RC_DT_Path_Orig)
      }

      rm(RC_Data_Orig)

      # +++++++++++++++++++++++++++++++++
      # Prepare plotting data: probability of occurrence
      # +++++++++++++++++++++++++++++++++

      Probabilities <- c(0.025, 0.5, 0.975) # quantiles to be calculated

      RC_Data_Prob <- purrr::map_dfr(
        .x = seq_len(length(Pred_Expected)),
        .f = function(Sample) {
          tibble::as_tibble(Pred_Expected[[Sample]]) %>%
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
            .f = ~{
              dplyr::reframe(
                .x,
                Pred = stats::quantile(Pred, Probabilities), Quantile = Probabilities,
                .by	= XVals)
            }),

          # Values at observed presence and absences
          Observed_PA = purrr::map(
            .x = Species,
            .f = ~tibble::tibble(
              XVals = Model$XData[, Variable], Pred = Model$Y[, .x])),

          # Positive trend probability
          PositiveTrendProb = purrr::map_dbl(
            .x = PlotData,
            .f = ~{
              dplyr::group_by(.x, SampleID) %>%
                dplyr::reframe(MM = dplyr::last(Pred) > dplyr::first(Pred)) %>%
                dplyr::pull(MM) %>%
                mean()
            })) %>%
        dplyr::mutate(
          Variable = Variable, Coords = Coords, NFV = NFV, .before = 1)

      # Save data
      IASDT.R::SaveAs(
        InObj = RC_Data_Prob, OutObj = paste0(RC_DT_Name,  "_Prob"),
        OutPath = RC_DT_Path_Prob)

      rm(Pred_Expected, RC_Data_Prob)

      # +++++++++++++++++++++++++++++++++
      # Prepare plotting data: Species richness
      # +++++++++++++++++++++++++++++++++

      # predicted species richness
      RC_Data_SR <- purrr::map(
        .x = seq_len(length(Pred_SR)),
        .f = ~tibble::tibble(
          XVals = XVals, SampleID = .x, SR = rowSums(Pred_SR[[.x]]))) %>%
        dplyr::bind_rows() %>%
        dplyr::arrange(XVals, SampleID)

      # Quantiles of species richness
      RC_Data_SR_Quant <- dplyr::reframe(
        RC_Data_SR, SR = stats::quantile(SR, Probabilities),
        Quantile = Probabilities, .by	= XVals)

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
        Variable = Variable, Coords = Coords, NFV = NFV,
        RC_Data_SR = RC_Data_SR,
        RC_Data_SR_Quant = RC_Data_SR_Quant,
        Observed_SR = Observed_SR,
        SR_PositiveTrendProb = SR_PositiveTrendProb)

      IASDT.R::SaveAs(
        InObj = RC_Data_SR, OutObj = paste0(RC_DT_Name,  "_SR"),
        OutPath = RC_DT_Path_SR)

      rm(RC_Data_SR, RC_Data_SR_Quant, Observed_SR, SR_PositiveTrendProb)
    }
    return(OutputTibble)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Prepare response curve data -------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare response curve data")

  Path_ResCurve <- dirname(dirname(Path_Model)) %>%
    file.path("Model_Postprocessing")
  Path_RespCurvDT <- file.path(Path_ResCurve, "RespCurv_DT")
  fs::dir_create(Path_RespCurvDT)

  # +++++++++++++++++++++++++++++++++
  # Extract names of the variables
  # +++++++++++++++++++++++++++++++++
  IASDT.R::CatTime("Extract names of the variables", Level = 1)
  ModelVars <- stringr::str_split(
    as.character(Model$XFormula)[2], "\\+", simplify = TRUE) %>%
    stringr::str_trim()

  # +++++++++++++++++++++++++++++++++
  # Predictions variants
  # +++++++++++++++++++++++++++++++++

  # Coords
  # Value of the `coordinates` argument of the `constructGradient` function
  # coordinates = "c" for mean of coordinates (default)
  # coordinates = "i" for infinite coordinates without effect of
  # spatial dependence

  # NFV
  # Value of the `non.focalVariables` argument of `constructGradient`
  # non.focalVariables = 1 sets the values of the non-focal variable to
  # the most likely value (defined as expected value for covariates, mode
  # for factors)
  # non.focalVariables = 2 sets the values of the non-focal variable to most
  # likely value, given the value of focal variable, based on a linear
  # relationship
  # non.focalVariables = 3 fixes to the value given

  ResCurvDT <- tidyr::expand_grid(
    Variable = ModelVars, Coords = c("c", "i"), NFV = c(1, 2)) %>%
    dplyr::mutate(
      RC_DT_Name = paste0("RC_", Variable, "_NFV_", NFV, "_Coords_", Coords),
      RC_DT_Path_Orig = file.path(
        Path_RespCurvDT, paste0(RC_DT_Name,  "_Orig.RData")),
      RC_DT_Path_Prob = file.path(
        Path_RespCurvDT, paste0(RC_DT_Name,  "_Prob.RData")),
      RC_DT_Path_SR = file.path(
        Path_RespCurvDT, paste0(RC_DT_Name,  "_SR.RData")),
      FileExists = purrr::pmap_lgl(
        .l = list(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR),
        .f = function(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) {
          c(RC_DT_Path_Orig, RC_DT_Path_Prob, RC_DT_Path_SR) %>%
            file.exists() %>%
            all()
        }))

  # +++++++++++++++++++++++++++++++++
  # Checking file existence
  # +++++++++++++++++++++++++++++++++

  MissingRows <- sum(!ResCurvDT$FileExists)

  if (MissingRows == 0) {
    IASDT.R::CatTime(
      "All response curve data files were already available on disk", Level = 1)
    ResCurvDT <- purrr::map(
      .x = seq_len(nrow(ResCurvDT)), .f = PrepRCData_Int) %>%
      dplyr::bind_rows()
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

    # +++++++++++++++++++++++++++++++++
    # Prepare working on parallel
    # +++++++++++++++++++++++++++++++++

    NCores <- max(min(NCores, MissingRows), 1)

    IASDT.R::CatTime(
      paste0("Prepare working on parallel, using ", NCores, " cores"),
      Level = 2)

    withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)
    
    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)

    # +++++++++++++++++++++++++++++++++
    # Prepare response curve data on parallel
    # +++++++++++++++++++++++++++++++++

    IASDT.R::CatTime("Prepare response curve data on parallel", Level = 2)

    ResCurvDT <- future.apply::future_lapply(
      X = seq_len(nrow(ResCurvDT)),
      FUN = PrepRCData_Int,
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = c("dplyr", "purrr", "tidyr"),
      future.globals = c("ResCurvDT", "Model")) %>%
      dplyr::bind_rows()

    future::plan("sequential")
  }

  IASDT.R::CatTime("Saving data to desk")
  save(ResCurvDT, file = file.path(Path_RespCurvDT, "ResCurvDT.RData"))

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  if (ReturnData) {
    return(ResCurvDT)
  } else {
    return(invisible(NULL))
  }
}
