## |------------------------------------------------------------------------| #
# Mod_CV_Fit ----
## |------------------------------------------------------------------------| #

#' Prepare Cross-Validated Hmsc Models for HPC Fitting
#'
#' This function prepares cross-validated Hmsc models for fitting using HPC. It
#' handles data preparation, model initialization, and generation of SLURM
#' commands.
#'
#' @param Model Either a path to a saved model file (character) or an Hmsc model
#'   object. If a path is provided, the model is loaded from the file.
#' @param ModelData Character. Path to input data used to fit the model.
#' @param CVName Character vector specifying the name of the column(s) in the
#'   model input data (see [IASDT.R::Mod_PrepData] and [IASDT.R::GetCV]) to be
#'   used to cross-validate the models. The function allows the possibility of
#'   using more than one way of assigning grid cells into cross-validation
#'   folders. If multiple names are provided, separate cross-validation models
#'   will be fitted for each column. Currently, there are three cross-validation
#'   strategies, created using the [IASDT.R::Mod_PrepData]: `CV_SAC`, `CV_Dist`,
#'   and `CV_Large` (see [IASDT.R::GetCV]).
#' @param partitions A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `CVName` argument. If the `partitions`
#'   vector is provided, the label used in the output files will be `CV_Custom`.
#' @param Path_CV The directory path where cross-validation models and outputs
#'   will be stored. If `Model` argument is a character vector, it will be
#'   estimated from the model path. If `Model` is an Hmsc model object, it has
#'   to be provided by the user, otherwise the function will give an error.
#' @param initPar a named list of parameter values used for initialization of
#'   MCMC states. See [Hmsc::computePredictedValues] for more information.
#'   Default: `NULL`.
#' @param JobName String specifying the name of the submitted job(s) for SLURM.
#'   Default: `CV_Models`.
#' @param updater a named list, specifying which conditional updaters should be
#'   omitted.  See [Hmsc::computePredictedValues] for more information. Defaults
#'   to `list(Gamma2 = FALSE, GammaEta = FALSE)` to disable the following
#'   warnings: `setting updater$Gamma2=FALSE due to specified phylogeny matrix`
#'   and `setting updater$GammaEta=FALSE: not implemented for spatial methods
#'   'GPP' and 'NNGP'`.
#' @param alignPost boolean flag indicating whether the posterior of each chains
#'   should be aligned. See [Hmsc::computePredictedValues] for more information.
#'   Default: `TRUE`.
#' @param ... Additional arguments passed to the [IASDT.R::Mod_SLURM] function.
#' @details The function copies part of the [Hmsc::computePredictedValues]
#'   function, which currently does not support performing cross-validation
#'   using Hmsc-HPC. Although it is possible to [make some
#'   changes](https://github.com/aniskhan25/hmsc-hpc/issues/14) to the
#'   [Hmsc::computePredictedValues] to make it possible to use Hmsc-HPC, this
#'   does not help as our intention is to use the Hmsc-HPC using GPU.
#'
#'   The function reads the following environment variables:
#'   - **`DP_R_Path_Python`** for Python path on LUMI.
#'   project on LUMI.
#' @author Ahmed El-Gabbas
#' @inheritParams Mod_SLURM
#' @inheritParams Mod_PrepData
#' @inheritParams Mod_Prep4HPC
#' @export
#' @name Mod_CV_Fit

Mod_CV_Fit <- function(
    Model = NULL, ModelData = NULL,
    CVName = c("CV_Dist", "CV_Large"), partitions = NULL,
    Path_CV = NULL, EnvFile = ".env", initPar = NULL, JobName = "CV_Models",
    updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    alignPost = TRUE, ToJSON = FALSE, FromHPC = TRUE, PrepSLURM = TRUE,
    MemPerCpu = NULL, Time = NULL, Path_Hmsc = NULL, Precision = 64, ...) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- Path_ModInit_rds <- CV <- partition <- ModName <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input parameters -----
  IASDT.R::CatTime("Check input parameters")

  NullVarsNames <- c(
    "Model", "ModelData", "Path_Hmsc", "MemPerCpu", "Time", "EnvFile")
  NullVars <- which(purrr::map_lgl(.x = NullVarsNames, .f = ~ is.null(get(.x))))

  if (length(NullVars) > 0) {
    NullVarsNames[NullVars]
    stop(
      paste0(
        paste0(NullVarsNames[NullVars], collapse = ", "),
        " cannot be missing."),
      call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    .x = AllArgs,
    .f = function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "Model", "ModelData", "JobName", "EnvFile",
    "Time", "MemPerCpu", "Path_Hmsc")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  if (!(Precision %in% c(32, 64))) {
    stop(
      paste0("Precision should be either of 32 or 64, not ", Precision),
      call. = FALSE)
  }

  rm(AllArgs, NullVarsNames, NullVars, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Loading model -----
  IASDT.R::CatTime("Loading model")

  if (inherits(Model, "character")) {
    if (!file.exists(Model)) {
      stop("Model path does not exist.", call. = FALSE)
    }

    if (is.null(Path_CV)) {
      Path_CV <- file.path(dirname(dirname(Model)), "Model_Fitting_CV")
    }

    Model_Full <- IASDT.R::LoadAs(Model)
  } else {
    if (is.null(Path_CV)) {
      stop(
        "Path_CV cannot be empty if the Model is provided as Hmsc object",
        call. = FALSE)
    }
    Model_Full <- Model
    rm(Model, envir = environment())
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating paths -----
  IASDT.R::CatTime("Creating paths")

  Path_Init <- file.path(Path_CV, "Model_Init")
  Path_Fitted <- file.path(Path_CV, "Model_Fitted")
  Path_Post <- file.path(Path_CV, "Model_Posterior")
  Dir_Pred <- file.path(Path_CV, "Model_Prediction")
  fs::dir_create(c(Path_CV, Path_Init, Path_Fitted, Path_Post, Dir_Pred))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Cross-validation partitions ----
  IASDT.R::CatTime("Cross-validation partitions")

  if (is.null(partitions)) {
    # if custom partitions is not provided, extract the CV column(s) available
    # in the modelling data

    CV_Data <- IASDT.R::LoadAs(ModelData)$DT_CV

    if (all(CVName %in% names(CV_Data))) {
      # Extract CV folds from the CV column(s)
      partitions <- purrr::map(CVName, ~ dplyr::pull(CV_Data, .x))
      names(partitions) <- stringr::str_remove(CVName, "^CV_")
    } else {
      # if any of the column names does not exist, stop the function
      MissingCV <- CVName[CVName %in% names(CV_Data) == FALSE]
      stop(
        paste0(
          "`partitions` was not defined (NULL) and column(s) for CV folds ",
          paste(MissingCV, collapse = " + "),
          " can not be found in species data"
        ),
        call. = FALSE
      )
    }
  } else {
    # If partitions is provided directly to the function, use "CV_Custom" as CV
    # name
    CVName <- "CV_Custom"
  }

  # Check the length of CV data equals the number of sampling units in the model
  if (any(sapply(partitions, length) != Model_Full$ny)) {
    stop("partitions parameter must be a vector of length ny", call. = FALSE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Verbose ----

  # If not provided explicitly as model argument, the verbose value should be
  # read directly from the model object. It seems that models fitted with
  # Hmsc-HPC have NULL verbose value; however, I did not test this. if the
  # verbose value stored in the model object is NULL, the following assigns a
  # value of 1000 to it.

  verbose <- Model_Full$verbose %||% 1000

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Number of chains ----
  NChains <- length(Model_Full$postList)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare cross-validated initial models -----
  IASDT.R::CatTime("Prepare cross-validated initial models")

  # Load coordinates
  Coords <- IASDT.R::LoadAs(ModelData)$DT_xy

  CV_DT <- tibble::tibble(
    partition = partitions, CVName = names(partitions)) %>%
    dplyr::mutate(
      nfolds = purrr::map_int(.x = partition, .f = ~ length(unique(.x))),

      # Loop over all cross-validation partitions used
      CV_Info = purrr::pmap(
        .l = list(partition, CVName, nfolds),
        .f = function(partition, CVName, nfolds) {

          # prepare data for each cross-validation strategy

          CV_DT0 <- purrr::map_dfr(
            .x = seq_len(nfolds),
            .f = function(k) {

              # The following is adapted from Hmsc::computePredictedValues()

              train <- (partition != k)
              val <- (partition == k)

              switch(
                class(Model_Full$X)[1L],
                matrix = {
                  XTrain <- Model_Full$X[train, , drop = FALSE]
                },
                list = {
                  XTrain <- purrr::map(
                    .x = Model_Full$X, .f = ~ .x[train, , drop = FALSE])
                })

              if (Model_Full$ncRRR > 0) {
                XRRRTrain <- Model_Full$XRRR[train, , drop = FALSE]
                XRRRVal <- Model_Full$XRRR[val, , drop = FALSE]
              } else {
                XRRRTrain <- XRRRVal <- NULL
              }

              dfPi <- droplevels(Model_Full$dfPi[train, , drop = FALSE])

              YTrain <- Model_Full$Y[train, , drop = FALSE]
              rownames(YTrain) <- rownames(dfPi)

              # Initial models
              Model_CV <- Hmsc::Hmsc(
                Y = YTrain, Loff = Model_Full$Loff[train, , drop = FALSE],
                X = XTrain, XRRR = XRRRTrain, ncRRR = Model_Full$ncRRR,
                XSelect = Model_Full$XSelect, distr = Model_Full$distr,
                studyDesign = dfPi, Tr = Model_Full$Tr, C = Model_Full$C,
                ranLevels = Model_Full$rL) %>%
                Hmsc::setPriors(
                  V0 = Model_Full$V0, f0 = Model_Full$f0,
                  mGamma = Model_Full$mGamma, UGamma = Model_Full$UGamma,
                  aSigma = Model_Full$aSigma, bSigma = Model_Full$bSigma,
                  rhopw = Model_Full$rhowp)

              Model_CV$XFormula <- Model_Full$XFormula
              Model_CV$XInterceptInd <- Model_Full$XInterceptInd

              # function to remove scaled attributes
              RemAttr <- function(x) {
                attr(x, "scaled:center") <- NULL
                attr(x, "scaled:scale") <- NULL
                return(x)
              }

              Model_CV$YScalePar <- Model_Full$YScalePar
              Model_CV$YScaled <- scale(
                Model_CV$Y, Model_CV$YScalePar[1, ], Model_CV$YScalePar[2, ])

              Model_CV$XScalePar <- Model_Full$XScalePar
              Model_CV$XScaled <- scale(
                Model_CV$X, Model_CV$XScalePar[1, ], Model_CV$XScalePar[2, ])

              if (Model_CV$ncRRR > 0) {
                Model_CV$XRRRScalePar <- Model_Full$XRRRScalePar
                Model_CV$XRRRScaled <- scale(
                  Model_CV$XRRR,
                  Model_CV$XRRRScalePar[1, ], Model_CV$XRRRScalePar[2, ])
                Model_CV$XRRRScaled <- RemAttr(Model_CV$XRRRScaled)
              }

              Model_CV$TrInterceptInd <- Model_Full$TrInterceptInd
              Model_CV$TrScalePar <- Model_Full$TrScalePar
              Model_CV$TrScaled <- scale(
                Model_CV$Tr, Model_CV$TrScalePar[1, ], Model_CV$TrScalePar[2, ])
              Model_CV$TrScaled <- RemAttr(Model_CV$TrScaled)


              # Save unfitted model
              Path_ModInit <- file.path(
                Path_Init, paste0("InitMod_", CVName, "_k", k, ".RData"))
              IASDT.R::SaveAs(
                InObj = Model_CV,
                OutObj = paste0("InitMod_", CVName, "_k", k),
                OutPath = Path_ModInit)


              # initiate sampling and save initial models to
              Model_CV <- Hmsc::sampleMcmc(
                hM = Model_CV, samples = Model_Full$samples,
                thin = Model_Full$thin, transient = Model_Full$transient,
                adaptNf = Model_Full$adaptNf, initPar = initPar,
                nChains = NChains, updater = updater, verbose = verbose,
                alignPost = alignPost, engine = "HPC")

              if (ToJSON) {
                Model_CV <- jsonify::to_json(Model_CV)
              }

              # Save model input as rds file
              Path_ModInit_rds <- file.path(
                Path_Init, paste0("InitMod_", CVName, "_k", k, ".rds"))
              saveRDS(Model_CV, file = Path_ModInit_rds)

              Path_ModFull <- dplyr::if_else(
                inherits(Model, "character"), Model, NA_character_)

              Path_ModFitted <- file.path(
                Path_Fitted, paste0("Model_", CVName, "_k", k, ".RData"))

              dfPi <- droplevels(Model_Full$dfPi[val, , drop = FALSE])
              valCoords <- Coords[val, ]
              rownames(dfPi) <- rownames(valCoords) <- dfPi$sample

              tibble::tibble(
                Path_ModFull = Path_ModFull,
                CV = k,
                Path_ModInit = Path_ModInit,
                Path_ModInit_rds = Path_ModInit_rds,
                Path_ModFitted = Path_ModFitted,
                val = list(val),
                valCoords = list(valCoords),
                XRRRVal = list(XRRRVal),
                dfPi = list(dfPi)) %>%
                return()
            }
          )
          return(CV_DT0)
        })) %>%
    dplyr::select(-"partition") %>%
    tidyr::unnest("CV_Info")

  # Prepare fitting command for each model chain
  CV_DT <- CV_DT %>%
    dplyr::mutate(ModName = paste0(CVName, CV), .after = CV) %>%
    dplyr::mutate(
      ModelPrep = purrr::pmap(
        .l = list(Path_ModInit_rds, CV, ModName),
        .f = function(Path_ModInit_rds, CV, ModName) {

          CV_Out <- purrr::map_dfr(
            .x = seq_len(NChains),
            .f = function(Chain) {

              # Path to save the posterior of the combination of CV and chain
              Path_Post <- file.path(
                Path_Post,
                paste0("Mod_", ModName, "_Ch", Chain, "_post.rds"))

              # Path to save the progress of model fitting
              Path_ModProg <- stringr::str_replace_all(
                Path_Post, "post.rds$", "Progress.txt")

              # Model fitting command
              Command_HPC <- paste0(
                "/usr/bin/time -v python3 -m hmsc.run_gibbs_sampler",
                " --input ", shQuote(Path_ModInit_rds),
                " --output ", shQuote(Path_Post),
                " --samples ", Model_Full$samples,
                " --transient ", Model_Full$transient,
                " --thin ", Model_Full$thin,
                " --verbose ", verbose,
                " --chain ", (Chain - 1),
                " --fp ", Precision,
                " >& ", shQuote(Path_ModProg))

              # data to be returned for each combination of CV and Chain
              tibble::tibble(
                Path_Post = Path_Post, Path_ModProg = Path_ModProg,
                Command_HPC = Command_HPC, NSamples = Model_Full$samples,
                Transient = Model_Full$transient, Thin = Model_Full$thin) %>%
                return()
            })

          return(dplyr::summarise_all(CV_Out, list))

        })) %>%
    tidyr::unnest("ModelPrep") %>%
    dplyr::mutate_at(
      .vars = c("NSamples", "Transient", "Thin"),
      .funs = ~ as.integer(unique(unlist(.x))))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model fitting commands -----
  IASDT.R::CatTime("Save model fitting commands")

  CommandFile <- file.path(Path_CV, "Commands2Fit.txt")
  f <- file(CommandFile, open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(unlist(CV_DT$Command_HPC), sep = "\n", append = FALSE, file = f)
  close(f)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare SLURM -----

  if (PrepSLURM) {
    IASDT.R::CatTime("Prepare SLURM script")

    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::Mod_SLURM(
      ModelDir = Path_CV, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC, Path_Hmsc = Path_Hmsc,
      Path_SLURM_Out = Path_Post, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save summary data to disk -----
  IASDT.R::CatTime("Save summary data to disk")
  save(CV_DT, file = file.path(Path_CV, "CV_DT.RData"))

  return(invisible(NULL))
}
