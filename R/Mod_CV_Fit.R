## |------------------------------------------------------------------------| #
# Mod_CV_Fit ----
## |------------------------------------------------------------------------| #

#' Prepare cross-Validated Hmsc models for HPC fitting
#'
#' This function prepares cross-validated Hmsc models for fitting using HPC. It
#' handles data preparation, model initialization, and generation of SLURM
#' commands.
#'
#' @param Model Character. Path to a saved model file (`*.qs2`).
#' @param CVName Character vector specifying the name of the column(s) in the
#'   model input data (see [IASDT.R::Mod_PrepData] and [IASDT.R::GetCV]) to be
#'   used to cross-validate the models. The function allows the possibility of
#'   using more than one way of assigning grid cells into cross-validation
#'   folders. If multiple names are provided, separate cross-validation models
#'   will be fitted for each column. Currently, there are three cross-validation
#'   strategies, created using the [IASDT.R::Mod_PrepData]: `CV_SAC`, `CV_Dist`,
#'   and `CV_Large` (see [IASDT.R::GetCV]). Defaults to `c("CV_Dist",
#'   "CV_Large")`.
#' @param Partitions A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `CVName` argument. If the `Partitions`
#'   vector is provided, the label used in the output files will be `CV_Custom`.
#' @param InitPar a named list of parameter values used for initialization of
#'   MCMC states. See [Hmsc::computePredictedValues] for more information.
#'   Default: `NULL`.
#' @param JobName String specifying the name of the submitted job(s) for SLURM.
#'   Default: `CV_Models`.
#' @param Updater a named list, specifying which conditional updaters should be
#'   omitted.  See [Hmsc::computePredictedValues] for more information. Defaults
#'   to `list(Gamma2 = FALSE, GammaEta = FALSE)` to disable the following
#'   warnings: `setting updater$Gamma2=FALSE due to specified phylogeny matrix`
#'   and `setting updater$GammaEta=FALSE: not implemented for spatial methods
#'   'GPP' and 'NNGP'`.
#' @param AlignPost boolean flag indicating whether the posterior of each chains
#'   should be aligned. See [Hmsc::computePredictedValues] for more information.
#'   Default: `TRUE`.
#' @param ... Additional arguments passed to the [IASDT.R::Mod_SLURM] function.
#' @details The function copies part of the [Hmsc::computePredictedValues]
#'   function, which currently does not support performing cross-validation
#'   using Hmsc-HPC.
#' @author Ahmed El-Gabbas
#' @inheritParams Mod_SLURM
#' @inheritParams Mod_PrepData
#' @inheritParams Mod_Prep4HPC
#' @export
#' @name Mod_CV_Fit

Mod_CV_Fit <- function(
    Model = NULL, CVName = c("CV_Dist", "CV_Large"), Partitions = NULL,
    EnvFile = ".env", InitPar = NULL, JobName = "CV_Models",
    Updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    AlignPost = TRUE, ToJSON = FALSE, FromHPC = TRUE, PrepSLURM = TRUE,
    MemPerCpu = NULL, Time = NULL, Path_Hmsc = NULL, Precision = 64, ...) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- Path_ModInit_rds <- CV <- ModName <- partition <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input parameters -----
  IASDT.R::CatTime("Check input parameters")

  NullVarsNames <- c("Model", "Path_Hmsc", "MemPerCpu", "Time", "EnvFile")
  NullVars <- which(purrr::map_lgl(.x = NullVarsNames, .f = ~ is.null(get(.x))))

  if (length(NullVars) > 0) {
    NullVarsNames[NullVars]
    stop(
      paste0(
        paste0(NullVarsNames[NullVars], collapse = ", "), " cannot be NULL"),
      call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    .x = AllArgs,
    .f = function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "Model", "JobName", "EnvFile", "Time", "MemPerCpu", "Path_Hmsc")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  if (!(Precision %in% c(32, 64))) {
    stop(
      paste0("Precision should be either of 32 or 64, not ", Precision),
      call. = FALSE)
  }

  if (!file.exists(Model)) {
    stop("Model path does not exist.", call. = FALSE)
  }

  rm(AllArgs, NullVarsNames, NullVars, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Loading model -----
  IASDT.R::CatTime("Loading model")

  # Path of the cross-validation folder
  Path_CV <- file.path(dirname(dirname(Model)), "Model_Fitting_CV")

  # Path of the model input data
  Path_ModelData <- list.files(
    path = dirname(dirname(Model)),
    pattern = "^ModDT_.+_subset.RData", full.names = TRUE)
  if (length(Path_ModelData) != 1) {
    stop(
      "There should be exactly one file matches model input data",
      call. = FALSE)
  }

  # Loading full model object
  Model_Full <- IASDT.R::LoadAs(Model)

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

  if (is.null(Partitions)) {
    # if custom partitions is not provided, extract the CV column(s) available
    # in the modelling data

    CV_Data <- IASDT.R::LoadAs(Path_ModelData)$DT_CV

    if (all(CVName %in% names(CV_Data))) {
      # Extract CV folds from the CV column(s)
      Partitions <- purrr::map(CVName, ~ dplyr::pull(CV_Data, .x))
      names(Partitions) <- stringr::str_remove(CVName, "^CV_")
    } else {
      # if any of the column names does not exist, stop the function
      MissingCV <- CVName[CVName %in% names(CV_Data) == FALSE]
      stop(
        paste0(
          "`Partitions` was not defined (NULL) and column(s) for CV folds ",
          paste(MissingCV, collapse = " + "),
          " can not be found in species data"
        ),
        call. = FALSE
      )
    }
  } else {
    # If Partitions is provided directly to the function, use "CV_Custom" as CV
    # name
    CVName <- "CV_Custom"
  }

  # Check the length of CV data equals the number of sampling units in the model
  if (any(sapply(Partitions, length) != Model_Full$ny)) {
    stop(
      paste0(
        "Partitions parameter must be a vector of the same length of the ",
        "sampling  units of the the full model"),
      call. = FALSE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Verbose ----

  # Verbose value should be extracted directly from the full model object. It
  # seems that models fitted with Hmsc-HPC have NULL verbose value (although not
  # tested). If the `verbose` value stored in the full model object is NULL, the
  # following assigns a value of 1000 to it (using the null coalescing operator
  # `%||%`).

  verbose <- Model_Full$verbose %||% 1000

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Number of chains ----
  NChains <- length(Model_Full$postList)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare cross-validated initial models -----
  IASDT.R::CatTime("Prepare cross-validated initial models")

  # Load coordinates from the full model
  Coords <- IASDT.R::LoadAs(Path_ModelData)$DT_xy

  CV_DT <- tibble::tibble(
    partition = Partitions, CVName = names(Partitions)) %>%
    dplyr::mutate(

      # number of cross-validation folders for each CV type
      nfolds = purrr::map_int(.x = partition, .f = ~ length(unique(.x))),

      # prepare data for each cross-validation strategy
      CV_Info = purrr::pmap(
        .l = list(partition, CVName, nfolds),
        .f = function(partition, CVName, nfolds) {

          IASDT.R::CatTime(paste0("Cross-validation Type: ", CVName), Level = 1)

          CV_DT0 <- purrr::map_dfr(
            .x = seq_len(nfolds),
            .f = function(k) {

              IASDT.R::CatTime(paste0("Fold ", k, "/", nfolds), Level = 2)

              # # |||||||||||||||||||||||
              # The following is adapted from Hmsc::computePredictedValues()
              # # |||||||||||||||||||||||

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


              # initiate sampling and save initial models
              Model_CV <- Hmsc::sampleMcmc(
                hM = Model_CV, samples = Model_Full$samples,
                thin = Model_Full$thin, transient = Model_Full$transient,
                adaptNf = Model_Full$adaptNf, initPar = InitPar,
                nChains = NChains, updater = Updater, verbose = verbose,
                alignPost = AlignPost, engine = "HPC")

              if (ToJSON) {
                Model_CV <- jsonify::to_json(Model_CV)
              }

              # Save model input as rds file
              Path_ModInit_rds <- file.path(
                Path_Init, paste0("InitMod_", CVName, "_k", k, ".rds"))
              saveRDS(Model_CV, file = Path_ModInit_rds)

              Path_ModFitted <- file.path(
                Path_Fitted, paste0("Model_", CVName, "_k", k, ".RData"))

              dfPi <- droplevels(Model_Full$dfPi[val, , drop = FALSE])
              valCoords <- Coords[val, ]
              rownames(dfPi) <- rownames(valCoords) <- dfPi$sample

              tibble::tibble(
                Path_ModFull = Model,
                CV = k,
                Path_ModInit = Path_ModInit,
                Path_ModInit_rds = Path_ModInit_rds,
                Path_ModFitted = Path_ModFitted,
                val = list(val),
                valCoords = list(valCoords),
                XRRRVal = list(XRRRVal),
                dfPi = list(dfPi)) %>%
                return()
            })
          return(CV_DT0)
        })) %>%
    dplyr::select(-"partition") %>%
    tidyr::unnest("CV_Info")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare model fitting commands -----
  IASDT.R::CatTime("Prepare model fitting commands")

  CV_DT <- CV_DT %>%
    dplyr::mutate(ModName = paste0(CVName, CV), .after = CV) %>%
    dplyr::mutate(
      ModelPrep = purrr::pmap(
        .l = list(Path_ModInit_rds, CV, ModName),
        .f = function(Path_ModInit_rds, CV, ModName) {

          CV_Out <- purrr::map_dfr(
            .x = seq_len(NChains),
            .f = function(Chain) {

              # Turn off scientific notation
              withr::local_options(list(scipen = 999))

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
                " --input ", IASDT.R::NormalizePath(Path_ModInit_rds),
                " --output ", IASDT.R::NormalizePath(Path_Post),
                " --samples ", Model_Full$samples,
                " --transient ", Model_Full$transient,
                " --thin ", Model_Full$thin,
                " --verbose ", verbose,
                " --chain ", (Chain - 1),
                " --fp ", Precision,
                " >& ", IASDT.R::NormalizePath(Path_ModProg))

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

  # Prepare SLURM script -----

  if (PrepSLURM) {
    IASDT.R::CatTime("Prepare SLURM script")

    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::Mod_SLURM(
      ModelDir = Path_CV, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC, Path_Hmsc = Path_Hmsc,
      Path_SLURM_Out = dirname(Path_Post), SLURM_Prefix = "CV_Bash_Fit", ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save summary data to disk -----
  IASDT.R::CatTime("Save summary data to disk")
  save(CV_DT, file = file.path(Path_CV, "CV_DT.RData"))

  return(invisible(NULL))
}
