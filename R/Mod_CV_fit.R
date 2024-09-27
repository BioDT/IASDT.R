## |------------------------------------------------------------------------| #
# Mod_CV_fit ----
## |------------------------------------------------------------------------| #

#' Prepare Cross-Validated Hmsc Models for HPC Fitting
#'
#' This function prepares initial cross-validated Hmsc models in R for fitting
#' by Hmsc-HPC. It includes data preparation, model initialization, and command
#' generation for running models on HPC.
#' @param Model Either a path to a saved model file (character) or an Hmsc model
#'   object. If a path is provided, the model is loaded from the file.
#' @param ModelData Character. Path to input data used to fit the model.
#' @param CVNames Character vector specifying the name of the column(s) in the
#'   model input data (see [IASDT.R::Mod_PrepData] and [IASDT.R::GetCV]) to be
#'   used to cross-validate the models. The function allows the possibility of
#'   using more than one way of assigning grid cells into cross-validation
#'   folders. If multiple names are provided, separate cross-validation models
#'   will be fitted for each column. Currently, there are three cross-validation
#'   strategies, created using the [IASDT.R::Mod_PrepData]: `CV_SAC`, `CV_Dist`,
#'   and `CV_Large` (see [IASDT.R::GetCV]).
#' @param partition A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `CVNames` argument. If the `partition`
#'   vector
#'   is provided, the label used in the output files will be `CV_Custom`.
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
#' @name Mod_CV_fit

Mod_CV_fit <- function(
    Model = NULL, ModelData = NULL,
    CVNames = c("CV_SAC", "CV_Dist", "CV_Large"), partition = NULL,
    Path_CV = NULL, EnvFile = ".env", initPar = NULL, JobName = "CV_Models",
    updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    alignPost = TRUE, ToJSON = FALSE, FromHPC = TRUE, PrepSLURM = TRUE,
    MemPerCpu = NULL, Time = NULL, Path_Hmsc = NULL, Precision = 64, ...) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- CV_DT0 <- NULL

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

  rm(AllArgs, NullVarsNames, NullVars)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Loading model -----
  IASDT.R::CatTime("Loading model")

  if (inherits(Model, "character")) {
    if (!file.exists(Model)) {
      stop("Model path does not exist.", call. = FALSE)
    }

    if (is.null(Path_CV)) {
      Path_CV <- file.path(dirname(dirname(Model)), "Model_Fitting_HPC_CV")
    }

    ModFull <- IASDT.R::LoadAs(Model)
  } else {
    if (is.null(Path_CV)) {
      stop(
        "Path_CV cannot be empty if the Model is provided as Hmsc object",
        call. = FALSE)
    }
    ModFull <- Model
    rm(Model)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating paths -----
  IASDT.R::CatTime("Creating paths")

  Path_Init <- file.path(Path_CV, "InitMod")
  fs::dir_create(c(Path_CV, Path_Init))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Cross-validation partitions ----
  IASDT.R::CatTime("Cross-validation partitions")

  if (is.null(partition)) {
    # if custom partition is not provided, extract the CV column(s) available in
    # the modelling data

    CV_Data <- IASDT.R::LoadAs(ModelData)$DT_CV

    if (all(CVNames %in% names(CV_Data))) {
      # Extract CV folds from the CV column(s)
      partition <- purrr::map(CVNames, ~ dplyr::pull(CV_Data, .x))
      names(partition) <- CVNames
    } else {
      # if any of the column names does not exist, stop the function
      MissingCV <- CVNames[CVNames %in% names(CV_Data) == FALSE]
      stop(
        paste0(
          "`partition` was not defined (NULL) and column(s) for CV folds ",
          paste(MissingCV, collapse = " + "),
          " can not be found in species data"
        ),
        call. = FALSE
      )
    }
  } else {
    # If partition is provided directly to the function, use "CV_Custom" as CV
    # name
    CVNames <- "CV_Custom"
  }

  # Check the length of CV data equals the number of sampling units in the model
  if (any(sapply(partition, length) != ModFull$ny)) {
    stop("partition parameter must be a vector of length ny", call. = FALSE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Verbose ----

  # If not provided explicitly as model argument, the verbose value should be
  # read directly from the model object. It seems that models fitted with
  # Hmsc-HPC have NULL verbose value; however, I did not test this. if the
  # verbose value stored in the model object is NULL, the following assign a
  # value of 1000 to it.
  if (is.null(ModFull$verbose)) {
    verbose <- 1000
  } else {
    verbose <- ModFull$verbose
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Number of chains ----

  NChains <- length(ModFull$postList)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare cross-validated initial models -----
  IASDT.R::CatTime("Prepare cross-validated initial models")

  Coords <- IASDT.R::LoadAs(ModelData)$DT_xy

  CV_DT <- tibble::tibble(partition = partition, CVNames = names(partition)) %>%
    dplyr::mutate(
      nfolds = purrr::map_int(.x = partition, .f = ~ length(unique(.x))),

      # Loop over all cross-validation partitions used
      CV_DT0 = purrr::pmap(
        .l = list(partition, CVNames, nfolds),
        .f = function(CV_Fold, CV_Name, CV_Folds) {
          # prepare data for each cross-validation strategy

          CV_DT0 <- purrr::map_dfr(
            .x = seq_len(CV_Folds),
            .f = function(k) {

              # The following is adapted from Hmsc::computePredictedValues()

              train <- (CV_Fold != k)
              val <- (CV_Fold == k)
              valCoords <- Coords[val, ]

              dfPi <- as.data.frame(matrix(NA, sum(train), ModFull$nr))
              colnames(dfPi) <- ModFull$rLNames

              for (r in seq_len(ModFull$nr)) {
                dfPi[, r] <- factor(ModFull$dfPi[train, r])
              }

              LoffVal <- ModFull$LoffVal[val, , drop = FALSE]

              switch(
                class(ModFull$X)[1L],
                matrix = {
                  XTrain <- ModFull$X[train, , drop = FALSE]
                  XVal <- ModFull$X[val, , drop = FALSE]

                },
                list = {
                  XTrain <- purrr::map(
                    .x = ModFull$X, .f = ~ .x[train, , drop = FALSE])
                  XVal <- lapply(ModFull$X, function(a) a[val, , drop = FALSE])
                })

              if (ModFull$ncRRR > 0) {
                XRRRTrain <- ModFull$XRRR[train, , drop = FALSE]
                XRRRVal <- ModFull$XRRR[val, , drop = FALSE]
              } else {
                XRRRTrain <- XRRRVal <- NULL
              }

              # Initial models
              Mod_CV <- Hmsc::Hmsc(
                Y = ModFull$Y[train, , drop = FALSE],
                Loff = ModFull$Loff[train, , drop = FALSE],
                X = XTrain, XRRR = XRRRTrain, ncRRR = ModFull$ncRRR,
                XSelect = ModFull$XSelect, distr = ModFull$distr,
                studyDesign = dfPi, Tr = ModFull$Tr, C = ModFull$C,
                ranLevels = ModFull$rL) %>%
                Hmsc::setPriors(
                  V0 = ModFull$V0, f0 = ModFull$f0,
                  mGamma = ModFull$mGamma, UGamma = ModFull$UGamma,
                  aSigma = ModFull$aSigma, bSigma = ModFull$bSigma,
                  rhopw = ModFull$rhowp)

              Mod_CV$YScalePar <- ModFull$YScalePar
              Mod_CV$YScaled <- (
                Mod_CV$Y -
                  matrix(
                    Mod_CV$YScalePar[1, ], Mod_CV$ny, Mod_CV$ns,
                    byrow = TRUE)) /
                matrix(
                  Mod_CV$YScalePar[2, ], Mod_CV$ny, Mod_CV$ns, byrow = TRUE)
              Mod_CV$XInterceptInd <- ModFull$XInterceptInd
              Mod_CV$XScalePar <- ModFull$XScalePar

              switch(
                class(ModFull$X)[1L],
                matrix = {
                  Mod_CV$XScaled <- (
                    Mod_CV$X - matrix(
                      Mod_CV$XScalePar[1, ], Mod_CV$ny, Mod_CV$ncNRRR,
                      byrow = TRUE)) /
                    matrix(
                      Mod_CV$XScalePar[2, ], Mod_CV$ny, Mod_CV$ncNRRR,
                      byrow = TRUE)
                },
                list = {
                  Mod_CV$XScaled <- list()
                  for (zz in seq_len(length(Mod_CV$X))) {
                    Mod_CV$XScaled[[zz]] <- (
                      Mod_CV$X[[zz]] - matrix(
                        Mod_CV$XScalePar[1, ], Mod_CV$ny, Mod_CV$ncNRRR,
                        byrow = TRUE)) /
                      matrix(Mod_CV$XScalePar[2, ], Mod_CV$ny, Mod_CV$ncNRRR,
                             byrow = TRUE)
                  }
                }
              )

              if (Mod_CV$ncRRR > 0) {
                Mod_CV$XRRRScalePar <- ModFull$XRRRScalePar
                Mod_CV$XRRRScaled <- (
                  Mod_CV$XRRR - matrix(
                    Mod_CV$XRRRScalePar[1, ], Mod_CV$ny, Mod_CV$ncORRR,
                    byrow = TRUE)) /
                  matrix(
                    Mod_CV$XRRRScalePar[2, ], Mod_CV$ny, Mod_CV$ncORRR,
                    byrow = TRUE)
              }

              Mod_CV$TrInterceptInd <- ModFull$TrInterceptInd
              Mod_CV$TrScalePar <- ModFull$TrScalePar
              Mod_CV$TrScaled <- (
                Mod_CV$Tr - matrix(
                  Mod_CV$TrScalePar[1, ], Mod_CV$ns,
                  Mod_CV$nt, byrow = TRUE)) /
                matrix(
                  Mod_CV$TrScalePar[2, ], Mod_CV$ns, Mod_CV$nt, byrow = TRUE)

              # Save unfitted model
              Path_ModInit <- file.path(
                Path_Init, paste0("InitMod_", CV_Name, "_k", k, ".RData"))
              IASDT.R::SaveAs(
                InObj = Mod_CV,
                OutObj = paste0("InitMod_", CV_Name, "_k", k),
                OutPath = Path_ModInit)

              # initiate sampling and save initial models to
              Mod_CV <- Hmsc::sampleMcmc(
                hM = Mod_CV, samples = ModFull$samples, thin = ModFull$thin,
                transient = ModFull$transient, adaptNf = ModFull$adaptNf,
                initPar = initPar, nChains = NChains, updater = updater,
                verbose = verbose, alignPost = alignPost, engine = "HPC")

              if (ToJSON) {
                Mod_CV <- jsonify::to_json(Mod_CV)
              }
              # Save model input as rds file
              Path_ModInit_rds <- file.path(
                Path_Init, paste0("InitMod_", CV_Name, "_k", k, ".rds"))
              saveRDS(Mod_CV, file = Path_ModInit_rds)

              # Prepare fitting command for each model chain
              CV_Out <- purrr::map_dfr(
                .x = seq_len(NChains),
                .f = function(Chain) {
                  # Path to save the posterior of the combination of CV and
                  # chain
                  Path_Post <- file.path(
                    Path_CV,
                    paste0("Mod_", CV_Name, "_k", k, "_Ch", Chain, "_post.rds"))

                  # Path to save the progress of model fitting
                  Path_ModProg <- stringr::str_replace_all(
                    Path_Post, "post.rds$", "Progress.txt")
                  
                  # Model fitting command
                  Command_HPC <- paste0(
                    "/usr/bin/time -v python3 -m hmsc.run_gibbs_sampler",
                    " --input ", shQuote(Path_ModInit_rds),
                    " --output ", shQuote(Path_Post),
                    " --samples ", ModFull$samples,
                    " --transient ", ModFull$transient,
                    " --thin ", ModFull$thin,
                    " --verbose ", verbose,
                    " --chain ", (Chain - 1),
                    " --fp ", Precision,
                    " >& ", shQuote(Path_ModProg))

                  FullModel <- dplyr::if_else(
                    inherits(Model, "character"), Model, NA_character_)

                  # data to be returned for each combination of CV and Chain
                  tibble::tibble(
                    FullModel = FullModel,
                    CV = k, Chain = Chain,
                    Path_ModInit = Path_ModInit,
                    Path_ModInit_rds = Path_ModInit_rds,
                    Path_Post = Path_Post, Path_ModProg = Path_ModProg,
                    Command_HPC = Command_HPC,
                    val = list(val), valCoords = list(valCoords),
                    LoffVal = list(LoffVal),
                    XVal = list(XVal), XRRRVal = list(XRRRVal)) %>%
                    return()
                })
              return(CV_Out)
            }
          )
          return(CV_DT0)
        })) %>%
    tidyr::unnest(CV_DT0)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model fitting commands -----
  IASDT.R::CatTime("Save model fitting commands")

  CommandFile <- file.path(Path_CV, "Commands2Fit.txt")
  f <- file(CommandFile, open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(CV_DT$Command_HPC, sep = "\n", append = FALSE, file = f)
  close(f)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare SLURM -----

  if (PrepSLURM) {
    IASDT.R::CatTime("Prepare SLURM script")

    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::Mod_SLURM(
      Path_Model = Path_CV, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC, Path_Hmsc = Path_Hmsc,
      Path_SLURM_Out = Path_CV, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save summary data to disk -----
  IASDT.R::CatTime("Save summary data to disk")
  save(CV_DT, file = file.path(Path_CV, "CV_DT.RData"))

  return(invisible(NULL))
}
