## |------------------------------------------------------------------------| #
# Mod_CV ----
## |------------------------------------------------------------------------| #

#' Prepare Cross-Validated Hmsc Models for HPC Fitting
#'
#' This function prepares initial cross-validated Hmsc models in R for fitting
#' by Hmsc-HPC. It includes data preparation, model initialization, and command
#' generation for running models on HPC.
#' @param Model Either a path to a saved model file (character) or an Hmsc model
#'   object. If a path is provided, the model is loaded from the file.
#' @param CVName Character vector specifying the name of the column(s) in the
#'   model input data (see [IASDT.R::Mod_PrepData] and [IASDT.R::GetCV]) to be
#'   used to cross-validate the models. The function allows the possibility of
#'   using more than one way of assigning grid cells into cross-validation
#'   folders. If multiple names are provided, separate cross-validation models
#'   will be fitted for each column. Currently, there are three cross-validation
#'   strategies, created using the [IASDT.R::Mod_PrepData]: `CV_SAC`, `CV_Dist`,
#'   and `CV_Median` (see [IASDT.R::GetCV]).
#' @param partition A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `CVName` argument. If the `partition` vector
#'   is provided, the label used in the output files will be `CV_Custom`.
#' @param Path_CV The directory path where cross-validation models and outputs
#'   will be stored. If `Model` argument is a character vector, it will be
#'   estimated from the model path. If `Model` is an Hmsc model object, it has
#'   to be provided by the user, otherwise the function will give an error.
#' @param EnvFile String specifying the path to read environment variables from,
#'   with a default value of `.env`.
#' @param initPar a named list of parameter values used for initialization of
#'   MCMC states. See [Hmsc::computePredictedValues] for more information. 
#'   Default: `NULL`.
#' @param PrepSLURM Logical indicating whether to prepare SLURM command files.
#'   If `TRUE` (default), the SLURM commands will be saved to disk using the
#'   [IASDT.R::Mod_SLURM] function.
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
#' @param ToJSON Logical indicating whether to convert unfitted models to JSON
#'   before saving to RDS file. Default: `FALSE`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param MemPerCpu String specifying the memory per CPU for the SLURM job. This
#'   value will be assigned to the `#SBATCH --mem-per-cpu=` SLURM argument.
#'   Example: `32G` to request 32 gigabyte. Only effective if `PrepSLURM =
#'   TRUE`.
#' @param Time String specifying the requested time for each job in the SLURM
#'   bash arrays. Example: `01:00:00` to request an hour. Only effective if
#'   `PrepSLURM = TRUE`.
#' @param Path_Hmsc String specifying the path for the Hmsc-HPC. This will be
#'   provided as the `Path_Hmsc` argument of the [IASDT.R::Mod_SLURM] function.
#' @param ... Additional arguments passed to the [IASDT.R::Mod_SLURM] function.
#' @param CheckPyPath Logical, whether to check the existence of the
#'   `Path_Python` directory extracted from reading `EnvFile` content. Default:
#'   `TRUE`.
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
#' @export
#' @name Mod_CV

Mod_CV <- function(
    Model = NULL, CVName = c("CV_SAC", "CV_Dist", "CV_Median"),
    partition = NULL, Path_CV = NULL, EnvFile = ".env", initPar = NULL,
    JobName = "CV_Models", updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    alignPost = TRUE, ToJSON = FALSE, FromHPC = TRUE, PrepSLURM = TRUE,
    MemPerCpu = NULL, Time = NULL, Path_Hmsc = NULL, CheckPyPath = TRUE, ...) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- CV_DT0 <- Path_Python <- NULL

  # # ++++++++++++++++++++++++++++++++++++
  # Check input parameters
  # # ++++++++++++++++++++++++++++++++++++
  if (is.null(Model) || is.null(Path_Hmsc) || is.null(MemPerCpu) ||
      is.null(Time) || is.null(EnvFile)) {
    stop("The following arguments cannot be empty: ",
         "Model, Path_Hmsc, MemPerCpu, Time, EnvFile")
  }

  AllArgs <- ls()
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

  rm(AllArgs)

  # # ++++++++++++++++++++++++++++++++++++
  # Loading environment variables
  # # ++++++++++++++++++++++++++++++++++++

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Python", "DP_R_Path_Python", CheckPyPath, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Python", "DP_R_Path_Python", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ++++++++++++++++++++++++++++++++++++
  # Loading model
  # # ++++++++++++++++++++++++++++++++++++
  if (inherits(Model, "character")) {
    # Ensure that model path does not contain the scratch path
    Path_CV <- file.path(dirname(dirname(Model)), "Model_Fitting_HPC_CV")
    ModFull <- IASDT.R::LoadAs(Model)
  } else {
    if (is.null(Path_CV)) {
      stop("Path_CV cannot be empty if the Model is provided as Hmsc object")
    }
    ModFull <- Model
    rm(Model)
  }

  # # ++++++++++++++++++++++++++++++++++++
  # Creating paths
  # # ++++++++++++++++++++++++++++++++++++
  Path_Init <- file.path(Path_CV, "InitMod")
  fs::dir_create(c(Path_CV, Path_Init))

  # # ++++++++++++++++++++++++++++++++++++
  # Cross-validation partitions
  # # ++++++++++++++++++++++++++++++++++++
  if (is.null(partition)) {
    # if custom partition is not provided, extract the CV column(s) available in
    # the modelling data
    if (all(CVName %in% names(ModFull$XData))) {
      # Extract CV folds from the CV column(s)
      partition <- purrr::map(CVName, ~dplyr::pull(ModFull$XData, .x))
      names(partition) <- CVName
    } else {
      # if any of the column names does not exist, stop the function
      stop("partition was not defined and column(s) for CV folds ",
           "can not be found in species data")
    }
  } else {
    # If partition is provided directly to the function, use "CV_Custom" as CV
    # name
    CVName <- "CV_Custom"
  }

  # Check the length of CV data equals the number of sampling units in the model
  if (any(sapply(partition, length) != ModFull$ny)) {
    stop("partition parameter must be a vector of length ny")
  }

  # # ++++++++++++++++++++++++++++++++++++
  # Verbose
  # # ++++++++++++++++++++++++++++++++++++

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

  # # ++++++++++++++++++++++++++++++++++++
  # Number of chains
  # # ++++++++++++++++++++++++++++++++++++
  nChains <- length(ModFull$postList)

  # # ++++++++++++++++++++++++++++++++++++
  # Prepare cross-validated initial models / model fitting commands / SLURM jobs
  # # ++++++++++++++++++++++++++++++++++++

  CV_DT <- tibble::tibble(partition = partition, CVName = names(partition)) %>%
    dplyr::mutate(
      nfolds = purrr::map_int(.x = partition, .f = ~length(unique(.x))),
      # Loop over all cross-validation partitions used
      CV_DT0 = purrr::pmap(
        .l = list(partition, CVName, nfolds),
        .f = function(CV_Fold, CV_Name, CV_Folds) {

          # prepare data for each cross-validation strategy

          CV_DT0 <- purrr::map_dfr(
            .x = seq_len(CV_Folds),
            .f = function(k) {

              # The following is copied from Hmsc::computePredictedValues()
              train <- (CV_Fold != k)
              val <- (CV_Fold == k)
              dfPi <- as.data.frame(matrix(NA, sum(train), ModFull$nr))
              colnames(dfPi) <- ModFull$rLNames

              for (r in seq_len(ModFull$nr)) {
                dfPi[, r] <- factor(ModFull$dfPi[train, r])
              }

              switch(
                class(ModFull$X)[1L],
                matrix = {
                  XTrain <- ModFull$X[train, , drop = FALSE]
                  XVal <- ModFull$X[val, , drop = FALSE]
                },
                list = {
                  XTrain <- purrr::map(
                    .x = ModFull$X, .f = ~.x[train, , drop = FALSE])
                  XVal <- purrr::map(
                    .x = ModFull$X, .f = ~.x[val, , drop = FALSE])
                }
              )

              if (ModFull$ncRRR > 0) {
                XRRRTrain <- ModFull$XRRR[train, , drop = FALSE]
                XRRRVal <- ModFull$XRRR[val, , drop = FALSE]
              } else {
                XRRRTrain <- NULL
                XRRRVal <- NULL
              }

              # Initial models
              ModCV <- Hmsc::Hmsc(
                Y = ModFull$Y[train, , drop = FALSE],
                Loff = ModFull$Loff[train, , drop = FALSE],
                X = XTrain, XRRR = XRRRTrain, ncRRR = ModFull$ncRRR,
                XSelect = ModFull$XSelect, distr = ModFull$distr,
                studyDesign = dfPi, Tr = ModFull$Tr, C = ModFull$C,
                ranLevels = ModFull$rL)

              ModCV <- Hmsc::setPriors(
                ModCV, V0 = ModFull$V0, f0 = ModFull$f0,
                mGamma = ModFull$mGamma, UGamma = ModFull$UGamma,
                aSigma = ModFull$aSigma, bSigma = ModFull$bSigma,
                rhopw = ModFull$rhowp)

              ModCV$YScalePar <- ModFull$YScalePar
              ModCV$YScaled <- (ModCV$Y - matrix(ModCV$YScalePar[1, ], ModCV$ny, ModCV$ns, byrow = TRUE)) / matrix(ModCV$YScalePar[2, ], ModCV$ny, ModCV$ns, byrow = TRUE)
              ModCV$XInterceptInd <- ModFull$XInterceptInd
              ModCV$XScalePar <- ModFull$XScalePar

              switch(
                class(ModFull$X)[1L],
                matrix = {
                  ModCV$XScaled <- (ModCV$X - matrix(ModCV$XScalePar[1, ], ModCV$ny, ModCV$ncNRRR, byrow = TRUE)) / matrix(ModCV$XScalePar[2, ], ModCV$ny, ModCV$ncNRRR, byrow = TRUE)
                },
                list = {
                  ModCV$XScaled <- list()
                  for (zz in seq_len(length(ModCV$X))) {
                    ModCV$XScaled[[zz]] <- (ModCV$X[[zz]] - matrix(ModCV$XScalePar[1, ], ModCV$ny, ModCV$ncNRRR, byrow = TRUE)) / matrix(ModCV$XScalePar[2, ], ModCV$ny, ModCV$ncNRRR, byrow = TRUE)
                  }
                }
              )

              if (ModCV$ncRRR > 0) {
                ModCV$XRRRScalePar <- ModFull$XRRRScalePar
                ModCV$XRRRScaled <- (ModCV$XRRR - matrix(ModCV$XRRRScalePar[1, ], ModCV$ny, ModCV$ncORRR, byrow = TRUE)) / matrix(ModCV$XRRRScalePar[2, ], ModCV$ny, ModCV$ncORRR, byrow = TRUE)
              }

              ModCV$TrInterceptInd <- ModFull$TrInterceptInd
              ModCV$TrScalePar <- ModFull$TrScalePar
              ModCV$TrScaled <- (ModCV$Tr - matrix(ModCV$TrScalePar[1, ], ModCV$ns, ModCV$nt, byrow = TRUE)) / matrix(ModCV$TrScalePar[2, ], ModCV$ns, ModCV$nt, byrow = TRUE)

              # initiate sampling and save initial models to
              ModCV <- Hmsc::sampleMcmc(
                hM = ModCV, samples = ModFull$samples, thin = ModFull$thin,
                transient = ModFull$transient, adaptNf = ModFull$adaptNf,
                initPar = initPar, nChains = nChains, updater = updater,
                verbose = verbose, alignPost = alignPost, engine = "HPC")

              if (ToJSON) {
                ModCV <- jsonify::to_json(ModCV)
              }

              # Save unfitted models to disk
              Path_ModInit <- file.path(
                Path_Init, paste0("Mod_", CV_Name, "_", k, ".rds"))
              saveRDS(ModCV, file = Path_ModInit)

              # Prepare fitting command for each model chain
              CV_Out <- purrr::map_dfr(
                .x = seq_len(nChains),
                .f = function(Chain) {

                  # Path to save the posterior of the combination of CV and chain
                  Path_Post <- file.path(
                    Path_CV,
                    paste0("Mod_", CV_Name, "_", k, "_Ch", Chain, "_post.rds"))

                  # Path to save the progress of model fitting
                  Path_ModProg <- stringr::str_replace_all(
                    Path_Post, "post.rds$", "Progress.txt")

                  # Model fitting command
                  Command_HPC <- paste0(
                    "export TF_CPP_MIN_LOG_LEVEL=3; /usr/bin/time -v ",
                    Path_Python,
                    " -m hmsc.run_gibbs_sampler",
                    " --input ", shQuote(Path_ModInit),
                    " --output ", shQuote(Path_Post),
                    " --samples ", ModFull$samples,
                    " --transient ", ModFull$transient,
                    " --thin ", ModFull$thin,
                    " --verbose ", verbose,
                    " --chain ", (Chain - 1),
                    " >& ", shQuote(Path_ModProg))

                  # data to be returned for each combination of CV and Chain
                  return(
                    tibble::tibble(
                      CV = k, Chain = Chain, Path_ModInit = Path_ModInit,
                      Path_Post = Path_Post, Path_ModProg = Path_ModProg,
                      Command_HPC = Command_HPC))
                })

              return(CV_Out)
            }
          )
          return(CV_DT0)
        }
      )) %>%
    dplyr::select(-partition) %>%
    tidyr::unnest(CV_DT0)

  # Save model fitting commands to "Commands2Fit.txt" file
  CommandFile <- file.path(Path_CV, "Commands2Fit.txt")
  f <- file(CommandFile, open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(CV_DT$Command_HPC, sep = "\n", append = FALSE, file = f)
  close(f)

  if (PrepSLURM) {
    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::Mod_SLURM(
      Path_Model = Path_CV, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC, Path_Hmsc = Path_Hmsc,
      Path_SLURM_Out = Path_CV, ...)
  }

  # Save summary data to disk
  save(CV_DT, file = file.path(Path_CV, "CV_DT.RData"))

  return(invisible(NULL))
}
