## |------------------------------------------------------------------------| #
# mod_CV_fit ----
## |------------------------------------------------------------------------| #

#' Prepare cross-validated Hmsc models for HPC fitting
#'
#' This function prepares cross-validated Hmsc models for fitting using HPC. It
#' handles data preparation, model initialization, and generation of SLURM
#' commands.
#'
#' @param path_model Character. Path to a saved model file (`*.qs2`).
#' @param CV_name Character vector. Column name(s) in the model input data to be
#'   used to cross-validate the models (see [mod_prepare_data] and
#'   [mod_CV_prepare]). The function allows the possibility of using more than
#'   one way of assigning grid cells into cross-validation folders. If multiple
#'   names are provided, separate cross-validation models will be fitted for
#'   each cross-validation type. Currently, there are three cross-validation
#'   strategies: `CV_SAC`, `CV_Dist`, and `CV_Large`. Defaults to `c("CV_Dist",
#'   "CV_Large")`.
#' @param Partitions A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `CV_name` argument. If the `Partitions`
#'   vector is provided, the label used in the output files will be `CV_Custom`.
#' @param init_par a named list of parameter values used for initialization of
#'   MCMC states. See [Hmsc::computePredictedValues] for more information.
#'   Default: `NULL`.
#' @param job_name Character. Name of the submitted job(s) for SLURM. Default:
#'   `CV_Models`.
#' @param updater named `list`. Which conditional updaters should be omitted?
#'   See [Hmsc::computePredictedValues] for more information. Defaults to
#'   `list(Gamma2 = FALSE, GammaEta = FALSE)` to disable the following warnings:
#'   `setting updater$Gamma2=FALSE due to specified phylogeny matrix` and
#'   `setting updater$GammaEta=FALSE: not implemented for spatial methods 'GPP'
#'   and 'NNGP'`.
#' @param align_posterior Logical. Whether the posterior of each chains should
#'   be aligned. See [Hmsc::computePredictedValues] for more information.
#'   Default: `TRUE`.
#' @param ... Additional arguments passed to the [IASDT.R::mod_SLURM] function.
#' @details The function copies part of the [Hmsc::computePredictedValues]
#'   function, which currently does not support performing cross-validation
#'   using Hmsc-HPC.
#' @author Ahmed El-Gabbas
#' @importFrom rlang %||%
#' @inheritParams mod_SLURM
#' @inheritParams mod_inputs
#' @export
#' @name mod_CV_fit

mod_CV_fit <- function(
    path_model = NULL, CV_name = c("CV_Dist", "CV_Large"), Partitions = NULL,
    env_file = ".env", init_par = NULL, job_name = "CV_Models",
    updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    align_posterior = TRUE, to_JSON = FALSE, SLURM_prepare = TRUE,
    memory_per_cpu = NULL, job_runtime = NULL, path_Hmsc = NULL,
    precision = 64, ...) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- Path_ModInit_rds <- CV <- ModName <- partition <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input parameters -----
  IASDT.R::cat_time("Check input parameters")

  NullVarsNames <- c(
    "path_model", "path_Hmsc", "memory_per_cpu", "job_runtime", "env_file")
  NullVars <- which(purrr::map_lgl(.x = NullVarsNames, .f = ~ is.null(get(.x))))

  if (length(NullVars) > 0) {
    NullVarsNames[NullVars]
    stop(
      toString(NullVarsNames[NullVars]), " cannot be NULL", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    .x = AllArgs,
    .f = function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "path_model", "job_name", "env_file",
    "job_runtime", "memory_per_cpu", "path_Hmsc")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = CharArgs, args_type = "character")

  # numeric arguments
  NumericArgs <- c("gpus_per_node", "cpus_per_task", "ntasks")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = NumericArgs, args_type = "numeric")

  if (!(precision %in% c(32, 64))) {
    stop(
      "precision should be either of 32 or 64, not ", precision, call. = FALSE)
  }

  if (!file.exists(path_model)) {
    stop("Model path does not exist.", call. = FALSE)
  }

  rm(AllArgs, NullVarsNames, NullVars, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Loading model -----
  IASDT.R::cat_time("Loading model")

  # Path of the cross-validation folder
  Path_CV <- IASDT.R::path(dirname(dirname(path_model)), "Model_Fitting_CV")

  # Path of the model input data
  Path_ModelData <- list.files(
    path = dirname(dirname(path_model)),
    pattern = "^ModDT_.*subset.RData", full.names = TRUE)
  if (length(Path_ModelData) != 1) {
    stop(
      "There should be exactly one file matches model input data",
      call. = FALSE)
  }

  # Loading full model object
  Model_Full <- IASDT.R::load_as(path_model)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating paths -----
  IASDT.R::cat_time("Creating paths")

  Path_Init <- IASDT.R::path(Path_CV, "Model_Init")
  Path_Fitted <- IASDT.R::path(Path_CV, "Model_Fitted")
  Path_Post <- IASDT.R::path(Path_CV, "Model_Posterior")
  Path_Log <- IASDT.R::path(Path_Post, "JobsLog")
  Dir_Pred <- IASDT.R::path(Path_CV, "Model_Prediction")
  fs::dir_create(c(Path_CV, Path_Init, Path_Fitted, Path_Post, Dir_Pred))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Cross-validation partitions ----
  IASDT.R::cat_time("Cross-validation partitions")

  if (is.null(Partitions)) {
    # if custom partitions is not provided, extract the CV column(s) available
    # in the modelling data

    CV_Data <- IASDT.R::load_as(Path_ModelData)$DT_CV

    if (!(all(CV_name %in% names(CV_Data)))) {
      # if any of the column names does not exist, stop the function
      MissingCV <- CV_name[isFALSE(CV_name %in% names(CV_Data))]
      stop(
        "`Partitions` was not defined (NULL) and column(s) for CV folds ",
        paste(MissingCV, collapse = " + "),
        " can not be found in species data", call. = FALSE
      )
    }

    # Extract CV folds from the CV column(s)
    Partitions <- purrr::map(CV_name, ~ dplyr::pull(CV_Data, .x))
    names(Partitions) <- stringr::str_remove(CV_name, "^CV_")

  } else {
    # If Partitions is provided directly to the function, use "CV_Custom" as CV
    # name
    CV_name <- "CV_Custom"
  }

  # Check the length of CV data equals the number of sampling units in the model
  if (any(purrr::map_int(Partitions, length) != Model_Full$ny)) {
    stop(
      "Partitions parameter must be a vector of the same length of the ",
      "sampling  units of the the full model", call. = FALSE)
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
  IASDT.R::cat_time("Prepare cross-validated initial models")

  # Load coordinates from the full model
  Coords <- IASDT.R::load_as(Path_ModelData)$DT_xy

  CV_DT <- tibble::tibble(
    partition = Partitions, CV_name = names(Partitions)) %>%
    dplyr::mutate(

      # number of cross-validation folders for each CV type
      nfolds = purrr::map_int(.x = partition, .f = ~ length(unique(.x))),

      # prepare data for each cross-validation strategy
      CV_Info = purrr::pmap(
        .l = list(partition, CV_name, nfolds),
        .f = function(partition, CV_name, nfolds) {

          IASDT.R::cat_time(
            paste0("Cross-validation Type: ", CV_name), level = 1)

          CV_DT0 <- purrr::map_dfr(
            .x = seq_len(nfolds),
            .f = function(k) {

              IASDT.R::cat_time(paste0("Fold ", k, "/", nfolds), level = 2)

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
              Path_ModInit <- IASDT.R::path(
                Path_Init, paste0("InitMod_", CV_name, "_k", k, ".RData"))
              IASDT.R::save_as(
                object = Model_CV,
                object_name = paste0("InitMod_", CV_name, "_k", k),
                out_path = Path_ModInit)


              # initiate sampling and save initial models
              Model_CV <- Hmsc::sampleMcmc(
                hM = Model_CV, samples = Model_Full$samples,
                thin = Model_Full$thin, transient = Model_Full$transient,
                adaptNf = Model_Full$adaptNf, initPar = init_par,
                nChains = NChains, updater = updater, verbose = verbose,
                alignPost = align_posterior, engine = "HPC")

              if (to_JSON) {
                Model_CV <- jsonify::to_json(Model_CV)
              }

              # Save model input as rds file
              Path_ModInit_rds <- IASDT.R::path(
                Path_Init, paste0("InitMod_", CV_name, "_k", k, ".rds"))
              saveRDS(Model_CV, file = Path_ModInit_rds)

              Path_ModFitted <- IASDT.R::path(
                Path_Fitted, paste0("Model_", CV_name, "_k", k, ".RData"))

              dfPi <- droplevels(Model_Full$dfPi[val, , drop = FALSE])
              valCoords <- Coords[val, ]
              rownames(dfPi) <- rownames(valCoords) <- dfPi$sample

              return(
                tibble::tibble(
                  Path_ModFull = path_model,
                  CV = k,
                  Path_ModInit = Path_ModInit,
                  Path_ModInit_rds = Path_ModInit_rds,
                  Path_ModFitted = Path_ModFitted,
                  val = list(val),
                  valCoords = list(valCoords),
                  XRRRVal = list(XRRRVal),
                  dfPi = list(dfPi)))
            })
          return(CV_DT0)
        })) %>%
    dplyr::select(-"partition") %>%
    tidyr::unnest("CV_Info")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare model fitting commands -----
  IASDT.R::cat_time("Prepare model fitting commands")

  CV_DT <- CV_DT %>%
    dplyr::mutate(ModName = paste0(CV_name, CV), .after = CV) %>%
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
              Post_File <- IASDT.R::path(
                Path_Post, paste0("Mod_", ModName, "_Ch", Chain, "_post.rds"))

              # Path to save the progress of model fitting
              ModProg_File <- stringr::str_replace_all(
                Post_File, "post.rds$", "Progress.txt")


              # Model fitting command
              Command_HPC <- paste0(
                "/usr/bin/time -v python3 -m hmsc.run_gibbs_sampler",
                " --input ", IASDT.R::normalize_path(Path_ModInit_rds),
                " --output ", IASDT.R::normalize_path(Post_File),
                " --samples ", Model_Full$samples,
                " --transient ", Model_Full$transient,
                " --thin ", Model_Full$thin,
                " --verbose ", verbose,
                " --chain ", (Chain - 1),
                " --fp ", precision,
                " >& ", IASDT.R::normalize_path(ModProg_File))

              # data to be returned for each combination of CV and Chain
              return(
                tibble::tibble(
                  Path_Post = Post_File, Path_ModProg = ModProg_File,
                  Command_HPC = Command_HPC, NSamples = Model_Full$samples,
                  Transient = Model_Full$transient, Thin = Model_Full$thin))
            })

          return(dplyr::summarise_all(CV_Out, list))

        })) %>%
    tidyr::unnest("ModelPrep") %>%
    dplyr::mutate_at(
      .vars = c("NSamples", "Transient", "Thin"),
      .funs = ~ as.integer(unique(unlist(.x))))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model fitting commands -----
  IASDT.R::cat_time("Save model fitting commands")

  purrr::walk(
    .x = CV_name,
    .f = ~{
      Curr_CV_Name <- stringr::str_remove_all(.x, "CV_")
      CommandFile <- IASDT.R::path(Path_CV, paste0("Commands2Fit_", .x, ".txt"))
      f <- file(CommandFile, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)

      CV_DT %>%
        dplyr::filter(CV_name == Curr_CV_Name) %>%
        dplyr::pull(Command_HPC) %>%
        unlist() %>%
        cat(sep = "\n", append = FALSE, file = f)
      close(f)
    })

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare SLURM script -----

  if (SLURM_prepare) {
    IASDT.R::cat_time("Prepare SLURM script")

    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::mod_SLURM(
      model_dir = Path_CV, job_name = job_name, memory_per_cpu = memory_per_cpu,
      job_runtime = job_runtime, env_file = env_file, path_Hmsc = path_Hmsc,
      SLURM_path_out = Path_Log, SLURM_prefix = "CV_Bash_Fit", ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save summary data to disk -----
  IASDT.R::cat_time("Save summary data to disk")
  save(CV_DT, file = IASDT.R::path(Path_CV, "CV_DT.RData"))

  return(invisible(NULL))
}
