## |------------------------------------------------------------------------| #
# Mod_Postprocess ----
## |------------------------------------------------------------------------| #

#' Model Pipeline for Hmsc Analysis
#'
#' This function sets up and runs an analysis pipeline for Hmsc models. It
#' includes steps for environment setup, loading packages, managing SLURM
#' refits, merging MCMC chains, convergence diagnostics, model summaries,
#' spatial predictions, response curve generation, and variance partitioning.
#' @param GPP_Dist Integer specifying the distance in *kilometers* between knots
#'   for GPP models.
#' @param Tree Character string specifying if phylogenetic tree was used in the
#'   model. Valid values are "Tree" or "NoTree". Default is "Tree".
#' @param Samples Integer specifying the value for the number of MCMC samples in
#'   the selected model. Defaults to 1000.
#' @param Thin Integer specifying the value for thinning in the selected model.
#'
#' @name Mod_Postprocess
#' @inheritParams Predict_Maps
#' @inheritParams Mod_CV_Fit
#' @inheritParams Merge_Chains
#' @inheritParams Mod_Prep4HPC
#' @inheritParams RespCurv_PrepData
#' @inheritParams Coda_to_tibble
#' @author Ahmed El-Gabbas
#' @export

Mod_Postprocess <- function(
    ModelDir = NULL, Hab_Abb = NULL, NCores = 8L, FromHPC = TRUE,
    EnvFile = ".env", Path_Hmsc = NULL, MemPerCpu = NULL, Time = NULL,
    FromJSON = FALSE, GPP_Dist = NULL, Tree = "Tree", Samples = 1000L,
    Thin = NULL, N_Grid = 50L, NOmega = 1000L, UseTF = TRUE, TF_Environ = NULL,
    TF_use_single = FALSE, LF_NCores = NCores, LF_Check = FALSE,
    LF_Temp_Cleanup = TRUE, Temp_Cleanup = TRUE,
    CC_Models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_Scenario = c("ssp126", "ssp370", "ssp585"),
    Pred_Clamp = TRUE, Fix_Efforts = "mean", Pred_NewSites = TRUE,
    CVName = c("CV_Dist", "CV_Large")) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  Ch1 <- function(Text) {
    IASDT.R::InfoChunk(
      paste0("\t", Text), Rep = 2, Char = "=", CharReps = 60, Red = TRUE,
      Bold = TRUE, Time = FALSE)
  }

  Ch2 <- function(Text) {
    IASDT.R::InfoChunk(
      paste0("\t", Text), Rep = 1, Char = "-", CharReps = 60, Red = TRUE,
      Bold = TRUE, Time = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check input arguments ----

  Hab_Abb <- as.character(Hab_Abb)

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Hab_Abb", "EnvFile", "ModelDir", "Tree", "Path_Hmsc"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("UseTF", "FromHPC", "Pred_Clamp", "Pred_NewSites", "FromJSON"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "NOmega", "N_Grid", "GPP_Dist", "Samples", "Thin"))
  rm(AllArgs, envir = environment())


  ValidHabAbbs <- c(as.character(0:3), "4a", "4b", "10", "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0(
        "Invalid Habitat abbreviation. Valid values are:\n >> ",
        paste0(ValidHabAbbs, collapse = ", ")),
      call. = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(
      paste0("Path for environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (!dir.exists(ModelDir)) {
    stop(paste0("Model directory: ", ModelDir, " was not found"), call. = FALSE)
  }

  if (Pred_Clamp && is.null(Fix_Efforts)) {
    stop(
      "`Fix_Efforts` can not be NULL when Clamping is implemented",
      call. = FALSE)
  }


  if (!all(
    CC_Models %in% c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"))) {
    stop(
      paste0(
        "Invalid climate models. Valid values are:\n >> ",
        paste0(
          c(
            "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
            "MRI-ESM2-0", "UKESM1-0-LL"), collapse = ", ")),
      call. = FALSE)
  }

  if (!all(CC_Scenario %in% c("ssp126", "ssp370", "ssp585"))) {
    stop(
      paste0(
        "Invalid climate scenarios. Valid values are:\n >> ",
        paste0(c("ssp126", "ssp370", "ssp585"), collapse = ", ")),
      call. = FALSE)
  }


  if (!(Tree %in% c("Tree", "NoTree"))) {
    stop(
      paste0(
        "Invalid value for Tree argument. Valid values ",
        "are: 'Tree' or 'NoTree'"),
      call. = FALSE)
  }

  if (!all(CVName %in% c("CV_Dist", "CV_Large"))) {
    stop(
      paste0(
        "Invalid value for CVName argument. Valid values ",
        "are: 'CV_Dist' or 'CV_Large'"),
      call. = FALSE)
  }


  LoadedPackages <- paste0(
    sort(IASDT.R::LoadedPackages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", IASDT.R::CurrOS(),
      "\n  >>> Model root: ", ModelDir,
      "\n  >>> NCores: ", NCores,
      "\n  >>> FromHPC: ", FromHPC,
      "\n  >>> EnvFile: ", EnvFile,
      "\n  >>> Path_Hmsc: ", Path_Hmsc,
      "\n  >>> SLURM MemPerCpu: ", MemPerCpu,
      "\n  >>> SLURM Time: ", Time,
      "\n  >>> NOmega: ", NOmega,
      "\n  >>> Hab_Abb: ", Hab_Abb,
      "\n  >>> UseTF: ", UseTF,
      "\n  >>> Python environment: ", TF_Environ,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"),
      "\n  >>> Loaded packages: \n", LoadedPackages))


  Temp_Dir <- file.path(ModelDir, "TEMP_Pred")

  ModelData <- list.files(
    path = ModelDir, full.names = TRUE,
    pattern = paste0("^ModDT_", Hab_Abb, "_.+_subset.RData"))

  if (length(ModelData) != 1) {
    stop("Model data was not found", call. = FALSE)
  }


  # # ..................................................................... ###
  # # ..................................................................... ###

  # Check unsuccessful models -----
  Ch1("Check unsuccessful models")

  IASDT.R::Mod_SLURM_Refit(
    ModelDir = ModelDir,
    JobName = stringr::str_remove(basename(ModelDir), "Mod_"),
    MemPerCpu = MemPerCpu, Time = Time, EnvFile = EnvFile, FromHPC = FromHPC,
    Path_Hmsc = Path_Hmsc)

  invisible(gc())

  # ****************************************************************

  Ch2("Merge chains and saving RData files")

  IASDT.R::Merge_Chains(
    ModelDir = ModelDir, NCores = NCores,
    FromHPC = FromHPC, FromJSON = FromJSON)

  invisible(gc())

  # ****************************************************************

  Ch2("Convergence of all model variants")

  IASDT.R::Convergence_Plot_All(
    ModelDir = ModelDir, maxOmega = NOmega, NCores = NCores,
    FromHPC = FromHPC, MarginType = "histogram")

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Path of selected model -----
  Ch1("Path of selected model")

  Path_Model <- file.path(
    ModelDir, "Model_Fitted",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Model.qs2"))

  Path_Coda <- file.path(
    ModelDir, "Model_Coda",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Coda.qs2"))

  cat(paste0("Path_Model:\n\t", Path_Model, "\nPath_Coda:\n\t", Path_Coda))

  if (!all(file.exists(Path_Model, Path_Coda))) {
    stop("Selected model files not found", call. = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Convergence ----
  Ch1("Convergence")

  ## Gelman_Plot -----
  Ch2("Gelman_Plot")

  IASDT.R::PlotGelman(
    Path_Coda = Path_Coda, Alpha = TRUE, Beta = TRUE, Omega = TRUE, Rho = TRUE,
    NOmega = NOmega, FromHPC = FromHPC, EnvFile = EnvFile)

  invisible(gc())

  # ****************************************************************

  ## Convergence plots ----
  Ch2("Convergence plots")

  IASDT.R::Convergence_Plot(
    Path_Coda = Path_Coda, Path_Model = Path_Model, EnvFile = EnvFile,
    FromHPC = FromHPC, NOmega = NOmega, NCores = NCores, NRC = c(2, 2),
    Beta_NRC = c(3, 3), MarginType = "histogram")

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Response curves -----
  Ch1("Response curves")

  ## Prepare data ------
  Ch2("Prepare data")

  IASDT.R::RespCurv_PrepData(
    Path_Model = Path_Model, N_Grid = N_Grid, NCores = NCores, UseTF = UseTF,
    TF_Environ = TF_Environ, TF_use_single = TF_use_single,
    LF_NCores = LF_NCores, LF_Temp_Cleanup = LF_Temp_Cleanup,
    LF_Check = LF_Check, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    Verbose = TRUE)

  invisible(gc())

  # ****************************************************************

  ## Plotting - species richness ------
  Ch2("Plotting - species richness")

  IASDT.R::RespCurv_PlotSR(ModelDir = ModelDir, Verbose = TRUE, NCores = NCores)
  invisible(gc())

  ## Plotting - species -----
  Ch2("Plotting - species")

  IASDT.R::RespCurv_PlotSp(
    ModelDir = ModelDir, NCores = NCores, EnvFile = EnvFile, FromHPC = FromHPC)

  invisible(gc())

  # ****************************************************************

  ## Plotting - all species together -------
  Ch2("Plotting - all species together")

  IASDT.R::RespCurv_PlotSpAll(ModelDir = ModelDir, NCores = NCores)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Model summary ------
  Ch1("Model summary")

  IASDT.R::Mod_Summary(
    Path_Coda = Path_Coda, EnvFile = EnvFile, FromHPC = FromHPC)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Plotting model parameters -----
  Ch1("Plotting model parameters")


  ## Omega -----
  Ch2("Plotting Omega parameter")

  IASDT.R::PlotOmegaGG(
    Path_Model = Path_Model, supportLevel = 0.95,
    PlotWidth = 22, PlotHeight = 20)

  invisible(gc())

  # ****************************************************************

  ## Beta ------
  Ch2("Plotting Beta parameters")

  IASDT.R::PlotBetaGG(
    Path_Model = Path_Model, supportLevel = 0.95,
    PlotWidth = 26, PlotHeight = 20)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Spatial predictions / evaluation ---------
  Ch1("Spatial predictions / evaluation")

  Model_Predictions <- IASDT.R::Predict_Maps(
    Path_Model = Path_Model, Hab_Abb = Hab_Abb, EnvFile = EnvFile,
    FromHPC = FromHPC, NCores = NCores, Pred_Clamp = Pred_Clamp,
    Fix_Efforts = Fix_Efforts, Pred_NewSites = Pred_NewSites, UseTF = UseTF,
    TF_Environ = TF_Environ, CC_Models = CC_Models, CC_Scenario = CC_Scenario,
    Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    TF_use_single = TF_use_single, LF_NCores = LF_NCores, LF_Check = LF_Check,
    LF_Temp_Cleanup = LF_Temp_Cleanup)

  rm(Model_Predictions)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Variance partitioning ------
  Ch1("Variance partitioning")

  IASDT.R::VarPar_Plot(
    Path_Model = Path_Model, EnvFile = EnvFile, FromHPC = FromHPC,
    UseTF = UseTF, TF_Environ = TF_Environ, NCores = 5,
    Fig_width = 30, Fig_height = 15)

  # # ..................................................................... ###
  # # ..................................................................... ###


  # Cross-validation -------

  Ch1("Prepare input data for cross-validation")
  IASDT.R::Mod_CV_Fit(
    Model = Path_Model, ModelData = ModelData, CVName = CVName,
    EnvFile = EnvFile, JobName = paste0("CV_", Hab_Abb), FromHPC = FromHPC,
    MemPerCpu = MemPerCpu, Time = Time, Path_Hmsc = Path_Hmsc)

  # # ..................................................................... ###
  # # ..................................................................... ###

  CatTime()
  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPostprocessing took ")

  # # ..................................................................... ###
  # # ..................................................................... ###

  return(invisible(NULL))

}


# # ========================================================================== #
# # ========================================================================== #


## |------------------------------------------------------------------------| #
# Mod_Postprocess_1_CPU ----
## |------------------------------------------------------------------------| #

#' Model Pipeline for Hmsc Analysis (using CPU)
#'
#' This function sets up and runs an analysis pipeline for Hmsc models. It
#' includes steps for environment setup, loading packages, managing SLURM
#' refits, merging MCMC chains, convergence diagnostics, model summaries, and
#' prepare input commands to use TensorFlow on GPU for predicting latent factor
#' for new spatial locations and response curves.
#' @param GPP_Dist Integer specifying the distance in *kilometers* between knots
#'   for GPP models.
#' @param Tree Character string specifying if phylogenetic tree was used in the
#'   model. Valid values are "Tree" or "NoTree". Default is "Tree".
#' @param Samples Integer specifying the value for the number of MCMC samples in
#'   the selected model. Defaults to 1000.
#' @param Thin Integer specifying the value for thinning in the selected model.
#' @param NCores_VP Integer specifying the number of cores to use for variance
#'   partitioning. Defaults to 6.
#' @name Mod_Postprocess_1_CPU
#' @inheritParams Predict_Maps
#' @inheritParams Mod_CV_Fit
#' @inheritParams Merge_Chains
#' @inheritParams Mod_Prep4HPC
#' @inheritParams RespCurv_PrepData
#' @inheritParams Coda_to_tibble
#' @author Ahmed El-Gabbas
#' @export

Mod_Postprocess_1_CPU <- function(
    ModelDir = NULL, Hab_Abb = NULL, NCores = 8L, FromHPC = TRUE,
    EnvFile = ".env", Path_Hmsc = NULL, MemPerCpu = NULL, Time = NULL,
    FromJSON = FALSE, GPP_Dist = NULL, Tree = "Tree", Samples = 1000L,
    Thin = NULL, NOmega = 1000L, CVName = c("CV_Dist", "CV_Large"),
    N_Grid = 50L, UseTF = TRUE, TF_use_single = FALSE, LF_NCores = NCores,
    LF_Temp_Cleanup = TRUE, LF_Check = FALSE, Temp_Cleanup = TRUE,
    TF_Environ = NULL, Pred_Clamp = TRUE, Fix_Efforts = "mean",
    Pred_NewSites = TRUE, NCores_VP = 6) {

  .StartTime <- lubridate::now(tzone = "CET")

  # ****************************************************************

  Ch1 <- function(Text) {
    IASDT.R::InfoChunk(
      paste0("\t", Text), Rep = 2, Char = "=", CharReps = 60, Red = TRUE,
      Bold = TRUE, Time = FALSE)
  }

  # ****************************************************************

  # Check input arguments ----

  Hab_Abb <- as.character(Hab_Abb)

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Hab_Abb", "EnvFile", "ModelDir", "Tree", "Path_Hmsc"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c(
      "FromHPC", "FromJSON", "Pred_NewSites", "Pred_Clamp", "Temp_Cleanup",
      "UseTF", "LF_Check", "LF_Temp_Cleanup", "TF_use_single"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c(
      "NCores", "NOmega", "GPP_Dist", "Samples", "Thin", "NCores_VP",
      "LF_NCores", "N_Grid"))
  rm(AllArgs, envir = environment())


  ValidHabAbbs <- c(as.character(0:3), "4a", "4b", "10", "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      paste0(
        "Invalid Habitat abbreviation. Valid values are:\n >> ",
        paste0(ValidHabAbbs, collapse = ", ")),
      call. = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(
      paste0("Path for environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (!dir.exists(ModelDir)) {
    stop(paste0("Model directory: ", ModelDir, " was not found"), call. = FALSE)
  }

  if (!(Tree %in% c("Tree", "NoTree"))) {
    stop(
      paste0(
        "Invalid value for Tree argument. Valid values ",
        "are: 'Tree' or 'NoTree'"),
      call. = FALSE)
  }

  if (!all(CVName %in% c("CV_Dist", "CV_Large"))) {
    stop(
      paste0(
        "Invalid value for CVName argument. Valid values ",
        "are: 'CV_Dist' or 'CV_Large'"),
      call. = FALSE)
  }

  # ****************************************************************

  LoadedPackages <- paste0(
    sort(IASDT.R::LoadedPackages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", IASDT.R::CurrOS(),
      "\n  >>> Model root: ", ModelDir,
      "\n  >>> NCores: ", NCores,
      "\n  >>> FromHPC: ", FromHPC,
      "\n  >>> EnvFile: ", EnvFile,
      "\n  >>> Path_Hmsc: ", Path_Hmsc,
      "\n  >>> SLURM MemPerCpu: ", MemPerCpu,
      "\n  >>> SLURM Time: ", Time,
      "\n  >>> NOmega: ", NOmega,
      "\n  >>> Hab_Abb: ", Hab_Abb,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"),
      "\n  >>> Loaded packages: \n", LoadedPackages))

  # ****************************************************************

  ModelData <- list.files(
    path = ModelDir, full.names = TRUE,
    pattern = paste0("^ModDT_", Hab_Abb, "_.+_subset.RData"))

  if (length(ModelData) != 1) {
    stop("Model data was not found", call. = FALSE)
  }

  # ****************************************************************

  # Check unsuccessful models -----
  Ch1("Check unsuccessful models")

  IASDT.R::Mod_SLURM_Refit(
    ModelDir = ModelDir,
    JobName = stringr::str_remove(basename(ModelDir), "Mod_"),
    MemPerCpu = MemPerCpu, Time = Time, EnvFile = EnvFile, FromHPC = FromHPC,
    Path_Hmsc = Path_Hmsc)

  invisible(gc())

  # ****************************************************************

  Ch1("Merge chains and saving RData files")

  IASDT.R::Merge_Chains(
    ModelDir = ModelDir, NCores = NCores,
    FromHPC = FromHPC, FromJSON = FromJSON)

  invisible(gc())

  # ****************************************************************

  Ch1("Convergence of all model variants")

  IASDT.R::Convergence_Plot_All(
    ModelDir = ModelDir, maxOmega = NOmega, NCores = NCores,
    FromHPC = FromHPC, MarginType = "histogram")

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Path of selected model -----
  Ch1("Path of selected model")

  Path_Model <- file.path(
    ModelDir, "Model_Fitted",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Model.qs2"))

  Path_Coda <- file.path(
    ModelDir, "Model_Coda",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Coda.qs2"))

  cat(
    paste0("Path_Model:\n\t", Path_Model, "\nPath_Coda:\n\t", Path_Coda, "\n"))

  if (!all(file.exists(Path_Model, Path_Coda))) {
    stop("Selected model files not found", call. = FALSE)
  }

  Temp_Dir <- file.path(ModelDir, "TEMP_Pred")

  # ****************************************************************

  # Gelman_Plot -----
  Ch1("Gelman_Plot")

  IASDT.R::PlotGelman(
    Path_Coda = Path_Coda, Alpha = TRUE, Beta = TRUE, Omega = TRUE, Rho = TRUE,
    NOmega = NOmega, FromHPC = FromHPC, EnvFile = EnvFile)

  invisible(gc())

  # ****************************************************************

  # Convergence plots ----
  Ch1("Convergence plots")

  IASDT.R::Convergence_Plot(
    Path_Coda = Path_Coda, Path_Model = Path_Model, EnvFile = EnvFile,
    FromHPC = FromHPC, NOmega = NOmega, NCores = NCores, NRC = c(2, 2),
    Beta_NRC = c(3, 3), MarginType = "histogram")

  invisible(gc())

  # ****************************************************************

  # Model summary ------
  Ch1("Model summary")

  IASDT.R::Mod_Summary(
    Path_Coda = Path_Coda, EnvFile = EnvFile, FromHPC = FromHPC)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Omega -----
  Ch1("Plotting model parameters - Omega")

  IASDT.R::PlotOmegaGG(
    Path_Model = Path_Model, supportLevel = 0.95,
    PlotWidth = 22, PlotHeight = 20)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Beta -----
  Ch1("Plotting model parameters - Beta")

  IASDT.R::PlotBetaGG(
    Path_Model = Path_Model, supportLevel = 0.95,
    PlotWidth = 26, PlotHeight = 20)

  invisible(gc())

  # ****************************************************************

  # Prepare input data for cross-validation -------

  Ch1("Prepare input data for cross-validation")
  IASDT.R::Mod_CV_Fit(
    Model = Path_Model, ModelData = ModelData, CVName = CVName,
    EnvFile = EnvFile, JobName = paste0("CV_", Hab_Abb), FromHPC = FromHPC,
    MemPerCpu = MemPerCpu, Time = Time, Path_Hmsc = Path_Hmsc)

  invisible(gc())

  # ****************************************************************

  # Prepare scripts for predicting latent factors of the response curves -------
  Ch1("Prepare scripts for predicting latent factors of the response curves")

  IASDT.R::RespCurv_PrepData(
    Path_Model = Path_Model, N_Grid = N_Grid, NCores = NCores, UseTF = UseTF,
    TF_Environ = TF_Environ, TF_use_single = TF_use_single,
    LF_NCores = LF_NCores, LF_Temp_Cleanup = LF_Temp_Cleanup,
    LF_Check = LF_Check, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    Verbose = TRUE, LF_Commands_Only = TRUE)

  invisible(gc())

  # ****************************************************************

  # Prepare scripts for predicting latent factors for new sampling units -------
  Ch1("Prepare scripts for predicting latent factors for new sampling units")

  IASDT.R::Predict_Maps(
    Path_Model = Path_Model, Hab_Abb = Hab_Abb, EnvFile = EnvFile,
    FromHPC = FromHPC, NCores = NCores, Pred_Clamp = Pred_Clamp,
    Fix_Efforts = Fix_Efforts, Pred_NewSites = Pred_NewSites, UseTF = UseTF,
    TF_Environ = TF_Environ, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    TF_use_single = TF_use_single, LF_NCores = LF_NCores, LF_Check = LF_Check,
    LF_Temp_Cleanup = LF_Temp_Cleanup, LF_Only = TRUE, LF_Commands_Only = TRUE)

  # ****************************************************************

  # Prepare scripts for computing variance partitioning -------
  Ch1("Prepare scripts for computing variance partitioning")

  IASDT.R::VarPar_Compute(
    Path_Model = Path_Model,
    NCores = NCores_VP, UseTF = UseTF, TF_Environ = TF_Environ,
    TF_use_single = TF_use_single, Temp_Cleanup = Temp_Cleanup,
    Chunk_size = 50L, Verbose = TRUE, VarParFile = "VarPar",
    VP_Commands_Only = TRUE)

  # ****************************************************************

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPostprocessing took ")

  return(invisible(NULL))
}


# # ========================================================================== #
# # ========================================================================== #


## |------------------------------------------------------------------------| #
# Mod_Prep_TF ----
## |------------------------------------------------------------------------| #

#' Prepare Batch Scripts for LF Prediction on GPU
#'
#' This function automates the preparation of batch scripts for predicting the
#' latent factors (LF) of response curves and at new sites using GPU resources.
#' It reads input files that match a specified pattern, merges their contents,
#' sorts the commands, and distributes them into a user-defined number of output
#' files. Each batch script is designed to be compatible with an HPC
#' environment, such as LUMI, with TensorFlow setup included. The function
#' limits the number of output files to a specified maximum, defaults to 210 for
#' compatibility with LUMI's job limits.
#'
#' @param Path Character. Directory containing input files with commands.
#' @param NumFiles Integer. Number of output batch files to create. Must be less
#'   than or equal to the maximum job limit of the HPC environment.
#' @param Pattern Character. Regex pattern to match input files.
#' @param WD Character. Optional. Working directory to be set in batch scripts.
#'   If If `NULL`, the working directory will not be changed.
#' @param WD Character. Working directory for batch files. If `NULL`, defaults
#'   to the current directory.
#' @param Path_Out Character. Directory to save output files. Default is
#'   `TF_BatchFiles`.
#' @note This function is designed specifically for the LUMI HPC environment. It
#'   assumes the `tensorflow` module is available and pre-configured with all
#'   necessary Python packages. On other HPC systems, users may need to modify
#'   the function to load a Python virtual environment or install the required
#'   dependencies for TensorFlow and related packages.
#' @return None. Writes batch files to `Path_Out`.
#' @author Ahmed El-Gabbas
#' @name Mod_Prep_TF
#' @export

Mod_Prep_TF <- function(
    Path = "datasets/processed/model_fitting", NumFiles = 210,
    Pattern = "LF_.+_Commands_.+.txt|LF_.+_Commands_.+txt|VP_.+Command.txt",
    WD = NULL, Path_Out = "TF_BatchFiles") {


  # ****************************************************************

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", 
    Args = c("Path", "Pattern", "Path_Out"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = "NumFiles")
  rm(AllArgs, envir = environment())


  if (!dir.exists(Path)) {
    stop("`Path` must be a valid directory.", call. = FALSE)
  }

  if (Pattern == "") {
    stop("`Pattern` must be a non-empty string.", call. = FALSE)
  }

  if (NumFiles <= 0) {
    stop("`NumFiles` must be a positive integer.", call. = FALSE)
  }

  # ****************************************************************

  Path_Out <- file.path(Path, Path_Out)
  fs::dir_create(Path_Out)

  # ****************************************************************

  # Prepare batch files ----

  # Basic commands for TensorFlow setup
  BasicCommands <- c(
    "#!/bin/bash\n",
    "# Load TensorFlow module and configure environment",
    "ml use /appl/local/csc/modulefiles",
    "ml tensorflow",
    "export TF_CPP_MIN_LOG_LEVEL=3",
    "export TF_ENABLE_ONEDNN_OPTS=0\n",
    "# Verify GPU availability",
    paste0(
      'python3 -c "import tensorflow as tf; ',
      'print(\\\"Num GPUs Available:\\\", ',
      'len(tf.config.list_physical_devices(\\\"GPU\\\")))"'),
    "")

  # Change working directory if specified
  if (!is.null(WD)) {
    WD <- normalizePath(WD, winslash = "/", mustWork = TRUE)
    BasicCommands <- c(
      BasicCommands, "# Change to working directory", paste0("cd ", WD), "")
  }

  # Find list of files matching the pattern
  LF_InFiles <- list.files(
    path = Path, pattern = Pattern, recursive = TRUE, full.names = TRUE)

  if (length(LF_InFiles) == 0) {
    stop(
      paste0("No files found matching the pattern `", Pattern, "` in ", Path),
      call. = FALSE)
  }

  cat(
    paste0(
      "Found ", length(LF_InFiles), " files matching the pattern `",
      Pattern, "`.\n"),
    paste0(" >>>  ", LF_InFiles, sep = "\n"))

  
  # Read and merge commands from input files
  LF_commands <- purrr::map(LF_InFiles, readr::read_lines) %>%
    unlist() %>%
    gtools::mixedsort()
  
  paste0(
    "Total number of commands to be executed: ",
    length(LF_commands), "\n") %>%
    cat()

  if (length(LF_commands) < NumFiles) {
    cat(
      paste0(
        "Fewer commands than the requested number of files.",
        "Adjusting `NumFiles`.\n"))
    NumFiles <- length(LF_commands)
  }

  cat(paste0("Splitting commands into ", NumFiles, " files:\n"))
  LF_commands <- IASDT.R::SplitVector(LF_commands, NSplit = NumFiles)

  # ****************************************************************

  # Write batch files

  purrr::walk(
    .x = seq_len(length(LF_commands)),
    .f = ~ {

      File <- file.path(
        Path_Out,
        paste0(
          "TF_Chunk_",
          stringr::str_pad(.x, pad = "0", width = nchar(NumFiles)), ".txt"))

      cat(paste0("  >>>  ", File, "\n"))

      readr::write_lines(x = BasicCommands, file = File, append = FALSE)
      readr::write_lines(
        x = paste0(
          "# ", length(LF_commands[[.x]]), " commands to be executed:"),
        file = File, append = TRUE)
      readr::write_lines(x = LF_commands[[.x]], file = File, append = TRUE)
      readr::write_lines(
        x = paste0("\n# Remove this file after execution\nrm -f ", File),
        file = File, append = TRUE)
      readr::write_lines(
        x = c(
          "\n#----------------------------------------",
          paste0("# ", lubridate::now(tzone = "CET")),
          "#----------------------------------------"),
        file = File, append = TRUE)
      return(invisible(NULL))
    })

  return(invisible(NULL))
}
