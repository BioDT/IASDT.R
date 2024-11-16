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

# NEED TO REVISE DOCUMENTATION

Mod_Postprocess <- function(
    ModelDir = NULL, Hab_Abb = NULL, NCores = 54, FromHPC = TRUE,
    EnvFile = ".env", Path_Hmsc = NULL, MemPerCpu = NULL, Time = NULL,
    FromJSON = FALSE, GPP_Dist = NULL, Tree = "Tree", Samples = 1000,
    Thin = NULL, N_Grid = 50, NOmega = 1000, UseTF = TRUE, TF_Environ = NULL,
    CC_Models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_Scenario = c("ssp126", "ssp370", "ssp585"),
    Pred_Clamp = TRUE, Fix_Efforts = "mean", Pred_NewSites = TRUE,
    CVName = c("CV_Dist", "CV_Large")) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  Ch1 <- function(Text) {
    IASDT.R::InfoChunk(
      paste0("\t", Text), Extra1 = 2, Extra2 = 2, Rep = 2,
      Char = "=", CharReps = 60)
  }

  Ch2 <- function(Text) {
    IASDT.R::InfoChunk(
      paste0("\t", Text), Extra1 = 1, Extra2 = 1, Rep = 1,
      Char = "-", CharReps = 60)
  }

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
    stop("`Fix_Efforts` can not be NULL when Clamping is implemented")
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
      "\n  >>> Loaded packages: \n", LoadedPackages
    ))

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
  IASDT.R::Merge_Chains(ModelDir = ModelDir, NCores = NCores, FromHPC = FromHPC)

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

  # ****************************************************************


  Temp_Dir <- file.path(ModelDir, "TEMP2Pred")

  ModelData <- list.files(
    path = ModelDir, pattern = paste0("^ModDT_", Hab_Abb, "_.+_subset.RData"))

  if (length(ModelData) != 0) {
    stop("Model data was not found", call. = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Path of selected model -----
  Ch1("Path of selected model")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Path_Model <- file.path(
    ModelDir, "Model_Fitted",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Model.qs"))
  Path_Coda <- file.path(
    ModelDir, "Model_Coda",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Coda.qs"))

  cat(paste0("Path_Model:\n\t", Path_Model, "\nPath_Coda:\n\t", Path_Coda))

  if (!all(file.exists(Path_Model, Path_Coda))) {
    stop("Selected model files not found", call. = FALSE)
  }


  # # ..................................................................... ###
  # # ..................................................................... ###


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Convergence ----
  Ch1("Convergence")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Response curves -----
  Ch1("Response curves")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Prepare data ------
  Ch2("Prepare data")

  IASDT.R::RespCurv_PrepData(
    Path_Model = Path_Model, N_Grid = N_Grid, NCores = NCores, UseTF = UseTF,
    TF_Environ = TF_Environ, Temp_Dir = Temp_Dir)

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


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Model summary ------
  Ch1("Model summary")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::Mod_Summary(
    Path_Coda = Path_Coda, EnvFile = EnvFile, FromHPC = FromHPC)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plotting model parameters -----
  Ch1("Plotting model parameters")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Spatial predictions / evaluation ---------
  Ch1("Spatial predictions / evaluation")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::Predict_Maps(
    Path_Model = Path_Model, Hab_Abb = Hab_Abb, EnvFile = EnvFile,
    FromHPC = FromHPC, NCores = NCores, Pred_Clamp = Pred_Clamp,
    Fix_Efforts = Fix_Efforts, Pred_NewSites = Pred_NewSites, UseTF = UseTF,
    TF_Environ = TF_Environ, CC_Models = CC_Models, CC_Scenario = CC_Scenario,
    Temp_Dir = Temp_Dir)

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Variance paritioning ------
  Ch1("Variance paritioning")
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::VarPar_Plot(
    Path_Model = Path_Model, EnvFile = EnvFile, FromHPC = FromHPC,
    UseTF = UseTF, TF_Environ = TF_Environ, NCores = 5,
    Fig_width = 30, Fig_height = 15)

  # # ..................................................................... ###
  # # ..................................................................... ###


  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Cross-validation -------
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  IASDT.R::InfoChunk("Prepare input data for cross-validation")

  IASDT.R::Mod_CV_Fit(
    Model = Path_Model, ModelData = ModelData, CVName = CVName,
    EnvFile = EnvFile, JobName = paste0("CV_", Hab_Abb), FromHPC = FromHPC,
    MemPerCpu = MemPerCpu, Time = Time, Path_Hmsc = Path_Hmsc)


  # # ..................................................................... ###
  # # ..................................................................... ###

  CatTime()
  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPostprocessing took ")

  return(invisible(NULL))

}
