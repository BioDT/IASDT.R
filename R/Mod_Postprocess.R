#' Model pipeline for post-processing fitted Hmsc models
#'
#'
#' These functions post-process fitted Hmsc models on both CPU and GPU. The
#' pipeline is under active development and may change in future updates.
#' Currently, there are three main functions in this script:
#' `Mod_Postprocess_1_CPU()`, `Mod_Prep_TF()`, and `Mod_Postprocess_2_CPU()`.
#' See details for more information.
#' @param ModelDir Character. Path to the root directory of the fitted model.
#' @param GPP_Dist Integer. Distance in *kilometers* between knots for the
#'   selected model.
#' @param Tree Character. Whether a phylogenetic tree was used in the selected
#'   model. Accepts "Tree" (default) or "NoTree".
#' @param Thin,Samples Integer. Thinning value and the number of MCMC samples of
#'   the selected model.
#' @param NCores_VP Integer. Number of cores to use for variance partitioning.
#'   Defaults to 3.
#' @param NumFiles Integer. Number of output batch files to create. Must be less
#'   than or equal to the maximum job limit of the HPC environment.
#' @param WD Character. Optionally sets the working directory in batch scripts
#'   to this path. If `NULL`, the directory remains unchanged.
#' @param Partition_Name Character. Name of the partition to submit the SLURM
#'   jobs to. Default is `small-g`.
#' @param LF_Time,VP_Time Character. Time limit for latent factor prediction and
#'   variance partitioning processing jobs, respectively. Default is `01:00:00`.
#' @param RC_NCores Integer. The number of cores to use for response curve
#'   prediction. Defaults to `8`.
#' @param PlotWidth_Omega,PlotHeight_Omega,PlotWidth_Beta,PlotHeight_Beta
#'   Integer. The width and height of the generated heatmaps of the Omega and
#'   Beta parameters in centimeters.
#' @rdname Mod_postprocessing
#' @name Mod_postprocessing
#' @order 1
#' @inheritParams Predict_Maps
#' @inheritParams Mod_CV_Fit
#' @inheritParams Mod_Merge_Chains
#' @inheritParams Mod_inputs
#' @inheritParams Response_curves
#' @inheritParams Coda_to_tibble
#' @author Ahmed El-Gabbas
#' @export
#' @details
#'
#' **Mod_Postprocess_1_CPU**
#'
#' This function performs the initial post-processing step for habitat-specific
#' fitted models, automating the following tasks:
#'
#' - check unsuccessful models: [Mod_SLURM_Refit]
#' - merge chains and save R objects (fitted model object and coda object) to
#' `qs2` or `RData` files: [Mod_Merge_Chains]
#' - visualize the convergence of all model variants fitted
#' [Convergence_Plot_All]
#' - visualize the convergence of selected model, including plotting
#' Gelman-Rubin-Brooks [PlotGelman] and [Convergence_Plot] for model convergence
#' diagnostics of the `rho`, `alpha`, `omega`, and `beta` parameters.
#' - extract and save model summary: [Mod_Summary]
#' - plotting model parameters: [Mod_Heatmap_Omega], [Mod_Heatmap_Beta]
#' - prepare data for cross-validation and fit initial cross-validated models:
#' [Mod_CV_Fit]
#' - Prepare scripts for GPU processing, including:
#'    - predicting latent factors of the response curves: [RespCurv_PrepData]
#'    - predicting latent factors for new sampling units: [Predict_Maps]
#'    - computing variance partitioning: [VarPar_Compute]
#'
#' <hr>
#'
#' **Mod_Prep_TF**
#'
#' After running `Mod_Postprocess_1_CPU` for all habitat types, this function
#' prepares batch scripts for GPU computations of all habitat types:
#'    - for <u>variance partitioning</u>, the function matches all files with
#' the pattern ` "VP_.+Command.txt"` (created by [VarPar_Compute] and merges
#' their contents into a single file (`TF_postprocess/VP_Commands.txt`). Then,
#' it prepares a SLURM script for variance partitioning computations
#' (`TF_postprocess/VP_SLURM.slurm`).
#'    - for <u>latent factor predictions</u>, the function matches all files
#' with the pattern `"^LF_NewSites_Commands_.+.txt|^LF_RC_Commands_.+txt"` and
#' split their contents into multiple scripts at the `TF_postprocess` directory
#' for processing as a batch job. The function prepares a SLURM script for
#' latent factor predictions (`LF_SLURM.slurm`).
#'
#' This function is tailored for the LUMI HPC environment and assumes that the
#' `tensorflow` module is installed and correctly configured with all required
#' Python packages. On other HPC systems, users may need to modify the function
#' to load a Python virtual environment or install the required dependencies for
#' TensorFlow and related packages.
#'
#' <hr>
#'
#' **Mod_Postprocess_2_CPU**
#'
#' This function continues running the analysis pipeline for post-processing
#' Hmsc by automating the following steps:
#' - process and visualize response curves: [Response_curves]
#' - predict habitat suitability across different climate options:
#' [Predict_Maps]
#' - plot species & SR predictions as JPEG: [Mod_Predict_Plot]
#' - plot latent factors as JPEG: [Mod_Plot_LF]
#' - process and visualize variance partitioning: [VarPar_Compute] and
#' [VarPar_Plot]
#' - compute and visualizing model internal evaluation (explanatory power):
#' [Mod_Eval_Plot]
#'  - initiate post-processing of fitted cross-validated models: prepare
#' commands for latent factor predictions on GPU --- **Ongoing**
#'
#' This function should be run after:
#' - completing `Mod_Postprocess_1_CPU` and `Mod_Prep_TF` on CPU,
#' - running `VP_SLURM.slurm` and `LF_SLURM.slurm` on GPU to process response
#' curves and latent factor predictions (both scripts are generated by
#' `Mod_Prep_TF`).
#' - submitting SLURM jobs for cross-validated model fitting.
#'

## |------------------------------------------------------------------------| #
# Mod_Postprocess_1_CPU ----
## |------------------------------------------------------------------------| #

Mod_Postprocess_1_CPU <- function(
    ModelDir = NULL, Hab_Abb = NULL, NCores = 8L, EnvFile = ".env",
    Path_Hmsc = NULL, MemPerCpu = NULL, Time = NULL, FromJSON = FALSE,
    GPP_Dist = NULL, Tree = "Tree", Samples = 1000L, Thin = NULL,
    NOmega = 1000L, CVName = c("CV_Dist", "CV_Large"), N_Grid = 50L,
    UseTF = TRUE, TF_use_single = FALSE, LF_NCores = NCores,
    LF_Temp_Cleanup = TRUE, LF_Check = FALSE, Temp_Cleanup = TRUE,
    TF_Environ = NULL, Pred_Clamp = TRUE, Fix_Efforts = "q90",
    Fix_Rivers = "q90", Pred_NewSites = TRUE, NCores_VP = 3,
    PlotWidth_Omega = 26, PlotHeight_Omega = 22.5,
    PlotWidth_Beta = 25, PlotHeight_Beta = 35) {

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
      "FromJSON", "Pred_NewSites", "Pred_Clamp", "Temp_Cleanup",
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
      "Invalid Habitat abbreviation. Valid values are:\n >> ",
      toString(ValidHabAbbs), call. = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(
      "Error: Environment file '", EnvFile, "' is invalid or does not exist.",
      call. = FALSE)
  }

  if (!dir.exists(ModelDir)) {
    stop("Model directory: ", ModelDir, " was not found", call. = FALSE)
  }

  if (!(Tree %in% c("Tree", "NoTree"))) {
    stop(
      "Invalid value for Tree argument. Valid values ",
      "are: 'Tree' or 'NoTree'", call. = FALSE)
  }

  if (!all(CVName %in% c("CV_Dist", "CV_Large"))) {
    stop(
      "Invalid value for CVName argument. Valid values ",
      "are: 'CV_Dist' or 'CV_Large'", call. = FALSE)
  }

  # ****************************************************************

  IASDT.R::RecordArgs(
    ExportPath = IASDT.R::Path(ModelDir, "Args_Mod_Postprocess_1_CPU.RData"))

  # ****************************************************************

  LoadedPackages <- paste(
    sort(IASDT.R::LoadedPackages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "\n  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", IASDT.R::CurrOS(),
      "\n  >>> Model root: ", ModelDir,
      "\n  >>> NCores: ", NCores,
      "\n  >>> EnvFile: ", EnvFile,
      "\n  >>> Path_Hmsc: ", Path_Hmsc,
      "\n  >>> SLURM MemPerCpu: ", MemPerCpu,
      "\n  >>> SLURM Time: ", Time,
      "\n  >>> NOmega: ", NOmega,
      "\n  >>> Hab_Abb: ", Hab_Abb,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"),
      "\n  >>> Loaded packages: \n", LoadedPackages, "\n"))

  # ****************************************************************

  # Check unsuccessful models -----
  Ch1("Check unsuccessful models")

  IASDT.R::Mod_SLURM_Refit(
    ModelDir = ModelDir,
    JobName = stringr::str_remove(basename(ModelDir), "Mod_"),
    MemPerCpu = MemPerCpu, Time = Time, EnvFile = EnvFile,
    Path_Hmsc = Path_Hmsc)

  invisible(gc())

  # ****************************************************************

  Ch1("Merge chains and saving RData files")

  IASDT.R::Mod_Merge_Chains(
    ModelDir = ModelDir, NCores = NCores, FromJSON = FromJSON)

  invisible(gc())

  # ****************************************************************

  Ch1("Convergence of all model variants")

  IASDT.R::Convergence_Plot_All(
    ModelDir = ModelDir, NOmega = NOmega, NCores = NCores,
    MarginType = "histogram")

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Path of selected model -----
  Ch1("Path of selected model")

  Path_Model <- IASDT.R::Path(
    ModelDir, "Model_Fitted",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Model.qs2"))

  Path_Coda <- IASDT.R::Path(
    ModelDir, "Model_Coda",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Coda.qs2"))

  cat(
    paste0("Path_Model:\n\t", Path_Model, "\nPath_Coda:\n\t", Path_Coda, "\n"))

  if (!all(file.exists(Path_Model, Path_Coda))) {
    stop("Selected model files not found", call. = FALSE)
  }

  Temp_Dir <- IASDT.R::Path(ModelDir, "TEMP_Pred")

  # ****************************************************************

  # Gelman_Plot -----
  Ch1("Gelman_Plot")

  IASDT.R::PlotGelman(
    Path_Coda = Path_Coda, Alpha = TRUE, Beta = TRUE, Omega = TRUE, Rho = TRUE,
    NOmega = NOmega, EnvFile = EnvFile)

  invisible(gc())

  # ****************************************************************

  # Convergence plots ----
  Ch1("Convergence plots")

  IASDT.R::Convergence_Plot(
    Path_Coda = Path_Coda, Path_Model = Path_Model, EnvFile = EnvFile,
    NOmega = NOmega, NCores = NCores, NRC = c(2, 2),
    Beta_NRC = c(3, 3), MarginType = "histogram")

  invisible(gc())

  # ****************************************************************

  # Model summary ------
  Ch1("Model summary")

  IASDT.R::Mod_Summary(Path_Coda = Path_Coda, EnvFile = EnvFile)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Omega -----
  Ch1("Plotting model parameters - Omega")

  IASDT.R::Mod_Heatmap_Omega(
    Path_Model = Path_Model, SupportLevel = 0.95,
    PlotWidth = PlotWidth_Omega, PlotHeight = PlotHeight_Omega)

  invisible(gc())

  # ****************************************************************

  # Plotting model parameters - Beta -----
  Ch1("Plotting model parameters - Beta")

  IASDT.R::Mod_Heatmap_Beta(
    Path_Model = Path_Model, SupportLevel = 0.95,
    PlotWidth = PlotWidth_Beta, PlotHeight = PlotHeight_Beta)

  invisible(gc())

  # ****************************************************************

  # Prepare input data for cross-validation -------

  Ch1("Prepare input data for cross-validation")
  IASDT.R::Mod_CV_Fit(
    Path_Model = Path_Model, CVName = CVName, EnvFile = EnvFile,
    JobName = paste0("CV_", Hab_Abb), MemPerCpu = MemPerCpu, Time = Time,
    Path_Hmsc = Path_Hmsc)

  invisible(gc())

  # ****************************************************************

  # latent factors of the response curves -------
  Ch1("Prepare scripts for predicting latent factors of the response curves")

  IASDT.R::RespCurv_PrepData(
    Path_Model = Path_Model, N_Grid = N_Grid, NCores = NCores, UseTF = UseTF,
    TF_Environ = TF_Environ, TF_use_single = TF_use_single,
    LF_NCores = LF_NCores, LF_Temp_Cleanup = LF_Temp_Cleanup,
    LF_Check = LF_Check, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    Verbose = TRUE, LF_Commands_Only = TRUE)

  invisible(gc())

  # ****************************************************************

  # latent factors for new sampling units -------
  Ch1("Prepare scripts for predicting latent factors for new sampling units")

  IASDT.R::Predict_Maps(
    Path_Model = Path_Model, Hab_Abb = Hab_Abb, EnvFile = EnvFile,
    NCores = NCores, Pred_Clamp = Pred_Clamp, Fix_Efforts = Fix_Efforts,
    Fix_Rivers = Fix_Rivers, Pred_NewSites = Pred_NewSites, UseTF = UseTF,
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
    InitTime = .StartTime, Prefix = "\nPost-processing using CPU took ")

  return(invisible(NULL))
}


# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# Mod_Prep_TF ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Mod_postprocessing
#' @name Mod_postprocessing
#' @order 2
#' @author Ahmed El-Gabbas

Mod_Prep_TF <- function(
    NumFiles = 210L, EnvFile = ".env", WD = NULL, Partition_Name = "small-g",
    LF_Time = "01:00:00", VP_Time = "01:30:00") {

  # ****************************************************************

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjectID <- Path_Model <- NULL

  # ****************************************************************

  # Check input arguments ----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("LF_Time", "VP_Time", "Partition_Name", "EnvFile"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NumFiles")
  rm(AllArgs, envir = environment())

  if (NumFiles <= 0) {
    stop("`NumFiles` must be a positive integer.", call. = FALSE)
  }

  # # Load environment variables, for project ID
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "ProjectID", "DP_R_LUMI_gpu", FALSE, FALSE,
    "Path_Model", "DP_R_Model_path", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  IASDT.R::RecordArgs(
    ExportPath = IASDT.R::Path(Path_Model, "Args_Mod_Prep_TF.RData"))

  # ****************************************************************

  # Path to store TF commands
  Path_TF <- IASDT.R::Path(Path_Model, "TF_postprocess")
  # Path to store log files
  Path_Log <- IASDT.R::NormalizePath(IASDT.R::Path(Path_TF, "log"))
  fs::dir_create(c(Path_TF, Path_Log))

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  # VARIANCE PARTITIONING ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  # Prepare post-processing data for calculating variance partitioning
  IASDT.R::CatTime(
    "Prepare post-processing data for calculating variance partitioning")

  # Create paths for VP SLURM script and commands
  Path_VP_SLURM <- IASDT.R::Path(Path_TF, "VP_SLURM.slurm")
  Path_VP_Commands <- IASDT.R::Path(Path_TF, "VP_Commands.txt")

  # Merge and organize TensorFlow commands for computing variance partitioning 
  IASDT.R::CatTime(
    "Merge and organize TF TensorFlow for computing variance partitioning",
    Level = 1, Time = FALSE)

  # Find list of files matching the pattern
  VP_InFiles <- list.files(
    path = Path_Model, pattern = "VP_.+Command.txt", recursive = TRUE,
    full.names = TRUE) %>%
    purrr::map(readr::read_lines) %>%
    unlist() %>%
    gtools::mixedsort()

  # Save all VP commands to single file for batch processing
  readr::write_lines(x = VP_InFiles, file = Path_VP_Commands, append = FALSE)

  # ****************************************************************

  # Prepare batch files ----
  IASDT.R::CatTime("Prepare batch files", Level = 1, Time = FALSE)

  VP_Commands <- c(
    "#!/bin/bash",
    "#SBATCH --job-name=VP_TF",
    "#SBATCH --ntasks=1",
    "#SBATCH --ntasks-per-node=1",
    paste0("#SBATCH --account=", ProjectID),
    "#SBATCH --cpus-per-task=1",
    "#SBATCH --mem-per-cpu=200G",
    "#SBATCH --gpus-per-node=1",
    paste0("#SBATCH --time=", VP_Time),
    paste0("#SBATCH --partition=", Partition_Name),
    paste0("#SBATCH --output=", IASDT.R::Path(Path_Log, "%x-%A-%a.out")),
    paste0("#SBATCH --error=", IASDT.R::Path(Path_Log, "%x-%A-%a.out")),
    paste0("#SBATCH --array=1-", length(VP_InFiles)),
    "\n# File containing commands to be executed",
    paste0("File=", Path_VP_Commands),
    "\n# Load TensorFlow module and configure environment",
    "ml use /appl/local/csc/modulefiles",
    "ml tensorflow",
    "export TF_CPP_MIN_LOG_LEVEL=3",
    "export TF_ENABLE_ONEDNN_OPTS=0\n",
    "# Verify GPU availability",
    paste0(
      'python3 -c "import tensorflow as tf; ',
      'print(\\\"Num GPUs Available:\\\", ',
      'len(tf.config.list_physical_devices(\\\"GPU\\\")))"'),
    "\n# Run array job",
    "head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash",
    "\necho End of program at `date`\n",
    paste0("# ", strrep("-", 50)),
    paste0(
      "# This script was created on: ",
      format(lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"),
      " CET"),
    paste0("# ", strrep("-", 50)))

  IASDT.R::CatTime(
    paste0("Writing SLURM script to: `", Path_VP_SLURM, "`"),
    Level = 2, Time = FALSE)

  readr::write_lines(x = VP_Commands, file = Path_VP_SLURM, append = FALSE)
  # Make the file executable
  Sys.chmod(Path_VP_SLURM, mode = "755")

  # ****************************************************************
  # ****************************************************************

  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #
  # LF PREDICTIONS ----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||| # #

  # Prepare post-processing data for LF predictions
  IASDT.R::CatTime("\nPrepare post-processing data for LF predictions")

  # Ensure that the total number of simultaneous jobs (LF + VP) = NumFiles; so
  # both can be run on the same time.
  NumFiles <- NumFiles - length(VP_InFiles)

  # Merge and organize TensorFlow commands for LF predictions ----
  IASDT.R::CatTime(
    paste0(
      "Merge and organize TensorFlow commands for LF predictions ",
      "into a maximum of ", NumFiles, " files"),
    Level = 1, Time = FALSE)

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
    WD <- IASDT.R::NormalizePath(WD, MustWork = TRUE)
    BasicCommands <- c(
      BasicCommands, "# Change to working directory", paste0("cd ", WD), "")
  }

  # Find list of files matching the pattern
  # Regex pattern to match input files
  LF_Pattern <- "^LF_NewSites_Commands_.+.txt|^LF_RC_Commands_.+txt"
  LF_InFiles <- list.files(
    path = Path_Model, pattern = LF_Pattern,
    recursive = TRUE, full.names = TRUE) %>%
    gtools::mixedsort()

  if (length(LF_InFiles) == 0) {
    stop(
      "No files found matching the pattern `", LF_Pattern,
      "` in ", Path_Model, call. = FALSE)
  }

  IASDT.R::CatTime(
    paste0(
      "Found ", length(LF_InFiles), " files matching the pattern `",
      LF_Pattern, "`"),
    Level = 2, Time = FALSE)
  purrr::walk(LF_InFiles, IASDT.R::CatTime, Level = 3, Time = FALSE)

  # Read and merge commands from input files
  LF_commands <- purrr::map(LF_InFiles, readr::read_lines) %>%
    unlist() %>%
    gtools::mixedsort()

  IASDT.R::CatTime(
    paste0(
      "Total number of commands to be executed: ", length(LF_commands)),
    Level = 2, Time = FALSE)

  if (length(LF_commands) < NumFiles) {
    IASDT.R::CatTime(
      paste0(
        "Fewer commands than the requested number of files. ",
        "Setting `NumFiles=", NumFiles, "`."),
      Level = 2, Time = FALSE)
    NumFiles <- length(LF_commands)
  }

  IASDT.R::CatTime(
    paste0("Splitting commands into ", NumFiles, " files"),
    Time = FALSE, Level = 2)
  LF_commands <- IASDT.R::SplitVector(LF_commands, NSplit = NumFiles)

  purrr::walk(
    .x = seq_len(length(LF_commands)),
    .f = ~ {

      File <- IASDT.R::Path(
        Path_TF,
        paste0(
          "TF_Chunk_",
          stringr::str_pad(.x, pad = "0", width = nchar(NumFiles)), ".txt"))

      readr::write_lines(x = BasicCommands, file = File, append = FALSE)
      readr::write_lines(
        x = paste0(
          "# ", length(LF_commands[[.x]]), " commands to be executed:"),
        file = File, append = TRUE)
      readr::write_lines(x = LF_commands[[.x]], file = File, append = TRUE)
      readr::write_lines(
        x = c(
          "\n#----------------------------------------",
          paste0("# ", lubridate::now(tzone = "CET")),
          "#----------------------------------------"),
        file = File, append = TRUE)
      return(invisible(NULL))
    })

  # ****************************************************************

  # Prepare LF batch files ----
  IASDT.R::CatTime("Prepare batch files", Level = 1, Time = FALSE)

  LF_slurm_script <- c(
    "#!/bin/bash",
    "#SBATCH --job-name=PP_LF",
    "#SBATCH --ntasks=1",
    "#SBATCH --ntasks-per-node=1",
    paste0("#SBATCH --account=", ProjectID),
    "#SBATCH --cpus-per-task=1",
    "#SBATCH --mem-per-cpu=200G",
    "#SBATCH --gpus-per-node=1",
    paste0("#SBATCH --time=", LF_Time),
    paste0("#SBATCH --partition=", Partition_Name),
    paste0("#SBATCH --output=", IASDT.R::Path(Path_Log, "%x-%A-%a.out")),
    paste0("#SBATCH --error=", IASDT.R::Path(Path_Log, "%x-%A-%a.out")),
    paste0("#SBATCH --array=1-", NumFiles),
    "",
    "# Define directories",
    paste0('OutputDir="', IASDT.R::Path(Path_Model, "TF_postprocess"), '"'),
    "",
    "# Find all the split files and sort them explicitly",
    paste0(
      'SplitFiles=($(find "$OutputDir" -type f ',
      '-name "TF_Chunk_*.txt" | sort -V))'),
    "",
    "# Check if files were found",
    "if [ ${#SplitFiles[@]} -eq 0 ]; then",
    '    echo "Error: No files matching TF_Chunk_*.txt found in $OutputDir"',
    "    exit 1",
    "fi",
    "",
    paste0("# Ensure no more than `, NumFiles, ` files are processed"),
    paste0("MaxFiles=", NumFiles),
    "if [ ${#SplitFiles[@]} -gt $MaxFiles ]; then",
    '    SplitFiles=("${SplitFiles[@]:0:$MaxFiles}")',
    paste0(
      '    echo "More than $MaxFiles files found, ',
      'limiting to the first $MaxFiles files."'),
    "fi",
    "",
    "# Get the index of the current task based on SLURM_ARRAY_TASK_ID",
    "TaskIndex=$((SLURM_ARRAY_TASK_ID - 1))",
    "",
    "# Validate TaskIndex",
    "if [ $TaskIndex -ge ${#SplitFiles[@]} ] || [ $TaskIndex -lt 0 ]; then",
    paste0(
      '    echo "Error: TaskIndex $TaskIndex is out of range. ',
      'Valid range: 0 to $((${#SplitFiles[@]} - 1))"'),
    "    exit 1",
    "fi",
    "",
    "# Get the specific split file to process based on the job array task ID",
    'SplitFile="${SplitFiles[$TaskIndex]}"',
    "",
    "# Verify the selected split file",
    'if [ -z "$SplitFile" ] || [ ! -f "$SplitFile" ]; then',
    '    echo "Error: File $SplitFile does not exist or is invalid."',
    "    exit 1",
    "fi",
    "",
    "# Processing file",
    'echo "Processing file: $SplitFile"',
    "",
    "# Run the selected split file",
    'bash "$SplitFile"',
    "\necho End of program at `date`\n",
    paste0("# ", strrep("-", 50)),
    paste0(
      "# This script was created on: ",
      format(lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"), " CET"),
    paste0("# ", strrep("-", 50)))

  Path_LF_SLURM <- IASDT.R::Path(Path_TF, "LF_SLURM.slurm")
  IASDT.R::CatTime(
    paste0("Writing SLURM script to: `", Path_LF_SLURM, "`"),
    Level = 2, Time = FALSE)
  # Write the content to a file
  readr::write_lines(LF_slurm_script, Path_LF_SLURM, append = FALSE)
  # Make the file executable
  Sys.chmod(Path_LF_SLURM, mode = "755")

  # ****************************************************************
  # ****************************************************************

  return(invisible(NULL))
}

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# Mod_Postprocess_2_CPU ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Mod_postprocessing
#' @name Mod_postprocessing
#' @order 3
#' @author Ahmed El-Gabbas

Mod_Postprocess_2_CPU <- function(
    ModelDir = NULL, Hab_Abb = NULL, NCores = 8L, EnvFile = ".env",
    GPP_Dist = NULL, Tree = "Tree", Samples = 1000L, Thin = NULL, UseTF = TRUE,
    TF_Environ = NULL, TF_use_single = FALSE, LF_NCores = NCores,
    LF_Check = FALSE, LF_Temp_Cleanup = TRUE, Temp_Cleanup = TRUE, N_Grid = 50L,
    CC_Models = c(
      "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
      "MRI-ESM2-0", "UKESM1-0-LL"),
    CC_Scenario = c("ssp126", "ssp370", "ssp585"),
    RC_NCores = 8L, Pred_Clamp = TRUE, Fix_Efforts = "q90", Fix_Rivers = "q90",
    Pred_NewSites = TRUE, Tar = TRUE) {

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
    Args = c("Hab_Abb", "EnvFile", "ModelDir", "Tree"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("UseTF", "Pred_Clamp", "Pred_NewSites"))

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NCores", "GPP_Dist", "Samples", "Thin", "N_Grid"))
  rm(AllArgs, envir = environment())

  ValidHabAbbs <- c(as.character(0:3), "4a", "4b", "10", "12a", "12b")
  if (!(Hab_Abb %in% ValidHabAbbs)) {
    stop(
      "Invalid Habitat abbreviation. Valid values are:\n >> ",
      toString(ValidHabAbbs), call. = FALSE)
  }

  if (!file.exists(EnvFile)) {
    stop(
      "Error: Environment file '", EnvFile, "' is invalid or does not exist.",
      call. = FALSE)
  }

  if (!dir.exists(ModelDir)) {
    stop(
      "Model directory: `", ModelDir, "' is invalid or does not exist.",
      call. = FALSE)
  }

  if (!(Tree %in% c("Tree", "NoTree"))) {
    stop(
      "Invalid value for Tree argument. Valid values ",
      "are: 'Tree' or 'NoTree'", call. = FALSE)
  }

  if (Pred_Clamp && is.null(Fix_Efforts)) {
    stop(
      "`Fix_Efforts` can not be NULL when Clamping is implemented",
      call. = FALSE)
  }

  ValidModels <- c(
    "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
    "MRI-ESM2-0", "UKESM1-0-LL")

  if (!all(CC_Models %in% ValidModels)) {
    stop(
      "Invalid climate models. Valid values are:\n >> ",
      toString(ValidModels), call. = FALSE)
  }

  if (!all(CC_Scenario %in% c("ssp126", "ssp370", "ssp585"))) {
    stop(
      "Invalid climate scenarios. Valid values are:\n >> ",
      toString(c("ssp126", "ssp370", "ssp585")), call. = FALSE)
  }

  # ****************************************************************

  IASDT.R::RecordArgs(
    ExportPath = IASDT.R::Path(ModelDir, "Args_Mod_Postprocess_2_CPU.RData"))

  # ****************************************************************

  LoadedPackages <- paste(sort(IASDT.R::LoadedPackages()), collapse = " + ") %>%
    stringr::str_wrap(width = 60, indent = 8, exdent = 8)
  cat(
    paste0(
      "\n  >>> Working directory: ", getwd(),
      "\n  >>> Operating system: ", IASDT.R::CurrOS(),
      "\n  >>> Model root: ", ModelDir,
      "\n  >>> NCores: ", NCores,
      "\n  >>> RC_NCores: ", RC_NCores,
      "\n  >>> LF_NCores: ", LF_NCores,
      "\n  >>> EnvFile: ", EnvFile,
      "\n  >>> Hab_Abb: ", Hab_Abb,
      "\n  >>> UseTF: ", UseTF,
      "\n  >>> Python environment: ", TF_Environ,
      "\n  >>> .libPaths(): \n",
      paste0("\t", .libPaths(), collapse = "\n"),
      "\n  >>> Loaded packages: \n", LoadedPackages, "\n"))

  # ****************************************************************

  Temp_Dir <- IASDT.R::Path(ModelDir, "TEMP_Pred")

  ModelData <- list.files(
    path = ModelDir, full.names = TRUE, pattern = "^ModDT_.*subset.RData")

  if (length(ModelData) != 1) {
    stop(
      "Expected one model data file, but found ", length(ModelData), ".",
      call. = FALSE)
  }

  # ****************************************************************

  # Path of selected model -----
  Ch1("Path of selected model")

  Path_Model <- IASDT.R::Path(
    ModelDir, "Model_Fitted",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Model.qs2"))

  Path_Coda <- IASDT.R::Path(
    ModelDir, "Model_Coda",
    paste0(
      "GPP", GPP_Dist, "_", Tree, "_samp", Samples, "_th", Thin, "_Coda.qs2"))

  cat(
    paste0("Path_Model:\n\t", Path_Model, "\nPath_Coda:\n\t", Path_Coda, "\n"))

  if (!all(file.exists(Path_Model, Path_Coda))) {
    stop("Selected model files not found", call. = FALSE)
  }

  # ****************************************************************

  # Prepare response curve data -----
  Ch1("Prepare response curve data")

  IASDT.R::RespCurv_PrepData(
    Path_Model = Path_Model, N_Grid = N_Grid, NCores = RC_NCores,
    UseTF = UseTF, TF_Environ = TF_Environ, TF_use_single = TF_use_single,
    LF_NCores = LF_NCores, LF_Temp_Cleanup = LF_Temp_Cleanup,
    LF_Check = LF_Check, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    Verbose = TRUE, LF_Commands_Only = FALSE,
    ReturnData = FALSE, Probabilities = c(0.025, 0.5, 0.975))

  invisible(gc())

  # ****************************************************************

  # Plotting response curves - species richness -----
  Ch1("Plotting response curves - species richness")

  IASDT.R::RespCurv_PlotSR(
    ModelDir = ModelDir, Verbose = TRUE, NCores = RC_NCores)

  invisible(gc())

  # ****************************************************************

  # ## Plotting response curves - species -----
  Ch1("Plotting response curves - species")
  IASDT.R::RespCurv_PlotSp(
    ModelDir = ModelDir, NCores = RC_NCores, EnvFile = EnvFile)

  invisible(gc())

  # ****************************************************************

  # ## Plotting - all species together -----
  Ch1("Plotting response curves - all species together")

  IASDT.R::RespCurv_PlotSpAll(ModelDir = ModelDir, NCores = RC_NCores)

  invisible(gc())

  # ****************************************************************

  # Predicting habitat suitability across different climate options -------
  Ch1("Predicting habitat suitability across different climate options")

  IASDT.R::Predict_Maps(
    Path_Model = Path_Model, Hab_Abb = Hab_Abb, EnvFile = EnvFile,
    NCores = NCores, Pred_Clamp = Pred_Clamp, Fix_Efforts = Fix_Efforts,
    Fix_Rivers = Fix_Rivers, Pred_NewSites = Pred_NewSites, UseTF = UseTF,
    TF_Environ = TF_Environ, TF_use_single = TF_use_single,
    LF_NCores = LF_NCores, LF_Check = LF_Check,
    LF_Temp_Cleanup = LF_Temp_Cleanup, LF_Only = FALSE,
    LF_Commands_Only = FALSE, Temp_Dir = Temp_Dir, Temp_Cleanup = Temp_Cleanup,
    Tar = Tar, CC_Models = CC_Models, CC_Scenario = CC_Scenario)

  invisible(gc())

  # ****************************************************************

  # Compute Variance partitioning ------
  Ch1("Compute Variance partitioning")

  IASDT.R::VarPar_Compute(
    Path_Model = Path_Model,
    NCores = 1, UseTF = UseTF, TF_Environ = TF_Environ,
    TF_use_single = TF_use_single, Temp_Cleanup = Temp_Cleanup,
    Chunk_size = 50L, Verbose = TRUE, VarParFile = "VarPar",
    VP_Commands_Only = FALSE)

  # ****************************************************************

  # Plot Variance partitioning ------
  Ch1("Plot Variance partitioning")

  IASDT.R::VarPar_Plot(
    Path_Model = Path_Model, EnvFile = EnvFile, VarParFile = "VarPar",
    UseTF = UseTF, TF_Environ = TF_Environ, NCores = NCores, Fig_width = 30,
    Fig_height = 15)

  # ****************************************************************

  # Plot species & SR predictions as JPEG ------
  Ch1("Plot species & SR predictions as JPEG")

  IASDT.R::Mod_Predict_Plot(
    ModelDir = ModelDir, EnvFile = EnvFile, NCores = NCores)

  # ****************************************************************

  # Plot latent factors as JPEG ------
  Ch1("Plot latent factors as JPEG")

  IASDT.R::Mod_Plot_LF(Path_Model = Path_Model, EnvFile = EnvFile)

  # ****************************************************************

  # Plot explanatory Power ------
  Ch1("Plot explanatory Power")

  IASDT.R::Mod_Eval_Plot(ModelDir = ModelDir, EnvFile = EnvFile)

  # ****************************************************************


  # Post-processing cross-validated models ------

  # ****************************************************************

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "\nPost-processing took ")

  return(invisible(NULL))
}
