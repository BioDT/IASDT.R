## |------------------------------------------------------------------------| #
# Mod_SLURM ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM scripts for Hmsc-HPC model fitting
#'
#' The `Mod_SLURM` function generates SLURM job submission scripts for fitting
#' Hmsc-HPC models in an HPC environment. Additionally, `Mod_SLURM_Refit`
#' creates SLURM scripts for refitting models that failed or were not previously
#' fitted.
#'
#' @param ModelDir Character. Path to the root directory of the fitted model.
#' @param JobName Character. Name of the submitted job(s).
#' @param CatJobInfo Logical. If `TRUE`, additional bash commands are included
#'   to print job-related information. Default: `TRUE`.
#' @param ntasks Integer. Number of tasks to allocate for the job (`#SBATCH
#'   --ntasks`). Default: 1.
#' @param CpusPerTask Integer. Number of CPU cores allocated per task (`#SBATCH
#'   --cpus-per-task`). Default: 1.
#' @param GpusPerNode Integer. Number of GPUs requested per node (`#SBATCH
#'   --gpus-per-node`). Default: 1.
#' @param MemPerCpu Character. Memory allocation per CPU core. Example: "32G"
#'   for 32 gigabytes. Required --- if not provided, the function throws an
#'   error.
#' @param Time Character. Maximum allowed runtime for the job. Example:
#'   "01:00:00" for one hour. Required --- if not provided, the function throws
#'   an error.
#' @param Partition Character. Name of the SLURM partition to submit the job to.
#'   Default: "small-g", for running the array jobs on the GPU.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param Path_Hmsc Character. Path to the Hmsc-HPC installation.
#' @param Command_Prefix Character.Prefix for the bash commands used in job
#'   execution. Default: "`Commands2Fit`".
#' @param SLURM_Prefix Character. Prefix for the generated SLURM script
#'   filenames.
#' @param Path_SLURM_Out Character. Directory where SLURM script(s) will be
#'   saved. If `NULL` (default), the function derives the path from `ModelDir`.
#' @param Refit_Prefix Character. Prefix for files containing commands to refit
#'   failed or incomplete models.
#' @inheritParams Mod_inputs
#' @name Mod_SLURM
#' @rdname Mod_SLURM
#' @order 1
#' @author Ahmed El-Gabbas
#' @return This function does not return a value. Instead, it generates and
#'   writes SLURM script files to disk for model fitting and refitting.
#' @export

Mod_SLURM <- function(
    ModelDir = NULL, JobName = NULL, CatJobInfo = TRUE, ntasks = 1L,
    CpusPerTask = 1L, GpusPerNode = 1L, MemPerCpu = NULL, Time = NULL,
    Partition = "small-g", EnvFile = ".env", Path_Hmsc = NULL,
    Command_Prefix = "Commands2Fit", SLURM_Prefix = "Bash_Fit",
    Path_SLURM_Out = NULL) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjNum <- Path_GPU_Check <- NULL

  if (is.null(ModelDir) || is.null(JobName) || is.null(MemPerCpu) ||
      is.null(Time) || is.null(Path_Hmsc)) {
    stop(
      "ModelDir, JobName, MemPerCpu, Time, and Path_Hmsc cannot be empty",
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (!file.exists(EnvFile)) {
    stop("Environment file not found:", EnvFile, call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "ProjNum", "DP_R_LUMI_gpu", FALSE, FALSE,
    "Path_GPU_Check", "DP_R_LUMI_gpu_check", FALSE, Sys.info()[1] != "Windows")

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "ModelDir", "JobName", "EnvFile", "Time", "MemPerCpu", "Partition",
    "Path_Hmsc", "ProjNum", "Path_GPU_Check",
    "Command_Prefix", "SLURM_Prefix")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (!fs::dir_exists(ModelDir)) {
    stop("Model directory does not exist:", ModelDir, call. = FALSE)
  }

  ListCommands <- list.files(
    ModelDir, pattern = Command_Prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)
  if (NCommandFiles == 0) {
    stop("The file containing the bash commands does not exist", call. = FALSE)
  }

  if (is.null(Path_SLURM_Out)) {
    # This folder was created in the Mod_Prep4HPC function
    Path_SLURM_Out <- IASDT.R::Path(ModelDir, "Model_Fitting_HPC", "JobsLog")
  }

  purrr::walk(
    .x = seq_len(NCommandFiles),
    .f = function(x) {

      if (NCommandFiles == 1) {
        OutFile <- paste0(SLURM_Prefix, ".slurm")
        JobName0 <- JobName
      } else {
        OutFile <- paste0(SLURM_Prefix, "_", x, ".slurm")
        JobName0 <- paste0(JobName, "_", x)
      }
      NJobs <- R.utils::countLines(ListCommands[x])[1]


      # create connection to SLURM file

      # This is better than using sink to have a platform independent file
      # (here, to maintain a linux-like new line ending)
      f <- file(IASDT.R::Path(ModelDir, OutFile), open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)

      # a wrapper function of cat with new line separator
      cat2 <- function(...) {
        # always cat to the file connection
        cat(sep = "\n", file = f, ...)
      }

      # Writing bash commands to text file

      cat2("#!/bin/bash\n")

      cat2("# -----------------------------------------------------------")
      cat2("# Job array configuration")
      cat2("# -----------------------------------------------------------")
      cat2(paste0("#SBATCH --job-name=", JobName0))
      cat2(paste0("#SBATCH --ntasks=", ntasks))
      cat2(paste0(
        "#SBATCH --output=", IASDT.R::Path(Path_SLURM_Out, "%x-%A-%a.out")))
      cat2(paste0(
        "#SBATCH --error=", IASDT.R::Path(Path_SLURM_Out, "%x-%A-%a.out")))
      cat2(paste0("#SBATCH --account=", ProjNum))
      cat2(paste0("#SBATCH --cpus-per-task=", CpusPerTask))
      cat2(paste0("#SBATCH --mem-per-cpu=", MemPerCpu))
      cat2(paste0("#SBATCH --gpus-per-node=", GpusPerNode))
      cat2(paste0("#SBATCH --time=", Time))
      cat2(paste0("#SBATCH --partition=", Partition))
      cat2(paste0("#SBATCH --array=1-", NJobs, "\n"))

      if (CatJobInfo) {
        cat2("# -----------------------------------------------------------")
        cat2("# Job info")
        cat2("# -----------------------------------------------------------")
        cat2('echo "Start time = $(date)"')
        cat2('echo "Submitting directory = "$SLURM_SUBMIT_DIR')
        cat2('echo "working directory = "$PWD')
        cat2('echo "Project name = "$SLURM_JOB_ACCOUNT')
        cat2('echo "Job id = "$SLURM_JOB_ID')
        cat2('echo "Job name = "$SLURM_JOB_NAME')
        cat2('echo "memory per CPU = "$SLURM_MEM_PER_CPU')
        cat2(paste0(
          'echo "The GPU IDs of GPUs in the job allocation ',
          '(if any) = "$SLURM_JOB_GPUS'))
        cat2('echo "Node running the job script = "$SLURMD_NODENAME')
        cat2(paste0(
          'echo "Process ID of the process started for the task',
          ' = "$SLURM_TASK_PID'))
        cat2('echo "Dependency = "$SLURM_JOB_DEPENDENCY')
        cat2('echo "Number of nodes assigned to a job = "$SLURM_NNODES')
        cat2('echo "Number of tasks requested by the job = "$SLURM_NTASKS')
        cat2('echo "Number of cpus per task = "$SLURM_CPUS_PER_TASK')
        cat2('echo "Number of tasks in the array = "$SLURM_ARRAY_TASK_COUNT')
        cat2(
          'echo "Array\'s maximum ID (index) number = "$SLURM_ARRAY_TASK_MIN')
        cat2(
          'echo "Array\'s minimum ID (index) number = "$SLURM_ARRAY_TASK_MAX\n')
      }

      cat2("# -----------------------------------------------------------")
      cat2("# File contains bash commands for model fitting")
      cat2("# -----------------------------------------------------------")
      cat2(paste0("File=", ListCommands[x], "\n"))

      cat2("# -----------------------------------------------------------")
      cat2("# Loading Hmsc-HPC")
      cat2("# -----------------------------------------------------------")
      cat2(paste0("source ", IASDT.R::Path(Path_Hmsc, "setup-env.sh"), "\n"))

      cat2("# -----------------------------------------------------------")
      cat2("# Check GPU")
      cat2("# -----------------------------------------------------------")
      cat2("export TF_CPP_MIN_LOG_LEVEL=3")
      cat2(paste0("PythonCheckGPU=", Path_GPU_Check, "\n"))

      cat2("# -----------------------------------------------------------")
      cat2("# Some checking")
      cat2("# -----------------------------------------------------------")
      cat2("Path_Python=$(which python3)")
      cat2(paste0(
        'echo -e "Some Checking:\\n  >>  Working directory',
        ": $PWD\\n  >>  Python path:       $Path_Python\\n  >>  ",
        'Checking GPU:      $(python3 $PythonCheckGPU)\\n"'))
      cat2("")

      cat2("# -----------------------------------------------------------")
      cat2("# Run array job")
      cat2("# -----------------------------------------------------------")
      cat2("head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash\n")

      cat2('echo "End of program at `date`"\n\n')

      cat2("# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
      cat2("# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
      cat2(paste0(
        "# This script was created on: ",
        format(lubridate::now(tzone = "CET"), format = "%Y-%m-%d %H:%M"),
        " CET"))
      cat2("# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
      cat2("# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")

      # close connection to the file
      close(f)

      # Print the command to submit the job
      cat(paste0("\t sbatch ", IASDT.R::Path(ModelDir, OutFile), "\n"))
  })

  return(invisible(NULL))
}
