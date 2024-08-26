## |------------------------------------------------------------------------| #
# Mod_SLURM ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM Files for Hmsc-HPC Model Fitting
#'
#' This function prepares and writes SLURM file(s) for submitting model fitting
#' jobs on an HPC environment. It dynamically generates bash scripts based on
#' the provided parameters, which are then used to submit jobs to the SLURM
#' workload manager.
#' @param Path_Model String. Path to the model files (without trailing slash).
#' @param JobName String. The name of the submitted job(s).
#' @param CatJobInfo Logical. Add bash lines to print information on the
#'   submitted job. Default: `TRUE`.
#' @param ntasks Integer. The value for the `#SBATCH --ntasks=` SLURM argument.
#'   Default: 1.
#' @param CpusPerTask Integer. The value for the `#SBATCH --cpus-per-task=`
#'   SLURM argument. Default: 1.
#' @param GpusPerNode Integer. The value for the `#SBATCH --gpus-per-node=`
#'   SLURM argument. Default: 1.
#' @param MemPerCpu String. The value for the `#SBATCH --mem-per-cpu=` SLURM
#'   argument. Example: `32G` to request 32 gigabytes.
#' @param Time String. The value for the requested time for each job in the bash
#'   arrays. Example: `01:00:00` to request an hour.
#' @param Partition String. The name of the partition. Default: `small-g`.
#' @param EnvFile String. Path to read the environment variables. Default value:
#'   `.env`.
#' @param FromHPC Logical. Indicates if the operation is being performed from an
#'   HPC environment. This adjusts file paths accordingly. Default: `TRUE`.
#' @param Path_Hmsc String. Path for the Hmsc-HPC.
#' @param Command_Prefix String. Prefix for the bash commands to be executed. Default: `Commands2Fit`.
#' @param SLURM_Prefix String. Prefix for the exported SLURM file.
#' @param Path_SLURM_Out String indicating the directory where the SLURM file(s)
#'   will be saved. Defaults to `NULL`, which means to identify the path from
#'   `Path_Model`.
#' @name Mod_SLURM
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but writes SLURM script files
#'   to the disk.
#' @details The function reads the following environment variables:
#'    - **`LUMI_ProjNum`** for the BioDT LUMI project number.
#'    - **`DP_R_Path_GPU_Check`** for the path of the python for reporting if
#'    the GPU was used in the running SLURM job.
#' @export

Mod_SLURM <- function(
    Path_Model = NULL, JobName = NULL, CatJobInfo = TRUE, ntasks = 1,
    CpusPerTask = 1, GpusPerNode = 1, MemPerCpu = NULL, Time = NULL,
    Partition = "small-g", EnvFile = ".env", FromHPC = TRUE,
    Path_Hmsc = NULL, Command_Prefix = "Commands2Fit",
    SLURM_Prefix = "Bash_Fit", Path_SLURM_Out = NULL) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjNum <- Path_GPU_Check <- NULL

  if (is.null(Path_Model) || is.null(JobName) || is.null(MemPerCpu) ||
      is.null(Time) || is.null(Path_Hmsc)) {
    stop(
      "Path_Model, JobName, MemPerCpu, Time, and Path_Hmsc cannot be empty",
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (!file.exists(EnvFile)) {
    stop(paste("Environment file not found:", EnvFile), call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "ProjNum", "LUMI_ProjNum", FALSE, FALSE,
    "Path_GPU_Check", "DP_R_Path_GPU_Check", FALSE, FromHPC)

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "Path_Model", "JobName", "EnvFile", "Time", "MemPerCpu", "Partition",
    "Path_Hmsc", "ProjNum", "Path_GPU_Check",
    "Command_Prefix", "SLURM_Prefix")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ListCommands <- list.files(
    Path_Model, pattern = Command_Prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)
  if (NCommandFiles == 0) {
    stop("The file containing the bash commands does not exist", call. = FALSE)
  }

  if (is.null(Path_SLURM_Out)) {
    # This folder was created in the Mod_Prep4HPC function
    Path_SLURM_Out <- file.path(Path_Model, "Model_Fitting_HPC", "JobsLog")
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
      f <- file(file.path(Path_Model, OutFile), open = "wb")
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
        "#SBATCH --output=", file.path(Path_SLURM_Out, "Job_%x-%A-%a.out")))
      cat2(paste0(
        "#SBATCH --error=", file.path(Path_SLURM_Out, "Job_%x-%A-%a.out")))
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
        cat2(paste0('echo "The GPU IDs of GPUs in the job allocation ',
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
      cat2(paste0("source ", file.path(Path_Hmsc, "setup-env.sh"), "\n"))

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
  })

  return(invisible(NULL))
}
