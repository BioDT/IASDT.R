## |------------------------------------------------------------------------| #
# Mod_SLURM ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM file for model fitting on HPC
#'
#' This function prepares and writes a SLURM file(s) for submitting model fitting jobs on an HPC environment. It dynamically generates bash scripts based on the provided parameters, which are then used to submit jobs to the SLURM workload manager.
#'
#' @param Path_Model String. Path to the model files (without trailing slash).
#' @param JobName String. The name of the submitted job(s).
#' @param CatJobInfo Logical. Add bash lines to print information on the submitted job. Default: `TRUE`.
#' @param ntasks Integer. The value for the `#SBATCH --ntasks=` SLURM argument. Default: 1.
#' @param CpusPerTask Integer. The value for the `#SBATCH --cpus-per-task=` SLURM argument. Default: 1.
#' @param GpusPerNode Integer. The value for the `#SBATCH --gpus-per-node=` SLURM argument. Default: 1.
#' @param MemPerCpu String. The value for the `#SBATCH --mem-per-cpu=` SLURM argument. Example: "32G" to request 32 gigabytes.
#' @param Time String. The value for the requested time for each job in the bash arrays. Example: "01:00:00" to request an hour.
#' @param Partition String. The name of the partition. Default: `small-g`.
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`.
#' @param FromHPC Logical. Indicates if the operation is being performed from an HPC environment. This adjusts file paths accordingly. Default: `TRUE`.
#' @param Path_Hmsc String. Path for the Hmsc-HPC.
#' @param Command_Prefix String. Prefix for the bash commands to be executed.
#' @param SLURM_Prefix String. Prefix for the exported SLURM file.
#' @param Path_SLURM_Out String indicating the directory where the SLURM file(s) will be saved. Defaults to `NULL`, which means to identify the path from `Path_Model`.
#' @name Mod_SLURM
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but writes SLURM script files to the disk.
#' @export

Mod_SLURM <- function(
    Path_Model = NULL, JobName = NULL, CatJobInfo = TRUE, ntasks = 1,
    CpusPerTask = 1, GpusPerNode = 1, MemPerCpu = NULL, Time = NULL,
    Partition = "small-g", EnvFile = ".env", FromHPC = TRUE,
    Path_Hmsc = NULL, Command_Prefix = "Commands_All",
    SLURM_Prefix = "Bash_Fit", Path_SLURM_Out = NULL) {

  if (is.null(Path_Model) || is.null(JobName) || is.null(MemPerCpu)) {
    stop("Path_Model, JobName and MemPerCpu cannot be empty")
  }

  if (is.null(Time) || is.null(Path_Hmsc)) {
    stop("(Time and Path_Hmsc cannot be empty")
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    ProjNum <- Sys.getenv("IASDT_Proj_Number")
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
    Path_GPU_Check <- Sys.getenv("DP_R_Mod_Path_GPU_Check")
  } else {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  # temporarily setting the working directory
  if (FromHPC) {
    InitialWD <- getwd()
    setwd(Path_Scratch)
    on.exit(setwd(InitialWD), add = TRUE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "Path_Model", "JobName", "EnvFile", "Time", "MemPerCpu", "Partition",
    "Path_Hmsc", "ProjNum", "Path_Scratch", "Path_GPU_Check",
    "Command_Prefix", "SLURM_Prefix")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs)
  invisible(gc())

  # Ensure that model path does not contain the scratch path
  Path_Model <- stringr::str_remove(Path_Model, paste0(Path_Scratch, "/"))

  ListCommands <- list.files(
    Path_Model, pattern = Command_Prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)

  if (NCommandFiles == 0) {
    stop("The file containing the bash commands does not exist")
  }

  if (is.null(Path_SLURM_Out)) {
    # This folder was created in the Mod_Prep4HPC function
    Path_SLURM_Out <- file.path(
      Path_Scratch, Path_Model, "Model_Fitting_HPC", "SLURM_Results")
  }

  purrr::walk(
    .x = seq_len(NCommandFiles),
    .f = function(x) {

      if (NCommandFiles == 1) {
        OutFile <- paste0(SLURM_Prefix, ".slurm")
      } else {
        OutFile <- paste0(SLURM_Prefix, "_", x, ".slurm")
      }
      NJobs <- R.utils::countLines(ListCommands[x])[1]


      # create connection to SLURM file
      # This is better than using sink to have a platform independent file (here, to maintain a linux-like new line ending)
      f <- file(file.path(Path_Model, OutFile), open = "wb")
      on.exit(tryCatch(close(f)), add = TRUE)

      # a wrapper function of cat with new line separator
      cat2 <- function(...) {
        # always cat to the file connection
        cat(sep = "\n", file = f, ...)
      }

      # Writing bash commands to text file

      cat2("#!/bin/bash\n")

      cat2("# -----------------------------------------------")
      cat2("# Job array configuration")
      cat2("# -----------------------------------------------")
      cat2(paste0("#SBATCH --job-name=", paste0(JobName, x)))
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

      cat2("# -----------------------------------------------")
      cat2("# Change working directory to scratch")
      cat2("# -----------------------------------------------")
      cat2(paste0("cd ", Path_Scratch, "\n"))


      if (CatJobInfo) {
        cat2("# -----------------------------------------------")
        cat2("# Job info")
        cat2("# -----------------------------------------------")
        cat2('echo "Start time = $(date)"')
        cat2('echo "Submitting directory = "$SLURM_SUBMIT_DIR')
        cat2('echo "working directory = "$PWD')
        cat2('echo "Project name = "$SLURM_JOB_ACCOUNT')
        cat2('echo "Job id = "$SLURM_JOB_ID')
        cat2('echo "Job name = "$SLURM_JOB_NAME')
        cat2('echo "memory per CPU = "$SLURM_MEM_PER_CPU')
        cat2('echo "The GPU IDs of GPUs in the job allocation (if any) = "$SLURM_JOB_GPUS')
        cat2('echo "Node running the job script = "$SLURMD_NODENAME')
        cat2('echo "Process ID of the process started for the task = "$SLURM_TASK_PID')
        cat2('echo "Dependency = "$SLURM_JOB_DEPENDENCY')
        cat2('echo "Number of nodes assigned to a job = "$SLURM_NNODES')
        cat2('echo "Number of tasks requested by the job = "$SLURM_NTASKS')
        cat2('echo "Number of cpus per task = "$SLURM_CPUS_PER_TASK')
        cat2('echo "Number of tasks in the array = "$SLURM_ARRAY_TASK_COUNT')
        cat2('echo "Array\'s maximum ID (index) number = "$SLURM_ARRAY_TASK_MIN')
        cat2('echo "Array\'s minimum ID (index) number = "$SLURM_ARRAY_TASK_MAX\n')
      }

      cat2("# -----------------------------------------------")
      cat2("# File contains bash commands for model fitting")
      cat2("# -----------------------------------------------")
      cat2(paste0("File=", ListCommands[x], "\n"))

      cat2("# -----------------------------------------------")
      cat2("# Loading Hmsc-HPC")
      cat2("# -----------------------------------------------")
      cat2(paste0("source ", file.path(Path_Hmsc, "setup-env.sh"), "\n"))

      cat2("# -----------------------------------------------")
      cat2("# CHECK GPU")
      cat2("# -----------------------------------------------")
      cat2("export TF_CPP_MIN_LOG_LEVEL=3")
      cat2(paste0("PythonCheckGPU=", Path_GPU_Check, "\n"))

      cat2("# -----------------------------------------------")
      cat2("# Some checking")
      cat2("# -----------------------------------------------")
      cat2("Path_Python=$(which python3)")
      cat2('echo -e "Some Checking:\\n  >>  Working directory: $PWD\\n  >>  Python path:       $Path_Python\\n  >>  Checking GPU:      $(python3 $PythonCheckGPU)\\n"')
      cat2("")

      cat2("# -----------------------------------------------")
      cat2("# Run array job")
      cat2("# -----------------------------------------------")
      cat2("head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash\n")

      cat2('echo "End of program at `date`"')

      # close connection to the file
      close(f)
    })

  return(invisible(NULL))
}
