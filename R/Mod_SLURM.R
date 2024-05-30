## |------------------------------------------------------------------------| #
# Mod_SLURM ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM file for model fitting on HPC
#'
#' Prepare SLURM file for model fitting on HPC
#'
#' @param Path_Model String. Path to the model files (without trailing slash)
#' @param JobName String. The name of the submitted job(s)
#' @param CatJobInfo Logical. Add bash lines to print information on the submitted job. Default: `TRUE`
#' @param ntasks Integer. The value for the `#SBATCH --ntasks=` SLURM argument. Default: 1
#' @param CpusPerTask Integer. The value for the `#SBATCH --cpus-per-task=` SLURM argument. Default: 1
#' @param GpusPerNode Integer. The value for the `#SBATCH --gpus-per-node=` SLURM argument. Default: 1
#' @param MemPerCpu String. The value for the `#SBATCH --mem-per-cpu=` SLURM argument. Example: "32G" to request 32 gigabyte
#' @param Time String. The value for the requested time for each job in the bash arrays. Example: "01:00:00" to request an hour.
#' @param Partition String. The name of the partition. Default: `small-g`
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param FromHPC Logical. Work from HPC? This is to adjust the file paths.
#' @param Command_Prefix String. Prefix for the bash commands to be executed
#' @param SLURM_Prefix String. Prefix for the exported SLURM file.
#' @name Mod_SLURM
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_SLURM <- function(
    Path_Model = NULL, JobName = NULL, CatJobInfo = TRUE, ntasks = 1,
    CpusPerTask = 1, GpusPerNode = 1, MemPerCpu = NULL, Time = NULL,
    Partition = "small-g", EnvFile = ".env", FromHPC = TRUE,
    Command_Prefix = "Commands_All", SLURM_Prefix = "Bash_Fit") {

  InitialWD <- getwd()

  # # |||||||||||||||||||||||||||||||||||
  # # cat2
  # # |||||||||||||||||||||||||||||||||||

  # a wrapper function of cat with new line separator
  cat2 <- function(...) cat(sep = "\n", ...)

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_Hmsc <- Sys.getenv("DP_R_Mod_Path_Hmsc")
    Path_Python <- Sys.getenv("DP_R_Mod_Path_Python")
    ProjNum <- Sys.getenv("IASDT_Proj_Number")
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
    Path_GPU_Check <- Sys.getenv("DP_R_Mod_Path_GPU_Check")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }

  # temporarily setting the working directory
  if (FromHPC) setwd(Path_Scratch)


  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "Path_Model", "JobName", "EnvFile", "Time", "MemPerCpu", "Partition",
    "Path_Hmsc", "Path_Python", "ProjNum", "Path_Scratch", "Path_GPU_Check",
    "Command_Prefix", "SLURM_Prefix")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  # numeric arguments
  NumericArgs <- c("GpusPerNode", "CpusPerTask", "ntasks")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs)
  invisible(gc())

  ListCommands <- list.files(
    Path_Model, pattern = Command_Prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)

  if (NCommandFiles == 0) {
    stop("The file containing the bash commands does not exist")
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

      sink(file = file.path(Path_Model, OutFile))
      Path_SLURM_Out <- file.path(
        Path_Scratch, Path_Model, "Model_Fitting_HPC", "SLURM_Results")
      fs::dir_create(Path_SLURM_Out)

      cat2("#!/bin/bash\n")

      cat2("# -----------------------------------------------")
      cat2("# Job array configuration")
      cat2("# -----------------------------------------------")
      cat2(paste0("#SBATCH --job-name=", paste0(JobName, x)))
      cat2(paste0("#SBATCH --ntasks=", ntasks))
      cat2(paste0("#SBATCH --output=",
                  file.path(Path_SLURM_Out, "Job_%x-%A-%a.out")))
      cat2(paste0("#SBATCH --error=",
                  file.path(Path_SLURM_Out, "Job_%x-%A-%a.out")))
      cat2(paste0("#SBATCH --account=", ProjNum))
      cat2(paste0("#SBATCH --cpus-per-task=", CpusPerTask))
      cat2(paste0("#SBATCH --mem-per-cpu=", MemPerCpu))
      cat2(paste0("#SBATCH --gpus-per-node=", GpusPerNode))
      cat2(paste0("#SBATCH --time=", Time))
      cat2(paste0("#SBATCH --partition=", Partition))
      cat2(paste0("#SBATCH --array=1-", NJobs, "\n"))

      if (CatJobInfo) {
        cat2("# -----------------------------------------------")
        cat2("# Job info")
        cat2("# -----------------------------------------------")
        cat2('echo "Start date = $(date)"')
        cat2('echo "working directory = "$SLURM_SUBMIT_DIR')
        cat2('echo "Account name = "$SLURM_JOB_ACCOUNT')
        cat2('echo "Job id = "$SLURM_JOB_ID')
        cat2('echo "Job name = "$SLURM_JOB_NAME')
        cat2('echo "--mem-per-cpu = "$SLURM_MEM_PER_CPU')
        cat2('echo "GPU IDs allocated to the job = "$SLURM_JOB_GPUS')
        cat2('echo "Node running the job script = "$SLURMD_NODENAME')
        cat2('echo "Process ID of the task = "$SLURM_TASK_PID')
        cat2('echo "Dependency = "$SLURM_JOB_DEPENDENCY')
        cat2('echo "Number of nodes allocated = "$SLURM_NNODES')
        cat2('echo "Number of tasks allocated = "$SLURM_NTASKS')
        cat2('echo "Number of cpus per task = "$SLURM_CPUS_PER_TASK')
        cat2('echo "# tasks in a job array = "$SLURM_ARRAY_TASK_COUNT')
        cat2('echo "Array\'s master job ID number = "$SLURM_ARRAY_JOB_ID')
        cat2('echo "Array\'s maximum ID (index) number = "$SLURM_ARRAY_TASK_MIN')
        cat2('echo "Array\'s minimum ID (index) number = "$SLURM_ARRAY_TASK_MAX')

        cat("\n\n")
      }

      cat2("# -----------------------------------------------")
      cat2("# Change working directory to scratch")
      cat2("# -----------------------------------------------")
      cat2(paste0("cd ", Path_Scratch, "\n"))

      cat2("# -----------------------------------------------")
      cat2("# File contains bash commands for model fitting")
      cat2("# -----------------------------------------------")
      cat2(paste0("File=", ListCommands[x], "\n"))

      cat2("# -----------------------------------------------")
      cat2("# Loading modules")
      cat2("# -----------------------------------------------")
      cat2("module use /appl/local/csc/modulefiles/")
      cat2("module load tensorflow/2.12\n")

      cat2("# -----------------------------------------------")
      cat2("# HMSC-HPC path")
      cat2("# -----------------------------------------------")
      cat2(paste0("Path_Hmsc=", Path_Hmsc, "\n"))

      cat2("# -----------------------------------------------")
      cat2("# Python path")
      cat2("# -----------------------------------------------")
      cat2("export PYTHONPATH=$Path_Hmsc:$PYTHONPATH")
      cat2("Path_Python=$(which python3)\n")

      cat2("# -----------------------------------------------")
      cat2("# CHECK GPU")
      cat2("# -----------------------------------------------")
      cat2("export TF_CPP_MIN_LOG_LEVEL=3")
      cat2(paste0("PythonCheckGPU=", Path_GPU_Check, "\n"))

      cat2("# -----------------------------------------------")
      cat2("# Some checking")
      cat2("# -----------------------------------------------")
      cat2('echo -e "Some Checking:\\n>>  Working directory: $PWD\\n  >>  Python path:       $Path_Python\\n  >>  Checking GPU:      $(python3 $PythonCheckGPU)\\n"')
      cat("\n")

      cat2("# -----------------------------------------------")
      cat2("# Run array job")
      cat2("# -----------------------------------------------")
      cat2("head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash\n")

      cat2('echo "End of program at `date`"')
      sink()
    })

  if (FromHPC) setwd(InitialWD)

  return(invisible(NULL))
}
