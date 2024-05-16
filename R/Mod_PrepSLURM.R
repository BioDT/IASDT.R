# |---------------------------------------------------| #
# PrepSLURM ----
# |---------------------------------------------------| #

#' PrepSLURM
#'
#' PrepSLURM
#' @param Path_Model Path_Model
#' @param CommandFile CommandFile
#' @param JobName JobName
#' @param CatJobInfo CatJobInfo
#' @param ntasks ntasks
#' @param CpusPerTask CpusPerTask
#' @param GpusPerNode GpusPerNode
#' @param MemPerCpu MemPerCpu
#' @param Time Time
#' @param Partition Partition
#' @param Path_EnvFile Path_EnvFile
#' @name PrepSLURM
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PrepSLURM <- function(
    Path_Model = NULL, CommandFile = "Commands_All.txt",
    JobName = NULL, CatJobInfo = TRUE,
    ntasks = 1, CpusPerTask = 1, GpusPerNode = 1, MemPerCpu = NULL,
    Time = NULL, Partition = "small-g", Path_EnvFile = ".env") {

  # # |||||||||||||||||||||||||||||||||||
  # # cat2
  # # |||||||||||||||||||||||||||||||||||

  # a wrapper function of cat with new line separator

  cat2 <- function(...) cat(sep = "\n", ...)

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (file.exists(Path_EnvFile)) {
    readRenviron(Path_EnvFile)
  } else {
    MSG <- paste0(
      "Path for environment variables: ", Path_EnvFile, " was not found")
    stop(MSG)
  }

  Path_Hmsc <- Sys.getenv("DP_R_Mod_Path_Hmsc")
  Path_Python <- Sys.getenv("DP_R_Mod_Path_Python")
  ProjNum <- Sys.getenv("IASDT_Proj_Number")
  Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
  Path_GPU_Check <- Sys.getenv("DP_R_Mod_Path_GPU_Check")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # character arguments
  MissingArgs <- list(
    Path_Model = Path_Model, JobName = JobName,
    Path_EnvFile = Path_EnvFile,
    Time = Time, MemPerCpu = MemPerCpu, Partition = Partition,
    Path_Hmsc = Path_Hmsc, Path_Python = Path_Python,
    ProjNum = ProjNum, Path_Scratch = Path_Scratch,
    Path_GPU_Check = Path_GPU_Check
  ) %>%
    purrr::map(~inherits(.x, "character") && nchar(.x) > 0) %>%
    purrr::discard(.p = isTRUE) %>%
    names() %>%
    sort()
  if (length(MissingArgs) > 0) {
    MSG <- paste0("The following argument(s) must be provided\n  >>  ",
                  paste0(MissingArgs, collapse = " | "))
    stop(MSG)
  }

  # numeric arguments
  ArgsInt <- list(
    GpusPerNode = GpusPerNode, CpusPerTask = CpusPerTask, ntasks = ntasks) %>%
    purrr::map(~inherits(.x, "numeric")) %>%
    purrr::discard(.p = isTRUE) %>%
    names() %>%
    sort()
  if (length(ArgsInt) > 0) {
    MSG <- paste0(
      "The following argument(s) must be numeric (integer)\n  >>  ", ArgsInt)
    stop(MSG)
  }

  CommandFile <- file.path(Path_Model, CommandFile)
  if (file.exists(CommandFile)) {
    NJobs <- R.utils::countLines(CommandFile)[1]
  } else {
    stop("The file containing the bas commands does not exist")
  }

  sink(file = file.path(Path_Model, "BashCommand.slurm"))


  cat2("#!/bin/bash\n")

  cat2("# -----------------------------------------------")
  cat2("# Job array configuration")
  cat2("# -----------------------------------------------")
  cat2(paste0("#SBATCH --job-name=", JobName))
  cat2(paste0("#SBATCH --ntasks=", ntasks))
  cat2(paste0("#SBATCH --output=", file.path(Path_Model, "Job_%x-%A-%a.out")))
  cat2(paste0("#SBATCH --error=", file.path(Path_Model, "Job_%x-%A-%a.out")))
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
    cat2('echo "--mem = "$SLURM_MEM_PER_NODE')
    cat2('echo "--mem-per-cpu = "$SLURM_MEM_PER_CPU')
    cat2('echo "GPU IDs allocated to the job = "$SLURM_JOB_GPUS')
    cat2('echo "Node running the job script = "$SLURMD_NODENAME')
    cat2('echo "Node list used = "$SLURM_NODELIST')
    cat2('echo "Process ID of the task = "$SLURM_TASK_PID')
    cat2('echo "Dependency = "$SLURM_JOB_DEPENDENCY')
    cat2('echo "# cpus per task = "$SLURM_CPUS_PER_TASK')
    cat2('echo "# Nodes allocated = "$SLURM_NNODES')
    cat2('echo "# tasks Allocated = "$SLURM_NTASKS')
    cat2('echo "# tasks in a job array = "$SLURM_ARRAY_TASK_COUNT')
    cat2('echo "Array\'s master job ID number = "$SLURM_ARRAY_JOB_ID')
    cat2('echo "Array\'s maximum ID (index) number = "$SLURM_ARRAY_TASK_MIN')
    cat2('echo "Array\'s minimum ID (index) number = "$SLURM_ARRAY_TASK_MAX')

    cat("\n")
  }

  cat2("# -----------------------------------------------")
  cat2("# Change working directory to scratch")
  cat2("# -----------------------------------------------")
  cat2(paste0("cd ", Path_Scratch, "\n"))

  cat2("# -----------------------------------------------")
  cat2("# File contains bash commands for model fitting")
  cat2("# -----------------------------------------------")
  cat2(paste0("File=", CommandFile, "\n"))

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
  cat2('echo -e "Some Checing:\\n  >>  Time now:          $(date -u)\\n  >>  Working directory: $PWD\\n  >>  Python path:       $Path_Python\\n  >>  Checking GPU:      $(python3 $PythonCheckGPU)\\n"')
  cat("\n")

  cat2("# -----------------------------------------------")
  cat2("# Run array job")
  cat2("# -----------------------------------------------")
  cat2("head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash\n")

  cat2('echo "End of program at `date`"')
  sink()
}
