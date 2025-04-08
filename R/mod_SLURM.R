## |------------------------------------------------------------------------| #
# mod_SLURM ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM scripts for Hmsc-HPC model fitting
#'
#' The `mod_SLURM` function generates SLURM job submission scripts for fitting
#' Hmsc-HPC models in an HPC environment. Additionally, `mod_SLURM_refit`
#' creates SLURM scripts for refitting models that failed or were not previously
#' fitted.
#'
#' @param model_dir Character. Path to the root directory of the fitted model.
#' @param job_name Character. Name of the submitted job(s).
#' @param cat_job_info Logical. If `TRUE`, additional bash commands are included
#'   to print job-related information. Default: `TRUE`.
#' @param ntasks Integer. Number of tasks to allocate for the job (`#SBATCH
#'   --ntasks`). Default: 1.
#' @param cpus_per_task Integer. Number of CPU cores allocated per task
#'   (`#SBATCH
#'   --cpus-per-task`). Default: 1.
#' @param gpus_per_node Integer. Number of GPUs requested per node (`#SBATCH
#'   --gpus-per-node`). Default: 1.
#' @param memory_per_cpu Character. Memory allocation per CPU core. Example:
#'   "32G" for 32 gigabytes. Defaults to "64G".
#' @param job_runtime Character. Maximum allowed runtime for the job. Example:
#'   "01:00:00" for one hour. Required --- if not provided, the function throws
#'   an error.
#' @param HPC_partition Character. Name of the SLURM partition to submit the job
#'   to. Default: "small-g", for running the array jobs on the GPU.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param path_Hmsc Character. Path to the Hmsc-HPC installation.
#' @param command_prefix Character.Prefix for the bash commands used in job
#'   execution. Default: "`Commands2Fit`".
#' @param SLURM_prefix Character. Prefix for the generated SLURM script
#'   filenames.
#' @param SLURM_path_out Character. Directory where SLURM script(s) will be
#'   saved. If `NULL` (default), the function derives the path from `model_dir`.
#' @param refit_prefix Character. Prefix for files containing commands to refit
#'   failed or incomplete models.
#' @inheritParams mod_inputs
#' @name mod_SLURM
#' @rdname mod_SLURM
#' @order 1
#' @author Ahmed El-Gabbas
#' @return This function does not return a value. Instead, it generates and
#'   writes SLURM script files to disk for model fitting and refitting.
#' @export

mod_SLURM <- function(
    model_dir = NULL, job_name = NULL, cat_job_info = TRUE, ntasks = 1L,
    cpus_per_task = 1L, gpus_per_node = 1L, memory_per_cpu = "64G",
    job_runtime = NULL, HPC_partition = "small-g", env_file = ".env",
    path_Hmsc = NULL, command_prefix = "Commands2Fit",
    SLURM_prefix = "Bash_Fit", SLURM_path_out = NULL) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjNum <- Path_GPU_Check <- NULL

  if (is.null(model_dir) || is.null(job_name) || is.null(memory_per_cpu) ||
      is.null(job_runtime) || is.null(path_Hmsc)) {
    stop(
      "`model_dir`, `job_name`, `memory_per_cpu`, `job_runtime`, ",
      "and `path_Hmsc` ", "cannot be empty",
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (!file.exists(env_file)) {
    stop("Environment file not found:", env_file, call. = FALSE)
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "ProjNum", "DP_R_LUMI_gpu", FALSE, FALSE,
    "Path_GPU_Check", "DP_R_LUMI_gpu_check", FALSE, Sys.info()[1] != "Windows")

  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  # character arguments
  CharArgs <- c(
    "model_dir", "job_name", "env_file", "job_runtime", "memory_per_cpu",
    "HPC_partition", "path_Hmsc", "ProjNum", "Path_GPU_Check",
    "command_prefix", "SLURM_prefix")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = CharArgs,
    args_type = "character")

  # numeric arguments
  NumericArgs <- c("gpus_per_node", "cpus_per_task", "ntasks")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = NumericArgs,
    args_type = "numeric")

  rm(AllArgs, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (!fs::dir_exists(model_dir)) {
    stop("Model directory does not exist:", model_dir, call. = FALSE)
  }

  ListCommands <- list.files(
    model_dir, pattern = command_prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)
  if (NCommandFiles == 0) {
    stop("The file containing the bash commands does not exist", call. = FALSE)
  }

  if (is.null(SLURM_path_out)) {
    # This folder was created in the mod_prepare_HPC function
    SLURM_path_out <- IASDT.R::path(model_dir, "Model_Fitting_HPC", "JobsLog")
  }

  purrr::walk(
    .x = seq_len(NCommandFiles),
    .f = function(x) {

      if (NCommandFiles == 1) {
        OutFile <- paste0(SLURM_prefix, ".slurm")
        JobName0 <- job_name
      } else {
        OutFile <- paste0(SLURM_prefix, "_", x, ".slurm")
        JobName0 <- paste0(job_name, "_", x)
      }
      NJobs <- R.utils::countLines(ListCommands[x])[1]


      # create connection to SLURM file

      # This is better than using sink to have a platform independent file
      # (here, to maintain a linux-like new line ending)
      f <- file(IASDT.R::path(model_dir, OutFile), open = "wb")
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
        "#SBATCH --output=", IASDT.R::path(SLURM_path_out, "%x-%A-%a.out")))
      cat2(paste0(
        "#SBATCH --error=", IASDT.R::path(SLURM_path_out, "%x-%A-%a.out")))
      cat2(paste0("#SBATCH --account=", ProjNum))
      cat2(paste0("#SBATCH --cpus-per-task=", cpus_per_task))
      cat2(paste0("#SBATCH --mem-per-cpu=", memory_per_cpu))
      cat2(paste0("#SBATCH --gpus-per-node=", gpus_per_node))
      cat2(paste0("#SBATCH --time=", job_runtime))
      cat2(paste0("#SBATCH --partition=", HPC_partition))
      cat2(paste0("#SBATCH --array=1-", NJobs, "\n"))

      if (cat_job_info) {
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
      cat2(paste0("source ", IASDT.R::path(path_Hmsc, "setup-env.sh"), "\n"))

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
      cat(paste0("\t sbatch ", IASDT.R::path(model_dir, OutFile), "\n"))
  })

  return(invisible(NULL))
}
