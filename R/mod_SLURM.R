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

  # character arguments
  ecokit::check_args(
    args_to_check = c(
      "model_dir", "job_name", "HPC_partition", "path_Hmsc",
      "command_prefix", "SLURM_prefix"),
    args_type = "character")
  # numeric arguments
  ecokit::check_args(
    args_to_check = c("gpus_per_node", "cpus_per_task", "ntasks"),
    args_type = "numeric")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ProjNum <- Path_GPU_Check <- NULL

  # Validate memory_per_cpu and job_runtime
  memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
  job_runtime <- .validate_slurm_runtime(job_runtime)

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "ProjNum", "DP_R_lumi_gpu", FALSE, FALSE,
    "Path_GPU_Check", "DP_R_lumi_gpu_check", FALSE, Sys.info()[1] != "Windows")

  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  ListCommands <- list.files(
    model_dir, pattern = command_prefix, full.names = TRUE)
  NCommandFiles <- length(ListCommands)
  if (NCommandFiles == 0) {
    ecokit::stop_ctx(
      "The file containing the bash commands does not exist",
      NCommandFiles = NCommandFiles, ListCommands = basename(ListCommands),
      include_backtrace = TRUE)
  }

  if (is.null(SLURM_path_out)) {
    # This folder was created in the mod_prepare_HPC function
    SLURM_path_out <- fs::path(model_dir, "Model_Fitting_HPC", "JobsLog")
  }

  purrr::walk(
    .x = seq_len(NCommandFiles),
    .f = function(x) {

      SLURM_suffix <- basename(ListCommands[x]) %>%
        stringr::str_remove_all(
          paste(c(command_prefix, "CV_", ".txt"), collapse = "|")) %>%
        stringr::str_remove_all("^_")

      if (nzchar(SLURM_suffix)) {
        OutFile <- c(SLURM_prefix, SLURM_suffix) %>%
          paste(collapse = "_") %>%
          paste0(".slurm")
        JobName0 <- stringr::str_remove_all(OutFile, "_Fit|_Bash|.slurm$")
      } else {
        if (NCommandFiles == 1) { # nolint: unnecessary_nesting_linter
          OutFile <- paste0(SLURM_prefix, ".slurm")
          JobName0 <- job_name
        } else {
          OutFile <- paste0(SLURM_prefix, "_", x, ".slurm")
          JobName0 <- paste0(job_name, "_", x)
        }
      }

      NJobs <- R.utils::countLines(ListCommands[x])[1]

      # create connection to SLURM file

      # This is better than using sink to have a platform independent file
      # (here, to maintain a linux-like new line ending)
      f <- file(fs::path(model_dir, OutFile), open = "wb")
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
        "#SBATCH --output=", fs::path(SLURM_path_out, "%x-%A-%a.out")))
      cat2(paste0(
        "#SBATCH --error=", fs::path(SLURM_path_out, "%x-%A-%a.out")))
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
      cat2(paste0("source ", fs::path(path_Hmsc, "setup-env.sh"), "\n"))

      cat2("# -----------------------------------------------------------")
      cat2("# Check GPU")
      cat2("# -----------------------------------------------------------")
      cat2("export TF_CPP_MIN_LOG_LEVEL=3")
      cat2(paste0("PythonCheckGPU=", Path_GPU_Check, "\n"))

      cat2("# -----------------------------------------------------------")
      cat2("# Some checking")
      cat2("# -----------------------------------------------------------")
      cat2("path_python=$(which python3)")
      cat2(paste0(
        'echo -e "Some Checking:\\n  >>  Working directory',
        ": $PWD\\n  >>  Python path:       $path_python\\n  >>  ",
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
      cat(paste0("\t sbatch ", fs::path(model_dir, OutFile), "\n"))
    })

  return(invisible(NULL))
}
