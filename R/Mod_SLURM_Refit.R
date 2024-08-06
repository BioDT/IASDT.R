## |------------------------------------------------------------------------| #
# Mod_SLURM_Refit ----
## |------------------------------------------------------------------------| #

#' Prepare SLURM file for failed / Unfitted Models
#'
#' This function prepares SLURM files for models that have not been successfully
#' fitted yet. It generates command files for refitting and creates
#' corresponding SLURM batch files to execute these commands on a HPC
#' environment.
#' @param Path_Model String. Path to the model files directory. Must not end
#'   with a slash.
#' @param MaxJobCounts Integer. Maximum number of batch jobs that can be
#'   submitted in a single SLURM file.
#' @param JobName String (optional). Name of the SLURM jobs. If not provided, a
#'   default name based on the model path is used.
#' @param MemPerCpu String. Memory per CPU allocation for the SLURM job.
#'   Example: `32G` for 32 gigabytes.
#' @param Time String. Duration for which the job should run. Example:
#'   `01:00:00` for one hour.
#' @param Partition String. The SLURM partition to submit the job to. Default is
#'   "small-g".
#' @param EnvFile String. Path to the environment variables file. Defaults to
#'   `.env`.
#' @param CatJobInfo Logical. Whether to include commands in the SLURM script to
#'   print job information. Defaults to TRUE.
#' @param ntasks Integer. Number of tasks to request (`#SBATCH --ntasks=`) for
#'   the SLURM job. Defaults to 1.
#' @param CpusPerTask Integer. Number of CPUs per task (`#SBATCH
#'   --cpus-per-task=`) to request. Defaults to 1.
#' @param GpusPerNode Integer. Number of GPUs per node to request (`#SBATCH
#'   --gpus-per-node=`). Defaults to 1.
#' @param FromHPC Logical. Indicates whether the function is being called from
#'   an HPC environment. Adjusts file paths accordingly.
#' @param Path_Hmsc String. Path to the `Hmsc-HPC` directory.
#' @param Refit_Prefix String. Prefix for the files that contain the commands to
#'   refit the models.
#' @param SLURM_Prefix String. Prefix for the generated SLURM batch files.
#' @name Mod_SLURM_Refit
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but creates command and SLURM
#'   batch files for refitting models.
#' @details The function reads the following environment variable:
#'    - **`LUMI_Scratch`** for the path of the scratch folder of the
#'    `BioDT` project on LUMI.
#' @export

Mod_SLURM_Refit <- function(
    Path_Model = NULL, MaxJobCounts = 210, JobName = NULL, MemPerCpu = NULL,
    Time = NULL, Partition = "small-g", EnvFile = ".env", CatJobInfo = TRUE,
    ntasks = 1, CpusPerTask = 1, GpusPerNode = 1, FromHPC = TRUE,
    Path_Hmsc = NULL,
    Refit_Prefix = "Commands2Refit", SLURM_Prefix = "Bash_Refit") {

  if (is.null(Path_Model) || is.null(MemPerCpu) || is.null(Time) ||
      is.null(Path_Hmsc)) {
    stop("Path_Model, MemPerCpu, Time and Path_Hmsc cannot be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Command_HPC <- Post_Path <- Path_Scratch <- NULL

  if (magrittr::not(file.exists(EnvFile))) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Scratch", "LUMI_Scratch", FromHPC, FALSE)

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # temporarily setting the working directory
  if (FromHPC) {
    InitialWD <- getwd()
    setwd(Path_Scratch)
    on.exit(setwd(InitialWD), add = TRUE)
  }

  if (is.null(JobName)) {
    JobName <- paste0(basename(Path_Model), "_RF")
  }

  # checking arguments
  AllArgs <- ls()
  AllArgs <- purrr::map(
    .x = AllArgs,
    .f = ~get(.x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Partition", "EnvFile", "Time", "JobName", "MemPerCpu",
             "Path_Model", "Path_Hmsc", "SLURM_Prefix", "Refit_Prefix"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("MaxJobCounts", "ntasks", "CpusPerTask", "GpusPerNode"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("CatJobInfo"), Type = "logical")
  rm(AllArgs)


  # remove temp files and incomplete RDs files
  Path_Model_Fit <- file.path(Path_Model, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)
  if (length(tempFiles) > 0) {
    IASDT.R::CatTime(
      paste0("There are ", length(tempFiles),
             " unsuccessful model variants to be removed"))
    tempFilesRDs <- stringr::str_replace_all(tempFiles, ".rds_temp$", ".rds")
    purrr::walk(
      .x = c(tempFilesRDs, tempFiles),
      .f = ~{
        if (file.exists(.x)) {
          file.remove(.x)
        }
      })
  }

  # List of unfitted model variants
  Commands2Refit <- file.path(Path_Model, "Model_Info.RData") %>%
    IASDT.R::LoadAs() %>%
    tidyr::unnest_longer(c(Post_Path, Command_HPC)) %>%
    dplyr::select(Post_Path, Command_HPC) %>%
    dplyr::filter(
      magrittr::not(file.exists(file.path(Path_Scratch, Post_Path)))) %>%
    dplyr::pull(Command_HPC) %>%
    purrr::set_names(NULL)

  NJobs <- length(Commands2Refit)

  if (NJobs > 0) {
    if (NJobs > MaxJobCounts) {
      NSplits <- ceiling((NJobs / MaxJobCounts))
      IDs <- IASDT.R::SplitVector(Vector = seq_len(NJobs), NSplit = NSplits)
    } else {
      NSplits <- 1
      IDs <- list(seq_len(NJobs))
    }

    purrr::walk(
      .x = seq_len(NSplits),
      .f = function(x) {

        CurrIDs <- IDs[[x]]

        if (NSplits == 1) {
          OutCommandFile <- paste0(Refit_Prefix, ".txt")
        } else {
          OutCommandFile <- paste0(Refit_Prefix, "_", x, ".txt")
        }

        # create connection to SLURM file
        # This is better than using sink to have a platform independent file
        # (here, to maintain a linux-like new line ending)
        f <- file(file.path(Path_Model, OutCommandFile), open = "wb")
        on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
        cat(Commands2Refit[CurrIDs], sep = "\n", append = FALSE, file = f)
        close(f)
      })

    IASDT.R::Mod_SLURM(
      Path_Model = Path_Model, JobName = paste0(basename(Path_Model), "_RF"),
      MemPerCpu = MemPerCpu, Time = Time, Partition = Partition,
      EnvFile = EnvFile, GpusPerNode = GpusPerNode,
      CatJobInfo = CatJobInfo, ntasks = ntasks, CpusPerTask = CpusPerTask,
      Path_Hmsc = Path_Hmsc,
      Command_Prefix = Refit_Prefix, SLURM_Prefix = SLURM_Prefix,
      FromHPC = FromHPC,
      Path_SLURM_Out = NULL)

    IASDT.R::CatTime(
      paste0(NJobs, " model variants (in ", NSplits,
             " slurm files) need to be re-fitted"))
  } else {
    IASDT.R::CatTime("All models were already fitted!")
  }

  return(invisible(NULL))
}
