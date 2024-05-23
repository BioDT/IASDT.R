# |---------------------------------------------------| #
# Mod_SLURM_Refit ----
# |---------------------------------------------------| #

#' Prepare SLURM file for unfitted models on HPC
#'
#' Prepare SLURM file for unfitted models on HPC
#'
#' @param Path_Model String. Path to the model files (without trailing slash)
#' @param MaxJobCounts Integer. Maximum number of jobs per slurm file
#' @param JobName String. The name of the submitted jobs
#' @param MemPerCpu String. The value for the `#SBATCH --mem-per-cpu=` SLURM argument. Example: "32G" to request 32 gigabyte
#' @param Time String. The value for the requested time for each job in the bash arrays. Example: "01:00:00" to request an hour.
#' @param Partition String. The name of the partition. Default: `small-g`
#' @param Path_EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param CatJobInfo Logical. Add bash lines to report some job information. Default: `TRUE`
#' @param ntasks Integer. The value for the `#SBATCH --ntasks=` SLURM argument. Default: 1
#' @param CpusPerTask Integer. The value for the `#SBATCH --cpus-per-task=` SLURM argument. Default: 1
#' @param GpusPerNode Integer. The value for the `#SBATCH --gpus-per-node=` SLURM argument. Default: 1
#' @param Refit_Prefix String. Prefix for the file containing commands to be re-fitted
#' @param SLURM_Prefix String. Prefix for the exported SLURM file
#' @name Mod_SLURM_Refit
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_SLURM_Refit <- function(
    Path_Model = NULL, MaxJobCounts = 210, JobName = NULL, MemPerCpu = NULL,
    Time = NULL, Partition = "small-g", Path_EnvFile = ".env", CatJobInfo = TRUE,
    ntasks = 1, CpusPerTask = 1, GpusPerNode = 1,
    Refit_Prefix = "Commands2Refit", SLURM_Prefix = "Bash_Refit") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Command_HPC <- Post_Path <- NULL

  if (file.exists(Path_EnvFile)) {
    readRenviron(Path_EnvFile)
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", Path_EnvFile, " was not found")
    stop(MSG)
  }

  if (is.null(JobName)) {
    JobName <- paste0(basename(Path_Model), "_RF")
  }

  # checking arguments
  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Partition", "Path_EnvFile", "Time", "JobName",
             "MemPerCpu", "Path_Model", "SLURM_Prefix", "Refit_Prefix"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("MaxJobCounts", "ntasks", "CpusPerTask", "GpusPerNode"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = c("CatJobInfo"), Type = "logical")

  rm(AllArgs)
  invisible(gc())

  # remove temp files and incomplete RDs files
  Path_Model_Fit <- file.path(Path_Model, "ModelFitting")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)
  if (length(tempFiles) > 0) {
    IASDT.R::CatTime(
      paste0("There are ", length(tempFiles),
             " unsuccessful model variants to be removed"))
    tempFilesRDs <- stringr::str_replace_all(tempFiles, ".rds_temp$", ".rds")
    purrr::walk(
      .x = c(tempFilesRDs, tempFiles),
      .f = ~{ if(file.exists(.x)) file.remove(.x) } )
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
      NSplits <- ceiling(NJobs / MaxJobCounts)
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
        cat(Commands2Refit[CurrIDs], sep = "\n", append = FALSE,
            file = file.path(Path_Model, OutCommandFile))
      })

    IASDT.R::Mod_SLURM(
      Path_Model = Path_Model, JobName = paste0(basename(Path_Model), "_RF"),
      MemPerCpu = MemPerCpu, Time = Time, Partition = Partition,
      Path_EnvFile = Path_EnvFile, GpusPerNode = GpusPerNode,
      CatJobInfo = CatJobInfo, ntasks = ntasks, CpusPerTask = CpusPerTask,
      Command_Prefix = Refit_Prefix, SLURM_Prefix = SLURM_Prefix)

    IASDT.R::CatTime(
      paste0(NJobs, " model variants (in ", NSplits,
             " slurm files) need to be re-fitted"))
  } else {
    IASDT.R::CatTime("All models were already fitted!")
  }
  return(invisible(NULL))
}

