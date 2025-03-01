## |------------------------------------------------------------------------| #
# Mod_SLURM_Refit ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Mod_SLURM
#' @rdname Mod_SLURM
#' @order 2
#' @export

Mod_SLURM_Refit <- function(
    ModelDir = NULL, NumArrayJobs = 210L, JobName = NULL, MemPerCpu = NULL,
    Time = NULL, Partition = "small-g", EnvFile = ".env", CatJobInfo = TRUE,
    ntasks = 1L, CpusPerTask = 1L, GpusPerNode = 1L, PrepSLURM = TRUE,
    Path_Hmsc = NULL, Refit_Prefix = "Commands2Refit",
    SLURM_Prefix = "Bash_Refit") {

  # # ..................................................................... ###

  NullVarsNames <- c("ModelDir", "MemPerCpu", "Time", "Path_Hmsc")
  NullVars <- which(purrr::map_lgl(.x = NullVarsNames, .f = ~ is.null(get(.x))))

  if (length(NullVars) > 0) {
    stop(
      toString(NullVarsNames[NullVars]), " cannot be missing.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Command_HPC <- Post_Path <- rL <- n <- Summ <- NULL

  # # ..................................................................... ###

  if (!file.exists(EnvFile)) {
    stop(
      "Path to environment variables: ", EnvFile, " was not found",
      call. = FALSE)
  }

  if (is.null(JobName)) {
    JobName <- paste0(basename(ModelDir), "_RF")
  }

  # # ..................................................................... ###

  # Checking arguments -----

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(.x = AllArgs, .f = get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c(
      "Partition", "EnvFile", "Time", "JobName", "MemPerCpu",
      "ModelDir", "Path_Hmsc", "SLURM_Prefix", "Refit_Prefix"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NumArrayJobs", "ntasks", "CpusPerTask", "GpusPerNode"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = "CatJobInfo", Type = "logical")
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files -----

  if (!fs::dir_exists(ModelDir)) {
    stop("Model directory does not exist:", ModelDir, call. = FALSE)
  }

  Path_Model_Fit <- IASDT.R::Path(ModelDir, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)
  if (length(tempFiles) > 0) {
    IASDT.R::CatTime(
      paste0(
        "There are ", length(tempFiles),
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

  # # ..................................................................... ###

  # List of unfitted model variants -----
  FailedModels <- IASDT.R::Path(ModelDir, "Model_Info.RData") %>%
    IASDT.R::LoadAs() %>%
    tidyr::unnest_longer(c(Post_Path, Command_HPC)) %>%
    dplyr::filter(!file.exists(Post_Path))

  Commands2Refit <- FailedModels %>%
    dplyr::pull(Command_HPC) %>%
    purrr::set_names(NULL)

  NJobs <- length(Commands2Refit)

  # # ..................................................................... ###

  # Print summary of unfitted model variants

  if (NJobs > 0) {
    Failed_rL <- dplyr::count(FailedModels, rL) %>%
      dplyr::mutate(
        Summ = paste0(
          "For Knots spaced by ", rL, " km, there are ", n,
          " failed model variants")) %>%
      dplyr::pull(Summ)

    IASDT.R::CatTime(
      paste0(
        NJobs, " model variants for ", length(unique(FailedModels$M_Name_Fit)),
        " models were not successful."))
    purrr::walk(Failed_rL, IASDT.R::CatTime, Level = 1)
  } else {
    IASDT.R::CatTime("All models were already fitted!")
  }

  # # ..................................................................... ###

  if (PrepSLURM && NJobs > 0) {

    IASDT.R::CatTime("Preparing SLURM script(s) for failed models")

    if (NJobs > NumArrayJobs) {
      NSplits <- ceiling((NJobs / NumArrayJobs))
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
        f <- file(IASDT.R::Path(ModelDir, OutCommandFile), open = "wb")
        on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
        cat(Commands2Refit[CurrIDs], sep = "\n", append = FALSE, file = f)
        close(f)
      })

    IASDT.R::Mod_SLURM(
      ModelDir = ModelDir, JobName = JobName,
      MemPerCpu = MemPerCpu, Time = Time, Partition = Partition,
      EnvFile = EnvFile, GpusPerNode = GpusPerNode,
      CatJobInfo = CatJobInfo, ntasks = ntasks, CpusPerTask = CpusPerTask,
      Path_Hmsc = Path_Hmsc,
      Command_Prefix = Refit_Prefix, SLURM_Prefix = SLURM_Prefix,
      Path_SLURM_Out = NULL)

    IASDT.R::CatTime(
      paste0(
        NJobs, " model variants (", NSplits,
        " SLURM files) need to be re-fitted"), Level = 1)
  }

  # # ..................................................................... ###

  return(invisible(NULL))
}
