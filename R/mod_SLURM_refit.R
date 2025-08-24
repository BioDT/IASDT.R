## |------------------------------------------------------------------------| #
# mod_SLURM_refit ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name mod_SLURM
#' @rdname mod_SLURM
#' @order 2
#' @export

mod_SLURM_refit <- function(
    model_dir = NULL, n_array_jobs = 210L, job_name = NULL,
    memory_per_cpu = "64G", job_runtime = NULL, HPC_partition = "small-g",
    env_file = ".env", cat_job_info = TRUE, ntasks = 1L, cpus_per_task = 1L,
    gpus_per_node = 1L, SLURM_prepare = TRUE, path_Hmsc = NULL,
    refit_prefix = "Commands2Refit", SLURM_prefix = "Bash_Refit") {

  # # ..................................................................... ###

  # Checking arguments -----

  ecokit::check_args(
    args_to_check = c(
      "HPC_partition", "model_dir", "path_Hmsc",
      "SLURM_prefix", "refit_prefix"),
    args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "n_array_jobs", "ntasks", "cpus_per_task", "gpus_per_node"),
    args_type = "numeric")
  ecokit::check_args(
    args_to_check = c("cat_job_info", "SLURM_prepare"), args_type = "logical")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Command_HPC <- Post_Path <- rL <- n <- Summ <- NULL

  # # ..................................................................... ###

  if (is.null(job_name)) {
    job_name <- paste0(basename(model_dir), "_RF")
  }

  # # ..................................................................... ###

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  # Validate memory_per_cpu and job_runtime
  memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
  job_runtime <- .validate_slurm_runtime(job_runtime)

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files -----

  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "Model directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  Path_Model_Fit <- fs::path(model_dir, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)
  if (length(tempFiles) > 0) {
    ecokit::cat_time(
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
  FailedModels <- fs::path(model_dir, "Model_Info.RData") %>%
    ecokit::load_as() %>%
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

    ecokit::cat_time(
      paste0(
        NJobs, " model variants for ", length(unique(FailedModels$M_Name_Fit)),
        " models were not successful."))
    purrr::walk(Failed_rL, ecokit::cat_time, level = 1L)
  } else {
    ecokit::cat_time("All models were already fitted!")
  }

  # # ..................................................................... ###

  if (SLURM_prepare && NJobs > 0) {

    ecokit::cat_time("Preparing SLURM script(s) for failed models")

    if (NJobs > n_array_jobs) {
      NSplits <- ceiling((NJobs / n_array_jobs))
      IDs <- ecokit::split_vector(vector = seq_len(NJobs), n_splits = NSplits)
    } else {
      NSplits <- 1
      IDs <- list(seq_len(NJobs))
    }

    purrr::walk(
      .x = seq_len(NSplits),
      .f = function(x) {

        CurrIDs <- IDs[[x]]

        if (NSplits == 1) {
          OutCommandFile <- paste0(refit_prefix, ".txt")
        } else {
          OutCommandFile <- paste0(refit_prefix, "_", x, ".txt")
        }

        # create connection to SLURM file
        # This is better than using sink to have a platform independent file
        # (here, to maintain a linux-like new line ending)
        f <- file(fs::path(model_dir, OutCommandFile), open = "wb")
        on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
        cat(Commands2Refit[CurrIDs], sep = "\n", append = FALSE, file = f)
        on.exit(close(f))
      })

    IASDT.R::mod_SLURM(
      model_dir = model_dir, job_name = job_name,
      memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
      HPC_partition = HPC_partition, env_file = env_file,
      gpus_per_node = gpus_per_node, cat_job_info = cat_job_info,
      ntasks = ntasks, cpus_per_task = cpus_per_task, path_Hmsc = path_Hmsc,
      command_prefix = refit_prefix, SLURM_prefix = SLURM_prefix,
      SLURM_path_out = NULL)

    ecokit::cat_time(
      paste0(
        NJobs, " model variants (", NSplits,
        " SLURM files) need to be re-fitted"), level = 1L)
  }

  # # ..................................................................... ###

  return(invisible(NULL))
}
