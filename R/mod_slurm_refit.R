## |------------------------------------------------------------------------| #
# mod_slurm_refit ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name mod_slurm
#' @rdname mod_slurm
#' @order 2
#' @export

mod_slurm_refit <- function(
    model_dir = NULL, n_array_jobs = 210L, job_name = NULL,
    memory_per_cpu = "64G", job_runtime = NULL, hpc_partition = "small-g",
    env_file = ".env", cat_job_info = TRUE, ntasks = 1L, cpus_per_task = 1L,
    gpus_per_node = 1L, slurm_prepare = TRUE, path_hmsc = NULL,
    refit_prefix = "commands_to_refit", slurm_prefix = "bash_refit") {

  # # ..................................................................... ###

  # Checking arguments -----

  ecokit::check_args(
    args_to_check = c(
      "hpc_partition", "model_dir", "path_hmsc",
      "slurm_prefix", "refit_prefix"),
    args_type = "character")
  ecokit::check_args(
    args_to_check = c(
      "n_array_jobs", "ntasks", "cpus_per_task", "gpus_per_node"),
    args_type = "numeric")
  ecokit::check_args(
    args_to_check = c("cat_job_info", "slurm_prepare"), args_type = "logical")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  command_hpc <- path_post <- rl <- n <- NULL

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

  path_model_fit <- fs::path(model_dir, "model_fitting_hpc")
  temp_files <- list.files(
    path = path_model_fit, pattern = ".rds_temp$", full.names = TRUE)
  if (length(temp_files) > 0) {
    ecokit::cat_time(
      paste0(
        "There are ", length(temp_files),
        " unsuccessful model variants to be removed"))
    temp_files_rds <- stringr::str_replace_all(temp_files, ".rds_temp$", ".rds")
    purrr::walk(
      .x = c(temp_files_rds, temp_files),
      .f = ~{
        if (file.exists(.x)) {
          file.remove(.x)
        }
      })
  }

  # # ..................................................................... ###

  # List of unfitted model variants -----
  failed_models <- fs::path(model_dir, "model_info.RData") %>%
    ecokit::load_as() %>%
    tidyr::unnest_longer(c(path_post, command_hpc)) %>%
    dplyr::filter(!file.exists(path_post))

  commands_to_refit <- failed_models %>%
    dplyr::pull(command_hpc) %>%
    purrr::set_names(NULL)

  n_jobs <- length(commands_to_refit)

  # # ..................................................................... ###

  # Print summary of unfitted model variants

  if (n_jobs > 0) {
    failed_rl <- dplyr::count(failed_models, rl) %>%
      dplyr::mutate(
        summary = paste0(
          "For Knots spaced by ", rl, " km, there are ", n,
          " failed model variants")) %>%
      dplyr::pull(summary)

    ecokit::cat_time(
      paste0(
        n_jobs, " model variants for ",
        length(unique(failed_models$m_name_fit)),
        " models were not successful."))
    purrr::walk(failed_rl, ecokit::cat_time, level = 1L)
  } else {
    ecokit::cat_time("All models were already fitted!")
  }

  # # ..................................................................... ###

  if (slurm_prepare && n_jobs > 0) {

    ecokit::cat_time("Preparing SLURM script(s) for failed models")

    if (n_jobs > n_array_jobs) {
      n_splits <- ceiling((n_jobs / n_array_jobs))
      ids <- ecokit::split_vector(vector = seq_len(n_jobs), n_splits = n_splits)
    } else {
      n_splits <- 1
      ids <- list(seq_len(n_jobs))
    }

    purrr::walk(
      .x = seq_len(n_splits),
      .f = function(x) {

        curr_ids <- ids[[x]]

        if (n_splits == 1) {
          out_command_file <- paste0(refit_prefix, ".txt")
        } else {
          out_command_file <- paste0(refit_prefix, "_", x, ".txt")
        }

        # create connection to SLURM file
        # This is better than using sink to have a platform independent file
        # (here, to maintain a linux-like new line ending)
        f <- file(fs::path(model_dir, out_command_file), open = "wb")
        on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
        cat(commands_to_refit[curr_ids], sep = "\n", append = FALSE, file = f)
        on.exit(close(f))
      })

    IASDT.R::mod_slurm(
      model_dir = model_dir, job_name = job_name,
      memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
      hpc_partition = hpc_partition, env_file = env_file,
      gpus_per_node = gpus_per_node, cat_job_info = cat_job_info,
      ntasks = ntasks, cpus_per_task = cpus_per_task, path_hmsc = path_hmsc,
      command_prefix = refit_prefix, slurm_prefix = slurm_prefix,
      slurm_path_out = NULL)

    ecokit::cat_time(
      paste0(
        n_jobs, " model variants (", n_splits,
        " SLURM files) need to be re-fitted"), level = 1L)
  }

  # # ..................................................................... ###

  return(invisible(NULL))
}
