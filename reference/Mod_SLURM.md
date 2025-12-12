# Prepare SLURM scripts for Hmsc-HPC model fitting

The `mod_slurm` function generates SLURM job submission scripts for
fitting Hmsc-HPC models in an HPC environment. Additionally,
`mod_slurm_refit` creates SLURM scripts for refitting models that failed
or were not previously fitted.

## Usage

``` r
mod_slurm(
  model_dir = NULL,
  job_name = NULL,
  cat_job_info = TRUE,
  ntasks = 1L,
  cpus_per_task = 1L,
  gpus_per_node = 1L,
  memory_per_cpu = "64G",
  job_runtime = NULL,
  hpc_partition = "small-g",
  env_file = ".env",
  path_hmsc = NULL,
  command_prefix = "commands_to_fit",
  slurm_prefix = "bash_fit",
  slurm_path_out = NULL
)

mod_slurm_refit(
  model_dir = NULL,
  n_array_jobs = 210L,
  job_name = NULL,
  memory_per_cpu = "64G",
  job_runtime = NULL,
  hpc_partition = "small-g",
  env_file = ".env",
  cat_job_info = TRUE,
  ntasks = 1L,
  cpus_per_task = 1L,
  gpus_per_node = 1L,
  slurm_prepare = TRUE,
  path_hmsc = NULL,
  refit_prefix = "commands_to_refit",
  slurm_prefix = "bash_refit"
)
```

## Arguments

- model_dir:

  Character. Path to the root directory of the fitted model.

- job_name:

  Character. Name of the submitted job(s).

- cat_job_info:

  Logical. If `TRUE`, additional bash commands are included to print
  job-related information. Default: `TRUE`.

- ntasks:

  Integer. Number of tasks to allocate for the job (`#SBATCH --ntasks`).
  Default: 1.

- cpus_per_task:

  Integer. Number of CPU cores allocated per task
  (`#SBATCH --cpus-per-task`). Default: 1.

- gpus_per_node:

  Integer. Number of GPUs requested per node
  (`#SBATCH --gpus-per-node`). Default: 1.

- memory_per_cpu:

  Character. Memory allocation per CPU core. Example: "32G" for 32
  gigabytes. Defaults to "64G".

- job_runtime:

  Character. Maximum allowed runtime for the job. Example: "01:00:00"
  for one hour. Required — if not provided, the function throws an
  error.

- hpc_partition:

  Character. Name of the SLURM partition to submit the job to. Default:
  "small-g", for running the array jobs on the GPU.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- path_hmsc:

  Character. Path to the Hmsc-HPC installation.

- command_prefix:

  Character.Prefix for the bash commands used in job execution. Default:
  "`commands_to_fit`".

- slurm_prefix:

  Character. Prefix for the generated SLURM script filenames.

- slurm_path_out:

  Character. Directory where SLURM script(s) will be saved. If `NULL`
  (default), the function derives the path from `model_dir`.

- n_array_jobs:

  Integer. Number of jobs per SLURM script file. In LUMI HPC, there is a
  limit of 210 submitted jobs per user for the `small-g` partition. This
  argument is used to split the jobs into multiple SLURM scripts if
  needed. Default: 210. See [LUMI
  documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
  for more details.

- slurm_prepare:

  Logical. Whether to prepare SLURM command files. If `TRUE` (default),
  the SLURM commands will be saved to disk using the mod_slurm function.

- refit_prefix:

  Character. Prefix for files containing commands to refit failed or
  incomplete models.

## Value

This function does not return a value. Instead, it generates and writes
SLURM script files to disk for model fitting and refitting.

## Author

Ahmed El-Gabbas
