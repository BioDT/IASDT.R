# Prepare cross-validated Hmsc models for HPC fitting

This function prepares cross-validated Hmsc models for fitting using
HPC. It handles data preparation, model initialisation, and generation
of SLURM commands.

## Usage

``` r
mod_cv_fit(
  path_model = NULL,
  cv_name = c("cv_dist", "cv_large"),
  partitions = NULL,
  env_file = ".env",
  init_par = NULL,
  job_name = "cv_models",
  updater = list(Gamma2 = FALSE, GammaEta = FALSE),
  align_posterior = TRUE,
  to_json = FALSE,
  slurm_prepare = TRUE,
  memory_per_cpu = NULL,
  job_runtime = NULL,
  path_hmsc = NULL,
  precision = 64,
  ...
)
```

## Arguments

- path_model:

  Character. Path to a saved model file (`*.qs2`).

- cv_name:

  Character vector. Column name(s) in the model input data to be used to
  cross-validate the models (see
  [mod_prepare_data](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  and
  [mod_cv_prepare](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md)).
  The function allows the possibility of using more than one way of
  assigning grid cells into cross-validation folders. If multiple names
  are provided, separate cross-validation models will be fitted for each
  cross-validation type. Currently, there are three cross-validation
  strategies: `cv_sac`, `cv_dist`, and `cv_large`. Defaults to
  `c("cv_dist", "cv_large")`.

- partitions:

  A vector for cross-validation created by
  [Hmsc::createPartition](https://rdrr.io/pkg/Hmsc/man/createPartition.html)
  or similar. Defaults to `NULL`, which means to use column name(s)
  provided in the `cv_name` argument. If the `partitions` vector is
  provided, the label used in the output files will be `cv_custom`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- init_par:

  a named list of parameter values used for initialisation of MCMC
  states. See
  [Hmsc::computePredictedValues](https://rdrr.io/pkg/Hmsc/man/computePredictedValues.html)
  for more information. Default: `NULL`.

- job_name:

  Character. Name of the submitted job(s) for SLURM. Default:
  `cv_models`.

- updater:

  named `list`. Which conditional updaters should be omitted? See
  [Hmsc::computePredictedValues](https://rdrr.io/pkg/Hmsc/man/computePredictedValues.html)
  for more information. Defaults to
  `list(Gamma2 = FALSE, GammaEta = FALSE)` to disable the following
  warnings:
  `setting updater$Gamma2=FALSE due to specified phylogeny matrix` and
  `setting updater$GammaEta=FALSE: not implemented for spatial methods 'GPP' and 'NNGP'`.

- align_posterior:

  Logical. Whether the posterior of each chains should be aligned. See
  [Hmsc::computePredictedValues](https://rdrr.io/pkg/Hmsc/man/computePredictedValues.html)
  for more information. Default: `TRUE`.

- to_json:

  Logical. Whether to convert unfitted models to JSON before saving to
  RDS file. Default: `FALSE`.

- slurm_prepare:

  Logical. Whether to prepare SLURM command files. If `TRUE` (default),
  the SLURM commands will be saved to disk using the
  [mod_slurm](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  function.

- memory_per_cpu:

  Character. Memory allocation per CPU core. Example: "32G" for 32
  gigabytes. Defaults to "64G".

- job_runtime:

  Character. Maximum allowed runtime for the job. Example: "01:00:00"
  for one hour. Required — if not provided, the function throws an
  error.

- path_hmsc:

  Character. Path to the Hmsc-HPC installation.

- precision:

  Integer. Must be either 32 or 64 (default). Defines the floating-point
  precision mode for `Hmsc-HPC` sampling (–fp 32 or –fp 64).

- ...:

  Additional arguments passed to the
  [mod_slurm](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  function.

## Details

The function copies part of the
[Hmsc::computePredictedValues](https://rdrr.io/pkg/Hmsc/man/computePredictedValues.html)
function, which currently does not support performing cross-validation
using Hmsc-HPC.

## Author

Ahmed El-Gabbas
