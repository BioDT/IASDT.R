# Prepare and plot response curve data for Hmsc models

The `rc_*()` functions process and visualise response curves for Hmsc
models. They support parallel computation and optionally return
processed data. There are four functions in this group:

- `rc_prepare_data()`: Prepares response curve data for analysis

- `rc_plot_species()`: Generates response curve plots for individual
  species

- `rc_plot_species_all()`: Generates response curves for all species
  together in a single plot

- `rc_plot_sr()`: Plots response curves for species richness.

## Usage

``` r
rc_prepare_data(
  path_model = NULL,
  n_grid = 50L,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1500L,
  return_data = FALSE,
  probabilities = c(0.025, 0.5, 0.975),
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  n_cores_lf = n_cores,
  lf_check = FALSE,
  lf_temp_cleanup = TRUE,
  lf_commands_only = FALSE,
  temp_dir = "temp_pred",
  temp_cleanup = TRUE,
  verbose = TRUE
)

rc_plot_species(
  model_dir = NULL,
  n_cores = 20,
  env_file = ".env",
  return_data = FALSE
)

rc_plot_species_all(
  model_dir = NULL,
  n_cores = 8L,
  return_data = FALSE,
  plotting_alpha = 0.3
)

rc_plot_sr(
  model_dir,
  verbose = TRUE,
  n_cores = 8L,
  future_max_size = 1000L,
  strategy = "multisession"
)
```

## Arguments

- path_model:

  Character. Path to the file containing the fitted Hmsc model.

- n_grid:

  Integer. Number of points along the gradient for continuous focal
  variables. Higher values result in smoother curves. Default: 50. See
  [Hmsc::constructGradient](https://rdrr.io/pkg/Hmsc/man/constructGradient.html)
  for details.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Defaults
  to 8L for all functions, except for `rc_plot_species`, in which it
  defaults to 20L.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

- future_max_size:

  Numeric. Maximum allowed total size (in megabytes) of global variables
  identified. See `future.globals.maxSize` argument of
  [future::future.options](https://future.futureverse.org/reference/zzz-future.options.html)
  for more details.

- return_data:

  Logical. If `TRUE`, the function returns processed data as an R
  object. Default: `FALSE`.

- probabilities:

  Numeric vector. Quantiles to calculate in response curve predictions.
  Default: `c(0.025, 0.5, 0.975)`. See
  [stats::quantile](https://rdrr.io/r/stats/quantile.html) for details.

- use_tf:

  Logical. Whether to use `TensorFlow` for calculations. Defaults to
  `TRUE`.

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- n_cores_lf:

  Integer. Number of cores to use for parallel processing of latent
  factor prediction. Defaults to 8L.

- lf_check:

  Logical. If `TRUE`, the function checks if the output files are
  already created and valid. If `FALSE`, the function will only check if
  the files exist without checking their integrity. Default is `FALSE`.

- lf_temp_cleanup:

  Logical. Whether to delete temporary files in the `temp_dir` directory
  after finishing the LF predictions.

- lf_commands_only:

  Logical. If `TRUE`, returns the command to run the Python script.
  Default is `FALSE`.

- temp_dir:

  Character. Path for temporary storage of intermediate files.

- temp_cleanup:

  Logical. Whether to clean up temporary files. Defaults to `TRUE`.

- verbose:

  Logical. Whether to print a message upon successful saving of files.
  Defaults to `FALSE`.

- model_dir:

  Character. Path to the root directory containing fitted models. The
  function reads data from the `response_curves_data` subdirectory,
  which is created by `rc_prepare_data`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- plotting_alpha:

  Numeric. Opacity level for response curve lines (0 = fully
  transparent, 1 = fully opaque). Default: 0.3.

## Author

Ahmed El-Gabbas
