# Fit Hmsc-HPC models on UFZ Windows Server

This function fits Hmsc models on a UFZ Windows Server. It reads model
configurations from a specified path, loads environment variables,
checks input arguments for validity, and executes model fitting in
parallel if required.

## Usage

``` r
mod_fit_windows(
  path_model = NULL,
  python_ve = NULL,
  n_cores = NULL,
  strategy = "multisession"
)
```

## Arguments

- path_model:

  Character. Path to the model files. This argument can not be empty.

- python_ve:

  Character. Path to a valid Python virtual environment. Defaults to
  `NULL`. This argument can not be empty.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

## Value

The function does not return anything but prints messages to the console
regarding the progress and completion of model fitting.

## Author

Ahmed El-Gabbas
