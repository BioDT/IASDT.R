# Plot model convergence of multiple modelling alternatives

This function generates and saves a series of diagnostic plots to assess
the convergence of Hmsc models across multiple modelling alternatives.
It checks model convergence using trace plots and Gelman-Rubin
diagnostics for key model parameters.

## Usage

``` r
convergence_plot_all(
  model_dir = NULL,
  n_omega = 1000L,
  margin_type = "histogram",
  spatial_model = TRUE,
  n_cores = NULL,
  strategy = "multisession",
  future_max_size = 1000L,
  n_rc_alpha = c(2L, 3L)
)
```

## Arguments

- model_dir:

  Character. Path to the root directory of the fitted model. The
  convergence outputs will be saved to the `model_convergence_all`
  subdirectory.

- n_omega:

  Integer. Number of species interactions sampled for Omega parameter
  diagnostics. Default: 1000L

- margin_type:

  Character. The type of marginal plot to add to the main plot. Valid
  options are "histogram" (default) or "density".

- spatial_model:

  Logical. Whether the model is a spatial model. If `TRUE` (default),
  the function will generate additional plots for the model's `Alpha`
  parameter.

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

- future_max_size:

  Numeric. Maximum allowed total size (in megabytes) of global variables
  identified. See `future.globals.maxSize` argument of
  [future::future.options](https://future.futureverse.org/reference/zzz-future.options.html)
  for more details.

- n_rc_alpha:

  Numeric vector of length 2. Number of rows and columns for the
  convergence plots of the `alpha` parameter. Default: `c(2L, 3L)`.

## Value

The function returns `invisible(NULL)` and does not return any value,
but saves a series of diagnostic plots in the specified path.

## Author

Ahmed El-Gabbas
