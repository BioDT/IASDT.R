# Plot model convergence of a selected model

The `convergence_plot()` function generates and saves convergence
diagnostics plots for the `rho`, `alpha`, `omega`, and `beta` parameters
in an Hmsc model. These plots help assess whether the MCMC chains have
reached stationarity. It supports parallel processing and can work with
models fitted on HPC environments.

## Usage

``` r
convergence_plot(
  path_coda = NULL,
  env_file = ".env",
  title = " ",
  n_omega = 1000L,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 2000L,
  n_rc = list(alpha = c(2L, 3L), omega = c(2L, 2L), beta = c(3L, 3L)),
  pages_per_file = 20L,
  chain_colors = NULL,
  margin_type = "histogram",
  spatial_model = TRUE
)

convergence_alpha(
  posterior = NULL,
  title = NULL,
  n_rc_alpha = c(2L, 3L),
  add_footer = TRUE,
  add_title = TRUE,
  chain_colors = NULL,
  margin_type = "histogram",
  n_chains = NULL,
  n_samples = NULL
)

convergence_rho(
  posterior = NULL,
  title = NULL,
  chain_colors = NULL,
  margin_type = "histogram",
  n_chains = NULL,
  n_samples = NULL
)

convergence_beta_ranges(model_dir = NULL, beta_data = NULL, n_chains = NULL)
```

## Arguments

- path_coda:

  Character. Path to a saved coda object containing MCMC samples.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- title:

  Character. title for **rho** and **alpha** convergence plots. Default:
  " "

- n_omega:

  Integer. Number of species interactions sampled for Omega parameter
  diagnostics. Default: 1000L

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  8.

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

- n_rc:

  List of 3 numeric vectors representing the number of rows and columns
  for grid layout of the convergence plots of alpha, omega, and beta
  parameters. .

- pages_per_file:

  Integer. Number of plots per page in the Omega parameter output.
  Default: 20L.

- chain_colors:

  Character vector. MCMC chain colours (optional). Default: `NULL`.

- margin_type:

  Character. The type of marginal plot to add to the main plot. Valid
  options are "histogram" (default) or "density".

- spatial_model:

  Logical. Whether the model is a spatial model. If `TRUE` (default),
  the function will generate additional plots for the model's `Alpha`
  parameter.

- posterior:

  `mcmc.list` or character. Either an MCMC object (`mcmc.list`)
  containing posterior samples, or a file path to a saved coda object.

- n_rc_alpha:

  Numeric vector of length 2. Number of rows and columns for the
  convergence plots of the `alpha` parameter. Default: `c(2L, 3L)`.

- add_footer:

  Logical. If `TRUE` (default), adds a footer with page numbers to each
  plot.

- add_title:

  Logical. If `TRUE` (default), adds the main title (`title`) to the
  plot.

- n_chains:

  Integer. Number of MCMC chains.

- n_samples:

  Integer. Number of MCMC samples.

- model_dir:

  Character. A path to the model directory.

- beta_data:

  Data frame. Beta parameter summary data frame.

## Details

`convergence_alpha()`, `convergence_rho()`, and
`convergence_beta_ranges` are internal functions and should not be
called directly. The `convergence_beta_ranges` plots the convergence
range of the each species beta parameters. It can be used to check if
any of the chains show convergence issues; i.e., showing exceptionally
high or low beta values.

## Author

Ahmed El-Gabbas
