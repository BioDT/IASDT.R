# Plot Gelman-Rubin-Brooks

The `plot_gelman_*()` functions generate plots visualising the evolution
of the Gelman-Rubin-Brooks shrink factor for different model parameters
as the number of iterations increases. These plots help assess whether
MCMC chains have converged to a common distribution. Each plot includes:
median (solid line) and 97.5^(th) percentile (dashed line) of the shrink
factor and a dashed horizontal line at 1.1, representing the common
convergence threshold. The primary function for users is
`plot_gelman()`, which internally calls:

- `plot_gelman_alpha()`: Plots shrink factor for the **Alpha** parameter

- `plot_gelman_beta()`: Plots shrink factor for the **Beta** parameters

- `plot_gelman_omega()`: Plots shrink factor for the **Omega** parameter

- `plot_gelman_rho()`: Plots shrink factor for the **Rho** parameter

## Usage

``` r
plot_gelman(
  path_coda = NULL,
  alpha = TRUE,
  beta = TRUE,
  omega = TRUE,
  rho = TRUE,
  n_omega = 1000L,
  plotting_alpha = 0.25,
  env_file = ".env"
)

plot_gelman_alpha(coda_object, plotting_alpha = 0.25)

plot_gelman_beta(coda_object, env_file = ".env", plotting_alpha = 0.25)

plot_gelman_omega(coda_object, n_omega = 1000L, plotting_alpha = 0.25)

plot_gelman_rho(coda_object)
```

## Arguments

- path_coda:

  Character. Path to a file containing the coda object, representing
  MCMC samples.

- alpha, beta, omega, rho:

  Logical. If `TRUE`, plots the Gelman-Rubin statistic for the
  respective model parameters (alpha, beta, omega, or rho). Default:
  `TRUE` for all parameters.

- n_omega:

  Integer. Number of species sampled for the omega parameter. Default:
  1000L.

- plotting_alpha:

  Numeric. Transparency level (alpha) for plot lines (0 = fully
  transparent, 1 = fully opaque). Default: 0.25.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- coda_object:

  `mcmc.list`. An MCMC sample object containing posterior distributions
  from an Hmsc model.

## Author

Ahmed El-Gabbas
