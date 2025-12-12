# Summary of Hmsc model parameters

This function provides a comprehensive summary of Hmsc model parameters,
including `Alpha`, `Beta`, `Rho`, and `Omega`. It processes the model's
output, performs statistical summaries, and optionally returns the
summarised data.

## Usage

``` r
mod_summary(
  path_coda = NULL,
  env_file = ".env",
  return_data = FALSE,
  spatial_model = TRUE
)
```

## Arguments

- path_coda:

  Character. Path to the `.qs2` / `.RData` file containing the coda
  object.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- return_data:

  Logical. Whether the summarised data should be returned as an R
  object. If `TRUE`, the function returns a list containing summaries of
  `Alpha`, `Beta`, `Rho`, and `Omega` parameters. The default value is
  `FALSE`, which means the function will not return any data but will
  save the summaries to a specified directory.

- spatial_model:

  Logical. Whether the model is a spatial model. If `TRUE` (default),
  the function will also process the `Alpha` parameter.

## Value

If `return_data` is `FALSE` (default), the function does not return
anything and saves the summaries to a directory. If `return_data` is
`TRUE`, it also returns the data as R object.

## Author

Ahmed El-Gabbas
