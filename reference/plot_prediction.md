# Plot species and level of invasion predictions as JPEG files using `ggplot2`

Generate predictions for species and habitat models and saves the output
as JPEG files.

## Usage

``` r
plot_prediction(
  model_dir = NULL,
  env_file = ".env",
  n_cores = 8L,
  is_cv_model = FALSE
)
```

## Arguments

- model_dir:

  Character. Path to the model directory containing predictions.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  8.

- is_cv_model:

  Logical. Whether the model is a cross-validated model (`TRUE`) or
  fitted with the full dataset (`FALSE`; default).

## Value

Saves prediction plots as JPEG files in the specified output directory.

## Author

Ahmed El-Gabbas
