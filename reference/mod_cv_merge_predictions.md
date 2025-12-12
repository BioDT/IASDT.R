# Merge Cross-Validation Predictions and Generate Summary Maps

Merges prediction results from cross-validated Hmsc models, processes
clamped and not-clamped model outputs, and generates summary maps and
anomaly visualizations for species invasion predictions. It saves
processed prediction data and map images to disk.

## Usage

``` r
mod_cv_merge_predictions(
  model_prefix = NULL,
  n_cv_folds = 4L,
  hab_abb = NULL,
  spatial_model = NULL,
  n_cores = 8L,
  clamp = c("clamp", "no_clamp")
)
```

## Arguments

- model_prefix:

  Character. Prefix for model output directories and files.

- n_cv_folds:

  Integer. Number of cross-validation folds. Default is 4L.

- hab_abb:

  Character. Habitat abbreviation to process.

- spatial_model:

  Logical. Whether to use a spatial model (TRUE) or non-spatial (FALSE).

- n_cores:

  Integer. Number of cores for parallel processing. Default is 8L.

- clamp:

  Character vector. Indicates which prediction types to process:
  "clamp", "no_clamp", or both.

## Author

Ahmed El-Gabbas
