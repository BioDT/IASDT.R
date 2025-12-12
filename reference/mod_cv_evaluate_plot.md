# Plot Evaluation Results for Cross-Validated Hmsc Models

This function evaluates cross-validation results for Hmsc models
(spatial and/or non-spatial), summarizes performance metrics (RMSE,
TjurR2, AUC, and Boyce), and generates diagnostic plots comparing model
types.

## Usage

``` r
mod_cv_evaluate_plot(
  model_prefix = NULL,
  hab_abb = NULL,
  n_cv_folds = 4,
  spatial_model = c("gpp", "nonspatial")
)
```

## Arguments

- model_prefix:

  Character. Prefix for model directory name.

- hab_abb:

  Character. Habitat abbreviation.

- n_cv_folds:

  Integer. Number of cross-validation folds (default: 4L).

- spatial_model:

  Character vector. Specifies which models to evaluate: "gpp",
  "nonspatial", or both.

## Author

Ahmed El-Gabbas
