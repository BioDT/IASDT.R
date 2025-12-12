# Cross-validation Model Evaluation and Plotting

performs model evaluation using cross-validation results, calculates
multiple metrics (AUC, Tjur R², Boyce index, RMSE), and generates
summary plots for explanatory and predictive power.

## Usage

``` r
mod_cv_evaluate(model_dir = NULL, cv_type = "cv_dist")
```

## Arguments

- model_dir:

  Character. Path to the root directory of the fitted model.

- cv_type:

  Character. Cross-validation type. One of `cv_dist` (default) or
  `cv_large`. See
  [`mod_cv_fit()`](https://biodt.github.io/IASDT.R/reference/Mod_CV_Fit.md)
  for more details.

## Value

Invisibly returns the path to the saved evaluation data file.

## Author

Ahmed El-Gabbas
