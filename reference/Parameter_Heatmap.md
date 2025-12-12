# Heatmaps for the `beta` and `omega` parameters of the Hmsc model

The `mod_heatmap_beta()` and `mod_heatmap_omega()` functions generate
heatmaps using `ggplot2` to visualise parameter estimates or posterior
support values for species' environmental responses (`beta` parameters,
which describes how species (*Y*) respond to various covariates (*X*);
see [Hmsc::plotBeta](https://rdrr.io/pkg/Hmsc/man/plotBeta.html)) and
residual associations (`omega` parameter), respectively.

## Usage

``` r
mod_heatmap_beta(
  path_model = NULL,
  support_level = 0.95,
  width = 25,
  height = 35
)

mod_heatmap_omega(
  path_model = NULL,
  support_level = 0.95,
  width = 26,
  height = 22.5
)
```

## Arguments

- path_model:

  Character. Path to the fitted `Hmsc` model object.

- support_level:

  Numeric. The posterior support threshold for determining which values
  are considered significant in the heatmap. Defaults to 0.95,
  indicating 95% posterior support. Values above this threshold (or
  below 1 - threshold for negative associations) are considered
  significant and will be plotted (see
  [Hmsc::plotBeta](https://rdrr.io/pkg/Hmsc/man/plotBeta.html)).

- width, height:

  Integer. The width and height of the generated heatmaps in
  centimetres. Defaults to 26×22.5 for `omega`; 25×35 for `beta`.

## Value

Both functions do not return a value but saves heatmap plots as JPEG
files in the `model_postprocessing/parameters_summary` subdirectory.

## Details

The functions exports three types of visualisations (see
[Hmsc::plotBeta](https://rdrr.io/pkg/Hmsc/man/plotBeta.html)):

- `mean`: posterior mean estimate,

- `support`: statistical support level, measured by the posterior
  probability for a positive or negative response,

- `sign`: indicates whether the response is positive, negative, or
  neither of these based on the chosen `support_level`.

For the `omega` parameter, the `mod_heatmap_omega()` function generates
two JPEG files: signs and mean values. While for the `beta` parameter,
the `mod_heatmap_beta()` function generates four JPEG files : support,
signs, mean values (including and excluding the intercept).

## Author

Ahmed El-Gabbas. The `mod_heatmap_beta()` function is adapted from
[Hmsc::plotBeta](https://rdrr.io/pkg/Hmsc/man/plotBeta.html)
