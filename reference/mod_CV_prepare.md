# Prepare spatial-block cross-validation folds for spatial analysis

This function assign modelling input data into spatial-block
cross-validation folds using three strategies (see below) using
[blockCV::cv_spatial](https://rdrr.io/pkg/blockCV/man/cv_spatial.html).
The function is planned to be used inside the
[mod_prepare_hpc](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
function.

## Usage

``` r
mod_cv_prepare(
  input_data = NULL,
  env_file = ".env",
  x_vars = NULL,
  cv_n_folds = 4L,
  cv_n_grids = 20L,
  cv_n_rows = 2L,
  cv_n_columns = 2L,
  cv_sac = FALSE,
  out_path = NULL
)
```

## Arguments

- input_data:

  `data.frame`. A data frame or tibble containing the input dataset.
  This data frame should include two columns for `x` and `y` coordinates
  as long as other columns matching the names of predictors listed in
  `x_vars` argument. This argument is mandatory and can not be empty.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- x_vars:

  Character vector. Variables to be used in the model. This argument is
  mandatory and can not be empty.

- cv_n_folds:

  Integer. Number of cross-validation folds. Default: 4L.

- cv_n_grids:

  Integer. Number of grid cells in both directions used in the `cv_dist`
  cross-validation strategy (see below). Default: 20L.

- cv_n_rows, cv_n_columns:

  Integer. Number of rows and columns used in the `cv_large`
  cross-validation strategy (see below), in which the study area is
  divided into large blocks given the provided `cv_n_rows` and
  `cv_n_columns` values. Both default to 2L which means to split the
  study area into four large blocks at the median latitude and
  longitude.

- cv_sac:

  Logical. Whether to use the spatial autocorrelation to determine the
  block size. Defaults to `FALSE`,

- out_path:

  Character. Path for directory to save the cross-validation results.
  This argument is mandatory and can not be empty.

## Value

The function returns a modified version of the input dataset with
additional numeric columns (integer) indicating the cross-validation
strategy used.

## Note

The function uses the following cross-validation strategies:

- `cv_dist` in which the size of spatial cross-validation blocks is
  determined by the `cv_n_grids` argument. The default `cv_n_grids`
  value is 20L, which means blocks of 20×20 grid cell each.

- `cv_large` which splits the study area into large blocks, as
  determined by the `cv_n_rows` and `cv_n_columns` arguments. if
  `cv_n_rows = cv_n_columns` = 2L (default), four large blocks will be
  used, split the study area at the median coordinates.

- `cv_sac` in which the size of the blocks is determined by the median
  spatial autocorrelation range in the predictor data (estimated using
  [blockCV::cv_spatial_autocor](https://rdrr.io/pkg/blockCV/man/cv_spatial_autocor.html)).
  This requires the availability of the `automap` R package. This
  strategy is currently skipped by default.

## Author

Ahmed El-Gabbas
