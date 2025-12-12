# Computes and visualise variance partitioning of Hmsc models

Computes and plots variance components with respect to given grouping of
fixed effects and levels of random effects. The The
`variance_partitioning_compute()` function inherits the main
functionality from the
[Hmsc::computeVariancePartitioning](https://rdrr.io/pkg/Hmsc/man/computeVariancePartitioning.html)
function, but with the added functionality of parallel computation and
using `TensorFlow`. The `variance_partitioning_plot()` function
generates plots for variance partitioning as JPEG files. It allows for
sorting the predictors and species; e.g., by the mean value per
predictor; and by original species order. It also plots the raw variance
partitioning (relative variance partitioning multiplied by the training
and testing (if supported) Tjur-R² value).

## Usage

``` r
variance_partitioning_compute(
  path_model,
  group = NULL,
  group_names = NULL,
  start = 1L,
  na.ignore = FALSE,
  n_cores = 8L,
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  temp_cleanup = TRUE,
  chunk_size = 50L,
  verbose = TRUE,
  vp_file = "varpar",
  vp_commands_only = FALSE,
  temp_dir = NULL
)

variance_partitioning_plot(
  path_model = NULL,
  env_file = ".env",
  vp_file = "varpar",
  use_tf = TRUE,
  tf_environ = NULL,
  n_cores = 1L,
  width = 30,
  height = 15,
  axis_text = 4,
  spatial_model = TRUE,
  is_cv_model = FALSE,
  temp_dir = NULL
)
```

## Arguments

- path_model:

  Character. Path to fitted `Hmsc` model object.

- group:

  vector of numeric values corresponding to group identifiers in
  groupnames. If the model was defined with `XData` and `XFormula`, the
  default is to use model terms.

- group_names:

  vector of names for each group of fixed effect. Should match `group`.
  If the model was defined with `XData` and `XFormula`, the default is
  to use the labels of model terms.

- start:

  index of first MCMC sample included. Default: `1L`.

- na.ignore:

  Logical. If `TRUE`, covariates are ignored for sites where the focal
  species is NA when computing variance-covariance matrices for each
  species.

- n_cores:

  Integer. Number of CPU cores to use for computing variance
  partitioning using `TensorFlow`. This is only effective when `use_tf`
  is `TRUE`. Default: `1`.

- use_tf:

  Logical. Whether to use `TensorFlow` for calculations. Defaults to
  `TRUE`.

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- temp_cleanup:

  Logical. Whether to delete temporary files after processing. Default:
  `TRUE`.

- chunk_size:

  Integer. Size of each chunk of samples to process in parallel. Only
  relevant for `TensorFlow`. Default: `50`.

- verbose:

  Logical. Whether to print progress messages. Default: `TRUE`.

- vp_file:

  Character. Name of the output file to save the results. Default:
  `varpar`.

- vp_commands_only:

  Logical. If `TRUE`, returns the commands to run the Python script.
  Default is `FALSE`. Only relevant when `use_tf` is `TRUE`.

- temp_dir:

  Character. Path to a temporary directory to store intermediate files.
  Default: `NULL`, which creates a temporary directory in the same
  parent directory as the model file.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- width, height:

  Numeric. Width and height of the output plot in centimetres. Default:
  `30` and `15`, respectively.

- axis_text:

  Numeric. Size of the axis text. Default: `4`.

- spatial_model:

  Logical. Whether the fitted model is a spatial model. Defaults to
  `TRUE`.

- is_cv_model:

  Logical. Whether the model is a cross-validated model (`TRUE`) or
  fitted with the full dataset (`FALSE`; default). If `TRUE`, the
  explanatory and predictive power of the model will be used to estimate
  the raw variance partitioning.

## Author

Ahmed El-Gabbas
