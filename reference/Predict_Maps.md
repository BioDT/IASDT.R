# Predict habitat suitability of `Hmsc` models

This package provides two functions for predicting habitat suitability
of `Hmsc` models in the `IASDT` framework. predict_maps generates
current and future habitat suitability maps (mean, sd, cov) from a full
Hmsc model fit. predict_maps_cv predicts and evaluates cross-validated
Hmsc models for current climate conditions. For more details, see the
respective function documentation and the details section below.

## Usage

``` r
predict_maps(
  path_model = NULL,
  hab_abb = NULL,
  env_file = ".env",
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  clamp_pred = TRUE,
  fix_efforts = "q90",
  fix_rivers = "q90",
  pred_new_sites = TRUE,
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  n_cores_lf = n_cores,
  lf_check = FALSE,
  lf_temp_cleanup = TRUE,
  lf_only = FALSE,
  lf_commands_only = FALSE,
  temp_dir = "temp_pred",
  temp_cleanup = TRUE,
  tar_predictions = TRUE,
  spatial_model = TRUE,
  n_cores_pred = n_cores,
  climate_models = c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0",
    "UKESM1-0-LL"),
  climate_scenario = c("ssp126", "ssp370", "ssp585")
)

predict_maps_cv(
  model_dir = NULL,
  cv_name = NULL,
  cv_fold = NULL,
  n_cores = 8L,
  strategy = "multisession",
  env_file = ".env",
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  n_cores_lf = n_cores,
  lf_check = FALSE,
  lf_temp_cleanup = TRUE,
  lf_only = FALSE,
  lf_commands_only = FALSE,
  temp_cleanup = TRUE
)
```

## Arguments

- path_model:

  Character. Path to fitted `Hmsc` model object.

- hab_abb:

  Character. Habitat abbreviation indicating the specific
  [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type.
  Valid values: `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`. See
  [Pysek et al.](https://doi.org/10.23855/preslia.2022.447) for details.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

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

- clamp_pred:

  Logical indicating whether to clamp the sampling efforts at a single
  value. If `TRUE` (default), the `fix_efforts` argument must be
  provided.

- fix_efforts:

  Numeric or character. When `clamp_pred = TRUE`, fixes the sampling
  efforts predictor at this value during predictions. If numeric, uses
  the value directly (on log₁₀ scale). If character, must be one of
  `identity` (i.e., do not fix), `median`, `mean`, `max`, or `q90` (90%
  quantile). Using `max` may reflect extreme sampling efforts from
  highly sampled locations, while `q90` captures high sampling areas
  without extremes. Required if `clamp_pred = TRUE`.

- fix_rivers:

  Numeric, character, or `NULL`. Similar to `fix_efforts`, but for the
  river length predictor. If `NULL`, the river length is not fixed.
  Default: `q90`.

- pred_new_sites:

  Logical. Whether to predict suitability at new sites. Default: `TRUE`.

- use_tf:

  Logical. Whether to use `TensorFlow` for calculations. Defaults to
  `TRUE`.

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- n_cores_lf:

  Integer. Number of cores to use for parallel processing of latent
  factor prediction. Defaults to 8L.

- lf_check:

  Logical. If `TRUE`, the function checks if the output files are
  already created and valid. If `FALSE`, the function will only check if
  the files exist without checking their integrity. Default is `FALSE`.

- lf_temp_cleanup:

  Logical. Whether to delete temporary files in the `temp_dir` directory
  after finishing the LF predictions.

- lf_only:

  Logical. Whether to predict only the latent factor. This is useful for
  distributing processing load between GPU and CPU. When
  `lf_only = TRUE`, latent factor prediction needs to be computed
  separately on GPU. When computations are finished on GPU, the function
  can later be rerun with `lf_only = FALSE` (default) to predict habitat
  suitability using the already-computed latent factor predictions.

- lf_commands_only:

  Logical. If `TRUE`, returns the command to run the Python script.
  Default is `FALSE`.

- temp_dir:

  Character. Path for temporary storage of intermediate files.

- temp_cleanup:

  Logical. Whether to clean up temporary files. Defaults to `TRUE`.

- tar_predictions:

  Logical. Whether to compress tiff files for predicted habitat
  suitability into a single `*.tar` file (without compression). Default:
  `TRUE`.

- spatial_model:

  Logical. Whether the fitted model is a spatial model. Defaults to
  `TRUE`.

- n_cores_pred:

  Integer. Number of cores to use for predicting species' habitat
  suitability.

- climate_models:

  Character vector. Climate models for future predictions. Available
  options are
  `c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")`
  (default).

- climate_scenario:

  Character vector. Climate scenarios for future predictions. Available
  options are: `c("ssp126", "ssp370", "ssp585")` (default).

- model_dir:

  Character. Path to the directory containing cross-validated models.

- cv_name:

  Character. Cross-validation strategy. Valid values are `cv_dist`,
  `cv_large`, or `cv_sac`.

- cv_fold:

  Integer. The cross-validation fold number.

## Details

- **`predict_maps`**: Generates habitat suitability maps for `Hmsc`
  models fitted on the full dataset, for both current and future climate
  options. It produces maps for mean, standard deviation (sd), and
  coefficient of variation (cov) of suitability for each species and
  overall species richness. It evaluate model's explanatory power using
  various metrics. For future predictions, it also generates anomaly
  maps (future - current). The function supports ensemble predictions
  across multiple climate models and prepares data for upload to the
  [OPeNDAP](http://opendap.biodt.eu/ias-pdt/) server for use in the
  `IASDT` [Shiny App](https://app.biodt.eu/).

- **`predict_maps_cv`**: Computes predictions for cross-validated `Hmsc`
  models using only the testing folds. It evaluates model performance
  (predictive power) with various metrics and plots evaluation results
  for predictive and explanatory power. Unlike `predict_maps`, this
  function does not perform clamping and does not generate future
  climate predictions.

## See also

predict_hmsc

## Author

Ahmed El-Gabbas
