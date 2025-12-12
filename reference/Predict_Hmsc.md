# Calculates predicted values from a fitted Hmsc model

This function modifies the `Hmsc:::predict.Hmsc` function.

## Usage

``` r
predict_hmsc(
  path_model,
  Loff = NULL,
  x_data = NULL,
  X = NULL,
  XRRRData = NULL,
  XRRR = NULL,
  gradient = NULL,
  Yc = NULL,
  mcmcStep = 1L,
  expected = TRUE,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  model_name = "train",
  temp_dir = "temp_pred",
  temp_cleanup = TRUE,
  prediction_type = NULL,
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  lf_out_file = NULL,
  lf_return = FALSE,
  lf_input_file = NULL,
  lf_only = FALSE,
  n_cores_lf = n_cores,
  lf_check = FALSE,
  lf_temp_cleanup = TRUE,
  lf_commands_only = FALSE,
  pred_directory = NULL,
  pred_pa = NULL,
  pred_xy = NULL,
  evaluate = FALSE,
  evaluation_name = NULL,
  evaluation_directory = "evaluation",
  verbose = TRUE,
  spatial_model = TRUE
)
```

## Arguments

- path_model:

  Character. Path to the saved model object.

- Loff:

  See
  [Hmsc::predict.Hmsc](https://rdrr.io/pkg/Hmsc/man/predict.Hmsc.html)
  for more details.

- x_data:

  `data.frame`. The unpreprocessed covariates for the predictions to be
  made. Works only if the `XFormula` argument was specified in the
  [Hmsc::Hmsc](https://rdrr.io/pkg/Hmsc/man/Hmsc.html) model constructor
  call. Requirements are similar to those in the `Hmsc` model
  constructor.

- X:

  `matrix`. Covariates for the predictions to be made. Only one of
  `x_data` and `X` arguments may be provided.

- XRRRData:

  `data.frame`. Covariates for reduced-rank regression.

- XRRR:

  `matrix`. Covariates for reduced-rank regression.

- gradient:

  An object returned by
  [Hmsc::constructGradient](https://rdrr.io/pkg/Hmsc/man/constructGradient.html).
  Providing `gradient` is an alternative for providing `x_data`,
  `studyDesign` and `ranLevels`. Cannot be used together with `Yc`.

- Yc:

  `matrix`. Outcomes that are assumed to be known for conditional
  predictions. Cannot be used together with `gradient`.

- mcmcStep:

  Integer. Number of extra `mcmc` steps used for updating the random
  effects. Defaults to 1L.

- expected:

  Logical. Whether to return the location parameter of the observation
  models or sample the values from those. Defaults to `TRUE`.

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

- model_name:

  Character. Prefix for temporary file names. Defaults to `NULL`, in
  which case no prefix is used.

- temp_dir:

  Character. Path for temporary storage of intermediate files.

- temp_cleanup:

  Logical. Whether to clean up temporary files. Defaults to `TRUE`.

- prediction_type:

  Character. Type of predictions to be made. If `NULL` (default),
  predictions are made for the latent factors. If `c`, predictions are
  made for response curves at mean coordinates. If `i`, predictions are
  made for response curves at infinite coordinates.

- use_tf:

  Logical. Whether to use `TensorFlow` for calculations. Defaults to
  `TRUE`.

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- lf_out_file:

  Character. Path to save the outputs. If `NULL` (default), the
  predicted latent factors are not saved to a file. This should end with
  either `*.qs2` or `*.RData`.

- lf_return:

  Logical. Whether the output should be returned. Defaults to `FALSE`.
  If `lf_out_file` is `NULL`, this parameter cannot be set to `FALSE`
  because the function needs to return the result if it is not saved to
  a file.

- lf_input_file:

  Character. File name where the latent factor predictions are saved. If
  `NULL` (default), latent factor predictions will be computed. If
  specified, latent factor predictions are read from this path. This
  allows to predicting the latent factors for new sites only once.

- lf_only:

  Logical. Whether to return the latent factor predictions only.
  Defaults to `FALSE`. This helps in predicting to new sites, allowing
  to predicting the latent factors only once, then the output can be
  loaded in other predictions when needed.

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

- lf_commands_only:

  Logical. If `TRUE`, returns the command to run the Python script.
  Default is `FALSE`.

- pred_directory:

  Character. Directory path indicating where the predictions will be
  saved. Defaults to `NULL`, which saves model predictions to
  "`model_prediction`" folder of the current working directory.

- pred_pa:

  `matrix`. Presence-absence data for evaluation. If `NULL` (default),
  the presence-absence data from the model object is used. This argument
  is used only when `evaluate` is `TRUE`.

- pred_xy:

  `matrix`. Coordinates to be added to predicted values. If `NULL`
  (default), the coordinates from the model object is used.

- evaluate:

  Logical. Whether to evaluate the model predictions. Defaults to
  `FALSE`.

- evaluation_name:

  Character. Name of the evaluation results. If `NULL`, the default name
  is used (`eval_[model_name].qs2`).

- evaluation_directory:

  Character. Directory where the evaluation results will be saved.
  Defaults to `evaluation`.

- verbose:

  Logical. Whether to print a message upon successful saving of files.
  Defaults to `FALSE`.

- spatial_model:

  Logical. Whether the fitted model is a spatial model. Defaults to
  `TRUE`.
