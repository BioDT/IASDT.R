# Draws samples from the conditional predictive distribution of latent factors

This function is optimized for speed using parallel processing and
optionally `TensorFlow` for matrix operations. This function is adapted
from
[Hmsc::predictLatentFactor](https://rdrr.io/pkg/Hmsc/man/predictLatentFactor.html)
with equivalent results to the original function when
`predictMean = TRUE`.

## Usage

``` r
predict_latent_factor(
  units_pred,
  units_model,
  post_eta,
  post_alpha,
  lf_rl,
  n_cores_lf = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  temp_dir = "temp_pred",
  lf_temp_cleanup = TRUE,
  model_name = NULL,
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  lf_out_file = NULL,
  lf_return = FALSE,
  lf_check = FALSE,
  lf_commands_only = FALSE,
  solve_max_attempts = 5L,
  solve_chunk_size = 50L,
  verbose = TRUE
)
```

## Arguments

- units_pred:

  a factor vector with random level units for which predictions are to
  be made

- units_model:

  a factor vector with random level units that are conditioned on

- post_eta:

  Character. Path of `post_eta`; a list containing samples of random
  factors at conditioned units

- post_alpha:

  a list containing samples of range (lengthscale) parameters for latent
  factors

- lf_rl:

  a HmscRandomLevel-class object that describes the random level
  structure

- n_cores_lf:

  Integer. Number of cores to use for parallel processing of latent
  factor prediction. Defaults to 8L.

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

- temp_dir:

  Character. Path for temporary storage of intermediate files.

- lf_temp_cleanup:

  Logical. Whether to delete temporary files in the `temp_dir` directory
  after finishing the LF predictions.

- model_name:

  Character. Prefix for temporary file names. Defaults to `NULL`, in
  which case no prefix is used.

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

- lf_check:

  Logical. If `TRUE`, the function checks if the output files are
  already created and valid. If `FALSE`, the function will only check if
  the files exist without checking their integrity. Default is `FALSE`.

- lf_commands_only:

  Logical. If `TRUE`, returns the command to run the Python script.
  Default is `FALSE`.

- solve_max_attempts:

  Integer. Maximum number of attempts to run solve and crossprod
  internal function
  [run_crossprod_solve](https://biodt.github.io/IASDT.R/reference/run_crossprod_solve.md).
  Default is 5L.

- solve_chunk_size:

  Integer. Chunk size for `solve_and_multiply` Python function. Default
  is 50L.

- verbose:

  Logical. If `TRUE`, logs detailed information during execution.
  Default is `TRUE`.

## Details

The function is expected to be faster than the original function in the
`Hmsc` package, especially when using `TensorFlow` for calculations and
when working in parallel.

The main difference is that this function:

- allow for parallel processing (`n_cores_lf` argument);

- when `TensorFlow` is used (`use_tf = TRUE`), matrix calculations are
  much faster, particularly when used on GPU. The following Python
  modules are needed: `numpy`, `tensorflow`, `rdata`, `xarray`, and
  `pandas`. To use `TensorFlow` under Windows, the argument `tf_environ`
  should be set to the path of a Python environment with `TensorFlow`
  installed;

- if `use_tf` is set to `FALSE`, the function uses `R` (supported by
  relatively faster `CPP` functions) in the calculations;

- `d11` and `d12` matrices are processed only once and saved to disk and
  called when needed.

## See also

[Hmsc::predictLatentFactor](https://rdrr.io/pkg/Hmsc/man/predictLatentFactor.html)
