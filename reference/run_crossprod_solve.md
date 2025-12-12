# run_crossprod_solve

Internal function to executes a Python script that performs matrix
computations using `TensorFlow` with provided inputs. Retries up to
three times if the output file validation fails.

## Usage

``` r
run_crossprod_solve(
  tf_environ,
  s1,
  s2,
  post_eta,
  path_out,
  denom,
  chunk_size = 1000L,
  threshold_mb = 2000L,
  tf_use_single = TRUE,
  verbose = TRUE,
  solve_chunk_size = 50L,
  solve_max_attempts = 5L,
  lf_commands_only = FALSE
)
```

## Arguments

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- s1:

  Character. Path to the input file containing s1 coordinates.

- s2:

  Character Path to the input file containing s2 coordinates.

- post_eta:

  Character. Path to the file containing the `post_eta` matrix data.

- path_out:

  Character. Path to rds file where the output results will be saved.

- denom:

  Numeric. The denominator value used in the computation.

- chunk_size:

  Numeric (Optional). Size of chunks to process at a time. Default is
  1000.

- threshold_mb:

  Numeric (Optional). Memory threshold (in MB) to manage processing.
  Default is 2000.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- verbose:

  Logical. If `TRUE`, logs detailed information during execution.
  Default is `TRUE`.

- solve_chunk_size:

  Integer. Chunk size for `solve_and_multiply` Python function. Default
  is 50L.

- solve_max_attempts:

  Integer. Maximum number of attempts to run solve and crossprod
  internal function run_crossprod_solve. Default is 5L.

- lf_commands_only:

  Logical. If `TRUE`, returns the command to run the Python script.
  Default is `FALSE`.

## Value

Returns the `path_out` if successful. Returns `NULL` if all attempts
fail.

## Details

- The function checks for the existence of required input files and the
  Python executable in the specified virtual environment.

- Executes the Python script using `system2`.

- Verifies the output file validity using
  [`ecokit::check_data`](https://elgabbas.github.io/ecokit/reference/check_data.html).
  Retries up to 3 times if the output is invalid.

- Generates detailed logs if `verbose` is set to `TRUE`.

## Author

Ahmed El-Gabbas
