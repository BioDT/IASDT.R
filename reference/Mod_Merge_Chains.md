# Merge model chains into `Hmsc` and `coda` objects

These functions merge posterior chains from multiple runs of `Hmsc`
models into unified `Hmsc` and `coda` objects, facilitating further
analysis. They check for missing or incomplete chains, optionally report
these issues, and save the processed results to disk. `mod_merge_chains`
handles regular models, while `mod_merge_chains_cv` is designed for
cross-validated models.

## Usage

``` r
mod_merge_chains(
  model_dir = NULL,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  model_info_name = NULL,
  print_incomplete = TRUE,
  from_json = FALSE,
  out_extension = "qs2"
)

mod_merge_chains_cv(
  model_dir = NULL,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  cv_names = c("cv_dist", "cv_large"),
  from_json = FALSE,
  out_extension = "qs2"
)
```

## Arguments

- model_dir:

  Character. Path to the root directory where the model was fitted. For
  `mod_merge_chains`, subdirectories `model_fitted` and `model_coda` are
  created within this path to store the merged `Hmsc` and `coda`
  objects, respectively. For `mod_merge_chains_cv`, merged objects are
  stored under `model_fitting_cv/model_fitted`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Defaults
  to 8L.

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

- model_info_name:

  Character. Name of the file (without extension) where updated model
  information is saved. If `NULL`, overwrites the existing
  `model_info.RData` file in `model_dir` directory. If specified,
  creates a new `.RData` file with this name in `model_dir` directory.
  Applies only to `mod_merge_chains`.

- print_incomplete:

  Logical. If `TRUE`, prints the names of model variants that failed to
  merge due to missing or incomplete chains. Defaults to `TRUE`.

- from_json:

  Logical. Whether to convert loaded models from JSON format before
  reading. Defaults to `FALSE`.

- out_extension:

  Character. File extension (without dot) for output files containing
  merged `Hmsc` and `coda` objects. Options are `qs2` (faster read/write
  via the `qs2` package) or `RData` (standard R format). Defaults to
  `qs2`.

- cv_names:

  Character vector. Names of cross-validation strategies to merge,
  matching those used during model setup. Defaults to
  `c("cv_dist", "cv_large")`. The names should be one of `cv_dist`,
  `cv_large`, or `cv_sac`. Applies only to `mod_merge_chains_cv`.

## Value

Both functions return `invisible(NULL)` and save processed model
information and merged objects to disk in the specified locations.

## Details

- `mod_merge_chains` merges posterior chains from multiple runs of an
  `Hmsc` model fitted without cross-validation. It checks for missing or
  incomplete chains, aligns posteriors (using `alignPost = TRUE`,
  falling back to `FALSE` if alignment fails), and saves a merged `Hmsc`
  object and a `coda` object for MCMC diagnostics. It also records
  fitting time and memory usage from progress files.

- `mod_merge_chains_cv` performs a similar merging process for
  cross-validated `Hmsc` models, processing each fold of the specified
  `cv_names` separately. It saves merged `Hmsc` objects per fold but
  does not generate `coda` objects.

## Author

Ahmed El-Gabbas
