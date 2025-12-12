# Process CHELSA Climate Data for the `IASDT`

Downloads, processes, and projects [CHELSA](https://chelsa-climate.org/)
climate data at the European scale for the Invasive Alien Species
Digital Twin (`IASDT`). Supports multiple climate scenarios, outputting
data in TIFF and NetCDF formats. Orchestrated by `chelsa_process()`,
with helper functions `chelsa_prepare()` and `chelsa_project()`.

## Usage

``` r
chelsa_process(
  env_file = ".env",
  n_cores = 8L,
  strategy = "multisession",
  download = FALSE,
  overwrite = FALSE,
  download_attempts = 10L,
  sleep = 5L,
  other_variables = "npp",
  download_n_cores = 4,
  compression_level = 5,
  overwrite_processed = FALSE
)

chelsa_prepare(
  env_file = ".env",
  download = FALSE,
  n_cores = 8L,
  strategy = "multisession",
  overwrite = FALSE,
  download_attempts = 10L,
  sleep = 5L,
  other_variables = "npp"
)

chelsa_project(metadata = NULL, env_file = ".env", compression_level = 5L)
```

## Arguments

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

- download:

  Logical. If `TRUE`, downloads CHELSA files. Default: `FALSE`.

- overwrite:

  Logical. If `TRUE`, re-downloads existing files. Default: `FALSE`.

- download_attempts:

  Integer. Maximum download retries. Default: `10`.

- sleep:

  Integer. Seconds to wait between download attempts. Default: `5`.

- other_variables:

  Character. Additional variables to process (e.g., `"npp"` for Net
  Primary Productivity alongside 19 bioclimatic variables bio1-bio19).
  Use `""` for bioclimatic only. See
  [chelsa_variables](https://biodt.github.io/IASDT.R/reference/CHELSA_variables.md)
  for details. Default: `"npp"`.

- download_n_cores:

  Integer. Number of CPU cores to use for parallel downloading of CHELSA
  data. Only valid if download = `TRUE`. Defaults to 4.

- compression_level:

  Integer. NetCDF compression level (1 = least, 9 = most). Default: `5`.

- overwrite_processed:

  Logical. If `TRUE`, overwrites processed files. Default: `FALSE`.

- metadata:

  Tibble. Single-row metadata for input files, prepared by
  `chelsa_prepare()`

## Note

- `chelsa_prepare()` and `chelsa_project()` are internal helpers, not
  for direct use.

- Processes 19 bioclimatic variables (bio1–bio19) plus optional
  variables (e.g., NPP) for 46 scenarios (1 current, 45 future).

- Time-intensive; depends on file size and compute resources.

## Functions details

- **`chelsa_process()`**: Main function; optionally downloads CHELSA
  data, processes it to the European scale and reference grid, and saves
  TIFF and NetCDF outputs for 46 climate scenarios.

- **`chelsa_prepare()`**: Extracts metadata from local URL text files
  and manages optional downloads.

- **`chelsa_project()`**: Projects data to the `IASDT` reference grid
  with optional transformations.

## Author

Ahmed El-Gabbas
