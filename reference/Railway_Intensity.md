# Calculate railway intensity based on `OpenStreetMap` data

This function downloads, processes, and analyses railway data extracted
from [OpenRailwayMap](https://www.openrailwaymap.org) available from
[OpenStreetMap Data Extracts](https://download.geofabrik.de/). It
supports parallel processing for faster execution and can calculate the
total length of railways and distance to the nearest railway for each
grid cell in Europe.

## Usage

``` r
railway_intensity(
  env_file = ".env",
  n_cores = 6L,
  strategy = "multisession",
  delete_processed = TRUE
)
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

- delete_processed:

  Logical indicating whether to delete the raw downloaded railways files
  after processing them. This helps to free large unnecessary file space
  (\> 55 GB). Defaults to `TRUE`.

## Value

`NULL`. Outputs processed files to the directories specified in the
environment file.

## Author

Ahmed El-Gabbas
