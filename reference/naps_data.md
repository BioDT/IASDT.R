# Process and map Naturalized Alien Plant Species (NAPS) data for the `IASDT`

Processes and visualises Naturalized Alien Plant Species (NAPS)
distribution data from GBIF, EASIN, and eLTER for the Invasive Alien
Species Digital Twin (`IASDT`). Merges pre-processed data, creates
presence-absence rasters, summarises distributions, and generates maps
using helper functions.

## Usage

``` r
naps_process(env_file = ".env", n_cores = 6L, strategy = "multisession")

naps_distribution(
  species = NULL,
  env_file = ".env",
  verbose = FALSE,
  dist_citizen = 100L
)

naps_plot(species = NULL, env_file = ".env")

naps_standardisation(env_file = ".env")
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  6.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

- species:

  Character. Species name for distribution mapping.

- verbose:

  Logical. If `TRUE`, prints progress messages. Default: `FALSE`.

- dist_citizen:

  Numeric. Distance in km for spatial filtering of citizen science data
  in GBIF or grid cells in countries at which species has not yet
  recognized as "naturalized". Default is `100L` km.

## Functions details

- **`naps_standardisation()`**: Load pre-prepared standardisation
  information for NAPS verbatim names and prepare a unique `ias_id` for
  each standardised NAPS. This function should be called **only once**
  per workflow version. The `ias_id` is determined based on sorting the
  standardised species names alphabetically within their taxonomic
  hierarchy (in the order of: class, order, family, and taxon_name). If
  a new standardisation file is used, new `ias_id` values will be
  generated, which will break reproducibility and prevent consistent
  data linkage across analyses.

- **`naps_process()`**: Merges pre-processed GBIF
  ([gbif_process](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)),
  EASIN
  ([easin_process](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)),
  and eLTER
  ([elter_process](https://biodt.github.io/IASDT.R/reference/eLTER_Process.md))
  data (run these first). Outputs `SpatRaster` distribution rasters,
  summary tables, and JPEG maps using `naps_distribution()` and
  `naps_plot()`.

- **`naps_distribution()`**: Generates presence-absence maps (`.RData`,
  `.tif`) for a species, including all grid cells in the study area and
  a set excluding cultivated/casual-only countries. Returns a file path
  to a tibble with presence counts (total, by source) and summary
  statistics for biogeographical regions

- **`naps_plot()`**: Creates JPEG distribution maps from GBIF, EASIN,
  and eLTER data using `ggplot2`.

## Author

Ahmed El-Gabbas
