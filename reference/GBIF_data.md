# Process GBIF occurrence data for the `IASDT`

Extracts, processes, and visualises occurrence data from the [Global
Biodiversity Information Facility (GBIF)](https://www.gbif.org) for the
Invasive Alien Species Digital Twin (`IASDT`). Orchestrated by
`gbif_process()`, it requests, downloads, cleans, chunks, and maps
species data using helper functions.

## Usage

``` r
gbif_process(
  env_file = ".env",
  r_environ = ".Renviron",
  n_cores = 6L,
  strategy = "multisession",
  request = TRUE,
  download = TRUE,
  split_chunks = TRUE,
  overwrite = FALSE,
  delete_chunks = TRUE,
  chunk_size = 50000L,
  boundaries = c(-30, 50, 25, 75),
  start_year = 1981L
)

gbif_download(
  env_file = ".env",
  r_environ = ".Renviron",
  request = TRUE,
  download = TRUE,
  split_chunks = TRUE,
  chunk_size = 50000L,
  boundaries = c(-30L, 50L, 25L, 75L),
  start_year = 1981L
)

gbif_read_chunk(
  chunk_file,
  env_file = ".env",
  max_uncertainty = 10L,
  start_year = 1981L,
  save_rdata = TRUE,
  return_data = FALSE,
  overwrite = FALSE
)

gbif_species_data(
  species = NULL,
  env_file = ".env",
  verbose = TRUE,
  plot_tag = NULL
)
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- r_environ:

  Character. Path to `.Renviron` file with GBIF credentials
  (`GBIF_EMAIL`, `GBIF_USER`, `GBIF_PWD`). Default: `".Renviron"`. The
  credentials must be in the format:

  - `GBIF_EMAIL=your_email`

  - `GBIF_USER=your_username`

  - `GBIF_PWD=your_password`

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

- request:

  Logical. If `TRUE` (default), requests GBIF data; otherwise, loads
  from disk.

- download:

  Logical. If `TRUE` (default), downloads and saves GBIF data.

- split_chunks:

  Logical. If `TRUE` (default), splits data into chunks for easier
  processing.

- overwrite:

  Logical. If `TRUE`, reprocesses existing `.RData` chunks. Default:
  `FALSE`. This helps to continue working on previously processed chunks
  if the previous try failed, e.g. due to memory issue.

- delete_chunks:

  Logical. If `TRUE` (default), deletes chunk files.

- chunk_size:

  Integer. Records per data chunk. Default: `50000`.

- boundaries:

  Numeric vector (length 4). GBIF data bounds (Left, Right, Bottom,
  Top). Default: `c(-30, 50, 25, 75)`.

- start_year:

  Integer. Earliest collection year to be included. Default is 1981.

- chunk_file:

  Character. Path of chunk file for processing.

- max_uncertainty:

  Numeric. Maximum spatial uncertainty in kilometres. Default: `10`.

- save_rdata:

  Logical. If `TRUE` (default), saves chunk data as `.RData`.

- return_data:

  If `TRUE`, returns chunk data; otherwise, `invisible(NULL)`. Default:
  `FALSE`.

- species:

  Character. Species name for processing.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

- plot_tag:

  Character. Tag for plot titles.

## Note

Relies on a static RDS file listing IAS species, GBIF keys, and
metadata, standardized by Marina Golivets (Feb 2024).

## Functions details

- **`gbif_process()`**: Orchestrates GBIF data requests, downloads,
  processing, and mapping. Saves `RData`, Excel, and JPEG summary files.

- **`gbif_download()`**: Requests and downloads GBIF data (if
  `download = TRUE`), using the specified criteria (taxa, coordinates,
  time period, and boundaries), splits into small chunks (if
  `split_chunks = TRUE`), and saves metadata. Returns `invisible(NULL)`.

- **`gbif_read_chunk()`**: Filters chunk data (spatial/temporal, e.g.,
  spatial uncertainty, collection year, coordinate precision, and
  taxonomic rank), select relevant columns, and saves as `.RData` (if
  `save_rdata = TRUE`) or returns it (if `return_data = TRUE`). Skips if
  `.RData` exists and `overwrite = FALSE`.

- **`gbif_species_data()`**: Converts species-specific data to `sf` and
  raster formats, generating distribution maps.

## Author

Ahmed El-Gabbas
