# Process GBIF sampling effort data for the `IASDT`

Downloads and processes GBIF sampling effort data for vascular plants in
Europe, supporting the Invasive Alien Species Digital Twin (`IASDT`).
Orchestrated by `efforts_process()`, it uses helper functions to
request, download, split, summarise, and visualise data at the Order
level. The functions prepares raster maps for the number of vascular
plant observations and species per grid cell.

## Usage

``` r
efforts_process(
  env_file = ".env",
  r_environ = ".Renviron",
  request = TRUE,
  download = TRUE,
  n_cores = 6L,
  strategy = "multisession",
  start_year = 1981L,
  boundaries = c(-30, 50, 25, 75),
  chunk_size = 100000L,
  delete_chunks = TRUE,
  delete_processed = TRUE
)

efforts_request(
  env_file = ".env",
  n_cores = 3L,
  strategy = "multisession",
  start_year = 1981L,
  r_environ = ".Renviron",
  boundaries = c(-30, 50, 25, 75)
)

efforts_download(n_cores = 6L, strategy = "multisession", env_file = ".env")

efforts_summarize(
  env_file = ".env",
  n_cores = 6L,
  strategy = "multisession",
  chunk_size = 100000L,
  delete_chunks = TRUE
)

efforts_split(path_zip = NULL, env_file = ".env", chunk_size = 100000L)

efforts_plot(env_file = ".env")
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

- request:

  Logical. If `TRUE` (default), requests GBIF data; otherwise, loads
  existing data.

- download:

  Logical. If `TRUE` (default), downloads and saves GBIF data;
  otherwise, skips download. Default: `TRUE`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  6, except for `efforts_request`, which defaults to 3 with a maximum of
  3.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

- start_year:

  Integer. Earliest year for GBIF records (matches CHELSA climate data).
  Default: `1981`.

- boundaries:

  Numeric vector (length 4). GBIF data bounds (Left, Right, Bottom,
  Top). Default: `c(-30, 50, 25, 75)`.

- chunk_size:

  Integer. Rows per chunk file. Default: `100000`.

- delete_chunks:

  Logical. If `TRUE` (default), deletes chunk files post-processing.

- delete_processed:

  Logical. If `TRUE` (default), removes raw GBIF files after processing
  (\>22 GB).

- path_zip:

  Character. Path to zip file with CSV for splitting.

## Note

- `efforts_process()` is the main entry point for processing sampling
  effort data.

- Time-intensive (\>9 hours on 6-core Windows PC; GBIF request ~5
  hours).

- Detects and processes only new/missing data by order.

## Functions details

- **`efforts_process()`**: Manages the workflow for requesting,
  downloading, processing, and plotting GBIF vascular plant data.

- **`efforts_request()`**: Requests GBIF data by order in parallel.
  Stores results to disk.

- **`efforts_download()`**: Downloads GBIF data, validates files, and
  loads existing data if available. Returns a dataframe
  (`efforts_all_requests`) with paths.

- **`efforts_split()`**: Splits zipped CSV data by order into chunks,
  saving each separately.

- **`efforts_summarize()`**: Processes and summarises data into `RData`
  and TIFF rasters.

- **`efforts_plot()`**: Plots observation efforts (raw and log10
  scales).

## References

Data source: <https://www.gbif.org>

## Author

Ahmed El-Gabbas
