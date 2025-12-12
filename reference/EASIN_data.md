# Process EASIN data for the `IASDT`

Extracts, processes, and visualises data from the [European Alien
Species Information Network (EASIN)](https://easin.jrc.ec.europa.eu/)
for the Invasive Alien Species Digital Twin (`IASDT`). Manages taxonomy,
occurrence data, and plots, handling API pagination and server limits.
Orchestrated by `easin_process()` with helpers `easin_taxonomy()`,
`easin_download()`, and `easin_plot()`.

## Usage

``` r
easin_process(
  extract_taxa = TRUE,
  extract_data = TRUE,
  n_download_attempts = 10L,
  n_cores = 6L,
  strategy = "multisession",
  sleep_time = 10L,
  n_search = 1000L,
  env_file = ".env",
  delete_chunks = TRUE,
  start_year = 1981L,
  plot = TRUE
)

easin_taxonomy(
  env_file = ".env",
  kingdom = "Plantae",
  phylum = "Tracheophyta",
  n_search = 100L
)

easin_download(
  species_key,
  timeout = 200,
  verbose = FALSE,
  env_file = ".env",
  n_search = 1000,
  n_attempts = 10,
  sleep_time = 5,
  delete_chunks = TRUE,
  return_data = FALSE
)

easin_plot(env_file = ".env")
```

## Arguments

- extract_taxa:

  Logical. If `TRUE`, extracts taxonomy using `easin_taxonomy()`.
  Default: `TRUE`.

- extract_data:

  Logical.If `TRUE`, downloads occurrence data with `easin_download()`.
  Default: `TRUE`.

- n_download_attempts:

  Integer. Retry attempts for downloads. Default: `10`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing.
  Default: 6. The maximum number of allowed cores are 8.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

- sleep_time:

  Integer. Number of seconds to pause between each data retrieval
  request to prevent overloading the server. Default: 5 second.

- n_search:

  Integer. Number of records to attempt to retrieve per request.
  Default: 1000, which is the current maximum allowed by the API.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- delete_chunks:

  Logical. Whether to delete temporary files for data chunks from the
  `file_parts` subdirectory. Defaults to `TRUE`.

- start_year:

  Integer. Earliest year for occurrence data (excludes earlier records).
  Default: `1981` (aligned with CHELSA climate data).

- plot:

  Logical. If `TRUE`, generates plots via `easin_plot()`. Default:
  `TRUE`.

- kingdom:

  Character. Taxonomic kingdom to query. Default: `"Plantae"`.

- phylum:

  Character. Taxonomic phylum within kingdom. Default: `"Tracheophyta"`

- species_key:

  Character. EASIN taxon ID for which data is to be retrieved. This
  parameter cannot be `NULL`.

- timeout:

  Integer. Download timeout in seconds. Default: `200`.

- verbose:

  Logical. If `TRUE`, prints progress messages. Default: `FALSE`.

- n_attempts:

  Integer. Max download attempts per chunk. Default: `10`.

- return_data:

  Logical. If `TRUE`, returns data as a dataframe; otherwise, saves to
  disk and returns `invisible(NULL)`. Default: `FALSE`.

## Note

Uses a static RDS file with EASIN-GBIF taxonomic standardization,
prepared by Marina Golivets (Feb 2024).

## Functions details

- **`easin_process()`**: Orchestrates taxonomy extraction, data
  downloads, and plotting for EASIN species data.

- **`easin_taxonomy()`**: Fetches taxonomy data in chunks via the easin
  API, filtered by kingdom and phylum. Returns a tibble.

- **`easin_download()`**: Downloads occurrence data for a given easin
  ID, handling pagination and pauses. Returns a dataframe if
  `return_data = TRUE`, else `invisible(NULL)`.

- **`easin_plot()`**: Creates summary plots (observations count, species
  count, distribution by partner) as JPEGs. Returns `invisible(NULL)`.

## Author

Ahmed El-Gabbas
