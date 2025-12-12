# Process Corine Land Cover (CLC) data for the `IASDT`

Processes [Corine Land Cover
(CLC)](https://land.copernicus.eu/pan-european/corine-land-cover/clc2018)
data for the Invasive Alien Species Digital Twin (`IASDT`). Calculates
percentage coverage and most common classes per grid cell at three CLC
levels, plus `eunis19` and `SynHab` habitat types. Prepares a reference
grid and optionally generates percentage coverage maps as JPEG.

## Usage

``` r
clc_process(
  env_file = ".env",
  min_land_percent = 15L,
  plot_clc = TRUE,
  n_cores_plot = 5L,
  strategy = "multisession"
)
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- min_land_percent:

  Numeric. Minimum land percentage per grid cell for the reference grid.
  Default: `15`.

- plot_clc:

  Logical. If `TRUE`, plots percentage coverage for CLC levels and
  habitat types. Default: `TRUE`.

- n_cores_plot:

  Integer. Number of CPU cores to use for parallel plotting. Default:
  5L.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

## Value

Returns `invisible(NULL)`; saves processed data and optional plots to
disk.

## References

- Data source:
  <https://land.copernicus.eu/pan-european/corine-land-cover/clc2018>

- Data citation:
  <https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac>

## Author

Ahmed El-Gabbas
