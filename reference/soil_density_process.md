# Retrieve and project soil bulk density data

Downloads and projects [SoilGrids](https://soilgrids.org)' bulk density
of the fine earth fraction (bdod) layers for specified depth intervals.
See [here](https://www.isric.org/explore/soilgrids/faq-soilgrids) for
more details. This function project the original data at different
depths intervals to the modelling reference grid.

## Usage

``` r
soil_density_process(depths = NULL, env_file = ".env", n_cores = 6L)
```

## Arguments

- depths:

  Character vector of depth intervals (e.g. `c("0-5","5-15")`). If
  `NULL` (default) all valid depths are processed. Valid `depths` are
  the SoilGrids BDOD standard horizons (cm): "0-5", "5-15", "15-30",
  "30-60", "60-100", and "100-200". Each depth string should omit the
  trailing `cm`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  6.

## Value

(Invisibly) the path to the `RData` file containing the processed soil
bulk density data.

## References

- SoilGrids: <https://soilgrids.org>

- SoilGrids250m 2.0 - [Bulk
  density](https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/713396f9-1687-11ea-a7c0-a0481ca9e724)

- Poggio *et al.* (2021): SoilGrids 2.0: producing soil information for
  the globe with quantified spatial uncertainty. Soil.
  [10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021)

## Author

Ahmed El-Gabbas
