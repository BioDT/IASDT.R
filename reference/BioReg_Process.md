# Process biogeographical regions dataset

Downloads and processes the Biogeographical Regions dataset (Europe
2016, v1) from the [European Environment
Agency](https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618).
This function extracts biogeographical region names per reference grid
cell for use in counting species presence across biogeographical
regions.

## Usage

``` r
bioreg_process(env_file = ".env")
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

## Value

Invisible `NULL`. Processed data is saved to disk as raster, vector, and
RData files.

## Details

- *Temporal coverage*: 2011-2015

- *Spatial coverage*: 28°E to 81°E, 31.27°W to 62°E

- *CRS*: EPSG:3035

- *file format*: shapefile (compressed in zip file)

- *Requirements*: `curl` (download) and `unzip` (extraction)

## Author

Ahmed El-Gabbas
