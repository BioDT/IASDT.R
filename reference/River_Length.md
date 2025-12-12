# Calculate the length of rivers in each Strahler order per grid cell

This function processes EU-Hydro River Network Database to calculate the
length of rivers in each Strahler number. The Strahler number is used as
an index for river network classification, with higher numbers
representing larger, more significant river segments. The function reads
and processes zip-compressed geographic data (GPKG files), extracts
relevant information about river segments, computes the length of rivers
for each Strahler order per grid cell, and outputs the results both as
raster files and RData objects. The calculated length represents the
total length of rivers in each Strahler number or larger (e.g., for
STRAHLER_5, the length of rivers with Strahler values of 5 or higher).

## Usage

``` r
river_length(env_file = ".env", cleanup = FALSE)
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- cleanup:

  Logical indicating whether to clean up temporary files from the
  Interim directory after finishing calculations. Default: `FALSE`.

## Value

`NULL`. The function outputs processed files to the specified
directories.

## Details

The data provides at pan-European level a photo-interpreted river
network, consistent of surface interpretation of water bodies (lakes and
wide rivers), and a drainage model (also called Drainage Network),
derived from EU-DEM, with catchments and drainage lines and nodes.

- **Data source**: EU-Hydro River Network Database v013 \|

- **Temporal extent**: 2006-2012; **Format**: Vector (GPKG); **Size**: 4
  GB

## References

- DOI: <https://doi.org/10.2909/393359a7-7ebd-4a52-80ac-1a18d5f3db9c>

- Download link:
  <https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database>

## Author

Ahmed El-Gabbas
