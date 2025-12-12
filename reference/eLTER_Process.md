# Process eLTER data for the `IASDT`

This function processes pre-cleaned and pre-standardized Integrated
European Long-Term Ecosystem, critical zone and socio-ecological
Research ([eLTER](https://elter-ri.eu/)) data for the Invasive Alien
Species Digital Twin (`IASDT`).

## Usage

``` r
elter_process(env_file = ".env", start_year = 1981)
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- start_year:

  Numeric. The starting year for the occurrence data. Only records from
  this year onward will be processed. Default is `1981`, which matches
  the year ranges of CHELSA current climate data.

## Value

Returns `NULL` invisibly after saving the processed data.

## Note

This function processes pre-cleaned vascular plants data from eLTER
sites, harmonized by Ahmed El-Gabbas. The original eLTER biodiversity
data were highly heterogeneous in format and structure, requiring
standardization and cleaning before use. Taxonomic standardization with
the GBIF backbone was performed by Marina Golivets (Feb. 2024).

## Author

Ahmed El-Gabbas
