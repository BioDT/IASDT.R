# Get species name or information of an `IASDT` species ID

This function retrieves detailed information on `IASDT` species list,
optionally filtered by a specific `IASDT` species ID (`species_id`).

## Usage

``` r
get_species_name(species_id = NULL, env_file = ".env")
```

## Arguments

- species_id:

  optional IASDT species ID for which detailed information is required.
  If not provided, the function returns the entire species list.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

## Value

A data frame containing species information. If a species ID
`species_id` is provided, it only returns species information for the
listed species, otherwise return the full list of IAS.

## Author

Ahmed El-Gabbas
