# Download and Process Topographic Wetness Index Data

Downloads, extracts, and processes the global topographic wetness index
data at 30 arc-second resolution (Title & Bemmels, 2018). The function
checks for existing processed data, downloads the required dataset if
necessary, extracts the relevant TIFF file, re-projects and masks it to
match a reference grid, and saves the processed raster in both TIFF and
`RData` formats.

## Usage

``` r
wetness_index_process(env_file = ".env")
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

## Value

(Invisibly) The path to the processed wetness index `RData` file.

## References

- Title & Bemmels (2018): ENVIREM: an expanded set of bioclimatic and
  topographic variables increases flexibility and improves performance
  of ecological niche modeling. Ecography.
  [10.1111/ecog.02880](https://doi.org/10.1111/ecog.02880)

- [ENVIREM](https://deepblue.lib.umich.edu/data/concern/generic_works/gt54kn05f):
  ENVIronmental Rasters for Ecological Modeling version 1.0

## Author

Ahmed El-Gabbas
