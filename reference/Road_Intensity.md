# Calculate road intensity per grid cell

This function downloads, processes, and analyses [GRIP global roads
data](https://www.globio.info/download-grip-dataset) ([Meijer et al.
2018](https://iopscience.iop.org/article/10.1088/1748-9326/aabd42/meta)).
The function calculates the total road lengths and the distance to the
nearest road per grid cell (for any road type and per road type).

## Usage

``` r
road_intensity(env_file = ".env")
```

## Arguments

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

## Value

`NULL`. The function outputs processed files to the specified
directories.

## Note

- The function downloads the most recent version of Global Roads
  Inventory Project (`GRIP`) data from the URL specified in the
  environment variable `DP_R_roads_url`. Original data format is a
  zipped file containing global road data in the form of `fgdb`
  (`EPSG:3246`).

- On LUMI HPC, loading the `libarchive` module is necessary to use the
  `archive` R package: `module load libarchive/3.6.2-cpeGNU-23.09`

- The distance to roads is calculated by determining the distance from
  each grid cell to the nearest grid cell that overlaps with a road (not
  to the nearest road line). Note that this is different from
  calculating the actual distance to the nearest road line, which is
  computationally intensive and not performed in this function.

## Author

Ahmed El-Gabbas
