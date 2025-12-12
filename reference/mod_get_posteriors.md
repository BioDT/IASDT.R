# Combines posteriors exported by `Hmsc-HPC` into an Hmsc object

This function converts posterior files exported by `Hmsc-HPC` into an
Hmsc object. It can either read the data directly from RDS files or
convert it from JSON format if specified.

## Usage

``` r
mod_get_posteriors(path_posterior = NULL, from_json = FALSE)
```

## Arguments

- path_posterior:

  Character vector. Path to the RDS files containing the exported
  posterior files. This argument is mandatory and cannot be empty.

- from_json:

  Logical. Whether the loaded models should be converted from `JSON`
  format. Defaults to `FALSE`, meaning the data will be read directly
  from RDS files without conversion.

## Value

Depending on the `from_json` parameter, returns an Hmsc object either
directly from the RDS files or after converting it from JSON format.

## Author

Ahmed El-Gabbas
