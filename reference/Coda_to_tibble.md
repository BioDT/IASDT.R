# Convert a Coda object to a tibble with specified parameter transformations

This function converts a Coda object (`mcmc.list` or `mcmc`) into a
tibble format, facilitating further analysis and visualisation. It
supports transformation for specific parameter types: `rho`, `alpha`,
`omega`, and `beta`.

## Usage

``` r
coda_to_tibble(
  coda_object = NULL,
  posterior_type = NULL,
  env_file = ".env",
  n_omega = 100
)
```

## Arguments

- coda_object:

  An object of class `mcmc.list` or `mcmc`, representing the MCMC
  output.

- posterior_type:

  Character. The parameter type to transform and extract. Must be one of
  `rho`, `alpha`, `omega`, or `beta`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- n_omega:

  Integer. The number of species to be sampled for the `Omega` parameter
  transformation. Defaults to 100.

## Value

A tibble containing the transformed parameters based on the specified
`posterior_type`. The structure of the returned tibble varies depending
on the `posterior_type` parameter.

## Author

Ahmed El-Gabbas

## Examples

``` r
ecokit::load_packages(Hmsc, coda, dplyr, data.table)

coda_object <- Hmsc::convertToCodaObject(Hmsc::TD$m)
ecokit::ht(
  IASDT.R::coda_to_tibble(
    coda_object = coda_object$Alpha[[1]], posterior_type = "Alpha"))
#>                Alpha alpha_num  factor  chain  iter value
#>               <fctr>    <fctr>  <fctr> <fctr> <int> <num>
#>   1: Alpha1[factor1]    Alpha1 factor1      1    51     0
#>   2: Alpha1[factor1]    Alpha1 factor1      1    52     0
#>   3: Alpha1[factor1]    Alpha1 factor1      1    53     0
#>   4: Alpha1[factor1]    Alpha1 factor1      1    54     0
#>   5: Alpha1[factor1]    Alpha1 factor1      1    55     0
#>  ---                                                     
#> 396: Alpha1[factor2]    Alpha1 factor2      2   146     0
#> 397: Alpha1[factor2]    Alpha1 factor2      2   147     0
#> 398: Alpha1[factor2]    Alpha1 factor2      2   148     0
#> 399: Alpha1[factor2]    Alpha1 factor2      2   149     0
#> 400: Alpha1[factor2]    Alpha1 factor2      2   150     0
```
