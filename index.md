# IASDT.R: Modelling the distribution of invasive alien plant species in Europe

## Overview

[`IASDT.R`](https://biodt.github.io/IASDT.R) supports the invasive alien
species (IAS) Digital Twin (`IASDT`), as part of the EU-funded
[BioDT](https://biodt.eu/) project. This R package provides functions
for processing input data (biotic and abiotic) and modelling the
distribution of naturalised alien plant species (NAPS) across Europe
using joint species distribution models.

- More information on the `BioDT` project can be found at this
  [link](https://biodt.eu/).
- For an overview on the `IASDT`, see Khan, El-Gabbas, et al. [![DOI:
  10.3897/rio.10.e124579](https://img.shields.io/badge/DOI:-10.3897/rio.10.e124579-blue)](https://doi.org/10.3897/rio.10.e124579).
- Documentation of the package functions can be found
  [here](https://biodt.github.io/IASDT.R/reference/index.html).
- Click [here](https://github.com/BioDT/uc-ias-workflows) for the
  workflow of the `IASDT`.
- The package also depends on
  [ecokit](https://elgabbas.github.io/ecokit/) for general ecological
  data analysis and visualisation.

![BioDT_logo](https://git.ufz.de/uploads/-/system/group/avatar/4444/biodt.png)

------------------------------------------------------------------------

## Installing the package

To install the most-recent development version of the package, use the
following command:

``` r
remotes::install_github(repo = "BioDT/IASDT.R", dependencies = TRUE)
```

To update the package, use:

``` r
remotes::update_packages("IASDT.R", dependencies = TRUE)
```

If you are using renv, update the package with:

``` r
renv::update("IASDT.R", prompt = FALSE)
```

------------------------------------------------------------------------

## Contribute to the package

Contributions, suggestions, and bug reports are welcome. Please make
changes on a new branch and submit a pull request or make an issue
[here](https://github.com/BioDT/IASDT.R/issues).

> For questions, please get in touch with
> [me](https://elgabbas.netlify.app/) at `elgabbas[at]outlook[dot]com`

------------------------------------------------------------------------

## Citation

If you use the `IASDT.R` package, please cite it as: [![DOI:
10.5281/zenodo.14834384](https://zenodo.org/badge/DOI/10.5281/zenodo.14834384.svg)](https://doi.org/10.5281/zenodo.14834384)

> El-Gabbas, A. (2025). **IASDT.R: Modelling the distribution of
> invasive alien plant species in Europe**.
> [10.5281/zenodo.14834384](https://doi.org/10.5281/zenodo.14834384), R
> package version 0.1.05;
> [https://github.com/BioDT/IASDT.R](https://github.com/BioDT/IASDT.R);
> [https://biodt.eu](https://biodt.eu).

Last update: 2025-12-12
