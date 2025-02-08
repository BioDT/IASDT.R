
# IASDT.R: Modelling the distribution of invasive alien plant species in Europe

<hr>

## Overview

[`IASDT.R`](https://biodt.github.io/IASDT.R/index.html) supports the
invasive alien species (IAS) prototype Digital Twin (`IAS-pDT`), as part
of the EU-funded [`BioDT`](https://biodt.eu/) project. This package
provides functions for processing input data (biotic and abiotic) and
modeling the distribution of invasive alien plant species across Europe
using joint species distribution models. By the end of the project, it
is planned that the package will enable end users to access model
outputs (e.g., prediction maps) directly from R.

- More information on the BioDT project can be found at this
  [`link`](https://biodt.eu/).

- For more information on the IAS-pDT, see Khan, El-Gabbas, et
  al. (2024)\[[`Link`](https://doi.org/10.3897/rio.10.e124579)\]

- The documentation of the package can be found
  [here](https://biodt.github.io/IASDT.R).

<center>
<img src="https://git.ufz.de/uploads/-/system/group/avatar/4444/biodt.png" alt="BioDT_logo" width="352">
</center>
<hr>

## Installing the package

The `IASDT.R` package is currently under development. To install the
package, use the following command:

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

<hr>

## Contribute to the package

Contributions, suggestions, and bug reports are welcome. Please make
changes on a new branch and submit a pull request or make an issue
[here](https://github.com/BioDT/IASDT.R/issues).

> For questions, please get in touch with
> [me](https://elgabbas.netlify.app/) at `ahmed.el-gabbas[at]ufz[dot]de`

<hr>

## Citation

If you use the `IASDT.R` package in your research, please cite it as:

> El-Gabbas A. (2025). IASDT.R: Modelling the distribution of invasive
> alien plant species in Europe. R package version 0.1.02. DOI:
> [10.5281/zenodo.14834385](https://doi.org/10.5281/zenodo.14834385).
> <https://github.com/BioDT/IASDT.R>; <https://biodt.eu>.

<span style="     color: grey !important;">Last update:
2025-02-07</span>
