
<!-- badges: start -->

[![R-CMD-check](https://github.com/BioDT/IASDT.R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BioDT/IASDT.R/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# IASDT.R: Modelling the distribution of invasive alien plant species in Europe

## Overview

[`IASDT.R`](https://biodt.github.io/IASDT.R) supports the invasive alien
species (IAS) Digital Twin (`IASDT`), as part of the EU-funded
<a href="https://biodt.eu/" target="_blank">BioDT</a> project. This R
package provides functions for processing input data (biotic and
abiotic) and modelling the distribution of invasive alien plant species
across Europe using joint species distribution models. By the end of the
project, it is planned that the package will enable end users to access
model outputs (e.g., prediction maps) directly from R.

- More information on the BioDT project can be found at this
  <a href="https://biodt.eu/" target="_blank">link</a>.
- For an overview on the IASDT, see Khan, El-Gabbas, et
  al. <a href="https://doi.org/10.3897/rio.10.e124579" target="_blank"><img src="https://img.shields.io/badge/DOI:-10.3897/rio.10.e124579-blue" alt="DOI: 10.3897/rio.10.e124579"/></a>.
- Documentation of the package functions can be found
  <a href="https://biodt.github.io/IASDT.R/reference/index.html" target="_blank">here</a>.
- Click <a href="https://github.com/BioDT/uc-ias-workflows">here</a> for
  the workflow of the `IASDT`.

<center>

<img src="https://git.ufz.de/uploads/-/system/group/avatar/4444/biodt.png" alt="BioDT_logo" width="400">
</center>

<hr>

## Installing the package

The `IASDT.R` package is currently under development. To install the
most-recent development version of the package, use the following
command:

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

If you use the `IASDT.R` package, please cite it as:
<a href="https://doi.org/10.5281/zenodo.14834384" target="_blank"><img role="button" tabindex="0" id="modal-858828210-trigger" aria-controls="modal-858828210" aria-expanded="false" class="doi-modal-trigger block m-0" src="https://zenodo.org/badge/DOI/10.5281/zenodo.14834384.svg" alt="DOI: 10.5281/zenodo.14834384"/></a>

> El-Gabbas, A. (2025). **IASDT.R: Modelling the distribution of
> invasive alien plant species in Europe**.
> <a href="https://doi.org/10.5281/zenodo.14834384" target="_blank">10.5281/zenodo.14834384</a>,
> R package version 0.1.05;
> <a href="https://github.com/BioDT/IASDT.R" target="_blank">https://github.com/BioDT/IASDT.R</a>;
> <a href="https://biodt.eu" target="_blank">https://biodt.eu</a>.

<span style="     color: grey !important;">Last update:
2025-07-15</span>
