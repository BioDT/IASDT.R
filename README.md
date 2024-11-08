
# IASDT.R R-package

The `IASDT.R` R package supports the invasive alien species (IAS)
prototype Digital Twin (`IAS-pDT`), as part of the EU-funded
[`BioDT`](https://biodt.eu/) project. This package provides functions
for processing input data (biotic and abiotic) and modeling the
distribution of IAS across Europe using joint species distribution
models. By the end of the project, it is planned that the package will
enable end users to access model outputs (e.g., prediction maps)
directly from R.

More information on the BioDT project can be found at this
[`link`](https://biodt.eu/). For more information on the IAS-pDT, see
[`Khan, El-Gabbas, et al. (2024)`](https://doi.org/10.3897/rio.10.e124579)
<br/>

<center>
<img
src="https://git.ufz.de/uploads/-/system/group/avatar/4444/biodt.png"
width="352" />
</center>
<hr>

## Installing the package

The `IASDT.R` package is currently in development and is private. To
install the package, use the following command:

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

## Contribute to the package

We welcome contributions and bug reports to the package. Please make
changes on a new branch and submit a pull request.

> For questions, please get in touch with
> [me](https://elgabbas.netlify.app/) at `ahmed.el-gabbas[at]ufz[dot]de`

<span style="     color: grey !important;">Last update:
2024-11-08</span>
