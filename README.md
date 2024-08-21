
# IASDT.R R-package

The `IASDT.R` is an R package that aims at supporting the invasive alien
species (IAS) prototype Digital Twin (`IAS-pDT`), as part of the
EU-funded [BioDT](https://biodt.eu/) project. These functions help in
the manipulation of IAS input data and the modelling the distribution of
IAS across Europe using joint species distribution models. It is planned
that the package will allow end users to access model outputs
(e.g. prediction maps) from within R by the end of 2024.

This package is currently hosted at a private
[GitLab](https://git.ufz.de/biodt/IASDT.R) repository, but will be moved
to the [BioDT GitHub Organisation](https://github.com/BioDT) soon.

More information on the BioDT project can be found at this
[link](https://biodt.eu/) <br/><br/>

For more information on the IAS-pDT, see [Taimur, El-Gabbas, et
al. (2024)](https://doi.org/10.3897/rio.10.e124579) <br/><br/>

<center>

<figure>
<img
src="https://git.ufz.de/uploads/-/system/group/avatar/4444/biodt.png"
width="500" alt="bioDT" />
<figcaption aria-hidden="true">bioDT</figcaption>
</figure>

</center>
<hr>

## Installing the package

The `IASDT.R` package is still at the development stage and currently
stores some helper function needed for the input data processing and
modelling of the IAS. The package is currently private (hosted at:
<https://git.ufz.de/>). This means that only users with access rights to
the git repository can use the package using a valid private access
token (PAT) authentication.

### Set PAT

To install the package, you need to set a valid gitlab PAT first. You
can create a PAT by visiting this
[link](https://git.ufz.de/biodt/IASDT.R/-/settings/access_tokens). After
creating the PAT, you need to add it to the `.Renviron` file. For more
information on the `.Renviron` file, click
[here](https://support.posit.co/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf).

In an interactive RStudio session, you add the PAT by using:

``` r
require(usethis)
require(remotes)
```

``` r
usethis::edit_r_environ()                             # global `.Renviron` file
# or alternatively
# usethis::edit_r_environ(scope = "project")          # project-specific settings
```

This opens the `.Renviron` file. You need to add the following line to
the file (replace `XXXXXXXXXX` with your valid PAT).

``` r
GITLAB_PAT=XXXXXXXXXX
```

Then you need to save the `.Renviron` file and ensure that it ends with
a newline.

RStudio should be restarted for making such changes effective. After you
restart RStudio, you can check the value of the `GITLAB_PAT` environment
variable by using the following command:

``` r
Sys.getenv("GITLAB_PAT")
# XXXXXXXXXX
```

### Install the package

To install the package you can use the following command

``` r
remotes::install_gitlab(
    repo = "BioDT/IASDT.R", host = "git.ufz.de", dependencies = TRUE)
```

Alternatively, you can provide the PAT directly to the
`remotes::install_gitlab` function; however, this is not recommended.

``` r
remotes::install_gitlab(
    repo = "BioDT/IASDT.R", host = "git.ufz.de", auth_token = "XXXXXXXXXX", 
    dependencies = TRUE)
```

### Update the package

To update the package, use the following

``` r
remotes::update_packages("IASDT.R", dependencies = TRUE)
```

or re-install it; you can use `force = TRUE` to force installation, even
if the remote state has not changed since the previous install.

``` r
remotes::install_gitlab(
    repo = "BioDT/IASDT.R", host = "git.ufz.de", dependencies = TRUE)
```

If you are using `renv`, you may use the following to update the package

``` r
renv::update("IASDT.R", prompt = FALSE)
```

### Contribute to the package

Your contribution to the package is welcome. Please make your changes to
a different branch and submit a pull request.

<hr>
<center>

For questions, please contact [me](https://elgabbas.netlify.app/) at:
`ahmed.el-gabbas[at]ufz[dot]de`

<span style="     color: grey !important;">Last update:
2024-08-21</span>

</center>
