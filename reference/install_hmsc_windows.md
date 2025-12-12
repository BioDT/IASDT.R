# Install Hmsc-HPC in a python virtual environment on Windows

This function sets up a Python virtual environment and installs the
[`Hmsc-HPC`](https://github.com/aniskhan25/hmsc-hpc) package from a
specified Git repository. It also installs `TensorFlow` and performs
checks to to verify the installation of Python, the virtual environment,
and the packages.

## Usage

``` r
install_hmsc_windows(
  path_python,
  path_ve,
  url_hmsc = "https://github.com/trossi/hmsc-hpc.git@simplify-io",
  url_rdata = "https://github.com/trossi/rdata.git@test-for-hmsc-v2",
  force = FALSE
)
```

## Arguments

- path_python:

  Character. Path to the Python executable.

- path_ve:

  Character. Path where the virtual environment will be created. This
  can not be an existing folder.

- url_hmsc:

  Character. URL of the Git repository (and branch name) of the
  `Hmsc-HPC` package.

- url_rdata:

  Character. URL of the Git repository (and branch name) of the `rdata`
  package. This is temporary to allow `rdata` package to write rds file.
  In the near future, this functionality will be pushed to the main
  branch of `rdata`.

- force:

  Logical. Whether to force the installation of `Hmsc-HPC`; i.e. using
  `--force-reinstall` as a suffix to the `pip install` command

## Value

The function performs installation steps and returns NULL invisibly.

## Details

The function performs the following steps:

- Checks if the virtual environment directory already exists and stops
  with an error if it does

- Verifies the Python version and installation

- Creates a new Python virtual environment

- Upgrades `pip` in the virtual environment

- Installs the `Hmsc-HPC` package from the provided Git URL

- Installs `TensorFlow` version 2.15

- Installs `rdata`

- Checks `TensorFlow` and `Hmsc-HPC` package installations

## Author

Ahmed El-Gabbas

## Examples

``` r
if (FALSE) { # \dontrun{
install_hmsc_windows(
   path_python = "C:/Python/Python310/python.exe",
   path_ve = "D:/Hmsc-HPC")
} # }
```
