## |------------------------------------------------------------------------| #
# install_hmsc_windows ----
## |------------------------------------------------------------------------| #

#' Install Hmsc-HPC in a python virtual environment on Windows
#'
#' This function sets up a Python virtual environment and installs the
#' [`Hmsc-HPC`](https://github.com/aniskhan25/hmsc-hpc) package from a specified
#' Git repository. It also installs `TensorFlow` and performs checks to to
#' verify the installation of Python, the virtual environment, and the packages.
#' @param path_python Character. Path to the Python executable.
#' @param path_ve Character. Path where the virtual environment will be created.
#'   This can not be an existing folder.
#' @param url_hmsc Character. URL of the Git repository (and branch name) of the
#'   `Hmsc-HPC` package.
#' @param url_rdata Character. URL of the Git repository (and branch name) of
#'   the `rdata` package. This is temporary to allow `rdata` package to write
#'   rds file. In the near future, this functionality will be pushed to the
#'   main branch of `rdata`.
#' @param force Logical. Whether to force the installation of `Hmsc-HPC`; i.e.
#'   using `--force-reinstall` as a suffix to the `pip install` command
#' @return The function performs installation steps and returns NULL invisibly.
#' @export
#' @name install_hmsc_windows
#' @author Ahmed El-Gabbas
#' @details The function performs the following steps:
#' - Checks if the virtual environment directory already exists and stops with
#'   an error if it does
#' - Verifies the Python version and installation
#' - Creates a new Python virtual environment
#' - Upgrades `pip` in the virtual environment
#' - Installs the `Hmsc-HPC` package from the provided Git URL
#' - Installs `TensorFlow` version 2.15
#' - Installs `rdata`
#' - Checks `TensorFlow` and `Hmsc-HPC` package installations
#' @examples
#' \dontrun{
#' install_hmsc_windows(
#'    path_python = "C:/Python/Python310/python.exe",
#'    path_ve = "D:/Hmsc-HPC")
#' }

install_hmsc_windows <- function(
    path_python, path_ve,
    url_hmsc = "https://github.com/trossi/hmsc-hpc.git@simplify-io",
    url_rdata = "https://github.com/trossi/rdata.git@test-for-hmsc-v2",
    force = FALSE) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check if the virtual environment directory already exists
  if (fs::dir_exists(path_ve)) {
    ecokit::stop_ctx(
      "Path to the virtual environment already exists", path_ve = path_ve,
      include_backtrace = TRUE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Python installation

  ecokit::info_chunk("Checking Python")
  python_version <- system2(
    path_python, "--version", stdout = TRUE, stderr = TRUE)

  if (!startsWith(python_version, "Python")) {
    ecokit::stop_ctx(
      "Python was not installed or not found at the provided path",
      path_python = path_python, include_backtrace = TRUE)
  }

  cat(paste0("  >>  ", python_version, "\n"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating the virtual environment

  ecokit::info_chunk("Creating virtual environment")
  system2(path_python, paste0("-m venv ", path_ve))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Python version in the virtual environment

  ecokit::info_chunk("Checking virtual environment")
  python <- fs::path(path_ve, "Scripts", "python.exe")
  python_version <- system2(python, "--version", stdout = TRUE, stderr = TRUE)

  if (!startsWith(python_version, "Python")) {
    ecokit::stop_ctx(
      "Python was not installed in the virtual environment",
      python_version = python_version, include_backtrace = TRUE)
  }

  cat(paste0("  >>  ", python_version, "\n"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Upgrading pip

  ecokit::info_chunk("Upgrading pip")
  system2(python, "-m pip install --upgrade pip")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Installing Hmsc-HPC package from Git

  ecokit::info_chunk("Installing Hmsc-HPC package from Git")

  if (force) {
    install_command <- paste0(
      "-m pip install git+", url_hmsc, " --force-reinstall")
  } else {
    install_command <- paste0("-m pip install git+", url_hmsc)
  }
  system2(python, install_command)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Install `TensorFlow` v2.15

  # This is to avoid the following error: ImportError: This version of
  # TensorFlow Probability requires TensorFlow version >= 2.15; Detected an
  # installation of version 2.13.1. Please upgrade TensorFlow to proceed.

  ecokit::info_chunk("Installing `TensorFlow` v2.15")
  system2(python, "-m pip install tensorflow==2.15")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking TensorFlow
  ecokit::info_chunk("Checking `TensorFlow`")
  system2(python, '-c "import tensorflow as tf; print(tf.constant(1))"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Install `rdata`
  ecokit::info_chunk("Installing rdata")
  install_command <- paste0("-m pip install git+", url_rdata)
  system2(python, install_command)

  # It is necessary to downgrade the version of numpy to work with tensorflow

  # ERROR: pip's dependency resolver does not currently take into account all
  # the packages that are installed. This behaviour is the source of the
  # following dependency conflicts. tensorflow-intel 2.15.0 requires
  # numpy<2.0.0, >=1.23.5, but you have numpy 2.0.2 which is incompatible.

  system2(python, '-m pip install numpy==1.26.4"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Hmsc
  ecokit::info_chunk("Checking Hmsc")
  system2(python, '-c "import hmsc"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  return(invisible(NULL))
}
